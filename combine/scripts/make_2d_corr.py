#!/usr/bin/env python3
from scipy.interpolate import griddata
import numpy as np
import argparse
import uproot
import boost_histogram as bh
import subprocess
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use([hep.style.CMS])
import ROOT

xsec_signal = 2.05
xsec_4top = 13.37

colors = ["#F7CE76", "#FFFFA6", "#FFFFFF"]
mycmap = mpl.colors.LinearSegmentedColormap.from_list("mycmap", colors)

class Data2D:
    def __init__(self, data_file, corr_file, n_points=100, x_range=[0.,3.], y_range=[0.,3.]):
        with uproot.open(data_file) as f:
            arr = f["limit"].arrays()
            self.y, self.x = arr.r_sig, arr.r_other
            self.dnll = arr.deltaNLL
        self.corr = np.nan
        if "Hesse" in corr_file:
            with uproot.open(corr_file) as f:
                self.corr = f['h_correlation'].to_boost()[bh.loc("r_sig"), bh.loc("r_other")]
        else:
            f =  ROOT.TFile.Open(str(corr_file))
            fit = f.Get("fit_mdf")
            self.corr = fit.correlation("r_sig", "r_other")
            f.Close()
        # except:
        #     pass
        self.n_points = n_points
        self.xlo, self.xhi = min(self.x), max(self.x)
        self.ylo, self.yhi = min(self.y), max(self.y)


    def get_grid(self):
        points = np.array([self.x, self.y]).transpose()
        grid_x, grid_y = np.mgrid[self.xlo:self.xhi:self.n_points*1j, self.ylo:self.yhi:self.n_points*1j]
        grid_vals = griddata(points, self.dnll, (grid_x,grid_y), "cubic")
        grid_vals = np.nan_to_num(grid_vals, nan=100)
        return (grid_x, grid_y, grid_vals)

def run_command(command, debug):
    if debug:
        print(command)
        print('-'*40)
    output = subprocess.run(command, shell=True, stdout = subprocess.PIPE,
                            stderr=subprocess.PIPE)
    if output.returncode != 0 or debug:
        print(output.stdout.decode())
        print(output.stderr.decode())
    if output.returncode != 0:
        exit()

def run_combine(cardname, year, other, sig, sig_range=[0,100], other_range=[0,10], hesse=False, debug=False, skip=False, blind=True):
    corr_command = '--robustHesse 1' if hesse else "--robustFit 1"
    corr_filename = "robustHesse" if hesse else "multidimfit"
    typ = 'exp' if blind else 'obs'
    out_corr, out_plot = f"corr_{typ}", f"plot_{typ}"
    outname = f'scan2d.{sig}-{other}.{year}'

    output_files = f'higgsCombine.{out_plot}.{outname}.MultiDimFit.mH125.root', f'{corr_filename}.{out_corr}.{outname}.root'

    if skip:
        return output_files

    sig_map = ""
    if isinstance(sig, list):
        for s in sig:
            sig_map += f'--PO "map=.*/{s}:r_sig[1,-200,200]" '
    else:
        sig_map = f'--PO "map=.*/{sig}:r_sig[1,-200,200]" '
    r_sig = ",".join(map(str, sig_range))
    r_other = ",".join(map(str, other_range))

    if 'txt' in cardname:
        run_command(f'text2workspace.py {cardname} -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel '+
                    f'{sig_map} --PO "map=.*/{other}:r_other[1,-200,200]" -o {cardname.replace("txt", "root")}', debug=debug)
        cardname = cardname.replace('txt', 'root')

    blind_args = '-t -1 --expectSignal=1' if blind else ""
    n = 30
    split = 30
    width = n**2//split
    cmd = []
    for i in range(split):
        start = i*width
        end = (i+1)*width
        cmd.append(f'combine -M MultiDimFit {cardname} -m 125 -n .{out_plot}.{outname}_{i} --gridPoints {n},{n} --algo grid '+
                   f'--cminDefaultMinimizerStrategy 0 -P r_sig -P r_other --setParameterRanges r_sig={r_sig}:r_other={r_other} ' +
                   f'--alignEdges 1 --firstPoint {start} --lastPoint {end} {blind_args}')
    print(cmd[0])

    if debug:
        stdout = None
        stderr = subprocess.STDOUT
    else:
        stdout = subprocess.DEVNULL
        stderr = subprocess.DEVNULL
    processes = [subprocess.Popen(i, shell=True, stdout=stdout, stderr=stderr) for i in cmd]
    for p in processes: p.wait()
    run_command(f'hadd -f higgsCombine.{out_plot}.{outname}.MultiDimFit.mH125.root higgsCombine.{out_plot}.{outname}_*.MultiDimFit.mH125.root',
                debug=debug)
    run_command(f'combine -M MultiDimFit {cardname} -m 125 -n .{out_corr}.{outname} {blind_args} ' +
                f'--cminDefaultMinimizerStrategy 0 -P r_sig -P r_other --setParameterRanges r_sig={r_sig}:r_other={r_other} ' +
                f'{corr_command} --saveFitResult --stepSize=0.01 --setRobustFitTolerance=1000 --setCrossingTolerance=0.001 ',
                debug=debug)

    return output_files

def make_2d(outfile, data, blind=False, use_xsec=True, **kwargs):
    cmap = mpl.colormaps.get_cmap('cividis_r')
    r_sig, r_other = (xsec_signal, xsec_4top) if use_xsec else (1, 1)

    fig, ax = plt.subplots(figsize=(13,13))

    grid_x, grid_y, grid_vals = data.get_grid()
    grid_x *= r_other
    grid_y *= r_sig
    min_idx = np.argmin(np.nan_to_num(grid_vals, nan=np.nanmax(grid_vals)).flatten())

    mesh = ax.pcolormesh(grid_x, grid_y, 2*grid_vals, cmap=mycmap, vmax=10)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.01)
    ax.contour(grid_x, grid_y, 2*grid_vals, [2.3], linewidths=3, colors='k')
    ax.contour(grid_x, grid_y, 2*grid_vals, [5.99], linewidths=3, colors='k',
               linestyles='dashed')
    sm_pt = ax.plot([r_other], [r_sig], markersize=15, color='r', marker='d', label='SM', ls='')
    bf_pt = ax.plot([grid_x.flatten()[min_idx]], [grid_y.flatten()[min_idx]], ls='',
                    markersize=15, color='k', marker='P', label='Best fit',)
    legend_elements = [
        sm_pt[0],
        plt.Line2D([0], [0], color='k', lw=3, label=r"$1\sigma$ CL"),
        plt.Line2D([0], [0], color='k', lw=3, label=r"$2\sigma$ CL", linestyle='dashed'),
        bf_pt[0],
    ]
    if use_xsec:
        ax.set_xlabel(r"$\sigma_{tttt}$ [fb]", loc='right')
        ax.set_ylabel(r"$\sigma_{ttt}$ [fb]", loc='top')
    else:
        ax.set_xlabel(r"$r_{TTTT}$", loc='right')
        ax.set_ylabel(r"$r_{SIG}$", loc='top')
    ax.text(0.05, 0.9, rf'$\rho = {data.corr:0.3f}$', transform=ax.transAxes, size='large')
    ax.legend(handles=legend_elements, facecolor='white', framealpha=1, edgecolor='k',
              frameon=True, fancybox=True)
    hep.cms.label(ax=ax, label="Preliminary", lumi="137.0", data=not blind)
    plt.savefig(outfile, bbox_inches='tight', dpi=300)
    print(data.corr)

def make_2d_expobs(outfile, exp, obs, use_xsec=True, **kwargs):
    cmap = mpl.colormaps.get_cmap('cividis_r')
    r_sig, r_other = (xsec_signal, xsec_4top) if use_xsec else (1, 1)

    fig, ax = plt.subplots(figsize=(13,13))

    exp_x, exp_y, exp_vals = exp.get_grid()
    exp_x *= r_other
    exp_y *= r_sig
    obs_x, obs_y, obs_vals = obs.get_grid()
    obs_x *= r_other
    obs_y *= r_sig
    min_exp = np.argmin(np.nan_to_num(exp_vals, nan=np.nanmax(exp_vals)).flatten())
    min_obs = np.argmin(np.nan_to_num(obs_vals, nan=np.nanmax(obs_vals)).flatten())
    print(np.max(obs_vals))

    mesh = ax.pcolormesh(obs_x, obs_y, 2*obs_vals, cmap=mycmap, vmax=10)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.01)
    cbar.set_label(r"Obs: $-2\Delta\ln(L)$")
    # ax.pcolormesh(exp_x, exp_y, 2*exp_vals, cmap=mycmap)
    ax.contour(exp_x, exp_y, 2*exp_vals, [2.3], colors='k', linewidths=3, linestyles='dashed')
    ax.contour(exp_x, exp_y, 2*exp_vals, [5.99], colors='r', linewidths=3,linestyles='dashed')
    ax.contour(obs_x, obs_y, 2*obs_vals, [2.3], colors='k', linewidths=3,)
    ax.contour(obs_x, obs_y, 2*obs_vals, [5.99], colors='r', linewidths=3,)
    sm_pt = ax.plot([r_other], [r_sig], markersize=15, color='r', marker='d', label='SM', ls='')
    bf_pt = ax.plot([obs_x.flatten()[min_obs]], [obs_y.flatten()[min_obs]], ls='',
                    markersize=15, color='k', marker='P', label='Best fit',)
    legend_elements = [
        sm_pt[0],
        plt.Line2D([0], [0], color='k', lw=3, label=r"$1\sigma$ Obs CL"),
        plt.Line2D([0], [0], color='r', lw=3, label=r"$2\sigma$ Obs CL"),
        bf_pt[0],
        plt.Line2D([0], [0], color='k', lw=3, label=r"$1\sigma$ Exp CL", linestyle='dashed'),
        plt.Line2D([0], [0], color='r', lw=3, label=r"$2\sigma$ Exp CL", linestyle='dashed'),
    ]
    if use_xsec:
        ax.set_xlabel(r"$\sigma_{tttt}$ [fb]", loc='right')
        ax.set_ylabel(r"$\sigma_{ttt}$ [fb]", loc='top')
    else:
        ax.set_xlabel(r"$r_{TTTT}$", loc='right')
        ax.set_ylabel(r"$r_{SIG}$", loc='top')
    # ax.text(0.05, 0.9, rf'$\rho = {data.corr:0.3f}$', transform=ax.transAxes, size='large')
    ax.legend(handles=legend_elements, facecolor='white', framealpha=1, edgecolor='k',
              frameon=True, fancybox=True, ncols=2)
    hep.cms.label(ax=ax, label="Preliminary", lumi="137.0", data=True)
    print(outfile)
    plt.savefig(outfile, bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="main", description="")
    parser.add_argument('-c', '--card', required=True)
    parser.add_argument('-o', '--other', default='TTTT')
    parser.add_argument('-s', '--signal', default='SIG')
    parser.add_argument('--hesse', action="store_true")
    parser.add_argument('-b', '--unblind', action="store_true")
    parser.add_argument('-sr', '--sig_range', nargs=2, type=float)
    parser.add_argument('-or', '--other_range', nargs=2, type=float)
    parser.add_argument("--skip", action="store_true")
    parser.add_argument("--debug", action='store_true')
    args = parser.parse_args()

    # Choose year. If fails to find, it defaults to the last value (might be a problem!)
    short_cardname = args.card[:args.card.index('.')]
    all_eras = ['2016pre', '2016post', '2017', '2018', 'all']
    for year in all_eras:
        if year in args.card:
            short_cardname = args.card[:args.card.index(year)-1]
            break

    data_file, corr_file = run_combine(
        args.card, year, args.other, sig=args.signal, sig_range=args.sig_range,
        other_range=args.other_range, debug=args.debug, hesse=args.hesse, skip=args.skip,
        blind=True
    )
    exp_data = Data2D(data_file, corr_file, y_range=args.sig_range, x_range=args.other_range)
    make_2d(f"correlation_exp_{short_cardname}_{args.signal}-{args.other}_{year}.png", exp_data)
    if args.unblind:
        data_file, corr_file = run_combine(
            args.card, year, args.other, sig=args.signal, sig_range=args.sig_range,
            other_range=args.other_range, debug=args.debug, hesse=args.hesse, skip=args.skip,
            blind=False
        )
        obs_data = Data2D(data_file, corr_file, y_range=args.sig_range, x_range=args.other_range)
        make_2d(f"correlation_obs_{short_cardname}_{args.signal}-{args.other}_{year}.png", obs_data, blind=True)
        make_2d_expobs(f"correlation_{short_cardname}_{args.signal}-{args.other}_{year}.png", exp_data, obs_data)
