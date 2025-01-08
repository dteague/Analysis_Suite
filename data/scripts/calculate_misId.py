#!/usr/bin/env python3
import argparse
from boost_histogram import axis
from scipy.optimize import minimize
import numpy as np
import pickle
# from matplotlib import colors

import analysis_suite.commons.configs as config
from analysis_suite.commons.histogram import Histogram
# from analysis_suite.commons.plot_utils import cms_label, plot, plot_colorbar
# from analysis_suite.commons.constants import lumi
# import analysis_suite.data.plotInfo.misId_fakerate as pinfo
from analysis_suite.plotting.plotter2 import Plotter
from analysis_suite.commons.user import workspace_area

cm_ptbins = axis.Variable([20, 30, 40, 55, 100, 150])
cm_etabins = axis.Variable([0.0, 1.1, 1.479, 1.653, 2.5])

latex_chan = {"Electron": "e", "Muon": "\mu",
              "EE": 'ee', "EM": 'e\mu', 'MM': '\mu\mu'}

formatter = {'extra_format': 'pdf',}

def plot_project(filename, tight, loose, axis_label, lumi, axis=None):
    with plot(filename, **formatter) as ax:
        if axis is None:
            eff_hist = Histogram.efficiency(tight, loose)
        else:
            eff_hist = Histogram.efficiency(tight.project(axis), loose.project(axis))
        eff_hist.plot_points(ax)
        ax.set_xlabel(axis_label)
        cms_label(ax, lumi=lumi)


def fr_plot(name, hist, year, **kwargs):
    with plot(name, subplot_info={'figsize': (18, 10)}, **formatter) as ax:
        xx = np.tile(hist.axes[0].edges, (len(hist.axes[1])+1, 1))
        yy = np.tile(hist.axes[1].edges, (len(hist.axes[0])+1, 1)).T
        mesh = ax.pcolormesh(xx, yy, hist.vals.T, norm=colors.LogNorm(vmin=1e-5, vmax=1e-2),
                             shading='flat', cmap="Blues", **kwargs)
        for j, y in enumerate(hist.axes[1].centers):
            for i, x in enumerate(hist.axes[0].centers):
                exp = int(np.floor(np.log10(hist.vals[i,j])))
                val_str = f'{hist.vals[i,j]:.2e}\n$\pm${hist.err[i,j]/10**exp:.2f}e{str(exp).zfill(3)}'
                text = ax.text(x, y, val_str, fontsize='xx-small', ha="center", va='center')
        ax.set_title("Electron Charge MisId Rate")
        ax.set_xlabel("$p_{T}(e)$ [GeV]")
        ax.set_ylabel("$|\eta(e)|$")
        cms_label(ax, year=year)
        plot_colorbar(mesh, ax, barpercent=2)

def get_fake_rate(part, fake_rate, idx):
    pt_axis, eta_axis = fake_rate.axes
    npt, neta = fake_rate.axes.size

    ptbin = np.digitize(part['pt', idx], pt_axis.edges) - 1
    ptbin = np.where(ptbin >= npt, npt-1, ptbin)
    etabin = np.digitize(part.abseta(idx), eta_axis.edges) - 1
    return fake_rate.vals.flatten()[etabin + neta*ptbin]

def scale_misId(vg, fake_rate):
    part = vg.TightElectron
    fr1 = get_fake_rate(part, fake_rate, 0)
    fr2 = get_fake_rate(part, fake_rate, 1)
    vg.scale = (fr1/(1-fr1)+fr2/(1-fr2))

def fit_template(data, flip):
    def chi2(factors, data, flip):
        mc = flip.hist*factors[0]
        tot_diff2 = np.where(data.vals < 0.1, 0, (data.vals - mc.view().value)**2)
        return np.sum(tot_diff2/(data.err+1e-6)**2)

    start_val = np.sum(data.vals)/np.sum(flip.vals)
    print(start_val)
    res = minimize(chi2, (start_val), args=(data, flip), method='Nelder-Mead')
    return res.x


def measurement(workdir, year, input_dir):
    ntuple = config.get_ntuple('charge_misId', 'measurement')
    plotter = Plotter(ntuple, year, bkg = ["DY_ht", "ttbar_lep", 'vv_inc'],
                      outdir=workdir/f'MR_{year}')
    hist_factory = HistGetter(ntuple, year=year)

    graphs = {
        'pt_lead':  GraphInfo('$p_{{T}}({}_{{lead}})$', axis.Regular(20, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', 0)),
        'pt_sub':   GraphInfo('$p_{{T}}({}_{{sub}})$', axis.Regular(20, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', 1)),
        'pt_all':   GraphInfo('$p_{{T}}(\ell)$', axis.Regular(20, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', -1)),
        'eta_lead': GraphInfo('$\eta({}_{{lead}})$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 0)),
        'eta_sub':  GraphInfo('$\eta({}_{{sub}})$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 1)),
        'eta_all':  GraphInfo('$\eta(\ell)$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', -1)),
        'mass':     GraphInfo('$M({})$', axis.Regular(30, 0, 400), lambda vg : (vg.dimass("TightLepton", 0, "TightLepton", 1), vg.scale)),
        'ht':       GraphInfo('HT', axis.Regular(30, 250, 1000), lambda vg : vg.get_hist("HT")),
        'met':      GraphInfo('MET', axis.Regular(30, 25, 250), lambda vg : vg.get_hist("Met")),
    }
    denom_fr = GraphInfo('pteta', (cm_ptbins, cm_etabins), lambda vg : vg['TightElectron'].get_hist2d('pt', 'abseta', -1)),
    num_fr = GraphInfo('pteta', (cm_ptbins, cm_etabins), lambda vg : fake_chargeMisId(vg)),

    chans = ['MM', 'EE', 'EM']
    for chan in chans:
        latex = latex_chan[chan]
        plotter.mask(lambda vg : vg["TightElectron"].num() == chan.count('E'))
        for graph in graphs:
            hists = hist_factory.get_hists(graph)
            plotter.plot_hists(hists)

    plotter.mask(lambda vg : vg["TightElectron"].num() > 0)
    plotter.fill_hists(graphs)
    for graph in graphs_1d:
        plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_e.png', chan="e\ell", region="$MR(e\ell)$", **formatter)

    denom_hist = hist_factory.get_hist(denom_fr)
    num_hist = hist_factory.get_hist(num_fr)

    fr = Histogram.efficiency(num_hist, denom_hist)
    fr_plot(plot_dir/f'fr_{year}', fr, year)

    plot_project(plot_dir/f'fr_pt.png', flip_fr, all_fr, "$p_{{T}}(e)$", lumi[year], axis=0)
    plot_project(plot_dir/f'fr_eta.png', flip_fr, all_fr, '$\eta(e)$', lumi[year], axis=1)

    # Dump fake rate
    with open(workdir/f"charge_misid_rate_{year}.pickle", "wb") as f:
        pickle.dump(fr, f)

    #################################
    plot_dir = workdir / f'MR_{year}'
    plot_dir.mkdir(exist_ok=True)

    bkg = ["DY_ht", "ttbar_lep", 'vv_inc']
    chans = ['MM', 'EE', 'EM']

    ntuple = config.get_ntuple('charge_misId', 'measurement')
    groups = ntuple.get_info().setup_groups(["data"] + bkg)
    filename = ntuple.get_filename(year=year, workdir=input_dir)

    graphs = pinfo.charge_misId['Measurement']
    graphs_1d = [graph for graph in graphs if graph.dim() == 1]

    plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
    plotter.set_groups(bkg=bkg)

    for chan in chans:
        latex = latex_chan[chan]
        plotter.mask(lambda vg : vg["TightElectron"].num() == chan.count('E'))
        plotter.fill_hists(graphs)
        for graph in graphs_1d:
            plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{chan}.png', chan=latex, region=f"$MR({latex})$", **formatter)

    plotter.mask(lambda vg : vg["TightElectron"].num() > 0)
    plotter.fill_hists(graphs)
    for graph in graphs_1d:
        plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_e.png', chan="e\ell", region="$MR(e\ell)$", **formatter)

    all_fr = plotter.get_sum(bkg, 'all_fr')
    flip_fr = plotter.get_sum(bkg, 'flip_fr')

    fr = Histogram.efficiency(flip_fr, all_fr)
    fr_plot(plot_dir/f'fr_{year}', fr, year)

    plot_project(plot_dir/f'fr_pt.png', flip_fr, all_fr, "$p_{{T}}(e)$", lumi[year], axis=0)
    plot_project(plot_dir/f'fr_eta.png', flip_fr, all_fr, '$\eta(e)$', lumi[year], axis=1)

    # Dump fake rate
    with open(workdir/f"charge_misid_rate_{year}.pickle", "wb") as f:
        pickle.dump(fr, f)


def closure(workdir, year, input_dir):
    plot_dir = workdir / f'CR_{year}'
    plot_dir.mkdir(exist_ok=True)

    mc_bkg = ["ttbar_lep", 'vv_inc', 'DY']
    graphs = pinfo.charge_misId['Closure']

    # # Opposite sign region
    # ntuple_os = config.get_ntuple('charge_misId', 'closure_os')
    # filename = ntuple_os.get_filename(year=year, workdir=input_dir)
    # plotter_os = Plotter(filename, groups, ntuple=ntuple_os, year=year)
    # plotter_os.cut(lambda vg : vg["Met"] < 50)
    # plotter_os.cut(lambda vg : vg["Nloose_Muon"]+vg["Nloose_Electron"]==2)
    # plotter_os.set_groups(bkg=mc_bkg)

    # plotter_os.fill_hists(graphs, ginfo)
    # for graph in graphs:
    #     plotter_os.plot_stack(graph.name, plot_dir/f'{graph.name}_OS.png', chan='ee', region="$OS({})$")

    # Same sign closure test


    ntuple_ss = config.get_ntuple('charge_misId', 'closure_ss')
    groups = ntuple_ss.get_info(keep_dd_data=True).setup_groups(["charge_flip", 'data'] + mc_bkg)
    filename = ntuple_ss.get_filename(year=year, workdir=input_dir)

    plotter_ss = Plotter(filename, groups, ntuple=ntuple_ss, year=year)
    plotter_ss.cut(lambda vg : vg["Met"] < 50)
    plotter_ss.set_groups(bkg=mc_bkg)

    # Scale Fake Rate
    with open(workdir/f"charge_misid_rate_{year}.pickle", "rb") as f:
        fake_rates = pickle.load(f)
    plotter_ss.scale(lambda vg : scale_misId(vg, fake_rates), groups='charge_flip')
    plotter_ss.fill_hists(graphs)

    # Get Scaling
    integrals = plotter_ss.get_integral()
    print(integrals)
    data_tot = integrals.pop("data")
    flip_tot = integrals.pop("charge_flip")
    mc_tot = sum(integrals.values())

    misID_scale = data_tot/flip_tot
    misID_mcscale = data_tot/mc_tot
    print("data/fake ratio", misID_scale)
    print("data/mc ratio", misID_mcscale)
    eta_hist = plotter_ss.get_hists("eta_all")
    fit_scale = fit_template(eta_hist['data'], eta_hist['charge_flip'])[0]
    print("template_fit", fit_scale)

    for graph in graphs:
        plotter_ss.plot_stack(graph.name, plot_dir/f'{graph.name}_SS_mc.png', chan='ee', region="$SS({})$", **formatter)

    plotter_ss.set_groups(bkg=mc_bkg, sig='charge_flip')
    plotter_ss.scale_hists('charge_flip', fit_scale)
    for graph in graphs:
        plotter_ss.plot_stack(graph.name, plot_dir/f'{graph.name}_SS_all.png', chan='ee', region="$SS({})$", **formatter)
    # Plot with just Data-Driven Background
    plotter_ss.set_groups(bkg=['charge_flip'])
    for graph in graphs:
        plotter_ss.plot_stack(graph.name, plot_dir/f'{graph.name}_SS_data.png', chan='ee', region="$SS({})$", **formatter)

    fake_rates *= fit_scale
    with open(workdir/f"fr_{year}.pickle", "wb") as f:
        pickle.dump(fake_rates, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-d', '--workdir', help="directory to run over. If nothing, use date",)
    parser.add_argument('-i', '--input_dir', default=None)
    parser.add_argument('-r', '--run', type=lambda x: [i.strip() for i in x.split(',')],
                        help="Regions to run through (sideband, measurement, closure)")
    args = parser.parse_args()

    workdir = workspace_area/ args.workdir / 'charge_misId'
    workdir.mkdir(exist_ok=True, parents=True)

    for year in args.years:
        if 'measurement' in args.run:
            print("Measurement")
            measurement(workdir, year, args.input_dir)
        if 'closure' in args.run:
            print("Closure")
            closure(workdir, year, args.input_dir)
