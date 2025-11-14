#!/usr/bin/env python3
import argparse
from boost_histogram import axis
from scipy.optimize import minimize
import numpy as np
import pickle

import analysis_suite.commons.configs as config
from analysis_suite.commons.histogram import Histogram
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
    chans = ['MM', 'EE', 'EM']
    ntuple = config.get_ntuple('charge_misId', 'measurement')
    plotter = Plotter(ntuple, year, bkg = ["DY_ht", "ttbar_lep", 'vv_inc'],
                      outdir=workdir/f'MR_{year}')
    hist_factory = HistGetter(ntuple, year=year)

    graphs = {
        'pt_lead':  GraphInfo('$p_{{T}}(\ell_{{lead}})$', axis.Regular(20, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', 0)),
        'pt_sub':   GraphInfo('$p_{{T}}(\ell_{{sub}})$', axis.Regular(20, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', 1)),
        'pt_all':   GraphInfo('$p_{{T}}(\ell)$', axis.Regular(20, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', -1)),
        'eta_lead': GraphInfo('$\eta(\ell_{{lead}})$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 0)),
        'eta_sub':  GraphInfo('$\eta(\ell_{{sub}})$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 1)),
        'eta_all':  GraphInfo('$\eta(\ell)$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', -1)),
        'mass':     GraphInfo('$M(\ell,\ell)$', axis.Regular(30, 0, 400), lambda vg : (vg.dimass("TightLepton", 0, "TightLepton", 1), vg.scale)),
        'ht':       GraphInfo('H_{{T}}', axis.Regular(30, 250, 1000), lambda vg : vg.get_hist("HT")),
        'met':      GraphInfo('MET', axis.Regular(30, 25, 250), lambda vg : vg.get_hist("Met")),
    }
    denom_fr = GraphInfo('pteta', (cm_ptbins, cm_etabins), lambda vg : vg['TightElectron'].get_hist2d('pt', 'abseta', -1)),
    num_fr = GraphInfo('pteta', (cm_ptbins, cm_etabins), lambda vg : fake_chargeMisId(vg)),

    for chan in chans:
        hist_factory.reset_mask()
        hist_factory.mask(lambda vg : vg["TightElectron"].num() == chan.count('E'))
        plotter.set_extra_text(chan)
        hists = hist_factory.get_hists(graphs)
        plotter.plot_hists(hists)

    # Any event with Electron
    hist_factory.reset_mask()
    hist_factory.mask(lambda vg : vg["TightElectron"].num() > 0)
    plotter.set_extra_text("anyE")
    hists = hist_factory.get_hists(graphs)
    plotter.plot_hists(hists)

    denom_hist = hist_factory.get_hist(denom_fr)
    num_hist = hist_factory.get_hist(num_fr)

    fr = Histogram.efficiency(num_hist, denom_hist)
    fr_plot(plot_dir/f'fr_{year}', fr, year)

    plot_project(plot_dir/f'fr_pt.png', flip_fr, all_fr, "$p_{{T}}(e)$", lumi[year], axis=0)
    plot_project(plot_dir/f'fr_eta.png', flip_fr, all_fr, '$\eta(e)$', lumi[year], axis=1)

    # Dump fake rate
    with open(workdir/f"charge_misid_rate_{year}.pickle", "wb") as f:
        pickle.dump(fr, f)

def closure(workdir, year, input_dir):
    chans = ['MM', 'EE', 'EM']
    ntuple = config.get_ntuple('charge_misId', 'closure_ss')
    plotter = Plotter(ntuple, year, bkg=['charge_flip'], outdir=workdir/f'CR_{year}')
    plotter_mc = Plotter(ntuple, year, data='charge_flip',
                         bkg=["DY", "ttbar_lep", 'vv_inc'], outdir=workdir/f'CR_mc_{year}')

    # Scale Fake Rate
    with open(workdir/f"charge_misid_rate_{year}.pickle", "rb") as f:
        fake_rates = pickle.load(f)

    graphs = {
        "pt_lead":   GraphInfo('$p_{{T}}({}_{{lead}})$', axis.Regular(30, 0, 150), lambda vg : vg['TightLepton'].get_hist('pt', 0)),
        "pt_sub":    GraphInfo('$p_{{T}}({}_{{sub}})$', axis.Regular(30, 0, 150), lambda vg : vg['TightLepton'].get_hist('pt', 1)),
        "pt_all":    GraphInfo('$p_{{T}}(\ell)$', axis.Regular(30, 0, 150), lambda vg : vg['TightLepton'].get_hist('pt', -1)),
        "eta_lead":  GraphInfo('$\eta({}_{{lead}})$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 0)),
        "eta_sub":   GraphInfo('$\eta({}_{{sub}})$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 1)),
        "eta_all":   GraphInfo('$\eta(\ell)$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', -1)),
        "mass":      GraphInfo('$M({})$', axis.Regular(16, 80, 100), lambda vg : (vg.dimass("TightLepton", 0, "TightLepton", 1), vg.scale)),
        "ht":        GraphInfo('HT', axis.Regular(30, 0, 250), lambda vg : vg.get_hist("HT")),
        "met":       GraphInfo('MET', axis.Regular(25, 0, 50), lambda vg : vg.get_hist("Met")),
        "metphi":    GraphInfo('$\phi(MET)$', axis.Regular(20, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
    }

    # Same sign closure test
    hist_factory = HistGetter(ntuple, year=year)
    hist_factory.scale(lambda vg : scale_misId(vg, fake_rates), groups='charge_flip')
    hists = hist_factory.fill_hists(graphs)
    plotter.plot_hists(hists)
    plotter_mc.plot_hists(hists)

    # Get Scaling
    data, flip, mc = 0, 0, 0
    for member, hist in next(iter(hists.values())).items():
        if member == "data":
            data += hist.integral()
        elif member == "charge_flip":
            flip += hist.integral()
        else:
            mc += hist.integral()
        print(member, hist.integral())

    print("data/fake ratio", data/flip)
    print("data/mc ratio", data/mc)
    eta_hist = hists("eta_all")
    fit_scale = fit_template(eta_hist['data'], eta_hist['charge_flip'])[0]
    print("template_fit", fit_scale)

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
