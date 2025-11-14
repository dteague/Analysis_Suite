#!/usr/bin/env python3
import argparse
from scipy.optimize import minimize
import numpy as np
import pickle
import awkward as ak
import boost_histogram.axis as axis
from datetime import datetime
import operator
from copy import copy

from analysis_suite.plotting.plotter import Plotter
from analysis_suite.plotting.hist_getter import GraphInfo, HistGetter
import analysis_suite.commons.configs as config
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import plot, plot_colorbar, ratio_plot
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.info import GroupInfo
from analysis_suite.commons.user import workspace_area
from scipy import stats

np_ptbins = {
    "Electron": axis.Variable([20, 25, 35, 45, 100]), # axis.Variable([20, 30, 45, 70]),
    "Muon": axis.Variable([20, 30, 45, 60, 100]), #axis.Variable([20, 35, 70])
}

np_etabins = {
    "Electron": axis.Variable([0.0, 0.8, 1.479, 2.5]),
    "Muon": axis.Variable([0.0, 1.2, 2.1, 2.4]),
}

def get_fake_rate(part, fake_rate, idx, name="Muon"):
    if isinstance(fake_rate, np.ndarray):
        pt_axis = np_ptbins[name]
        eta_axis = np_etabins[name]
        npt, neta = pt_axis.size, eta_axis.size
    else:
        pt_axis, eta_axis = fake_rate.axes
        npt, neta = fake_rate.axes.size
        fake_rate = fake_rate.values()
    ptbin = np.digitize(part['pt', idx], pt_axis.edges) - 1
    ptbin = np.where(ptbin >= npt, npt-1, ptbin)
    etabin = np.digitize(part.abseta(idx), eta_axis.edges) - 1
    return fake_rate.flatten()[etabin + neta*ptbin]


formatter = {'extra_format': 'pdf',}

latex_chan = {"Electron": "e", "Muon": "\mu",
              "EE": 'ee', "EM": 'e\mu', 'MM': '\mu\mu'}
loose_name = {"Muon": "muon", "Electron": "elec"}
bjet_cuts = {
    "none": {"2016": -1, "2017": -1, "2018": -1},
    "loose": {"2016": 0.05, "2017": 0.0532, "2018": 0.0490},
    "medium": {"2016": 0.25, "2017": 0.3040, "2018": 0.2783},
    "tight": {"2016": 0.64, "2017": 0.7476, "2018": 0.7100},
}
trig_cuts = {"Muon": 20, "Electron": 25}

close_btag_cut = 0.1


def plot_project(filename, tight, loose, axis, axis_label, lumi, data):
    with plot(filename, lumi=lumi, data=data, **formatter) as ax:
        eff_hist = Histogram.efficiency(tight.project(axis), loose.project(axis))
        eff_hist.plot_points(ax)
        ax.set_xlabel(axis_label)
        ax.set_ylabel("Fake Rate")

def plot_fr_diff(filename, data_ewk, qcd, axis_label, lumi):
    data_ewk.set_plot_details("Data-EWK", 'k')
    qcd.set_plot_details("QCD", 'r')

    ratio = Histogram(*qcd.axes, color="black")
    ratio += data_ewk/qcd

    with ratio_plot(filename, axis_label, qcd.get_xrange(), pad_label='Fake Rate',
                    subpad_label='(Data-EWK)/QCD', lumi=lumi, data=True,
                    **formatter) as ax:
        pad, subpad = ax
        data_ewk.plot_points(pad, capsize=1)
        qcd.plot_points(pad, capsize=1)
        ratio.plot_points(subpad, data=True)


def fit_template(data, qcd, ewk):
    def chi2(factors, data, qcd, ewk):
        mc = qcd.vals*factors[0] + ewk.vals*factors[1] + 1e-6
        tot_diff2 = (data.vals - mc)**2
        tot_diff2[tot_diff2 < 0] = 0.
        return np.sum(tot_diff2/mc)

    start_val = np.sum(data.vals)/(np.sum(qcd.vals+ewk.vals))
    res = minimize(chi2, (start_val, start_val), args=(data, qcd, ewk), method='Nelder-Mead',)
    return res.x

def scale_template(vg, chan, scales):
    pt_bins = np_ptbins[chan].edges
    eta_bins = np_etabins[chan].edges
    npt, neta = len(pt_bins)-1, len(eta_bins)-1

    ptbin = np.digitize(vg[chan]['pt', 0], pt_bins) - 1
    ptbin = np.where(ptbin >= npt, npt-1, ptbin)
    etabin = np.digitize(vg[chan]['abseta', 0], eta_bins) - 1
    vg.scale = (scales.flatten()[etabin + neta*ptbin], vg[chan].num() == 1)

def scale_fake(vg, chan, fakerate, i=0):
    part = vg[f'Fake{chan}']
    fr = get_fake_rate(part, fakerate, i, chan)
    vg.scale = (fr/(1-fr), part.num() > i)

def flip_fake(vg):
    vg.scale = (-1, (vg['FakeMuon'].num() + vg['FakeElectron'].num() == 2))

def fr_plot(name, ratio, chan, year, hasData, **kwargs):
    from matplotlib import ticker
    part = latex_chan[chan]
    with plot(name, lumi=lumi[year], data=hasData, **formatter) as ax:
        mesh = ratio.plot_2d(ax, **kwargs)
        ax.set_xlabel(f"$p_{{T}}({part})$ [GeV]")
        ax.set_ylabel(f"$|\eta({part})|$")
        plot_colorbar(mesh, ax)
        ax.set_xscale('log')
        ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
        ax.get_xaxis().set_minor_formatter(ticker.ScalarFormatter())

def fix_negative(data_ewk, qcd):
    axes = data_ewk.axes
    for x in range(axes.size[0]):
        for y in range(axes.size[1]):
            if data_ewk[x, y].value <= 0 or np.isnan(data_ewk[x, y].value):
                data_ewk[x, y] = qcd[x, y]


def plot_ptcorr(filename, pt, bins, lumi, chan, fact=None):
    with plot(filename, lumi=lumi, data=True, **formatter) as ax:
        print_spot= (max(pt) + min(pt))/2
        width = (bins[1]-bins[0])/2

        ax.hist(bins[:-1], bins, weights=pt, histtype='step')
        # ax.errorbar(bins[:-1]+width, pt, yerr=err, ls="", capsize=5, ecolor='blue')
        ax.set_xlim(-1, 1)
        ax.set_ylim(0.9*min(pt))
        ax.plot([0.4, 0.4], [0, max(pt)*1.5], color='red')

        fact_text = r"$\frac{p_{T}}{p_{T}^{ratio}}$" if fact is None else f"{fact}"+r"$\times\frac{p_{T}}{p_{T}^{ratio}}$"
        ax.text(-0.35, print_spot, fact_text)
        ax.text(0.65, print_spot, r"$p_{T}$")
        ax.set_xlabel("$disc_{TTH}$")
        ax.set_ylabel(f"$mean(p_{{T}}({chan}))$")


def plot_ptcorr_2d(filename, vals, x, y, year):
    with plot(filename, lumi=lumi[year], data=True, **formatter) as ax:
        mesh = ax.pcolormesh(x, y, vals, shading='flat')
        plot_colorbar(mesh, ax)
        ax.set_xscale('log')

def chi2(h1, h2):
    s1, s2 = h1.sumw2, h2.sumw2
    W1, W2 = h1.integral(), h2.integral()
    return np.sum((W1*h2.vals-W2*h2.vals)**2/(W1**2*s2 + W2**2*s1 + 1e-5))

#--------------------------------------------------------------------------

def sideband(workdir, year, input_dir, fake_qcd):
    chans = ['Electron', 'Muon']
    ntuple = config.get_ntuple('fake_rate', 'sideband')
    plotter = Plotter(ntuple, year, bkg='all', outdir=workdir / f'SB_pteta_{year}')

    fake_name = '_fake' if fake_qcd else ''
    mc_scale_factors = {chan: {reg: list() for reg in ['ewk', 'qcd']} for chan in chans}
    round5 = lambda x : np.round(0.3*x+14)*5

    graphs = {
        "mt":     GraphInfo('$M_{{T}}({chan}_{{tight}})$', axis.Regular(25, 0, 250),
                  lambda vg : vg['AllLepton'].get_hist('mt', 0)),
        'mt_fix': GraphInfo('$M_{{T}}({chan}_{{tight}})$', axis.Regular(12, 0, 120),
                  lambda vg : vg['AllLepton'].get_hist('mt_fix', 0)),
        'pt':     GraphInfo('$p_{{T}}({chan}_{{tight}})$', axis.Regular(30, 0, 150),
                  lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        'met':    GraphInfo('$MET_{{Puppi}}$', axis.Regular(30, 0, 150),
                  lambda vg : vg.get_hist('Met')),
    }

    hist_factory = HistGetter(ntuple, year, workdir=workdir)
    if fake_qcd:
        hist_factory.cut(lambda vg : vg['FakeLepton'].num() == 1, groups=['qcd'])
    else:
        hist_factory.cut(lambda vg : vg['TightLepton'].num() == 1, groups=['qcd'])

    for chan in chans:
        pt_bins = np_ptbins[chan].edges
        eta_bins = np_etabins[chan].edges
        for pt_lo, pt_hi in zip(pt_bins, pt_bins[1:]):
            for eta_lo, eta_hi in zip(eta_bins, eta_bins[1:]):
                hist_factory.reset_mask()
                plotter.set_extra_text(f'{chan}_{int(pt_lo)}_{eta_lo}')
                print(chan, f'{pt_lo}-{pt_hi}', f'{eta_lo}-{eta_hi}')
                hist_factory.mask(lambda vg : vg[chan].num() == 1)
                hist_factory.mask(lambda vg : vg[chan]['abseta', 0] >= eta_lo)
                hist_factory.mask(lambda vg : vg[chan]['abseta', 0] < eta_hi)
                hist_factory.mask(lambda vg : vg[chan]['pt', 0] >= pt_lo)
                if pt_hi != pt_bins[-1]:
                    hist_factory.mask(lambda vg : vg[chan]['pt', 0] < pt_hi)
                graphs['mt'].bin_tuple = axis.Regular(15, 0, round5(pt_lo))
                hists = hist_factory.get_hists(graphs, chan=latex_chan[chan])

                # calculate templated fit
                qcd_f, ewk_f = fit_template(hists['mt']['data'], hists['mt']['qcd'], hists['mt']["ewk"])
                for hist in hists.values():
                    hist['ewk'].scale(ewk_f)
                    hist['qcd'].scale(qcd_f)

                mc_scale_factors[chan]['ewk'].append(ewk_f)
                mc_scale_factors[chan]['qcd'].append(qcd_f)

                plotter.plot_hists(hists)

        for name, val in mc_scale_factors[chan].items():
            print(name, val)
            mc_scale_factors[chan][name] = np.reshape(val, (len(pt_bins)-1, len(eta_bins)-1))

        scales = mc_scale_factors[chan]
        plot_ptcorr_2d(workdir / f'SB_pteta_{year}'/f'qcd_fact_{chan}{fake_name}.png', scales['qcd'], eta_bins, pt_bins, year)
        plot_ptcorr_2d(workdir / f'SB_pteta_{year}'/f'ewk_fact_{chan}{fake_name}.png', scales['ewk'], eta_bins, pt_bins, year)

    pickle.dump(mc_scale_factors, open(workdir/f"mc_scales{fake_name}_{year}.pickle", "wb"))


def measurement(workdir, year, input_dir, fake_qcd):
    chans = ['Electron', 'Muon']
    ntuple = config.get_ntuple('fake_rate', 'measurement')
    plotter = Plotter(ntuple, year, bkg='all', outdir=workdir / f'MR_{year}')

    # Load MC scale factors
    with open(workdir/f"mc_scales_{year}.pickle", "rb") as f:
        mc_scale_factors = pickle.load(f)
    if fake_qcd:
        with open(workdir/f"mc_scales_fake_{year}.pickle", "rb") as f:
            tmp = pickle.load(f)
            for chan in chans:
                mc_scale_factors[chan]['ewk'] = tmp[chan]['ewk']
    fake_name = '_fake' if fake_qcd else ""

    graphs = {
        'pt':     GraphInfo('$p_{{T}}({chan})$', axis.Regular(20, 0, 100), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        'mt':     GraphInfo('$M_{{T}}({chan})$', axis.Regular(30, 0, 150), lambda vg : vg['AllLepton'].get_hist('mt', 0)),
        'met':    GraphInfo('$MET_{{Puppi}}$', axis.Regular(20, 0, 50), lambda vg : vg.get_hist("Met")),
        'ht':     GraphInfo('$H_T$', axis.Regular(30, 0, 300), lambda vg : vg.get_hist("HT")),
        'metphi': GraphInfo('$\phi(MET)$', axis.Regular(20, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
        'njets':  GraphInfo('$N_j$', axis.Regular(6, 0, 6), lambda vg : (vg.Jets.num(), vg.scale)),
        'iso':    GraphInfo('$iso({chan})$', axis.Regular(20, 0, 0.4), lambda vg : vg['AllLepton'].get_hist('iso', 0)),
        'j_btag': GraphInfo("", axis.Regular(25, 0, 1), lambda vg : vg["Jets"].get_hist("discriminator", 0)),
    }

    hist_factory = HistGetter(ntuple, year, workdir=workdir, inputdir=input_dir)

    for chan in chans:
        for group, member, df in hist_factory.df_iter():
            if group not in ['ewk', 'qcd']:
                continue
            scale_template(df, chan, mc_scale_factors[chan][group])

    if year == '2016':
        bcuts = {"Muon": "loose", "Electron": "loose"}
    if year == '2017':
        bcuts = {"Muon": "loose", "Electron": "loose"}
    elif year == '2018':
        bcuts = {"Muon": "loose", "Electron": "loose"}

    fake_rates = {chan: dict() for chan in chans}

    for chan in chans:
        hist_factory.reset_mask()
        print(chan)
        latex = latex_chan[chan]
        np_bins = (np_ptbins[chan], np_etabins[chan])
        fr_graph = {
            'tight': GraphInfo('', np_bins, lambda vg : vg['TightLepton'].get_hist2d('pt', 'abseta', 0)),
            'loose': GraphInfo('', np_bins, lambda vg : vg['AllLepton'].get_hist2d('pt', 'abseta', 0)),
        }

        hist_factory.mask(lambda vg : vg[chan].num() == 1)
        if bcuts[chan] != 'none':
            hist_factory.mask(lambda vg: vg[f"N_b{bcuts[chan]}"] > 0)
            # hist_factory.mask(lambda vg: ak.count_nonzero(vg.Jets['discriminator', -1] > bjet_cuts[bcuts[chan]][year], axis=-1) > 0)

        hists = hist_factory.get_hists(graphs, chan=latex)
        plotter.set_extra_text(f'{chan}')
        plotter.plot_hists(hists)

        fr = hist_factory.get_hists(fr_graph)
        for name, subhist in fr.items():
            fr[name]['data_ewk'] = subhist['data'] - subhist['ewk']

        fr_pt = dict()
        for typ in ['qcd', 'ewk', 'data_ewk']:
            fake_rates[chan][typ] = Histogram.efficiency(fr['tight'][typ], fr['loose'][typ])
            fr_pt[typ] = Histogram.efficiency(fr['tight'][typ].project(0), fr['loose'][typ].project(0))

        fix_negative(fake_rates[chan]['data_ewk'], fake_rates[chan]['qcd'])
        plot_fr_diff(workdir/f'MR_{year}'/f'fr_{chan}_diff_pt{fake_name}.pdf', fr_pt['data_ewk'], fr_pt['qcd'], f"$p_{{T}}({latex})$", lumi[year])

        for key, fr_ in fake_rates[chan].items():
            plot_project(workdir/f'MR_{year}'/f'fr_{key}_{chan}_pt{fake_name}.pdf', fr['tight'][key], fr['loose'][key], 0, f"$p_{{T}}({latex})$", lumi[year], 'data' in key)
            plot_project(workdir/f'MR_{year}'/f'fr_{key}_{chan}_eta{fake_name}.pdf', fr['tight'][key], fr['loose'][key], 1, f"$p_{{T}}({latex})$", lumi[year], 'data' in key)
            fr_plot(workdir/f'MR_{year}'/f'fr_{key}_{chan}{fake_name}.pdf', fr_, chan, year, 'data' in key, vmin=0, vmax=0.4)
            print(key, "fr", np.array2string(fr_.vals, separator=','))

        hist_factory.mask(lambda vg : vg[f'Tight{chan}'].num() == 1)
        hists = hist_factory.get_hists(graphs, chan=latex)
        plotter.set_extra_text(f"tight_{chan}")
        plotter.plot_hists(hists)

    with open(workdir/f"fr_{year}_new.pickle", "wb") as f:
        pickle.dump(fake_rates, f)

def closure_tf(workdir, year, input_dir, fake_qcd):
    chans = ['MM', 'EE', 'EM']
    ntuple = config.get_ntuple('fake_rate', 'closure_tf')
    plotter = Plotter(ntuple, year, bkg=['ttbar_lep', 'wjet_ht'], outdir=workdir / f'CR_TF_{year}')

    graphs = {
        'pt1':    GraphInfo('$p_{{T}}(\ell_{{1}})$', axis.Regular(20, 0, 250), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        'pt2':    GraphInfo('$p_{{T}}(\ell_{{2}})$', axis.Regular(24, 0, 120), lambda vg : vg['AllLepton'].get_hist('pt', 1)),
        'eta1':   GraphInfo('$\eta(\ell_{{1}})$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),
        'eta2':   GraphInfo('$\eta(\ell_{{2}})$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 1)),
        'pt':     GraphInfo('$p_{{T}}(\ell)$', axis.Regular(20, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', -1)),
        'eta':    GraphInfo('$\eta(\ell)$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),
        'met':    GraphInfo('$MET$', axis.Regular(20, 0, 300), lambda vg : vg.get_hist("Met")),
        'ht':     GraphInfo('$H_T$', axis.Regular(20, 0, 750), lambda vg : vg.get_hist("HT")),
        'ptsum':     GraphInfo('$H_T$', axis.Regular(20, 0, 750), lambda vg : (ht(vg), vg.scale)),
        'metphi': GraphInfo('$\phi(MET)$', axis.Regular(20, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
        # 'nbjets': GraphInfo(r"$N_{{b}}$", axis.Regular(4, 0, 4), lambda vg: vg.get_hist("N_bmedium")),
        # 'nbloose': GraphInfo(r"$N_{{b}}$", axis.Regular(5, 0, 5), lambda vg: vg.get_hist("N_bloose")),
        'jetpt':     GraphInfo('$p_{{T}}(j)$', axis.Regular(20, 0, 200), lambda vg : vg['Jets'].get_hist('pt', -1)),
        'njets':  GraphInfo('$N_{{j}}$', axis.Regular(6, 0, 6), lambda vg : (vg.Jets.num(), vg.scale)),
    }
    hist_factory = HistGetter(ntuple, year, workdir=workdir, inputdir=input_dir)

    for chan in chans:
        hist_factory.reset_mask()
        hist_factory.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
        plotter.set_extra_text(f'{chan}')
        print('\n', chan)

        hists = hist_factory.get_hists(graphs)
        data, np_mc, mc = 0, 0, 0
        for member, hist in hists['njets'].items():
            if member == "data":
                data += hist.integral()
            else:
                mc += hist.integral()
            print(member, hist.integral())

        print("diff:", abs(mc-data)/mc)

        plotter.plot_hists(hists)

def closure(workdir, year, input_dir, fake_qcd):
    chans = ['MM', 'EE', 'EM']
    ntuple = config.get_ntuple('fake_rate', 'closure_tt')
    plotter = Plotter(ntuple, year, data='nonprompt', bkg=['ttbar_lep', 'wjet_ht'], outdir=workdir / f'CR_{year}')
    plotter_mc = Plotter(ntuple, year, data='nonprompt_mc', bkg=['ttbar_lep', 'wjet_ht'], outdir=workdir / f'CR_{year}')
    plotter_mc.set_extra_text("MC")

    # Load MC scale factors
    with open(workdir/f"fr_{year}_new.pickle", "rb") as f:
        fake_rates = pickle.load(f)

    fake_type = 'data_ewk'

    graphs = {
        # 'pt1':    GraphInfo('$p_{{T}}(\ell_{{1}})$', axis.Regular(20, 0, 250), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        # 'pt2':    GraphInfo('$p_{{T}}(\ell_{{2}})$', axis.Regular(24, 0, 120), lambda vg : vg['AllLepton'].get_hist('pt', 1)),
        # 'eta1':   GraphInfo('$\eta(\ell_{{1}})$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),
        # 'eta2':   GraphInfo('$\eta(\ell_{{2}})$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 1)),
        'pt':     GraphInfo('$p_{{T}}(\ell)$', axis.Regular(20, 0, 300), lambda vg : vg['AllLepton'].get_hist('pt', -1)),
        'eta':    GraphInfo('$\eta(\ell)$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),
        'met':    GraphInfo('$MET$', axis.Regular(20, 0, 300), lambda vg : vg.get_hist("Met")),
        'ht':     GraphInfo('$H_T$', axis.Regular(20, 0, 750), lambda vg : vg.get_hist("HT")),
        # 'metphi': GraphInfo('$\phi(MET)$', axis.Regular(20, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
        # 'nbjets': GraphInfo(r"$N_{{b}}$", axis.Regular(4, 0, 4), lambda vg: vg.get_hist("N_bmedium")),
        # 'nbloose': GraphInfo(r"$N_{{b}}$", axis.Regular(5, 0, 5), lambda vg: vg.get_hist("N_bloose")),
        # 'jetpt':     GraphInfo('$p_{{T}}(j)$', axis.Regular(20, 0, 200), lambda vg : vg['Jets'].get_hist('pt', -1)),
        'njets':  GraphInfo('$N_{{j}}$', axis.Regular(6, 0, 6), lambda vg : (vg.Jets.num(), vg.scale)),
    }

    hist_factory = HistGetter(ntuple, year, workdir=workdir, inputdir=input_dir)

    for group, member, df in hist_factory.df_iter():
        if 'nonprompt' not in group:
            continue
        for chan in ["Muon", 'Electron']:
            scale_fake(df, chan, fake_rates[chan][fake_type])
            scale_fake(df, chan, fake_rates[chan][fake_type], i=1)
        flip_fake(df)

    for chan in chans:
        hist_factory.reset_mask()
        hist_factory.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
        plotter.set_extra_text(f'{fake_type}_{chan}')
        plotter_mc.set_extra_text(f'{fake_type}_{chan}_MC')
        print('\n', chan)

        hists = hist_factory.get_hists(graphs)

        data, np_mc, mc = 0, 0, 0
        for member, hist in hists['njets'].items():
            if member == "nonprompt":
                data += hist.integral()
            elif member == "nonprompt_mc":
                np_mc += hist.integral()
            else:
                mc += hist.integral()
            print(member, hist.integral())

        print("diff:", abs(mc-data)/mc)
        print("diff np:", abs(np_mc-mc)/mc)

        def get_chi2(hists):
            np_h = hists['nonprompt']
            np_mc_h = hists['nonprompt_mc']
            mc_h = hists['wjet_ht'] + hists['ttbar_lep']
            dof = len(np_h.vals)-1
            x2_data = chi2(np_h, mc_h)
            x2_mc = chi2(mc_h, np_mc_h)
            print('data vs mc: ', x2_data, 1-stats.chi2.cdf(x2_data, dof))
            print('mc vs mc: ', x2_mc, 1-stats.chi2.cdf(x2_mc, dof))

        for var in ['pt', 'eta']:
            print(var)
            print('--------------')
            get_chi2(hists[var])

        for name, hist in hists.items():
            syst = (0.3*hist['nonprompt'].values())**2
            syst_mc = (0.3*hist['nonprompt_mc'].values())**2

            plotter.plot_hist(name, hist, syst_err=syst)
            plotter_mc.plot_hist(name, hist, syst_err=syst)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-d', '--workdir', default=datetime.now().strftime("%m%d"),
                        help="directory to run over. If nothing, use date",)
    parser.add_argument('-i', '--input_dir', default=None)
    parser.add_argument('-r', '--run', type=lambda x: [i.strip() for i in x.split(',')],
                        help="Regions to run through (sideband, measurement, closure, dy)")
    parser.add_argument('--fakeqcd', action="store_true")
    args = parser.parse_args()

    workdir = workspace_area / args.workdir / 'fake_rate'
    workdir.mkdir(exist_ok=True, parents=True)

    for year in args.years:
        if 'sideband' in args.run:
            print("Sideband")
            sideband(workdir, year, args.input_dir, args.fakeqcd)

        if 'measurement' in args.run:
            print("Measurement")
            measurement(workdir,  year, args.input_dir, args.fakeqcd)
        if 'closure' in args.run:
            print("Closure")
            closure(workdir, year, args.input_dir, args.fakeqcd)
        if 'closure_tf' in args.run:
            print("Closure")
            closure_tf(workdir, year, args.input_dir, args.fakeqcd)
