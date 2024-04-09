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

from analysis_suite.plotting.plotter import GraphInfo
import analysis_suite.commons.configs as config
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import hep, plot, plot_colorbar
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.info import GroupInfo
import analysis_suite.data.plotInfo.nonprompt_fakerate as pinfo
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.commons.user import workspace_area
from scipy import stats

from analysis_suite.combine.hist_writer import HistWriter
from analysis_suite.combine.card_maker import Card_Maker

def get_fake_rate(part, fake_rate, idx, name="Muon"):
    if isinstance(fake_rate, np.ndarray):
        pt_axis = pinfo.np_ptbins[name]
        eta_axis = pinfo.np_etabins[name]
        npt, neta = pt_axis.size, eta_axis.size
    else:
        pt_axis, eta_axis = fake_rate.axes
        npt, neta = fake_rate.axes.size
        fake_rate = fake_rate.valutes()
    ptbin = np.digitize(part['pt', idx], pt_axis.edges) - 1
    ptbin = np.where(ptbin >= npt, npt-1, ptbin)
    etabin = np.digitize(part.abseta(idx), eta_axis.edges) - 1
    return fake_rate.flatten()[etabin + neta*ptbin]


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

def eff_err2(top, bot, bot_w2):
    alf = (1-0.682689492137)/2
    pas = top*bot/(bot_w2+1e-6)+1
    tot = (bot-top)*bot/(bot_w2+1e-6)+1

    lo,hi = beta.ppf(alf, pas, tot), beta.ppf(1-alf, pas, tot)
    eff = pas/tot
    return ((eff-lo)**2 + (hi-eff)**2)/2


def plot_project(filename, tight, loose, axis, axis_label, lumi):
    with plot(filename) as ax:
        eff_hist = Histogram.efficiency(tight.project(axis), loose.project(axis))
        eff_hist.plot_points(ax)
        ax.set_xlabel(axis_label)
        hep.cms.label(ax=ax, lumi=lumi)

def plot_fr_diff(filename, data_ewk, qcd, axis_label, lumi):
    data_ewk.name = "Data-EWK"
    qcd.name = "QCD"
    qcd.color = 'r'

    with plot(filename) as ax:
        data_ewk.plot_points(ax, capsize=1)
        qcd.plot_points(ax, capsize=1)
        ax.legend()
        ax.set_xlabel(axis_label)
        hep.cms.label(ax=ax, lumi=lumi)

def fit_template(data, qcd, ewk, fit_type="chi2"):
    def log_gamma(val):
        return val*np.log(val)-val+0.5*np.log(2*np.pi*val)

    def log_poisson(obs, mean):
        return (mean-obs)-obs*np.log(mean/obs)+0.5*np.log(2*np.pi*obs)

    def log_gaussian(f, err_ratio):
        return 0.5*((f-1)/err_ratio)**2+0.5*np.log(2*np.pi*err_ratio**2)

    def log_norm(f, val, err):
        return 0.5*((np.log(f*val)-val)/err)**2 + np.log(f*val*err)

    def likelihood(factors, data, qcd, ewk):
        mc = factors[0]*qcd.vals + factors[1]*ewk.vals
        base = np.sum(log_poisson(data.vals, mc))
        err_qcd = log_gaussian(factors[0], factors[2])
        err_ewk = log_gaussian(factors[1], factors[3])
        return base+err_qcd+err_ewk

    def chi2(factors, data, qcd, ewk):
        mc = qcd.hist*factors[0] + ewk.hist*factors[1]
        tot_diff2 = np.where(data.vals < 0.1, 0, (data.vals - mc.view().value)**2)
        return np.sum(tot_diff2/(mc.view().value+1e-6))

    start_val = np.sum(data.vals)/(np.sum(qcd.vals+ewk.vals))
    if fit_type == "chi2":
        res = minimize(chi2, (start_val, start_val), args=(data, qcd, ewk), method='Nelder-Mead',)
    elif fit_type == "ml":
        res = minimize(likelihood, (1, 1, 0.1, 0.1), args=(data, qcd, ewk), method='Nelder-Mead')
    elif fit_type == "signal":
        res = minimize(chi2_single, (start_val), args=(data, qcd, ewk), method='Nelder-Mead')
    return res.x

def fit_fakes(data, mc):
    def chi2(factors, data, mc):
        new_data = np.dot(factors, data)
        return np.sum((new_data - mc)**2/(mc+1e-6))

    total_data = np.sum(data)
    start_val = np.full(len(data), np.sum(mc)/total_data)
    res = minimize(chi2, start_val, args=(data, mc), method='Nelder-Mead')
    print(res.x)
    print(res.x/(1+res.x))

def scale_trigger(vg, chan, trig_scale, ptCut):
    pt = vg["AllLepton"]['rawPt', 0]
    num = vg[chan].num() == 1
    vg.scale = (trig_scale[0], (pt < ptCut)*num)
    vg.scale = (trig_scale[1], (pt >= ptCut)*num)

def scale_trigger_pt(vg, chan, trig_scale):
    pt_bins = pinfo.np_ptbins[chan].edges
    pt_bins[-1] = 1000
    pt = vg["AllLepton"]['rawPt', 0]
    num = vg[chan].num() == 1
    for i, pt_bot in enumerate(pt_bins[:-1]):
        vg.scale = (trig_scale[i], (pt >= pt_bot)*(pt < pt_bins[i+1])*num)

def scale_template(vg, chan, scales):
    pt_bins = pinfo.np_ptbins[chan].edges
    eta_bins = pinfo.np_etabins[chan].edges
    npt, neta = len(pt_bins)-1, len(eta_bins)-1

    ptbin = np.digitize(vg[chan]['pt', 0], pt_bins) - 1
    ptbin = np.where(ptbin >= npt, npt-1, ptbin)
    etabin = np.digitize(vg[chan]['abseta', 0], eta_bins) - 1
    vg.scale = (scales.flatten()[etabin + neta*ptbin], vg[chan].num() == 1)

def scale_fake(vg, chan, fakerate, i=0):
    part = vg[f'Fake{chan}']
    fr = get_fake_rate(part, fakerate, i, chan)
    vg.scale = (fr/(1-fr), part.num() > i)

def scale_bjet(vg, chan, wp):
    i = ['loose', 'medium', 'tight'].index(wp)
    vg.scale = vg['bjet_scale'][:, i][vg[chan].num() == 1], vg[chan].num() == 1

def scale_hlt(vg):
    pt = vg["AllLepton"]['pt', 0]
    hlt_scale = np.where(
        (pt < trig_cuts['Muon'])*(vg['Muon'].num() == 1)+(pt < trig_cuts['Electron'])*(vg['Electron'].num() == 1),
        vg['hlt_loPt_prescale'], vg['hlt_hiPt_prescale'])
    vg.scale = hlt_scale


def flip_fake(vg):
    vg.scale = (-1, (vg['FakeMuon'].num() == 2)+(vg['FakeElectron'].num() == 2))
    # vg.scale = (-1, (vg['FakeLepton'].num() == 2))

def fr_plot(name, ratio, chan, **kwargs):
    part = latex_chan[chan]
    with plot(name) as ax:
        mesh = ratio.plot_2d(ax, **kwargs)
        ax.set_xlabel(f"$p_{{T}}({part})$ [GeV]")
        ax.set_ylabel(f"$|\eta({part})|$")
        plot_colorbar(mesh, ax)

def fix_negative(data_ewk, qcd):
    axes = data_ewk.axes
    for x in range(axes.size[0]):
        for y in range(axes.size[1]):
            if data_ewk.hist[x, y].value <= 0 or np.isnan(data_ewk.hist[x, y].value):
                data_ewk.hist[x, y] = qcd.hist[x, y]


def plot_ptcorr(filename, pt, bins, lumi, chan, fact=None):
    with plot(filename) as ax:
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
        hep.cms.label(ax=ax, lumi=lumi, data=True)

def plot_ptcorr_2d(filename, vals, x, y):
    with plot(filename) as ax:
        # xx = np.tile(mva_bin, (len(pt_bin)-1, 1))
        # yy = np.tile(pt_bin, (len(mva_bin)-1, 1)).T
        # print(xx, yy, vals)
        mesh = ax.pcolormesh(x, y, vals, shading='flat')
        # mesh = ax.contourf(mva_bin[:-1], pt_bin[:-1], vals, 100)
        plot_colorbar(mesh, ax)

#--------------------------------------------------------------------------

def trigger_turnon(workdir, year, input_dir):
    plot_dir = workdir / f'trigger_eff'
    plot_dir.mkdir(exist_ok=True)

    mc = ['qcd', 'ewk']
    ntuple = config.get_ntuple('fake_rate', 'measurement')
    groups = ntuple.get_info().setup_groups(["data",] + mc)
    filename = ntuple.get_filename(year=year, workdir=input_dir)
    chans = ['Electron','Muon']

    graph = {
        "Electron": {
            "low": GraphInfo("rawpt", '$p_{{T}}({})$', axis.Regular(15, 10, 25), lambda vg : vg['AllLepton'].get_hist('rawPt', 0)),
            "high": GraphInfo("rawpt", '$p_{{T}}({})$', axis.Regular(15, 25, 40), lambda vg : vg['AllLepton'].get_hist('rawPt', 0)),
        },
        "Muon": {
            "low": GraphInfo("rawpt", '$p_{{T}}({})$', axis.Regular(15, 5, 20), lambda vg : vg['AllLepton'].get_hist('rawPt', 0)),
            "high": GraphInfo("rawpt", '$p_{{T}}({})$', axis.Regular(15, 20, 35), lambda vg : vg['AllLepton'].get_hist('rawPt', 0)),
        }
    }

    plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
    plotter.set_groups(bkg=mc)
    plotter.cut(lambda vg : vg["N_loose_mu"]+vg["N_loose_el"]==0)

    for chan in chans:
        latex = latex_chan[chan]
        for reg, op in {"low": operator.le, "high": operator.gt}.items():
            plotter.mask(lambda vg : vg[chan].num() == 1)
            plotter.mask(lambda vg : op(vg[chan]['rawPt', 0], trig_cuts[chan]), clear=False)

            plotter.fill_hists(graph[chan][reg], flow=False)
            plotter.plot_stack("rawpt", plot_dir/f'rawpt_{reg}_{chan}_{year}.png', chan=latex_chan[chan])

def percentages(workdir, year, input_dir):
    flav_split = GraphInfo("flav_split", "flavor", axis.Regular(23, 0, 23), lambda vg: vg['FakeLepton'].get_hist("jet_flav", 0))

    mc = ['qcd', 'ewk']
    ntuple = config.get_ntuple('fake_rate', 'measurement')
    groups = ntuple.get_info().setup_groups(mc)
    filename = ntuple.get_filename(year=year, workdir=input_dir)
    chans = ['Electron','Muon']

    measure_plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
    for chan in chans:
        for reg in ['ewk', 'qcd']:
            plotter.scale(lambda vg : scale_template(vg, chan, mc_scale_factors[chan][reg]), groups=reg)
    # cut_lowPt(measure_plotter)
    measure_plotter.set_groups(bkg=mc)

    # for name, df in {name: df for name, dfs in measure_plotter.dfs.items() for df in dfs.values()}.items():
    #     # print(name, sum(df.scale), sum(df['bjet_scale'][:,0])/len(df['bjet_scale'][:,0]), sum(df.scale*df['bjet_scale'][:,0]))
    #     df.scale = df['bjet_scale'][:,2]

    output = {chan : {bname: np.zeros(4) for bname in bjet_cuts.keys()} for chan in chans} # [all, b, c]
    for bname, cuts in bjet_cuts.items():
        if bname != "none": continue
        for chan in chans:
            out = output[chan][bname]
            measure_plotter.mask(lambda vg : vg[chan].num() == 1)
            measure_plotter.mask(lambda vg : vg["Jets"]["discriminator", 0] > cuts[year], clear=False)

            measure_plotter.fill_hists(flav_split)

            for name, hist in measure_plotter.get_hists("flav_split").items():
                out[0] += np.sum(hist.vals) # all
                out[1] += hist.vals[5] # bjets
                out[2] += hist.vals[4] # cjets
                out[3] += hist.vals[3] # light or unknown

            print(out)
            print(chan, bname, # out,
                  np.round(100*out[1:]/out[0], 1))
    return
    # TF setup
    ntuple = config.get_ntuple('fake_rate', 'closure_tf')
    filename = ntuple.get_filename(year=year, workdir=input_dir)
    groups = ntuple.get_info().setup_groups(["ttbar_lep", "wjet_ht",])
    plotter_tf = Plotter(filename, groups, ntuple=ntuple, year=year)
    cut_lowPt(plotter_tf)
    plotter_tf.cut(lambda vg : vg["N_loose_mu"]+vg["N_loose_el"]==0)

    output["Electron"]["closure"] = np.zeros(4)
    output["Muon"]["closure"] = np.zeros(4)
    for chan in chans:
        out = output[chan]['closure']
        plotter_tf.mask(lambda vg : vg[f'Fake{chan}'].num() == 1)
        plotter_tf.fill_hists(flav_split)
        for name, hist in plotter_tf.get_hists("flav_split").items():
            out[0] += np.sum(hist.vals) # all
            out[1] += hist.vals[5] # bjets
            out[2] += hist.vals[4] # cjets
            out[3] += hist.vals[3] # light or unknown

    for chan, bs in output.items():
        for bname, vals in bs.items():
            print(chan, bname, vals, vals/vals[0])



def sideband(workdir, year, input_dir, fake_qcd):
    plot_dir = workdir / f'SB_pteta_{year}'
    plot_dir.mkdir(exist_ok=True)

    chans = ['Electron', 'Muon']

    mc = ["qcd", "ewk"]
    ntuple = config.get_ntuple('fake_rate', 'sideband')
    groups = ntuple.get_info().setup_groups(["data"]+mc)
    filename = ntuple.get_filename(year=year, workdir=input_dir)

    plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
    plotter.set_groups(bkg=mc)

    # graphs = pinfo.nonprompt['SideBand_none']
    graphs = pinfo.nonprompt['SideBand']


    for chan in chans:
        plotter.mask_part(chan, 'pt',  lambda var : var > 20)

    plotter.cut(lambda vg : vg['TightLepton'].num() == 1, groups=['data', "ewk"])
    if fake_qcd:
        plotter.cut(lambda vg : vg['FakeLepton'].num() == 1, groups=['qcd'])
    else:
        plotter.cut(lambda vg : vg['TightLepton'].num() == 1, groups=['qcd'])

    fake_name = '_fake' if fake_qcd else ''
    mc_scale_factors = {chan: {reg: list() for reg in ['ewk', 'qcd',]} for chan in chans}

    round5 = lambda x : np.round(0.3*x+14)*5

    for chan in chans:
        pt_bins = pinfo.np_ptbins[chan].edges
        eta_bins = pinfo.np_etabins[chan].edges
        # pt_bins[-1] = 1000
        for i in range(len(pt_bins)-1):
            pt_lo, pt_hi = pt_bins[i], pt_bins[i+1]
            for eta_lo, eta_hi in zip(eta_bins[:-1], eta_bins[1:]):
                print(chan, f'{pt_lo}-{pt_hi}', f'{eta_lo}-{eta_hi}')
                plotter.mask(lambda vg : vg[chan].num() == 1)
                plotter.mask(lambda vg : vg[chan]['pt', 0] >= pt_lo, clear=False)
                plotter.mask(lambda vg : vg[chan]['pt', 0] < pt_hi, clear=False)
                plotter.mask(lambda vg : vg[chan]['abseta', 0] >= eta_lo, clear=False)
                plotter.mask(lambda vg : vg[chan]['abseta', 0] < eta_hi, clear=False)

                graphs[0].bin_tuple = axis.Regular(15, 0, round5(pt_lo))
                reg = "lo" if pt_hi < 30 else "hi"
                # plotter.fill_hists(graphs, 0, f"hlt_{reg}Pt_prescale")
                plotter.fill_hists(graphs)
                tot_val = {key: hist.integral() for key, hist in plotter.get_hists('mt').items()}

                # calculate templated fit
                tightmt = plotter.get_hists(f'mt')
                qcd_f, ewk_f = fit_template(tightmt['data'], tightmt['qcd'], tightmt["ewk"], fit_type="chi2")
                print(qcd_f, ewk_f)
                plotter.scale_hists('ewk', ewk_f)
                plotter.scale_hists('qcd', qcd_f)

                mc_scale_factors[chan]['ewk'].append(ewk_f)
                mc_scale_factors[chan]['qcd'].append(qcd_f)

                for graph in graphs:
                    plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{chan}_{int(pt_lo)}_{eta_lo}{fake_name}.png', chan=latex_chan[chan], region='$SB[{}]$')

        for name, val in mc_scale_factors[chan].items():
            print(name, val)
            mc_scale_factors[chan][name] = np.reshape(val, (len(pt_bins)-1, len(eta_bins)-1))

        scales = mc_scale_factors[chan]
        plot_ptcorr_2d(plot_dir/f'qcd_fact_{chan}{fake_name}.png', scales['qcd'], eta_bins, pt_bins)
        plot_ptcorr_2d(plot_dir/f'ewk_fact_{chan}{fake_name}.png', scales['ewk'], eta_bins, pt_bins)

    pickle.dump(mc_scale_factors, open(workdir/f"mc_scales{fake_name}_{year}.pickle", "wb"))



def mt_split(workdir, year, input_dir):
    plot_dir = workdir / f'mt_split_{year}'
    plot_dir.mkdir(exist_ok=True)

    mc = ['qcd', "ewk"]
    ntuple = config.get_ntuple('fake_rate', 'measurement')
    groups = ntuple.get_info().setup_groups( ["data",] + mc)
    filename = ntuple.get_filename(year=year, workdir=input_dir)
    chans = ['Electron','Muon']

    mt = [
        GraphInfo("mt", '$M_{{T}}({})$', axis.Regular(30, 0, 150), lambda vg : vg['AllLepton'].get_hist('mt', 0)),
        GraphInfo("mt_tight", '$M_{{T}}({})$', axis.Regular(30, 0, 150), lambda vg : vg['TightLepton'].get_hist('mt', 0)),
        GraphInfo("mt_fake", '$M_{{T}}({})$', axis.Regular(30, 0, 150), lambda vg : vg['FakeLepton'].get_hist('mt', 0)),
        GraphInfo("lepjet_mass", '$M_{{\ell,j}}$', axis.Regular(40, 0, 200), lambda vg : (vg.mass('AllLepton', 0, "Jets", 0), vg.scale)),
        GraphInfo("btag", "flavor", axis.Regular(80, 0, 0.4), lambda vg: vg['AllLepton'].get_hist("jet_btag", 0)),
        GraphInfo("btag_tight", "flavor", axis.Regular(80, 0, 0.4), lambda vg: vg['TightLepton'].get_hist("jet_btag", 0)),
        GraphInfo("btag_jet", "flavor", axis.Regular(80, 0, 1), lambda vg: vg['Jets'].get_hist('discriminator', 0)),
        GraphInfo("mvaTTH", "mvaTTH", axis.Regular(80, -1, 1), lambda vg: vg['FakeLepton'].get_hist('mvaTTH', 0)),
        GraphInfo("ptRatio", "ptRatio", axis.Regular(24, 0, 1.2), lambda vg: vg['FakeLepton'].get_hist("ptRatio", 0)),
        # GraphInfo("ptRatio2", "ptRatio2", axis.Regular(24, 0, 1.2), lambda vg: vg['FakeLepton'].get_hist("ptRatio2", 0)),
        GraphInfo("pt", '$M_{{T}}({})$', axis.Regular(40, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("pt_tight", '$M_{{T}}({})$', axis.Regular(40, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', 0)),
        GraphInfo("pt_fake", '$M_{{T}}({})$', axis.Regular(40, 0, 200), lambda vg : vg['FakeLepton'].get_hist('pt', 0)),
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(20, 0, 50), lambda vg : vg.get_hist("Met")),
        GraphInfo("ht", '$H_T$', axis.Regular(30, 0, 300), lambda vg : vg.get_hist("HT")),
        GraphInfo("dr", '$H_T$', axis.Regular(30, 0, 6), lambda vg : (vg.dr('AllLepton', 0, "Jets", 0), vg.scale)),
        GraphInfo("metphi", '$\phi(MET)$', axis.Regular(20, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
    ]

    # Load MC scale factors
    with open(workdir/f"mc_scales_{year}.pickle", "rb") as f:
        normal_factors = pickle.load(f)

    with open(workdir/f"mc_scales_fake_{year}.pickle", "rb") as f:
        fake_factors = pickle.load(f)

    plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
    plotter.set_groups(bkg=mc)

    for chan in chans:
        plotter.mask_part(chan, 'pt',  lambda var : var > 20)

    # plotter.mask_part('AllLepton', 'pt', lambda var : var > 20)
    # plotter.mask_part('Electron', 'mvaTTH', lambda var : var > -0.95)
    # plotter.mask_part('FakeLepton', 'ptRatio', lambda var : var > 0.6)
    plotter.cut(lambda vg: vg["AllLepton"].num() == 1)


    if year == '2016':
        bcuts = {"Muon": "loose", "Electron": "none"}
    if year == '2017':
        bcuts = {"Muon": "loose", "Electron": "loose"}
    elif year == '2018':
        bcuts = {"Muon": "loose", "Electron": "loose"}

    round5 = lambda x : np.round((x+7)/3)*5

    output_hist = {chan: {} for chan in chans}

    output = dict()
    for chan in chans:
        output[chan] = {i: list() for i in ['data_ewk', 'de_refit', 'de_fake', 'qcd']}

    for chan in chans:
        pt_bins = pinfo.np_ptbins[chan].edges
        eta_bins = pinfo.np_etabins[chan].edges

        pteta_loose = Histogram("", pinfo.np_ptbins[chan], pinfo.np_etabins[chan])
        pteta_tight = Histogram("", pinfo.np_ptbins[chan], pinfo.np_etabins[chan])
        pteta_qcd_loose = Histogram("", pinfo.np_ptbins[chan], pinfo.np_etabins[chan])
        pteta_qcd_tight = Histogram("", pinfo.np_ptbins[chan], pinfo.np_etabins[chan])
        latex = latex_chan[chan]
        for i in range(len(pt_bins)-1):
            pt_lo, pt_hi = pt_bins[i], pt_bins[i+1]
            for j in range(len(eta_bins)-1):
                eta_lo, eta_hi = eta_bins[j], eta_bins[j+1]

                print(chan, f'{pt_lo}-{pt_hi}', j, i)
                plotter.mask(lambda vg : vg[chan].num() == 1)
                plotter.mask(lambda vg : vg["Jets"]["discriminator", 0] > bjet_cuts[bcuts[chan]][year], clear=False)
                plotter.mask(lambda vg : vg[chan]['abseta', 0] >= eta_lo, clear=False)
                plotter.mask(lambda vg : vg[chan]['abseta', 0] < eta_hi, clear=False)
                plotter.mask(lambda vg : vg[chan]['pt', 0] >= pt_lo, clear=False)
                # if pt_hi != pt_bins[-1]:
                #     plotter.mask(lambda vg : vg[chan]['pt', 0] < pt_hi, clear=False)
                plotter.mask(lambda vg : vg[chan]['pt', 0] < pt_hi, clear=False)
                if pt_lo >= 50:
                    plotter.mask(lambda vg: vg[chan]['mt', 0] < 50, clear=False)

                mt[0].bin_tuple = axis.Regular(10, 0, round5(pt_lo))
                mt[1].bin_tuple = axis.Regular(10, 0, round5(pt_lo))
                mt[2].bin_tuple = axis.Regular(10, 0, round5(pt_lo))

                plotter.fill_hists(mt)

                loosemt = plotter.get_hists(f'mt')
                tightmt = plotter.get_hists(f'mt_tight')
                q_fact, e_fact = fit_template(tightmt['data'], tightmt['qcd'], tightmt["ewk"], fit_type="chi2")
                print(q_fact, e_fact)
                d_l, d_t = np.sum(loosemt['data'].vals), np.sum(tightmt['data'].vals)
                dl_err, dt_err = np.sum(loosemt['data'].sumw2), np.sum(tightmt['data'].sumw2)
                e_l, e_t = np.sum(loosemt['ewk'].vals), np.sum(tightmt['ewk'].vals)
                el_err, et_err = np.sum(loosemt['ewk'].sumw2), np.sum(tightmt['ewk'].sumw2)
                q_l, q_t = np.sum(loosemt['qcd'].vals), np.sum(tightmt['qcd'].vals)
                ql_err, qt_err = np.sum(loosemt['qcd'].sumw2), np.sum(tightmt['qcd'].sumw2)

                factor = normal_factors[chan]['ewk'][i,j]
                method1 = (d_t-factor*e_t)/(d_l-factor*e_l)
                pteta_loose.hist[i, j] = (d_l - factor*e_l, dl_err+factor**2*el_err)
                pteta_tight.hist[i, j] = (d_t - factor*e_t, dt_err+factor**2*et_err)
                pteta_qcd_loose.hist[i, j] = (q_l, ql_err)
                pteta_qcd_tight.hist[i, j] = (q_t, qt_err)

                factor = fake_factors[chan]['ewk'][i,j]
                method2 = (d_t-factor*e_t)/(d_l-factor*e_l)

                method3 = q_t/q_l

                factor = e_fact
                method4 = (d_t-factor*e_t)/(d_l-factor*e_l)

                output[chan]['data_ewk'].append(method1)
                output[chan]['de_fake'].append(method2)
                output[chan]['qcd'].append(method3)
                output[chan]['de_refit'].append(method4)

                print(f'{method1:0.4f}', f'{method2:0.4f}', f'{method3:0.4f}', f'{method4:0.4f}')
                plotter.scale_hists('qcd', normal_factors[chan]['qcd'][i,j])
                plotter.scale_hists('ewk', normal_factors[chan]['ewk'][i,j])

                for graph in mt:
                    plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{pt_lo}_{eta_lo}_{chan}.png', chan=latex, region=f'$MR[{latex}]$')

        fr_qcd = Histogram.efficiency(pteta_qcd_tight, pteta_qcd_loose)
        fr_ewk = Histogram.efficiency(pteta_tight, pteta_loose)

        fix_negative(fr_ewk, fr_qcd)

        output_hist[chan]['qcd'] = fr_qcd
        output_hist[chan]['ewk'] = fr_ewk

        fr_pt_qcd = Histogram.efficiency(pteta_qcd_tight.project(0), pteta_qcd_loose.project(0))
        fr_pt_ewk = Histogram.efficiency(pteta_tight.project(0), pteta_loose.project(0))

        plot_fr_diff(plot_dir/f'fr_{chan}_pt.png', fr_pt_ewk, fr_pt_qcd, f"$p_{{T}}({latex})$", lumi[year])

        plot_project(plot_dir/f'fr_qcd_{chan}_pt.png', pteta_qcd_tight, pteta_qcd_loose, 0, f"$p_{{T}}({latex})$", lumi[year])
        plot_project(plot_dir/f'fr_qcd_{chan}_eta.png', pteta_qcd_tight, pteta_qcd_loose, 1, f"$p_{{T}}({latex})$", lumi[year])
        plot_project(plot_dir/f'fr_ewk_{chan}_pt.png', pteta_tight, pteta_loose, 0, f"$p_{{T}}({latex})$", lumi[year])
        plot_project(plot_dir/f'fr_ewk_{chan}_eta.png', pteta_tight, pteta_loose, 1, f"$p_{{T}}({latex})$", lumi[year])

        fr_plot(plot_dir/f'fr_qcd_{chan}.png', fr_qcd, chan, vmin=0, vmax=0.4)
        fr_plot(plot_dir/f'fr_ewk_{chan}.png', fr_ewk, chan, vmin=0, vmax=0.4)
        print("fr", np.array2string(fr_ewk.vals, separator=','))


    with open(workdir/f"fr_test_{year}.pickle", "wb") as f:
        for chan in chans:
            for name, hist in output[chan].items():
                output[chan][name] = np.reshape(hist, (len(pt_bins)-1, len(eta_bins)-1))
        pickle.dump(output, f)

    with open(workdir/f"fr_{year}.pickle", "wb") as f:
        pickle.dump(output_hist, f)

def closure_test(workdir, year, input_dir):
    plot_dir = workdir / f'CR_test_{year}'
    plot_dir.mkdir(exist_ok=True)

    mc_list = ["ttbar_lep", "wjet_ht",]
    chans = ['MM', 'EE', 'EM']
    graphs = pinfo.nonprompt['Closure']

    # Load MC scale factors
    fake_rates = dict()
    with open(workdir/f"fr_test_{year}.pickle", "rb") as f:
        fake_rates = pickle.load(f)
        # fake_rates['Electron']['data_ewk'][-1, 0] = 0.14
        # fake_rates['Electron']['data_ewk'][-1, 1] = 0.14
        # fake_rates['Electron']['data_ewk'][-1, 2] = 0.14
        # fake_rates['Electron']['data_ewk'][-1] = fake_rates['Electron']['qcd'][-1]


    graphs = [
        GraphInfo("pt", '$p_{{T}}(\ell)$', axis.Regular(20, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("eta", '$\eta(\ell)$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),
        GraphInfo("ht", '$H_T$', axis.Regular(20, 0, 750), lambda vg : vg.get_hist("HT")),
    ]

    ntuple_tt = config.get_ntuple('fake_rate', 'closure_tt')
    groups = ntuple_tt.get_info().setup_groups(['nonprompt', 'nonprompt_mc'] + mc_list)
    filename = ntuple_tt.get_filename(year=year, workdir=input_dir)
    for fake_type in fake_rates['Electron'].keys():
        if fake_type != "data_ewk": continue
        print()
        print(fake_type)
        plotter_tt = Plotter(filename, groups, ntuple=ntuple_tt, year=year)

        for chan in ['Electron', "Muon"]:
            plotter_tt.mask_part(chan, 'pt',  lambda var : var > 20)
        plotter_tt.cut(lambda vg : vg["AllLepton"].num() == 2)

        for chan in ["Muon", 'Electron']:
            print(fake_rates[chan][fake_type])
            plotter_tt.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan][fake_type]), groups=['nonprompt', 'nonprompt_mc'])
            plotter_tt.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan][fake_type], i=1), groups=['nonprompt', 'nonprompt_mc'])
        plotter_tt.scale(lambda vg : flip_fake(vg), groups=['nonprompt', 'nonprompt_mc'])

        for chan in chans:
            print(chan)
            plotter_tt.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
            data, np_mc, mc = 0, 0, 0
            for vg, mem, group in plotter_tt.getters(keys=True):
                if group == "nonprompt":
                    data += sum(vg.scale)
                elif group == "nonprompt_mc":
                    np_mc += sum(vg.scale)
                else:
                    mc += sum(vg.scale)
            print('data', data)
            print('np_mc', np_mc)
            print("mc", mc)
            print("diff:", abs(mc-data)/mc)
            print("diff np:", abs(np_mc-mc)/mc)
            def get_chi2(plotter, var, np_name):
                hists = plotter.get_hists(var)
                np_h = hists[np_name].vals
                mc_h = hists['wjet_ht'].vals + hists['ttbar_lep'].vals
                mc_err = hists['wjet_ht'].sumw2 + hists['ttbar_lep'].sumw2
                pchi2 = np.sum((mc_h-np_h)**2/(mc_h+1e-5))
                print(var, pchi2/(len(np_h)-1), stats.chi2.sf(pchi2, len(np_h)-1))

            plotter_tt.set_groups(bkg=mc_list, data='nonprompt')
            plotter_tt.fill_hists(graphs)
            # for graph in graphs:
            #     plotter_tt.plot_stack(graph.name, plot_dir/f'{graph.name}_{chan}_{fake_type}.png', chan=latex_chan[chan], region='$CR[{}]$')
            get_chi2(plotter_tt, 'pt', 'nonprompt')
            get_chi2(plotter_tt, 'eta', 'nonprompt')
            get_chi2(plotter_tt, 'ht', 'nonprompt')

            plotter_tt.set_groups(bkg=mc_list, data='nonprompt_mc')
            plotter_tt.fill_hists(graphs)
            # for graph in graphs:
            #     plotter_tt.plot_stack(graph.name, plot_dir/f'{graph.name}_fromMC_{chan}_{fake_type}.png', chan=latex_chan[chan], region='$CR[{}]$')

            get_chi2(plotter_tt, 'pt', 'nonprompt_mc')
            get_chi2(plotter_tt, 'eta', 'nonprompt_mc')
            get_chi2(plotter_tt, 'ht', 'nonprompt_mc')


def closure_dy(workdir, year, input_dir):
    plot_dir = workdir / f'CR_DY_{year}'
    plot_dir.mkdir(exist_ok=True)

    chans = ['Muon', 'Electron',]
    graphs = pinfo.nonprompt['DY_closure']

    # Load MC scale factors
    fake_rates = dict()
    with open(workdir/f"fr_test_{year}.pickle", "rb") as f:
        fake_rates = pickle.load(f)

    ntuple = config.get_ntuple('fake_rate', 'dy_tight')
    groups = ntuple.get_info().setup_groups(['data', "VV", 'nonprompt', 'DY_J', 'nonprompt_mc'])
    filename = ntuple.get_filename(year=year, workdir=input_dir)
    for fake_type in fake_rates['Electron'].keys():
        if fake_type != "data_ewk": continue
        print()
        print(fake_type)
        plotter = Plotter(filename, groups, ntuple=ntuple, year=year)

        for chan in ['Electron', "Muon"]:
            plotter.mask_part(chan, 'pt',  lambda var : var > 20)
        plotter.cut(lambda vg : vg["AllLepton"].num() == 3)
        plotter.cut(lambda vg : vg["Met"] < 75)

        for chan in ["Muon", 'Electron']:
            print(fake_rates[chan][fake_type])
            plotter.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan][fake_type]), groups='nonprompt')
            plotter.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan][fake_type]), groups='nonprompt_mc')

        for chan in chans:
            print(chan)
            plotter.mask(lambda vg : vg[chan].num() == 1)

            def get_chi2(plotter, var, np_name):
                hists = plotter.get_hists(var)
                np_h = hists[np_name].vals
                mc_h = hists['DY_J'].vals
                print((np_h-mc_h)**2/(mc_h+1e-5))
                pchi2 = np.sum((mc_h-np_h)**2/(mc_h+1e-5))
                print(var, pchi2/(len(np_h)-1), stats.chi2.sf(pchi2, len(np_h)-1))

            # nonprompt+VV vs data
            plotter.set_groups(bkg=['nonprompt', 'VV'], data='data')
            plotter.fill_hists(graphs, chan)
            # for graph in graphs:
            #     plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{chan}_{fake_type}_data.png', chan=latex_chan[chan])

            # MC vs MC
            plotter.set_groups(bkg=['DY_J'], data='nonprompt_mc')
            plotter.fill_hists(graphs, chan)
            get_chi2(plotter, 'mass_z', 'nonprompt_mc')
            get_chi2(plotter, 'met', 'nonprompt_mc')
            get_chi2(plotter, 'mt', 'nonprompt_mc')
            # for graph in graphs:
            #     plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{chan}_{fake_type}_mc.png', chan=latex_chan[chan])

def measurement(workdir, year, input_dir, fake_qcd):
    plot_dir = workdir / f'MR_{year}'
    plot_dir.mkdir(exist_ok=True)

    mc = ['qcd' , 'ewk']
    chans = ['Electron','Muon']

    ntuple = config.get_ntuple('fake_rate', 'measurement')
    groups = ntuple.get_info().setup_groups(["data",] + mc)
    filename = ntuple.get_filename(year=year, workdir=input_dir)

    graphs = pinfo.nonprompt['Measurement_bjet']

    # Load MC scale factors
    with open(workdir/f"mc_scales_{year}.pickle", "rb") as f:
        mc_scale_factors = pickle.load(f)

    if fake_qcd:
        with open(workdir/f"mc_scales_fake_{year}.pickle", "rb") as f:
            tmp = pickle.load(f)
            for chan in chans:
                mc_scale_factors[chan]['ewk'] = tmp[chan]['ewk']# *tmp[chan]['data_mc_scale']
    fake_name = '_fake' if fake_qcd else ""

    plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
    # plotter.cut(lambda vg : (vg["AllLepton"]["mt", 0] < 55))
    plotter.set_groups(bkg=mc)

    for chan in chans:
        for reg in mc:
            plotter.scale(lambda vg : scale_template(vg, chan, mc_scale_factors[chan][reg]), groups=reg)

    plotter.mask_part('AllLepton', 'pt', lambda var : var > 20)
    # plotter.mask_part('FakeLepton', 'ptRatio', lambda var : var > 0.6)
    plotter.cut(lambda vg: vg["AllLepton"].num() == 1)
    plotter.cut(lambda vg: (vg['AllLepton']['pt', 0] < 50) + (vg['AllLepton']['mt', 0] < 50))


    if year == '2016':
        bcuts = {"Muon": "loose", "Electron": "none"}
    if year == '2017':
        bcuts = {"Muon": "loose", "Electron": "loose"}
    elif year == '2018':
        bcuts = {"Muon": "loose", "Electron": "loose"}

    fake_rates = {chan: dict() for chan in chans}
    fr_eta_bins, fr_pt_bins = pinfo.np_etabins, pinfo.np_ptbins

    for chan in chans:
        print(chan)
        latex = latex_chan[chan]
        np_bins = (fr_pt_bins[chan], fr_eta_bins[chan])
        # np_bins = (pinfo.pt_binning, fr_eta_bins[chan])
        fr_graph = [
            GraphInfo("tightfr", '', np_bins, lambda vg, chan : vg['TightLepton'].get_hist2d('pt', 'abseta', 0)),
            GraphInfo("loosefr", '', np_bins, lambda vg, chan : vg['AllLepton'].get_hist2d('pt', 'abseta', 0)),
        ]

        plotter.mask(lambda vg : vg[chan].num() == 1)
        plotter.mask(lambda vg : vg["Jets"]["discriminator", 0] > bjet_cuts[bcuts[chan]][year], clear=False)

        # plotter.fill_hists(graphs)
        # for graph in graphs:
        #     plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{chan}{fake_name}.png', chan=latex_chan[chan], region='$MR[{}]$')


        plotter.fill_hists(fr_graph, chan)
        loosefr = plotter.get_hists('loosefr')
        tightfr = plotter.get_hists('tightfr')

        loosefr['data_ewk'] = loosefr['data'] - loosefr['ewk']
        tightfr['data_ewk'] = tightfr['data'] - tightfr['ewk']


        fake_rates[chan]['qcd'] = Histogram.efficiency(tightfr['qcd'], loosefr['qcd'])
        fake_rates[chan]['data_ewk'] = Histogram.efficiency(tightfr['data_ewk'], loosefr['data_ewk'])

        fr_pt_qcd = Histogram.efficiency(tightfr['qcd'].project(0), loosefr['qcd'].project(0))
        fr_pt_ewk = Histogram.efficiency(tightfr['data_ewk'].project(0), loosefr['data_ewk'].project(0))

        plot_fr_diff(plot_dir/f'fr_{chan}_pt{fake_name}.png', fr_pt_ewk, fr_pt_qcd, f"$p_{{T}}({latex})$", lumi[year])

        for key, fr_ in fake_rates[chan].items():
            plot_project(plot_dir/f'fr_{key}_{chan}_pt{fake_name}.png', tightfr[key], loosefr[key], 0, f"$p_{{T}}({latex})$", lumi[year])
            plot_project(plot_dir/f'fr_{key}_{chan}_eta{fake_name}.png', tightfr[key], loosefr[key], 1, f"$p_{{T}}({latex})$", lumi[year])
            fr_plot(plot_dir/f'fr_{key}_{chan}{fake_name}.png', fr_, chan, vmin=0, vmax=0.4)
            print(key, "fr", np.array2string(fr_.vals, separator=','))

    # with open(workdir/f"fr_{year}{fake_name}.pickle", "wb") as f:
    #     pickle.dump(fake_rates, f)


def closure(workdir, year, input_dir, fake_qcd):
    plot_dir = workdir / f'CR_{year}'
    plot_dir.mkdir(exist_ok=True)

    mc_list = ["ttbar_lep", "wjet_ht",]
    chans = ['MM', 'EE', 'EM']
    graphs = pinfo.nonprompt['Closure']

    fake_name = '_fake' if fake_qcd else ""

    # Load MC scale factors
    fake_rates = dict()
    with open(workdir/f"fr_{year}{fake_name}.pickle", "rb") as f:
        fake_rates = pickle.load(f)

    # TF Setup
    ntuple_tf = config.get_ntuple('fake_rate', 'closure_tf')
    groups = ntuple_tf.get_info().setup_groups(['data'])
    filename = ntuple_tf.get_filename(year=year, workdir=input_dir)

    plotter_tf = Plotter(filename, groups, ntuple=ntuple_tf, year=year)
    fake_graphs = pinfo.nonprompt['FakeClosure']
    # plotter_tf.set_groups(bkg=['data'])
    for chan in ['Electron']:
        plotter_tf.mask(lambda vg : vg[f'Fake{chan}'].num() == 1)

        plotter_tf.fill_hists(fake_graphs)
        for graph in fake_graphs:
            plotter_tf.plot_stack(graph.name, plot_dir/f'{graph.name}_TF_{chan}.png')

    # # for chan in ["Muon", 'Electron']:
    # #     print(chan)
    # #     plotter_tf.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan]['qcd'].hist), groups='data')


    # for chan in ['EE', 'EM', 'ME', 'MM']:
    #     plotter_tf.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
    #     fake = "Muon" if chan[1] == 'M' else 'Electron'
    #     plotter_tf.fill_hists(fake_graphs)

    #     for graph in fake_graphs:
    #         plotter_tf.plot_stack(graph.name, plot_dir/f'{graph.name}_TF_{chan}.png')

    # exit()

    # TT Setup
    ntuple_tt = config.get_ntuple('fake_rate', 'closure_tt')
    groups = ntuple_tt.get_info().setup_groups(['nonprompt'] + mc_list)
    filename = ntuple_tt.get_filename(year=year, workdir=input_dir)
    plotter_tt = Plotter(filename, groups, ntuple=ntuple_tt, year=year)
    plotter_tt.set_groups(bkg=mc_list, data='nonprompt')

    fake_type = 'qcd'
    for chan in ['Electron', "Muon"]:
        plotter_tt.mask_part(chan, 'pt',  lambda var : var > 20)
        # plotter_tt.mask_part(f'Fake{chan}', 'ptRatio', lambda var : var > 0.6)
    plotter_tt.cut(lambda vg : vg["AllLepton"].num() == 2)

    for chan in ["Muon", 'Electron']:
        print(chan)
        plotter_tt.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan][fake_type].hist), groups='nonprompt')
        plotter_tt.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan][fake_type].hist, i=1), groups='nonprompt')
    plotter_tt.scale(lambda vg : flip_fake(vg), groups='nonprompt')

    for chan in chans:
        print(chan)
        plotter_tt.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
        data, mc = 0, 0
        for member, dfs in plotter_tt.dfs.items():
            for tree, df in dfs.items():
                if member == "nonprompt":
                    data += sum(df.scale)
                    print(tree, sum(df.scale))
                else:
                    mc += sum(df.scale)
        print("mc", mc)
        print("diff:", abs(mc-data)/mc)
    exit()
    for chan in chans:
        print(chan)
        plotter_tt.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
        plotter_tt.fill_hists(graphs)

        def get_chi2(plotter, var):
            hists = plotter.get_hists(var)
            np_h = hists['nonprompt'].vals
            mc_h = hists['wjet_ht'].vals + hists['ttbar_lep'].vals
            chi2 = (mc_h-np_h)**2/(mc_h+1e-5)
            print(var, np.sum(chi2), len(np_h))

        get_chi2(plotter_tt, 'pt')
        get_chi2(plotter_tt, 'eta')
        print()
        for graph in graphs:
            plotter_tt.plot_stack(graph.name, plot_dir/f'{graph.name}_TT_{fake_type}_{chan}{fake_name}.png', chan=latex_chan[chan], region="$TT({})$")


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
        if 'trigger' in args.run:
            print("Trigger Turn On")
            trigger_turnon(workdir, year, args.input_dir)

        if 'percentage' in args.run:
            print("BJet Percentage")
            percentages(workdir, year, args.input_dir)

        if 'sideband' in args.run:
            print("Sideband")
            sideband(workdir, year, args.input_dir, args.fakeqcd)

        if 'measurement' in args.run:
            print("Measurement")
            measurement(workdir,  year, args.input_dir, args.fakeqcd)

        if 'mt_split' in args.run:
            print("MT Split")
            mt_split(workdir, year, args.input_dir)

        if 'closure' in args.run:
            print("Closure")
            closure(workdir, year, args.input_dir, args.fakeqcd)

        if 'closure_test' in args.run:
            print("Closure")
            closure_test(workdir, year, args.input_dir)

        if 'closure_dy' in args.run:
            print("Closure")
            closure_dy(workdir, year, args.input_dir)
