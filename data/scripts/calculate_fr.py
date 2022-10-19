#!/usr/bin/env python3
import argparse
from scipy.optimize import minimize
import numpy as np
import pickle
import awkward as ak

import analysis_suite.commons.configs as config
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import hep, plot, plot_colorbar
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.info import GroupInfo
import analysis_suite.data.plotInfo.nonprompt_fakerate as pinfo
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.commons.user import workspace_area
from datetime import datetime

def get_fake_rate(part, fake_rate, idx, name="Muon"):
    pt_axis, eta_axis = fake_rate.axes
    npt, neta = fake_rate.axes.size

    ptbin = np.digitize(part['pt', idx], pt_axis.edges) - 1
    ptbin = np.where(ptbin >= npt, npt-1, ptbin)
    etabin = np.digitize(part.abseta(idx), eta_axis.edges) - 1
    return fake_rate.values().flatten()[etabin + neta*ptbin]


latex_chan = {"Electron": "e", "Muon": "\mu",
              "EE": 'ee', "EM": 'e\mu', 'MM': '\mu\mu'}

trig_scale = {
    "2016pre": {
        # "Electron" : [0.8670, 1.0171], # Ele8 & Ele 17
        "Electron" : [0.8670, 1.0125], # Ele8 & Ele23
        # "Electron" : [0.9158, 1.0125], # Ele12 & Ele23
        "Muon" : [0.7023, 0.9768],
    },
    "2016post": {
        "Electron" : [0.8670, 1.0125],
        "Muon" : [0.7023, 0.9768],
    },
    "2017": {
        "Electron" : [0.9540, 1.1406],
        "Muon" : [0.7123, 1.0358],
    },
    "2018": {
        "Electron" : [0.9491, 1.1130],
        "Muon" : [0.9074, 1.1205],
    },
}

def plot_project(filename, tight, loose, axis, axis_label, lumi):
    with plot(filename) as ax:
        eff_hist = Histogram.efficiency(tight.project(axis), loose.project(axis))
        eff_hist.plot_points(ax)
        ax.set_xlabel(axis_label)
        hep.cms.label(ax=ax, lumi=lumi)


def fit_template(data, qcd, ewk, fit_type="chi2"):
    def log_gamma(val):
        return val*np.log(val)-val+0.5*np.log(2*np.pi/val)

    def likelihood(factors, data, qcd, ewk):
        mc = factors[0]*qcd + factors[1]*ewk
        return np.sum(data) + np.sum(log_gamma(mc) - mc*np.log(data))

    def chi2(factors, data, qcd, ewk):
        mc = qcd.hist*factors[0] + ewk.hist*factors[1]
        tot_diff2 = (data.vals - mc.view().value)**2
        return np.sum(tot_diff2/data.err)
    # return np.sum(tot_diff2/mc.err)

    start_val = np.sum(data.vals)/(np.sum(qcd.vals+ewk.vals))
    if fit_type == "chi2":
        res = minimize(chi2, (start_val, start_val), args=(data, qcd, ewk), method='Nelder-Mead')
    elif fit_type == "ml":
        res = minimize(likelihood, (start_val, start_val), args=(data.vals, qcd.vals, ewk.vals), method='Nelder-Mead')
        print(res)
    return res.x

def scale_trigger(vg, chan, trig_scale, ptCut):
    pt = vg["AllLepton"]['pt', 0, True]
    num = vg[chan].num() == 1

    if f'Tight{chan}/rawPt' in vg.branches:
        pt = vg["AllLepton"]['rawPt', 0]
    else:
        ptratio = vg["AllLepton"]['ptRatio', 0, True]
        ptrel = vg["AllLepton"]['ptRel', 0, True]
        iso = vg["AllLepton"]['iso', 0, True]
        pt = np.where((ptrel > 8.)*(iso > 0.1), pt/(0.9+iso), pt)
        pt = np.where((ptrel < 8.)*(ptratio < 0.75), pt*ptratio/0.75, pt)

    print("before", chan, sum(vg.scale))
    vg.scale = (trig_scale[0], (pt < ptCut)*num)
    vg.scale = (trig_scale[1], (pt >= ptCut)*num)
    print("after", chan, sum(vg.scale))

def scale_fake(vg, chan, fakerate, i=0):
    part = vg[f'Fake{chan}']

    # fake_lead = vg[f'AllLepton']['pt', 0] == part['pt', 0]
    # lead_ratio = ak.count_nonzero(fake_lead)/len(fake_lead)
    # print(len(fake_lead), lead_ratio)
    # print(np.histogram(part['rawPt', 0][fake_lead], np.linspace(0, 200, 21))[0])
    # print(np.histogram(vg['AllLepton']['rawPt', 1][fake_lead], np.linspace(0, 200, 21))[0])

    # print(np.histogram(vg[f'AllLepton']['pt', -1], np.linspace(0, 200, 21))[0])
    # print(np.histogram(part['pt', -1], np.linspace(0, 200, 21))[0])
    fr = get_fake_rate(part, fakerate, i, chan)
    # fr = 0.05 if chan == "Muon" else 0.2
    vg.scale = (fr/(1-fr), part.num() > i)

def flip_fake(vg):
    vg.scale = (-1, vg['FakeLepton'].num() == 2)

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
            if data_ewk.hist[x, y].value < 0:
                data_ewk.hist[x, y] = qcd.hist[x, y]


def sideband(workdir, ginfo, year, input_dir):
    plot_dir = workdir / f'SB_{year}'
    plot_dir.mkdir(exist_ok=True)

    groups = ginfo.setup_groups(["data", "qcd", "ewk",])
    ntuple = config.get_ntuple('fake_rate', 'sideband')
    chans = ['Electron',
             'Muon']

    graphs = pinfo.nonprompt['SideBand']

    plotter = Plotter(ntuple.get_file(year=year, workdir=input_dir), groups, ntuple=ntuple, year=year)
    # exit()
    plotter.set_groups(bkg=['qcd', 'ewk'])
    plotter.scale(lambda vg : scale_trigger(vg, "Muon", trig_scale[year]["Muon"], 20), groups='data')
    plotter.scale(lambda vg : scale_trigger(vg, "Electron", trig_scale[year]["Electron"], 25), groups='data')

    mc_scale_factors = {chan: dict() for chan in chans}
    for chan in chans:
        plotter.mask(lambda vg : vg[f'Tight{chan}'].num() == 1)
        plotter.mask(lambda vg : vg[f'Tight{chan}']['pt', 0, True] > 20, clear=False)

        plotter.fill_hists(graphs, ginfo)
        scale_factor = mc_scale_factors[chan]

        for graph in graphs:
            plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{chan}.png', chan=latex_chan[chan])

        # calculate templated fit
        tightmt = plotter.get_hists('tightmt')
        qcd_f, ewk_f = fit_template(tightmt['data'], tightmt['qcd'], tightmt["ewk"])
        scale_factor['ewk'] = ewk_f
        scale_factor['qcd'] = qcd_f
        print(chan, qcd_f, ewk_f)

        # Scale template fit stuff
        plotter.scale_hists('ewk', ewk_f)
        plotter.scale_hists('qcd', qcd_f)

        for graph in graphs:
            plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_scaled_ratio_{chan}.png', chan=latex_chan[chan])

    # Dump MC scale factors
    with open(workdir/f"mc_scales_{year}.pickle", "wb") as f:
        pickle.dump(mc_scale_factors, f)


def measurement(workdir, ginfo, year, input_dir):
    plot_dir = workdir / f'MR_{year}'
    plot_dir.mkdir(exist_ok=True)

    groups = ginfo.setup_groups(["data", "qcd", "ewk",])
    ntuple = config.get_ntuple('fake_rate', 'measurement')
    chans = ['Electron',
             'Muon']

    graphs = pinfo.nonprompt['Measurement']
    graphs_1d = [graph for graph in graphs if graph.dim() == 1]
    fr_eta_bins = pinfo.np_etabins
    fr_pt_bins = pinfo.np_ptbins

    # Load MC scale factors
    with open(workdir/f"mc_scales_{year}.pickle", "rb") as f:
        mc_scale_factors = pickle.load(f)

    # mc_scale_factors = {
    #     "Muon": {'ewk': 0.9, "qcd": 0.7},
    #     "Electron": {'ewk': 0.9, "qcd": 0.7},
    # }

    plotter = Plotter(ntuple.get_file(year=year, workdir=input_dir), groups, ntuple=ntuple, year=year)
    plotter.set_groups(bkg=['ewk', 'qcd'])
    plotter.scale(lambda vg : scale_trigger(vg, "Muon", trig_scale[year]["Muon"], 20), groups='data')
    plotter.scale(lambda vg : scale_trigger(vg, "Electron", trig_scale[year]["Electron"], 25), groups='data')
    plotter.mask(lambda vg: vg.Muon.num() == 1)
    pt_cone = np.array([])
    scale = np.array([])
    mva = np.array([])
    for group, df in plotter.dfs.items():
        # if group != "data":
        #     continue
        for tree_name, tree in df.items():
            # scale = np.concatenate([scale, tree.scale])
            # mva = np.concatenate([mva, tree.Muon["mvaTTH", 0]])
            # pt_cone = np.concatenate([pt_cone, tree.Muon["rawPt", 0]/tree.Muon["ptRatio", 0]])
            tree.arr["FakeMuon/pt"] = tree.arr["FakeMuon/pt"]*1.2
    # nbins = 30
    # mva_bin = np.linspace(-1, 1, nbins+1)
    # pt_binned = list()
    # for i in range(nbins):
    #     mask = (mva >= mva_bin[i]) * (mva < mva_bin[i+1])
    #     pt_binned.append(np.sum(pt_cone[mask]*scale[mask])/np.sum(scale[mask]))
    # print(pt_binned)
    # exit()

    fr_opt = {"cmap":'jet', "vmin": 0.0, "vmax":0.5}

    fake_rates = {chan: dict() for chan in chans}
    for chan in chans:
        latex = latex_chan[chan]
        for i, graph in enumerate(graphs):
            if graph.name == "fr":
                graphs[i].bin_tuple = (fr_pt_bins[chan], fr_eta_bins[chan])
                break

        plotter.mask(lambda vg : vg[chan].num() == 1)
        plotter.mask(lambda vg : vg[chan]['pt', 0] > 15, clear=False)
        # plotter.mask(lambda vg : vg[chan]['mvaTTH', 0] < -0.9, clear=False)
        plotter.fill_hists(graphs, ginfo)

        # scale = 0.9
        # Scale template fit stuff
        for group, sf in mc_scale_factors[chan].items():
            plotter.scale_hists(group, sf)
        # plotter.scale_hists('ewk', scale)
        loosefr = plotter.get_hists('fr')
        for graph in graphs_1d:
            plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_loose_{chan}.png', chan=latex, region=f'$MR[{latex}]$')

        # Tight
        plotter.mask(lambda vg : vg[f'Tight{chan}'].num() == 1, clear=False)
        plotter.fill_hists(graphs, ginfo)
        for group, sf in mc_scale_factors[chan].items():
            plotter.scale_hists(group, sf)
        # plotter.scale_hists('ewk', scale)
        tightfr = plotter.get_hists('fr')
        for graph in graphs_1d:
            plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_tight_{chan}.png', chan=latex, region=f'$MR[{latex}]$')

        # Fake Rate
        tightfr['data_ewk'] = tightfr['data'] - tightfr['ewk']
        loosefr['data_ewk'] = loosefr['data'] - loosefr['ewk']
        fix_negative(tightfr['data_ewk'], tightfr['qcd'])
        for key in tightfr.keys():
            fake_rates[chan][key] = Histogram.efficiency(tightfr[key], loosefr[key])
            print(key, fake_rates[chan][key].vals)
            fr_plot(plot_dir/f"fr_{key}_{chan}.png", fake_rates[chan][key], chan, **fr_opt)

        plot_project(plot_dir/f'fr_pt_{chan}.png', tightfr['data_ewk'], loosefr['data_ewk'], 0, f"$p_{{T}}({latex})$", lumi[year])
        plot_project(plot_dir/f'fr_eta_{chan}.png', tightfr['data_ewk'], loosefr['data_ewk'], 1, f'$\eta({latex})$', lumi[year])

    # Dump Fake rates
    if fake_rates["Muon"] and fake_rates["Electron"]:
        with open(workdir/f"fr_{year}.pickle", "wb") as f:
            pickle.dump(fake_rates, f)


def closure(workdir, ginfo, year, input_dir):
    plot_dir = workdir / f'CR_{year}'
    plot_dir.mkdir(exist_ok=True)

    groups = ginfo.setup_groups(["data", "ttbar_lep", "wjet_ht",
                                 'nonprompt'])
    chans = ['MM', 'EE', 'EM']
    graphs = pinfo.nonprompt['Closure']

    # # TF setup
    # ntuple_tf = config.get_ntuple('fake_rate', 'closure_tf')
    # plotter_tf = Plotter(ntuple_tf.get_file(year=year, workdir=input_dir), groups, ntuple=ntuple_tf, year=year)
    # plotter_tf.cut(lambda vg : (vg["AllLepton"]['rawPt', 0] > 25) + (vg["AllLepton"]['rawPt', 1] > 25))
    # plotter_tf.cut(lambda vg : (vg["AllLepton"]['rawPt', 0] > 12) * (vg["AllLepton"]['rawPt', 1] > 12))

    # plotter_tf.set_groups(bkg=["ttbar_lep", "wjet_ht"])
    # for chan in chans:
    #     plotter_tf.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
    #     print(chan)
    #     mc = 0
    #     for member, dfs in plotter_tf.dfs.items():
    #         for tree, df in dfs.items():
    #             if len(df.scale) < 10:
    #                 continue
    #             if member == "data":
    #                 print(tree, sum(df.scale))
    #                 print("data", min(df["FakeLepton"]['rawPt', 0]))
    #             else:
    #                 mc += sum(df.scale)
    #                 print(member, min(df["FakeLepton"]['rawPt', 0]))
    #     print("mc", mc)

    #     plotter_tf.fill_hists(graphs, ginfo)
    #     for graph in graphs:
    #         plotter_tf.plot_stack(graph.name, plot_dir/f'{graph.name}_TF_{chan}.png', chan=latex_chan[chan], region="$TF({})$")

    # exit()
    # Load MC scale factors
    with open(workdir/f"fr_{year}.pickle", "rb") as f:
        fake_rates = pickle.load(f)

    # TT Setup
    ntuple_tt = config.get_ntuple('fake_rate', 'closure_tt')
    plotter_tt = Plotter(ntuple_tt.get_file(year=year, workdir=input_dir), groups, ntuple=ntuple_tt, year=year)
    plotter_tt.set_groups(bkg=["ttbar_lep", "wjet_ht"], data='nonprompt')

    for chan in chans:
        print(chan)
        plotter_tt.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
        data = 0
        mc = 0
        for member, dfs in plotter_tt.dfs.items():
            for tree, df in dfs.items():
                if member == "nonprompt":
                    print(tree, sum(df.scale))
                else:
                    mc += sum(df.scale)
        print("mc", mc)

    #     print(chan)
    #     print(plotter_tt.get_integral())

    # plotter_tt.mask(lambda vg : vg["Muon"].num() == 2)
            # print(member, len(df.scale), sum(df.scale), np.sqrt(sum(df.scale*df.scale)))
    plotter_tt.mask(lambda vg: True)
    for chan in ["Muon", 'Electron']:
        plotter_tt.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan]['qcd'].hist), groups='nonprompt')
        plotter_tt.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan]['qcd'].hist, i=1), groups='nonprompt')
    plotter_tt.scale(lambda vg : flip_fake(vg), groups='nonprompt')


    for chan in chans:
        plotter_tt.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
        plotter_tt.fill_hists(graphs, ginfo)
        for graph in graphs:
            plotter_tt.plot_stack(graph.name, plot_dir/f'{graph.name}_TT_{chan}.png', chan=latex_chan[chan], region="$TT({})$")


if __name__ == "__main__":
    workdir = workspace_area/'fake_rate'

    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-d', '--workdir', help="directory to run over. If nothing, use date",
                        # choices=[d.stem for d in workdir.glob('*/') if d.is_dir()]
                        )
    parser.add_argument('-i', '--input_dir', required=True)
    parser.add_argument('-r', '--run', type=lambda x: [i.strip() for i in x.split(',')],
                        help="Regions to run through (sideband, measurement, closure)")
    args = parser.parse_args()

    workdir /= datetime.now().strftime("%m%d") if args.workdir is None else args.workdir
    workdir.mkdir(exist_ok=True)

    color_by_group = {
        "data": "black",
        "qcd": "grey",
        "qcd_em": "grey",
        "qcd_mu": "grey",
        "ewk": "orange",
        "wjet_ht": "olive",
        "wjets": "olive",
        "ttbar_lep": "royalblue",
        'nonprompt': 'black'
    }
    ginfo = GroupInfo(color_by_group)
    for year in args.years:
        if 'sideband' in args.run:
            print("Side Band")
            sideband(workdir, ginfo, year, args.input_dir)
        if 'measurement' in args.run:
            print("Measurement")
            measurement(workdir, ginfo, year, args.input_dir)
        if 'closure' in args.run:
            print("Closure")
            closure(workdir, ginfo, year, args.input_dir)
