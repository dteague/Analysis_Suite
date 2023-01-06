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
loose_name = {"Muon": "muon", "Electron": "elec"}

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
        tot_diff2 = np.where(data.vals < 0.1, 0, (data.vals - mc.view().value)**2)
        return np.sum(tot_diff2/(data.err+1e-6)**2)

    def chi2_single(factor, data, qcd, ewk):
        mc = qcd.hist*factors[0] + ewk.hist
        tot_diff2 = np.where(data.vals < 0.1, 0, (data.vals - mc.view().value)**2)
        return np.sum(tot_diff2/(data.err+1e-6)**2)


    start_val = np.sum(data.vals)/(np.sum(qcd.vals+ewk.vals))
    print(start_val)
    if fit_type == "chi2":
        res = minimize(chi2, (start_val, start_val), args=(data, qcd, ewk), method='Nelder-Mead')
    elif fit_type == "ml":
        res = minimize(likelihood, (start_val, start_val), args=(data.vals, qcd.vals, ewk.vals), method='Nelder-Mead')
        print(res)
    elif fit_type == "signal":
        res = minimize(chi2_single, (start_val), args=(data, qcd, ewk), method='Nelder-Mead')
    return res.x


def scale_trigger(vg, chan, trig_scale, ptCut):
    pt = vg["AllLepton"]['rawPt', 0]
    num = vg[chan].num() == 1
    vg.scale = (trig_scale[0], (pt < ptCut)*num)
    vg.scale = (trig_scale[1], (pt >= ptCut)*num)

def scale_fake(vg, chan, fakerate, i=0):
    part = vg[f'Fake{chan}']
    fr = get_fake_rate(part, fakerate, i, chan)
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

def conecorrection(workdir, ginfo, year, input_dir):
    groups = ginfo.setup_groups(# ["qcd", "ewk",]
                                ['data']
    )
    ntuple = config.get_ntuple('fake_rate', 'sideband')
    chans = ['Electron', 'Muon']
    filename = ntuple.get_filename(year=year, workdir=input_dir)

    mva_bins = np.linspace(-1., 1, 31)

    nbins = len(mva_bins) -1
    mva_cut_bin = np.where((mva_bins > 0.4)==1)[0][0]

    graphs = pinfo.nonprompt['SideBand']

    plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
    plotter.set_groups(bkg=['qcd', 'ewk'])
    fact = {"Electron": 0.85, "Muon": 0.725}
    # fact = {"Electron": 1, "Muon": 1}

    for chan in chans:
        print(' '*8 + f'"{chan}": {{')
        jet_pt = np.zeros(nbins)
        raw_pt = np.zeros(nbins)
        nom_pt =  np.zeros(nbins)
        wgt =  np.zeros(nbins)

        plotter.mask(lambda vg : vg[chan].num() == 1)
        for member, dfs in plotter.dfs.items():
            for tree, df in dfs.items():
                pt = df[chan]['rawPt', 0]
                jetpt = df[chan]['rawPt', 0]/df[chan]['ptRatio', 0]*fact[chan]
                nompt = df[chan]['pt', 0]
                mva = df[chan]['mvaTTH', 0]
                weight = df.scale
                for i in range(nbins):
                    mask = (mva > mva_bins[i]) * (mva < mva_bins[i+1])
                    jet_mask = mask * (jetpt > 15)
                    jet_pt[i] += np.sum(jetpt[jet_mask]*weight[jet_mask])
                    raw_pt[i] += np.sum(pt[mask]*weight[mask])
                    # nom_pt[i] += np.sum(nompt[mask]*weight[mask])
                    wgt[i] += np.sum(weight[mask])

        print(" "*12 + f"'jet': np.array({np.array2string(jet_pt/wgt, separator=',')}),")
        print(" "*12 + f"'raw': np.array({np.array2string(raw_pt/wgt, separator=',')}),")
        print(" "*12 + f"'wgt': np.array({np.array2string(wgt, separator=',')}),")
        print(" "*8 + '}, ')

def sideband(workdir, ginfo, year, input_dir):
    import operator

    groups = ginfo.setup_groups(["data", "qcd", "ewk",])
    ntuple = config.get_ntuple('fake_rate', 'sideband')
    filename = ntuple.get_filename(year=year, workdir=input_dir)
    plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
    # plotter = Plotter(filename, groups, ntuple=ntuple, year="2017_e")
    plotter.set_groups(bkg=['qcd', 'ewk'])

    plotter.mask_part("AllLepton", "iso", lambda var: var < 0.2)

    chans = ['Electron', 'Muon']
    trig_cuts = {"Muon": 20, "Electron": 25}
    mc_scale_factors = {chan: {reg: list() for reg in ['ewk', 'qcd', 'qcd_fake']} for chan in chans}
    qcd_amount = {chan: {"low": None, 'high': None} for chan in chans}
    fake_qcd = False

    if fake_qcd:
        plot_dir = workdir / f'SB_fake_{year}'
    else:
        plot_dir = workdir / f'SB_tight_{year}'
    plot_dir.mkdir(exist_ok=True)

    # for chan in chans:
    #     for reg, op in {"low": operator.le, "high": operator.gt}.items():
    #         plotter.mask(lambda vg : vg[f'Tight{chan}'].num() == 1)
    #         plotter.mask(lambda vg : op(vg[chan]['rawPt', 0], trig_cuts[chan]), clear=False)
    #         plotter.mask(lambda vg : vg["LooseLepton"].num() == 0, clear=False)
    #         qcd = 0
    #         data = 0
    #         ewk = 0
    #         for member, dfs in plotter.dfs.items():
    #             for tree, df in dfs.items():
    #                 if 'qcd' in member:
    #                     qcd += sum(df.scale)
    #                 elif 'data' in member:
    #                     data += sum(df.scale)
    #                 else:
    #                     ewk += sum(df.scale)
    #         qcd_amount[chan][reg] = (data, qcd, ewk)
    #         print(chan, reg, data, qcd, ewk)
    # plotter.mask(lambda vg : True)

    plotter.cut(lambda vg : vg['TightLepton'].num() == 1, groups=['data', 'ewk',])
    if fake_qcd:
        plotter.cut(lambda vg : vg['FakeLepton'].num() == 1, groups=['qcd'])
    else:
        plotter.cut(lambda vg : vg['TightLepton'].num() == 1, groups=['qcd'])

    for chan in chans:
        for reg, op in {"low": operator.le, "high": operator.gt}.items():
            graphs = [graph for graph in pinfo.nonprompt['SideBand']
                      if ("mt" not in graph.name and "pt" not in graph.name) or reg in graph.name]

            plotter.mask(lambda vg : vg[chan].num() == 1)
            plotter.mask(lambda vg : op(vg[chan]['rawPt', 0], trig_cuts[chan]), clear=False)


            plotter.fill_hists(graphs, ginfo)
            for graph in graphs:
                plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{reg}_{chan}.png', chan=latex_chan[chan], region='$SB[{}]$')

            plotter.mask(lambda vg : vg["LooseLepton"].num() == 0, clear=False)

            plotter.fill_hists(graphs, ginfo)
            for graph in graphs:
                plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{reg}_{chan}_veto.png', chan=latex_chan[chan], region='$SB[{}]$')

            # calculate templated fit
            tightmt = plotter.get_hists(f'tight_{reg}_mt')

            qcd_f, ewk_f = fit_template(tightmt['data'], tightmt['qcd'], tightmt["ewk"])
            mc_scale_factors[chan]['ewk'].append(ewk_f)
            if fake_qcd:
                data, qcd, ewk = qcd_amount[chan][reg]
                print("qcd fake: ", (data-ewk_f*ewk)/qcd)
                mc_scale_factors[chan]['qcd'].append((data-ewk_f*ewk)/qcd)
                mc_scale_factors[chan]['qcd_fake'].append(qcd_f)
            else:
                mc_scale_factors[chan]['qcd'].append(qcd_f)
            print(chan, qcd_f, ewk_f)

            # Scale template fit stuff
            plotter.scale_hists('ewk', ewk_f)
            plotter.scale_hists('qcd', qcd_f)

            for graph in graphs:
                plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{reg}_{chan}_scaled.png', chan=latex_chan[chan], region='$SB[{}]$')

    # Dump MC scale factors
    if fake_qcd:
        scale_file = f"mc_scales_fake_{year}.pickle"
    else:
        scale_file = f"mc_scales_{year}.pickle"

    with open(workdir/scale_file, "wb") as f:
        pickle.dump(mc_scale_factors, f)

def mt_split(workdir, ginfo, year, input_dir):
    plot_dir = workdir / f'mt_split_{year}'
    plot_dir.mkdir(exist_ok=True)

    mc = ['qcd', 'ewk']
    groups = ginfo.setup_groups(["data",] + mc)
    ntuple = config.get_ntuple('fake_rate', 'measurement')
    filename = ntuple.get_filename(year=year, workdir=input_dir)
    chans = ['Electron','Muon']
    mt = [pinfo.mt]

    # Load MC scale factors
    with open(workdir/f"mc_scales_{year}.pickle", "rb") as f:
        mc_scale_factors = pickle.load(f)

    plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
    
    plotter.set_groups(bkg=mc)
    for reg in mc:
        plotter.scale(lambda vg : scale_trigger(vg, "Muon", mc_scale_factors["Muon"][reg], 20), groups=reg)
        plotter.scale(lambda vg : scale_trigger(vg, "Electron", mc_scale_factors["Electron"][reg], 25), groups=reg)
    trig_cuts = {"Muon": 20, "Electron": 25}
    fr_opt = {"cmap":'jet', "vmin": 0.0, "vmax":0.5}

    pt_bins = [15, 20, 25, 35, 45, 60, 1000]
    for chan in chans:
        latex = latex_chan[chan]
        for i in range(len(pt_bins)-1):
            print(chan, pt_bins[i])
            plotter.mask(lambda vg : vg[chan].num() == 1)
            plotter.mask(lambda vg : vg[chan]['pt', 0] > pt_bins[i], clear=False)
            plotter.mask(lambda vg : vg[chan]['pt', 0] < pt_bins[i+1], clear=False)
            plotter.fill_hists(mt, ginfo)

            loosemt = plotter.get_hists('mt')
            pt_range = f'{pt_bins[i]}-{pt_bins[i+1]}'
            plotter.plot_stack("mt", plot_dir/f'mt_tight_{pt_range}_{chan}.png', chan=latex, region=f'$MR[{latex}]$')

            plotter.mask(lambda vg : vg[f'Tight{chan}'].num() == 1, clear=False)
            plotter.fill_hists(mt, ginfo)
            tightmt = plotter.get_hists('mt')
            plotter.plot_stack("mt", plot_dir/f'mt_tight_{pt_range}_{chan}.png', chan=latex, region=f'$MR[{latex}]$')

            loosemt['data_ewk'] = loosemt['data'] - loosemt['ewk']
            tightmt['data_ewk'] = tightmt['data'] - tightmt['ewk']
            for key in tightmt.keys():
                if key == "ewk" or key == "data":
                    continue
                print(key, np.cumsum(tightmt[key].vals)/np.cumsum(loosemt[key].vals))


def measurement(workdir, ginfo, year, input_dir):
    plot_dir = workdir / f'MR_{year}'
    plot_dir.mkdir(exist_ok=True)

    mc = ['qcd', 'ewk']
    groups = ginfo.setup_groups(["data",] + mc)
    ntuple = config.get_ntuple('fake_rate', 'measurement')
    filename = ntuple.get_filename(year=year, workdir=input_dir)
    chans = ['Electron','Muon']
    graphs = pinfo.nonprompt['Measurement']
    graphs_1d = [graph for graph in graphs if graph.dim() == 1]
    fr_eta_bins = pinfo.np_etabins
    fr_pt_bins = pinfo.np_ptbins

    # Load MC scale factors
    with open(workdir/f"mc_scales_{year}.pickle", "rb") as f:
        mc_scale_factors = pickle.load(f)
        print(mc_scale_factors)

    plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
    plotter.mask_part("AllLepton", "iso", lambda var: var < 0.2)
    plotter.mask_part("LooseLepton", "iso", lambda var: var < 0.2)
    plotter.cut(lambda vg : vg["AllLepton"].num() >= 1)
    plotter.set_groups(bkg=['ewk', 'qcd'])
    for reg in mc:
        plotter.scale(lambda vg : scale_trigger(vg, "Muon", mc_scale_factors["Muon"][reg], 20), groups=reg)
        plotter.scale(lambda vg : scale_trigger(vg, "Electron", mc_scale_factors["Electron"][reg], 25), groups=reg)
    trig_cuts = {"Muon": 20, "Electron": 25}
    fr_opt = {"cmap":'jet', "vmin": 0.0, "vmax":0.5}

    fake_rates = {chan: dict() for chan in chans}
    for chan in chans:
        print(chan)
        latex = latex_chan[chan]
        for i, graph in enumerate(graphs):
            if graph.name == "fr":
                graphs[i].bin_tuple = (fr_pt_bins[chan], fr_eta_bins[chan])
                break

        plotter.mask(lambda vg : vg[chan].num() == 1)
        plotter.mask(lambda vg : vg["LooseLepton"].num() == 0, clear=False)
        plotter.mask(lambda vg : vg[chan]['pt', 0] > 15, clear=False)
        plotter.mask(lambda vg : (vg[chan]['pt', 0] < 35) + (vg[chan]['mt', 0] < 55), clear=False)

        plotter.fill_hists(graphs, ginfo)
        loosefr = plotter.get_hists('fr')
        loosefr['data_ewk'] = loosefr['data'] - loosefr['ewk']
        for graph in graphs_1d:
            plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_loose_{chan}.png', chan=latex, region=f'$MR[{latex}]$')

        # Tight
        plotter.mask(lambda vg : vg[f'Tight{chan}'].num() == 1, clear=False)
        plotter.fill_hists(graphs, ginfo)
        tightfr = plotter.get_hists('fr')
        tightfr['data_ewk'] = tightfr['data'] - tightfr['ewk']
        fix_negative(tightfr['data_ewk'], tightfr['qcd'])
        for graph in graphs_1d:
            plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_tight_{chan}.png', chan=latex, region=f'$MR[{latex}]$')

        for key in tightfr.keys():
            fake_rates[chan][key] = Histogram.efficiency(tightfr[key], loosefr[key])
            if key == "qcd" or key == "data_ewk":
                print(key, "fr", np.array2string(Histogram.efficiency(tightfr[key], loosefr[key]).vals, separator=','))
                fr_plot(plot_dir/f"fr_{key}_{chan}.png", fake_rates[chan][key], chan, **fr_opt)
                plot_project(plot_dir/f'fr_{key}_{chan}_pt.png', tightfr[key], loosefr[key], 0, f"$p_{{T}}({latex})$", lumi[year])
                plot_project(plot_dir/f'fr_{key}_{chan}_eta.png', tightfr[key], loosefr[key], 1, f'$\eta({latex})$', lumi[year])

    # Dump Fake rates
    if "Muon" in fake_rates and "Electron" in fake_rates:
        with open(workdir/f"fr_{year}.pickle", "wb") as f:
            pickle.dump(fake_rates, f)


def closure(workdir, ginfo, year, input_dir):
    plot_dir = workdir / f'CR_{year}'
    plot_dir.mkdir(exist_ok=True)

    mc_list = ["ttbar_lep", "wjet_ht",]
    groups = ginfo.setup_groups(["data", 'nonprompt'] + mc_list)
    chans = ['MM', 'EE', 'EM']
    graphs = pinfo.nonprompt['Closure']

    # # TF setup
    # ntuple_tf = config.get_ntuple('fake_rate', 'closure_tf')
    # filename = ntuple_tf.get_filename(year=year, workdir=input_dir)
    # plotter_tf = Plotter(filename, groups, ntuple=ntuple_tf, year=year)


    # # for chan in chans:
    # #     print(chan)
    # #     for member, dfs in plotter_tf.dfs.items():
    # #         for tree, df in dfs.items():
    # #             print(member, np.mean(df["LooseMuon"].num()))
    # #             df["LooseMuon"].mask_part("iso", lambda part : part < 0.1)
    # #             print(member, np.mean(df["LooseMuon"].num()))


    #     # plotter_tf.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
    #     # plotter_tf.mask(lambda vg : vg["LooseLepton"].num() == 2, clear=False)
    #     # data = 0
    #     # mc = 0
    #     # for member, dfs in plotter_tf.dfs.items():
    #     #     for tree, df in dfs.items():
    #     #         if member == "data":
    #     #             data += sum(df.scale)
    #     #             print(tree, sum(df.scale))

    #     #         else:
    #     #             mc += sum(df.scale)
    #     # print("mc", mc)
    #     # print("diff:", abs(mc-data)/mc)
    # exit()

    # plotter_tf.set_groups(bkg=mc_list)
    # for chan in chans:
    #     plotter_tf.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
    #     plotter_tf.mask(lambda vg : vg["LooseLepton"].num() == 2, clear=False)

    #     plotter_tf.fill_hists(graphs, ginfo)
    #     for graph in graphs:
    #         plotter_tf.plot_stack(graph.name, plot_dir/f'{graph.name}_TF_{chan}.png', chan=latex_chan[chan], region="$TF({})$")


    # Load MC scale factors
    with open(workdir/f"fr_{year}.pickle", "rb") as f:
        fake_rates = pickle.load(f)

    # TT Setup
    ntuple_tt = config.get_ntuple('fake_rate', 'closure_tt')
    filename = ntuple_tt.get_filename(year=year, workdir=input_dir)
    plotter_tt = Plotter(filename, groups, ntuple=ntuple_tt, year=year)
    # plotter_tt.mask_part("AllLepton", "iso", lambda var: var < 0.3)
    plotter_tt.set_groups(bkg=mc_list, data='nonprompt')

    plotter_tt.mask(lambda vg: True)
    for chan in ["Muon", 'Electron']:
        plotter_tt.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan]['data_ewk'].hist), groups='nonprompt')
        plotter_tt.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan]['data_ewk'].hist, i=1), groups='nonprompt')
    plotter_tt.scale(lambda vg : flip_fake(vg), groups='nonprompt')

    for chan in chans:
        print(chan)
        plotter_tt.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
        plotter_tt.mask(lambda vg : vg["LooseLepton"].num() == 2, clear=False)
        data = 0
        mc = 0
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
        plotter_tt.mask(lambda vg : vg["Muon"].num() == chan.count('M'))
        plotter_tt.mask(lambda vg : vg["LooseLepton"].num() == 2, clear=False)

        plotter_tt.fill_hists(graphs, ginfo)
        for graph in graphs:
            plotter_tt.plot_stack(graph.name, plot_dir/f'{graph.name}_TT_{chan}.png', chan=latex_chan[chan], region="$TT({})$")


def dy_closure(workdir, ginfo, year, input_dir):
    plot_dir = workdir / f'DY_CR_{year}'
    plot_dir.mkdir(exist_ok=True)

    bkgs = ["DY", "ttbar", "VV"]
    groups = ginfo.setup_groups(["data", 'nonprompt']+bkgs)
    chans = ['Electron', "Muon"]
    op_chan = {"Electron": "Muon", "Muon": "Electron"}
    graphs = pinfo.nonprompt['DY_closure']

    # TF setup
    ntuple_fake = config.get_ntuple('fake_rate', 'dy_fake')
    filename = ntuple_fake.get_filename(year=year, workdir=input_dir)
    plotter_fake = Plotter(filename, groups, ntuple=ntuple_fake, year=year)

    plotter_fake.set_groups(bkg=bkgs)
    for chan in chans:
        region = '${0}{0}-{1}'.format(latex_chan[op_chan[chan]], latex_chan[chan]) + "_{{fake}}$"
        plotter_fake.mask(lambda vg : vg[chan].num() == 1)
        plotter_fake.mask(lambda vg : vg["LooseLepton"].num() == 3, clear=False)
        plotter_fake.fill_hists(graphs, ginfo, chan)
        for graph in graphs:
            plotter_fake.plot_stack(graph.name, plot_dir/f'{graph.name}_fake_{chan}.png', chan=latex_chan[chan], region=region)


    # Load MC scale factors
    with open(workdir/f"fr_{year}.pickle", "rb") as f:
        fake_rates = pickle.load(f)

    # TT Setup
    ntuple_tight = config.get_ntuple('fake_rate', 'dy_tight')
    filename = ntuple_tight.get_filename(year=year, workdir=input_dir)
    plotter_tight = Plotter(filename, groups, ntuple=ntuple_tight, year=year)

    for chan in ["Muon", 'Electron']:
        plotter_tight.scale(lambda vg : scale_fake(vg, chan, fake_rates[chan]['data_ewk'].hist), groups='nonprompt')
    print(plotter_tight.get_integral())

    for chan in chans:
        plotter_tight.mask(lambda vg : vg[chan].num() == 1)
        plotter_fake.mask(lambda vg : vg["LooseLepton"].num() == 3, clear=False)
        plotter_tight.fill_hists(graphs, ginfo, chan)
        print(plotter_tight.get_integral())
        for graph in graphs:
            latex = latex_chan[chan] if "z" not in graph.name else latex_chan[op_chan[chan]]
            region = '${0}{0}-{1}$'.format(latex_chan[op_chan[chan]], latex_chan[chan])
            plotter_tight.set_groups(bkg=['DY'], data='nonprompt')
            plotter_tight.plot_stack(graph.name, plot_dir/f'{graph.name}_tight_dy_{chan}.png', chan=latex, region=region)

            plotter_tight.set_groups(bkg=['VV', 'nonprompt'], data='data')
            plotter_tight.plot_stack(graph.name, plot_dir/f'{graph.name}_tight_all_{chan}.png', chan=latex, region=region)

            plotter_tight.set_groups(bkg=['nonprompt'], data='data')
            plotter_tight.hists[graph.name]['data'] -= plotter_tight.hists[graph.name]['VV']
            plotter_tight.plot_stack(graph.name, plot_dir/f'{graph.name}_tight_data_{chan}.png', chan=latex, region=region)


if __name__ == "__main__":
    workdir = workspace_area/'fake_rate'

    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-d', '--workdir', help="directory to run over. If nothing, use date",)
    parser.add_argument('-i', '--input_dir', default=None)
    parser.add_argument('-r', '--run', type=lambda x: [i.strip() for i in x.split(',')],
                        help="Regions to run through (sideband, measurement, closure, dy)")
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
        'nonprompt': 'grey',
        "DY_ht": "goldenrod",
        "DY": "goldenrod",
        "VV": 'mediumorchid',

        "ttbar_2l2n": 'blue',
        "ttbar_semilep": 'mediumblue',
        "ttbar_hadronic": 'cornflowerblue',
    }
    ginfo = GroupInfo(color_by_group)
    for year in args.years:
        if 'sideband' in args.run:
            print("Side Band")
            sideband(workdir, ginfo, year, args.input_dir)
        if 'mt_split' in args.run:
            print("mt_split")
            mt_split(workdir, ginfo, year, args.input_dir)

        if 'measurement' in args.run:
            print("Measurement")
            measurement(workdir, ginfo, year, args.input_dir)
        if 'closure' in args.run:
            print("Closure")
            closure(workdir, ginfo, year, args.input_dir)
        if 'dy' in args.run:
            print("DY Closure")
            dy_closure(workdir, ginfo, year, args.input_dir)
        if 'cone' in args.run:
            print("Cone Correction")
            conecorrection(workdir, ginfo, year, args.input_dir)
