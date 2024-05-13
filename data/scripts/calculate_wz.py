#!/usr/bin/env python3
import argparse
import numpy as np
import boost_histogram.axis as axis
import pickle
import itertools

import analysis_suite.commons.user as user
from analysis_suite.plotting.plotter import GraphInfo
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.commons.histogram import Histogram
import analysis_suite.commons.configs as config
from analysis_suite.commons.plot_utils import hep, plot, plot_colorbar
from analysis_suite.commons.constants import lumi

def scale_wz(vg, fakerate):
    njets = vg.Jets.num()
    jet_bin = np.digitize(njets, np.arange(1, 6)) - 1
    vg.scale = fakerate[jet_bin]

def measurement(workdir, year, input_dir):
    plot_dir = workdir / f'WZ_control_{year}'
    plot_dir.mkdir(exist_ok=True)

    njet_bins = axis.Regular(5, 1, 6)
    graphs = [
        GraphInfo('njets', '$N_{j}$', axis.Regular(7, 0, 7), lambda vg: (vg['Jets'].num(), vg.scale)),
        GraphInfo('njets', '$N_{j}$', njet_bins, lambda vg: (vg['Jets'].num(), vg.scale)),
        GraphInfo('nbjets', '$N_{b}$', axis.Regular(5, 0, 5), lambda vg: (vg['NBjets_medium'], vg.scale)),
        GraphInfo('nbjets_noscale', '$N_{b}$', axis.Regular(5, 0, 5), lambda vg: (vg['NBjets_medium'], vg.scale/vg["NBjets_wgt_medium"])),
        GraphInfo('ntightbjets', '$N_{b}$', axis.Regular(5, 0, 5), lambda vg: (vg['NBjets_tight'], vg.scale*vg["NBjets_wgt_tight"]/vg["NBjets_wgt_medium"])),
        GraphInfo('nloosebjets', '$N_{b}$', axis.Regular(5, 0, 5), lambda vg: (vg['NBjets_loose'], vg.scale*vg["NBjets_wgt_loose"]/vg["NBjets_wgt_medium"])),
        GraphInfo('zmass', '$M_{Z}$', axis.Regular(15, 75, 105), lambda vg: vg.get_hist("Zmass")),
        GraphInfo('met', '$p_{T}^{miss}$ (GeV)', axis.Regular(20, 0, 300), lambda vg: vg.get_hist("Met")),
        GraphInfo('l1Pt', '$p_{T}(l_{1})$ (GeV)', axis.Regular(20, 0, 500), lambda vg: vg['TightLepton'].get_hist('pt', 0)),
        GraphInfo('l2Pt', '$p_{T}(l_{2})$ (GeV)', axis.Regular(20, 0, 300), lambda vg: vg['TightLepton'].get_hist('pt', 0)),
        GraphInfo('l3Pt', '$p_{T}(l_{3})$ (GeV)', axis.Regular(20, 0, 200), lambda vg: vg['TightLepton'].get_hist('pt', 0)),
        GraphInfo('nmu', '$N_{\mu}$', axis.Regular(4, 0, 4), lambda vg: (vg["TightMuon"].num(), vg.scale)),

        GraphInfo('hlt_dilep', '$HLT_{2\ell}$', axis.Regular(2, 0, 2), lambda vg: (vg["HLT_dilepton"], vg.scale)),
        GraphInfo('hlt_trilep', '$HLT_{3\ell}$', axis.Regular(2, 0, 2), lambda vg: (vg["HLT_trilepton"], vg.scale)),
        GraphInfo('hlt_dilep_ht', '$HLT_{2\ell+HT}$', axis.Regular(2, 0, 2), lambda vg: (vg["HLT_dilepton_HT"], vg.scale)),

    ]

    ntuple = config.get_ntuple('WZ')
    ginfo = ntuple.get_info()
    mc = [i for i in ginfo.group2color.keys() if i != 'data']
    filename = ntuple.get_filename(year=year, workdir=input_dir)
    plotter = Plotter(filename, ginfo.setup_groups(), ntuple=ntuple, year=year)
    plotter.set_groups(bkg=mc)

    def trig_test(name, isMC=True):
        trig_by_mu = np.zeros((6, 4))
        for member, dfs in plotter.dfs.items():
            if (member == 'data') is isMC:
                continue
            for tree, df in dfs.items():
                scale = df.scale
                trig_info = df[name]
                n_trig = len(trig_info[0])
                for mu in range(4):
                    trig_mask = np.zeros(len(scale), dtype=bool)
                    for trig in range(n_trig):
                        trig_mask = np.any((trig_mask, trig_info[:,trig]), axis=0)
                        trig_by_mu[trig, mu] += np.sum(scale[(df["TightMuon"].num() == mu)*(trig_info[:,trig])])
                    trig_by_mu[n_trig, mu] += np.sum(scale[(df["TightMuon"].num() == mu)*(trig_mask)])
            # print(member)
            # print(trig_by_mu)
        print(trig_by_mu)

    def trig_pass(name, mapping, isMC=True):
        trig_by_mu = np.zeros((4, 2))
        for member, dfs in plotter.dfs.items():
            if (member == 'data') is isMC:
                continue
            for tree, df in dfs.items():
                scale = df.scale
                trig_info = df[name]
                n_trig = len(trig_info[0])
                for mu, trigs in mapping.items():
                    mu_mask = df["TightMuon"].num() == mu
                    trig_mask = np.zeros(len(scale), dtype=bool)
                    for trig in trigs:
                        trig_mask = np.any((trig_mask, trig_info[:,trig]), axis=0)
                    trig_by_mu[mu, 0] += np.sum(scale[mu_mask])
                    trig_by_mu[mu, 1] += np.sum(scale[mu_mask*trig_mask])
            # print(member)
            # print(trig_by_mu)
        print(trig_by_mu)


    dilept_map = {0: [0], 1: [0, 1, 2], 2: [1, 2, 3], 3: [3] }
    dilept_ht_map = {0: [0], 1: [0, 1], 2: [1, 2], 3: [2] }
    trilepton_map = {0: [0], 1: [1], 2: [2], 3: [3] }

    # trig_test('hlt_ind_dilepton')
    # trig_test('hlt_ind_dilepton', isMC=False)
    # print()
    # trig_test('hlt_ind_dilepton_HT')
    # trig_test('hlt_ind_dilepton_HT', isMC=False)
    # print()
    # trig_test('hlt_ind_trilepton')
    # trig_test('hlt_ind_trilepton', isMC=False)
    # print()
    # print("-"*40)
    # trig_pass('hlt_ind_dilepton', dilept_map)
    # trig_pass('hlt_ind_dilepton', dilept_map, isMC=False)
    # print()
    # trig_pass('hlt_ind_dilepton_HT', dilept_ht_map)
    # trig_pass('hlt_ind_dilepton_HT', dilept_ht_map, isMC=False)
    # print()
    # trig_pass('hlt_ind_trilepton', trilepton_map)
    # trig_pass('hlt_ind_trilepton', trilepton_map, isMC=False)
    # print()
    # exit()

    plotter.fill_hists(graphs)
    for graph in graphs:
        plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{year}.png')

    # for name, hist in plotter.get_hists("nbjets").items():
    #     print(name, hist.vals)


    # return
    # for name, hist in plotter.get_hists("nmu").items():
    #     print(name, hist.vals)

    data, mc, wz = np.zeros((3, 5))
    data_err, mc_err, wz_err = np.zeros((3, 5))

    for name, hist in plotter.get_hists("njets").items():
        if name == 'data':
            data += hist.vals
            data_err += hist.sumw2
        elif name == 'wz':
            wz += hist.vals
            wz_err += hist.sumw2
        else:
            mc += hist.vals
            mc_err += hist.sumw2

    scale = (data-mc)/wz
    scale_err = data_err/wz**2
    scale_hist = Histogram("", njet_bins)
    scale_hist.hist.values()[:] = scale
    scale_hist.hist.variances()[:] = scale_err
    print(scale)
    print(scale_err)
    print()

    with plot(plot_dir/f'wz_scale_{year}.png') as ax:
        scale_hist.plot_points(ax)
        for i, sc in enumerate(scale):
            print(i, sc)
            ax.annotate(f'${sc:0.3f}\pm{np.sqrt(scale_err[i]):0.3f}$', (i+1, sc),
                        xytext=(i+1.1, sc*1.1), fontsize=15)
        ax.set_xlabel(r'$N_{j}$')
        hep.cms.label(ax=ax, lumi=lumi[year])

    # Dump MC scale factors
    scale_file = workdir/f"wz_scale_factor.pickle"
    if scale_file.exists():
        with open(scale_file, "rb") as f:
            wz_scales = pickle.load(f)
    else:
        wz_scales = dict()
    wz_scales[year] = scale_hist

    with open(scale_file, "wb") as f:
        wz_scales = pickle.dump(wz_scales, f)

    plotter.scale(lambda vg: scale_wz(vg, scale), 'wz')
    plotter.fill_hists(graphs)
    for graph in graphs:
        plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_scaled_{year}.png')

def ttz_test(workdir, year, input_dir):
    plot_dir = workdir / f'ttz_control_{year}'
    plot_dir.mkdir(exist_ok=True)

    scale_file = workdir/f"wz_scale_factor.pickle"
    with open(scale_file, "rb") as f:
        wz_scales = pickle.load(f)

    graphs = [
        GraphInfo('njets', '$N_{j}$', axis.Regular(7, 0, 7), lambda vg: (vg['Jets'].num(), vg.scale)),
        # GraphInfo('njets', '$N_{j}$', njet_bins, lambda vg: (vg['Jets'].num(), vg.scale)),
        GraphInfo('nbjets', '$N_{b}$', axis.Regular(5, 0, 5), lambda vg: (vg['NBjets_medium'], vg.scale)),
        GraphInfo('ht', '$H_{T}$', axis.Regular(20, 250, 1250), lambda vg: (vg['HT'], vg.scale)),
        GraphInfo('ntightbjets', '$N_{b}$', axis.Regular(5, 0, 5), lambda vg: (vg['NBjets_tight'], vg.scale*vg["NBjets_wgt_tight"]/vg["NBjets_wgt_medium"])),
        GraphInfo('nloosebjets', '$N_{b}$', axis.Regular(5, 0, 5), lambda vg: (vg['NBjets_loose'], vg.scale*vg["NBjets_wgt_loose"]/vg["NBjets_wgt_medium"])),
        GraphInfo('zmass', '$M_{Z}$', axis.Regular(15, 75, 105), lambda vg: vg.get_hist("Zmass")),
        GraphInfo('met', '$p_{T}^{miss}$ (GeV)', axis.Regular(20, 0, 300), lambda vg: vg.get_hist("Met")),
        GraphInfo('l1Pt', '$p_{T}(l_{1})$ (GeV)', axis.Regular(20, 0, 500), lambda vg: vg['TightLepton'].get_hist('pt', 0)),
        GraphInfo('l2Pt', '$p_{T}(l_{2})$ (GeV)', axis.Regular(20, 0, 300), lambda vg: vg['TightLepton'].get_hist('pt', 0)),
        GraphInfo('l3Pt', '$p_{T}(l_{3})$ (GeV)', axis.Regular(20, 0, 200), lambda vg: vg['TightLepton'].get_hist('pt', 0)),
        GraphInfo('nmu', '$N_{\mu}$', axis.Regular(4, 0, 4), lambda vg: (vg["TightMuon"].num(), vg.scale)),
    ]

    ntuple = config.get_ntuple('ttzCR')
    ginfo = ntuple.get_info()
    mc = [i for i in ginfo.group2color.keys() if i != 'data']
    filename = ntuple.get_filename(year=year, workdir=input_dir)
    plotter = Plotter(filename, ginfo.setup_groups(), ntuple=ntuple, year=year)
    plotter.set_groups(bkg=mc)

    plotter.fill_hists(graphs)
    for graph in graphs:
        plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_unscaled_{year}.png')

    plotter.scale(lambda vg: scale_wz(vg, wz_scales[year].vals), 'wz')
    plotter.fill_hists(graphs)
    for graph in graphs:
        plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{year}.png')




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016pre", "2016post", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area / x,
                        help="Working Directory")
    parser.add_argument('-i', '--input_dir', default=None)
    parser.add_argument('-r', '--run', type=lambda x: [i.strip() for i in x.split(',')],
                        help="Regions to run through (sideband, measurement, closure, dy)")
    args = parser.parse_args()

    for year in args.years:
        if 'measurement' in args.run:
            measurement(args.workdir, year, args.input_dir)
        if 'ttz' in args.run:
            ttz_test(args.workdir, year, args.input_dir)
