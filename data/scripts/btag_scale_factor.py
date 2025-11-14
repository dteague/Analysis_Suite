#!/usr/bin/env python3
import argparse
import boost_histogram.axis as axis
import uproot
import numpy as np
import collections.abc
import multiprocessing as mp
from copy import deepcopy
import warnings
import pickle

warnings.simplefilter("ignore", UserWarning)

from analysis_suite.commons.user import workspace_area
from analysis_suite.commons.constants import all_eras
from analysis_suite.commons.configs import get_ntuple
from analysis_suite.plotting.hist_getter import GraphInfo, HistGetter
from analysis_suite.commons.plot_utils import nonratio_plot, cms_label

from analysis_suite.commons.constants import lumi

def plot_graphs(workdir, year):
    graphs_with = {
        "ht": GraphInfo(r"$H_T$", axis.Regular(25, 0, 1500), lambda vg : vg.get_hist('HT')),
        'met': GraphInfo(r"Met", axis.Regular(25, 0, 500), lambda vg : vg.get_hist('Met')),
        'nelec': GraphInfo('$N_e$', axis.Regular(4, 0, 4), lambda vg: (vg.TightElectron.num(), vg.scale)),
        'nmuon': GraphInfo('$N_\mu$', axis.Regular(4, 0, 4), lambda vg: (vg.TightMuon.num(), vg.scale)),
        # 'l1Pt': GraphInfo('$p_T(l_{{1}})$ (GeV)', axis.Regular(25, 0, 500), lambda vg: vg.TightLepton.get_hist('pt', 0)),
        # 'l2Pt': GraphInfo('$p_T(l_{{2}})$ (GeV)', axis.Regular(25, 0, 200), lambda vg: vg.TightLepton.get_hist('pt', 1)),
        # 'j1Pt': GraphInfo('$p_T(j_{{1}})$ (GeV)', axis.Regular(25, 0, 500), lambda vg: vg.Jets.get_hist('pt', 0)),
        # 'j2Pt': GraphInfo('$p_T(j_{{2}})$ (GeV)', axis.Regular(25, 0, 300), lambda vg: vg.Jets.get_hist('pt', 1)),
        # 'b1Pt': GraphInfo('$p_T(b_{{2}})$ (GeV)', axis.Regular(25, 0, 300), lambda vg: vg.BJets.get_hist('pt', 1)),
        'njet': GraphInfo(r"$N_j$", axis.Regular(9, 0, 9), lambda vg : (vg.Jets.num(), vg.scale)),
        'nbjet': GraphInfo(r"$N_b$", axis.Regular(5, 0, 5), lambda vg: vg.get_hist('NBjets_medium')),
        'nloose': GraphInfo(r"$N_{{loose}}$", axis.Regular(3, 0, 3), lambda vg: (vg['N_loose_el']+vg['N_loose_mu'], vg.scale)),
        'Zmass': GraphInfo(r"$N_{{loose}}$", axis.Regular(15, 75, 105), lambda vg: vg.get_hist('Zmass')),
        # 'eta_j': GraphInfo('eta', axis.Regular(25, -2.4, 2.4), lambda vg: vg.Jets.get_hist('eta', -1)),
        # 'pt_j': GraphInfo('pt', axis.Regular(25, 0, 500), lambda vg: vg.Jets.get_hist('pt', -1)),
        # 'phi_j': GraphInfo('phi', axis.Regular(25, -3.14, 3.14), lambda vg: vg.Jets.get_hist('phi', -1)),
        # 'eta_l': GraphInfo('eta', axis.Regular(25, -2.4, 2.4), lambda vg: vg.TightLepton.get_hist('eta', -1)),
        # 'pt_l': GraphInfo('pt', axis.Regular(25, 0, 500), lambda vg: vg.TightLepton.get_hist('pt', -1)),
        # 'phi_l': GraphInfo('phi', axis.Regular(25, -3.14, 3.14), lambda vg: vg.TightLepton.get_hist('phi', -1)),
        # 'mt': GraphInfo(r"$M_T$", axis.Regular(25, 0, 250), lambda vg: vg.TightLepton.get_hist('mt', -1)),
    }
    graphs_without = {
        "ht": GraphInfo(r"$H_T$", axis.Regular(25, 0, 1500), lambda vg : (vg['HT'], vg.get_sf(vg.syst_name)*vg['wgt_nobtag'])),
        'met': GraphInfo(r"Met", axis.Regular(25, 0, 500), lambda vg : (vg['Met'], vg.get_sf(vg.syst_name)*vg['wgt_nobtag'])),
        'nelec': GraphInfo('$N_e$', axis.Regular(4, 0, 4), lambda vg: (vg.TightElectron.num(), vg.get_sf(vg.syst_name)*vg['wgt_nobtag'])),
        'nmuon': GraphInfo('$N_\mu$', axis.Regular(4, 0, 4), lambda vg: (vg.TightMuon.num(), vg.get_sf(vg.syst_name)*vg['wgt_nobtag'])),
        # 'l1Pt': GraphInfo('$p_T(l_{{1}})$ (GeV)', axis.Regular(25, 0, 500), lambda vg: vg.TightLepton.get_hist('pt', 0)),
        # 'l2Pt': GraphInfo('$p_T(l_{{2}})$ (GeV)', axis.Regular(25, 0, 200), lambda vg: vg.TightLepton.get_hist('pt', 1)),
        # 'j1Pt': GraphInfo('$p_T(j_{{1}})$ (GeV)', axis.Regular(25, 0, 500), lambda vg: vg.Jets.get_hist('pt', 0)),
        # 'j2Pt': GraphInfo('$p_T(j_{{2}})$ (GeV)', axis.Regular(25, 0, 300), lambda vg: vg.Jets.get_hist('pt', 1)),
        # 'b1Pt': GraphInfo('$p_T(b_{{2}})$ (GeV)', axis.Regular(25, 0, 300), lambda vg: vg.BJets.get_hist('pt', 1)),
        'njet': GraphInfo(r"$N_j$", axis.Regular(9, 0, 9), lambda vg : (vg.Jets.num(), vg.get_sf(vg.syst_name)*vg['wgt_nobtag'])),
        'nbjet': GraphInfo(r"$N_b$", axis.Regular(5, 0, 5), lambda vg: (vg["NBjets_medium"], vg.get_sf(vg.syst_name)*vg['wgt_nobtag'])),
        'nloose': GraphInfo(r"$N_{{loose}}$", axis.Regular(3, 0, 3), lambda vg: (vg['N_loose_el']+vg['N_loose_mu'], vg.get_sf(vg.syst_name)*vg['wgt_nobtag'])),
        'Zmass': GraphInfo(r"$N_{{loose}}$", axis.Regular(15, 75, 105), lambda vg: (vg['Zmass'], vg.get_sf(vg.syst_name)*vg['wgt_nobtag'])),
        # 'eta_j': GraphInfo('eta', axis.Regular(25, -2.4, 2.4), lambda vg: vg.Jets.get_hist('eta', -1)),
        # 'pt_j': GraphInfo('pt', axis.Regular(25, 0, 500), lambda vg: vg.Jets.get_hist('pt', -1)),
        # 'phi_j': GraphInfo('phi', axis.Regular(25, -3.14, 3.14), lambda vg: vg.Jets.get_hist('phi', -1)),
        # 'eta_l': GraphInfo('eta', axis.Regular(25, -2.4, 2.4), lambda vg: vg.TightLepton.get_hist('eta', -1)),
        # 'pt_l': GraphInfo('pt', axis.Regular(25, 0, 500), lambda vg: vg.TightLepton.get_hist('pt', -1)),
        # 'phi_l': GraphInfo('phi', axis.Regular(25, -3.14, 3.14), lambda vg: vg.TightLepton.get_hist('phi', -1)),
        # 'mt': GraphInfo(r"$M_T$", axis.Regular(25, 0, 250), lambda vg: vg.TightLepton.get_hist('mt', -1)),
    }
    ntuple = get_ntuple('btag')
    ginfo = ntuple.get_info()
    plot_dir = workdir/ 'btag_scales'
    plot_dir.mkdir(exist_ok=True, parents=True)
    hist_factory = HistGetter(ntuple, year, workdir=workdir, scales=['btag_jetbinned'])
    hists_with = hist_factory.get_hists(graphs_with)
    hists_without = hist_factory.get_hists(graphs_without)
    for name, graph in graphs_without.items():
        h_with = hists_with[name]
        h_without = hists_without[name]
        for group in h_with.keys():
            scale = h_without[group].integral()/h_with[group].integral()
            print(group, name, scale)
            print(scale*h_with[group].vals)
            print(scale*h_with[group].vals/h_without[group].vals)

def process(workdir, year, nlep):
    ntuple = get_ntuple('btag')
    ginfo = ntuple.get_info()
    bins = axis.Regular(5, 2, 7)
    plot_dir = workdir/ 'btag_scales'
    plot_dir.mkdir(exist_ok=True, parents=True)
    graph = {
        'with': GraphInfo("$N_j$", bins, lambda vg: (vg.Jets.num(), vg.scale)),
        'without': GraphInfo("$N_j$", bins, lambda vg: (vg.Jets.num(), vg.get_sf(vg.syst_name)*vg['wgt_nobtag']))
    }
    out_jet = {group: dict() for group in ginfo.get_groups()}
    out_int = {group: dict() for group in ginfo.get_groups()}

    hist_factory = HistGetter(ntuple, year)
    if nlep == 2:
        hist_factory.cut(lambda vg: vg["TightLepton"].num() == 2)
    else:
        hist_factory.cut(lambda vg: vg["TightLepton"].num() > 2)

    for syst in hist_factory.systs:
        if syst != 'Nominal' and "BJet" not in syst and "Jet_JE" not in syst:
            continue
        hist_factory.reset_syst(syst)
        hists = hist_factory.get_hists(graph, fix_negative=True)
        print(syst)
        with_btag = hists['with']
        without_btag = hists['without']
        for group in ginfo.get_groups():
            out_jet[group][syst] = without_btag[group]/with_btag[group]
            print(syst, without_btag[group].integral(), with_btag[group].integral())
            print(without_btag[group].integral()/with_btag[group].integral())
            out_int[group][syst] = without_btag[group].integral()/with_btag[group].integral()
            # print(group, without_btag[group].integral()/with_btag[group].integral())
            # print((without_btag[group]/with_btag[group]).vals)

        # continue
        if syst != "Nominal":
            continue
        with nonratio_plot(plot_dir/f"bjet_scales_{year}.pdf", "$N_{j}$", bins.edges,
                           pad_label="$w_{w/o\ btag}/w_{w/\ btag}$", zero_bot=False) as ax:
            for group, hist in out_jet.items():
                hist = hist['Nominal']
                hist.color = ginfo.get_color(group)
                hist.plot_label = ginfo.get_legend_name(group)
                hist.plot_points(ax)
            ax.legend()
            cms_label(ax, year=year)

    # exit()
    # Dump MC scale factors
    scale_file = workdir/f"btag_scales_lep{nlep}.pkl"
    if scale_file.exists():
        with open(scale_file, "rb") as f:
            scales = pickle.load(f)
    else:
        scales = dict()
    scales[year] = out_int
    with open(scale_file, "wb") as f:
        pickle.dump(scales, f)

    # Jet binned
    scale_file = workdir/f"btag_scales_njet_lep{nlep}.pkl"
    if scale_file.exists():
        with open(scale_file, "rb") as f:
            scales = pickle.load(f)
    else:
        scales = dict()
    scales[year] = out_jet
    with open(scale_file, "wb") as f:
        pickle.dump(scales, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : all_eras if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-d', '--workdir', required=True,
                        help="directory to run over. If nothing, use date",)
    args = parser.parse_args()

    workdir = workspace_area / args.workdir
    workdir.mkdir(exist_ok=True, parents=True)

    for year in args.years:
        print(year)
        process(workdir, year, nlep=2)
        process(workdir, year, nlep=3)
        # plot_graphs(workdir, year)
