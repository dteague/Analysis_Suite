#!/usr/bin/env python3
import argparse
import boost_histogram.axis as axis
import uproot
import numpy as np

import analysis_suite.commons.user as user
from analysis_suite.commons.constants import all_eras
from analysis_suite.commons.info import GroupInfo
from analysis_suite.commons.configs import get_list_systs, get_ntuple_info, get_graph
from analysis_suite.combine.card_maker import Card_Maker
from analysis_suite.combine.hist_writer import HistWriter
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.flatten.flatgetter import FlatGetter
from analysis_suite.combine.combine_wrapper import runCombine
from analysis_suite.plotting.plotter import GraphInfo
from analysis_suite.data.systs import systematics, get_shape_systs, dummy
from analysis_suite.commons.histogram import Histogram

def order_list(groups, sig):
    glist = list(groups.keys())
    if glist[0] != sig:
        idx = glist.index(sig)
        glist[0], glist[idx] = glist[idx], glist[0]
    return glist


def get_sig(vg):
    j, b = np.array(vg["NJets"]), vg['NmediumBJets']
    masker = (j >= 8) + (b >= 5) + (j >= 7)*(b >= 3)
    return np.where(masker, 1.01, vg["Signal"]), vg.scale

# def get_systs(avail_systs):
#     syst_objs = list()
#     for syst in systematics:
#         if f'{syst.name}_up' in avail_systs and f'{syst.name}_down' in avail_systs:
#             syst_objs.append(syst)
#         if syst.syst_type == "lnN":
#             syst_objs.append(syst)
#     return syst_objs

def get_systs(filename):
    with uproot.open(filename) as f:
        pass

def fill_group(syst_hists, f, year, group, members, graph, mask=None):
    for member in members:
        fg = FlatGetter(f, member)
        if not fg:
            continue
        if mask is not None:
            fg.mask = mask
        for syst, final_name in get_shape_systs(year):
            if not fg.includes_systs(syst):
                continue
            if final_name not in syst_hists:
                syst_hists[final_name] = dict()
            if group not in syst_hists[final_name]:
                syst_hists[final_name][group] = Histogram("", *graph.bins())
            fg.set_syst(syst)
            vals, weight = fg.get_graph(graph)
            # scale = 1000.
            # if member == "tttw" or member == "tttj" or member == "tttt":
            #     weight = scale*weight
            syst_hists[final_name][group].fill(vals, weight=weight, member=member)
    return syst_hists

def make_card(workdir, combinedir, year, infile, ntupleName, region, graph, unblind=False, **kwargs):
    remove = ['nonprompt_mc']
    if region == "Multi":
        remove.append('charge_flip')
    ginfo = get_ntuple_info(ntupleName, remove=remove)
    groups = ginfo.setup_groups()
    # groups['4top'] = groups.pop('tttt')
    graph_name = graph.name

    if not kwargs.get("skip", True):
        file_systs = get_list_systs(workdir, 'combine', ['all'])
        syst_hists = {}
        for file_syst in file_systs:
            with uproot.open(workdir/year/infile.format(file_syst)) as f:
                for group, members in groups.items():
                    mask = kwargs.get("mask", None)
                    syst_hists = fill_group(syst_hists, f, year, group, members, graph, mask)

        with HistWriter(combinedir / f'{graph_name}_{year}_{region}.root') as writer:
            for syst, hists in syst_hists.items():
                writer.add_syst(hists, syst=syst, blind=not unblind)

    # signals = ['ttt', 'tttt']
    # signals = ['tttj', "tttw", "4top",]
    signals = ['ttt']

    group_list = [g for g in groups.keys() if g not in signals]
    nosyst = kwargs.get("no_systs", False)
    with Card_Maker(combinedir, year, region, signals, group_list, graph_name, nosyst=nosyst) as card:
        card.write_preamble()
        if not nosyst:
            card.write_systematics(systematics)
        else:
            card.write_systematics(dummy)
        card.add_stats()
        for rate_param in kwargs.get("rate_params", []):
            card.add_rateParam(rate_param)
    if nosyst:
        region += "_nosyst"
    return f"{region}={graph_name}_{year}_{region}_card.txt "


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', "--nominal", action="store_true", help="Run with no systematics (only nominal)")
    parser.add_argument("-y", "--years", required=True, type=lambda x : all_eras if x == "all" else x.split(','),
                        help="Year to use")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area / x,
                        help="Working Directory")
    parser.add_argument("-t", '--extra_text', default="")
    parser.add_argument('-s', '--skip_hist', action="store_true")
    parser.add_argument('-ns', '--no_syst', action='store_true')
    args = parser.parse_args()
    combine_dir = args.workdir/"combine"/args.extra_text
    combine_dir.mkdir(exist_ok=True, parents=True)
    runCombine.work_dir = combine_dir

    signal_graph = GraphInfo('signal', '', axis.Variable([0.0, 0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.94, 1.0, 1.1]), lambda x : get_sig(x))
    ttz_graph = GraphInfo('ZMass', '', axis.Regular(15, 75, 105), "ZMass")
    # rate_params = ['ttz']
    rate_params = []
    all_command = 'combineCards.py '

    for year in args.years:
        print("Start:", year)
        combine_cmd = "combineCards.py "
        combine_cmd += make_card(args.workdir, combine_dir, year, 'test_{}_signal.root', 'signal', 'Dilepton', signal_graph,
                                 mask=lambda vg : vg['NMuons']+vg['NElectrons'] == 2, rate_params=rate_params, skip=args.skip_hist, no_systs=args.no_syst)
        combine_cmd += make_card(args.workdir, combine_dir, year, 'test_{}_signal.root', 'signal', 'Multi', signal_graph,
                                 mask=lambda vg : vg['NMuons']+vg['NElectrons'] > 2, rate_params=rate_params, skip=args.skip_hist, no_systs=args.no_syst)
        # combine_cmd += make_card(args.workdir, combine_dir, year, 'processed_{}_ttzCR.root', 'ttzCR', 'ttz', ttz_graph,
        #                          rate_params=rate_params, skip=args.skip_hist)
        final_card = f"final_{year}" + ("_nosyst" if args.no_syst else "") + "_card.txt"
        combine_cmd += f'> {final_card}'
        all_command += f'era{year}={final_card} '
        runCombine(combine_cmd)


    if args.years == all_eras:
        final_card = f"final_all" + ("_nosyst" if args.no_syst else "") + "_card.txt"
        all_command += f"> {final_card}"
        runCombine(all_command)
