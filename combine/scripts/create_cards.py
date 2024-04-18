#!/usr/bin/env python3
import argparse
import boost_histogram.axis as axis
import uproot
import numpy as np
import collections.abc
import multiprocessing as mp

import analysis_suite.commons.user as user
from analysis_suite.commons.constants import all_eras
from analysis_suite.commons.configs import get_list_systs, get_ntuple_info, get_graph
from analysis_suite.combine.card_maker import Card_Maker
from analysis_suite.combine.hist_writer import HistWriter
from analysis_suite.flatten.flatgetter import FlatGetter
from analysis_suite.combine.combine_wrapper import runCombine
from analysis_suite.plotting.plotter import GraphInfo
from analysis_suite.data.systs import systematics, get_shape_systs, dummy, get_syst_name
from analysis_suite.commons.histogram import Histogram

masks = {}

def deep_update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = deep_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def order_list(groups, sig):
    glist = list(groups.keys())
    if glist[0] != sig:
        idx = glist.index(sig)
        glist[0], glist[idx] = glist[idx], glist[0]
    return glist

def get_systs(filename):
    with uproot.open(filename) as f:
        pass

def clean_negatives(syst_hists):
    for systs, hists in syst_hists.items():
        for group, hist in hists.items():
            if np.sum(hist.vals) < 0:
                hist.hist.view().value = np.abs(hist.vals)

def fill_group(f, year, group, members, graph, mask=None):
    output = {}
    for member in members:
        fg = FlatGetter(f, member)
        if not fg:
            # print(member, 'not found!')
            continue
        if mask is not None:
            fg.mask = mask
        for syst in f[f'weights/{member}'].keys():
            if syst == 'index':
                continue
            final_name = get_syst_name(syst, year)
            fg.set_syst(syst)
            if abs(sum(fg.scale)) < 1e-5:
                continue
            if final_name not in output:
                output[final_name] = {group: Histogram("", *graph.bins())}
            vals, weight = fg.get_graph(graph)
            output[final_name][group].fill(vals, weight=weight, member=member)
    return output


def make_card(indir, filename, outdir, year, groups, region, graph, rate_params, unblind):
    graph_name = graph.name
    signals = ['ttt']
    mask = masks.get(region, None)

    syst_hists = {}
    for infile in (indir/year).glob(filename):
        with uproot.open(infile) as f:
            for group, members in groups.items():
                syst_hists = deep_update(syst_hists, fill_group(f, year, group, members, graph, mask))
    clean_negatives(syst_hists)

    used_syst_names = set()
    with HistWriter(outdir / f'{graph_name}_{year}_{region}.root') as writer:
        for syst, hists in syst_hists.items():
            final_name = syst[:syst.rfind('_')]
            if final_name not in used_syst_names:
                used_syst_names.add(final_name)
            writer.add_syst(hists, syst=syst, blind=not unblind)

    def keep_systs(x):
        return x.syst_type != 'shape' or x.dan_name in used_syst_names

    group_list = [g for g in syst_hists['Nominal'].keys() if g not in signals]
    with Card_Maker(outdir, year, region, signals, group_list, graph_name) as card:
        card.write_preamble()
        card.write_systematics(list(filter(keep_systs, systematics)))
        card.add_stats()
        for rp in rate_params:
            card.add_rateParam(rp)

    return f"{region}={graph_name}_{year}_{region}_card.txt"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', "--nominal", action="store_true", help="Run with no systematics (only nominal)")
    parser.add_argument("-y", "--years", required=True, type=lambda x : all_eras if x == "all" else x.split(','),
                        help="Year to use")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area / x,
                        help="Working Directory")
    parser.add_argument("-t", '--extra_text', default="")
    parser.add_argument("-u", '--unblind', action='store_true')
    args = parser.parse_args()
    combine_dir = args.workdir/"combine"/args.extra_text
    combine_dir.mkdir(exist_ok=True, parents=True)
    runCombine.work_dir = combine_dir

    sig_bins = [0.0, 0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.94, 1.0]
    signal_graph = GraphInfo('sig_3top', '', axis.Variable(sig_bins), '3top_sig')
    top4_graph = GraphInfo('sig_4top', '', axis.Regular(1, 0, 1), '4top_sig')
    ttz_graph = GraphInfo('HT', '', axis.Regular(20, 250, 750), "HT")

    rate_params = ['ttz']
    all_command = 'combineCards.py '

    groups = get_ntuple_info('signal', remove=['nonprompt_mc']).setup_groups()

    masks['Dilepton'] = lambda vg : vg['NMuons']+vg['NElectrons'] == 2
    masks['Multi'] = lambda vg : vg['NMuons']+vg['NElectrons'] == 3
    masks['ttttCR'] = lambda vg : vg['4top_sig'] > 0.97

    for year in args.years:
        print("Start:", year)
        combine_cmd = "combineCards.py "
        inputs = [
            (args.workdir/'second_train', 'test*', combine_dir, year, groups, 'Dilepton',
             signal_graph, rate_params, args.unblind),
            (args.workdir/'second_train', 'test*', combine_dir, year, groups, 'Multi',
             signal_graph, rate_params, args.unblind),
            (args.workdir/'first_train', 'test*', combine_dir, year, groups, 'ttttCR',
             top4_graph, rate_params, args.unblind),
            (args.workdir, 'processed*ttzCR.root', combine_dir, year, groups, 'ttzCR',
             ttz_graph, rate_params, args.unblind),
        ]

        with mp.Pool(len(inputs)) as pool:
            combine_txt = pool.starmap(make_card, inputs)
        combine_cmd += " ".join(combine_txt)

        final_card = f"final_{year}_card.txt"
        combine_cmd += f'> {final_card}'
        all_command += f'era{year}={final_card} '
        runCombine(combine_cmd)


    if args.years == all_eras:
        print("Combining all cards")
        final_card = f"final_all_card.txt"
        all_command += f"> {final_card}"
        runCombine(all_command)
