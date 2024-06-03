#!/usr/bin/env python3
import argparse
import boost_histogram.axis as axis
import uproot
import numpy as np
import collections.abc
import multiprocessing as mp
from copy import deepcopy
import warnings

import analysis_suite.commons.user as user
from analysis_suite.commons.constants import all_eras
from analysis_suite.commons.configs import get_list_systs, get_ntuple_info, get_graph, get_inputs
from analysis_suite.combine.card_maker import Card_Maker
from analysis_suite.combine.hist_writer import HistWriter
from analysis_suite.flatten.flatgetter import FlatGetter
from analysis_suite.combine.combine_wrapper import runCombine
from analysis_suite.plotting.plotter import GraphInfo
from analysis_suite.data.systs import systematics, get_shape_systs, dummy
from analysis_suite.commons.histogram import Histogram

from analysis_suite.commons.constants import lumi
from analysis_suite.commons.plot_utils import plot, nonratio_plot, ratio_plot, hep
from analysis_suite.plotting.stack import Stack

masks = {}
graphs = {}
ginfo = None
signals = ['tttj', 'tttw']

def plot_stack(hists, outfile, graph, year, region, normed=False):
    signal = Histogram("ttt", graph.bin_tuple, color='crimson', name=ginfo.get_legend_name('ttt'))
    data = Histogram("data", graph.bin_tuple, color='black')
    stack = Stack(*graph.bins())

    for group, hist in hists.items():
        hist.set_plot_details(ginfo)
        if group in signals and region != 'ttzCR':
            signal += hist
        elif group == 'data':
            data += hist
        else:
            stack += hist

    plotter = ratio_plot if data else nonratio_plot
    with plotter(outfile, graph.axis_name, stack.get_xrange(), normed=normed) as ax:
        ratio = Histogram("Ratio", graph.bin_tuple, color="black")
        band = Histogram("Ratio", graph.bin_tuple, color="plum")
        error = Histogram("Stat Errors", graph.bin_tuple, color="plum")

        if data:
            ratio += data/stack
        band += stack/stack
        for hist in stack.stack:
            error += hist

        if data:
            pad, subpad = ax
        else:
            pad, subpad = ax, None

        #upper pad
        if region != 'ttzCR' and region != 'ttttCR':
            rounder = 50
            sig_scale = np.max(stack.vals)/np.max(signal.vals)
            sig_scale = int(sig_scale//rounder*rounder)
            signal.scale(sig_scale, changeName=True, forPlot=True)

        stack.plot_stack(pad, normed=normed)
        signal.plot_points(pad, normed=normed)
        data.plot_points(pad, normed=normed)
        error.plot_band(pad, normed=normed)

        # ratio pad
        ratio.plot_points(subpad)
        band.plot_band(subpad)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hep.cms.label(ax=pad, lumi=lumi[year], data=data)


def get_syst_name(systname, group, year):
    if systname == 'Nominal':
        return 'Nominal'
    updown_break = systname.rfind('_')
    if updown_break == -1:
        print('ERROR: Systematic does not have _up/_down in name --', systname)
        return None
    updown = systname[updown_break:]
    rawsyst = systname[:updown_break]
    for syst in systematics:
        if syst.name == rawsyst and syst.good_syst(group, year):
            return syst.get_name(year)+updown

    return None


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
            continue
        if mask is not None:
            fg.mask = mask
        for syst in f[f'weights/{member}'].keys():
            if syst == 'index':
                continue
            final_name = get_syst_name(syst, group, year)
            if final_name is None:
                continue
            fg.set_syst(syst)
            if abs(sum(fg.scale)) < 1e-5:
                continue
            if final_name not in output:
                output[final_name] = {group: Histogram(group, *graph.bins())}
            vals, weight = fg.get_graph(graph)
            output[final_name][group].fill(vals, weight=weight, member=member)
    return output


def make_hists(file_list, outdir, year, region, unblind):
    groups = ginfo.setup_groups()
    mask = masks.get(region, None)
    graph = graphs.get(region, None)
    graph_name = graph.name
    out_signals = [ginfo.get_combine_name(s) for s in signals]
    group_list = list()

    syst_hists = {}
    for infile in file_list:
        # if "Nominal" not in str(infile):
        #     continue
        with uproot.open(infile) as f:
            for group, members in groups.items():
                syst_hists = deep_update(syst_hists, fill_group(f, year, group, members, graph, mask))

    clean_negatives(syst_hists)

    plot_stack(syst_hists['Nominal'], outdir/'plots'/f'{graph_name}_{year}_{region}.png', graph, year, region)
    plot_stack(syst_hists['Nominal'], outdir/'plots'/f'{graph_name}_{year}_{region}_normed.png', graph, year, region, normed=True)

    with HistWriter(outdir / f'{graph_name}_{year}_{region}.root') as writer:
        for syst, hists in syst_hists.items():
            final_name = syst[:syst.rfind('_')]
            writer.add_syst(hists, ginfo, syst=syst, blind=not unblind)
    for g in syst_hists['Nominal'].keys():
        if g not in signals and 'data' not in g:
            group_list.append(ginfo.get_combine_name(g))
    print(f'Finished {region}: {year}')


def make_card(outdir, year, region, graph, rate_params, nosyst):
    groups = ginfo.setup_groups()
    graph_name = graph.name
    out_signals = [ginfo.get_combine_name(s) for s in signals]
    group_list = list()

    with uproot.open(outdir / f'{graph_name}_{year}_{region}.root') as f:
        used_syst_names = list()
        for d, cls in f.classnames().items():
            if '/' in d:
                continue
            elif 'TH1' in cls and d[:-2] not in out_signals and 'data' not in d:
                group_list.append(d[:-2])
            elif 'Up' in d:
                used_syst_names.append(d[:-4])

    def keep_systs(x):
        return x.syst_type != 'shape' or x.get_name(year) in used_syst_names

    with Card_Maker(outdir, year, region, out_signals, group_list, graph_name, nosyst) as card:
        card.write_preamble()
        if nosyst:
            card.write_systematics(dummy, ginfo)
        else:
            card.write_systematics(list(filter(keep_systs, systematics)), ginfo)
            # card.add_stats()
        for rp in rate_params:
            card.add_rateParam(ginfo.get_combine_name(rp))

    return f"{region}={graph_name}_{year}_{region}{'_nosyst' if args.no_systs else ''}_card.txt"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', "--nominal", action="store_true", help="Run with no systematics (only nominal)")
    parser.add_argument("-y", "--years", required=True, type=lambda x : all_eras if x == "all" else x.split(','),
                        help="Year to use")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area / x,
                        help="Working Directory")
    parser.add_argument("-t", '--extra_text', default="")
    parser.add_argument("-u", '--unblind', action='store_true')
    parser.add_argument("-j", '--cores', default=1, type=int)
    parser.add_argument('-ns', '--no_systs', action='store_true')
    parser.add_argument('--skip', action='store_true')
    args = parser.parse_args()
    combine_dir = args.workdir/"combine"/args.extra_text
    combine_dir.mkdir(exist_ok=True, parents=True)
    (combine_dir/'plots').mkdir(exist_ok=True, parents=True)
    runCombine.work_dir = combine_dir

    all_command = 'combineCards.py '

    ginfo = get_ntuple_info('signal', remove=['nonprompt_mc', 'ttt'], add={'data': 'black'})

    combine_info = get_inputs(args.workdir, 'combine_info')
    rate_params = combine_info.rate_params

    if not args.skip:
        hist_inputs = []
        for region, info in combine_info.regions.items():
            masks[region] = info['mask']
            graphs[region] = info['graph']
            for year in args.years:
                infiles = list((args.workdir/info['dir']/year).glob(info['glob']))
                hist_inputs.append(( infiles, combine_dir, year, region, args.unblind))

        if args.cores == 1:
            for input in hist_inputs:
                make_hists(*input)
        else:
            with mp.Pool(args.cores) as pool:
                pool.starmap(make_hists, hist_inputs)

    for year in args.years:
        combine_cmd = "combineCards.py"
        final_card = f"final_{year}{'_nosyst' if args.no_systs else ''}_card.txt"

        combine_txt = list()
        for region, info in combine_info.regions.items():
            combine_cmd += " " + make_card(combine_dir, year, region, info['graph'], rate_params, args.no_systs)

        combine_cmd += f' > {final_card}'
        all_command += f'era{year}={final_card} '
        runCombine(combine_cmd)

        comb_card = f"brown_wisc_{year}_card.txt"
        brown_year = year[2:4]
        if year == '2016pre':
            brown_year += "APV"
        brown_card = user.analysis_area/'daniel_cards'/f'workspace_{brown_year}.txt'

        runCombine(f'combineCards.py brown={brown_card} wisc={final_card} > {comb_card}')


    if args.years == all_eras:
        print("Combining all cards")
        final_card = f"final_all_card.txt"
        all_command += f"> {final_card}"
        runCombine(all_command)
