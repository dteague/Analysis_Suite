#!/usr/bin/env python3
import argparse
import boost_histogram.axis as axis
import uproot
import numpy as np
import multiprocessing as mp
import warnings

warnings.simplefilter("ignore", UserWarning)
from statsmodels.nonparametric.smoothers_lowess import lowess

import analysis_suite.commons.user as user
from analysis_suite.commons.constants import all_eras
from analysis_suite.commons.configs import get_ntuple, get_ntuple_info, get_inputs
from analysis_suite.combine.card_maker import Card_Maker
from analysis_suite.combine.hist_writer import HistWriter
from analysis_suite.combine.combine_wrapper import runCombine
from analysis_suite.plotting.hist_getter import GraphInfo, HistGetter
from analysis_suite.data.systs import systematics, dummy
from analysis_suite.combine.systematics import use_lowess
from analysis_suite.commons.histogram import Histogram
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.plotting.LogFile import LogFile

from analysis_suite.commons.constants import lumi

def deep_update(mapping, *updating_mappings):
    updated_mapping = mapping.copy()
    for updating_mapping in updating_mappings:
        for k, v in updating_mapping.items():
            if k in updated_mapping and isinstance(updated_mapping[k], dict) and isinstance(v, dict):
                updated_mapping[k] = deep_update(updated_mapping[k], v)
            else:
                updated_mapping[k] = v
    return updated_mapping

signals = ['ttt_nlo']

dilep_bins = [0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.9]
multi_bins = [0, 0.24, 0.36, 0.48, 0.6, 0.9]
# multi_bins = [0.0, 0.3, 0.45, 0.6, 0.94]
# multi_bins = np.linspace(0, 1, 26)
ht_bins = [250, 275, 300, 350, 400, 450, 550, 700, 850]
# ht_bins = np.linspace(250, 1000, 31)
dilep_graph = GraphInfo(r'$BDT_{{3top}}$', axis.Variable(dilep_bins), 'BDT')
multi_graph = GraphInfo(r'$BDT_{{3top}}$',axis.Variable(multi_bins), 'BDT')
top4_graph = GraphInfo(r'$BDT_{{4top}}$',axis.Regular(1, 0, 1), 'BDT')
ttz_graph = GraphInfo(r'$H_{{T}}$ (GeV)', axis.Variable(ht_bins), "HT")
njet_graph = GraphInfo(r'$N_j$', axis.Regular(6, 2, 8), "NJets")
ttz_valid_graph = GraphInfo(r'$BDT_{{3top}}$',axis.Regular(15,0,1), 'BDT')

masks = {
    'Dilepton': lambda vg : (vg['NLeps'] == 2),
    'Multi': lambda vg : (vg['NLeps'] > 2),
}
all_graphs = {
    # 'Dilepton': {"combine": dilep_graph},
    # 'Multi': {"combine": multi_graph},
    # 'ttzCR': {"combine": ttz_graph},
    # 'ttttCR': {"combine": top4_graph},
    'ttz4Valid': {"combine": ttz_valid_graph},
    'ttz3Valid': {"combine": ttz_valid_graph},
}


def split_hists(hists, systname, year, chan):
    if systname != 'Nominal':
        updown_break = systname.rfind('_')
        rawsyst, updown = systname[:updown_break], systname[updown_break+1:]
    output = {}
    for group, hist in hists.items():
        if systname == "Nominal":
            yield systname, group, hist
        else:
            for syst in systematics:
                if syst.name == rawsyst and syst.good_syst(group, year, chan):
                    out_name = syst.get_name(year, with_lowess=False)
                    yield f'{out_name}_{updown}', group, hist

def get_hists(infile, systs, workdir, ntuple_name, region, year, unblind):
    ntuple = get_ntuple(*ntuple_name)
    ntuple.remove_group("nonprompt_mc")
    if not unblind:
        ntuple.remove_group('data')
    mask = masks.get(region, None)
    graphs = all_graphs.get(region, None)
    syst_hists = {name: {} for name in graphs.keys()}

    # Hack to get things working
    if region == 'ttttCR':
        cut = get_inputs(workdir, 'params').cut
        mask = lambda vg: vg['BDT'] > cut[year]

    hist_factory = HistGetter(ntuple, year, region=region, filename=infile,
                                workdir=workdir, scales=['btag_jetlep', 'wz', 'theory_rescale'], mask=mask)
    systs = hist_factory.systs if systs is None else systs
    for syst in systs:
        if "data" in str(infile) and syst != "Nominal":
            continue
        if "Nominal" not in infile.name and syst not in infile.name:
            continue
        hist_factory.reset_syst(syst)
        all_hists = hist_factory.get_hists(graphs, fix_negative=True)
        for graph_name, hists in all_hists.items():
            for systname, group, hist in split_hists(hists, syst, year, region):
                if systname not in syst_hists[graph_name]:
                    syst_hists[graph_name][systname] = {}
                syst_hists[graph_name][systname][group] = hist
    return {f'{region}-{year}': syst_hists}

def smooth_hist(nom, up, down, frac=0.67, it=5, symm=False):
    centers = nom.axis.centers
    if symm:
        up_ratio = 1+(up.vals-down.vals)/(2*nom.vals+1e-5)
        down_ratio = 1+(down.vals-up.vals)/(2*nom.vals+1e-5)
    else:
        up_ratio = (up.vals+1e-5)/(nom.vals+1e-5)
        down_ratio = (down.vals+1e-5)/(nom.vals+1e-5)

    if len(centers) > 2:
        up_ratio_lowess = lowess(up_ratio, centers, frac=frac, it=it).T[1]
        down_ratio_lowess = lowess(down_ratio, centers, frac=frac, it=it).T[1]
    else:
        up_ratio_lowess = up_ratio
        down_ratio_lowess = down_ratio

    up_lowess = Histogram(nom.axis)
    up_lowess.set_data(up_ratio_lowess*nom.vals, up.variances())
    down_lowess = Histogram(nom.axis)
    down_lowess.set_data(down_ratio_lowess*nom.vals, down.variances())

    if symm:
        up.set_data(up_ratio*nom.vals)
        down.set_data(down_ratio*nom.vals)
    return up_lowess, down_lowess


def write_hists(graph_name, syst_hists, region, year, outdir, unblind):
    nom_hist = syst_hists.pop('Nominal')
    ginfo = get_ntuple_info('signal')

    for syst in systematics:
        if not use_lowess(syst.name):
            continue
        unsmooth, smooth = syst.get_name(year, with_lowess=False), syst.get_name(year)
        if f'{unsmooth}_up' not in syst_hists:
            continue
        up_unsmooth = syst_hists[f'{unsmooth}_up']
        down_unsmooth = syst_hists[f'{unsmooth}_down']
        syst_hists[smooth+'_up'] = {}
        syst_hists[smooth+'_down'] = {}
        for group in up_unsmooth.keys():
            up, down = smooth_hist(nom_hist[group], up_unsmooth[group], down_unsmooth[group], symm=True, frac=0.67)
            syst_hists[smooth+'_up'][group] = up
            syst_hists[smooth+'_down'][group] = down

    with HistWriter(outdir / f'{graph_name}_{year}_{region}.root') as writer:
        writer.add_syst(nom_hist, ginfo, syst="Nominal", blind=not unblind)
        for syst, hists in syst_hists.items():
            if "JEC" not in syst or "JER" not in syst:
                syst = syst.replace('LOWESS', "")
            writer.add_syst(hists, ginfo, syst=syst, blind=not unblind)

    logger = LogFile(graph_name, lumi[year], graph_name)
    for group, hist in nom_hist.items():
        logger.add_breakdown(group, hist)
        if group in signals:
            logger.add_mc(group, hist, 'signal')
        elif group == 'data':
            logger.add_data(hist)
        else:
            logger.add_mc(group, hist, "bkg")
    logger.write_out(outdir/'plots', f'{graph_name}_{region}_{year}')
    print(f'Finished {region}: {year}')

def make_card(outdir, year, region, graph_name, rate_params, nosyst):
    ginfo = get_ntuple_info('signal')
    groups = ginfo.setup_groups()
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
            card.add_stats()
        for rp in rate_params:
            card.add_rateParam(ginfo.get_combine_name(rp))

    if nosyst:
        region += "_nosyst"
    return outdir/f'{graph_name}_{year}_{region}_card.txt'
    # return f"{region}={graph_name}_{year}_{region}{'_nosyst' if args.no_systs else ''}_card.txt"

def plot_prefit(ntuple, workdir, datacard, graph):
    combine_dir = workdir/'tmp'
    combine_dir.mkdir(exist_ok=True, parents=True)
    # shutils.copy(datacard, combine_dir)
    ginfo = ntuple.get_info()
    _, region, year, graph_name = [word[::-1] for word in datacard.stem[::-1].split('_', 3)]
    bins = graph.bins()[0]
    axis_name = graph.axis_name
    nbins = bins.size

    command = 'combine -M FitDiagnostics --skipSBFit --skipBOnlyFit --saveShapes --saveWithUncertainties'
    fit_file = combine_dir/f'fitDiagnostics.{graph_name}.{year}.root'
    if not fit_file.exists() or datacard.stat().st_mtime > fit_file.stat().st_mtime:
        runCombine(f"{command} -d {datacard} -n .{graph_name}.{year}", output=False, workdir=combine_dir)

    # ntuple = get_ntuple('signal')
    plotter = Plotter(ntuple, year, sig='ttt_nlo', bkg='all', outdir=workdir, from_ntuple=False)

    with uproot.open(fit_file) as f:
        f = f['shapes_prefit']
        hists = {}
        syst_err = None
        for reg in f.keys(cycle=False, recursive=False):
            for combine_group, hist in f[reg].items(cycle=False):
                if combine_group == 'total_background':
                    if syst_err is None:
                        syst_err = hist.variances()
                    else:
                        syst_err += hist.variances()
                    continue
                elif "total" in combine_group:
                    continue
                elif "TGraph" in repr(hist):
                    vals, sumw2 = hist.values()[1], hist.errors('mean')[1]
                    sumw2 = sumw2**2
                else:
                    vals, sumw2 = hist.values(), hist.variances()
                group = ginfo.get_group_from_combine(combine_group)
                tmp_hist = Histogram(bins, axis_name=axis_name)
                tmp_hist.set_data(val=vals[:nbins], sumw2=sumw2[:nbins])
                if group in hists:
                    hists[group] += tmp_hist
                else:
                    hists[group] = tmp_hist
        plotter.set_extra_text(f'{region}_prefit')
        plotter.plot_hist(graph_name, hists, syst_err=syst_err)

ntuple_infos = {
    # "Dilepton": ('signal', "dilep_ntuple"),
    # "Multi": ('signal', "multi_ntuple"),
    # "ttttCR": ('signal', "info"),
    # "ttzCR": ('ttzCR', "info"),
    "ttz4Valid": ('ttzCR', "info"),
    "ttz3Valid": ('ttzCR', "info"),
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", type=lambda x : x.split(','), help="Year to use")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area / x,
                        help="Working Directory")
    parser.add_argument("-n", '--ntuple', default="signal")
    parser.add_argument("-f", '--flat', action='store_true')
    parser.add_argument('-ns', '--no_systs', action='store_true')
    parser.add_argument('-p', '--prefit', action='store_true')
    parser.add_argument("-u", '--unblind', action='store_true')
    parser.add_argument('--skip', action='store_true')
    parser.add_argument("-j", '--cores', default=1, type=int)
    parser.add_argument("-t", '--extra_text', default="")
    args = parser.parse_args()

    combine_dir = args.workdir/"combine"/args.extra_text
    combine_dir.mkdir(exist_ok=True, parents=True)
    (combine_dir/'plots').mkdir(exist_ok=True, parents=True)
    runCombine.work_dir = combine_dir

    combine_info = get_inputs(args.workdir, 'combine_info')
    hist_inputs = []

    if not args.skip:
        nom_systs = [syst.name for syst in systematics if syst.syst_type == "shape" and "Jet_JE" not in syst.name]
        nom_systs = np.vstack((np.char.add(nom_systs, "_up"), np.char.add(nom_systs, "_down")))
        nom_systs = np.unique(np.concatenate((nom_systs.flatten(),  ["Nominal"])))

        for region, info in combine_info.regions.items():
            ntuple_name = ntuple_infos[region]
            for year in args.years:
                file_dir = args.workdir/info['dir']
                if 'year_split' in info and info['year_split']:
                    file_dir /= year
                for filename in file_dir.glob(info['glob']):
                    if 'Nominal' in filename.name and args.cores > 1:
                        for syst in nom_systs:
                            hist_inputs.append((filename, [syst], args.workdir, ntuple_name, region, year, args.unblind))
                    else:
                        hist_inputs.append((filename, None, args.workdir, ntuple_name, region, year, args.unblind))

        if args.cores == 1:
            file_hists = [get_hists(*input) for input in hist_inputs]
        else:
            with mp.Pool(args.cores) as pool:
                file_hists = pool.starmap(get_hists, hist_inputs)

        all_hists = dict()
        for dicts in file_hists:
            all_hists = deep_update(all_hists, dicts)

        print(all_hists)
        for region, graph_hists in all_hists.items():
            region, year = region.split("-")
            for graph_name, hists in graph_hists.items():
                write_hists(graph_name, hists, region, year, combine_dir, args.unblind)

    # Make the cards
    rate_params = combine_info.rate_params
    all_command = 'combineCards.py '
    for year in args.years:
        combine_cmd = "combineCards.py"
        for region, info in combine_info.regions.items():
            ntuple = get_ntuple(*ntuple_infos[region])
            for graph_name, graph in all_graphs[region].items():
                card = make_card(combine_dir, year, region, graph_name, rate_params, args.no_systs)
                if args.prefit:
                    plot_prefit(ntuple, combine_dir/'plots', card, graph)
                if graph_name == "combine":
                    combine_cmd += f" {region}={card}"

        if '=' in combine_cmd:
            final_card = f"final_{year}{'_nosyst' if args.no_systs else ''}_card.txt"
            combine_cmd += f' > {final_card}'
            all_command += f'era{year}={final_card} '
            runCombine(combine_cmd)

    if args.years == all_eras and "=" in all_command:
        print("Combining all cards")
        full_command = all_command

        final_card = f"final_all_card.txt"
        all_command += f" > {final_card}"
        runCombine(all_command)

        final_brown_card = 'brown_wisc_all_card.txt'
        brown_card = user.analysis_area/'daniel_cards'/f'workspace.txt'
        full_command += f'brown={brown_card} '
        full_command += f"> {final_brown_card}"
        runCombine(full_command)
