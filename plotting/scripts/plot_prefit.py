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
from prettytable import PrettyTable
import awkward as ak

warnings.simplefilter("ignore", UserWarning)
from statsmodels.nonparametric.smoothers_lowess import lowess

import analysis_suite.commons.user as user
from analysis_suite.commons.constants import all_eras
from analysis_suite.commons.configs import get_ntuple, get_ntuple_info, get_inputs
from analysis_suite.combine.card_maker import Card_Maker
from analysis_suite.combine.hist_writer import HistWriter
from analysis_suite.combine.combine_wrapper import runCombine
from analysis_suite.plotting.hist_getter import GraphInfo, HistGetter
from analysis_suite.data.systs import systematics
from analysis_suite.commons.histogram import Histogram
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.plotting.LogFile import LogFile

from analysis_suite.commons.constants import lumi

region_mask = {
     "Dilepton": lambda vg: vg["NMuons"] + vg["NElectrons"] == 2,
     "Multi": lambda vg: vg["NMuons"] + vg["NElectrons"] > 2
}

def write_table(hists, ginfo, outfile):
    def get_str(hist):
        val, err = hist.get_int_err(sqrt_err=True)
        return rf'${val:0.2f}\pm{err:0.2f}$'

    table = PrettyTable(["Process", 'Prefit'])
    # Signal
    sig_name = 'ttt_nlo'
    table.add_row([ginfo.get_legend_name(sig_name), get_str(hists[sig_name])])
    table.add_divider()
    # Background
    total = np.array([0., 0.])
    for group, hist in hists.items():
        if group in [sig_name, 'data']:
            continue
        total += hist.get_int_err(sqrt_err=True)
        table.add_row([ginfo.get_legend_name(group), get_str(hist)])
    table.add_divider()
    # Total
    row = 'total_background'
    table.add_row(['Total Bkg.', rf'${total[0]:0.2f}\pm{total[1]:0.2f}$'])
    table.add_divider()
    # Data
    data = int(hists['data'].sum(flow=True).value)
    table.add_row(["Data", f'{data}'])
    # Signal
    with open(outfile, 'w') as output:
        output.write(table.get_latex_string()+'\n')

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

def fix_jec(nom_hists, syst_hists):
    for syst in systematics:
        if "Jet_JE" not in syst.name:
            continue
        unsmooth, smooth = syst.get_name(year, with_lowess=False), syst.get_name(year)
        if f'{unsmooth}_up' not in syst_hists:
            continue
        up_unsmooth = syst_hists.pop(f'{unsmooth}_up')
        down_unsmooth = syst_hists.pop(f'{unsmooth}_down')
        syst_hists[smooth+'_up'], syst_hists[smooth+'_down'] = {}, {}
        for group in up_unsmooth.keys():
            up, down = smooth_hist(nom_hists[group], up_unsmooth[group], down_unsmooth[group], symm=True)
            syst_hists[smooth+'_up'][group] = up
            syst_hists[smooth+'_down'][group] = down
    return syst_hists

def split_hists(hists, systname, year, chan):
    if systname != 'Nominal':
        updown_break = systname.rfind('_')
        rawsyst, updown = systname[:updown_break], systname[updown_break+1:]
    output = {}
    for group, hist in hists.items():
        if systname == "Nominal":
            yield systname, group, hist
            continue
        for syst in systematics:
            if syst.name == rawsyst and syst.good_syst(group, year, chan):
                out_name = syst.get_name(year, with_lowess=False)
                yield f'{out_name}_{updown}', group, hist

# Read in files/hists
def read_histograms(workdir, ntuple, year, graphs, infile, region):
    syst_hists = {name: dict() for name in graphs}
    ntuple = get_ntuple(*ntuple)
    ntuple.remove_group('nonprompt_mc')

    mask = region_mask.get(region, None)
    hist_factory = HistGetter(ntuple, year, filename=infile, workdir=workdir,
                              scales=['btag_jetlep', 'wz'],
                              mask = mask,
                              )
    nom_vals = {}
    group_nom = {}
    for group, member, df in hist_factory.df_iter():
        if member == 'data':
            nom_vals[group] = np.sum(df.scale)
        else:
            nom_vals[member] = np.sum(df.scale)
        if group not in group_nom:
            group_nom[group] = 0.
        group_nom[group] += np.sum(df.scale)

    for syst in hist_factory.systs:
        print(f"Processing: {syst}")
        hist_factory.reset_syst(syst)
        hists = hist_factory.get_hists(graphs, fix_negative=True)
        for graph_name, subhists in hists.items():
            for systname, group, hist in split_hists(subhists, syst, year, region):
                if systname not in syst_hists[graph_name]:
                    syst_hists[graph_name][systname] = {}
                syst_hists[graph_name][systname][group] = hist
    return syst_hists

def make_hist(ntuple, workdir, outdir, year, region, graphs, useFlat=False, cores=1):
    outfiles = []
    all_hists = {name: dict() for name in graphs}
    if useFlat:
        infiles = (workdir/year).glob(f'processed_*_signal.root')
    else:
        infiles = [None]

    if cores == 1 or not useFlat:
        outputs = [read_histograms(workdir, ntuple, year, graphs, f, region) for f in infiles]
    else:
        inputs = []
        for infile in infiles:
            inputs.append((workdir, ntuple, year, graphs, infile, region))
        with mp.Pool(cores) as pool:
            outputs = pool.starmap(read_histograms, inputs)
    for output in outputs:
        for name, graph in output.items():
            all_hists[name].update(graph)

    ginfo = get_ntuple_info(*ntuple)
    for graph_name, graph in graphs.items():
        nom_hist = all_hists[graph_name].pop('Nominal')
        syst_hists = fix_jec(nom_hist, all_hists[graph_name])

        outfile = outdir / f'{graph_name}_{year}_{region}.root'
        outfiles.append(outfile)
        with HistWriter(outfile) as writer:
            writer.add_syst(nom_hist, ginfo, syst="Nominal", blind=False)
            for syst, hists in syst_hists.items():
                writer.add_syst(hists, ginfo, syst=syst, blind=False)
    print("Finished Hist making")
    return outfiles


def make_card(workdir, rootfile):
    combine_info = get_inputs(workdir, 'combine_info')
    ginfo = get_ntuple_info('signal')
    outdir = rootfile.parent
    region, year, graph_name = [word[::-1] for word in rootfile.stem[::-1].split('_', 2)]
    print(graph_name, year, region, "card")

    group_list = list()
    used_syst_names = list()
    with uproot.open(rootfile) as f:
        for d, cls in f.classnames(cycle=False, recursive=False).items():
            if 'TH1' in cls and d not in ['data_obs', 'SIG']:
                group_list.append(d)
            elif 'Up' in d:
                used_syst_names.append(d[:-2])

    def keep_systs(x):
        return x.syst_type != 'shape' or x.get_name(year) in used_syst_names

    with Card_Maker(outdir, year, region, ["SIG"], group_list, graph_name) as card:
        card.write_preamble()
        card.write_systematics(list(filter(keep_systs, systematics)), ginfo)
        card.add_stats()
        for rp in combine_info.rate_params:
            card.add_rateParam(ginfo.get_combine_name(rp))

    return outdir/(rootfile.stem+"_card.txt")

def plot_prefit(ntuple, workdir, datacard, bins, axis_name):
    outdir = datacard.parent
    ginfo = ntuple.get_info()
    _, region, year, graph_name = [word[::-1] for word in datacard.stem[::-1].split('_', 3)]
    nbins = bins.size

    command = 'combine -M FitDiagnostics --skipSBFit --skipBOnlyFit --saveShapes --saveWithUncertainties'
    fit_file = outdir/f'fitDiagnostics.{graph_name}.{year}.root'
    # if not fit_file.exists():
    runCombine(f"{command} -d {datacard} -n .{graph_name}.{year}", output=False)

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
        logfile = workdir/f'{graph_name}_{region}_{year}.log'
        write_table(hists, ginfo, logfile)
        plotter.set_extra_text(f'{region}')
        plotter.plot_hist(graph_name, hists, syst_err=syst_err)

def starplot(workdir, rootfile, ntuple_name, write_dir, bins, axis_name):
    ntuple = get_ntuple(*ntuple_name)
    ntuple.remove_group("nonprompt_mc")
    datacard = make_card(workdir, rootfile)
    plot_prefit(ntuple, write_dir, datacard, bins, axis_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", type=lambda x : x.split(','), help="Year to use")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area / x,
                        help="Working Directory")
    parser.add_argument("-n", '--ntuple', default="signal")
    parser.add_argument("-f", '--flat', action='store_true')
    parser.add_argument('-a', '--all', action='store_true')
    parser.add_argument('-r', '--region', required=True)
    parser.add_argument('-t', '--extra', default="")
    parser.add_argument('-ne', '--ntuple_extra', default="info")
    parser.add_argument("-j", '--cores', default=1, type=int)
    args = parser.parse_args()
    years = args.years if args.years is not None else []

    graphs = {
        'njet': GraphInfo(r"$N_j$", axis.Regular(11, 0, 11), 'NJets'),
        'nbjet': GraphInfo(r"$N_b$", axis.Regular(6, 0, 6), "NmediumBJets"),
        "ht": GraphInfo(r"$H_T$", axis.Regular(25, 0, 1200), "HT"),
        'met': GraphInfo(r"Met", axis.Regular(25, 0, 400), 'Met'),
        'l1Pt': GraphInfo('$p_T(l_{{1}})$ (GeV)', axis.Regular(25, 0, 350), 'l1Pt'),
        'l2Pt': GraphInfo('$p_T(l_{{2}})$ (GeV)', axis.Regular(25, 0, 150), 'l2Pt'),
    }
    # graphs = {
    #     "ht": GraphInfo(r"$H_T$", axis.Regular(25, 0, 1000), lambda vg : vg.get_hist('HT')),
    #     'met': GraphInfo(r"Met", axis.Regular(35, 0, 350), lambda vg : vg.get_hist('Met')),
    #     'njet': GraphInfo(r"$N_j$", axis.Regular(9, 0, 9), lambda vg : (vg.Jets.num(), vg.scale)),
    #     'nbjet': GraphInfo(r"$N_b$", axis.Regular(5, 0, 5), lambda vg: vg.get_hist('NBjets_medium')),
    #     'pt_1': GraphInfo('pt', axis.Regular(25, 0, 500), lambda vg: vg.TightLepton.get_hist('pt', 0)),
    #     'pt_2': GraphInfo('pt', axis.Regular(25, 0, 500), lambda vg: vg.TightLepton.get_hist('pt', 1)),

    #     # 'nmuon': GraphInfo('$N_\mu$', axis.Regular(4, 0, 4), lambda vg: (vg.TightMuon.num(), vg.scale)),
    #     # 'elec_pt': GraphInfo('$p_T(e)$ (GeV)', axis.Regular(25, 0, 200), lambda vg: vg.TightElectron.get_hist('pt', -1)),
    #     # 'muon_pt': GraphInfo('$p_T(\mu)$ (GeV)', axis.Regular(25, 0, 200), lambda vg: vg.TightMuon.get_hist('pt', -1)),
    #     # 'njet': GraphInfo(r"$N_j$", axis.Regular(9, 0, 9), lambda vg : (vg.Jets.num(), vg.scale)),

    #     # 'eta_j': GraphInfo('eta', axis.Regular(25, -2.4, 2.4), lambda vg: vg.Jets.get_hist('eta', -1)),
    #     # 'pt_j': GraphInfo('pt', axis.Regular(25, 0, 500), lambda vg: vg.Jets.get_hist('pt', -1)),
    #     # 'phi_j': GraphInfo('phi', axis.Regular(25, -3.14, 3.14), lambda vg: vg.Jets.get_hist('phi', -1)),
    #     # 'eta_b': GraphInfo('eta', axis.Regular(25, -2.4, 2.4), lambda vg: vg.BJets.get_hist('eta', -1)),
    #     # 'pt_b': GraphInfo('pt', axis.Regular(25, 0, 500), lambda vg: vg.BJets.get_hist('pt', -1)),
    #     # 'phi_b': GraphInfo('phi', axis.Regular(25, -3.14, 3.14), lambda vg: vg.BJets.get_hist('phi', -1)),
    #     # 'eta_l': GraphInfo('eta', axis.Regular(25, -2.4, 2.4), lambda vg: vg.TightLepton.get_hist('eta', -1)),
    #     # 'pt_l': GraphInfo('pt', axis.Regular(25, 0, 500), lambda vg: vg.TightLepton.get_hist('pt', -1)),
    #     # 'phi_l': GraphInfo('phi', axis.Regular(25, -3.14, 3.14), lambda vg: vg.TightLepton.get_hist('phi', -1)),
    # }

    region = args.region

    ntuple = (args.ntuple, args.ntuple_extra)
    write_dir_tmp = args.workdir/args.extra
    for year in years:
        if not year:
            continue
        write_dir = write_dir_tmp/year
        combine_dir = write_dir/"tmp"
        combine_dir.mkdir(exist_ok=True, parents=True)
        runCombine.work_dir = combine_dir
        rootfiles = [combine_dir/f'{name}_{year}_{region}.root' for name in graphs]
        make_hist(ntuple, args.workdir, combine_dir, year, region, graphs, args.flat, args.cores)

        input_files = []
        for rootfile, graph in zip(rootfiles, graphs.values()):
            bins = graph.bins()[0]
            axis_name = graph.axis_name
            input_files.append(
                (args.workdir, rootfile, ntuple, write_dir, bins, axis_name)
            )
        if args.cores == 1:
            for inputs in input_files:
                starplot(*inputs)
        else:
            with mp.Pool(args.cores) as pool:
                pool.starmap(starplot, input_files)

    # Combine all the cards together and run
    if args.all:
        write_dir = write_dir_tmp/'all'
        combine_dir = write_dir/"tmp"
        combine_dir.mkdir(exist_ok=True, parents=True)
        runCombine.work_dir = combine_dir
        ntuple = get_ntuple(*ntuple)
        ntuple.remove_group("nonprompt_mc")
        for name, graph in graphs.items():
            print(name)
            func = f'combineCards.py'
            for year in all_eras:
                card = write_dir_tmp/year/'tmp'/f'{name}_{year}_{region}_card.txt'
                func += f' yr{year}={card}'
            datacard = combine_dir/f'{name}_all_{region}_card.txt'
            runCombine(func+f' > {datacard}')
            plot_prefit(ntuple, write_dir, datacard, graph.bins()[0], graph.axis_name)


