#!/usr/bin/env python3
import argparse
import uproot
import numpy as np
from pathlib import Path
from prettytable import PrettyTable
import boost_histogram.axis as axis

import warnings
warnings.simplefilter("ignore", UserWarning)

import analysis_suite.commons.user as user
from analysis_suite.commons.configs import get_ntuple, get_inputs
from analysis_suite.combine.combine_wrapper import runCombine
from analysis_suite.commons.histogram import Histogram
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.plotting.LogFile import LogFile
from analysis_suite.plotting.hist_getter import GraphInfo

from analysis_suite.commons.constants import all_eras

bkg_names = {
    "SIG": 'ttt',
    'NONPROMPT': 'Nonprompt',
    'TTW': r'$\ttbar+W$',
    'TTZ': r'$\ttbar+Z$',
    'TTH': r'$\ttbar+H$',
    'TTTT': r'$\ttbar\ttbar$',
    'RARE': 'Rare',
    'XG': r'$X+\gamma$',
    'CHARGE_FLIP': "Charge Flip",
    "total_background": "Total Bkg.",
}

extra_graphs = {
    'hot0_muon': GraphInfo("DNN", axis.Regular(25, 0.48, 0.98), ""),
    'hot1_muon':GraphInfo("DNN", axis.Regular(25, 0.48, 0.98), ""),
    'hot0_elec':GraphInfo("DNN", axis.Regular(25, 0.48, 0.98), ""),
    'hot1_elec':GraphInfo("DNN", axis.Regular(25, 0.48, 0.98), ""),
}


def write_table(infile, regions, outfile):
    prefit_data = {}
    postfit_data = {}
    data = 0
    for region in regions:
        prefit = infile['shapes_prefit'][region]
        postfit = infile['shapes_fit_s'][region]
        for fname, outname in bkg_names.items():
            if fname not in prefit or sum(prefit[fname].values()) < 1e-3:
                continue
            if outname not in prefit_data:
                prefit_data[outname] = np.zeros(2)
                postfit_data[outname] = np.zeros(2)
            prefit_data[outname] += np.array([sum(prefit[fname].values()), sum(prefit[fname].variances())])
            postfit_data[outname] += np.array([sum(postfit[fname].values()), sum(postfit[fname].variances())])
        data += int(sum(prefit['data'].values()[1]))


    table = PrettyTable(["Process", 'Prefit', 'Postfit', 'Post/Prefit'])
    for key, pre_val in prefit_data.items():
        post_val = postfit_data[key]
        def get_str(hist):
            return rf'${hist[0]:0.2f}\pm{np.sqrt(hist[1]):0.2f}$'
        table.add_row([key, get_str(pre_val), get_str(post_val), f'${post_val[0]/pre_val[0]:0.2f}$'])
    table.add_divider()
    table.add_row(["Data", data, data, ""])

    with open(outfile, 'w') as f:
        f.write(table.get_latex_string()+'\n')



def plot(workdir, year, graph_info, extra='wisc', prefit=True, combined=False):
    if prefit:
        tree = "shapes_prefit"
        name = 'prefit'
    else:
        tree = "shapes_fit_s"
        name = 'postfit'
    first_year = year

    command = 'combineTool.py -M FitDiagnostics --saveToys --rMax 150 --rMin -150 --saveShapes --saveWithUncertainties --cminDefaultMinimizerStrategy 0 -v 1'
    # May need to change based on needs
    # datacard = f'final_{year}_card.root'
    datacard = f"3top_Run2_card.root"
    fitFile = f'fitDiagnostics.{year}.{extra}.root'
    if not (workdir/fitFile).exists():
        runCombine(f"{command} -d {datacard} -n .{year}.{extra}", output=True)
    print('done')

    with uproot.open(workdir/fitFile) as f:
        hists = {}
        syst_err = {}
        if not prefit and "final" in datacard:
            regs = list(f[tree].keys(recursive=False, cycle=False))
            sub_regs = np.unique([r.split("_")[-1] for r in regs])
            for reg in sub_regs:
                logfile = workdir/f'yield_{reg}_{year}.log'
                write_table(f, [r for r in regs if reg in r], logfile)
        for reg, r_dir in f[tree].items(recursive=False, cycle=False):
            if 'brown' in reg or "TTTX" in reg:
                year_dict = {'16': '2016post', '16APV': '2016pre', '17': '2017', '18': '2018'}
                year = year_dict[reg.split('_')[-1]]
                reg_name = reg.split('_')[-3]
                region = f"hot{reg_name[9]}_"
                region += "elec" if reg_name[2] == "E" else "muon"
                ntuple = get_ntuple('brown')
            else:
                reg_split = reg.split("_")
                if len(reg_split) == 1:
                    region = reg
                else:
                    year = reg_split[-2].replace("yr", "")
                    region = reg_split[-1]
                ntuple = get_ntuple('signal')
                if 'era' in year:
                    year = year[3:]

            if combined:
                year = 'all'
            if year not in hists:
                hists[year] = {}
                syst_err[year] = {}
            if region not in hists[year]:
                hists[year][region] = {}
                syst_err[year][region] = None
            ginfo = ntuple.get_info()

            graph = graph_info[region]
            bins = graph.bins()[0]
            nbins = bins.size
            axis_name = graph.axis_name

            # print(region, year)
            def get_hist(group, hist):
                if combine_group == 'total_background':
                    return hist.variances()[:nbins]
                elif "total" in combine_group:
                    return None
                elif "TGraph" in repr(hist):
                    vals, sumw2 = hist.values()[1], hist.errors('mean')[1]
                    sumw2 = sumw2**2
                else:
                    vals, sumw2 = hist.values(), hist.variances()
                tmp_hist = Histogram(bins, axis_name=axis_name)
                tmp_hist.set_data(val=vals[:nbins], sumw2=sumw2[:nbins])
                return tmp_hist

            for combine_group, h in r_dir.items(cycle=False):

                hist = get_hist(combine_group, h)
                group = ginfo.get_group_from_combine(combine_group)
                if hist is None:
                    continue
                elif isinstance(hist, Histogram):
                    if group in hists[year][region]:
                        hists[year][region][group] += hist
                    else:
                        hists[year][region][group] = hist
                else:
                    if syst_err[year][region] is None:
                        syst_err[year][region] = hist
                    else:
                        syst_err[year][region] += hist

        for year, year_dict in hists.items():
            for region, plot_hists in year_dict.items():
                if "hot" in region:
                    ntuple = get_ntuple('brown')
                else:
                    ntuple = get_ntuple('signal')
                sig = None
                if prefit:
                    sig = 'ttt_nlo'
                plotter = Plotter(ntuple, year, sig=sig, bkg='all',
                                outdir=workdir/'comb_plots', from_ntuple=False)
                plotter.set_year(year)
                plotter.set_extra_text(region)
                text = False
                log = False
                rtop, rbot = 1.5, 0.5
                if 'hot' in region:
                    log = True
                    lep = 'e' if 'elec' in region else r'\mu'
                    text = fr'${lep}, N_{{HOT}}={region[3]}$'
                    rtop, rbot = 1.3, 0.7
                else:
                    reg_dict = {
                        'Dilepton': r"$N_{\ell}=2$", "Multi": r"$N_{\ell}\geq3$",
                        "ttzCR": "ttZ CR", "ttttCR": "tttt CR", 'ttz3Valid': "ttZ CR",
                        'ttz4Valid': "ttZ CR"
                    }
                    text = reg_dict[region]

                print(name, year, region)
                plotter.plot_hist(name, plot_hists, syst_err=syst_err[year][region],
                                  log=log, text=text, ratio_top=rtop, ratio_bot=rbot)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True, type=lambda x: x.split(','),
                        help="Year to use")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area / x,
                        help="Working Directory")
    parser.add_argument("-t", '--extra_text', default="")
    parser.add_argument("--pre", action='store_true')
    parser.add_argument("-a", '--combine', action='store_true')
    parser.add_argument("-e", '--extra', help='Extra name in the output of the files', default='wisc')
    args = parser.parse_args()

    graph_info = {n: v['graph'] for n, v in get_inputs(args.workdir, 'combine_info').regions.items()}
    graph_info.update(extra_graphs)
    combine_dir = Path('.')

    for year in args.years:
        if args.pre:
            plot(combine_dir, year, graph_info, args.extra, args.pre, combined=args.combine)
        plot(combine_dir, year, graph_info, args.extra, False, combined=args.combine)
