#!/usr/bin/env python3
import argparse
import numpy as np
import boost_histogram.axis as axis
import pickle
import warnings
from copy import copy
warnings.simplefilter("ignore", UserWarning)

from analysis_suite.plotting.plotter import Plotter
from analysis_suite.plotting.hist_getter import GraphInfo, HistGetter
import analysis_suite.commons.configs as config
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import plot, plot_colorbar, cms_label, ratio_plot
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.info import GroupInfo
from analysis_suite.commons.user import workspace_area

def scale_wz(vg, fakerate):
    njets = vg.Jets.num()
    jet_bin = np.digitize(njets, np.arange(2, 6)) - 1
    vg.scale = fakerate[jet_bin]

def measurement(workdir, year, input_dir):
    chans = ['Electron', 'Muon']
    ntuple = config.get_ntuple('WZ')
    plotter = Plotter(ntuple, year, bkg='all', outdir=workdir/f'WZ_CR/MR_{year}')

    njet_bins = axis.Regular(4, 2, 6)
    graphs = {
        'scale':   GraphInfo('$N_j$', njet_bins, lambda vg: (vg['Jets'].num(), vg.scale)),

        'njets':   GraphInfo('$N_j$', axis.Regular(7, 0, 7), lambda vg: (vg['Jets'].num(), vg.scale)),
        # 'nbjets':  GraphInfo('$N_b$', axis.Regular(5, 0, 5), lambda vg: vg.get_hist('NBjets_medium')),
        'zmass':   GraphInfo('$M_Z$', axis.Regular(20, 75, 105), lambda vg: vg.get_hist("Zmass")),
        'mt':      GraphInfo('$M_T$', axis.Regular(20, 0, 400), lambda vg: vg.TightLepton.get_hist("mt", -1)),
        'met':     GraphInfo('$p_T^{{miss}}$ (GeV)', axis.Regular(20, 0, 400), lambda vg: vg.get_hist("Met")),
        'l1Pt':    GraphInfo('$p_T(l_1)$ (GeV)', axis.Regular(20, 0, 500), lambda vg: vg['TightLepton'].get_hist('pt', 0)),
        'l2Pt':    GraphInfo('$p_T(l_2)$ (GeV)', axis.Regular(20, 0, 300), lambda vg: vg['TightLepton'].get_hist('pt', 1)),
        'l3Pt':    GraphInfo('$p_T(l_3)$ (GeV)', axis.Regular(20, 0, 200), lambda vg: vg['TightLepton'].get_hist('pt', 2)),
        'nmuon':   GraphInfo(r'$N_{{\mu}}$', axis.Regular(4, 0, 4), lambda vg: (vg["TightMuon"].num(), vg.scale)),
    }


    hist_factory = HistGetter(ntuple, year, workdir=workdir, scales=['btag_jetbinned'])
    hists = hist_factory.get_hists(graphs)
    print(hists['mt']['wz'].vals)
    plotter.plot_hists(hists)

    njet_hist = hists['scale']
    mc = None
    for group, hist in njet_hist.items():
        if group not in ['data', 'wz']:
            if mc is None:
                mc = copy(hist)
            else:
                mc += hist
    scale_hist = (njet_hist['data']-mc)/njet_hist['wz']

    with plot(workdir/'WZ_CR'/f'wz_scale_{year}.png') as ax:
        scale_hist.plot_points(ax)
        i = 0
        for val, err in zip(scale_hist.vals, scale_hist.err):
            print(val, err)
            ax.annotate(f'${val:0.3f}\pm{err:0.3f}$', (i+1, val),
                        xytext=(i+1.1, val*1.1), fontsize=15)
            i += 1
        ax.set_xlabel(r'$N_{j}$')
        cms_label(ax, year=year, hasData=True)

    # Dump MC scale factors
    scale_file = workdir/f"wz_scale_factor.pickle"
    if scale_file.exists():
        with open(scale_file, "rb") as f:
            wz_scales = pickle.load(f)
    else:
        wz_scales = dict()
    wz_scales[year] = scale_hist
    with open(scale_file, 'wb') as f:
        pickle.dump(wz_scales, f)

    hist_factory = HistGetter(ntuple, year, workdir=workdir, scales=['btag_jetbinned', 'wz'])
    hists = hist_factory.get_hists(graphs)
    print(hists['mt']['wz'].vals)
    plotter.set_extra_text("scaled")
    plotter.plot_hists(hists)


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
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : workspace_area / x,
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
