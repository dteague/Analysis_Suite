#!/usr/bin/env python3
import argparse
import numpy as np
import pickle
import awkward as ak
import boost_histogram.axis as axis
from datetime import datetime
import operator
from copy import copy
from matplotlib.colors import LogNorm

from analysis_suite.plotting.plotter import GraphInfo
import analysis_suite.commons.configs as config
from analysis_suite.commons.plot_utils import hep, plot, plot_colorbar
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.info import GroupInfo
import analysis_suite.data.plotInfo.nonprompt_fakerate as pinfo
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.commons.user import workspace_area

latex_chan = {"Electron": "e", "Muon": "\mu",
              "EE": 'ee', "EM": 'e\mu', 'MM': '\mu\mu'}

def fr_plot(name, ratio, lumi):
    vals = ratio.vals
    with plot(name) as ax:
        xx, yy = np.meshgrid(ratio.axes[0].edges, ratio.axes[1].edges)
        mesh = ax.pcolormesh(xx, yy, vals.T, norm=LogNorm(), shading='auto')
        xx, yy = np.meshgrid(ratio.axes[0].centers, ratio.axes[1].centers)
        ax.contour(xx, yy, vals.T, np.logspace(2, int(np.log10(np.max(vals))), 3))
        ax.plot([-1, 1], [-1, 1], color='red')
        ax.set_xlabel("$Disc_{old}$")
        ax.set_ylabel("$Disc_{new}$")
        plot_colorbar(mesh, ax)
        hep.cms.label(ax=ax, lumi=lumi)

def compare_mva(name, old, new, lumi):
    old.name = "$disc_{old}$"
    old.color = 'r'
    new.name = "$disc_{new}$"
    new.color = 'b'
    with plot(name) as ax:
        old.plot_shape(ax)
        new.plot_shape(ax)
        ax.set_xlabel("$Disc_{tth}$")
        ax.set_ylabel("Normalized Events (A.U)")
        ax.legend()
        hep.cms.label(ax=ax, lumi=lumi)

def test_mva(workdir, year, input_dir):
    plot_dir = workdir / f'mva_new'
    plot_dir.mkdir(exist_ok=True)

    chan = "Electron"

    mc = ['qcd',  'ewk']

    ntuple = config.get_ntuple('fake_rate', 'measurement')
    groups = ntuple.get_info().setup_groups(mc)
    filename = ntuple.get_filename(year=year, workdir=input_dir)

    plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
    plotter.set_groups(sig='ewk', bkg=['qcd'])
    plotter.cut(lambda vg : vg[chan].num() == 1)
    graphs = [
        GraphInfo(f"mva", '', (axis.Regular(24, -1, 1), axis.Regular(24, -1, 1)), lambda vg : vg[chan].get_hist2d('mvaTTH', 'new_mvaTTH', 0)),
        GraphInfo(f"mva_old", '$Disc_{{old}}$', axis.Regular(40, -1, 1), lambda vg : vg[chan].get_hist('mvaTTH', 0)),
        GraphInfo(f"mva_new", '$Disc_{{new}}$', axis.Regular(40, -1, 1), lambda vg : vg[chan].get_hist('new_mvaTTH', 0)),
    ]
    plotter.fill_hists(graphs)
    plotter.sig_colors[0] = "EWK"
    plotter.bkg_colors[0] = "QCD"

    def print_info(dic):
        for name, hist in dic.items():
            print(name)
            edges = hist.axis.edges[:-1]
            print(edges[edges > 0.35])
            print((np.cumsum(hist.vals[::-1])[::-1]/np.sum(hist.vals))[edges > 0.35])
            print()

    # print("old")
    # print_info(mva_old)
    # print("new")
    # print_info(mva_new)
    # exit()
    for graph in graphs:
        if graph.dim() != 1:
            continue
        plotter.plot_shape(graph.name, plot_dir/f'{graph.name}_{chan}_{year}.png', chan='e')


    mva = plotter.get_hists('mva')
    mva_old = plotter.get_hists('mva_old')
    mva_new = plotter.get_hists('mva_new')

    compare_mva(plot_dir/f'mva_{chan}_compare_qcd_{year}.png', mva_old['qcd'], mva_new['qcd'], lumi[year])
    compare_mva(plot_dir/f'mva_{chan}_compare_ewk_{year}.png', mva_old['ewk'], mva_new['ewk'], lumi[year])
    exit()
    fr_plot(plot_dir/f"mva_{chan}_qcd_{year}.png", mva['qcd'], lumi[year])
    fr_plot(plot_dir/f"mva_{chan}_ewk_{year}.png", mva['ewk'], lumi[year])

if __name__ == "__main__":
    workdir = workspace_area

    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-d', '--workdir', default=datetime.now().strftime("%m%d"),
                        help="directory to run over. If nothing, use date")
    parser.add_argument('-i', '--input_dir', default=None)
    args = parser.parse_args()

    workdir = workspace_area/args.workdir
    workdir.mkdir(exist_ok=True)

    for year in args.years:
        test_mva(workdir, year, args.input_dir)
