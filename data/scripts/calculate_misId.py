#!/usr/bin/env python3
import boost_histogram as bh
import argparse
from scipy.optimize import minimize
import numpy as np
import awkward as ak;
import pickle
from pathlib import Path
from matplotlib import colors

import analysis_suite.commons.configs as config
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import hep, plot, plot_colorbar
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.info import GroupInfo
import analysis_suite.data.plotInfo.misId_fakerate as pinfo
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.commons.user import workspace_area
from datetime import datetime

latex_chan = {"Electron": "e", "Muon": "\mu",
              "EE": 'ee', "EM": 'e\mu', 'MM': '\mu\mu'}

def plot_project(filename, tight, loose, axis_label, lumi, axis=None):
    with plot(filename) as ax:
        if axis is None:
            eff_hist = Histogram.efficiency(tight, loose)
        else:
            eff_hist = Histogram.efficiency(tight.project(axis), loose.project(axis))
        eff_hist.plot_points(ax)
        ax.set_xlabel(axis_label)
        hep.cms.label(ax=ax, lumi=lumi)

def fr_plot(name, hist, year, **kwargs):
    with plot(name, figsize=(18, 10)) as ax:
        xx = np.tile(hist.axes[0].edges, (len(hist.axes[1])+1, 1))
        yy = np.tile(hist.axes[1].edges, (len(hist.axes[0])+1, 1)).T
        mesh = ax.pcolormesh(xx, yy, hist.vals.T, norm=colors.LogNorm(vmin=1e-5, vmax=1e-2),
                             shading='flat', cmap="Blues", **kwargs)
        for j, y in enumerate(hist.axes[1].centers):
            for i, x in enumerate(hist.axes[0].centers):
                exp = int(np.floor(np.log10(hist.vals[i,j])))
                val_str = f'{hist.vals[i,j]:.2e}\n$\pm${hist.err[i,j]/10**exp:.2f}e{str(exp).zfill(3)}'
                text = ax.text(x, y, val_str, fontsize='xx-small', ha="center", va='center')
        ax.set_title("Electron Charge MisId Rate")
        ax.set_xlabel("$p_{T}(e)$ [GeV]")
        ax.set_ylabel("$|\eta(e)|$")
        hep.cms.label(ax=ax, year=year)
        plot_colorbar(mesh, ax, barpercent=2)

def get_fake_rate(part, fake_rate, idx, name="Muon"):
    pt_axis, eta_axis = fake_rate.axes
    npt, neta = fake_rate.axes.size

    ptbin = np.digitize(part['pt', idx], pt_axis.edges) - 1
    ptbin = np.where(ptbin >= npt, npt-1, ptbin)
    etabin = np.digitize(part.abseta(idx), eta_axis.edges) - 1
    return fake_rate.values().flatten()[etabin + neta*ptbin]

def scale_misId(vg, fake_rate):
    part = vg.TightElectron
    fr1 = get_fake_rate(part, fake_rate, 0)
    fr2 = get_fake_rate(part, fake_rate, 1)
    print(sum((fr1/(1-fr1)+fr2/(1-fr2))))
    vg.scale = (fr1/(1-fr1)+fr2/(1-fr2))*0.86

def measurement(workdir, ginfo, year, input_dir):
    plot_dir = workdir / f'MR_{year}'
    plot_dir.mkdir(exist_ok=True)

    bkg = ["DY_ht", "ttbar_lep", 'VV']
    groups = ginfo.setup_groups(["data"] + bkg)
    ntuple = config.get_ntuple('charge_misId', 'measurement')
    chans = ['MM', 'EE', 'EM']

    graphs = pinfo.charge_misId['Measurement']
    graphs_1d = [graph for graph in graphs if graph.dim() == 1]

    plotter = Plotter(ntuple.get_file(year=year, workdir=input_dir), groups, ntuple=ntuple, year=year)
    plotter.cut(lambda vg : vg["Met"] > 25)
    plotter.set_groups(bkg=bkg)
    # for chan in chans:
    #     latex = latex_chan[chan]
    #     plotter.mask(lambda vg : vg["TightElectron"].num() == chan.count('E'))
    #     plotter.fill_hists(graphs, ginfo)
    #     for graph in graphs_1d:
    #         plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_{chan}.png', chan=latex, region=f"$MR({latex})$")

    plotter.mask(lambda vg : vg["TightElectron"].num() > 0)
    plotter.fill_hists(graphs, ginfo)
    # for graph in graphs_1d:
    #     plotter.plot_stack(graph.name, plot_dir/f'{graph.name}_e.png', chan="e\ell", region="$MR(e\ell)$")

    all_fr = plotter.get_sum(bkg, 'all_fr')
    flip_fr = plotter.get_sum(bkg, 'flip_fr')
    all_pt = plotter.get_sum(bkg, 'pt_allq')
    flip_pt = plotter.get_sum(bkg, 'pt_flipq')

    fr = Histogram.efficiency(flip_fr, all_fr)
    fr_pt = Histogram.efficiency(flip_pt, all_pt)
    print(fr.vals)
    print(fr_pt.vals)
    fr_plot(plot_dir/f'fr_{year}', fr, year)
    plot_project(plot_dir/f'fr_pt_tight.png', flip_pt, all_pt, "$p_{{T}}(e)$", lumi[year])

    plot_project(plot_dir/f'fr_pt.png', flip_fr, all_fr, "$p_{{T}}(e)$", lumi[year], axis=0)
    plot_project(plot_dir/f'fr_eta.png', flip_fr, all_fr, '$\eta(e)$', lumi[year], axis=1)

    # Dump fake rate
    with open(workdir/f"charge_misid_rate_{year}.pickle", "wb") as f:
        pickle.dump(fr.hist, f)


def closure(workdir, ginfo, year, input_dir):
    plot_dir = workdir / f'CR_{year}'
    plot_dir.mkdir(exist_ok=True)

    mc_bkg = [# "DY_ht",
              "ttbar_lep", 'VV', 'DY']
    groups = ginfo.setup_groups([ "data", "charge_flip"]+ mc_bkg)
    graphs = pinfo.charge_misId['Closure']

    # # TF setup
    # ntuple_os = config.get_ntuple('charge_misId', 'closure_os')
    # plotter_os = Plotter(ntuple_os.get_file(year=year, workdir=input_dir), groups, ntuple=ntuple_os, year=year)
    # # plotter_os.cut(lambda vg : vg["LHE_HT"] > 70, "DYm50_amc")
    # plotter_os.set_groups(bkg=["DY", "ttbar_lep", 'VV'])
    # plotter_os.fill_hists(graphs, ginfo)
    # for graph in graphs:
    #     plotter_os.plot_stack(graph.name, plot_dir/f'{graph.name}_OS_allMet.png', chan='ee', region="$OS({})$")

    # plotter_os.mask(lambda vg : vg["Met"] < 50, clear=False)
    # plotter_os.fill_hists(graphs, ginfo)
    # for graph in graphs:
    #     plotter_os.plot_stack(graph.name, plot_dir/f'{graph.name}_OS.png', chan='ee', region="$OS({})$")
    

    ntuple_ss = config.get_ntuple('charge_misId', 'closure_ss')
    plotter_ss = Plotter(ntuple_ss.get_file(year=year, workdir=input_dir), groups, ntuple=ntuple_ss, year=year)
    plotter_ss.cut(lambda vg : vg["Met"] < 50)
    # plotter_ss.cut(lambda vg : vg["LHE_HT"] < 70, "DYm50_amc")
    plotter_ss.set_groups(bkg=["DY", "ttbar_lep", 'VV'])
    # plotter_ss.fill_hists(graphs, ginfo)
    # for graph in graphs:
    #     plotter_ss.plot_stack(graph.name, plot_dir/f'{graph.name}_SS_allMet_mc.png', chan='ee', region="$SS({})$")

    # plotter_ss.fill_hists(graphs, ginfo)
    # for graph in graphs:
    #     plotter_ss.plot_stack(graph.name, plot_dir/f'{graph.name}_SS_mc.png', chan='ee', region="$SS({})$")
    # exit()

    with open(workdir/f"charge_misid_rate_{year}.pickle", "rb") as f:
        fake_rates = pickle.load(f)
        print(fake_rates)
    plotter_ss.scale(lambda vg : scale_misId(vg, fake_rates), groups='charge_flip')
    print(plotter_ss.get_integral())
    plotter_ss.set_groups(bkg=mc_bkg, sig='charge_flip')
    plotter_ss.fill_hists(graphs, ginfo)
    # for graph in graphs:
    #     plotter_ss.plot_stack(graph.name, plot_dir/f'{graph.name}_SS_all.png', chan='ee', region="$SS({})$")
    # Plot with just Data-Driven Background
    plotter_ss.set_groups(bkg=['charge_flip'])
    for graph in graphs:
        plotter_ss.plot_stack(graph.name, plot_dir/f'{graph.name}_SS_data.png', chan='ee', region="$SS({})$")


if __name__ == "__main__":
    workdir = workspace_area/'charge_misId'

    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-d', '--workdir', help="directory to run over. If nothing, use date",)
    parser.add_argument('-i', '--input_dir', required=True)
    parser.add_argument('-r', '--run', type=lambda x: [i.strip() for i in x.split(',')],
                        help="Regions to run through (sideband, measurement, closure)")
    args = parser.parse_args()

    workdir /= datetime.now().strftime("%m%d") if args.workdir is None else args.workdir
    workdir.mkdir(exist_ok=True)

    color_by_group = {
        "data": "black",
        "DY_ht": "goldenrod",
        "DY": "goldenrod",
        "ttbar_lep": 'royalblue',
        "VV": 'mediumorchid',
        'charge_flip': 'seagreen',
    }
    ginfo = GroupInfo(color_by_group)
    for year in args.years:
        if 'measurement' in args.run:
            print("Measurement")
            measurement(workdir, ginfo, year, args.input_dir)
        if 'closure' in args.run:
            print("Closure")
            closure(workdir, ginfo, year, args.input_dir)
