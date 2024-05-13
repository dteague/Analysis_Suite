#!/usr/bin/env python3
import numpy as np
import boost_histogram.axis as axis
import warnings

from analysis_suite.plotting.plotter import Plotter
from analysis_suite.commons.info import NtupleInfo
from analysis_suite.plotting.plotter import GraphInfo
import analysis_suite.commons.user as user
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.plot_utils import plot, nonratio_plot, ratio_plot, hep


ntuple = NtupleInfo(
    filename = user.hdfs_area/'workspace/endcap_study/{year}/{workdir}/',
    trees = ["Signal"],
    color_by_group = {
        "ttt": "crimson",
        'nonprompt_mc': 'gray',
        "xg": "indigo",
        "ttw": "olivedrab",
        "tth": "goldenrod",
        "ttz": "steelblue",
        "ttXY": "teal",
        "rare_nowz": "deeppink",
        "wz": 'slateblue',
        "tttt": "tomato",
    }
)

graphs = [
    GraphInfo('dxy', '$d_{{xy}}$', axis.Regular(40, 0, 0.1), lambda vg: vg["Electrons"].get_hist('dxy', -1)),
    GraphInfo('dz', '$d_{{z}}$', axis.Regular(40, 0, 0.2), lambda vg: vg["Electrons"].get_hist('dz', -1)),
    GraphInfo('lostHits', 'lostHits', axis.Regular(3, 0, 3), lambda vg: vg["Electrons"].get_hist('lostHits', -1)),
    GraphInfo('sip3d', '$\sigma_{IP}$', axis.Regular(40, 0, 10), lambda vg: vg["Electrons"].get_hist('sip3d', -1)),
    GraphInfo('mvaTTH', '$Disc_{TTH}$', axis.Regular(40, -1, 1), lambda vg: vg["Electrons"].get_hist('mvaTTH', -1)),
]

year = '2018'
output = user.workspace_area/'endcap'
groups = ntuple.get_info().setup_groups()
filename = ntuple.get_filename(year=year)
plotter = Plotter(filename, groups, ntuple=ntuple, year=year)
plotter.cut(lambda vg : vg['Electrons'].num()>0)

signal = 'ttt'
plotter.set_groups(sig=signal, bkg=[g for g in groups if g != signal])

def mask(plotter, truth=True, eta_split=None):
    plotter.reset_part("Electrons")
    plotter.mask_part('Electrons', "convVeto", lambda var : var==1)
    plotter.mask_part('Electrons', "dxy", lambda var : var < 0.05)
    plotter.mask_part('Electrons', "dz", lambda var : var < 0.1)
    plotter.mask_part('Electrons', "lostHits", lambda var : var == 0)
    # plotter.mask_part('Electrons', "sip3d", lambda var : var < 4)
    if truth is None:
        pass
    elif truth:
        plotter.mask_part("Electrons", 'truth', lambda var: var == 1)
    else:
        plotter.mask_part("Electrons", 'truth', lambda var: var != 1)
    if eta_split == "EB":
        plotter.mask_part("Electrons", 'eta', lambda var : abs(var) < 1.5)
    elif eta_split == "EE":
        plotter.mask_part("Electrons", 'eta', lambda var : abs(var) > 1.5)

def plot_stack(plotter, name, outfile, sig=None):
    graph = plotter.graphs[name]
    signal = plotter.hists[name].get(sig, None)
    stack = plotter.make_stack(name)
    axis_name = graph.axis_name
    if name in ['mvaTTH']:
        ratio_func = lambda x : np.cumsum(x[::-1])[::-1]/np.sum(x)
    else:
        ratio_func = lambda x : np.cumsum(x)/np.sum(x)

    with ratio_plot(outfile, axis_name, stack.get_xrange(), ratio_top=1.1, subpad_label="Cut Eff.") as ax:
        pad, subpad = ax

        error = Histogram("Stat Errors", graph.bin_tuple, color="plum")
        for hist in stack.stack:
            error += hist

        #upper pad
        if signal:
            rounder = 50
            sig_scale = int(np.max(stack.vals)/np.max(signal.vals)//rounder*rounder)
            signal.scale(sig_scale, changeName=True, forPlot=True)

        stack.plot_stack(pad)
        if signal is not None:
            signal.plot_points(pad)
        error.plot_band(pad)

        # ratio pad
        ratio = Histogram("Ratio", graph.bin_tuple, color="black")
        ratio.hist.values()[:] = ratio_func(error.vals)
        ratio.plot_points(subpad)
        if signal:
            ratio_sig = Histogram("Ratio", graph.bin_tuple, color="red")
            ratio_sig.hist.values()[:] = ratio_func(signal.vals)
            ratio_sig.plot_points(subpad)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hep.cms.label(ax=pad, lumi=plotter.lumi)

def plot_eff(plotter, graph, outfile, eta_split, sig=None):
    name = graph.name
    axis_name = graph.axis_name
    if name in ['mvaTTH']:
        ratio_func = lambda x : np.cumsum(x[::-1])[::-1]/np.sum(x)
    else:
        ratio_func = lambda x : np.cumsum(x)/np.sum(x)
    mask(plotter, truth=False, eta_split=eta_split)
    plotter.fill_hists([graph])
    false_hist = plotter.get_sum(groups.keys(), graph)
    false_hist.set_plot_details(["False Rate", 'red'])
    false_hist.hist.values()[:] = ratio_func(false_hist.vals)
    false_hist.hist.variances()[:] = np.zeros(len(false_hist.vals))

    mask(plotter, truth=True, eta_split=eta_split)
    plotter.fill_hists([graph])
    true_hist = plotter.get_sum(groups.keys(), graph)
    true_hist.set_plot_details(["True Rate", 'black'])
    true_hist.hist.values()[:] = ratio_func(true_hist.vals)
    true_hist.hist.variances()[:] = np.zeros(len(false_hist.vals))

    with ratio_plot(outfile, axis_name, false_hist.get_xrange(), ratio_bot=1, ratio_top=3) as ax:
        pad, subpad = ax

        false_hist.plot_points(pad)
        true_hist.plot_points(pad)

        # ratio pad
        ratio = Histogram("True/False", graph.bin_tuple, color="black")
        ratio += true_hist/false_hist
        ratio.plot_points(subpad)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hep.cms.label(ax=pad, lumi=plotter.lumi)


def plot_pt(plotter, outfile, eta_split, sig=None):

    plotter.reset_part("Electrons")
    plotter.mask_part('Electrons', "convVeto", lambda var : var==1)
    plotter.mask_part('Electrons', "dxy", lambda var : var < 0.05)
    plotter.mask_part('Electrons', "dz", lambda var : var < 0.1)
    plotter.mask_part('Electrons', "lostHits", lambda var : var == 0)
    plotter.mask_part('Electrons', "sip3d", lambda var : var < 4)

    name = graph.name
    axis_name = graph.axis_name
    if name in ['mvaTTH']:
        ratio_func = lambda x : np.cumsum(x[::-1])[::-1]/np.sum(x)
    else:
        ratio_func = lambda x : np.cumsum(x)/np.sum(x)
    mask(plotter, truth=False, eta_split=eta_split)
    plotter.fill_hists([graph])
    false_hist = plotter.get_sum(groups.keys(), graph)
    false_hist.set_plot_details(["False Rate", 'red'])
    false_hist.hist.values()[:] = ratio_func(false_hist.vals)
    false_hist.hist.variances()[:] = np.zeros(len(false_hist.vals))

    mask(plotter, truth=True, eta_split=eta_split)
    plotter.fill_hists([graph])
    true_hist = plotter.get_sum(groups.keys(), graph)
    true_hist.set_plot_details(["True Rate", 'black'])
    true_hist.hist.values()[:] = ratio_func(true_hist.vals)
    true_hist.hist.variances()[:] = np.zeros(len(false_hist.vals))

    with ratio_plot(outfile, axis_name, false_hist.get_xrange(), ratio_bot=1, ratio_top=3) as ax:
        pad, subpad = ax

        false_hist.plot_points(pad)
        true_hist.plot_points(pad)

        # ratio pad
        ratio = Histogram("True/False", graph.bin_tuple, color="black")
        ratio += true_hist/false_hist
        ratio.plot_points(subpad)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hep.cms.label(ax=pad, lumi=plotter.lumi)

nbins = 40
pt_tot = np.zeros(nbins)
wgt_tot = np.zeros(nbins)
for vg in plotter.getters():
    mva = vg.Electrons['mvaTTH', -1]
    pt = vg.Electrons['pt', -1]
    scale = vg.Electrons.scale(-1)
    for i, mva_cut in enumerate(np.linspace(-1, 1, nbins, endpoint=False)):
        mask = mva > mva_cut
        pt_tot[i] += np.sum(pt[mask]*scale[mask])
        wgt_tot[i] += np.sum(scale[mask])

bins = np.linspace(-1, 1, nbins+1)
with nonratio_plot(output/'mean_pt.png', "Cut Applied to $Disc_{TTH}$", bins, pad_label='$mean(p_{T})$') as pad:
    pad.hist(x=bins[:-1], weights=pt_tot/wgt_tot, bins=bins, histtype='step', linewidth=1.5)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hep.cms.label(ax=pad, lumi=plotter.lumi)

# graph_2d = [
#     GraphInfo('mva_2d', '', (axis.Regular(40, -1, 1), axis.Regular(9, 15, 60)),
#               lambda vg: vg.Electrons.get_hist2d('mvaTTH', 'pt', -1))
# ]
# mask(plotter, truth=None)
# plotter.fill_hists(graph_2d)
# hists_tot = plotter.get_hists('mva_2d')

# tot_hist = Histogram("total", *graph_2d[0].bins())
# for key, h in hists_tot.items():
#     tot_hist += h

# mask(plotter, truth=True)
# plotter.fill_hists(graph_2d)
# hists_true = plotter.get_hists('mva_2d')
# true_hist = Histogram("true", *graph_2d[0].bins())
# for key, h in hists_true.items():
#     true_hist += h

# print((true_hist / tot_hist).vals)

# for graph in graphs:
#     plot_eff(plotter, graph, output/f'{graph.name}_EE_eff.png', "EE")
#     plot_eff(plotter, graph, output/f'{graph.name}_EB_eff.png', "EB")
# # exit()

# version = 'none'

# # Fake Rates
# mask(plotter, truth=False)
# plotter.fill_hists(graphs)
# for graph in graphs:
#     plot_stack(plotter, graph.name, output/f'{graph.name}_fake_{version}.png')

# mask(plotter, truth=False, eta_split="EB")
# plotter.fill_hists(graphs)
# for graph in graphs:
#     plot_stack(plotter, graph.name, output/f'{graph.name}_fake_EB_{version}.png')

# mask(plotter, truth=False, eta_split="EE")
# plotter.fill_hists(graphs)
# for graph in graphs:
#     plot_stack(plotter, graph.name, output/f'{graph.name}_fake_EE_{version}.png')

# # True Rates
# mask(plotter, truth=True)
# plotter.fill_hists(graphs)
# for graph in graphs:
#     plot_stack(plotter, graph.name, output/f'{graph.name}_true_{version}.png', sig='ttt')

# mask(plotter, truth=True, eta_split="EB")
# plotter.fill_hists(graphs)
# for graph in graphs:
#     plot_stack(plotter, graph.name, output/f'{graph.name}_true_EB_{version}.png', sig='ttt')

# mask(plotter, truth=True, eta_split="EE")
# plotter.fill_hists(graphs)
# for graph in graphs:
#     plot_stack(plotter, graph.name, output/f'{graph.name}_true_EE_{version}.png', sig='ttt')
