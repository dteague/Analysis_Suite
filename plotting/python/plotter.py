#!/usr/bin/env python3
import uproot
import numpy as np
from dataclasses import dataclass
import awkward as ak
from typing import ClassVar
from sklearn.metrics import roc_curve
from scipy.stats import kurtosis
from copy import copy
from prettytable import PrettyTable
from pathlib import Path

import analysis_suite.flatten as getter
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.constants import lumi
import analysis_suite.commons.configs as config
from analysis_suite.commons.info import fileInfo
from analysis_suite.commons.plot_utils import plot, nonratio_plot, ratio_plot, hep

from .utils import likelihood_sig, get_syst_index
from .stack import Stack

@dataclass
class GraphInfo:
    name: str = ""
    axis_name: str = ""
    bin_tuple: object = None
    func: object = None
    output: str = 'stack'
    info: ClassVar[object] = None
    cuts: object = None

    def bins(self):
        if isinstance(self.bin_tuple, tuple):
            return self.bin_tuple
        else:
            return (self.bin_tuple,)

    def dim(self):
        return len(self.bins())

    def edges(self, i=0):
        if isinstance(self.bin_tuple, tuple):
            return self.bin_tuple[0].edges
        else:
            return self.bin_tuple.edges


class Plotter:
    fom_label = {'likely': '\sqrt{-2\ln{\lambda(0)}}',
                 's/sqrtb': 'S/\sqrt{B}'}
    sig_colors = ['Signal', '#ef8a62'] #     '#f1a340' '#7fbf7b' '#d8b365' '#ef8a62'
    bkg_colors = ['Background', "#67a9cf"] # '#998ec3' '#af8dc3' '#5ab4ac' "#67a9cf"

    def __init__(self, filename, groups, year, ntuple=None, cuts=None, **kwargs):
        self.dfs = dict()
        self.sig = ""
        self.bkg = []
        self.data = 'data'
        self.groups = copy(groups)
        self.graphs = dict()
        self.lumi = lumi[year]
        self.workdir = Path('.')

        if ntuple is None:
            self.setup_flat(filename, cuts)
        else:
            self.setup_ntuple(filename, ntuple, **kwargs)
        self.clean_groups(groups)


    def set_groups(self, sig="", bkg=None, data=None):
        self.sig = sig
        self.bkg = list() if bkg is None else [b for b in bkg if b in self.groups]
        if data is not None:
            self.data = data

    def clean_groups(self, groups):
        for group, members in groups.items():
            new_members = [member for member in members if member in self.dfs]
            self.groups[group] = new_members
        

    def get_hists(self, name, group=None):
        if group is None:
            return self.hists[name]
        else:
            return self.hists[name][group]


    def setup_ntuple(self, filenames, ntuple, systName="Nominal"):
        root_files = filenames.glob("*.root") if filenames.is_dir() else [filenames]
        for root_file in root_files:
            syst = get_syst_index(root_file, systName)
            if syst == -1:
                continue
            for group, members in self.groups.items():
                for member in members:
                    xsec = 1. if fileInfo.is_data(member) else fileInfo.get_xsec(member)*self.lumi*1000
                    for tree in ntuple.trees:
                        if ntuple.ignore(tree, group):
                            continue
                        vg = getter.NtupleGetter(root_file, tree, member, xsec, syst, systName)
                        if not vg:
                            continue
                        name = ntuple.get_change(tree, member)
                        if name not in self.dfs:
                            self.dfs[name] = dict()
                        ntuple.setup_branches(vg)
                        ntuple.apply_cut(vg)
                        self.dfs[name][tree] = vg

    def setup_flat(self, filename, cuts):
        remove = list()
        with uproot.open(filename) as f:
            for group, members in self.groups.items():
                for member in members:
                    fg = getter.FlatGetter(f, member)
                    if not fg:
                        continue
                    fg.cut(cuts)
                    self.dfs[member] = fg

    def get_integral(self):
        output = {keys: 0. for keys in self.dfs.keys()}
        for key, vg in self.dfs.items():
            if isinstance(vg, dict):
                for subvg in vg.values():
                    output[key] += sum(subvg.scale)
            else:
                output[key] += sum(vg.scale)
        return output

    def fill_hists(self, graphs, ginfo, *args, subset=None):
        if isinstance(graphs, dict):
            self.graphs = graphs
        elif isinstance(graphs, list):
            self.graphs = {graph.name: graph for graph in graphs}
        else:
            self.graphs = {graphs.name: graphs}
        self.hists = {name: dict() for name in self.graphs.keys() }
        for name, graph in self.graphs.items():
            if subset is not None and name not in subset:
                continue
            for group, members in self.groups.items():
                self.hists[name][group] = Histogram(group, *graph.bins(), group_info=ginfo)
                for member in members:
                    vals, weight = self.get_array(member, graph, *args)
                    self.hists[name][group].fill(*vals, weight=weight, member=member)

    def mask(self, mask, clear=True, groups=None):
        for key, vg in self.dfs.items():
            if groups is not None and key not in groups:
                continue
            if isinstance(vg, dict):
                for subvg in vg.values():
                    if clear:
                        subvg.clear_mask()
                    subvg.mask = mask
            else:
                if clear:
                    vg.clear_mask()
                vg.mask = mask

    def cut(self, mask, groups=None):
        def _cut(vg, mask):
            if isinstance(vg, dict):
                for key, subvg in vg.items():
                    subvg.cut(mask)
            else:
                vg.cut(mask)

        if groups is None:
            for vg in self.dfs.values():
                _cut(vg, mask)
        elif isinstance(groups, str):
            _cut(self.dfs[groups], mask)
        else:
            for group in groups:
                for member in self.groups[group]:
                    _cut(self.dfs[member], mask)


    def scale(self, scaler, groups=None):
        def _scale(vg, scaler):
            if isinstance(vg, dict):
                for key, subvg in vg.items():
                    scaler(subvg)
            else:
                scaler(vg)

        if groups is None:
            for vg in self.dfs.values():
                _scale(vg, scaler)
        elif isinstance(groups, str):
            for member in self.groups[groups]:
                _scale(self.dfs[member], scaler)
        else:
            for group in groups:
                for member in self.groups[group]:
                    _scale(self.dfs[member], scaler)


    def scale_hists(self, group, scale, scale_df=False):
        for name, hist in self.hists.items():
            hist[group].scale(scale, changeName=True)
        if not scale_df:
            return

        for member in self.groups[group]:
            if isinstance(self.dfs[member], dict):
                for subdf in self.dfs[member].values():
                    subdf.scale = scale
            else:
                self.dfs[member].scale = scale


    def get_array(self, member, graph, *args):
        if isinstance(graph, str):
            graph = self.graphs[graph]

        if isinstance(self.dfs[member], dict):
            vals = [[] for _ in range(graph.dim())]
            scales = np.array([])
            for tree, vg in self.dfs[member].items():
                v, s = vg.get_graph(graph, *args)
                if len(s) == 0:
                    continue
                if graph.dim() == 1:
                    vals = [np.concatenate((val, v)) for i, val in enumerate(vals)]
                else:
                    vals = [np.concatenate((val, v[i])) for i, val in enumerate(vals)]
                scales = np.concatenate((scales, s))
            return vals, scales
        else:
            vals, scales = self.dfs[member].get_graph(graph, *args)
            return [vals], scales


    def get_fom(self, graph, bins, sig_type="likely"):
        s_tot = np.zeros(len(bins))
        b_tot = np.zeros(len(bins))
        for group, members in self.groups.items():
            for member in members:
                values, scales = self.get_array(member, graph)
                if group in self.sig:
                    b_tot += np.array([np.sum(scales[values > val]) for val in bins])
                elif group in self.bkg:
                    s_tot += np.array([np.sum(scales[values > val]) for val in bins])
        if sig_type == "likely":
            return likelihood_sig(s_tot, b_tot)
        elif sig_type == "s/sqrtb":
            return config.asymptotic_sig(s_tot, b_tot)
        else:
            raise AttributeError(f"Significance type {sig_type} not allowed!")


    def plot_from_graph(self, graph, outfile, **kwargs):
        if graph.output == "fom":
            self.plot_fom(graph.name, outfile, **kwargs)
        elif graph.output == "shape":
            self.plot_shape(graph.name, outfile, **kwargs)
        elif graph.output == 'roc':
            self.plot_roc(outfile)
        elif graph.output == 'stack':
            self.plot_stack(graph.name, outfile, **kwargs)
        else:
            raise AttributeError(f"Not correct Output type: graph.output=={graph.output}")
        

    def plot_fom(self, name, outfile, bins=None, sig_type='likely', chan=None, **kwargs):
        graph = self.graphs[name]
        if bins is None:
            bins = graph.bin_tuple.edges
        discrete = np.all(np.unique(np.diff(bins)) == 1)

        with nonratio_plot(self.workdir/outfile, graph.axis_name.format(chan), graph.bin_tuple.edges, legend=False) as ax:
            sig = self.get_sum(self.sig, graph, self.sig_colors)
            bkg = self.get_sum(self.bkg, graph, self.bkg_colors)
            sig.plot_shape(ax)
            bkg.plot_shape(ax)
            ax.set_ylabel("Normalized Events (A.U)")

            # Plot FOM details
            ax2 = ax.twinx()
            fom = self.get_fom(graph, bins, sig_type)
            maxbin = bins[np.argmax(fom)]
            full_label = f'${self.fom_label[sig_type]}={max(fom):.3f}$\n cut={maxbin:.2f}'
            if discrete:
                ax2.hist(bins=bins, x=bins, weights=fom, label=full_label, histtype='step', edgecolor='k')
            else:
                ax2.plot(bins, fom, label=full_label, color = 'k')
            ax2.plot(np.linspace(bins[0], bins[-1], 5), [max(fom)]*5, linestyle=':', color='k')
            ax2.set_ylabel("FOM", horizontalalignment='right', y=1.0)
            _, top_lim = ax2.get_ylim()
            ax2.set_ylim(top=top_lim*0.975)

            # Make Legend
            lines, labels = ax2.get_legend_handles_labels()
            lines2, labels2 = ax.get_legend_handles_labels()
            ax2.legend(lines + lines2, labels + labels2)

        return (fom, maxbin)


    def plot_shape(self, name, outfile, chan=None, **kwargs):
        graph = self.graphs[name]
        axis_name = graph.axis_name if chan is None else graph.axis_name.format(chan)
        with nonratio_plot(self.workdir/outfile, axis_name, graph.bin_tuple.edges) as ax:
            sig = self.get_sum(self.sig, graph, self.sig_colors)
            bkg = self.get_sum(self.bkg, graph, self.bkg_colors)
            sig.plot_shape(ax)
            bkg.plot_shape(ax)
            ax.set_ylabel("Normalized Events (A.U)")


    def plot_roc(self):
        with nonratio_plot(self.workdir/'roc.png', "False Positive Rate", [0., 1.]) as ax:
            pred, truth, weight = np.empty((3,0))
            for group, members in self.groups.items():
                for member in members:
                    df = self.dfs[member]
                    if group in self.sig:
                        truth = np.append(truth, np.ones(len(df)))
                    elif group in self.bkg:
                        truth = np.append(truth, np.zeros(len(df)))
                    weight = np.append(weight, df.scale)
                    pred = np.append(pred, df['Signal'])
            fpr, tpr, _ = roc_curve(truth, pred, sample_weight=weight)
            auc = np.trapz(tpr, fpr)

            ax.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
            ax.plot(np.linspace(0, 1, 5), np.linspace(0, 1, 5), linestyle=':')
            ax.set_ylabel("True Positive Rate", horizontalalignment='right', y=1.0)


    def plot_stack(self, name, outfile, chan=None, **kwargs):
        graph = self.graphs[name]
        signal = self.hists[name].get(self.sig, Histogram(""))
        data = copy(self.hists[name].get(self.data, Histogram("")))
        data.color = 'k'
        stack = self.make_stack(name)
        plotter = ratio_plot if data else nonratio_plot

        axis_name = graph.axis_name if chan is None else graph.axis_name.format(chan)
        with plotter(self.workdir/outfile, axis_name, stack.get_xrange(), **kwargs) as ax:
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
            stack.plot_stack(pad)
            signal.plot_points(pad)
            data.plot_points(pad)
            error.plot_band(pad)

            if (region := kwargs.get('region', False)):
                (subpad if subpad else pad).text(graph.edges()[0], -0.7, region.format(chan))

            # ratio pad
            ratio.plot_points(subpad)
            band.plot_band(subpad)

            hep.cms.label(ax=pad, lumi=self.lumi, data=data)


    def get_sum(self, groups, graph, details=None):
        if isinstance(graph, str):
            graph = self.graphs[graph]
        all_hist = Histogram("", *graph.bins())
        if isinstance(groups, str):
            all_hist += self.hists[graph.name][groups]
        else:
            for group in groups:
                if group not in self.hists[graph.name]:
                    continue
                all_hist += self.hists[graph.name][group]
        if details is not None:
            all_hist.set_plot_details(details)
        return all_hist


    def print_info(self, varname):
        info = list()
        graph = self.graphs[varname]
        for group, members in self.groups.items():
            for member in members:
                vals, weight = self.get_array(member, graph)
                info.append({
                    "name" : member,
                    "group" : group,
                    "mean": np.mean(vals),
                    "std":  np.std(vals),
                    "kurtosis": kurtosis(vals),
                    'nRaw': len(vals),
                    'nEvents': np.sum(weight)
                })

        table = PrettyTable(["Name", "Group", "Events", "Raw", "Mean", "Kurtosis"])
        info = sorted(info, reverse=True, key=lambda x: x["mean"])

        for arr in info:
            table.add_row([arr['name'], arr['group'], config.sig_fig(arr['nEvents']), arr['nRaw'],
                          f'{arr["mean"]:.2f}+-{arr["std"]:.2f}', config.sig_fig(arr['kurtosis'])])
        print(table.get_string())


    def make_stack(self, name):
        stack = Stack(*self.graphs[name].bins())
        for group in self.bkg:
            stack += self.hists[name][group]
        return stack



# for key, hist in file["shapes_fit_s/yr2018"].items():
#     key = key[:key.index(";")]
#     # Need to deal with Data
#     if "data" in key:
#         x, y = hist.member("fX"), hist.member("fY")
#         data += data.fill(x, weight=y)
#     if "TH1" not in repr(hist):
#         continue
#     hist = hist.to_boost()
#     if "total" in key:
#         continue
#     elif key == signalName:
#         signal += hist


#     groupHist = Histogram(key, hist.axes[0])
#     groupHist.set_plot_details(group_info)
#     groupHist += hist
#     stacker += groupHist
