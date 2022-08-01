#!/usr/bin/env python3
import uproot
import numpy as np
from dataclasses import dataclass
import awkward as ak
from typing import ClassVar
from sklearn.metrics import roc_curve
from scipy.stats import kurtosis

from analysis_suite.Variable_Creator.data_processor import DataProcessor
from analysis_suite.commons.histogram import Histogram
import analysis_suite.commons.configs as config
from analysis_suite.commons.plot_utils import plot, nonratio_plot, ratio_plot
from analysis_suite.commons.fake_rate_helper import make_stack

@dataclass
class GraphInfo:
    name: str
    axis_name: str
    bin_tuple: object
    func: object
    info: ClassVar[object] = None
    lumi: ClassVar[float] = None
    cuts: object = None

    def bins(self):
        if isinstance(self.bin_tuple, tuple):
            return self.bin_tuple
        else:
            return (self.bin_tuple,)

class Plotter:
    fom_label = {'likely': '\sqrt{-2\ln{\lambda(0)}}',
                 's/sqrtb': 'S/\sqrt{B}'}
    sig_colors = ['Signal', '#ef8a62'] #     '#f1a340' '#7fbf7b' '#d8b365' '#ef8a62'
    bkg_colors = ['Background', "#67a9cf"] # '#998ec3' '#af8dc3' '#5ab4ac' "#67a9cf"

    def __init__(self, filename, groups, cuts=None, **kwargs):
        self.dfs = dict()
        self.sig = ""
        self.bkg = []
        self.groups = groups
        self.readtype = 'flat'
        self.graphs = dict()

        if self.readtype == 'flat':
            self.setup_flat(filename, groups, cuts)
        elif self.readtype == 'ntuple':
            self.setup_ntuple(filename, groups, cuts, **kwargs)

    def set_groups(self, sig, bkg):
        self.sig = sig
        self.bkg = bkg

    def get_hists(self, name):
        return self.hists[name]

    def setup_ntuple(self, filenames, groups, cut=None, year="2016", trees='Signal_Dilepton', systName="Nominal"):
        data = DataProcessor(filenames, PlotInfo.lumi[year], systName, cut=cut)
        if isinstance(trees, list):
            for tree in trees:
                data.process_year(filenames, trees)
        else:
            data.process_year(filenames, trees)
        for group, members in groups.items():
            self.dfs[group] = pd.DataFrame()
            for member in members:
                if member not in data.final_set:
                    continue
                self.dfs[group] = pd.concat([self.dfs[group], data.final_set[member]])
            if not len(self.dfs[group]):
                del self.dfs[group]


    def setup_flat(self, filename, groups, cuts):
        with uproot.open(filename) as f:
           for group, members in groups.items():
               self.dfs[group] = ak.Array([])
               for member in members:
                   if member not in f:
                       continue
                   self.dfs[group] = ak.concatenate((self.dfs[group], self.apply_cut(f[member].arrays(), cuts)))
               if not len(self.dfs[group]):
                   del self.dfs[group]


    def fill_hists(self, graphs, ginfo=None, subset=None):
        self.graphs = graphs if isinstance(graphs, dict) else {graphs.name: graphs}
        self.hists = {name: dict() for name in self.graphs.keys() }
        for group in self.groups:
            for name, graph in self.graphs.items():
                if subset is not None and name not in subset:
                    continue
                self.hists[name][group] = Histogram(group, *graph.bins())
                if self.readtype == 'flat' and group in self.dfs:
                    self.hists[name][group].fill(np.nan_to_num(self.dfs[group][graph.func], nan=-1000),
                                                 weight=self.dfs[group]['scale_factor'])
                elif self.readtype == 'ntuple':
                    pass
                else:
                    del self.hists[name][group]
                if ginfo is not None:
                    self.hists[name][group].set_plot_details(ginfo)

    def apply_cut(self, df, cut):
        if cut is None:
            return df
        elif isinstance(cut, list):
            for cut_ in cut:
                df = apply_cut(df, cut_)
            return df
        else:
            if len(cutter := cut.split("<")) > 1:
                return df[df[cutter[0]] < float(cutter[1])]
            elif len(cutter := cut.split(">")) > 1:
                return df[df[cutter[0]] > float(cutter[1])]
            elif len(cutter := cut.split("==")) > 1:
                return df[df[cutter[0]] == float(cutter[1])]
            else:
                raise Exception(f'{cut} is not formatted correctly')


    def get_fom(self, var, bins, sig_type="likely"):
        s_tot = np.array([np.sum(self.dfs[self.sig]["scale_factor"][self.dfs[self.sig][var] > val]) for val in bins])
        b_tot = np.zeros(len(bins))
        for bkg in self.bkg:
            if bkg not in self.dfs: continue
            b_tot += np.array([np.sum(self.dfs[bkg]["scale_factor"][self.dfs[bkg][var] > val]) for val in bins])
        if sig_type == "likely":
            return config.likelihood_sig(s_tot, b_tot)
        elif sig_type == "s/sqrtb":
            return config.asymptotic_sig(s_tot, b_tot)
        else:
            raise AttributeError(f"Significance type {sig_type} not allowed!")


    def plot_fom(self, name, bins=None, sig_type='likely'):
        graph = self.graphs[name]
        if bins is None:
            bins = graph.bin_tuple.edges
        discrete = np.unique(np.diff(bins)) == 1

        with plot(f'test.png') as ax:
            sig = self.get_sum(self.sig, graph, self.sig_colors)
            bkg = self.get_sum(self.bkg, graph, self.bkg_colors)
            sig.plot_shape(ax)
            bkg.plot_shape(ax)
            ax.set_ylabel("Normalized Events (A.U)")
            ax.set_xlabel(graph.axis_name, horizontalalignment='right', x=1.0)

            # Plot FOM details
            ax2 = ax.twinx()
            fom = self.get_fom(graph.func, bins, sig_type)
            maxbin = bins[np.argmax(fom)]
            full_label = f'${self.fom_label[sig_type]}={max(fom):.3f}$\n cut={maxbin:.2f}'
            if discrete:
                ax2.hist(bins=bins, x=bins, weights=fom, label=full_label, histtype='step', edgecolor='k')
            else:
                ax2.plot(bins, fom, label=full_label, color = 'k')
            ax2.plot(np.linspace(bins[0], bins[-1], 5), [max(fom)]*5, linestyle=':', color='k')
            ax2.set_ylabel("FOM", horizontalalignment='right', y=1.0)

            # Make Legend
            lines, labels = ax2.get_legend_handles_labels()
            lines2, labels2 = ax.get_legend_handles_labels()
            ax2.legend(lines + lines2, labels + labels2)

        return (fom, maxbin)

    def plot_shape(self, name):
        graph = self.graphs[name]
        with nonratio_plot(f'test.png', graph.axis_name, graph.bin_tuple.edges) as ax:
            sig = self.get_sum(self.sig, graph, self.sig_colors)
            bkg = self.get_sum(self.bkg, graph, self.bkg_colors)
            sig.plot_shape(ax)
            bkg.plot_shape(ax)
            ax.set_ylabel("Normalized Events (A.U)")

    def plot_roc(self):
        with nonratio_plot(f'roc.png', "False Positive Rate", [0., 1.]) as ax:
            pred, truth, weight = np.array([]), np.array([]), np.array([])
            for group, df in self.dfs.items():
                if group in self.sig:
                    truth = np.append(truth, np.ones(len(df)))
                else:
                    truth = np.append(truth, np.zeros(len(df)))
                weight = np.append(weight, df['scale_factor'])
                pred = np.append(pred, df['Signal'])
            fpr, tpr, _ = roc_curve(truth, pred, sample_weight=weight)
            auc = np.trapz(tpr, fpr)

            ax.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
            ax.plot(np.linspace(0, 1, 5), np.linspace(0, 1, 5), linestyle=':')
            ax.set_ylabel("True Positive Rate", horizontalalignment='right', y=1.0)

    def plot_stack(self, name):
        graph = self.graphs[name]
        signal = self.hists[graph.name][self.sig]
        stack = make_stack({key: val for key, val in self.hists[graph.name].items() if key != self.sig})

        with ratio_plot('test.png', graph.axis_name, stack.get_xrange()) as ax:
            ratio = Histogram("Ratio", graph.bin_tuple, color="black")
            band = Histogram("Ratio", graph.bin_tuple, color="plum")
            error = Histogram("Stat Errors", graph.bin_tuple, color="plum")

            ratio += signal/stack
            band += stack/stack
            for hist in stack.stack:
                error += hist

            pad, subpad = ax

            #upper pad
            stack.plot_stack(pad)
            signal.plot_points(pad)
            error.plot_band(pad)

            # ratio pad
            ratio.plot_points(subpad)
            band.plot_band(subpad)

            hep.cms.label(ax=pad, lumi=lumi, data=data)


    def get_sum(self, groups, graph, details=None):
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
        print("| name      | events | raw events | mean+-std | kurtosis |")
        info = list()
        for group, df in self.dfs.items():
            var = df[varname]

            info.append({
                "name": group,
                "mean": np.mean(var),
                "std":  np.std(var),
                "kurtosis": kurtosis(var),
                'nRaw': len(var),
                'nEvents': np.sum(df.scale_factor)
            })
        info = sorted(info, reverse=True, key=lambda x: x["mean"])
        print("-"*50)
        for arr in info:
            print(f'|{arr["name"]:10} | {arr["nEvents"]:.2f}| {arr["nRaw"]} | {arr["mean"]:.2f}+-{arr["std"]:.2f} | {arr["kurtosis"]:0.2f} |')
        print("-"*50)