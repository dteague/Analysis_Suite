#!/usr/bin/env python3
import numpy as np
from sklearn.metrics import roc_curve
from pathlib import Path
# from prettytable import PrettyTable
import warnings

from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.constants import lumi
import analysis_suite.commons.configs as config
from analysis_suite.commons.plot_utils import plot, nonratio_plot, ratio_plot, cms_label

from .utils import likelihood_sig, get_syst_index
from .stack import Stack

def get_a_dictval(hists):
    return next(iter(hists.values()))

class Plotter:
    fom_label = {'likely': r'\sqrt{-2\ln{\lambda(0)}}',
                 's/sqrtb': r'S/\sqrt{B}',
                 's/b': 'S/B'}
    sig_colors = ['Signal', '#ef8a62'] #     '#f1a340' '#7fbf7b' '#d8b365' '#ef8a62'
    bkg_colors = ['Background', "#67a9cf"] # '#998ec3' '#af8dc3' '#5ab4ac' "#67a9cf"

    def __init__(self, ntuple, year, sig="", bkg=None, data='data', out_format='png', outdir="."):
        # Setup groups
        self.ntuple = ntuple
        self.ginfo = ntuple.get_info(keep_dd_data=True)
        self.sig = sig
        if bkg is None:
            bkg = []
        elif bkg == "all":
            self.bkg = [g for g in self.ginfo.get_groups() if g not in [sig, data]]
        else:
            self.bkg = bkg
        self.data = data
        self.groups = self.ginfo.setup_groups([self.sig, self.data]+self.bkg)
        # Setup output name
        self.ft = out_format
        self.outdir = Path(outdir)
        self.outdir.mkdir(exist_ok=True)
        self.extra = str(year)
        self.year = year
        self.lumi = lumi[year]

    def set_hist_groups(self, hists):
        output = {}
        for member, subdict in hists.items():
            for tree, hist in subdict.items():
                group = self.get_group(member, tree)
                if group is None:
                    continue
                if group not in output:
                    color = self.ginfo.get_color(group)
                    name = self.ginfo.get_legend_name(group)
                    axis_name = hist.axis_name
                    output[group] = Histogram(hist.axis, color=color, name=name, axis_name=axis_name)
                output[group] += hist
        return output

    def set_extra_text(self, extra):
        self.extra = f"{year}_{extra}"

    def get_group(self, member, tree):
        for group, members in self.groups.items():
            if member in members:
                if self.ntuple.pass_group(tree, group):
                    return group
        return None

    def make_stack(self, hists):
        stack = Stack()
        for group in self.bkg:
            stack += self.hists[name][group]
        return stack

    def plot_hist(self, name, hist, **kwargs):
        self.plot_hists({name: hist}, **kwargs)

    def plot_hists(self, hist_dict, plot_type="stack", **kwargs):
        for name, hists in hist_dict.items():
            group_hists = self.set_hist_groups(hists)
            # if graph.output == "fom":
            #     self.plot_fom(graph.name, outfile, **kwargs)
            # elif graph.output == "shape":
            #     self.plot_shape(graph.name, outfile, **kwargs)
            # elif graph.output == 'roc':
            #     self.plot_roc(outfile)
            if plot_type == 'stack':
                self.plot_stack(group_hists, name, **kwargs)
            else:
                raise AttributeError(f"Not correct Output type: graph.output=={graph.output}")

    def add_total_systs(self, stack, nom_hists, hist_systs):
        for syst, syst_dict in hist_systs.items():
            for group, syst_shift in syst_dict.items():
                if isinstance(syst_shift, float):
                    stack.add_syst(syst_shift*nom_hists[group])
                else:
                    stack.add_syst(self.set_hist_groups(syst_shift)[group])

    def add_group_systs(self, hist, hist_group, hist_systs):
        for syst, syst_dict in hist_systs.items():
            for group, syst_shift in syst_dict.items():
                if hist_group != group:
                    continue
                if isinstance(syst_shift, float):
                    hist.add_syst(syst_shift*hist)
                else:
                    hist.add_syst(self.set_hist_groups(syst_shift)[group])

    def plot_stack(self, hists, name, **kwargs):
        a_hist = get_a_dictval(hists)
        axis = a_hist.axes
        plot_syst = 'syst_hists' in kwargs

        data = hists.get(self.data, Histogram())
        signal = hists.get(self.sig, Histogram())
        make_ratio = bool(data)

        stack = Stack(*axis)
        for bkg in self.bkg:
            if bkg in hists:
                stack += hists[bkg]
        if plot_syst:
            self.add_total_systs(stack, hists, kwargs['syst_hists'])
            self.add_group_systs(signal, self.sig, kwargs['syst_hists'])

        # Lower plot
        plot_func = ratio_plot if make_ratio else nonratio_plot
        if make_ratio:
            ratio = Histogram(*axis, color="black")
            band = stack.get_self_ratio()
            ratio += data/stack

        file_name = f"{self.outdir}/{name}_{self.extra}.{self.ft}"
        with plot_func(file_name, a_hist.axis_name, a_hist.get_xrange(), **kwargs) as ax:
            if make_ratio:
                pad, subpad = ax
            else:
                pad, subpad = ax, None

            stack.plot_stack(pad)
            signal.plot_hist(pad)
            data.plot_points(pad)
            # Ratio
            ratio.plot_points(subpad)
            # Error bands
            if plot_syst:
                stack.plot_error_band(pad, color='coral', hatch=r'\\', systs=True)
                band.plot_error_band(subpad, color='coral', hatch=r'\\', systs=True)
            else:
                stack.plot_error_band(pad)
                band.plot_error_band(subpad)

            if (region := kwargs.get('region', False)) and data:
                (subpad if subpad else pad).text(graph.edges()[0], -0.7, region.format(chan))

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cms_label(pad, lumi=self.lumi, hasData=data)

    # def set_groups(self, sig="", bkg=None, data=None):
    #     self.sig = sig
    #     self.bkg = list() if bkg is None else [b for b in bkg if b in self.groups]
    #     if data is not None:
    #         self.data = data


    # def clean_groups(self, groups):
    #     for group, members in groups.items():
    #         new_members = [member for member in members if member in self.dfs]
    #         self.groups[group] = new_members


    # def get_hists(self, name, group=None):
    #     if group is None:
    #         return self.hists[name]
    #     else:
    #         return self.hists[name][group]


    # def setup_ntuple(self, filenames, ntuple, systName="Nominal", **kwargs):
    #     root_files = filenames.glob("*.root") if filenames.is_dir() else [filenames]
    #     self.ginfo = ntuple.get_info()
    #     for root_file in root_files:
    #         for group, members in self.groups.items():
    #             for member in members:
    #                 xsec = 1. if fileInfo.is_data(member) else fileInfo.get_xsec(member)*self.lumi*1000
    #                 # if member not in uproot.open(root_file):
    #                 #     continue # nice for debugging
    #                 for tree in ntuple.trees:
    #                      if not ntuple.pass_group(tree, group):
    #                         continue
    #                     vg = NtupleGetter(root_file, tree, member, xsec, systName=systName)
    #                     if not vg:
    #                         continue
    #                     ntuple.setup_branches(vg)
    #                     ntuple.apply_part_cut(vg)
    #                     ntuple.apply_cut(vg)
    #                     # To avoid data-driven bkg being written out as 'data'
    #                     out_name = group if self.ginfo.is_data_driven(group) else member

    #                     if vg:
    #                         if out_name not in self.group_by_mem:
    #                             self.group_by_mem[out_name] = {}
    #                         self.group_by_mem[out_name][tree] = group
    #                         if out_name not in self.dfs:
    #                             self.dfs[out_name] = dict()
    #                         self.dfs[out_name][tree] = vg


    # def setup_flat(self, filename, cuts):
    #     region = filename.stem.split('_')[-1]
    #     with uproot.open(filename) as f:
    #         for group, members in self.groups.items():
    #             for member in members:
    #                 fg = FlatGetter(f, member)
    #                 if not fg:
    #                     continue
    #                 fg.cut(cuts)
    #                 if member not in self.group_by_mem:
    #                     self.group_by_mem[member] = {}
    #                 self.group_by_mem[member][region] = group
    #                 self.dfs[member] = {region: fg}

    # def get_integral(self):
    #     output = {}
    #     for vg, key, group in self.getters(keys=True):
    #         if group not in output:
    #             output[group] = 0.
    #         output[group] += sum(vg.scale)
    #     return output


    # def fill_hists(self, graphs, *args, subset=None, **kwargs):
    #     if isinstance(graphs, dict):
    #         self.graphs = graphs
    #     elif isinstance(graphs, list):
    #         self.graphs = {graph.name: graph for graph in graphs}
    #     else:
    #         self.graphs = {graphs.name: graphs}
    #     self.hists = {name: dict() for name in self.graphs.keys() }
    #     for name, graph in self.graphs.items():
    #         if subset is not None and name not in subset:
    #             continue
    #         self.hists[name] = {g: Histogram(g, *graph.bins(), group_info=self.ginfo) for g in self.groups.keys()}
    #         for vg, member, group in self.getters(keys=True):
    #             vals, weight = vg.get_graph(graph, *args)
    #             if graph.dim() == 1:
    #                 self.hists[name][group].fill(vals, weight=weight, member=member, **kwargs)
    #             else:
    #                 self.hists[name][group].fill(*vals, weight=weight, member=member, **kwargs)

    # def mask(self, mask, clear=True, groups=None):
    #     for vg in self.getters(groups):
    #         if clear:
    #             vg.clear_mask()
    #         vg.mask = mask


    # def cut(self, mask, groups=None):
    #     for vg in self.getters(groups):
    #         vg.cut(mask)


    # def getters(self, groups=None, members=None, keys=False):
    #     if isinstance(groups, str):
    #         groups = [groups]
    #     if isinstance(members, str):
    #         members = [members]

    #     for mem, vg in self.dfs.items():
    #         if members is not None and mem not in members:
    #             continue
    #         for subkey, subvg in vg.items():
    #             group = self.group_by_mem[mem][subkey]
    #             if groups is not None and group not in groups:
    #                 continue
    #             if keys:
    #                 yield subvg, mem, group
    #             else:
    #                 yield subvg

    # def mask_part(self, part, var, func):
    #     for vg in self.getters():
    #         vg[part].mask_part(var, func)
    #         vg.reset()

    # def reset_part(self, part):
    #     for vg in self.getters():
    #         vg[part].clear_mask()

    # def scale(self, scaler, groups=None, inputs=None):
    #     for vg in self.getters(groups):
    #         if inputs is not None:
    #             scaler(vg, *inputs)
    #         else:
    #             scaler(vg)

    # def scale_hists(self, groups, scale, scale_df=False):
    #     if isinstance(groups, str):
    #         groups = [groups]
    #     for name, hist in self.hists.items():
    #         for group in groups:
    #             hist[group].scale(scale, changeName=True)
    #     if scale_df:
    #         for vg in self.getters(group):
    #             vg.scale = scale


    # def get_fom(self, graph, bins, sig_type="likely"):
    #     s_tot = np.zeros(len(bins))
    #     b_tot = np.zeros(len(bins))
    #     for df, member, group in self.getters(keys=True):
    #         values, scales = vg.get_graph(graph, *args)
    #         if group in self.bkg:
    #             b_tot += np.array([np.sum(scales[values > val]) for val in bins])
    #         elif group in self.sig:
    #             s_tot += np.array([np.sum(scales[values > val]) for val in bins])
    #     if sig_type == "likely":
    #         return likelihood_sig(s_tot, b_tot)
    #     elif sig_type == "s/sqrtb":
    #         return config.asymptotic_sig(s_tot, b_tot)
    #     else:
    #         raise AttributeError(f"Significance type {sig_type} not allowed!")


    # def plot_from_graph(self, graph, outfile, **kwargs):
    #     if graph.output == "fom":
    #         self.plot_fom(graph.name, outfile, **kwargs)
    #     elif graph.output == "shape":
    #         self.plot_shape(graph.name, outfile, **kwargs)
    #     elif graph.output == 'roc':
    #         self.plot_roc(outfile)
    #     elif graph.output == 'stack':
    #         self.plot_stack(graph.name, outfile, **kwargs)
    #     else:
    #         raise AttributeError(f"Not correct Output type: graph.output=={graph.output}")


    # def plot_fom(self, name, outfile, bins=None, sig_type='likely', chan=None, **kwargs):
    #     graph = self.graphs[name]
    #     if bins is None:
    #         bins = graph.bin_tuple.edges
    #     discrete = np.all(np.unique(np.diff(bins)) == 1)
    #     axis_name = graph.axis_name if chan is None else graph.axis_name.format(chan)

    #     with nonratio_plot(self.workdir/outfile, axis_name, graph.bin_tuple.edges, legend=False) as ax:
    #         sig = self.get_sum(self.sig, graph, self.sig_colors)
    #         bkg = self.get_sum(self.bkg, graph, self.bkg_colors)
    #         sig.plot_shape(ax)
    #         bkg.plot_shape(ax)
    #         ax.set_ylabel("Normalized Events (A.U)")

    #         # Plot FOM details
    #         ax2 = ax.twinx()
    #         fom = self.get_fom(graph, bins, sig_type)
    #         maxbin = bins[np.argmax(fom)]
    #         full_label = f'${self.fom_label[sig_type]}={max(fom):.3f}$\n cut={maxbin:.2f}'
    #         if discrete:
    #             ax2.hist(bins=bins, x=bins, weights=fom, label=full_label, histtype='step', edgecolor='k')
    #         else:
    #             ax2.plot(bins, fom, label=full_label, color = 'k')
    #         ax2.plot(np.linspace(bins[0], bins[-1], 5), [max(fom)]*5, linestyle=':', color='k')
    #         ax2.set_ylabel("FOM", horizontalalignment='right', y=1.0)
    #         _, top_lim = ax2.get_ylim()
    #         ax2.set_ylim(top=top_lim*0.975)

    #         # Make Legend
    #         lines, labels = ax2.get_legend_handles_labels()
    #         lines2, labels2 = ax.get_legend_handles_labels()
    #         ax2.legend(lines + lines2, labels + labels2)

    #     return (fom, maxbin)


    # def plot_shape(self, name, outfile, chan=None, **kwargs):
    #     graph = self.graphs[name]
    #     axis_name = graph.axis_name if chan is None else graph.axis_name.format(chan)
    #     with nonratio_plot(self.workdir/outfile, axis_name, graph.bin_tuple.edges) as ax:
    #         sig = self.get_sum(self.sig, graph, self.sig_colors)
    #         bkg = self.get_sum(self.bkg, graph, self.bkg_colors)
    #         sig.plot_shape(ax)
    #         bkg.plot_shape(ax)
    #         ax.set_ylabel("Normalized Events (A.U)")


    # def plot_roc(self):
    #     with nonratio_plot(self.workdir/'roc.png', "False Positive Rate", [0., 1.]) as ax:
    #         pred, truth, weight = np.empty((3,0))
    #         for df, member, group in self.getters(keys=True):
    #             if group in self.sig:
    #                 truth = np.append(truth, np.ones(len(df)))
    #             elif group in self.bkg:
    #                 truth = np.append(truth, np.zeros(len(df)))
    #             weight = np.append(weight, df.scale)
    #             pred = np.append(pred, df['Signal'])
    #         fpr, tpr, _ = roc_curve(truth, pred, sample_weight=weight)
    #         auc = np.trapz(tpr, fpr)

    #         ax.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
    #         ax.plot(np.linspace(0, 1, 5), np.linspace(0, 1, 5), linestyle=':')
    #         ax.set_ylabel("True Positive Rate", horizontalalignment='right', y=1.0)




    # def get_sum(self, groups, graph, details=None):
    #     if isinstance(graph, str):
    #         graph = self.graphs[graph]
    #     all_hist = Histogram("", *graph.bins())
    #     if isinstance(groups, str):
    #         all_hist += self.hists[graph.name][groups]
    #     else:
    #         for group in groups:
    #             if group not in self.hists[graph.name]:
    #                 continue
    #             all_hist += self.hists[graph.name][group]
    #     if details is not None:
    #         all_hist.set_plot_details(details)
    #     return all_hist


    # def print_info(self, varname):
    #     info = list()
    #     graph = self.graphs[varname]
    #     for df, member, group in self.getters(keys=True):
    #         vals, weight = df.get_graph(graph)
    #         info.append({
    #             "name" : member,
    #             "group" : group,
    #             "mean": np.mean(vals),
    #             "std":  np.std(vals),
    #             'nRaw': len(vals),
    #             'nEvents': np.sum(weight)
    #         })

    #     table = PrettyTable(["Name", "Group", "Events", "Raw", "Mean", "Kurtosis"])
    #     info = sorted(info, reverse=True, key=lambda x: x["mean"])

    #     for arr in info:
    #         table.add_row([arr['name'], arr['group'], config.sig_fig(arr['nEvents']), arr['nRaw'],
    #                       f'{arr["mean"]:.2f}+-{arr["std"]:.2f}', config.sig_fig(arr['kurtosis'])])
    #     print(table.get_string())


