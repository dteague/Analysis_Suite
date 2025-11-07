#!/usr/bin/env python3
import uproot
import numpy as np
from dataclasses import dataclass
from importlib import import_module
from typing import Callable

from analysis_suite.flatten import NtupleGetter, FlatGetter
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.info import fileInfo
from analysis_suite.commons.user import analysis_area

@dataclass
class GraphInfo:
    axis_name: str = ""
    bin_tuple: object = None
    func: object = None

    def bins(self):
        if isinstance(self.bin_tuple, tuple):
            return self.bin_tuple
        else:
            return (self.bin_tuple,)

    def dim(self):
        return len(self.bins())

    def edges(self, i=0):
        if isinstance(self.bin_tuple, tuple):
            return self.bin_tuple[i].edges
        else:
            return self.bin_tuple.edges

    def get_axis_name(self, *args, **kwargs):
        return self.axis_name.format(*args, **kwargs)

class HistGetter:
    def __init__(self, ntuple_info, year, cores=1, **kwargs):
        self.dfs = dict()
        self.lumi = lumi[year]
        self.year = year
        self.cores = cores
        self.ntuple = ntuple_info
        self.systs = []
        self.scales = []
        self.scaled = []
        self.workdir = kwargs.get('workdir', analysis_area/"data")
        self.limit_samples = kwargs.get('limit', True)
        if kwargs.get('scales', False):
            for scale in kwargs['scales']:
                self.scales.append(import_module(f"analysis_suite.data.scales.{scale}").scale)
        filename = kwargs.get('filename', None)
        if filename is not None:
            self.setup_flat(filename, cuts=kwargs.get('mask', None))
        else:
            filename = ntuple_info.get_filename(year, **kwargs)
            self.root_files = []
            if filename.is_dir():
                for root_file in filename.glob("*root"):
                    self.setup_ntuple(root_file, **kwargs)
            else:
                self.setup_ntuple(filename, **kwargs)

    def df_iter(self, members=None):
        for member, df in self.dfs.items():
            if members is not None and member not in members:
                continue
            if isinstance(df, dict):
                for tree, subdf in df.items():
                    group = self.ntuple.get_group_name(member, tree)
                    if group is None:
                        # print(group, member, "PROBLEM!")
                        continue
                    yield group, member, subdf
            else:
                group = self.ntuple.get_group_name(member, None, False)
                if group is None:
                    # print(group, member, "PROBLEM!")
                    continue
                yield group, member, df

    def setup_ntuple(self, root_file, systName="Nominal", **kwargs):
        f = uproot.open(root_file,
                        decompression_executor=uproot.ThreadPoolExecutor(self.cores)
        )
        self.root_files.append(f)
        members = [m for m in f.keys(recursive=False, cycle=False)]
        for member in members:
            xsec = 1. if fileInfo.is_data(member) else fileInfo.get_xsec(member)*self.lumi*1000
            for tree in self.ntuple.trees:
                if self.limit_samples and self.ntuple.get_group_name(member, tree) is None:
                    continue
                vg = NtupleGetter(f, tree, member, xsec, systName=systName, cuts=self.ntuple.cut)
                if not vg.tree:
                    continue
                self.ntuple.setup_branches(vg)
                vg.set_systematic(systName)

                self.systs = np.unique(np.concatenate([self.systs, vg.all_systs]))
                group = self.ntuple.get_group_name(member, tree)
                self._internal_scale(vg, group, member)
                if member not in self.dfs:
                    self.dfs[member] = dict()
                self.dfs[member][tree] = vg

    def setup_flat(self, filename, cuts=None):
        with uproot.open(filename) as f:
            if self.year in f:
                f = f[self.year]
            members = [m for m in f.keys(recursive=False, cycle=False)]
            for member in members:
                if self.limit_samples and self.ntuple.get_group_name(member, None, False) is None:
                    continue
                fg = FlatGetter(f, member)
                if not fg:
                    continue
                self.systs = np.unique(np.concatenate([self.systs, fg.list_systs()]))
                group = self.ntuple.get_group_name(member, None, False)
                self._internal_scale(fg, group, member)
                fg.cut(cuts)
                self.dfs[member] = fg

    def cut(self, mask, groups=None):
        for group, member, df in self.df_iter(self.dfs.keys()):
            if groups is not None and group not in groups:
                continue
            df.cut(mask)

    def mask(self, mask, groups=None):
        for group, member, df in self.df_iter(self.dfs.keys()):
            if groups is not None and group not in groups:
                continue
            df.mask = mask


    def _internal_scale(self, df, group, member):
        if self.ntuple.get_info().is_data_driven(group) or group == 'data':
            return
        for scale in self.scales:
            # print(group, member, df.syst_name)
            scale(df, self.workdir, group, member, self.year, df.syst_name)
            self.scaled.append(df.syst_name)

    def scale(self, scale, groups=None):
        if isinstance(groups, str):
            groups = [groups]
        for group, member, df in self.df_iter(self.dfs.keys()):
            if groups is not None and group not in groups:
                continue
            if isinstance(scale, Callable):
                df.scale = scale(df)
            else:
                df.scale = scale

    def reset_mask(self):
        for group, member, df in self.df_iter(self.dfs.keys()):
            df.reset()

    def reset(self):
        for _, _, df in self.df_iter():
            df.set_systematic("Nominal")

    def reset_syst(self, systName, members=None):
        for group, member, df in self.df_iter(members):
            df.set_systematic(systName)
            if not df.correct_syst:
                continue
            # if df.syst_name is None or df.syst_name in self.scaled:
            #     continue
            self._internal_scale(df, group, member)

    def get_hist(self, graph, *args, **kwargs):
        hists = {}
        members = kwargs.get("members", self.dfs.keys())
        axis_name = graph.get_axis_name(**kwargs)
        for group, member, df in self.df_iter(members):
            vals, weight = df.get_graph(graph, *args)
            # if "wjet" in group and member != "wjets_ht100-200":
            #     continue
            # if "wjet" in group or 'ttbar' in group:
            #     print(group, member)
            #     print(np.histogram(vals, weights=weight, bins=np.linspace(0,750,21))[0])
            if len(vals) == 0:
                continue
            if group not in hists:
                hists[group] = Histogram(*graph.bins(), axis_name=axis_name)
            if graph.dim() == 1:
                hists[group].fill(vals, weight=weight, member=member, flow=kwargs.get('flow', False))
            else:
                hists[group].fill(*vals, weight=weight, member=member, flow=kwargs.get('flow', False))
        if kwargs.get('fix_negative', False):
            for group, hist in hists.items():
                if min(hist.vals) < 0:
                    hist.values()[:] = np.abs(hist.vals)
        return hists

    def get_hists(self, graphs, *args, **kwargs):
        return {
            key: self.get_hist(graph, *args, **kwargs) for key, graph in graphs.items()
        }

    def get_syst_hist(self, graph, members, systName, *args, **kwargs):
        hists = {}
        self.reset_syst(systName+"_up", members)
        hist_up = self.get_hist(graph, *args, members=members, **kwargs)
        self.reset_syst(systName+"_down", members)
        hist_down = self.get_hist(graph, *args, members=members, **kwargs)

        for member in members:
            if member not in hist_up:
                continue
            if member not in hists:
                hists[member] = {}
            for tree in hist_up[member].keys():
                up, down = hist_up[member][tree], hist_down[member][tree]
                hists[member][tree] = (up-down)/2
        return hists


    def setup_systs(self, syst_hists, systs, chan):
        out_hists = []
        groups = nom_hists.keys()
        for syst in systs:
            if syst.syst_type == 'lnN':
                out_hists.append({g: syst.get_lnN_error(self.year, chan) for g in groups
                             if syst.good_syst(group, self.year, chan)})
            else:
                up = syst_hists[syst.get_name(year)+"_up"]
                down = syst_hists[syst.get_name(year)+"_down"]
                out_hists.append((up-down)/2)

        return out_hists
