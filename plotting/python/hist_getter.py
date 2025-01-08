#!/usr/bin/env python3
import uproot
from dataclasses import dataclass

from analysis_suite.flatten import NtupleGetter, FlatGetter
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.info import fileInfo

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
        if ntuple_info is None:
            self.setup_flat(filename, cuts)
        else:
            filename = ntuple_info.get_filename(year)
            self.root_files = []
            if filename.is_dir():
                for root_file in filename.glob("*root"):
                    self.setup_ntuple(root_file, ntuple_info, **kwargs)
            else:
                self.setup_ntuple(filename, ntuple_info, **kwargs)
            self.ginfo = ntuple_info.get_info()

    def setup_ntuple(self, root_file, ntuple, systName="Nominal", **kwargs):
        f = uproot.open(root_file,
                        decompression_executor=uproot.ThreadPoolExecutor(self.cores)
        )
        self.root_files.append(f)
        members = [m for m in f.keys(recursive=False, cycle=False)]
        for member in members:
            xsec = 1. if fileInfo.is_data(member) else fileInfo.get_xsec(member)*self.lumi*1000
            for tree in ntuple.trees:
                vg = NtupleGetter(f, tree, member, xsec, systName=systName)
                if not vg:
                    continue
                ntuple.setup_branches(vg)
                ntuple.apply_part_cut(vg)
                ntuple.apply_cut(vg)
                if not vg:
                    continue

                if member not in self.dfs:
                    self.dfs[member] = dict()
                self.dfs[member][tree] = vg

    def setup_flat(self, filename, cuts):
        region = filename.stem.split('_')[-1]
        with uproot.open(filename) as f:
            for group, members in self.groups.items():
                for member in members:
                    fg = FlatGetter(f, member)
                    if not fg:
                        continue
                    fg.cut(cuts)
                    if member not in self.group_by_mem:
                        self.group_by_mem[member] = {}
                    self.group_by_mem[member][region] = group
                    self.dfs[member] = {region: fg}

    def reset(self):
        for subdf in self.dfs.values():
            for df in subdf.values():
                df.set_systematic("Nominal")

    def reset_syst(self, systName, members=None):
        if members is None:
            members = self.dfs.keys()

        for member in members:
            if member not in self.dfs:
                continue
            for tree, df in self.dfs[member].items():
                df.set_systematic(systName)

    def get_hist(self, graph, *args, **kwargs):
        members = kwargs.get("members", self.dfs.keys())
        hists = {}
        axis_name = graph.get_axis_name(**kwargs)
        for member in members:
            if member not in self.dfs:
                continue
            hists[member] = {}
            for tree, df in self.dfs[member].items():
                vals, weight = df.get_graph(graph, *args)
                hists[member][tree] = Histogram(*graph.bins(), axis_name=axis_name)
                if graph.dim() == 1:
                    hists[member][tree].fill(vals, weight=weight)
                else:
                    hists[member][tree].fill(*vals, weight=weight)
        return hists

    def load_variables(self, variables):
        for member, subdf in self.dfs.items():
            for tree, df in subdf.items():
                for var in variables:
                    if "/" in var:
                        part, var = var.split("/")
                        df[part][var, 0]
                    else:
                        df[var]

    def get_hists(self, graphs, *args):
        return {
            key: self.get_hist(graph, *args) for key, graph in graphs.items()
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

    def get_all_systs(self, graph, groups, systs, chan):
        out_hists = {}
        for syst in systs:
            for group, members in groups.items():
                if not syst.good_syst(group, self.year, chan):
                    continue
                if syst.name not in out_hists:
                    out_hists[syst.name] = {}
                if syst.syst_type == 'lnN':
                    out_hists[syst.name][group] = syst.get_lnN_error(self.year, chan)
                else:
                    out_hists[syst.name][group] = self.get_syst_hist(graph, members, syst.name)
        return out_hists
