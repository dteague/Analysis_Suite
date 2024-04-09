#!/usr/bin/env python3
import numpy as np
import re
from dataclasses import dataclass, field
from typing import  Callable
from pathlib import Path
import logging

from analysis_suite.data.FileInfo import info as finfo
from analysis_suite.data.PlotGroups import info as ginfo

class GroupInfo:
    def __init__(self, group2color=None, **kwargs):
        self.group2color = group2color if group2color is not None else {}
        self.group2MemberMap = self.get_memberMap()

    def get_groups(self):
        return list(self.group2color.keys())

    def get_legend_name(self, group):
        if group in ginfo:
            return ginfo[group]["Name"]
        else:
            return group

    def get_color(self, group):
        if group not in self.group2color:
            logging.warning("No color info given, default to black")
            return 'black'
        return self.group2color[group]

    def get_memberMap(self):
        final = dict()
        for key, info in ginfo.items():
            members = info["Members"]
            if "Composite" in info and info["Composite"]:
                tmpMembers = list()
                for mem in members:
                    if mem in ginfo:
                        tmpMembers += ginfo[mem]["Members"]
                    else:
                        tmpMembers.append(mem)
                members = tmpMembers
            if "DataDriven" in info and info['DataDriven']:
                members = [key]
            final[key] = members
        return final

    def get_members(self, group):
        return self.group2MemberMap[group]

    def setup_groups(self, groups=None):
        if groups is None:
            groups = self.group2color.keys()
        return {group : self.group2MemberMap[group] for group in groups}

    def setup_members(self, groups=None):
        groups = self.setup_groups(groups)
        if not groups.values():
            return list()
        return np.concatenate(list(groups.values()))

    def is_data_driven(self, group):
        if "DataDriven" not in ginfo[group]:
            return False
        return ginfo[group]["DataDriven"]
    
class FileInfo:
    def __init__(self):
        self.dasNames = {key: info["DAS"] for key, info in finfo.items()}

    def get_group(self, splitname):
        if isinstance(splitname, str) and splitname in self.dasNames:
            return self.dasNames[splitname]
        elif self.is_data(splitname):
            return 'data'

        sample_name = next(filter(lambda x: "13TeV" in x, splitname), None)
        for name, reName in self.dasNames.items():
            if re.match(reName, sample_name) is not None:
                return name
        return None

    def get_info(self, alias):
        return finfo[alias]

    def get_xsec(self, group):
        if self.is_data(group):
            return 1.
        info = self.get_info(group)
        scale = info['cross_section']
        if "kfactor" in info:
            scale *= info["kfactor"]
        return scale

    def is_data(self, group):
        if group == 'data':
            return True
        elif isinstance(group, str) and  group not in finfo:
            return True
        for split in map(str.lower, group):
            if "data" in split:
                return True
        return False

fileInfo = FileInfo()

@dataclass
class NtupleInfo:
    filename: str
    trees: list
    color_by_group: dict
    cut : Callable[[object], bool] = None
    branches: list = None
    changes: dict = field(default_factory=dict)
    ignore: dict = field(default_factory=dict)
    part_cut: list = None
    tree_groups: dict = field(default_factory=dict)
    group_trees: dict = field(default_factory=dict)

    def get_info(self, remove=None, add=None):
        color_by_group = self.color_by_group
        if not isinstance(remove, list):
            remove = [remove]
        if remove is not None:
            color_by_group = {group: color for group, color in color_by_group.items() if group not in remove}
        if add is not None:
            color_by_group = {**clr_by_group, **add}
        return GroupInfo(color_by_group)

    def get_file(self, **kwargs):
        return Path(str(self.filename).format(**kwargs))

    def get_filename(self, year, workdir=None):
        if workdir is None:
            path = Path(str(self.filename).format(year=year, workdir=""))
            workdir = max([int(d.name) for d in path.glob("*") if d.name.isnumeric()])
            logging.info(f"Getting from workdir {workdir}")
        return Path(str(self.filename).format(year=year, workdir=workdir))

    def set_groups_trees(self, trees, groups):
        if isinstance(trees, str):
            trees = [trees]
        if isinstance(groups, str):
            groups = [groups]
        for tree in trees:
            self.tree_groups[tree] = groups
        for group in groups:
            self.group_trees[group] = tree

    def pass_group(self, tree, group):
        if tree in self.tree_groups:
            return group in self.tree_groups[tree]
        elif group in self.group_trees:
            return tree in self.group_trees[group]
        else:
            return True

    def apply_part_cut(self, vg):
        if self.part_cut is None:
            return
        for cut in self.part_cut:
            vg[cut[0]].mask_part(cut[1], cut[2])


    def apply_cut(self, vg, *args):
        def cut_vg(vg, cut):
            vg.cut(cut)
            for name, part in vg.parts.items():
                part.reset()

        if self.cut is None:
            return
        elif isinstance(self.cut, list):
            for cut in self.cut:
                cut_vg(vg, cut)
        else:
            cut_vg(vg, self.cut)

    def setup_branches(self, vg):
        if self.branches is None:
            return
        if isinstance(self.branches, list):
            for func in self.branches:
                func(vg)
        else:
            self.branches(vg)
