#!/usr/bin/env python3
import numpy as np
import re
from dataclasses import dataclass, field
from typing import  Callable
from pathlib import Path

from analysis_suite.data.FileInfo import info as finfo
from analysis_suite.data.PlotGroups import info as ginfo

class GroupInfo:
    def __init__(self, group2color=None, **kwargs):
        self.group2color = group2color if group2color is not None else {}
        self.group2MemberMap = self.get_memberMap()

    def get_legend_name(self, group):
        return ginfo[group]["Name"]

    def get_color(self, group):
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
            final[key] = members
        return final

    def get_members(self, group):
        return self.group2MemberMap[group]

    def setup_groups(self, groups=None):
        if groups is None:
            groups = self.group2color.keys()
        return dict(filter(lambda x: x[0] in groups, self.group2MemberMap.items()))

    def setup_members(self, groups=None):
        groups = self.setup_groups(groups)
        if not groups.values():
            return list()
        return np.concatenate(list(groups.values()))

    
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
    region: str
    cut : Callable[[object], bool] = None
    branches: list = None
    changes: dict = field(default_factory=dict)
    ignores: dict = field(default_factory=dict)

    def get_file(self, **kwargs):
        return Path(str(self.filename).format(**kwargs))

    def get_filename(self, year, workdir=None):
        if workdir is None:
            path = Path(str(self.filename).format(year=year, workdir=""))
            workdir = max([int(d.name) for d in path.glob("*") if d.name.isnumeric()])
            print(f"Getting from workdir {workdir}")
        return Path(str(self.filename).format(year=year, workdir=workdir))

    def add_change(self, tree, changes):
        self.changes[tree] = changes

    def add_ignore(self, tree, ignores):
        if isinstance(ignores, str):
            ignores = [ignores]
        self.ignores[tree] = ignores

    def ignore(self, tree, group):
        if tree in self.ignores:
            return group in self.ignores[tree]
        return False

    def get_change(self, tree, member):
        if tree in self.changes and member in self.changes[tree]:
            return self.changes[tree][member]
        return member

    def apply_cut(self, vg, *args):
        def cut_vg(vg, cut):
            vg.cut(self.cut(vg))
            for name, part in vg.parts.items():
                part.reset()

        if self.cut is None:
            return
        elif isinstance(self.cut, list):
            for cut in self.cut:
                cut_vg(cut)
        else:
            cut_vg(self.cut)





    def setup_branches(self, vg):
        if self.branches is None:
            return
        if isinstance(self.branches, list):
            for func in self.branches:
                func(vg)
        else:
            self.branches(vg)
