#!/usr/bin/env python3
import numpy as np
import os
import glob
import imp
import subprocess
from pathlib import Path

class BasicInfo:
    def __init__(self, analysis="", selection="", lumi=140000,
               mcPath = "data/FileInfo/montecarlo/montecarlo_2016.py", **kwargs):
        self.analysis = analysis
        self.selection = selection
        self.lumi = lumi
        self.mcInfo = self.readAllInfo(mcPath)

    def readAllInfo(self, file_path):
        info = {}
        for info_file in glob.glob(file_path):
            file_info = self.readInfo(info_file)
            if file_info:
                info.update(file_info)
        return info

    def readInfo(self, file_path):
        if ".py" not in file_path[-3:] and ".json" not in file_path[-5:]:
            if os.path.isfile(file_path + ".py"):
                file_path = file_path + ".py"
            elif os.path.isfile(file_path + ".json"):
                file_path = file_path + ".json"
            else:
                return
        if ".py" in file_path[-3:]:
            file_info = imp.load_source("info_file", file_path)
            info = file_info.info
        else:
            info = self.readJson(file_path)
        return info

    def readJson(self, json_file_name):
        json_info = {}
        with open(json_file_name) as json_file:
            try:
                json_info = json.load(json_file)
            except ValueError as err:
                print("Error reading JSON file {}. The error message was:"
                      .format(json_file_name))
                print(err)
        return json_info

    def get_xsec(self, group):
        scale = self.mcInfo[group]['cross_section']
        if "kfactor" in self.mcInfo[group]:
            scale *= self.mcInfo[group]["kfactor"]
        return scale

class PlotInfo(BasicInfo):
    def __init__(self, plotInfo, **kwargs):
        super().__init__(**kwargs)
        self.plotSpecs = self.readAllInfo(plotInfo)

    def __getattr__(self, name):
        return {h: (vals[name] if name in vals else None) for h, vals in self.plotSpecs.items()}

    def __getitem__(self, key):
        return self.plotSpecs[key]

    def get_hists(self):
        return self.plotSpecs.keys()

    def get_binning(self, histname):
        import boost_histogram as bh
        bins = np.array(self.Binning[histname], dtype=float)
        if self.Discrete[histname]:
            bins[1:] = bins[1:] - 0.5
        return bh.axis.Regular(int(bins[0]), *bins[1:])


class FileInfo(BasicInfo):
    def __init__(self, group2color={}, **kwargs):
        super().__init__(**kwargs)
        self.group2color = group2color
        self.fileInfo = self.readAllInfo(
            "data/FileInfo/{}/{}.py".format(self.analysis,self.selection))
        self.groupInfo = self.readAllInfo(
            "data/PlotGroups/{}.py".format(self.analysis))
        if group2color:
            self.group2MemberMap = {key: item["Members"] for key, item in
                                    self.groupInfo.items() if key in self.group2color}
        else:
            self.group2MemberMap = {key: item["Members"]
                                    for key, item in self.groupInfo.items()}

    def get_file_with_size(self, group_list=None):
        return_dict = dict()
        for key, item in self.fileInfo.items():
            files = item["file_path"]
            if isinstance(files, str):
                dirs = Path(files)
                return_dict[key] = [(str(f), f.stat().st_size) for f in dirs.glob("**/*.root")]
            else:
                return_dict[key] = files

        return return_dict
        # if group_list is None:
        #     return {key: item["file_path"] for key, item in self.fileInfo.items()}
        # else:
        #     return_dict = dict()
        #     for group in group_list:
        #         if group in self.group2MemberMap:
        #             return_dict.update({sample: self.fileInfo[sample]["file_path"]
        #                                 for sample in self.group2MemberMap[group]})
        #         else:
        #             return_dict[group] = self.fileInfo[group]["file_path"]
        #     return return_dict


    def get_file_dict(self, group_list=None):
        if group_list is None:
            return {key: item["file_path"] for key, item in self.fileInfo.items()}
        else:
            return_dict = dict()
            for group in group_list:
                if group in self.group2MemberMap:
                    return_dict.update({sample: self.fileInfo[sample]["file_path"]
                                        for sample in self.group2MemberMap[group]})
                else:
                    return_dict[group] = self.fileInfo[group]["file_path"]
            return return_dict

    def get_file_dict_with_xsec(self, group_list=None):
        out_dict = dict()
        print(group_list)
        print(self.get_file_dict(group_list))
        for group, files in self.get_file_dict(group_list).items():
            out_dict[group] = (files, self.get_xsec(group))
        return out_dict
    
    def get_xrootd(self, dirpath):
        return subprocess.run(["hdfs", "dfs", "-find", dirpath, "-name", "*root"],
                       capture_output=True).stdout.split()

    def get_color(self, group):
        return self.group2color[group]

    def get_legend_name(self, group):
        return self.groupInfo[group]['Name']