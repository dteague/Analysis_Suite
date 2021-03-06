#!/usr/bin/env python3

import numpy as np
import pandas as pd
import logging
from pandas.api.types import is_numeric_dtype
from .vargetter import VarGetter
from analysis_suite.commons.configs import setup_pandas
import awkward1 as ak
from pathlib import Path
import uproot4
import uproot as upwrite
import json
import re

from sklearn.model_selection import train_test_split

class DataProcessor:
    def __init__(self, use_vars, groupDict, systName="Nominal"):
        """Constructor method
        """
        self.split_ratio = 1/3.
        self.max_events = 2000
        self.max_events_scaled = self.max_events/self.split_ratio
        self.sample_name_map = dict()

        self.group_dict = groupDict.copy()
        self.train_groups = set(sum(self.group_dict.values(), []))
        self.systName = systName

        self.use_vars = use_vars
        self._include_vars = list(use_vars.keys())
        self._drop_vars = ["groupName", "scale_factor"]
        self._all_vars = self._include_vars + self._drop_vars

    def get_final_dict(self, directory):
        arr_dict = dict()
        path = Path(directory)
        root_files = path.rglob("*.root") if path.is_dir() else [path]
        for root_file in root_files:
            groups = list()
            syst = 0
            with uproot4.open(root_file) as f:
                groups = [key.strip(";1") for key in f.keys() if "/" not in key]
                for i, syst_tnamed in enumerate(f[groups[0]]["Systematics"]):
                    if syst_tnamed.member("fName") == self.systName:
                        syst = i
                        break
            for group in groups:
                if group not in arr_dict:
                    arr_dict[group] = VarGetter(root_file, group, syst)
                else:
                    arr_dict[group] += VarGetter(root_file, group, syst)
        return arr_dict

    def process_year(self, infile, outdir):
        train_set, test_set = setup_pandas(self.use_vars, self._all_vars)

        # Process input file
        arr_dict = self.get_final_dict(infile)
        allGroups = set(arr_dict.keys())
        self.group_dict["NotTrained"] = list(allGroups-self.train_groups)
        classID_dict = {"Signal": 1, "NotTrained": 0, "Background": 0}

        for group, samples in self.group_dict.items():
            class_id = classID_dict[group]
            for sample in samples:
                if sample not in arr_dict:
                    logging.warning(f'Could not found sample {sample}')
                    continue
                if not len(arr_dict[sample]):
                    logging.warning(f'Sample {sample} has no events in it!')
                    continue

                if sample not in self.sample_name_map:
                    self.sample_name_map[sample] = len(self.sample_name_map)

                arr = arr_dict[sample]

                df_dict = {varname: func.apply(arr) for varname, func in self.use_vars.items()}
                df_dict["scale_factor"] = ak.to_numpy(arr.scale)
                df = pd.DataFrame.from_dict(df_dict)
                df["classID"] = class_id
                df["groupName"] = self.sample_name_map[sample]

                if group == "NotTrained" or len(arr) < 10:
                    test_set = pd.concat([df.reset_index(drop=True), test_set], sort=True)
                    continue
                
                split_ratio = self.split_ratio if len(df) < self.max_events_scaled \
                    else self.max_events
                train, test = train_test_split(df, train_size=split_ratio,
                                               random_state=12345)
                test["scale_factor"] *= len(df)/len(test)
                train["scale_factor"] *= len(df)/len(train)

                test_set = pd.concat([test.reset_index(drop=True), test_set], sort=True)
                train_set = pd.concat([train.reset_index(drop=True), train_set], sort=True)

        self._write_out(outdir / f'test_{self.systName}.root', test_set)
        self._write_out(outdir / f'train_{self.systName}.root', train_set)

    def _write_out(self, outfile, workSet):
        """**Write out pandas file as a compressed pickle file

        Args:
          outfile(string): Name of file to write
          workSet(pandas.DataFrame): DataFrame of variables to write out
          prediction(pandas.DataFrame): DataFrame of BDT predictions

        """
        workSet["groupName"] = workSet["groupName"].astype("int")
        keepList = [key for key in workSet.columns if is_numeric_dtype(workSet[key])]
        branches = {key: np.int32 if key[0] == "N" else  np.float32 for key in keepList}
        with upwrite.recreate(outfile) as f:
            f["sample_map"] = json.dumps(self.sample_name_map)
            for group in self.sample_name_map.keys():
                groupNum = self.sample_name_map[group]
                groupSet = workSet[workSet.groupName == groupNum][keepList]
                if len(groupSet) == 0:
                    continue
                f[group] = upwrite.newtree(branches)
                f[group].extend(groupSet.to_dict('list'))
