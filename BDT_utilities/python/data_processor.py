#!/usr/bin/env python3

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import xgboost as xgb
from analysis_suite.commons import VarGetter
import awkward1 as ak
from pathlib import Path
import uproot4
import uproot as upwrite
import json

from sklearn.model_selection import train_test_split

class DataProcessor:
    def __init__(self, use_vars, groupDict):
        """Constructor method
        """
        self.split_ratio = 1/3.
        self.max_events = 300000000
        self.sample_name_map = dict()

        self.group_dict = groupDict
        self.train_groups = set(sum(self.group_dict.values(), []))

        self.use_vars = use_vars
        self._include_vars = list(use_vars.keys())
        self._drop_vars = ["groupName", "scale_factor"]
        self._all_vars = self._include_vars + self._drop_vars

        self.train_all = pd.DataFrame(columns=self._all_vars)
        for key, func in self.use_vars.items():
            dtype = "int" if "num" in func else 'float'
            self.train_all[key] = self.train_all[key].astype(dtype)

    def get_final_dict(self, directory):
        arr_dict = dict()
        path = Path(directory)
        root_files = path.rglob("*.root") if path.is_dir() else [path]
        for root_file in root_files:
            groups = list()
            with uproot4.open(root_file) as f:
                groups = [key.strip(";1") for key in f.keys() if "/" not in key]
            for group in groups:
                if group not in arr_dict:
                    arr_dict[group] = VarGetter(root_file, group)
                else:
                    arr_dict[group] += VarGetter(root_file, group)

        return arr_dict

    def setup_test_train(self, arr, group, sample, class_id, noTrain):
        df_dict = dict()
        for varname, func in self.use_vars.items():
            df_dict[varname] = ak.to_numpy(eval("arr.{}".format(func)))
        df_dict["scale_factor"] = ak.to_numpy(arr.scale)

        df = pd.DataFrame.from_dict(df_dict)
        # df = self._cut_frame(df)
        df["classID"] = class_id

        print(df["groupName"])
        if noTrain:
            return df, None

        split_ratio = self.split_ratio if len(df) < self.max_events else self.max_events
        return train_test_split(df, train_size=split_ratio, random_state=12345)



    def process_year(self, year, outdir):
        # Setup dataframes to be used
        train_set = pd.DataFrame(columns=self._all_vars)
        test_set = pd.DataFrame(columns=self._all_vars)
        for key, func in self.use_vars.items():
            dtype = "int" if "num" in func else 'float'
            train_set[key] = train_set[key].astype(dtype)
            test_set[key] = test_set[key].astype(dtype)

        # Process input file
        arr_dict = self.get_final_dict("result_{}.root".format(year))
        allGroups = set(arr_dict.keys())
        self.group_dict["NotTrained"] = list(allGroups-self.train_groups)
        classID_dict = {"Signal": 1, "NotTrained": -1, "Background": 0}
        
        for group, samples in self.group_dict.items():
            class_id = classID_dict[group]
            for sample in samples:
                noTrain = False
                if sample not in arr_dict:
                    print("Could not found sample {}".format(sample))
                    continue
                if len(arr_dict[sample].scale) == 0:
                    print("Sample {} has no events in it!".format(sample))
                    continue

                noTrain = class_id < 0 or len(arr_dict[sample].scale) < 10

                df_dict = dict()
                arr = arr_dict[sample]
                for varname, func in self.use_vars.items():
                    df_dict[varname] = ak.to_numpy(eval("arr.{}".format(func)))
                df_dict["scale_factor"] = ak.to_numpy(arr.scale)

                df = pd.DataFrame.from_dict(df_dict)
                # df = self._cut_frame(df)
                df["classID"] = class_id
                if sample not in self.sample_name_map:
                    self.sample_name_map[sample] = len(self.sample_name_map)
                df["groupName"] = self.sample_name_map[sample]
                
                if noTrain:
                    test_set = pd.concat([df.reset_index(drop=True), test_set], sort=True)
                    continue
                
                split_ratio = self.split_ratio if len(df) < self.max_events \
                    else self.max_events
                train, test = train_test_split(df, train_size=split_ratio,
                                               random_state=12345)
                test_set = pd.concat([test.reset_index(drop=True), test_set], sort=True)
                train_set = pd.concat([train.reset_index(drop=True), train_set], sort=True)

        self.train_all = pd.concat([train_set.reset_index(drop=True), self.train_all], sort=True)
        self._write_out("{}/{}/test.root".format(outdir, year), test_set)
        self._write_out("{}/{}/train.root".format(outdir, year), train_set)

    def write_train(self, outdir):
        self._write_out("{}/train.root".format(outdir), self.train_all)


    def _write_out(self, outfile, workSet):
        """**Write out pandas file as a compressed pickle file

        Args:
          outfile(string): Name of file to write
          workSet(pandas.DataFrame): DataFrame of variables to write out
          prediction(pandas.DataFrame): DataFrame of BDT predictions

        """
        workSet["groupName"] = workSet["groupName"].astype("int")
        keepList = [key for key in workSet.columns if is_numeric_dtype(workSet[key])]
        branches = {key: workSet[key].dtype for key in keepList}
        with upwrite.recreate(outfile) as f:
            f["sample_map"] = json.dumps(self.sample_name_map)
            for group in self.sample_name_map.keys():
                groupNum = self.sample_name_map[group]
                groupSet = workSet[workSet.groupName == groupNum][keepList]
                f[group] = upwrite.newtree(branches)
                f[group].extend(groupSet.to_dict('list'))
