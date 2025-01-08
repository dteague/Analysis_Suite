#!/usr/bin/env python3

"""
.. module:: XGBoostMaker
   :synopsis: Takes in ROOT file to run a BDT training over it using XGBoost
.. moduleauthor:: Dylan Teague
"""
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
from pathlib import Path
from random import randint
import uproot
import operator
import logging
from matplotlib import colors as clr

from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.model_selection import train_test_split

from analysis_suite.plotting.utils import likelihood_sig
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.plot_utils import plot, plot_colorbar, cms_label

pd.options.mode.chained_assignment = None

formatter = {'extra_format': 'pdf',}

def generic_plot_setup(ax, year, bins=None):
    if bins is None:
        ax.set_xlim(0., 1.)
    else:
        ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel("$disc_{BDT}$")
    ax.set_ylabel("A.U.")
    ax.legend()
    cms_label(ax, year=year)

def setup_pandas(all_vars):
    df_set = pd.DataFrame(columns = all_vars)
    for key in all_vars:
        if key[0] == "N" or key == "classID":
            dtype = 'int'
        else:
            dtype = 'float'
        df_set[key] = df_set[key].astype(dtype)
    return df_set

class MLHolder:
    """Wrapper for XGBoost training. Takes an uproot input, a list of
    groups to do a multiclass training, as well as a cut string if
    needed and trains the data. After it is done, the results can be
    outputed to be piped into MVAPlotter

    Args:
      split_ratio(float): Ratio of test events for train test splitting
      group_names(list): List of the names of the different groups
      pred_train(dict): Dictionary of group name to BDT associated with it for train set
      pred_test(dict): Dictionary of group name to BDT associated with it for test set
      train_set(pandas.DataFrame): DataFrame of the training events
      test_set(pandas.DataFrame): DataFrame of the testing events
      cuts(list): List of ROOT style cuts to apply
      param(dict): Variables used in the training

    """
    def __init__(self, use_vars, sig, all_samples, region="signal", systName="Nominal", **kwargs):
        """Constructor method
        """
        # Group information
        self.signal = sig
        self.nonTrained = kwargs.get("nonTrained", ['nonprompt', 'charge_flip', 'data'])
        self.onlyTrained = kwargs.get("onlyTrained", ['nonprompt_mc'])
        self.samples = sorted(all_samples)
        self.sample_map = {val: i for i, val in enumerate(self.samples)}

        self.outfile_info = f'{systName}_{region}'
        nonTrain_vars = ["scale_factor", "split_weight"]
        derived_vars = ["classID", "sampleName", "train_weight"]
        self.use_vars = use_vars
        self._file_vars = use_vars + nonTrain_vars
        self._drop_vars = nonTrain_vars + derived_vars
        self.all_vars = self._file_vars + derived_vars

        self.split_ratio = 0.3
        self.validation_ratio = 0.15
        # self.max_events = 1000 #*5/2
        # self.max_train_events = self.max_events*self.split_ratio
        self.max_train_events = 300000

        self.min_train_events = 100
        self.random_state = 12345#randint(0, 2**32-1)#12345

        self.train_set = setup_pandas(self.all_vars)
        self.validation_set = setup_pandas(self.all_vars)
        self.test_sets = dict()
        self.test_weights = dict()

        self.pred_test = dict()
        self.pred_train = dict()
        self.pred_validation = dict()
        self.cuts = list()
        self.best_iter = 0
        self.minmax = [1, 0]
        self.outdir = None

        # self.group_names = {
        #     'ttt': ['tttw_nlo', 'tttj_nlo', 'tttw_m', 'tttw_p', 'tttj_m', 'tttj_p'],
        #     'ttX': ['ttz', 'ttw', 'tth', 'ttz_m1-10'],
        #     'dd': ['nonprompt', 'charge_flip']+list(self.group_dict["OnlyTrain"]),
        #     'tttt': ['tttt']
        # }
        # used_groups = [i for sub in self.group_names.values() for i in sub]
        # if 'Background' in self.group_dict:
        #     self.group_names['other'] =  [key for key in self.group_dict['Background'] if key not in used_groups]


    def __bool__(self):
        return bool(self.test_sets)

    def set_outdir(self, outdir):
        outdir.mkdir(exist_ok=True, parents=True)
        self.outdir = outdir

    def should_train(self, df, sample):
        enough_events = len(df)*self.split_ratio > self.min_train_events
        trainable_class = sample not in self.nonTrained
        return enough_events and trainable_class

    def read_in_files(self, indir, year, typ='test'):
        test_set= setup_pandas(self.all_vars)
        test_weights = dict()

        with uproot.open(indir/year/f'{typ}_{self.outfile_info}.root') as f:
            for df, sample, weights in self.read_in_dataframe(f):
                test_set = pd.concat([df, test_set], ignore_index=True)
                test_weights[sample] = weights
        self.test_sets[year] = test_set
        self.test_weights[year] = test_weights


    def read_in_train_files(self, indir):
        # Training files
        with uproot.open(indir/f'train_{self.outfile_info}.root') as f:
            for df, sample, _ in self.read_in_dataframe(f):
                self.train_set = pd.concat([df, self.train_set], ignore_index=True)
        # Validation Files
        with uproot.open(indir/f'validation_{self.outfile_info}.root') as f:
            for df, sample, _ in self.read_in_dataframe(f):
                self.validation_set = pd.concat([df, self.validation_set], ignore_index=True)


    def read_in_dataframe(self, root_file, create=False):
        for sample in self.samples:
            classID = 1 if sample in self.signal else 0
            if sample not in root_file:
                logging.debug(f"{sample} not found")
                continue
            # print(sample, "loaded")
            df = root_file[sample].arrays(library="pd")
            mask = self._cut_mask(df)
            if np.count_nonzero(mask) == 0:
                continue
            if create:
                df = df[self._file_vars][mask]
            else:
                df = df[self.all_vars][mask]
            logging.debug(f'{sample}, {len(df)}, {sum(df.scale_factor)}, {sum(df.scale_factor)/len(df)}')
            df.loc[:, "classID"] = classID
            df.loc[:, "sampleName"] = self.sample_map[sample]
            df.loc[:, "train_weight"] = sum(df.scale_factor) / len(df)
            df.loc[:, "split_weight"] = np.ones(len(df))
            if "weights" not in root_file:
                all_weights = None
            else:
                all_weights = root_file[f'weights/{sample}'].arrays(library='pd')[mask]
            yield df, sample, all_weights

    def setup_weights(self, filename):
        self.total_yield = 0.
        with uproot.open(filename) as f:
            # Setup estimation of number of events (do based on ttW events)
            sample_wgt = f['ttw']['scale_factor'].array()
            self.total_trained_evt = len(sample_wgt)/0.2*self.split_ratio
            # Setup total yield amount
            for sample in self.samples:
                if sample not in f or sample in self.nonTrained:
                    continue
                sample_wgt = f[sample]["scale_factor"].array()
                if len(sample_wgt) < self.min_train_events/self.split_ratio:
                    continue
                self.total_yield += np.sum(sample_wgt)


    def setup_year(self, directory, year="2018", split=True):
        """**Fill the dataframes with all info in the input files**

        This grabs all the variable information about each sample,
        does some preliminary weighting and splits the data into the
        test and train set (based on `self.split_ratio`)

        Args:
            directory(string): Path to directory where root files are kept
        """
        print(year)
        self.test_sets[year] = setup_pandas(self.all_vars)
        self.test_weights[year] = dict()
        infile = directory / year / f'processed_{self.outfile_info}.root'
        self.setup_weights(infile)
        with uproot.open(infile) as f:
            for df, sample, weights in self.read_in_dataframe(f):
                test, train, validation = self.setup_split(df, sample, split)
                total_scale = np.sum(df.scale_factor)
                exp_train_evt = self.total_trained_evt*(total_scale/self.total_yield)
                print(sample, len(df), total_scale, f'{int(exp_train_evt)}')
                print("    ", len(test), len(train), len(validation), f'{1.-len(test)/(len(df)+1e-5):0.3f}')
                if not test.empty:
                    self.test_weights[year][sample] = pd.DataFrame()
                    for key, col in weights.items():
                        if key == 'index':
                            continue
                        self.test_weights[year][sample].insert(
                            0, key, col.loc[test.index]*np.sum(col)/np.sum(col[test.index]))
                    self.test_sets[year] = pd.concat([test, self.test_sets[year]], ignore_index=True)
                if not train.empty:
                    self.train_set = pd.concat([train, self.train_set], ignore_index=True)
                if not validation.empty:
                    self.validation_set = pd.concat([validation, self.validation_set], ignore_index=True)
        print()

    def setup_split(self, df, sample, split=True):
        train = setup_pandas(self.all_vars)
        test = setup_pandas(self.all_vars)
        validation = setup_pandas(self.all_vars)
        total_scale = np.sum(df.scale_factor)

        if not split or not self.should_train(df, sample):
            return df, train, validation

        exp_train_evt = self.total_trained_evt*(total_scale/self.total_yield)
        split_ratio = self.split_ratio

        # 3top and 4top
        if "ttt" in sample:
            exp_train_evt = split_ratio*len(df)

        # Only train, i.e. use all the events possible
        if sample in self.onlyTrained:
            print(exp_train_evt)
            if exp_train_evt <= self.min_train_events:
                pass
            elif exp_train_evt > len(df):
                df.loc[:, "split_weight"] *= exp_train_evt/len(df)
                train = df
            else:
                _, train = self.split(df, exp_train_evt/len(df), total_scale)
        # Do Train
        else:
            if exp_train_evt <= len(df)*split_ratio and exp_train_evt > self.min_train_events:
                split_ratio = exp_train_evt/len(df)
            elif exp_train_evt <= self.min_train_events:
                split_ratio = self.min_train_events/len(df)
                df.loc[:, "split_weight"] *= exp_train_evt/self.min_train_events
            else:
                df.loc[:, "split_weight"] *= exp_train_evt/(split_ratio*len(df))
            test, train = self.split(df, split_ratio, total_scale)

        if self.validation_ratio > 0. and len(train)*self.validation_ratio > 1:
            train, validation = self.split(train, self.validation_ratio, total_scale)

        return test, train, validation


    def split(self, workset, split_ratio, target_scale):
        train, test = train_test_split(workset, train_size=split_ratio, random_state=self.random_state)
        test.loc[:, ["scale_factor", "train_weight"]] *= target_scale/np.sum(test.scale_factor)
        train.loc[:, ["scale_factor", "train_weight"]] *= target_scale/np.sum(train.scale_factor)
        return test, train


    def apply_model(self, year, skip_train=False):
        def get_pred(use_set, directory):
            return self.predict(use_set, directory).T[1]

        self.pred_test[year] = get_pred(self.test_sets[year], self.outdir)
        if not skip_train:
            self.pred_train = get_pred(self.train_set, self.outdir)
            self.pred_validation = get_pred(self.validation_set, self.outdir)


    def get_auc(self):
        is_sig = np.array([], dtype='bool')
        pred = np.array([])
        scale = np.array([])
        for year, df in self.test_sets.items():
            is_sig= np.concatenate([is_sig, df.classID == 1], dtype='bool')
            scale = np.concatenate([scale, df.scale_factor])
            pred = np.concatenate([pred, self.pred_test[year]])

        nbins = 100
        bins = np.linspace(0, 1+1/nbins, nbins+2)

        s_hist = np.histogram(pred[is_sig], bins, weights=scale[is_sig])[0]
        b_hist = np.histogram(pred[~is_sig], bins, weights=scale[~is_sig])[0]
        tp = np.cumsum(s_hist[::-1])/(np.sum(s_hist)+1e-5)
        fp = np.cumsum(b_hist[::-1])/(np.sum(b_hist)+1e-5)

        delta = fp[1:]-fp[:-1]
        trap = (tp[1:]+tp[:-1])/2
        auc = np.sum(delta*trap)
        return auc

    def get_fom(self):
        is_sig = np.array([], dtype='bool')
        pred = np.array([])
        scale = np.array([])
        for year, df in self.test_sets.items():
            is_sig= np.concatenate([is_sig, df.classID == 1], dtype='bool')
            scale = np.concatenate([scale, df.scale_factor])
            pred = np.concatenate([pred, self.pred_test[year]])

        nbins = 20
        bins = np.linspace(0, 1, nbins+1)

        s = np.histogram(pred[is_sig], bins, weights=scale[is_sig])[0]
        b = np.histogram(pred[~is_sig], bins, weights=scale[~is_sig])[0]

        with np.errstate(invalid='ignore'):
            fom = np.sqrt(2*np.sum((s + b)*np.log(1 + s/(b+1e-5)) - s))
        return fom

    # Private Functions
    def _cut_mask(self, frame):
        """**Reduce frame using root style cut string**

        Args:
          frame(pandas.DataFrame): DataFrame to cut on

        """
        mask = np.ones(len(frame), dtype=bool)
        for cut in self.cuts:
            if cut.find("<") != -1:
                cutter= (operator.lt, *cut.split("<"))
            elif cut.find(">") != -1:
                cutter= (operator.gt, *cut.split(">"))
            elif cut.find("==") != -1:
                cutter= (operator.eq, *cut.split("=="))
            else:
                raise Exception(f'{cut} is not formatted correctly')
            mask *= cutter[0](frame[cutter[1]], float(cutter[2]))
        return mask

    def add_cut(self, cut_string):
        self.cuts.append(cut_string)


    def output_files(self):
        for year in self.test_sets:
            (self.outdir/year).mkdir(exist_ok=True, parents=True)
            self._output(self.test_sets[year], self.outdir / year / f"test_{self.outfile_info}.root",
                         self.test_weights[year])
        self._output(self.train_set, self.outdir / f"train_{self.outfile_info}.root")
        self._output(self.validation_set, self.outdir / f"validation_{self.outfile_info}.root")

    def _output(self, workSet, outfile, weights=None):
        """**Write out pandas file as a compressed pickle file

        Args:
          workSet(pandas.DataFrame): DataFrame of variables to write out
          outfile(string): Name of file to write
        """
        keepList = [key for key in workSet.columns if is_numeric_dtype(workSet[key])]
        samples = np.unique(workSet.sampleName)
        with uproot.recreate(outfile) as f:
            for sample, value in self.sample_map.items():
                if value not in samples:
                    continue
                f[sample] = workSet[workSet.sampleName == value][keepList]
                if weights is not None:
                    f[f"weights/{sample}"] = weights[sample]

    def output_bdt(self, outfile, sig_out='Signal'):
        with uproot.recreate(outfile) as f:
            for sample, value in self.sample_map.items():
                # Test set
                for year, pred in self.pred_test.items():
                    mask = self.test_sets[year].sampleName == value
                    if np.count_nonzero(mask) > 0:
                        f[f'{year}/{sample}'] = {"BDT": pred[mask]}
                # Training set
                mask = self.train_set.sampleName == value
                if np.count_nonzero(mask) > 0:
                    f[f'train/{sample}'] = {"BDT": self.pred_train[mask]}
                # Validation set
                mask = self.validation_set.sampleName == value
                if np.count_nonzero(mask) > 0:
                    f[f'validation/{sample}'] = {"BDT": self.pred_validation[mask]}
