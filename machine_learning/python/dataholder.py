#!/usr/bin/env python3

"""
.. module:: XGBoostMaker
   :synopsis: Takes in ROOT file to run a BDT training over it using XGBoost
.. moduleauthor:: Dylan Teague
"""
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import uproot
import logging
from copy import copy

from sklearn.model_selection import train_test_split
from analysis_suite.data.PlotGroups import info as ginfo

from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import plot, cms_label

pd.options.mode.chained_assignment = None

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
        self.onlyTrained = kwargs.get("onlyTrained", ginfo['nonprompt_mc']['Members'])
        self.samples = sorted(all_samples)
        self.sample_map = {val: i for i, val in enumerate(self.samples)}

        self.outfile_info = f'{systName}_{region}'
        nonTrain_vars = ["scale_factor"]
        derived_vars = ["classID", "sampleName", "train_weight", "split_weight"]
        self.use_vars = use_vars
        self._file_vars = use_vars + nonTrain_vars
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
        self.best_iter = 0
        self.minmax = [1, 0]
        self.outdir = None


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

    def read_in_bdt(self, infile, variable, usevar=False, onlyTest=False):
        if usevar:
            self.usevars.append(variable)

        with uproot.open(infile) as f:
            for year, test in self.test_sets.items():
                test_bdt = np.array([])
                for sample in self.samples:
                    if sample in f[year]:
                        test_bdt = np.concatenate([f[f'{year}/{sample}']['BDT'].array(), test_bdt])
                test.insert(0, variable, test_bdt)

            if onlyTest:
                return
            train_bdt = np.array([])
            validation_bdt = np.array([])
            for sample in self.samples:
                if sample in f['train']:
                    train_bdt = np.concatenate([f[f'train/{sample}']['BDT'].array(), train_bdt])
                if sample in f['validation']:
                    validation_bdt = np.concatenate([f[f'validation/{sample}']['BDT'].array(), validation_bdt])
            self.train_set.insert(0, variable, train_bdt)
            self.validation_set.insert(0, variable, validation_bdt)


    def mask(self, func):
        for year, test in self.test_sets.items():
            pass_mask = func(test, year).to_numpy()
            for sample, df in self.test_weights[year].items():
                sample_mask = pass_mask[test.sampleName == self.sample_map[sample]]
                weight = self.test_weights[year][sample]
                self.test_weights[year][sample] = weight[sample_mask]
            self.test_sets[year] = test[pass_mask]
        if not self.train_set.empty:
            self.train_set = self.train_set[func(self.train_set, '2018')]
            self.validation_set = self.validation_set[func(self.validation_set, '2018')]


    def read_in_dataframe(self, root_file, create=False):
        for sample in self.samples:
            classID = 1 if sample in self.signal else 0
            if sample not in root_file:
                logging.debug(f"{sample} not found")
                continue
            # print(sample, "loaded")
            df = root_file[sample].arrays(library="pd")
            if create:
                df.loc[:, "split_weight"] = np.ones(len(df))
            df.loc[:, "classID"] = classID
            df.loc[:, "sampleName"] = self.sample_map[sample]
            df.loc[:, "train_weight"] = sum(df.scale_factor) / len(df)
            if "weights" not in root_file:
                all_weights = None
            else:
                all_weights = root_file[f'weights/{sample}'].arrays(library='pd')
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
            for df, sample, weights in self.read_in_dataframe(f, create=True):
                test, train, validation = self.setup_split(df, sample, split)
                total_scale = np.sum(df.scale_factor)
                exp_train_evt = self.total_trained_evt*(total_scale/self.total_yield)
                # print(sample, len(df), total_scale, f'{int(exp_train_evt)}')
                # print("    ", len(test), len(train), len(validation), f'{1.-len(test)/(len(df)+1e-5):0.3f}')
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
        # print()

    def setup_split(self, df, sample, split=True):
        train = setup_pandas(self.all_vars)
        test = setup_pandas(self.all_vars)
        validation = setup_pandas(self.all_vars)

        total_scale = np.sum(df.scale_factor)
        exp_train_evt = self.total_trained_evt*(total_scale/self.total_yield)
        split_ratio = self.split_ratio

        # 3top and 4top
        # Kinda hacky solution, but ok for now
        if sample == 'tttt':
            exp_train_evt = split_ratio*len(df)
        elif "ttt" in sample:
            exp_train_evt = split_ratio*len(df)
            # roughly 3x tttw_p/tttw_m than tttj_m and 6x for tttj_p
            # doesn't have to be perfect, just to get roughly alright
            if sample == 'tttj_m':
                df.loc[:, "split_weight"] *= 8./15
            elif sample == 'tttj_p':
                df.loc[:, "split_weight"] *= 4./15
            else:
                df.loc[:, "split_weight"] *= 8./5


        # Only train, i.e. use all the events possible
        if sample in self.onlyTrained:
            if exp_train_evt <= self.min_train_events:
                pass
            elif exp_train_evt > len(df):
                df.loc[:, "split_weight"] *= exp_train_evt/len(df)
                train = df
            else:
                _, train = self.split(df, exp_train_evt/len(df), total_scale)
        elif sample in self.nonTrained:
            print("not trained", sample)
            return df, train, validation
        elif len(df)*self.split_ratio < self.min_train_events:
            print("not enough", sample)
            return df, train, validation
        else:    # Do Train
            if exp_train_evt <= len(df)*split_ratio and exp_train_evt > self.min_train_events:
                split_ratio = exp_train_evt/len(df)
            elif exp_train_evt <= self.min_train_events:
                split_ratio = self.min_train_events/len(df)
                df.loc[:, "split_weight"] *= exp_train_evt/self.min_train_events
            else:
                df.loc[:, "split_weight"] *= exp_train_evt/(split_ratio*len(df))
            test, train = self.split(df, split_ratio, total_scale)

        # Setup validation set
        if self.validation_ratio > 0. and len(train)*self.validation_ratio > 1:
            train, validation = self.split(train, self.validation_ratio, total_scale)

        return test, train, validation


    def split(self, workset, split_ratio, target_scale):
        train, test = train_test_split(workset, train_size=split_ratio, random_state=self.random_state)
        test.loc[:, ["scale_factor", "train_weight"]] *= target_scale/np.sum(test.scale_factor)
        train.loc[:, ["scale_factor", "train_weight"]] *= target_scale/np.sum(train.scale_factor)
        return test, train


    def apply_model(self, years, skip_train=False):
        def get_pred(use_set, directory):
            return self.predict(use_set, directory).T[1]
        for year in years:
            self.pred_test[year] = get_pred(self.test_sets[year], self.outdir)
        if not skip_train:
            self.pred_train = get_pred(self.train_set, self.outdir)
            self.pred_validation = get_pred(self.validation_set, self.outdir)


    def get_hist(self, nbins, year=None, useTrain=False):
        output = {}
        if useTrain:
            workset = self.train_set
            pred = self.pred_train
        else:
            workset = self.test_sets[year]
            pred = self.pred_test[year]
        for sample, value in self.sample_map.items():
            mask = workset.sampleName == value
            if np.count_nonzero(mask) == 0:
                continue
            output[sample] = Histogram(*(nbins, 0, 1), axis_name="BDT")
            output[sample].fill(pred[mask], weight=workset[mask].scale_factor)
        return output

    def get_test_hist(self, bins, hist_year=None):
        is_sig = np.array([], dtype='bool')
        pred = np.array([])
        scale = np.array([])
        for year, df in self.test_sets.items():
            if hist_year is not None and year != hist_year:
                continue
            is_sig= np.concatenate([is_sig, df.classID == 1], dtype='bool')
            scale = np.concatenate([scale, df.scale_factor])
            pred = np.concatenate([pred, self.pred_test[year]])
        s_hist = np.histogram(pred[is_sig], bins, weights=scale[is_sig])[0]
        b_hist = np.histogram(pred[~is_sig], bins, weights=scale[~is_sig])[0]
        return s_hist, b_hist

    def get_train_hist(self, bins):
        is_sig = self.train_set.classID == 1
        pred = self.pred_train
        scale = self.train_set.scale_factor
        s_hist = np.histogram(pred[is_sig], bins, weights=scale[is_sig])[0]
        b_hist = np.histogram(pred[~is_sig], bins, weights=scale[~is_sig])[0]
        return s_hist, b_hist

    def get_auc(self, auc_year=None, get_train=False, return_roc=False):
        nbins = 100
        bins = np.linspace(0, 1+1/nbins, nbins+2)

        if get_train:
            s_hist, b_hist = self.get_train_hist(bins)
        else:
            s_hist, b_hist = self.get_test_hist(bins, auc_year)
        tp = np.cumsum(s_hist[::-1])/(np.sum(s_hist)+1e-5)
        fp = np.cumsum(b_hist[::-1])/(np.sum(b_hist)+1e-5)

        delta = fp[1:]-fp[:-1]
        trap = (tp[1:]+tp[:-1])/2
        auc = np.sum(delta*trap)
        if return_roc:
            return fp, tp, auc
        else:
            return auc

    def get_fom(self, get_train=False, fom_year=None):
        nbins = 20
        bins = np.linspace(0, 1, nbins+1)
        if get_train:
            s, b= self.get_train_hist(bins)
        else:
            s, b= self.get_test_hist(bins, fom_year)

        with np.errstate(invalid='ignore'):
            s = abs(s)
            b = abs(b)
            fom = np.sqrt(2*np.sum((s + b)*np.log(1 + s/(b+1e-5)) - s))
        return fom

    def plot_overtrain(self, plot_year=None):
        nbins = 15
        bins = np.linspace(0, 1, nbins+1)

        test_s, test_b= self.get_train_hist(bins)
        train_s, train_b= self.get_test_hist(bins, plot_year)

        year = 'all' if plot_year is None else plot_year

        kw_hist = {"alpha": 0.3, "hatch": '///', "histtype": "stepfilled"}
        kw_err = {"markersize": 4, "fmt": "o"}
        with plot(self.outdir/f"overtrain_{year}.pdf") as ax:
            ax.hist(x=bins[:-1], bins=bins, weights=test_s/np.sum(test_s), color='r',
                    label="Signal (test)", **kw_hist)
            ax.hist(x=bins[:-1], bins=bins, weights=test_b/np.sum(test_b), color='b',
                    label="Background (test)", **kw_hist)
            ax.errorbar(x=bins[:-1]+1/(2*nbins), xerr=1/(2*nbins), color='r',
                        y=train_s/np.sum(train_s), label="Signal (train)", **kw_err)
            ax.errorbar(x=bins[:-1]+1/(2*nbins), xerr=1/(2*nbins), color='b',
                        y=train_b/np.sum(train_b), label="Background (train)", **kw_err)
            ax.set_xlim(0., 1.)
            ax.set_xlabel("$disc_{BDT}$")
            ax.set_ylabel("A.U.")
            ax.legend()
            cms_label(ax, year=year)

    def plot_roc(self, year):
        fp_test, tp_test, auc_test = self.get_auc(year, return_roc=True)
        fp_train, tp_train, auc_train= self.get_auc(year, get_train=True, return_roc=True)

        with plot(self.outdir/f"roc_{year}.pdf") as ax:
            ax.plot(fp_test, tp_test, linewidth=3, label=f"Test set: AUC={auc_test:0.3f}")
            ax.plot(fp_train, tp_train, linewidth=3, label=f"Train set: AUC={auc_train:0.3f}")
            ax.plot([0, 1], [0, 1], linestyle='dashed')
            ax.legend()
            ax.set_xlim(0., 1.)
            ax.set_ylim(0., 1.)
            ax.set_xlabel("False Positive Rate")
            ax.set_ylabel("True Positive Rate")
            cms_label(ax, year=year)
            print(auc_train, auc_test)

    def plot_train_breakdown(self, year, signal='Signal', use_test=True, bins=None):
        name = 'test' if use_test else 'train'
        if use_test:
            pred = self.pred_test[year][signal]
            df = self.test_sets[year]
        else:
            pred = self.pred_train[signal]
            df = self.train_set


        nbins = 50
        if bins is None:
            bins = np.linspace(0, 1, nbins+1)
            plotname = f'{name}_breakdown_{year}.pdf'
        else:
            plotname = f'{name}_breakdown_{year}_{bins[0]}-{bins[-1]}.png'
        # bins = np.concatenate([[0], np.linspace(self.minmax[0], self.minmax[1], nbins+1), [1]])
        def get_hist(pred, df, group):
            mask = self.get_mask(df, self.group_names[group])
            return np.histogram(pred[mask], bins, weights=df.scale_factor[mask])[0]

        with plot(self.outdir/plotname) as ax:
            ttt_hist = get_hist(pred, df, 'ttt')
            tttt_hist = get_hist(pred, df, 'tttt')
            ttX_hist = get_hist(pred, df, 'ttX')
            dd_hist = get_hist(pred, df, 'dd')
            other_hist = get_hist(pred, df, 'other')
            stack = np.array([other_hist, dd_hist, ttX_hist])
            colors = ['blueviolet', 'cornflowerblue', 'olivedrab']

            ax.hist(x=bins[:-1], bins=bins, weights=ttt_hist/np.sum(ttt_hist), color='r', label="ttt", linewidth=3, histtype='step')
            ax.hist(x=bins[:-1], bins=bins, weights=tttt_hist/np.sum(tttt_hist), color='orange', label="4top", linewidth=3, histtype='step')
            if np.sum(stack) != 0:
                n, bins, patches = ax.hist(
                    weights=stack.T/np.sum(stack), bins=bins, x=np.tile(bins[:-1], (len(stack), 1)).T,
                    label=['Other', "Data-Driven", "ttX"], histtype='stepfilled', stacked=True,
                    color=colors
                )
                for p in patches:
                    p.set(ec='#000000')

            generic_plot_setup(ax, year, bins=bins)


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
                    test = self.test_sets[year]
                    mask = test.sampleName == value
                    if np.count_nonzero(mask) > 0:
                        f[f'{year}/{sample}'] = {
                            "BDT": pred[mask],
                            'scale_factor': self.test_sets[year].scale_factor[mask],
                            'NLeps': test['NMuons'][mask] + test['NElectrons'][mask],
                            'NJets': test['NJets'][mask],
                        }
                        f[f'{year}/weights/{sample}'] = self.test_weights[year][sample]
                # Training set
                mask = self.train_set.sampleName == value
                if np.count_nonzero(mask) > 0:
                    f[f'train/{sample}'] = {
                        "BDT": self.pred_train[mask],
                        'scale_factor': self.train_set.scale_factor[mask]
                    }
                # Validation set
                mask = self.validation_set.sampleName == value
                if np.count_nonzero(mask) > 0:
                    f[f'validation/{sample}'] = {
                        "BDT": self.pred_validation[mask],
                        'scale_factor': self.validation_set.scale_factor[mask]
                    }

    def output_all(self, outfile, sig_out='Signal'):
        with uproot.recreate(outfile) as f:
            for year, pred in self.pred_test.items():
                test = copy(self.test_sets[year])
                test.insert(0, 'BDT', pred)
                for sample, value in self.sample_map.items():
                    if sample != 'data':
                        continue
                    mask = test.sampleName == value
                    if np.count_nonzero(mask) > 0:
                        f[f'{year}/{sample}'] = test[mask]
                        f[f'{year}/weights/{sample}'] = self.test_weights[year][sample]
