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
from analysis_suite.commons.plot_utils import plot, hep, plot_colorbar

pd.options.mode.chained_assignment = None

def generic_plot_setup(ax, year, bins=None):
    if bins is None:
        ax.set_xlim(0., 1.)
    else:
        ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel("$disc_{BDT}$")
    ax.set_ylabel("A.U.")
    ax.legend()
    hep.cms.label(ax=ax, lumi=lumi[year], label="Preliminary")

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
    def __init__(self, use_vars, groupDict, region="signal", systName="Nominal", **kwargs):
        """Constructor method
        """
        self.classID_by_className = {"Signal": 1, "Background": 0, "NotTrained": 0, "OnlyTrain": 0, "4top": 2}
        self.group_dict = groupDict
        samples = list()
        for val_list in self.group_dict.values():
            samples += list(val_list)
        self.sample_map = {val: i for i, val in enumerate(sorted(samples))}
        self.systName = systName
        self.region = region

        nonTrain_vars = ["scale_factor"]
        derived_vars = ["classID", "sampleName", "train_weight", "split_weight"]
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

        self.group_names = {
            'ttt': ['tttw', 'tttj'],
            'ttX': ['ttz', 'ttw', 'tth', 'ttz_m1-10'],
            'dd': ['nonprompt', 'charge_flip']+list(self.group_dict["OnlyTrain"]),
            'tttt': ['tttt']
        }
        used_groups = [i for sub in self.group_names.values() for i in sub]
        if 'Background' in self.group_dict:
            self.group_names['other'] =  [key for key in self.group_dict['Background'] if key not in used_groups]


    def __bool__(self):
        return bool(self.test_sets)

    def set_outdir(self, outdir):
        self.outdir = outdir

    def should_train(self, df, className):
        enough_events = len(df)*self.split_ratio > self.min_train_events
        trainable_class = className != "NotTrained"
        is_SR_nominal = self.region == "signal" and self.systName == "Nominal"
        return enough_events and trainable_class and is_SR_nominal

    def read_in_file(self, directory, year="2018", typ='test', mask=None):
        test_set= setup_pandas(self.all_vars)
        test_weights = dict()

        infile = directory/year/f'{typ}_{self.systName}_{self.region}.root'
        for className, df, sample, weights in self.read_in_dataframe(infile):
            test_set = pd.concat([df, test_set], ignore_index=True)
            test_weights[sample] = weights
        self.test_sets[year] = test_set
        self.test_weights[year] = test_weights

        # Try to get train sets
        infile = directory/'train_files'/f'train_{self.systName}_{self.region}.root'
        if infile.exists():
            for className, df, sample, weights in self.read_in_dataframe(infile):
                self.train_set = pd.concat([df, self.train_set], ignore_index=True)
        infile = directory/'train_files'/f'validation_{self.systName}_{self.region}.root'
        if infile.exists():
            for className, df, sample, weights in self.read_in_dataframe(infile):
                self.validation_set = pd.concat([df, self.validation_set], ignore_index=True)


    def total_class(self, infile, className):
        total_wgt = 0.
        with uproot.open(infile) as f:
            for sample in self.group_dict[className]:
                if sample not in f:
                    continue
                sample_wgt = f[sample]["scale_factor"].array()
                if len(sample_wgt) < self.min_train_events/self.split_ratio:
                    continue
                # print(sample, np.sum(sample_wgt), len(sample_wgt))
                total_wgt += np.sum(sample_wgt)
        return total_wgt


    def read_in_dataframe(self, infile):
        with uproot.open(infile) as f:
            for className, samples in self.group_dict.items():
                for sample in samples:
                    if sample not in f:
                        logging.debug(f"{sample} not found")
                        continue
                    # print(sample, "loaded")
                    df = f[sample].arrays(library="pd")
                    mask = self._cut_mask(df)
                    df = df[self._file_vars][mask]
                    if df.empty:
                        continue
                    logging.debug(f'{sample}, {len(df)}, {sum(df.scale_factor)}, {sum(df.scale_factor)/len(df)}')
                    df.loc[:, "classID"] = self.classID_by_className[className]
                    df.loc[:, "sampleName"] = self.sample_map[sample]
                    df.loc[:, "train_weight"] = sum(df.scale_factor) / len(df)
                    df.loc[:, "split_weight"] = np.ones(len(df))
                    if "weights" not in f:
                        all_weights = None
                    else:
                        all_weights = f[f'weights/{sample}'].arrays(library='pd')[mask]
                    yield className, df, sample, all_weights

    def setup_weights(self, filename):
        self.train_total = {className: self.total_class(filename, className) for className in self.group_dict.keys()}
        self.train_total["OnlyTrain"] += self.train_total["Background"]
        self.train_total["Background"] = self.train_total["OnlyTrain"]
        with uproot.open(filename) as f:
            sample_wgt = f['ttw']['scale_factor'].array()
            self.total_train = len(sample_wgt)/0.2*self.split_ratio


    def setup_year(self, directory, year="2018", split=True):
        """**Fill the dataframes with all info in the input files**

        This grabs all the variable information about each sample,
        does some preliminary weighting and splits the data into the
        test and train set (based on `self.split_ratio`)

        Args:
            directory(string): Path to directory where root files are kept
        """
        self.test_sets[year] = setup_pandas(self.all_vars)
        self.test_weights[year] = dict()
        infile = directory / year / f'processed_{self.systName}_{self.region}.root'
        self.setup_weights(infile)
        for className, df, sample, weights in self.read_in_dataframe(infile):
            test, train, validation = self.setup_split(df, className, sample, split)
            # print(sample, len(df), len(test), len(train), len(validation))
            if not test.empty:
                self.test_weights[year][sample] = weights.loc[test.index]*len(df)/len(test)
                self.test_sets[year] = pd.concat([test, self.test_sets[year]], ignore_index=True)
                # print(sample, np.sum(test.scale_factor), np.sum(test.split_weight))
            if not train.empty:
                self.train_set = pd.concat([train, self.train_set], ignore_index=True)
                # print(sample, "train", np.sum(train.scale_factor), np.sum(train.split_weight))
            if not validation.empty:
                self.validation_set = pd.concat([validation, self.validation_set], ignore_index=True)
        exit()

    def setup_split(self, df, className, sample, split=True):
        train = setup_pandas(self.all_vars)
        test = setup_pandas(self.all_vars)
        validation = setup_pandas(self.all_vars)
        total_scale = np.sum(df.scale_factor)

        if not split:
            return df, train, validation

        exp_train_evt = total_scale/self.train_total[className]*self.total_train
        split_ratio = self.split_ratio

        if sample == "tttj" or sample == "tttw":
            exp_train_evt = split_ratio*len(df)

        # Only train, i.e. use all the events possible
        if className == "OnlyTrain":
            if exp_train_evt <= self.min_train_events:
                pass
            elif exp_train_evt > len(df):
                df.loc[:, "split_weight"] *= exp_train_evt/len(df)
                train = df
            else:
                _, train = self.split(df, exp_train_evt/len(df), total_scale)
        # Dont Train (not enough events)
        elif not self.should_train(df, className):
            test = df
        # Do Train
        else:
            if exp_train_evt < len(df)*split_ratio and exp_train_evt > self.min_train_events:
                split_ratio = exp_train_evt/len(df)
            elif exp_train_evt <= self.min_train_events:
                split_ratio = self.min_train_events/len(df)
                df.loc[:, "split_weight"] *= exp_train_evt/self.min_train_events
            else:
                df.loc[:, "split_weight"] *= exp_train_evt/(split_ratio*len(df))
            test, train = self.split(df, split_ratio, total_scale)
        if self.validation_ratio > 0. and len(train)*self.validation_ratio > 1:
            train, validation = self.split(train, self.validation_ratio, total_scale)

        print(sample,  np.sum(df.scale_factor), np.sum(test.scale_factor), np.sum(train.scale_factor), np.sum(validation.scale_factor))
        if len(train) > 0 and len(test) > 0:
            print("NJets", np.average(df.NJets, weights=df.scale_factor), np.average(test.NJets, weights=test.scale_factor), np.average(train.NJets, weights=train.scale_factor))
            print("NlooseBJets", np.average(df.NlooseBJets, weights=df.scale_factor), np.average(test.NlooseBJets, weights=test.scale_factor), np.average(train.NlooseBJets, weights=train.scale_factor))
            print("HT", np.average(df.HT, weights=df.scale_factor), np.average(test.HT, weights=test.scale_factor), np.average(train.HT, weights=train.scale_factor))
            print("j1Pt", np.average(df.j1Pt, weights=df.scale_factor), np.average(test.j1Pt, weights=test.scale_factor), np.average(train.j1Pt, weights=train.scale_factor))
        print()
        return test, train, validation

        
    def split(self, workset, split_ratio, target_scale):
        train, test = train_test_split(workset, train_size=split_ratio, random_state=self.random_state)
        # test.loc[:, ["scale_factor", "train_weight"]] *= target_scale/np.sum(test.scale_factor)
        # train.loc[:, ["scale_factor", "train_weight"]] *= target_scale/np.sum(train.scale_factor)
        test.loc[:, ["scale_factor", "train_weight"]] *= len(workset)/len(test)
        train.loc[:, ["scale_factor", "train_weight"]] *= len(workset)/len(train)

        return test, train


    def apply_model(self, year, skip_train=False):
        def get_pred(use_set, directory):
            unique_labels = np.unique(use_set.classID.astype(int))
            pred = self.predict(use_set, directory)
            return {grp: pred.T[i] for grp, i in self.classID_by_className.items() if i in unique_labels}
        self.pred_test[year] = get_pred(self.test_sets[year], self.outdir)
        if not skip_train:
            self.pred_train = get_pred(self.train_set, self.outdir)
            self.pred_validation = get_pred(self.validation_set, self.outdir)


    def get_stats(self, year, cut):
        truth_vals = self.test_sets[year].classID.astype(int)
        pred_mask = self.pred_test[year]["Signal"] > cut
        tn, fp, fn, tp = confusion_matrix(truth_vals, pred_mask).ravel()
        precision = tp / (tp+fp)
        recall = tp / (tp + fn)
        f1_score = 2*(precision*recall)/(precision+recall)
        s = (tp+fn)/len(truth_vals)
        p = (tp+fp)/len(truth_vals)

        cut_set = self.test_sets[year][pred_mask]

        sig = np.sum(cut_set[cut_set.classID == 1].scale_factor)
        bkg = np.sum(cut_set[cut_set.classID != 1].scale_factor)
        fom = sig/np.sqrt(sig+bkg+1e-5)

        matthew_coef = (tp/len(truth_vals)-s*p)/np.sqrt(p*s*(1-p)*(1-s))
        print(f'Cut {cut:0.3f} for year {year}: {precision:0.3f} {recall:0.3f} {f1_score:0.3f} {matthew_coef:0.3f} {fom:0.3f}')


    def get_mask(self, df, groups):
        num = [self.sample_map[key] for key in groups if key in self.sample_map]
        if len(num) == 0:
            return np.full(len(df), False)
        return np.any([df.sampleName == i for i in num], axis=0)

    def get_corr(self, typ='train', year=None):
        import matplotlib.pyplot as plt

        if typ == 'train':
            work_set = self.train_set[self.use_vars]
            weights = abs(self.train_set.scale_factor.to_numpy())
            filename = 'corr_train.png'
        elif typ == 'test' and year is not None:
            work_set = self.test_sets[year][self.use_vars]
            weights = abs(self.test_sets[year].scale_factor.to_numpy())
            filename = f'corr_test_{year}.png'

        cov = np.cov(work_set.to_numpy(), rowvar=False, aweights=weights)
        stddev = np.sqrt(np.diag(cov))
        corr = cov/(stddev[:, None]*stddev[None, :])

        with plot(self.outdir/filename) as ax:
            nVars = corr.shape[0]
            x = np.linspace(0, nVars, nVars+1)
            xx = np.tile(x, (nVars+1, 1))
            yy = np.tile(x, (nVars+1, 1)).T

            color_plot = ax.pcolormesh(xx, yy, corr, shading='flat')
            plot_colorbar(color_plot, ax)
            ax.minorticks_off()
            ax.tick_params(direction='out', length=0)
            ax.set_xticks(np.arange(nVars)+0.5)
            ax.set_xticklabels(work_set.columns, rotation=70, ha='right', rotation_mode="anchor", fontsize=10)
            ax.set_yticks(np.arange(nVars)+0.5)
            ax.set_yticklabels(work_set.columns, fontsize=10)

        #     for i in range(work_set.shape[1]):
        #         for j in range(work_set.shape[1]):
        #             text = ax.text(j, i, corr[i, j], ha="center", va="center", color="w")
            ax.set_title("Correlation")

    def get_auc(self, year, signal="Signal"):
        classID = self.classID_by_className[signal]
        df = self.test_sets[year]
        pred = self.pred_test[year][signal]
        is_sig = df.classID == classID

        nbins = 100
        bins = np.linspace(0, 1+1/nbins, nbins+2)

        s_hist = np.histogram(pred[is_sig], bins, weights=df.scale_factor[is_sig])[0]
        b_hist = np.histogram(pred[~is_sig], bins, weights=df.scale_factor[~is_sig])[0]
        tp = np.cumsum(s_hist[::-1])/(np.sum(s_hist)+1e-5)
        fp = np.cumsum(b_hist[::-1])/(np.sum(b_hist)+1e-5)

        delta = fp[1:]-fp[:-1]
        trap = (tp[1:]+tp[:-1])/2
        auc = np.sum(delta*trap)
        return auc


    def roc_curve(self, year, signal="Signal"):
        classID = self.classID_by_className[signal]
        test = self.test_sets[year]
        pred_test = self.pred_test[year][signal]
        is_sig_test = test.classID == classID
        train = self.train_set
        pred_train = self.pred_train[signal]
        is_sig_train = train.classID == classID

        def get_roc(pred, df, is_sig, group=None):
            nbins = 100
            bins = np.linspace(0, 1+1/nbins, nbins+2)

            s_hist = np.histogram(pred[is_sig], bins, weights=df.scale_factor[is_sig])[0]
            if group is None:
                b_hist = np.histogram(pred[~is_sig], bins, weights=df.scale_factor[~is_sig])[0]
            else:
                mask = self.get_mask(df, self.group_names[group])
                b_hist = np.histogram(pred[mask], bins, weights=df.scale_factor[mask])[0]
            tp = np.cumsum(s_hist[::-1])/(np.sum(s_hist)+1e-5)
            fp = np.cumsum(b_hist[::-1])/(np.sum(b_hist)+1e-5)

            delta = fp[1:]-fp[:-1]
            trap = (tp[1:]+tp[:-1])/2
            auc = np.sum(delta*trap)
            return tp, fp, auc

        tp_test, fp_test, auc_test = get_roc(pred_test, test, is_sig_test)
        tp_train, fp_train, auc_train = get_roc(pred_train, train, is_sig_train)

        tp_ttt, fp_ttt, auc_ttt = get_roc(pred_test, test, is_sig_test, 'ttt')
        tp_4top, fp_4top, auc_4top = get_roc(pred_test, test, is_sig_test, 'tttt')
        tp_ttX, fp_ttX, auc_ttX = get_roc(pred_test, test, is_sig_test, 'ttX')
        tp_dd, fp_dd, auc_dd = get_roc(pred_test, test, is_sig_test, 'dd')
        tp_other, fp_other, auc_other = get_roc(pred_test, test, is_sig_test, 'other')

        with plot(self.outdir/f"roc_{year}.png") as ax:
            ax.plot(fp_test, tp_test, linewidth=3, label=f"Test set: AUC={auc_test:0.3f}")
            ax.plot(fp_train, tp_train, linewidth=3, label=f"Train set: AUC={auc_train:0.3f}")
            ax.plot([0, 1], [0, 1], linestyle='dashed')
            ax.legend()
            ax.set_xlim(0., 1.)
            ax.set_ylim(0., 1.)
            ax.set_xlabel("False Positive Rate")
            ax.set_ylabel("True Positive Rate")
            hep.cms.label(ax=ax, lumi=lumi[year], label="Preliminary")
            print(auc_train, auc_test)

        with plot(self.outdir/f"roc_sep_{year}.png") as ax:
            ax.plot(fp_test, tp_test, linewidth=3, label=f"Total set: AUC={auc_test:0.3f}")
            ax.plot([0, 1], [0, 1], linestyle='dashed', color='k')

            if abs(auc_ttt-0.5) > 1e-3:
                ax.plot(fp_ttt, tp_ttt, linewidth=3, label=f"ttt: AUC={auc_ttt:0.3f}", linestyle='dashed')
            if abs(auc_4top-0.5) > 1e-3:
                ax.plot(fp_4top, tp_4top, linewidth=3, label=f"4-Top: AUC={auc_4top:0.3f}", linestyle='dashed')
            ax.plot(fp_ttX, tp_ttX, linewidth=3, label=f"ttX: AUC={auc_ttX:0.3f}", linestyle='dashed')
            ax.plot(fp_dd, tp_dd, linewidth=3, label=f"Data-Driven Bkg: AUC={auc_dd:0.3f}", linestyle='dashed')
            ax.plot(fp_other, tp_other, linewidth=3, label=f"Other Bkgs: AUC={auc_other:0.3f}", linestyle='dashed')

            ax.legend()
            ax.set_xlim(0., 1.)
            ax.set_ylim(0., 1.)
            ax.set_xlabel("False Positive Rate")
            ax.set_ylabel("True Positive Rate")
            hep.cms.label(ax=ax, lumi=lumi[year], label="Preliminary")

    def get_fom(self, year, signal='Signal'):
        test = self.test_sets[year]
        pred_test = self.pred_test[year][signal]
        classID = self.classID_by_className[signal]

        nbins = 20
        bins = np.linspace(0, 1, nbins+1)

        def get_hist(pred, df, sig=False):
            mask = df.classID==classID if sig else df.classID!=classID
            return np.histogram(pred[mask], bins, weights=df.scale_factor[mask])[0]

        s = get_hist(pred_test, test, sig=True)
        b = get_hist(pred_test, test)


        with np.errstate(invalid='ignore'):
            fom = np.sqrt(2*np.sum((s+b)*np.log(np.where(b > 1e-5, 1+s/b, 1))-s))
        return fom


    def plot_overtrain(self, year, signal='Signal'):
        test = self.test_sets[year]
        pred_test = self.pred_test[year][signal]
        pred_train = self.pred_train[signal]
        train = self.train_set
        classID = self.classID_by_className[signal]

        nbins = 15
        bins = np.linspace(0, 1, nbins+1)

        def get_hist(pred, df, sig=False):
            mask = df.classID==classID if sig else df.classID!=classID
            return np.histogram(pred[mask], bins, weights=df.scale_factor[mask])[0]

        train_s = get_hist(pred_train, train, sig=True)
        train_b = get_hist(pred_train, train)
        test_s = get_hist(pred_test, test, sig=True)
        test_b = get_hist(pred_test, test)

        def fom(s, b):
            with np.errstate(invalid='ignore'):
                fom_val = np.sqrt(2*np.sum((s+b)*np.log(np.where(b > 1e-5, 1+s/b, 1))-s))
            return fom_val

        print(f"Train FOM: {fom(train_s, train_b)}")
        print(f"Test FOM: {fom(test_s, test_b)}")

        kw_hist = {"alpha": 0.3, "hatch": '///', "histtype": "stepfilled"}
        kw_err = {"markersize": 4, "fmt": "o"}
        with plot(self.outdir/f"overtrain_{year}.png") as ax:
            ax.hist(x=bins[:-1], bins=bins, weights=test_s/np.sum(test_s), color='r', label="Signal (test)", **kw_hist)
            ax.hist(x=bins[:-1], bins=bins, weights=test_b/np.sum(test_b), color='b', label="Background (test)", **kw_hist)
            ax.errorbar(x=bins[:-1]+1/(2*nbins), xerr=1/(2*nbins), y=train_s/np.sum(train_s), color='r', label="Signal (train)", **kw_err)
            ax.errorbar(x=bins[:-1]+1/(2*nbins), xerr=1/(2*nbins), y=train_b/np.sum(train_b), color='b', label="Background (train)", **kw_err)
            generic_plot_setup(ax, year)


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
            plotname = f'{name}_breakdown_{year}.png'
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
                # Apply patch to edge colors
                dark = 0.3
                for p, c in zip(patches, colors):
                    ec =  [i - dark if i > dark else 0.0 for i in clr.to_rgb(c)]
                    if isinstance(p, list):
                        p[0].set_ec(ec)
                    else:
                        p.set(ec=ec)

            generic_plot_setup(ax, year, bins=bins)


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


    def output(self, year, sig_out='Signal', signal="Signal"):
        workSet = self.test_sets[year]
        if sig_out not in workSet:
            workSet.insert(0, sig_out, self.pred_test[year][signal])
        self._output(workSet, self.outdir / year / f"test_{self.systName}_{self.region}.root", self.test_weights[year])


    def output_train(self, sig_out='Signal', signal='Signal'):
        train_out = self.outdir / 'train_files'
        train_out.mkdir(exist_ok=True, parents=True)
        if sig_out not in self.train_set:
            self.train_set.insert(0, sig_out, self.pred_train[signal])
        self._output(self.train_set, train_out / f"train_{self.systName}_{self.region}.root")
        if sig_out not in self.validation_set:
            self.validation_set.insert(0, sig_out, self.pred_validation[signal])
        self._output(self.validation_set, train_out / f"validation_{self.systName}_{self.region}.root")


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
