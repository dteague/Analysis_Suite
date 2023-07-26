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

from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.model_selection import train_test_split

from analysis_suite.plotting.utils import likelihood_sig
from analysis_suite.commons.constants import lumi

pd.options.mode.chained_assignment = None

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
        self.classID_by_className = {"Signal": 1, "Background": 0, "NotTrained": 0, "OnlyTrain": 0}
        self.group_dict = groupDict

        self.sample_map = dict()
        self.systName = systName
        self.region = region

        nonTrain_vars = ["scale_factor"]
        derived_vars = ["classID", "sampleName", "train_weight", "split_weight"]
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
        self.auc = dict()
        self.fom = dict()

        self.best_iter = 0


    def __bool__(self):
        return bool(self.test_sets)
    
    def should_train(self, df, className):
        enough_events = len(df)*self.split_ratio > self.min_train_events
        trainable_class = className != "NotTrained"
        is_SR_nominal = self.region == "signal" and self.systName == "Nominal"
        return enough_events and trainable_class and is_SR_nominal

    def read_in_file(self, directory, year="2018", rerun=False):
        test_set= setup_pandas(self.all_vars)
        test_weights = dict()
        if (directory/year/f'test_{self.systName}_{self.region}.root').exists() and not rerun:
            return
        infile = directory/year/f'processed_{self.systName}_{self.region}.root'
        for className, df, sample, weights in self.read_in_dataframe(infile):
            test_set = pd.concat([df, test_set], ignore_index=True)
            test_weights[sample] = weights
        self.test_sets[year] = test_set
        self.test_weights[year] = test_weights

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

    def get_total_events(self, infile):
        with uproot.open(infile) as f:
            sample_wgt = f['ttw']['scale_factor'].array()
            return len(sample_wgt)/0.2*self.split_ratio


    def read_in_dataframe(self, infile):
        with uproot.open(infile) as f:
            allSet = set([name[:name.index(";")] for name in f.keys() if "/" not in name])
            self.update_sample_map(allSet)
            for className, samples in self.group_dict.items():
                for sample in samples:
                    if sample not in f:
                        print(f"{sample} not found")
                        continue
                    # print(sample, "loaded")
                    df = f[sample].arrays(self._file_vars, library="pd")
                    logging.debug(f'{sample}, {len(df)}, {sum(df.scale_factor)}, {sum(df.scale_factor)/len(df)}')
                    df.loc[:, "classID"] = self.classID_by_className[className]
                    df.loc[:, "sampleName"] = self.sample_map[sample]
                    df.loc[:, "train_weight"] = sum(df.scale_factor) / len(df)
                    df.loc[:, "split_weight"] = np.ones(len(df))
                    all_weights = f[f'weights/{sample}'].arrays(library='pd')

                    yield className, df, sample, all_weights


    def setup_year(self, directory, year="2018", save_train=False):
        """**Fill the dataframes with all info in the input files**

        This grabs all the variable information about each sample,
        does some preliminary weighting and splits the data into the
        test and train set (based on `self.split_ratio`)

        Args:
            directory(string): Path to directory where root files are kept
        """
        train_set= setup_pandas(self.all_vars)
        test_set= setup_pandas(self.all_vars)
        test_weights = dict()
        validation_set= setup_pandas(self.all_vars)
        infile = directory / year / f'processed_{self.systName}_{self.region}.root'
        train_total = {className: self.total_class(infile, className) for className in self.group_dict.keys()}
        train_total["OnlyTrain"] = train_total["OnlyTrain"] + train_total["Background"]
        train_total["Background"] = train_total["OnlyTrain"]
        total_train = self.get_total_events(infile)

        for className, df, sample, weights in self.read_in_dataframe(infile):
            split_ratio = self.split_ratio
            train_percent = np.sum(df.scale_factor)/train_total[className]*total_train/len(df)
            # print(f"{sample}: {np.sum(df.scale_factor)/train_total[className]:0.3f}")
            if className == "OnlyTrain":
                if len(df) < self.min_train_events:
                    continue
                elif train_percent > 1:
                    df.loc[:, "split_weight"] *= train_percent
                    train = df
                elif train_percent*len(df) <= self.min_train_events:
                    continue
                else:
                    _, train = self.split(df, train_percent)
                # print(sample, len(train), ": not used", len(df)-len(train))
            elif not self.should_train(df, className):
                test_set = pd.concat([df, test_set], ignore_index=True)
                test_weights[sample] =  weights
                print(sample, "skip")
                continue
            else:
                if train_percent < split_ratio and train_percent*len(df) > self.min_train_events:
                    split_ratio = train_percent
                elif train_percent*len(df) <= self.min_train_events:
                    split_ratio = self.min_train_events/len(df)
                    df.loc[:, "split_weight"] *= train_percent/split_ratio
                else:
                    df.loc[:, "split_weight"] *= train_percent/split_ratio
                test, train = self.split(df, split_ratio)
                test_weights[sample] =  weights.loc[test.index]*len(df)/len(test)
                test_set = pd.concat([test, test_set], ignore_index=True)
                # print(sample, len(train), len(test))
            if self.validation_ratio > 0. and len(train)*self.validation_ratio > 1:
                train, validation = self.split(train, self.validation_ratio)
                validation_set = pd.concat([validation, validation_set], ignore_index=True)
            train_set = pd.concat([train, train_set], ignore_index=True)

        if save_train:
            self._output(train_set, directory / year / f"train_{self.systName}.root")
            self._output(validation_set, directory / year / f"validation_{self.systName}.root")

        self.train_set = pd.concat([self.class_reweight(train_set),
                                    self.train_set], ignore_index=True)
        self.validation_set = pd.concat([validation_set, self.validation_set],
                                        ignore_index=True)
        self.test_sets[year] = test_set
        self.test_weights[year] = test_weights

    def split(self, workset, split_ratio):
        train, test = train_test_split(workset, train_size=split_ratio, random_state=self.random_state)
        test.loc[:, ["scale_factor", "train_weight"]] *= len(workset)/len(test)
        train.loc[:, ["scale_factor", "train_weight"]] *= len(workset)/len(train)
        return test, train


    def update_sample_map(self, allSet):
        for sample in (allSet - set(self.sample_map)):
            self.sample_map[sample] = len(self.sample_map)

    def class_reweight(self, workset):
        for className, classID in self.classID_by_className.items():
            # print(className, classID)
            class_mask = workset["classID"] == classID
            class_set = workset[class_mask]
            if not len(class_set):
                continue
            scale = len(class_set)/sum(class_set.train_weight)
            # print(scale)
            workset.loc[class_mask, "train_weight"] *= scale
        return workset


    def apply_model(self, directory, year, get_auc=False):
        use_set = self.test_sets[year]
        pred = self.predict(use_set, directory)
        weights = use_set.scale_factor
        labels = use_set.classID.astype(int)

        self.pred_test[year] = {grp: pred.T[i] for grp, i in self.classID_by_className.items()}

        if get_auc:
            self.auc[year] = roc_auc_score(labels, pred.T[1], sample_weight = abs(weights))

            self.fom[year] = 0
            fom_bins = np.linspace(0, 1, 101)
            sig = np.cumsum(np.histogram(self.pred_test[year]["Signal"][labels==1], bins=fom_bins,
                                     weights=weights[labels==1])[0][::-1])[::-1]
            tot = np.cumsum(np.histogram(self.pred_test[year]["Signal"], bins=fom_bins,
                                         weights=weights)[0][::-1])[::-1]
            self.fom[year] = max(sig/np.sqrt(tot))
            fom_maxbin = np.argmax(sig/np.sqrt(tot))

            print(f'AUC for year {year}: {self.auc[year]}')
            print(f'FOM (cut) for year {year}: {self.fom[year]:0.3f} at val {fom_bins[fom_maxbin]:0.2f}')

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
        bkg = np.sum(cut_set[cut_set.classID == 0].scale_factor)
        fom = sig/np.sqrt(sig+bkg)

        matthew_coef = (tp/len(truth_vals)-s*p)/np.sqrt(p*s*(1-p)*(1-s))
        print(f'Cut {cut:0.3f} for year {year}: {precision:0.3f} {recall:0.3f} {f1_score:0.3f} {matthew_coef:0.3f} {fom:0.3f}')

    def roc_curve(self, directory, year):

        def get_roc(dataset, pred, scale):
            nbins = 100
            bins = np.linspace(0, 1+1/nbins, nbins+2)
            s_hist = np.histogram(pred[dataset.classID==1], bins, weights=scale[dataset.classID == 1])[0]
            b_hist = np.histogram(pred[dataset.classID==0], bins, weights=scale[dataset.classID == 0])[0]
            tp = np.cumsum(s_hist[::-1])/np.sum(s_hist)
            fp = np.cumsum(b_hist[::-1])/np.sum(b_hist)

            delta = fp[1:]-fp[:-1]
            trap = (tp[1:]+tp[:-1])/2
            auc = np.sum(delta*trap)
            return tp, fp, auc

        test = self.test_sets[year]
        nt_nums = [self.sample_map[key] for key in self.group_dict['NotTrained'] if key in self.sample_map]
        nontrain_mask = np.any([test.sampleName == num for num in nt_nums], axis=0)
        pred_test = self.predict(test, directory).T[1]
        scale_test = test.scale_factor
        tp_test, fp_test, auc_test = get_roc(test, pred_test, scale_test)

        ttt_names = ['tttw', 'tttj']
        ttX_names = ['ttz', 'ttw', 'tth', 'ttz_m1-10']
        data_drive_names = ['nonprompt', 'charge_flip']
        other_names = [key for key in self.group_dict['Background'] if key not in ['tttt']+ttX_names]

        tttt_nums = [self.sample_map[key] for key in ['tttt']+ttt_names if key in self.sample_map]
        ttX_nums = [self.sample_map[key] for key in ttX_names+ttt_names if key in self.sample_map]
        data_drive_nums = [self.sample_map[key] for key in data_drive_names+ttt_names if key in self.sample_map]
        other_nums = [self.sample_map[key] for key in other_names+ttt_names if key in self.sample_map]

        tttt_mask = np.any([test.sampleName == num for num in tttt_nums], axis=0)
        ttX_mask = np.any([test.sampleName == num for num in ttX_nums], axis=0)
        data_drive_mask = np.any([test.sampleName == num for num in data_drive_nums], axis=0)
        other_mask = np.any([test.sampleName == num for num in other_nums], axis=0)

        tp_4top, fp_4top, auc_4top = get_roc(test[tttt_mask], pred_test[tttt_mask], scale_test[tttt_mask])
        tp_ttX, fp_ttX, auc_ttX = get_roc(test[ttX_mask], pred_test[ttX_mask], scale_test[ttX_mask])
        tp_dd, fp_dd, auc_dd = get_roc(test[data_drive_mask], pred_test[data_drive_mask], scale_test[data_drive_mask])
        tp_other, fp_other, auc_other = get_roc(test[other_mask], pred_test[other_mask], scale_test[other_mask])

        print(auc_test, auc_4top, auc_ttX, auc_dd, auc_other)
        # exit()
        # ttz_sig_mask = test.sampleName


        ot_nums = [self.sample_map[key] for key in self.group_dict['OnlyTrain'] if key in self.sample_map]
        train = self.train_set
        if len(ot_nums):
            train = train[np.all([train.sampleName != num for num in ot_nums], axis=0)]
        train = pd.concat((train, test[nontrain_mask]))
        pred_train = self.predict(train, directory).T[1]
        scale_train = train.scale_factor
        tp_train, fp_train, auc_train = get_roc(train, pred_train, scale_train)

        from analysis_suite.commons.plot_utils import plot, hep
        with plot(directory/f"roc_{year}.png") as ax:
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

        with plot(directory/f"roc_sep_{year}.png") as ax:
            ax.plot(fp_test, tp_test, linewidth=3, label=f"Total set: AUC={auc_test:0.3f}")
            ax.plot([0, 1], [0, 1], linestyle='dashed')

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




    def overtrain_test(self, directory, year):
        train = self.train_set
        test = self.test_sets[year]
        ot_nums = [self.sample_map[key] for key in self.group_dict['OnlyTrain'] if key in self.sample_map]
        if len(ot_nums):
            train = train[np.all([train.sampleName != num for num in ot_nums], axis=0)]
        nt_nums = [self.sample_map[key] for key in self.group_dict['NotTrained'] if key in self.sample_map]
        nontrain_mask = np.any([test.sampleName == num for num in nt_nums], axis=0)

        train_pred = self.predict(train, directory).T[1]
        test_pred = self.predict(test, directory).T[1]

        nbins = 15
        bins = np.linspace(0, 1, nbins+1)

        def get_hist(pred, scale, sampleNames, names):
            nums = [self.sample_map[key] for key in names if key in self.sample_map]
            mask = np.any([sampleNames == num for num in nums], axis=0)
            return np.histogram(pred[mask], bins, weights=scale[mask])[0]

        train_b = np.histogram(train_pred[train.classID==0], bins, weights=train.scale_factor[train.classID==0])[0]
        nontrain = np.histogram(test_pred[nontrain_mask], bins, weights=test.scale_factor[nontrain_mask])[0]
        train_b = (nontrain+train_b)
        train_s = np.histogram(train_pred[train.classID==1], bins, weights=train.scale_factor[train.classID==1])[0]

        test_b = np.histogram(test_pred[test.classID==0], bins, weights=test.scale_factor[test.classID==0])[0]
        test_s = np.histogram(test_pred[test.classID==1], bins, weights=test.scale_factor[test.classID==1])[0]

        kw_hist = {"alpha": 0.3, "hatch": '///', "histtype": "stepfilled"}
        kw_err = {"markersize": 4, "fmt": "o"}

        fom = lambda s, b: np.sqrt(2*np.sum((s+b)*np.log(1+s/(b+1e-5))-s))
        print(f"Train FOM: {fom(train_s, train_b)}")
        print(f"Test FOM: {fom(test_s, test_b)}")


        from analysis_suite.commons.plot_utils import plot, hep
        with plot(directory/f"overtrain_{year}.png") as ax:
            ax.hist(x=bins[:-1], bins=bins, weights=test_s/np.sum(test_s), color='r', label="Signal (test)", **kw_hist)
            ax.hist(x=bins[:-1], bins=bins, weights=test_b/np.sum(test_b), color='b', label="Background (test)", **kw_hist)
            ax.errorbar(x=bins[:-1]+1/(2*nbins), xerr=1/(2*nbins), y=train_s/np.sum(train_s), color='r', label="Signal (train)", **kw_err)
            ax.errorbar(x=bins[:-1]+1/(2*nbins), xerr=1/(2*nbins), y=train_b/np.sum(train_b), color='b', label="Background (train)", **kw_err)
            ax.set_xlim(0., 1.)
            ax.set_xlabel("$disc_{BDT}$")
            ax.set_ylabel("A.U.")
            ax.legend()
            hep.cms.label(ax=ax, lumi=lumi[year], label="Preliminary")

        with plot(directory/f'train_breakdown_{year}.png') as ax:
            tttt_hist = get_hist(test_pred, test.scale_factor, test.sampleName, ['tttt'])
            ttX_hist = get_hist(test_pred, test.scale_factor, test.sampleName, ['ttz', 'ttw', 'tth', 'ttz_m1-10'])
            dd_hist = get_hist(test_pred, test.scale_factor, test.sampleName, ['nonprompt', 'charge_flip'])
            other_hist = test_b - (tttt_hist + ttX_hist + dd_hist)
            ax.hist(x=bins[:-1], bins=bins, weights=test_s/np.sum(test_s), color='r', label="Signal", linewidth=3, histtype='step')
            stack = np.array([other_hist, dd_hist, ttX_hist])/np.sum(test_b)
            print(np.sum(stack))
            n, bins, patches = ax.hist(
                weights=stack.T, bins=bins, x=np.tile(bins[:-1], (len(stack), 1)).T,
                label=['Other', "Data-Driven", "ttX"], histtype='stepfilled', stacked=True,
                color=['blueviolet', 'cornflowerblue', 'olivedrab'],
            )
            ax.hist(x=bins[:-1], bins=bins, weights=tttt_hist/np.sum(tttt_hist), color='orange', label="tttt", linewidth=3, histtype='step')
            ax.set_xlim(0., 1.)
            ax.set_xlabel("$disc_{BDT}$")
            ax.set_ylabel("A.U.")
            ax.legend()
            hep.cms.label(ax=ax, lumi=lumi[year], label="Preliminary")

    # Private Functions
    def _cut_mask(self, frame):
        """**Reduce frame using root style cut string**

        Args:
          frame(pandas.DataFrame): DataFrame to cut on

        """
        mask = np.ones(len(frame))
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
        self.cuts = cut_string


    def output(self, outdir, year):
        workSet = self.test_sets[year]
        for key, arr in self.pred_test[year].items():
            workSet.insert(0, key, arr)
        self._output(workSet, outdir / year / f"test_{self.systName}_{self.region}.root", self.test_weights[year])


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
