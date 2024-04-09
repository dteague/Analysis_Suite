#!/usr/bin/env python3
import logging
import numpy as np
import os
from importlib import import_module
from copy import copy
import pandas as pd
import uproot

from analysis_suite.commons.configs import get_list_systs, get_inputs, get_ntuple_info
import analysis_suite.commons.user as user

def setup(cli_args):
    inputs = get_inputs(cli_args.workdir)
    os.environ["NUMEXPR_MAX_THREADS"] = "8"

    argList = list()
    if cli_args.train:
        allSysts = ["Nominal"]
    else:
        allSysts = get_list_systs(cli_args.workdir, cli_args.tool, cli_args.systs)
    for syst in allSysts:
        argList.append((inputs.groups, inputs.usevars, cli_args.workdir, cli_args.model, cli_args.train,
                        cli_args.years, cli_args.region, syst, cli_args.save, cli_args.rerun))
    return argList

def get_set(typ, indir, model, year='2018', cut=0.99):
    extra = "_signal" if typ == "test" else ""
    final_df = pd.DataFrame()
    filename = f"{typ}_Nominal_signal.root" if typ == 'test' else f"{typ}_Nominal.root"
    file_dir = indir/year if typ == 'test' else indir/'train_files'

    with uproot.open(file_dir/filename) as f:
        for key, cls in f.classnames().items():
            if "/" in key or cls != "TTree": continue
            df = f[key].arrays(library='pd')
            classID = np.unique(df['classID'])[0]
            # print(key, np.sum(df.split_weight), len(df), np.unique(df.split_weight), np.sum(df.scale_factor), np.sum(df.split_weight)/np.sum(df.scale_factor))
            if "4top" in key:
                if typ != "test": continue
                df['classID'] = 0 if classID == 1 else 1
            elif "tttj" in key or 'tttw' in key:
                df['classID'] = 0 if classID == 1 else 1
            if "4top_sig" not in df.columns:
                print("HERE", key)
                pred = model.predict(df, indir).T[1]
                df['4top_sig'] = pred
            # print(key, np.unique(df["sampleName"]), len(df))
            final_df = pd.concat([final_df, df])
    final_df = final_df[final_df["4top_sig"] < cut]
    return final_df


def get_weights(indir, year):
    weights = dict()
    with uproot.open(indir/year/"test_Nominal_signal.root") as f:
        wgt_dir = f['weights']
        weights = {samp[:-2]: wgt_dir[samp].arrays(library='pd') for samp in wgt_dir}
    return weights


def run(groups, usevars, workdir, model_type, train, years, region, systName, save_train, rerun):
    module = import_module(f'.{model_type}', "analysis_suite.machine_learning")
    maker = getattr(module, next(filter(lambda x : "Maker" in x, dir(module))))

    remove_group = ['nonprompt_mc'] if region != "signal" else []
    ginfo = get_ntuple_info(region, remove=remove_group)
    groupDict = dict()
    for cls, samples in groups.items():
        if cls == "Signal":
            samples = ['4top']
        elif cls == "Background":
            samples = samples+['ttt']
        groupDict[cls] = ginfo.setup_members(filter(lambda s: s in ginfo.group2color, samples))

    # First Training
    model = maker(usevars, groupDict, region=region, systName=systName)
    first_params = {
        'colsample_bytree': 0.5283131479037887, 'eta': 0.037489501191992,
        'eval_metric': 'logloss', 'gamma': 5.439613578244588, 'max_depth': 3,
        'min_child_weight': 1.0, 'n_estimators': 500, 'subsample': 0.9
    }
    model.update_params(first_params)
    output_first = workdir/'first_train'
    model.set_outdir(output_first)
    for year in years:
        (output_first/year).mkdir(exist_ok=True, parents=True)
        model.setup_year(workdir, year)
    model.train()
    model.plot_training_progress()
    for year in years:
        model.apply_model(year, get_auc=True)
        model.roc_curve(year)
        model.plot_overtrain(year)
        model.plot_train_breakdown(year)
        model.plot_train_breakdown(year, use_test=False)
        model.output(output_first, year, '4top_sig')
    model.output_train(output_first, '4top_sig')

    print(f"Training for syst {systName}")

    groupDict = dict()
    for cls, samples in groups.items():
        if cls == "Signal":
            samples = ['ttt']
        elif cls == "NotTrained":
            samples = samples + ['4top']
        groupDict[cls] = ginfo.setup_members(filter(lambda s: s in ginfo.group2color, samples))
    model = maker(usevars, groupDict, region=region, systName=systName)
    second_params = {
        'colsample_bytree': 0.5283131479037887, 'eta': 0.037489501191992,
        'eval_metric': 'logloss', 'gamma': 5.439613578244588, 'max_depth': 3,
        'min_child_weight': 1.0, 'n_estimators': 500, 'subsample': 0.9
    }
    model.update_params(second_params)
    output_second = workdir/'second_train'
    model.set_outdir(output_second)
    cut = 0.93
    for year in years:
        (output_second/year).mkdir(exist_ok=True, parents=True)
        model.test_sets[year] = get_set("test", output_first, model, year=year, cut=cut)
        model.test_weights[year] = get_weights(output_first, year)
    model.train_set = get_set("train", output_first, model, cut=cut)
    model.validation_set = get_set("validation", output_first, model, cut=cut)
    model.train()
    model.plot_training_progress()
    for year in years:
        model.apply_model(year, get_auc=True)
        model.roc_curve(year)
        model.plot_overtrain(year)
        model.plot_train_breakdown(year)
        model.plot_train_breakdown(year, use_test=False)
        model.output(output_first, year, '3top_sig')
    model.output_train(output_first, '3top_sig')


# def run(groupDict, usevars, workdir, model, train, years, region, systName, save_train, rerun):
#     module = import_module(f'.{model}', "analysis_suite.machine_learning")
#     maker = getattr(module, next(filter(lambda x : "Maker" in x, dir(module))))
#     ml_runner = maker(usevars, groupDict, region=region, systName=systName)
#     ml_runner.set_outdir(workdir)

#     print(f"Training for syst {systName}")

#     if train:
#         params = None
#         if model == 'XGBoost':
#             # params =  {'colsample_bytree': 0.7796308478596582, 'eta': 0.024579670267464545,
#             #           'eval_metric': "logloss", 'gamma': 8.60453282531311, 'max_depth': 4,
#             #           'min_child_weight': 5.0, 'n_estimators': 425, 'subsample': 0.65}
#             # params = {'colsample_bytree': 0.7796308478596582, 'eta': 0.024579670267464545, 'eval_metric': "error", 'gamma': 8.60453282531311, 'max_depth': 4, 'min_child_weight': 5.0, 'n_estimators': 425, 'subsample': 0.5317988260961465}
#             # Three Top full RunII
#             params = {'colsample_bytree': 0.5283131479037887, 'eta': 0.037489501191992,
#                       'eval_metric': 'logloss', 'gamma': 5.439613578244588, 'max_depth': 5,
#                       'min_child_weight': 4.0, 'n_estimators': 250, 'subsample': 0.7009163217893063,
#                       }
#             # params = {'colsample_bytree': 0.6, 'eta': 0.024579670267464545, 'eval_metric': "logloss",
#             #           'gamma': 8.60453282531311, 'max_depth': 3, 'min_child_weight': 5.0,
#             #           'n_estimators': 1500, 'subsample': 0.65}
#             # 4top training
#             # params = {'colsample_bytree': 0.8, 'eta': 0.15, 'eval_metric': "mae",
#             #           'gamma': 2, 'max_depth': 2, 'min_child_weight': 6.0,
#             #           'n_estimators': 1000, 'subsample': 0.9
#             #           }
#             # best fom
#             # params = {"max_depth": 1, "colsample_bytree": 1.0, "min_child_weight": 3.1622776601683795e-05,
#             #           "subsample": 1.0, "eta": 0.05, "eval_metric": "rmse"}
#              # best likelihood
#             # params = {"max_depth": 4, "colsample_bytree": 1.0, "min_child_weight": 0.03162277660168379,
#             #           "subsample": 0.5, "eta": 0.05, "eval_metric": "logloss"}
#         elif model == 'CutBased':
#             pass
#             # mvaRunner.add_cut(mva_params.cuts) # Only used in CutBased for now
#         if params is not None:
#             ml_runner.update_params(params)
#         # Setup test/train/validate sets
#         for year in years:
#             ml_runner.setup_year(workdir, year)
#         ml_runner.train()
#         ml_runner.get_importance()
#         # ml_runner.plot_training_progress()
#         for year in years:
#             ml_runner.apply_model(year, get_auc=True)
#             ml_runner.roc_curve(year)
#             ml_runner.plot_overtrain(year)
#             ml_runner.plot_train_breakdown(year)
#             ml_runner.plot_train_breakdown(year, use_test=False)
#     for year in years:
#         if not train:
#             ml_runner.read_in_file(workdir, year, rerun)
#             if not ml_runner:
#                 return
#             else:
#                 ml_runner.apply_model(year, get_auc=(systName=="Nominal"))
#         ml_runner.output(workdir, year)
def cleanup(cli_args):
    return
