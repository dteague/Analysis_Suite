#!/usr/bin/env python3
import logging
import numpy as np
import os
from importlib import import_module

from analysis_suite.commons.configs import get_list_systs, get_inputs, get_ntuple_info
import analysis_suite.commons.user as user

def setup(cli_args):
    inputs = get_inputs(cli_args.workdir)
    remove_group = inputs.remove if cli_args.region != "Signal" else []
    ginfo = get_ntuple_info(cli_args.region, remove=remove_group)
    os.environ["NUMEXPR_MAX_THREADS"] = "8"

    argList = list()
    groupDict = dict()
    if cli_args.train:
        allSysts = ["Nominal"]
        for cls, samples in inputs.groups.items():
            groupDict[cls] = ginfo.setup_members(filter(lambda s: s in ginfo.group2color, samples))
    else:
        allSysts = get_list_systs(cli_args.workdir, cli_args.tool, cli_args.systs)
        groupDict['Signal'] = ginfo.setup_members(['ttt'])
        groupDict['Background'] = [mem for mem in ginfo.setup_members() if mem not in groupDict['Signal']]

    for syst in allSysts:
        argList.append((groupDict, inputs.usevars, cli_args.workdir, cli_args.model, cli_args.train,
                        cli_args.years, cli_args.region, syst, cli_args.save, cli_args.rerun))

    return argList


def run(groupDict, usevars, workdir, model, train, years, region, systName, save_train, rerun):
    module = import_module(f'.{model}', "analysis_suite.machine_learning")
    maker = getattr(module, next(filter(lambda x : "Maker" in x, dir(module))))
    ml_runner = maker(usevars, groupDict, region=region, systName=systName)

    print(f"Training for syst {systName}")

    if train:
        params = None
        if model == 'XGBoost':
            # params =  {'colsample_bytree': 0.7796308478596582, 'eta': 0.024579670267464545,
            #           'eval_metric': "logloss", 'gamma': 8.60453282531311, 'max_depth': 4,
            #           'min_child_weight': 5.0, 'n_estimators': 425, 'subsample': 0.65}
            # params = {'colsample_bytree': 0.7796308478596582, 'eta': 0.024579670267464545, 'eval_metric': "error", 'gamma': 8.60453282531311, 'max_depth': 4, 'min_child_weight': 5.0, 'n_estimators': 425, 'subsample': 0.5317988260961465}
            # Three Top full RunII
            params = {'colsample_bytree': 0.5283131479037887, 'eta': 0.037489501191992,
                      'eval_metric': 'logloss', 'gamma': 5.439613578244588, 'max_depth': 5,
                      'min_child_weight': 4.0, 'n_estimators': 950, 'subsample': 0.7009163217893063,
                      }
            # params = {'colsample_bytree': 0.6, 'eta': 0.024579670267464545, 'eval_metric': "logloss",
            #           'gamma': 8.60453282531311, 'max_depth': 3, 'min_child_weight': 5.0,
            #           'n_estimators': 1500, 'subsample': 0.65}
            # 4top training
            # params = {'colsample_bytree': 0.8, 'eta': 0.15, 'eval_metric': "mae",
            #           'gamma': 2, 'max_depth': 2, 'min_child_weight': 6.0,
            #           'n_estimators': 1000, 'subsample': 0.9
            #           }
            # best fom
            # params = {"max_depth": 1, "colsample_bytree": 1.0, "min_child_weight": 3.1622776601683795e-05,
            #           "subsample": 1.0, "eta": 0.05, "eval_metric": "rmse"}
             # best likelihood
            # params = {"max_depth": 4, "colsample_bytree": 1.0, "min_child_weight": 0.03162277660168379,
            #           "subsample": 0.5, "eta": 0.05, "eval_metric": "logloss"}
        elif model == 'CutBased':
            pass
            # mvaRunner.add_cut(mva_params.cuts) # Only used in CutBased for now
        if params is not None:
            ml_runner.update_params(params)
        # Setup test/train/validate sets
        for year in years:
            ml_runner.setup_year(workdir, year, save_train)
        ml_runner.train(workdir)
        ml_runner.get_importance(workdir)
        ml_runner.plot_training_progress(workdir)
        for year in years:
            ml_runner.roc_curve(workdir, year)
            ml_runner.overtrain_test(workdir, year)

    for year in years:
        if not train:
            ml_runner.read_in_file(workdir, year, rerun)
            if not ml_runner: return
        # Apply to test sets and save
        ml_runner.apply_model(workdir, year, get_auc=(systName=="Nominal"))
        ml_runner.output(workdir, year)

def cleanup(cli_args):
    return
