#!/usr/bin/env python3
import pandas as pd
import numpy as np
import xgboost as xgb
import argparse
import pprint
import pickle
from sklearn.metrics import accuracy_score
from hyperopt import STATUS_OK, Trials, fmin, hp, tpe

from analysis_suite.commons.configs import get_list_systs, get_inputs, get_ntuple_info
import analysis_suite.commons.user as user
from analysis_suite.commons.constants import all_eras
from analysis_suite.machine_learning.XGBoost import XGBoostMaker


space = {
    'max_depth': hp.quniform("max_depth", 1, 5, 1),
    'gamma': hp.uniform ('gamma', 1,9),
    "eta": hp.loguniform("eta", np.log(0.001), np.log(0.5)),
    "eval_metric": hp.choice("eval_metric", ["logloss", "rmse", "mae",  "auc"]), # error
    "colsample_bytree": hp.uniform("colsample_bytree", 0.5, 1),
    'min_child_weight' : hp.quniform('min_child_weight', 0, 10, 1),
    "subsample": hp.uniform("subsample", 0.50, 1.0),
    'n_estimators': hp.quniform('n_estimators', 100, 1000, 50),
    'seed': 0,
}

model = None

def objective(space):
    space['max_depth'] = int(space['max_depth'])
    space['n_estimators'] = int(space['n_estimators'])
    model.update_params(space)
    model.train(None)

    total_fom = 0
    min_auc = 1.
    for year in all_eras:
        model.apply_model(year)
        total_fom += model.get_fom(year)**2
        min_auc = min(min_auc, model.get_auc(year))
    total_fom = np.sqrt(total_fom)
    space['fom'] = total_fom
    space['auc'] = min_auc
    return {
        'loss': -total_fom,
        'status': STATUS_OK,
        'auc': min_auc,
        "fom": total_fom,
    }

def get_groups(groups, ginfo, additions):
    groupDict = dict()
    for cls, samples in groups.items():
        groupDict[cls] = ginfo.setup_members(filter(lambda s: s in ginfo.group2color, samples))
        if cls in additions:
            add_groups = ginfo.setup_members(filter(lambda s: s in ginfo.group2color, additions[cls]))
            groupDict[cls] = np.concatenate((add_groups, groupDict[cls]))
    return groupDict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="main", description="")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area/x,
                        help="Working Directory")
    parser.add_argument('-s', '--signal', default='ttt')
    parser.add_argument('-n', '--number_calls', default=100, type=int)
    args = parser.parse_args()
    # user parameters

    # Derived Variables
    if args.signal == 'ttt':
        group_add = {"Signal": ['ttt'], 'NotTrained': ['4top']}
        outdir = args.workdir / 'hyperopt_ttt'
    else:
        group_add = {"Signal": ['4top'], 'Background': ['ttt']}
        outdir = args.workdir / 'hyperopt_4top'
    outdir.mkdir(exist_ok=True, parents=True)

    logfile = outdir/"trials.pkl"
    hyper_file = outdir/"hyperParams.py"

    ginfo = get_ntuple_info('signal')
    inputs = get_inputs(args.workdir)
    groupDict = get_groups(inputs.groups, ginfo, group_add)
    model = XGBoostMaker(inputs.usevars, groupDict, region='signal')
    model.set_outdir(outdir)
    if args.signal == 'ttt':
        model.add_cut(f'4top_sig<0.925')
        for year in all_eras:
            # model.setup_year(args.workdir, year)
            model.read_in_file(args.workdir/'first_train', year)

    trials = pickle.load(open(logfile, "rb")) if logfile.exists() else Trials()

    best_hyperparams = fmin(fn = objective,
                            space = space,
                            algo = tpe.suggest,
                            max_evals = args.number_calls,
                            trials = trials)

    pickle.dump(trials, open(logfile, "wb"))
    print(best_hyperparams)

    # Save used parameters to file
    with open( hyper_file, "w" ) as f:
        varName = "hyperParams="
        f.write(varName)
        f.write(pprint.pformat(best_hyperparams, indent=len(varName)))
        f.write("\n")
        print( "[OK ] Parameters saved to dataset folder." )

