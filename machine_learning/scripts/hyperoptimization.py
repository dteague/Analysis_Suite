#!/usr/bin/env python3
import warnings
warnings.filterwarnings('ignore')
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
    # 'n_estimators': hp.quniform('n_estimators', 1, 5, 1),
    'seed': 0,
}

model = None

def objective(space):
    space['max_depth'] = int(space['max_depth'])
    space['n_estimators'] = int(space['n_estimators'])
    model.update_params(space)
    model.train(None)

    for year in all_eras:
        model.apply_model(year)

    space['fom'] = model.get_fom()
    space['auc'] = model.get_auc()
    return {
        'loss': -space['fom'],
        'status': STATUS_OK,
        'auc': space['auc'],
        "fom": space['fom'],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="main", description="")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area/x,
                        help="Working Directory")
    parser.add_argument('-s', '--signal', default='ttt')
    parser.add_argument('-n', '--number_calls', default=100, type=int)
    args = parser.parse_args()
    # user parameters

    # Derived Variables
    inputs = get_inputs(args.workdir)
    input_files = args.workdir/'split_files'

    ginfo = get_ntuple_info('signal')
    samples = ginfo.setup_members()
    if args.signal == 'ttt':
        signal = ginfo.get_members('ttt_nlo')
        outdir = args.workdir / 'hyperopt_ttt'
    else:
        outdir = args.workdir / 'hyperopt_4top'
    outdir.mkdir(exist_ok=True, parents=True)
    logfile = outdir/"trials.pkl"
    hyper_file = outdir/"hyperParams.py"

    model = XGBoostMaker(inputs.usevars, ginfo.get_members('ttt_nlo'),
                         samples)
    model.set_outdir(outdir)
    for year in all_eras:
        model.read_in_files(input_files, year)
    model.read_in_train_files(input_files)

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
