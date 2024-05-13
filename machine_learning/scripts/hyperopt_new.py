#!/usr/bin/env python3
import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.metrics import accuracy_score
from hyperopt import STATUS_OK, Trials, fmin, hp, tpe

from analysis_suite.commons.info import GroupInfo
from analysis_suite.commons.configs import get_list_systs, get_inputs
import analysis_suite.commons.user as user
from analysis_suite.machine_learning.XGBoost import XGBoostMaker as maker

space = {
    'max_depth': hp.quniform("max_depth", 1, 5, 1),
    'gamma': hp.uniform ('gamma', 1,9),
    "eta": hp.loguniform("eta", np.log(0.001), np.log(0.5)),
    "eval_metric": hp.choice("eval_metric", ["logloss", "rmse", "mae",  "auc"]), # error
    "colsample_bytree": hp.uniform("colsample_bytree", 0.5, 1),
    'min_child_weight' : hp.quniform('min_child_weight', 0, 10, 1),
    "subsample": hp.uniform("subsample", 0.50, 1.0),
    'n_estimators': hp.quniform("n_estimators", 100, 1000, 50),
    'seed': 0,
}

workdir = user.workspace_area/"June_fullrun"
inputs = get_inputs(workdir)
ginfo = GroupInfo(inputs.color_by_group)
groupDict = {cls_name: ginfo.setup_members(samples) for cls_name, samples in inputs.groups.items()}
ml_runner = maker(inputs.usevars, groupDict)
for year in ['2016pre', '2016post', '2017', '2018']:
    ml_runner.setup_year(workdir, year)


def objective(space):
    space['max_depth'] = int(space['max_depth'])
    space['n_estimators'] = int(space['n_estimators'])
    x_train = ml_runner.train_set.drop(ml_runner._drop_vars, axis=1)
    y_train = ml_runner.train_set.classID
    w_train = ml_runner.train_set.split_weight.copy()
    s_train = abs(ml_runner.train_set.scale_factor.copy())

    x_test = ml_runner.validation_set.drop(ml_runner._drop_vars, axis=1)
    y_test = ml_runner.validation_set.classID
    w_test = ml_runner.validation_set.split_weight.copy()
    s_test = abs(ml_runner.validation_set.scale_factor.copy())

    clf=xgb.XGBClassifier(**space, use_label_encoder=False)
    evaluation = [(x_train, y_train), (x_test, y_test)]

    def fom_metric(y_pred, dtrain):
        prev = fom_metric.prev
        nbins = 20
        bins = np.linspace(0, 1, nbins+1)

        y_true = dtrain.get_label()
        y_pred = 1/(1+np.exp(-y_pred))
        weight = dtrain.get_weight()
        b = np.histogram(y_pred[y_true==0], bins, weights=weight[y_true==0])[0]
        s = np.histogram(y_pred[y_true==1], bins, weights=weight[y_true==1])[0]
        fom = -np.sqrt(2*np.sum((s+b)*np.log(1+s/(b+1e-5))-s))
        # print(min(y_pred), max(y_pred))

        diff = abs((fom-prev)/fom)
        if diff > 0.2 and prev < -0.09:
            fom = prev
        fom_metric.prev = fom
        return 'fom', fom
    fom_metric.prev = 0

    clf.fit(x_train, y_train, sample_weight=w_train, eval_metric=fom_metric,
            eval_set=evaluation, sample_weight_eval_set=[w_train, w_test],
            early_stopping_rounds=100, verbose=False)

    pred = clf.predict_proba(x_test).T[0]
    accuracy = accuracy_score(y_test, pred>0.5, sample_weight=w_test)
    bins = np.linspace(0, 1, 16)
    s = np.histogram(pred[y_test==1], bins, weights=s_test[y_test==1])[0]
    b = np.histogram(pred[y_test==0], bins, weights=s_test[y_test==0])[0]
    fom = np.sqrt(2*np.sum((s+b)*np.log(1+s/(b+1e-5))-s))
    space['fom'] = fom
    space['accuracy'] = accuracy
    write_line("hyperopt.log", space)
    return {'loss': -fom, 'status': STATUS_OK }

def write_line(filename, data):
    with open(filename, "a") as logfile:
        logfile.write(f'{data},\n')



if __name__ == "__main__":
    # user parameters
    work_year = "2017"
    extraOptions = {"verbose": False}
    number_calls = 500

    # # Derived Variables
    # pNames = list(hyperParams.keys())
    # logfile_name = cli_args.workdir / "minimize_results.csv"
    # hyper_file = cli_args.workdir / "hyperParams.py"
    # group_info = GroupInfo(**vars(cli_args))
    # groupDict = getGroupDict(mva_params.groups, group_info)
    # opt_space = create_opt_space(hyperParams)

    trials = Trials()

    best_hyperparams = fmin(fn = objective,
                            space = space,
                            algo = tpe.suggest,
                            max_evals = number_calls,
                            trials = trials)
    print(best_hyperparams)

    # # Save used parameters to file
    # with open( hyper_file, "w" ) as f:
    #     varName = "hyperParams="
    #     f.write(varName)
    #     f.write(pprint.pformat(hyperParams, indent=len(varName)))
    #     f.write("\n")
    #     print( "[OK ] Parameters saved to dataset folder." )

    # # Start the logfile
    # with open(logfile_name,'w') as f: pass # clear file
    # write_line(logfile_name, list(hyperParams), ["iters", "AUC"])

    # # Setup mva
    # mvaRunner = get_mva_runner(cli_args.train)(mva_params.usevars, groupDict)
    # for year in cli_args.years:
    #     mvaRunner.setup_year(cli_args.workdir, year)

    # # Perform the optimization
    # @skopt.utils.use_named_args(opt_space)
    # def objective(**X):
    #     return run_mva(mvaRunner, **X)

    # res_gp = skopt.gp_minimize(
    #     func = objective,
    #     dimensions = opt_space,
    #     n_calls = number_calls,
    #     n_random_starts = number_calls,
    #     verbose = True
    # )

    # result = {key: res_gp.x[i] for i, key in enumerate(pNames)}
    # print("TTT DNN Hyper Parameter Optimization Parameters")
    # print(f"Static and Parameter Space stored in: {hyper_file}")
    # print("Optimized Parameters:")
    # for key, val in result.items():
    #     print(f"    {key}: {val}")
    # # Report results
    # with open(hyper_file, "a") as f:
    #     varName = "best_params="
    #     f.write(varName)
    #     f.write(pprint.pformat(result, indent=len(varName)))
    #     f.write("\n")
