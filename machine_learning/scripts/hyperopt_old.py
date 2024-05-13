#!/usr/bin/env python3
import skopt
from collections import OrderedDict
import numpy as np
import pprint
import csv

from analysis_suite.commons.info import GroupInfo
from analysis_suite.commons.configs import get_cli, getGroupDict
import analysis_suite.data.inputs as mva_params
from analysis_suite.BDT_utilities.job_main import get_mva_runner

import sys
sys.argv.insert(1, "mva")


nVars = len(mva_params.usevars)
hyperParams = OrderedDict({
    "hidden_layers": np.array([1, 2]),
    "initial_nodes": np.array([nVars, 2*nVars ]),
    "node_pattern": np.array(["static", "dynamic"]),
    "batch_power": np.array([8, 11]), # 8, 9, 10, 11
    "learning_rate": np.array([0.001, 0.005, 0.01, 0.05]),
    "regulator": np.array(["dropout", "none"]),
    "activation": np.array(["softplus", "elu"]),
    "epochs": np.array([100, 200, 300, 400, 500]),
})
# sqrtNVar = int(np.sqrt(nVars))
# hyperParams = OrderedDict({
#     "max_depth": np.array([1, 2]),
#     'colsample_bytree': np.array([0.5, 0.75, 1.0]),
#     'min_child_weight': np.geomspace(1e-6, 1, 5),
#     'subsample': np.linspace(0.5, 1, 5),
#     'eta': np.array([0.05, 0.1, 0.2, 0.3]),
#     "eval_metric": np.array(["logloss", "rmse", "auc", "error"]),
#     # "objective": np.array(["reg:squaredlogerror", "binary:logistic"]),
# })




def write_line(filename, params, other):
    with open(filename, "a") as logfile:
        writer = csv.writer(logfile)
        if not isinstance(other, list):
            other = list(other)
        writer.writerow(params + other)


# Determine optimization space
def create_opt_space(params):
    opt_space = []
    for key, vals in params.items():
        if vals.dtype.type is np.str_ or len(vals) != 2:
            opt_space.append(skopt.space.Categorical(vals, name=key))
        else:
            opt_space.append(skopt.space.Integer(*vals, name=key))
    return opt_space


# Objective function
def run_mva(mvaRunner, **X):
    print("\n>> Configuration:")
    print(X)

    params = X.update(extraOptions)
    mvaRunner.update_params(X)
    mvaRunner.train(cli_args.workdir)
    mvaRunner.apply_model(cli_args.workdir, "2016")
    fom = mvaRunner.fom["2016"]
    auc = mvaRunner.auc["2016"]
    likely, cut = mvaRunner.approx_likelihood("NBJets", np.linspace(0, 6, 7), "2016")

    print(f"\n>> Obtained ROC-Integral value: {auc}")
    print(f">> Obtained FOM value: {fom}")
    print(f">> Obtained Likelihood value: {likely} at cut {cut}")

    write_line(logfile_name, [X[i] for i in pNames], [mvaRunner.best_iter, fom, cut, likely])

    return -fom


if __name__ == "__main__":
    # user parameters
    cli_args = get_cli()
    work_year = "2016"
    extraOptions = {"verbose": False}
    number_calls = 200

    # Derived Variables
    pNames = list(hyperParams.keys())
    logfile_name = cli_args.workdir / "minimize_results.csv"
    hyper_file = cli_args.workdir / "hyperParams.py"
    group_info = GroupInfo(**vars(cli_args))
    groupDict = getGroupDict(mva_params.groups, group_info)
    opt_space = create_opt_space(hyperParams)

    # Save used parameters to file
    with open( hyper_file, "w" ) as f:
        varName = "hyperParams="
        f.write(varName)
        f.write(pprint.pformat(hyperParams, indent=len(varName)))
        f.write("\n")
        print( "[OK ] Parameters saved to dataset folder." )

    # Start the logfile
    with open(logfile_name,'w') as f: pass # clear file
    write_line(logfile_name, list(hyperParams), ["iters", "AUC"])

    # Setup mva
    mvaRunner = get_mva_runner(cli_args.train)(mva_params.usevars, groupDict)
    for year in cli_args.years:
        mvaRunner.setup_year(cli_args.workdir, year)

    # Perform the optimization
    @skopt.utils.use_named_args(opt_space)
    def objective(**X):
        return run_mva(mvaRunner, **X)

    res_gp = skopt.gp_minimize(
        func = objective,
        dimensions = opt_space,
        n_calls = number_calls,
        n_random_starts = number_calls,
        verbose = True
    )

    result = {key: res_gp.x[i] for i, key in enumerate(pNames)}
    print("TTT DNN Hyper Parameter Optimization Parameters")
    print(f"Static and Parameter Space stored in: {hyper_file}")
    print("Optimized Parameters:")
    for key, val in result.items():
        print(f"    {key}: {val}")
    # Report results
    with open(hyper_file, "a") as f:
        varName = "best_params="
        f.write(varName)
        f.write(pprint.pformat(result, indent=len(varName)))
        f.write("\n")
