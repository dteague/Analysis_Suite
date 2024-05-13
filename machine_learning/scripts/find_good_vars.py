#!/usr/bin/env python3
import numpy as np
import pandas as pd
import argparse

from analysis_suite.commons.configs import get_list_systs, get_inputs, get_ntuple_info
from analysis_suite.machine_learning.XGBoost import XGBoostMaker
import analysis_suite.commons.user as user


def get_groups(groups, ginfo, additions):
    groupDict = dict()
    for cls, samples in groups.items():
        groupDict[cls] = ginfo.setup_members(filter(lambda s: s in ginfo.group2color, samples))
        if cls in additions:
            add_groups = ginfo.setup_members(filter(lambda s: s in ginfo.group2color, additions[cls]))
            groupDict[cls] = np.concatenate((add_groups, groupDict[cls]))
    return groupDict


def run_train(workdir, years):
    inputs = get_inputs(workdir)
    usevars = inputs.usevars
    groups = inputs.groups

    params = get_inputs(workdir, 'params')

    ginfo = get_ntuple_info('signal')

    # First Training
    groupDict = get_groups(groups, ginfo, {"Signal": ['ttt'], 'Background': ['4top']})
    output_first = workdir/'optimizing'

    model = XGBoostMaker(usevars, groupDict, region='signal')
    model.update_params(params.params_first)
    model.set_outdir(output_first)
    for year in years:
        (output_first/year).mkdir(exist_ok=True, parents=True)
        model.setup_year(workdir, year)

    model.train()
    for imp in ['weight', 'gain', 'cover', 'total_gain', 'total_cover']:
        model.get_importance(imp)
    model.plot_training_progress()
    model.get_corr()
    for year in years:
        model.get_corr(typ='test', year=year)
        model.apply_model(year, get_auc=True)
        model.roc_curve(year)
        model.plot_overtrain(year)
        model.plot_train_breakdown(year)
        model.plot_train_breakdown(year, use_test=False)

        # model.output_train(params.signal_first)

    # for year in years:
    #     model.output(year, params.signal_first)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016pre", "2016post", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area / x,
                        help="Working Directory")
    args = parser.parse_args()

    run_train(args.workdir, args.years)
