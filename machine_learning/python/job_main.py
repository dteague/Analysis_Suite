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
                        cli_args.years, cli_args.region, syst, cli_args.save))
    return argList


def get_groups(groups, ginfo, additions):
    groupDict = dict()
    for cls, samples in groups.items():
        groupDict[cls] = ginfo.setup_members(filter(lambda s: s in ginfo.group2color, samples))
        if cls in additions:
            add_groups = ginfo.setup_members(filter(lambda s: s in ginfo.group2color, additions[cls]))
            groupDict[cls] = np.concatenate((add_groups, groupDict[cls]))
    return groupDict


def training(model, years):
    model.train()
    model.plot_training_progress()
    for year in years:
        model.apply_model(year)
        model.roc_curve(year)
        model.plot_overtrain(year)
        model.plot_train_breakdown(year)
        model.plot_train_breakdown(year, bins=np.linspace(0.8, 1, 41))
        model.plot_train_breakdown(year, use_test=False)


def run(groups, usevars, workdir, model_type, train, years, region, systName, save_train):
    isNom = systName == "Nominal"
    if not train and isNom:
        print("Can only Train Nominal, not apply model: skipping")
        return

    module = import_module(f'.{model_type}', "analysis_suite.machine_learning")
    maker = getattr(module, next(filter(lambda x : "Maker" in x, dir(module))))
    params = get_inputs(workdir, 'params')

    remove_group = ['nonprompt_mc'] if region != "signal" else []
    ginfo = get_ntuple_info(region, remove=remove_group)

    # First Training
    groupDict = get_groups(groups, ginfo, {"Signal": ['4top'], 'Background': ['ttt']})
    output_first = workdir/'first_train'

    model = maker(usevars, groupDict, region=region, systName=systName)
    model.update_params(params.params_first)
    model.set_outdir(output_first)
    for year in years:
        (output_first/year).mkdir(exist_ok=True, parents=True)
        model.setup_year(workdir, year, isNom)

    if train:
        training(model, years)
        model.output_train(params.signal_first)

    for year in years:
        if not train:
            model.apply_model(year, skip_train=True)
        model.output(year, params.signal_first)

    print(f"Training for syst {systName}")

    # Second Training
    groupDict = get_groups(groups, ginfo, {"Signal": ['ttt'], 'NotTrained': ['4top']})
    output_second = workdir/'second_train'

    model2 = maker(usevars, groupDict, region=region, systName=systName)
    model2.update_params(params.params_second)
    model2.set_outdir(output_second)
    model2.add_cut(f'{params.signal_first}<{params.cut}')
    for year in years:
        (output_second/year).mkdir(exist_ok=True, parents=True)
        model2.read_in_file(output_first, year)

    if train:
         training(model2, years)

    for year in years:
        if not train:
            model2.apply_model(year, skip_train=True)
        model2.output(year, params.signal_second)



def cleanup(cli_args):
    return
