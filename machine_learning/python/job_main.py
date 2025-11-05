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

    if cli_args.train and cli_args.region != "signal":
        print("Training not allowed for non signal region")
        exit()

    # Setup splitting (if needed)
    if cli_args.setup_files and cli_args.train:
        print("Reproducing training files split")
        setup_split_files(cli_args.workdir, cli_args.model, cli_args.years,
                          inputs.usevars)

    argList = list()
    if cli_args.train:
        allSysts = ["Nominal"]
    else:
        allSysts = get_list_systs(cli_args.workdir, cli_args.tool, cli_args.systs)
    for syst in allSysts:
        argList.append((inputs.usevars, cli_args.workdir, cli_args.model, cli_args.train,
                        cli_args.years, cli_args.region, syst))
    return argList


def setup_split_files(workdir, model, years, variables):
    module = import_module(f'.{model}', "analysis_suite.machine_learning")
    maker = getattr(module, next(filter(lambda x : "Maker" in x, dir(module))))
    params = get_inputs(workdir, 'params')

    ginfo = get_ntuple_info('signal')
    samples = ginfo.setup_members()

    model = maker(variables, list(), samples)
    model.set_outdir(workdir/'split_files')
    for year in years:
        model.setup_year(workdir, year)
    model.output_files()


def training(model, years):
    model.train()
    model.plot_training_progress()
    for year in years:
        model.apply_model(year)
    print('auc', model.get_auc())
    print('fom', model.get_fom())
    #     model.roc_curve(year)
    #     model.plot_overtrain(year=year)
    #     model.plot_train_breakdown(year)
    #     model.plot_train_breakdown(year, bins=np.linspace(0.8, 1, 41))
    #     model.plot_train_breakdown(year, use_test=False)
    # model.plot_overtrain()

def run(usevars, workdir, model_type, train, years, region, systName):
    isNom = systName == "Nominal"
    isSignal = region == "signal"
    if not train and isNom and isSignal:
        print("Can only Train Nominal, not apply model: skipping")
        return
    elif train and not isSignal:
        print("No training for non signal region")
        return
    elif not isNom and not isSignal:
        print("Can only apply model Nominal for non-signal region")
        return

    input_files = workdir/'split_files'

    module = import_module(f'.{model_type}', "analysis_suite.machine_learning")
    maker = getattr(module, next(filter(lambda x : "Maker" in x, dir(module))))
    params = get_inputs(workdir, 'params')
    ginfo = get_ntuple_info(region)
    samples = ginfo.setup_members()

    # First Training
    output_first = workdir/'first_train'
    bdt_var = params.signal_first
    bdt_file = input_files/f'bdt_{bdt_var}.root'

    signal = ginfo.get_members('4top')
    background = ginfo.get_members('ttt_nlo')
    nontrain = [s for s in samples if s not in signal+background]
    model = maker(usevars, signal, samples, region=region, nonTrained=nontrain,
                  systName=systName)
    model.update_params(params.params_first)
    model.set_outdir(output_first)
    for year in years:
        model.read_in_files(input_files, year)

    if train:
        model.read_in_train_files(input_files)
        training(model, years)
    for year in years:
        if not train:
            model.apply_model(year, skip_train=True)
    model.output_bdt(bdt_file, bdt_var)

    print(f"Training for syst {systName}")
    if params.cut > 1:
        return

    # Second Training
    output_second= workdir/'second_train'
    bdt_var2 = params.signal_second
    nontrained = ['nonprompt', 'charge_flip', 'data', '4top']
    model = maker(usevars, ginfo.get_members('ttt_nlo'), samples, region=region,
                  systName=systName, nonTrained=nontrained)
    model.update_params(params.params_second)
    model.set_outdir(output_second)

    for year in years:
        model.read_in_files(input_files, year)
    if train:
        model.read_in_train_files(input_files)
    model.read_in_bdt(bdt_file, bdt_var, onlyTest=not train)
    model.mask(lambda df: df[bdt_var] < params.cut)

    if train:
        training(model, years)
    for year in years:
        if not train:
            model.apply_model(year, skip_train=True)
    model.output_bdt(input_files/f'bdt_{bdt_var2}.root', bdt_var2)


def cleanup(cli_args):
    return
