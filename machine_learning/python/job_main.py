#!/usr/bin/env python3
import logging
import numpy as np
import os
from importlib import import_module
from copy import copy
import pandas as pd
import uproot

from analysis_suite.commons.configs import get_list_systs, get_inputs, get_ntuple_info, get_ntuple
import analysis_suite.commons.user as user
from analysis_suite.plotting.plotter import Plotter

def setup(cli_args):
    inputs = get_inputs(cli_args.workdir)
    os.environ["NUMEXPR_MAX_THREADS"] = "8"

    if cli_args.train and cli_args.ntuple != "signal":
        print("Training not allowed for non signal region")
        exit()

    # Setup splitting (if needed)
    if cli_args.setup_files or not (cli_args.workdir/'split_files').exists():
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
                        cli_args.years, cli_args.ntuple, syst))
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


def plot(model, plotter, years):
    # Plots
    model.plot_training_progress()
    plotter.set_year('all')
    hists = plotter.fill_hist_groups(model.get_hist(25, useTrain=True))
    plotter.plot_hist('BDT_train', hists, plot_type='stack_shape',
                      plot_sigs=["4top", 'ttt_nlo'])
    for year in years:
        plotter.set_year(year)
        hists = plotter.fill_hist_groups(model.get_hist(25, year))
        plotter.plot_hist('BDT_train', hists, plot_type='stack_shape',
                          plot_sigs=["4top", 'ttt_nlo'])
        model.plot_roc(year)
        model.plot_overtrain(year)
    print('auc', model.get_auc())
    print('fom', model.get_fom())

def run(usevars, workdir, model_type, train, years, ntuple, systName, blind=True):
    params = get_inputs(workdir, 'params')
    ginfo = get_ntuple_info(ntuple)
    samples = ginfo.setup_members()

    isNom = systName == "Nominal"
    isSignal = ntuple == "signal"
    out_syst = systName
    if not train and isNom and isSignal:
        print("Can only Train Nominal, not apply model: skipping")
        return
    elif train and not isSignal:
        print("No training for non signal region")
        return
    # elif not isNom and not isSignal:
    #     print("Can only apply model Nominal for non-signal region")
    #     return

    module = import_module(f'.{model_type}', "analysis_suite.machine_learning")
    maker = getattr(module, next(filter(lambda x : "Maker" in x, dir(module))))

    # First Training
    output_first = workdir/'first_train'
    bdt_var = params.signal_first
    bdt_dir = workdir/'split_files'
    bdt_file = bdt_dir/f'bdt_{bdt_var}_{out_syst}_{ntuple}.root'
    plotter = Plotter(get_ntuple('bdt'), 'all', bkg='all', outdir=output_first)
    if train or (isNom and isSignal):
        in_type = 'test'
        input_files = workdir/'split_files'
    else:
        in_type = 'processed'
        input_files = workdir

    signal = ginfo.get_members('4top')
    model = maker(usevars, signal, samples, region=ntuple,
                  systName=systName)
    model.update_params(params.params_first)
    model.set_outdir(output_first)
    for year in years:
        model.read_in_files(input_files, year, typ=in_type)

    if train:
        model.read_in_train_files(input_files)
        model.train()
    model.apply_model(years, skip_train=not train)

    if train:
        plot(model, plotter, years)
    # model.output_all(workdir/f'all_{bdt_var}_Nominal_{ntuple}.root')
    model.output_bdt(bdt_file, bdt_var)

    print(f"Training for syst {systName}")

    # Second Training
    output_second= workdir/'second_train'
    bdt_var2 = params.signal_second
    nontrained = ['nonprompt', 'charge_flip', 'data', '4top']
    plotter.sig = 'ttt_nlo'
    plotter.outdir = output_second

    model = maker(usevars, ginfo.get_members('ttt_nlo'), samples, region=ntuple,
                  systName=systName, nonTrained=nontrained)
    model.update_params(params.params_second)
    model.set_outdir(output_second)

    for year in years:
        model.read_in_files(input_files, year, typ=in_type)
    if train:
        model.read_in_train_files(input_files)
    if ntuple == 'signal':
        model.read_in_bdt(bdt_file, bdt_var, onlyTest=not train)
        print(bdt_var, params.cut[year])
        model.mask(lambda df, year: df[bdt_var] < params.cut[year])

    if train:
        model.read_in_train_files(input_files)
        model.train()
    model.apply_model(years, skip_train=not train)

    # Plots
    if train:
        plot(model, plotter, years)
    # model.output_all(workdir/f'all_{bdt_var2}_Nominal_{ntuple}.root', bdt_var2)
    model.output_bdt(bdt_dir/f'bdt_{bdt_var2}_{out_syst}_{ntuple}.root', bdt_var2)


def cleanup(cli_args):
    return
