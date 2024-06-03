#!/usr/bin/env python3
import pandas as pd
import logging
import uproot

from analysis_suite.commons.configs import get_list_systs, get_inputs, get_ntuple
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.commons.constants import lumi
import analysis_suite.commons.user as user

# from .data_processor import DataProcessor

def setup(cli_args):
    argList = list()

    ntuple = get_ntuple(cli_args.ntuple)
    for year in cli_args.years:
        outdir = cli_args.workdir /year
        outdir.mkdir(exist_ok=True, parents=True)
        allSysts = get_list_systs(ntuple.get_filename(year=year), cli_args.tool, cli_args.systs)
        for syst in allSysts:
            argList.append((cli_args.workdir, outdir, cli_args.ntuple, year, syst, cli_args.inputs))
    return argList


def run(workdir, outdir, tupleName, year, syst, inputs):
    inputs = get_inputs(workdir)
    ntuple = get_ntuple(tupleName)
    groups = ntuple.get_info(keep_dd_data=True).setup_groups()
    plotter = Plotter(ntuple.get_filename(year=year),
                      groups, ntuple=ntuple, systName=syst, year=year)
    final_set = dict()
    final_weight = dict()

    allvars = inputs.allvar
    if "CR" in tupleName:
        allvars = {
            "HT": inputs.allvar['HT'],
            "NJets": inputs.allvar['NJets'],
        }

    for member, vgs in plotter.dfs.items():
        for tree, vg in vgs.items():
            df_dict = {varname: func(vg) for varname, func in allvars.items()}
            df_dict["scale_factor"] = vg.scale
            df = pd.DataFrame.from_dict(df_dict)
            df = df.astype({col: int for col in df.columns if col[0] == 'N'})
            if member not in final_set:
                final_set[member] = df
                final_weight[member] = pd.DataFrame.from_dict(vg.get_all_weights())
            else:
                final_set[member] = pd.concat([final_set[member], df], ignore_index=True)
                final_weight[member] = pd.concat([final_weight[member], pd.DataFrame.from_dict(vg.get_all_weights())], ignore_index=True)
            # print(f"Finished setting up {member} in tree {tree}")

    with uproot.recreate(outdir / f'processed_{syst}_{tupleName}.root') as f:
        for group, df in final_set.items():
            if not len(df):
                continue
            f[group] = df
            f[f"weights/{group}"] = final_weight[group]
    print(f"Finished {syst} for region {tupleName}")

def cleanup(cli_args):
    pass
