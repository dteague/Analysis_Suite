#!/usr/bin/env python3
import pandas as pd
import numpy as np
import logging
import uproot
import subprocess

from analysis_suite.commons.configs import get_list_systs, get_inputs, get_ntuple
from analysis_suite.commons.constants import lumi
import analysis_suite.commons.user as user

from analysis_suite.flatten import NtupleGetter
from analysis_suite.commons.info import fileInfo


def setup(cli_args):
    argList = list()

    ntuple = get_ntuple(cli_args.ntuple)
    groups = []
    for year in cli_args.years:
        outdir = cli_args.workdir/year
        outdir.mkdir(exist_ok=True, parents=True)
        allSysts = get_list_systs(ntuple.get_filename(year=year), cli_args.tool, cli_args.systs)
        for syst in allSysts:
            argList.append((cli_args.workdir, cli_args.ntuple, year, syst))
    return argList


def run(workdir, tupleName, year, syst):
    print(year, syst, "start")
    output = workdir / year/ f'processed_{syst}_{tupleName}.root'
    # if output.exists():
    #     return
    ntuple = get_ntuple(tupleName)
    filename = ntuple.get_filename(year)
    inputs = get_inputs(workdir)
    final_set = {}
    final_weight = {}
    final_ratio = {}

    allvars = inputs.allvar
    # if "CR" in tupleName and syst != 'Nominal':
    #     allvars = {
    #         "HT": inputs.allvar['HT'],
    #         "NJets": inputs.allvar['NJets'],
    #         "NLeps": lambda vg: vg['TightLepton'].num(),
    #     }

    for root_file in filename.glob("*root"):
        f = uproot.open(root_file)
        members = [m for m in f.keys(recursive=False, cycle=False)]
        for member in members:
            xsec = 1. if fileInfo.is_data(member) else fileInfo.get_xsec(member)*lumi[year]*1000
            for tree in ntuple.trees:
                group = ntuple.get_group_name(member, tree)
                outname = group if member == 'data' else member
                if group is None:
                    continue
                if group == "nonprompt_mc" and syst != 'Nominal':
                    continue # Skip since only used in training
                vg = NtupleGetter(f, tree, member, xsec, systName=syst, cuts=ntuple.cut)
                if not vg.tree or not vg.correct_syst:
                    continue
                ntuple.setup_branches(vg)
                vg.set_systematic(syst)
                if not vg or not vg.correct_syst:
                    continue

                # print(outname, group, len(vg))
                df_dict = {}
                for varname, func in allvars.items():
                    df_dict[varname] = func(vg)
                df_dict["scale_factor"] = vg.scale
                df = pd.DataFrame.from_dict(df_dict)
                df = df.astype({col: int for col in df.columns if col[0] == 'N'})
                if outname not in final_set:
                    final_set[outname] = pd.DataFrame()
                    final_weight[outname] = pd.DataFrame()
                    final_ratio[outname] = pd.DataFrame()

                final_set[outname] = pd.concat([final_set[outname], df], ignore_index=True)
                final_weight[outname] = pd.concat([final_weight[outname], pd.DataFrame.from_dict(vg.get_all_weights())],
                                                  ignore_index=True)
                if syst != "Nominal":
                    final_ratio[outname] = pd.concat([final_ratio[outname], pd.DataFrame.from_dict({syst:vg.scale/vg.get_nom()})], ignore_index=True)

    with uproot.recreate(output) as f:
        for outname in final_set:
            f[outname] = final_set[outname]
            f[f"weights/{outname}"] = final_weight[outname]
            if syst != "Nominal":
                f[f"ratio/{outname}"] = final_ratio[outname]

            # f[f"/{outname}"] = final_weight[outname]
    print(year, syst, "finished")

def cleanup(cli_args):
    pass
