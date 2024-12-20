#!/usr/bin/env python3
import pandas as pd
import logging
import uproot
import subprocess

from analysis_suite.commons.configs import get_list_systs, get_inputs, get_ntuple
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.commons.constants import lumi
import analysis_suite.commons.user as user

# from .data_processor import DataProcessor


def setup(cli_args):
    argList = list()

    ntuple = get_ntuple(cli_args.ntuple)
    groups = []
    for g, members in ntuple.get_info(keep_dd_data=True).setup_groups().items():
        for m in members:
            groups.append([g, m])

    for year in cli_args.years:
        outdir = cli_args.workdir /year
        outdir.mkdir(exist_ok=True, parents=True)
        allSysts = get_list_systs(ntuple.get_filename(year=year), cli_args.tool, cli_args.systs)
        for syst in allSysts:
            # with uproot.recreate(outdir / f'processed_{syst}_{cli_args.ntuple}.root') as f:
            #     pass

            for group in groups:
                argList.append((cli_args.workdir, outdir, cli_args.ntuple, year, syst, cli_args.inputs, {group[0]: [group[1]]}))
    return argList


def run(workdir, outdir, tupleName, year, syst, inputs, group):
    use_group = next(iter(group))
    use_mem = group[use_group][0]
    if use_mem == 'data' and use_group != 'data':
        use_mem = use_group
    if (outdir / f'{syst}_{use_mem}.root').exists():
        return

    inputs = get_inputs(workdir)
    ntuple = get_ntuple(tupleName)
    plotter = Plotter(ntuple.get_filename(year=year),
                      group, ntuple=ntuple, systName=syst, year=year)
    final_set = pd.DataFrame()
    final_weight = pd.DataFrame()

    allvars = inputs.allvar
    if "CR" in tupleName and syst != "Nominal":
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
            final_set = pd.concat([final_set, df], ignore_index=True)
            final_weight = pd.concat([final_weight, pd.DataFrame.from_dict(vg.get_all_weights())], ignore_index=True)

    if final_set.empty:
        print(f"{use_mem} skipped")
        return

    with uproot.recreate(outdir / f'{syst}_{use_mem}.root') as f:
        f[use_mem] = final_set
        f[f"weights/{use_mem}"] = final_weight
    print(f"Finished {use_mem} for {syst} for region {tupleName}")


def cleanup(cli_args):
    ntuple = get_ntuple(cli_args.ntuple)
    for year in cli_args.years:
        outdir = cli_args.workdir /year
        allSysts = get_list_systs(ntuple.get_filename(year=year), cli_args.tool, cli_args.systs)
        for syst in allSysts:
            output_file = outdir / f'processed_{syst}_{cli_args.ntuple}.root'
            input_files = ""
            for f in (outdir.glob(f"{syst}*root")):
                input_files += f" {f}"
            # print(f"hadd -f -v 1 -j {cli_args.j} {output_file} {input_files}")
            subprocess.call(f"hadd -f -v 1 -j 5 {output_file} {input_files}", shell=True)
            for f in (outdir.glob(f"{syst}*root")):
                f.unlink()
