#!/usr/bin/env python3
import numpy as np
import uproot
import sys
from importlib import import_module

from .user import analysis_area
if str(analysis_area) not in sys.path:
    sys.path.append(str(analysis_area))

def get_inputs(workdir, filename='inputs'):
    return import_module(f'.{filename}', f'workspace.{workdir.stem}')

def get_ntuple(name, obj='info'):
    module = import_module(f'.{name}', 'ntuple_info')
    return getattr(module, obj)

def get_ntuple_info(name, obj='info', **kwargs):
    module = import_module(f'.{name}', 'ntuple_info')
    return getattr(module, obj).get_info(**kwargs)

def get_list_systs(infile, tool, systs=["all"], **kwargs):
    allSysts = list()

    if tool == 'flatten':
        from analysis_suite.data.systs import get_change_systs
        allSysts = np.concatenate((['Nominal'], get_change_systs()))
    elif tool == 'mva':
        allSysts = ["_".join(f.stem.split('_')[1:-1]) for f in infile.glob("**/processed*root")]
        allSysts = np.unique(allSysts)
    elif tool == 'combine':
        allSysts = ["_".join(f.stem.split('_')[1:-1]) for f in infile.glob("**/test*root")]
        allSysts = np.unique(allSysts)

    if systs != ['all']:
        finalSysts = list()
        for syst in systs:
            if f'{syst}_up' in allSysts and f'{syst}_down' in allSysts:
                finalSysts += [f'{syst}_up', f'{syst}_down']
            elif syst == "Nominal":
                finalSysts.append("Nominal")
        allSysts = finalSysts
    return allSysts
