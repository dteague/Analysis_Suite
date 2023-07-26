#!/usr/bin/env python3
import numpy as np
import uproot
import sys
from importlib import import_module

from .user import analysis_area
if str(analysis_area) not in sys.path:
    sys.path.append(str(analysis_area))

def get_graph(name, graphs):
    for graph in graphs:
        if graph.name == name:
            return graph
    raise Exception(f"Graph ({name}) has not been found")
    
def get_inputs(workdir):
    return import_module('.inputs', f'workspace.{workdir.stem}')

def get_ntuple(name, obj='info'):
    module = import_module(f'.{name}', 'ntuple_info')
    return getattr(module, obj)

def get_ntuple_info(name, obj='info', **kwargs):
    module = import_module(f'.{name}', 'ntuple_info')
    return getattr(module, obj).get_info(**kwargs)

def get_list_systs(infile, tool, systs=["all"], **kwargs):
    allSysts = list()

    if tool == 'flatten':
        if infile.is_dir():
            all_files = list(infile.glob('*root'))
        else:
            all_files = [infile]

        for file_ in all_files:
            with uproot.open(file_) as f:
                for group in [k for k in f.keys() if "/" not in k]:
                    all_systs = [name.member("fName") for name in f[group]['Systematics']]
                    good_systs = {}
                    for name in f[group]["Syst_Index"]:
                        novelNum = int(name.member("fTitle"))
                        if novelNum not in good_systs:
                            good_systs[novelNum] = all_systs[int(name.member("fName"))]
                    allSysts = np.concatenate((allSysts, list(good_systs.values())))
                    allSysts = np.unique(allSysts)
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

def sig_fig(x, p=3):
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags

@np.vectorize
def asymptotic_sig(s, b):
    return s/np.sqrt(b+1e-5)

