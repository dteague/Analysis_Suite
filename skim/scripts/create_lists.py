#!/usr/bin/env python3
import subprocess
import numpy as np
from pathlib import Path
import argparse
import itertools
import re

from analysis_suite.commons.info import GroupInfo, fileInfo
from analysis_suite.commons.user import analysis_area, xrd_tag

condition_dict = {
    "data" : {
        "2016" : [
            "Run2016B-ver2*UL2016_MiniAODv2_NanoAODv9-v*",
            "Run2016C*UL2016_MiniAODv2_NanoAODv9-v*",
            "Run2016D*UL2016_MiniAODv2_NanoAODv9-v*",
            "Run2016E*UL2016_MiniAODv2_NanoAODv9-v*",
            "Run2016F*UL2016_MiniAODv2_NanoAODv9-v*",
            "Run2016G*UL2016_MiniAODv2_NanoAODv9-v*",
            "Run2016H*UL2016_MiniAODv2_NanoAODv9-v*",
        ],
        "2016pre" : [
            "Run2016B-ver2*UL2016_MiniAODv2_NanoAODv9-v*",
            "Run2016*-HIPM_UL2016_MiniAODv2_NanoAODv9-v*",
        ],
        "2016post" : [
            "Run2016*-UL2016_MiniAODv2_NanoAODv9-v*",
        ],
        "2017" : [
            "Run2017B-UL2017_MiniAODv2_NanoAODv9-v*",
            "Run2017C-UL2017_MiniAODv2_NanoAODv9-v*",
            "Run2017D-UL2017_MiniAODv2_NanoAODv9-v*",
            "Run2017E-UL2017_MiniAODv2_NanoAODv9-v*",
            "Run2017F-UL2017_MiniAODv2_NanoAODv9-v*",
        ],
        "2018": [
            "Run2018A-UL2018_MiniAODv2_NanoAODv9-v*",
            "Run2018B-UL2018_MiniAODv2_NanoAODv9-v*",
            "Run2018C-UL2018_MiniAODv2_NanoAODv9-v*",
            "Run2018D-UL2018_MiniAODv2_NanoAODv9-v*",
        ],
    },
    "mc" : {
        "2016pre":  ["RunIISummer20UL16NanoAODAPVv9*-106X_mcRun2_asymptotic_preVFP*v*"],
        "2016post": ["RunIISummer20UL16NanoAODv9*-106X_mcRun2_asymptotic_v*"],
        "2017":     ["RunIISummer20UL17NanoAODv9*-106X_mc2017_realistic_v*"],
        "2018":     ["RunIISummer20UL18NanoAODv9*-106X_upgrade2018_realistic_v*"],
    },

}

dataset_dict = {
    "data" : {
        "MM" : {
            '2016' : "DoubleMuon",
            '2017' : "DoubleMuon",
            '2018' : "DoubleMuon",
        },
        "EM" : {
            '2016' : "MuonEG",
            '2017' : "MuonEG",
            '2018' : "MuonEG",
        },
        "EE" : {
            '2016' : "DoubleEG",
            '2017' : "DoubleEG",
            '2018' : "EGamma",
        },
        "E" : {
            '2016' : 'DoubleEG',
            '2017' : 'SingleElectron',
            '2018' : 'EGamma',
        },
        "Met": {
            "2016": "MET",
            "2017": "MET",
            "2018": "MET",
        }
    },

    "mc" : {
        'ttt': [
            'TTTJ_TuneCP5*',
            'TTTW_TuneCP5*',
        ],
        'tttt' : [
            'TTTT_TuneCP5*',
        ],
        'ttX' : [
            'TTWJetsToLNu_TuneCP5_13TeV*',
            'TTZToLLNuNu_M-10_TuneCP5_13TeV*',
            'TTZToLL_M-1to10_TuneCP5*',
            'ttHToNonbb_M125_TuneCP5*',
            'ttHTobb_M125_TuneCP5*',
        ],
        'ttXY' : [
            'TTHH_TuneCP5*',
            'TTWH_TuneCP5*',
            'TTWW_TuneCP5*',
            'TTWZ_TuneCP5*',
            'TTZH_TuneCP5*',
            'TTZZ_TuneCP5*',
        ],
        'xg' : [
            'TGJets_TuneCP5*',
            'TTGamma_Dilept_TuneCP5_13TeV*',
            'TTGamma_Hadronic_TuneCP5*',
            'TTGamma_SingleLept_TuneCP5_13TeV*',
            'WGToLNuG_TuneCP5*',
            'WWG_TuneCP5*',
            'WZG_TuneCP5*',
            'ZGToLLG_01J_5f_TuneCP5*',
        ],
        'vv_inc' : [
            'WW_TuneCP5*',
            'WZ_TuneCP5*',
            'ZZ_TuneCP5*',
        ],
        'vvv' : [
            'WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8',
            'WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8',
            'WZZ_TuneCP5_13TeV-amcatnlo-pythia8',
            'ZZZ_TuneCP5_13TeV-amcatnlo-pythia8',
        ],
        'qcd' : [
            "QCD_Pt-20_MuEnrichedPt15_Tune*",
            "QCD_Pt-*_MuEnrichedPt5_Tune*",
            "QCD_Pt-*_EMEnriched_Tune*",
            "QCD_Pt_*_bcToE_Tune*",
        ],
        'ttbar' : [
            'TT_TuneCH3_13TeV-powheg-herwig7',
        ],
        'ttbar_lep' : [
            'TTTo2L2Nu_TuneCP5_13TeV*',
            'TTToSemiLeptonic_TuneCP5_13TeV*',
            'TTToHadronic_TuneCP5_13TeV*',
        ],
        'dy' : [
            'DYJetsToLL_M-10to50_TuneCP5*',
            'DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8',
        ],
        'dy_ht' : [
            "DYJetsToLL_M-50_HT-*_TuneCP5*13TeV-madgraphMLM-pythia8",
        ],
        'wjet' : [
            'WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8',
        ],
        'wjet_ht' : [
            'WJetsToLNu_HT-*_TuneCP5_13TeV-madgraphMLM-pythia8'
        ],
        'other' : [
            'GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8',
            'ST_tW_Dilept_5f_DR_TuneCP5*',
            'tZq_ll_4f_ckm_NLO_TuneCP5*',
        ],
        'vv' : [
            'WZTo3LNu_Tune*',
            'WWTo2L2Nu_TuneCP5_DoubleScattering_*',
            'VHToNonbb_M125_Tune*',
            'ZZTo4L_Tune*',
            'WpWpJJ_EWK-QCD_Tune*',
            "SSWW_Tune*",
        ],
        # VV
        # To be done
    },
}


def get_dataset(data_type, groups, year):
    year = year if "2016" not in year else "2016"
    final_set = list()
    for group in groups:
        # first check if all exist
        if group not in dataset_dict[data_type]:
            raise Exception(f'{group} not found in dataset! Please add or amend. Avaliable are {dataset_dict[data_type].keys()}')
        if isinstance(dataset_dict[data_type][group], list):
            final_set += dataset_dict[data_type][group]
        else:
            final_set.append(dataset_dict[data_type][group][year])
    return final_set
    
def get_condition(data_type, year):
    if year not in condition_dict[data_type]:
        raise Exception(f'{year} not found in conditions!')
    return condition_dict[data_type][year]

def get_DAS(data_type, groups, year):
    nano_type = "NANOAOD" if data_type == 'data' else "NANOAODSIM"
    conditions = get_condition(data_type, year)
    final_das = list()

    for dataset in get_dataset(data_type, groups, year):
        for condition in conditions:
            final_das.append(f'/{dataset}/{condition}/{nano_type}')
    return final_das

# DAS functions

def get_nevents(dataset):
    query = f'dasgoclient -query="file dataset={dataset} | grep file.nevents"'
    file_sizes = subprocess.check_output(query, shell=True).decode().split()
    return sum(np.array(file_sizes, dtype=int))

def get_files(dataset):
    query = f'dasgoclient -query="file dataset={dataset}"'
    return subprocess.check_output(query, shell=True).decode().split()

def get_fullDAS(dataset):
    query = f"dasgoclient -query='dataset={dataset}'"
    return subprocess.check_output(query, shell=True).decode().split()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", help="output filename")
    parser.add_argument('-g', "--groups", required=True, help="Groups to analyze", type=lambda x : x.split(','))
    parser.add_argument("-y", "--years", required=True, type=lambda x : [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-t', "--type", choices=['data', 'mc'], required=True,
                        help="Type of dataset to use")
    parser.add_argument('-r', '--run', choices = ['output', 'print', 'events'], required=True,
                        help="Print details about the files, don't create lists")
    args = parser.parse_args()

    runfile_dir = analysis_area / "runfiles"
    runfile_dir.mkdir(exist_ok=True)


    for year in args.years:
        if args.run == 'output':
            if args.filename is None:
                raise Exception("No output filename given!!")
            output = list()

        datasets = get_DAS(args.type, args.groups, year)
        for dataset in datasets:
            for das in get_fullDAS(dataset):
                if args.run == 'print':
                    print(das)
                elif args.run == 'events':
                    print(das, get_nevents(das))
                elif args.run == 'output':
                    output += get_files(das)

        if args.run == 'output':
            with open(runfile_dir / f'{args.filename}_{args.type}_{year}.dat', "w") as f:
                for data_file in output:
                    f.write(f'{xrd_tag}{data_file}\n')
