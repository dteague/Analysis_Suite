#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user
from copy import deepcopy

trees = ["Signal_Dilepton", 'Signal_Multi',
         "Nonprompt_Dilepton", "Nonprompt_Multi", "OS_Charge_MisId", ]

info = NtupleInfo(
    filename = user.hdfs_area / 'workspace/signal_region_nlo/{year}/{workdir}/',
    trees = trees,
    chan="Dilepton",
    cut=[
        lambda vg : vg['passZVeto'] == 1,
        lambda vg : vg.Jets.num() >= 2,
        lambda vg : vg["NBjets_medium"] >= 1,
        lambda vg : vg["Met"] > 50,
        # lambda vg : vg.BJets.num() > 0,
        lambda vg : vg["HT"] > 250,
        #
    ],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ],
    color_by_group = {
        'ttt_nlo': 'crimson',
        "4top": "tomato",
        'data_driven': 'gray',
        'ttX': 'blueviolet',
        'other': 'olivedrab'
    }
)
