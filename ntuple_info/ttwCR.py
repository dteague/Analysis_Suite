#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

trees = ["Signal_Dilepton", 'Signal_Multi',
         "Nonprompt_Dilepton", "Nonprompt_Multi", "OS_Charge_MisId", ]

info = NtupleInfo(
    filename = user.hdfs_area / 'workspace/signal_region/{year}/{workdir}/',
    trees = trees,
    cut=[
        lambda vg : vg["TightLepton"].num() == 2,
        lambda vg : vg['passZVeto'] == 1,
        lambda vg : vg["NBjets_medium"] >= 1,
        lambda vg : vg.Jets.num() >= 2,
        # lambda vg : vg.BJets.num() > 0,
        # lambda vg : vg["HT"] > 300,
        # lambda vg : vg["Met"] > 75
    ],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ],
    color_by_group = {
        "ttt": "crimson",
        'nonprompt': 'gray',
        "xg": "indigo",
        "ttw": "olivedrab",
        "tth": "goldenrod",
        "ttz": "steelblue",
        "ttXY": "teal",
        "rare": "deeppink",
        "tttt": "tomato",
        'charge_flip': 'mediumseagreen',
        'data': 'black',
    }
)
info.set_groups_trees(["Nonprompt_Dilepton", "Nonprompt_Multi"], ['nonprompt', 'nonprompt_mc'])
info.set_groups_trees(["OS_Charge_MisId"], ['charge_flip'])
