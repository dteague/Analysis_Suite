#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

trees = ["Signal_Dilepton", 'Signal_Multi',
         "Nonprompt_Dilepton", "Nonprompt_Multi", "OS_Charge_MisId", ]

info = NtupleInfo(
    filename = user.hdfs_area / 'workspace/signal_region/{year}/{workdir}/',
    trees = trees,
    cut=[
        lambda vg : vg['passZVeto'] == 1,
        # lambda vg : vg["N_bloose"] >= 2,
        lambda vg : vg.Jets.num() >= 2,
        lambda vg : vg["NBjets_medium"] >= 1,
        lambda vg : vg["Met"] > 50,
        # lambda vg : vg.BJets.num() > 0,
        # lambda vg : vg["HT"] > 300,

    ],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ],
    color_by_group = {
        "ttt": "crimson",
        'nonprompt': 'gray',
        'nonprompt_mc': 'blue',
        "xg": "indigo",
        "ttw": "olivedrab",
        "tth": "goldenrod",
        "ttz": "steelblue",
        "ttXY": "teal",
        "rare": "deeppink",
        "tttt": "tomato",
        # 'charge_flip': 'mediumseagreen',
    }
)

info.add_change('Nonprompt_Dilepton', {'nonprompt': 'data'})
info.add_change('Nonprompt_Multi', {'nonprompt': 'data'})
info.add_change('OS_Charge_MisId', {'charge_flip': 'data'})

info.set_ignore_trees("nonprompt_mc", ["Signal_Dilepton", "Signal_Multi", "OS_Charge_MisId"])
