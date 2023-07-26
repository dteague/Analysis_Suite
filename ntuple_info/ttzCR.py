#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

trees = ['Signal_Multi', "Nonprompt_Multi", "OS_Charge_MisId", ]

info = NtupleInfo(
    filename = user.hdfs_area / 'workspace/signal_region/{year}/{workdir}/',
    trees = trees,
    cut=[
        lambda vg : vg['passZVeto'] == 0,
        lambda vg : vg["TightLepton"].num() == 3,
        lambda vg : vg["N_loose_mu"]+vg["N_loose_el"] == 0,
        lambda vg : vg['NBjets_medium'] >= 1,
        lambda vg : vg.Jets.num() >= 2,
        # lambda vg : vg['HLT_trilepton'],
        lambda vg : vg["HLT_dilepton"] + vg['HLT_trilepton'],

        # lambda vg : vg["HLT_channel"],
        lambda vg : vg["Met"] > 50,
        # lambda vg : vg["HT"] > 300,
        # lambda vg : vg.TightLepton['pt', 1] > 20,
        # lambda vg : vg.TightLepton['pt', 2, -1] > 20,
    ],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron")
    ],
    color_by_group = {
        "ttt": "crimson",
        'nonprompt': 'gray',
        "xg": "indigo",
        "ttw": "olivedrab",
        "tth": "goldenrod",
        "ttz": "steelblue",
        "ttXY": "teal",
        # "rare_nowz": "deeppink",
        # "wz": 'slateblue',
        'rare': 'blue',
        "tttt": "tomato",
        # 'charge_flip': 'mediumseagreen',
        'data': 'black',
    }
)

info.add_change('Nonprompt_Multi', {'nonprompt': 'data'})
info.add_change('OS_Charge_MisId', {'charge_flip': 'data'})
