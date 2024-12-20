#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user
from copy import deepcopy

trees = ["Signal_Dilepton", 'Signal_Multi',
         "Nonprompt_Dilepton", "Nonprompt_Multi", "OS_Charge_MisId", ]

info = NtupleInfo(
    filename = user.hdfs_area / 'workspace/signal_region_nlo/{year}/{workdir}/',
    trees = trees,
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
        # "tttj_nlo": "firebrick",
        # "tttw_nlo": "indianred",
        "tttj_lo": "firebrick",
        "tttw_lo": "indianred",
        # 'ttt_lo': 'crimson',
        'nonprompt': 'gray',
        # 'ttt_nlo': 'crimson',
        'nonprompt_mc': 'blue',
        "xg": "indigo",
        "ttw": "olivedrab",
        "tth": "goldenrod",
        "ttz": "steelblue",
        "ttXY": "teal",
        "rare": "deeppink",
        "4top": "tomato",
        'charge_flip': 'mediumseagreen',
    }
)
info.set_groups_trees(["Nonprompt_Dilepton", "Nonprompt_Multi"], ['nonprompt', 'nonprompt_mc'])
info.set_groups_trees(["OS_Charge_MisId"], ['charge_flip'])

dilep_ntuple = deepcopy(info)
dilep_ntuple.cut.append(lambda vg: vg["TightLepton"].num() == 2)

multi_ntuple = deepcopy(info)
multi_ntuple.cut.append(lambda vg: vg["TightLepton"].num() > 2)
