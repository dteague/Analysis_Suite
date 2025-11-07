#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user
from copy import deepcopy
import awkward as ak

trees = ["Signal_Dilepton", 'Signal_Multi',
         "Nonprompt_Dilepton", "Nonprompt_Multi", "OS_Charge_MisId", ]

def any_zmass(vg):
    masses = ak.concatenate([vg['dilepton_masses'], vg['os_masses']], axis=1)
    # masses = vg['os_masses']
    return ak.any(abs(masses-91.188)<15, axis=1)

# trees = ['Signal_Multi',"Nonprompt_Multi"]

info = NtupleInfo(
    filename = user.hdfs_area / 'workspace/signal_region_nlo/{year}/{workdir}/',
    trees = trees,
    chan="Dilepton",
    cut=[
        lambda vg : vg['passZVeto'] == 1,
        lambda vg : vg.Jets.num() >= 2,
        lambda vg : vg["NBjets_medium"] >= 1,
        lambda vg : vg["Met"] > 50,
        lambda vg : vg["hasVetoJet"] == 0,
        lambda vg : vg['passZVeto'] == 1,
        # lambda vg : ~any_zmass(vg),
        # lambda vg : ~ss_zmass(vg),
        lambda vg : vg["HT"] > 250,
        #
    ],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ],
    color_by_group = {
        'data': 'k',
        'nonprompt': 'grey',
        'nonprompt_mc': 'grey',
        "xg": "purple",
        "ttw": "tan",
        "tth": "darkgrey",
        "ttz": "orange",
        "4top": "skyblue",
        'charge_flip': 'red',
        "rare": "brown",
        'ttt_nlo': 'k',
    }
)
info.set_groups_trees(["Nonprompt_Dilepton", "Nonprompt_Multi"], ['nonprompt', 'nonprompt_mc'])
info.set_groups_trees(["OS_Charge_MisId"], ['charge_flip'])


dilep_ntuple = NtupleInfo(
    filename = user.hdfs_area / 'workspace/signal_region_nlo/{year}/{workdir}/',
    trees = ["Signal_Dilepton", "Nonprompt_Dilepton", "OS_Charge_MisId", ],
    chan="Dilepton",
    cut=[
        lambda vg: vg["TightLepton"].num() == 2,
        lambda vg : vg['passZVeto'] == 1,
        lambda vg : vg.Jets.num() >= 2,
        lambda vg : vg["NBjets_medium"] >= 1,
        lambda vg : vg["hasVetoJet"] == 0,
        lambda vg : vg["Met"] > 50,
        lambda vg : vg["HT"] > 250,
        #
        # lambda vg: vg["TightLepton"].num() == 2,
        # lambda vg : vg['passZVeto'] == 1,
        # # lambda vg : vg.Jets.num() >= 2,
        # lambda vg : vg["NBjets_medium"] >= 1,
        # lambda vg : vg["Met"] > 50,
        # # lambda vg : vg.BJets.num() > 0,
        # lambda vg : vg["HT"] > 250,
        #
    ],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ],
    color_by_group = {
        'data': 'k',
        'nonprompt': 'grey',
        'nonprompt_mc': 'grey',
        "xg": "purple",
        "ttw": "tan",
        "tth": "darkgrey",
        "ttz": "orange",
        "4top": "skyblue",
        'charge_flip': 'red',
        "rare": "brown",
        'ttt_nlo': 'k',
    }
)
dilep_ntuple.set_groups_trees(["Nonprompt_Dilepton"], ['nonprompt', 'nonprompt_mc'])
dilep_ntuple.set_groups_trees(["OS_Charge_MisId"], ['charge_flip'])

multi_ntuple = NtupleInfo(
    filename = user.hdfs_area / 'workspace/signal_region_nlo/{year}/{workdir}/',
    trees = ['Signal_Multi', "Nonprompt_Multi"],
    chan="Dilepton",
    cut=[
        # lambda vg: vg["TightLepton"].num() > 2,
        lambda vg : vg['passZVeto'] == 1,
        lambda vg : vg.Jets.num() >= 2,
        lambda vg : vg["NBjets_medium"] >= 1,
        lambda vg : vg["hasVetoJet"] == 0,
        lambda vg : vg["Met"] > 50,
        lambda vg : vg["HT"] > 250,
        # lambda vg : vg["Met"] > 25,
        # lambda vg : vg["HT"] > 150,

    ],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ],
    color_by_group = {
        'data': 'k',
        'nonprompt': 'grey',
        'nonprompt_mc': 'grey',
        "xg": "purple",
        "ttw": "tan",
        "tth": "darkgrey",
        "ttz": "orange",
        "4top": "skyblue",
        "rare": "brown",
        'ttt_nlo': 'k',
    }
)
multi_ntuple.set_groups_trees(["Nonprompt_Multi"], ['nonprompt', 'nonprompt_mc'])
# multi_ntuple.set_groups_trees(["Nonprompt_Multi"], ['nonprompt'])
