#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

trees = ['Signal_Multi', "Nonprompt_Multi"]

info = NtupleInfo(
    filename = user.hdfs_area / 'workspace/signal_region_nlo/{year}/{workdir}/',
    trees = trees,
    cut=[
        # lambda vg : vg['passZVeto'] == 0,
        lambda vg : vg["TightLepton"].num() == 3,#
        lambda vg : vg["hasVetoJet"] == 0,
        lambda vg : vg.Jets.num() >= 2,
        # lambda vg : vg["N_loose_mu"]+vg["N_loose_el"] == 0,
        lambda vg : vg['NBjets_medium'] == 0,
    ],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron")
    ],
    color_by_group = {
        'data': 'k',
        'nonprompt': 'grey',
        "xg": "purple",
        'ttX': 'orange',
        "4top": "skyblue",
        "rare_nowz": "brown",
        'wz':   'red',
        'ttt_nlo': 'k',
    }
)
info.set_groups_trees(["Nonprompt_Multi"], ['nonprompt', 'nonprompt_mc'])
