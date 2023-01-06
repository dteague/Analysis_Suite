#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

trees = ["Nonprompt_FakeRate", "OS_Charge_MisId", "Signal_Dilepton", 'Signal_Multi']

info = NtupleInfo(
    filename = user.hdfs_area / 'workspace/signal_region/{year}/{workdir}/',
    trees = trees,
    region = 'Signal',
    cut=lambda vg : vg['passZVeto'] == 1,
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ]
)

info.add_change('Nonprompt_FakeRate', {'nonprompt': 'data'})
info.add_change('OS_Charge_MisId', {'charge_flip': 'data'})
