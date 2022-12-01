#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

measurement = NtupleInfo(
    filename = user.hdfs_workspace / 'charge_misId/{year}/{workdir}/',
    trees = ['OS_MR'],
    branches = [
        lambda vg : vg.mergeParticles("Electron", "TightElectron"),
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ]
)

closure_os = NtupleInfo(
    filename = user.hdfs_workspace / 'charge_misId/{year}/{workdir}/',
    trees = ['OS_CR'],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ]
)

closure_ss = NtupleInfo(
    filename = user.hdfs_workspace / 'charge_misId/{year}/{workdir}/',
    trees = ['SS', 'OS_CR'],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ]
)
closure_ss.add_change('OS_CR', {'data': 'charge_flip'})
closure_ss.set_exclusive('OS_CR', "charge_flip")
