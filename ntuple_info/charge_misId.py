#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

color_by_group = {
        "data": "black",
        "DY_ht": "goldenrod",
        "DY": "goldenrod",
        "ttbar_lep": 'royalblue',
        "VV": 'mediumorchid',
        'charge_flip': 'seagreen',
    }

measurement = NtupleInfo(
    filename = user.hdfs_workspace / 'charge_misId/{year}/{workdir}/',
    trees = ['OS_MR'],
    cut=[
        lambda vg : vg["Met"] > 25,
    ],
    branches = [
        lambda vg : vg.mergeParticles("Electron", "TightElectron"),
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ],
    color_by_group = color_by_group,
)

closure_os = NtupleInfo(
    filename = user.hdfs_workspace / 'charge_misId/{year}/{workdir}/',
    trees = ['OS_CR'],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ],
    color_by_group = color_by_group,
)

closure_ss = NtupleInfo(
    filename = user.hdfs_workspace / 'charge_misId/{year}/{workdir}/',
    trees = ['SS', 'OS_CR'],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ],
    color_by_group = color_by_group,
)
closure_ss.add_change('OS_CR', {'charge_flip': 'data'})
closure_ss.set_ignore_trees("ttbar_lep", ["OS_CR"])
closure_ss.set_ignore_trees("VV", ["OS_CR"])
closure_ss.set_ignore_trees("DY", ["OS_CR"])
