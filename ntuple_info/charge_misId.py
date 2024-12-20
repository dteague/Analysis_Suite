#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

color_by_group = {
        "data": "black",
        "DY_ht": "goldenrod",
        "DY": "goldenrod",
        "ttbar_lep": 'royalblue',
        "vv_inc": 'mediumorchid',
        'charge_flip': 'seagreen',
    }

measurement = NtupleInfo(
    filename = user.hdfs_workspace / 'charge_misId/{year}/{workdir}/',
    trees = ['OS_MR'],
    part_cut=[['TightLepton', 'pt', lambda var: var > 20]],
    cut=[
        lambda vg : vg["Met"] > 25,
        lambda vg : vg['TightLepton'].num() >= 2,
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
    part_cut=[['TightLepton', 'pt', lambda var: var > 20]],
    cut=[
        lambda vg : vg['TightLepton'].num() >= 2,
    ],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ],
    color_by_group = color_by_group,
)

closure_ss = NtupleInfo(
    filename = user.hdfs_workspace / 'charge_misId/{year}/{workdir}/',
    trees = ['SS', 'OS_CR'],
    part_cut=[['TightLepton', 'pt', lambda var: var > 20]],
    cut=[
        lambda vg : vg['TightLepton'].num() >= 2,
    ],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
    ],
    color_by_group = color_by_group,
)
closure_ss.set_groups_trees("OS_CR", ['charge_flip'])
