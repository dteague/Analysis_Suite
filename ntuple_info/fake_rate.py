#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

color_by_group = {
    "data": "black",
    "qcd": "grey",
    "qcd_em": "lightgrey",
    "qcd_mu": "darkgrey",
    "ewk": "orange",
    "wjet_ht": "olive",
    "wjets": "olive",
    "ttbar_lep": "royalblue",
    "ttbar": "royalblue",
    'nonprompt': 'grey',
    'nonprompt_mc': 'grey',
    "DY_ht": "goldenrod",
    "DY": "goldenrod",
    "DY_J": "goldenrod",
    "vv_inc": 'mediumorchid',
    "wz": 'mediumorchid',

    "ttbar_2l2n": 'blue',
    "ttbar_semilep": 'mediumblue',
    "ttbar_hadronic": 'cornflowerblue',
}

ptratio = 0.0
btag = 0.25

pt_correction = NtupleInfo(
    filename = user.hdfs_workspace / 'singleLep/{year}/{workdir}/',
    trees = ["Measurement"],
    branches = [
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("FakeLepton", "FakeMuon", "FakeElectron"),
    ],
    color_by_group = color_by_group,
)


sideband = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate_mcm/{year}/{workdir}/',
    trees = ['SideBand'],
    part_cut=[['AllLepton', 'pt', lambda var: var > 20]],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("FakeLepton", "FakeMuon", "FakeElectron"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("AllLepton", "FakeMuon", "FakeElectron", "TightMuon", "TightElectron",),
    ],
    color_by_group = color_by_group,
)

measurement = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate_mcm/{year}/{workdir}/',
    trees = ['Measurement'],
    part_cut=[['AllLepton', 'pt', lambda var: var > 20],
              ['FakeElectron', 'jet_btag', lambda var: var < btag],
              ['FakeElectron', 'ptRatio', lambda var: var > ptratio]
              ],
    cut=[lambda vg: vg['AllLepton'].num()==1],
    branches = [
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("FakeLepton", "FakeMuon", "FakeElectron"),
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "FakeMuon", "FakeElectron", "TightMuon", "TightElectron"),
    ],
    color_by_group = color_by_group,
)

closure_tf = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate_closure/{year}/{workdir}/',
    trees = ['Closure_TF'],
    part_cut=[['AllLepton', 'pt', lambda var: var > 20]],
    cut = [
        lambda vg : vg.Jets.num() >= 1,
        lambda vg : vg["HT"] > 125,
        lambda vg : vg['AllLepton'].num()==2,
    ],
    branches = [
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("FakeLepton", "FakeMuon", "FakeElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "TightMuon", "TightElectron", "FakeMuon", "FakeElectron",),
    ],
    color_by_group = color_by_group,
)


closure_tt = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate_closure/{year}/{workdir}/',
    trees = ['Closure_FF', 'Closure_TF', 'Closure_TT'],
    part_cut=[['AllLepton', 'pt', lambda var: var > 20],
              ['FakeElectron', 'jet_btag', lambda var: var < btag],
              ['FakeElectron', 'ptRatio', lambda var: var > ptratio]
              ],
    cut = [
        lambda vg : vg.Jets.num() >= 1,
        lambda vg : vg["HT"] > 125,
        lambda vg : vg['AllLepton'].num()==2,
    ],
    branches = [
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("FakeLepton", "FakeMuon", "FakeElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "TightMuon", "TightElectron", "FakeMuon", "FakeElectron",),
    ],
    color_by_group = color_by_group,
)
closure_tt.set_groups_trees(['Closure_FF', "Closure_TF"], ['nonprompt', 'nonprompt_mc'])


dy_tight = NtupleInfo(
    filename = user.hdfs_workspace / 'fr_closure/{year}/{workdir}/',
    trees = ["DY_Fake", "DY_Tight"],
    branches = [
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("AllLepton", "FakeMuon", "TightMuon", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
    ],
    color_by_group = color_by_group,
)
dy_tight.set_groups_trees("DY_Fake", ['nonprompt', 'nonprompt_mc'])
