#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

color_by_group = {
    "data": "black",
    "qcd": "grey",
    "qcd_em": "grey",
    "qcd_mu": "grey",
    "ewk": "orange",
    "wjet_ht": "olive",
    "wjets": "olive",
    "ttbar_lep": "royalblue",
    'nonprompt': 'grey',
    "DY_ht": "goldenrod",
    "DY": "goldenrod",
    "VV": 'mediumorchid',

        "ttbar_2l2n": 'blue',
    "ttbar_semilep": 'mediumblue',
    "ttbar_hadronic": 'cornflowerblue',
}

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
    filename = user.hdfs_workspace / 'fake_rate/{year}/{workdir}/',
    trees = ['SideBand'],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("FakeLepton", "FakeMuon", "FakeElectron"),
        # lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("AllLepton", "FakeMuon", "FakeElectron", "TightMuon", "TightElectron", # "LooseMuon", "LooseElectron"
                                 ),
    ],
    color_by_group = color_by_group,
)

measurement = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate/{year}/{workdir}/',
    trees = ['Measurement'],
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
    branches = [
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("FakeLepton", "FakeMuon", "FakeElectron"),
        # lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "TightMuon", "TightElectron", "FakeMuon", "FakeElectron", # "LooseMuon", "LooseElectron"
                                 ),
    ],
    color_by_group = color_by_group,
)

closure_tt = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate_closure/{year}/{workdir}/',
    trees = ['Closure_FF', 'Closure_TF', 'Closure_TT'],
    branches = [
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("FakeLepton", "FakeMuon", "FakeElectron"),
        # lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "TightMuon", "TightElectron", "FakeMuon", "FakeElectron", # "LooseMuon", "LooseElectron"
                                 ),
    ],
    color_by_group = color_by_group,
)
closure_tt.add_change('Closure_FF', {'nonprompt': 'data'})
closure_tt.add_change('Closure_TF', {'nonprompt': 'data'})
closure_tt.set_ignore_trees("ttbar_lep", ["Closure_FF", "Closure_TF"])
closure_tt.set_ignore_trees("wjet_ht", ["Closure_FF", "Closure_TF"])

dy_fake = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate_closure/{year}/{workdir}/',
    trees = ["DY_Fake",],
    branches = [
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        # lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "FakeMuon", "TightMuon", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
    ],
    color_by_group = color_by_group,
)



dy_tight = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate_closure/{year}/{workdir}/',
    trees = ["DY_Fake", "DY_Tight"],
    branches = [
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        # lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "FakeMuon", "TightMuon", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
    ],
    color_by_group = color_by_group,
)
dy_tight.add_change('DY_Fake', {'nonprompt': 'data'})
