#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

sideband = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate/{year}/{workdir}/',
    trees = ['SideBand'],
    branches = [
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("FakeLepton", "FakeMuon", "FakeElectron"),
        lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("AllLepton", "FakeMuon", "FakeElectron", "TightMuon", "TightElectron", "LooseMuon", "LooseElectron"),
    ]
)

measurement = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate/{year}/{workdir}/',
    trees = ['Measurement'],
    branches = [
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "FakeMuon", "FakeElectron", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
    ]
)

closure_tf = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate_closure/{year}/{workdir}/',
    trees = ['Closure_TF'],
    branches = [
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("FakeLepton", "FakeMuon", "FakeElectron"),
        lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "TightMuon", "TightElectron", "FakeMuon", "FakeElectron", "LooseMuon", "LooseElectron"),
    ]
)

closure_tt = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate_closure/{year}/{workdir}/',
    trees = ['Closure_FF', 'Closure_TF', 'Closure_TT'],
    branches = [
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("TightLepton", "TightMuon", "TightElectron"),
        lambda vg : vg.mergeParticles("FakeLepton", "FakeMuon", "FakeElectron"),
        lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "TightMuon", "TightElectron", "FakeMuon", "FakeElectron", "LooseMuon", "LooseElectron"),
    ]
)
closure_tt.add_change('Closure_FF', {'nonprompt': 'data'})
closure_tt.set_exclusive('Closure_FF', 'data')
closure_tt.add_change('Closure_TF', {'nonprompt': 'data'})
closure_tt.set_exclusive('Closure_TF', 'data')



dy_fake = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate_closure/{year}/{workdir}/',
    trees = ["DY_Fake",],
    branches = [
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "FakeMuon", "TightMuon", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
    ]
)



dy_tight = NtupleInfo(
    filename = user.hdfs_workspace / 'fake_rate_closure/{year}/{workdir}/',
    trees = ["DY_Fake", "DY_Tight"],
    branches = [
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("LooseLepton", "LooseMuon", "LooseElectron"),
        lambda vg : vg.mergeParticles("AllLepton", "FakeMuon", "TightMuon", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
    ]
)
dy_tight.add_change('DY_Fake', {'nonprompt': 'data'})
dy_tight.set_exclusive('DY_Fake', "data")
