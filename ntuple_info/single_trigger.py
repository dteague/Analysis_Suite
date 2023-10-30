#!/usr/bin/env python3

from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

measurement = NtupleInfo(
    filename = user.hdfs_workspace / 'lepton_trigger/{year}/{workdir}/',
    trees = ['Measurement'],
    branches = [
        lambda vg : vg.mergeParticles("Electron", "FakeElectron", "TightElectron"),
        lambda vg : vg.mergeParticles("Muon", "FakeMuon", "TightMuon"),
        lambda vg : vg.mergeParticles("AllLepton", "FakeMuon", "FakeElectron", "TightMuon", "TightElectron"),
    ]
)
