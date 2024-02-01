#!/usr/bin/env python3
from analysis_suite.commons.info import GroupInfo
import numpy as np
from mt2 import mt2

pad = -1

def mt2_l(vg):
    l1_x = vg.TightLepton['pt', 0]*np.cos(vg.TightLepton['phi', 0])
    l1_y = vg.TightLepton['pt', 0]*np.sin(vg.TightLepton['phi', 0])
    l2_x = vg.TightLepton['pt', 1]*np.cos(vg.TightLepton['phi', 1])
    l2_y = vg.TightLepton['pt', 1]*np.sin(vg.TightLepton['phi', 1])
    mx = vg["Met"]*np.cos(vg["Met_phi"])
    my = vg["Met"]*np.sin(vg["Met_phi"])
    return mt2(np.abs(vg.TightLepton['mass', 0]), l1_x, l1_y,
               np.abs(vg.TightLepton['mass', 1]), l2_x, l2_y,
               mx, my, 0., 0.)


# Variables used in Training
allvar = {
    "NJets" :           lambda vg : vg.Jets.num(),
    # "NResolvedTops":    lambda vg : vg.Tops.num(),
    "NlooseBJets":      lambda vg : vg['NBjets_loose'],
    "NmediumBJets":     lambda vg : vg['NBjets_medium'],
    "NtightBJets":      lambda vg : vg['NBjets_tight'],

    "NlooseMuons":      lambda vg : vg['N_loose_mu'],
    "NlooseElectrons":  lambda vg : vg['N_loose_el'],
    "NMuons":           lambda vg : vg.TightMuon.num(),
    "NElectrons":       lambda vg : vg.TightElectron.num(),

    "HT":               lambda vg : vg['HT'],
    "HT_b":             lambda vg : vg['HT_b'],
    "Met":              lambda vg : vg['Met'],
    "centrality":       lambda vg : vg['Centrality'],
    # "ZMass":            lambda vg : vg['Zmass'],

    "j1Pt":             lambda vg : vg.Jets['pt', 0, pad],
    "j2Pt":             lambda vg : vg.Jets['pt', 1, pad],
    "j3Pt":             lambda vg : vg.Jets['pt', 2, pad],
    "j4Pt":             lambda vg : vg.Jets['pt', 3, pad],
    "j5Pt":             lambda vg : vg.Jets['pt', 4, pad],
    "j6Pt":             lambda vg : vg.Jets['pt', 5, pad],
    # "j7Pt":             lambda vg : vg.Jets['pt', 6, pad],
    # "j8Pt":             lambda vg : vg.Jets['pt', 7, pad],

    "b1Pt":             lambda vg : vg.BJets['pt', 0, pad],
    "b2Pt":             lambda vg : vg.BJets['pt', 1, pad],
    "b3Pt":             lambda vg : vg.BJets['pt', 2, pad],
    "b4Pt":             lambda vg : vg.BJets['pt', 3, pad],

    "l1Pt":              lambda vg : vg.TightLepton["pt", 0],
    "l2Pt":              lambda vg : vg.TightLepton["pt", 1],

    "lep_mass" :         lambda vg : vg.mass("TightLepton", 0, "TightLepton", 1),
    # "lep_dphi" :         lambda vg : vg.dphi("TightLepton", 0, "TightLepton", 1),
    # "lep_deta" :         lambda vg : vg.TightLepton['eta', 0] - vg.TightLepton['eta', 1],
    # "lep_dr" :           lambda vg : vg.dr("TightLepton", 0, "TightLepton", 1),
    "lep_mt":            lambda vg : vg.dipart_mt("TightLepton", 0, "TightLepton", 1),
    "lep_cosTheta" :          lambda vg : vg.cosDtheta("TightLepton", 0, "TightLepton", 1),

    "jet_dr" :           lambda vg : vg.dr("Jets", 0, "Jets", 1),
    "jet_mass" :         lambda vg : vg.mass("Jets", 0, "Jets", 1),
    # "jet_dphi" :         lambda vg : vg.dphi("Jets", 0, "Jets", 1),
    # "jet_deta" :         lambda vg : vg.Jets['eta', 0] - vg.Jets['eta', 1],
    "jet_mt":            lambda vg : vg.dipart_mt("Jets", 0, "Jets", 1),
    "jet_cosTheta" :     lambda vg : vg.cosDtheta("Jets", 0, "Jets", 1),

    # "JetLep1_Cos" :     lambda vg : vg.cosDtheta("TightLepton", 0, "Jets", 0),
    # "JetLep2_Cos" :     lambda vg : vg.cosDtheta("TightLepton", 1, "Jets", 0),

    "mT_1":             lambda vg : vg.TightLepton['mt', 0],
    "mT_2":             lambda vg : vg.TightLepton['mt', 1],
    "mT2_l" :           lambda vg : mt2_l(vg),
    # "cosdphi_1":        lambda vg : np.cos(vg.TightLepton['phi', 0] - vg["Met_phi"]),
    # "cosdphi_2":        lambda vg : np.cos(vg.TightLepton['phi', 1] - vg["Met_phi"]),
}


# Vars to actually use in training
usevars = list(allvar.keys())

# Samples and the groups they are a part of
groups = {
    "Signal": ["ttt"],
    "Background": [
        "4top",
        "ttw", "ttz", "tth",
        "ttXY",
        "rare", "xg",
    ],
    "NotTrained": [
        "nonprompt", "charge_flip",
        'data'],
    "OnlyTrain": ['nonprompt_mc',]
}


remove = ['nonprompt_mc']

def wz_scale(vg, year):
    import numpy as np
    import pickle

    scale_file = workdir/f"wz_scale_factor.pickle"
    with open(scale_file, "rb") as f:
        scales = pickle.load(f)
    jet_bin = np.digitize(vg["NJets"], np.arange(1, 6)) - 1
    vg.scale = scales[year].vals

# scale = {}
scale = {'wz': lambda vg: wz_scale(vg)}

# Variables needed in code for things to work
assert "allvar" in locals()
assert 'usevars' in locals()
assert 'groups' in locals()
