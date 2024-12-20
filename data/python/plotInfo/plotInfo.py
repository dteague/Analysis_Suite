#!/usr/bin/env python3
import numpy as np
import boost_histogram.axis as axis
import awkward as ak
from copy import deepcopy

from analysis_suite.plotting.plotter import GraphInfo

main_plots = [
    GraphInfo('njets', '$N_{j}$', axis.Regular(12, 0, 12), 'NJets'),
    GraphInfo('nloosebjets', '$N_{looseb}$', axis.Regular(8, 0, 8), 'NlooseBJets'),
    GraphInfo('nmediumbjets', '$N_{mediumb}$', axis.Regular(8, 0, 8), 'NmediumBJets'),
    GraphInfo('ntightbjets', '$N_{tightb}$', axis.Regular(7, 0, 7), 'NtightBJets'),
    GraphInfo('ht', '$H_{T}$ (GeV)', axis.Regular(20, 0, 1500), 'HT'),
    GraphInfo('ht_b', '$H_{T}(b)$ (GeV)', axis.Regular(20, 0, 1200), 'HT_b'),
    GraphInfo('met', '$p_{T}^{miss}$ (GeV)', axis.Regular(20, 0, 500), 'Met'),
    GraphInfo('centrality', '$H_{T}/E$', axis.Regular(20, 0, 1), 'centrality'),
    GraphInfo('l1Pt', '$p_{T}(l_{1})$ (GeV)', axis.Regular(20, 0, 500), 'l1Pt'),
    GraphInfo('l2Pt', '$p_{T}(l_{2})$ (GeV)', axis.Regular(20, 0, 200), 'l2Pt'),

]

plots = [
    GraphInfo('njets', '$N_{j}$', axis.Regular(12, 0, 12), 'NJets'),
    # GraphInfo('nboostedtops', '$N_{t}$', axis.Regular(4, 0, 4), 'NResolvedTops'),
    GraphInfo('nloosebjets', '$N_{looseb}$', axis.Regular(8, 0, 8), 'NlooseBJets'),
    GraphInfo('nmediumbjets', '$N_{mediumb}$', axis.Regular(8, 0, 8), 'NmediumBJets'),
    GraphInfo('ntightbjets', '$N_{tightb}$', axis.Regular(7, 0, 7), 'NtightBJets'),

    GraphInfo('nlooseMus', '$N_{loose\mu}$', axis.Regular(4, 0, 4), 'NlooseMuons'),
    GraphInfo('nlooseEls', '$N_{loosee}$', axis.Regular(4, 0, 4), 'NlooseElectrons'),
    GraphInfo('nmuons', '$N_{\mu}$', axis.Regular(4, 0, 4), 'NMuons'),
    GraphInfo('nelectrons', '$N_{e}$', axis.Regular(4, 0, 4), 'NElectrons'),

    GraphInfo('ht', '$H_{T}$ (GeV)', axis.Regular(20, 0, 1500), 'HT'),
    GraphInfo('ht_b', '$H_{T}(b)$ (GeV)', axis.Regular(20, 0, 1200), 'HT_b'),
    GraphInfo('met', '$p_{T}^{miss}$ (GeV)', axis.Regular(20, 0, 500), 'Met'),
    GraphInfo('centrality', '$H_{T}/E$', axis.Regular(20, 0, 1), 'centrality'),

    GraphInfo('l1Pt', '$p_{T}(l_{1})$ (GeV)', axis.Regular(20, 0.1, 500), 'l1Pt'),
    GraphInfo('l2Pt', '$p_{T}(l_{2})$ (GeV)', axis.Regular(20, 0.1, 200), 'l2Pt'),

    GraphInfo('lep_mass', '$M_{\ell\ell}$ (GeV)', axis.Regular(20, 0, 500), 'lep_mass'),
    # GraphInfo('lep_dphi', '$\Delta\phi_{\ell\ell}$', axis.Regular(25, -np.pi, np.pi), 'lep_dphi'),
    # GraphInfo('lep_deta', '$\Delta\eta_{\ell\ell}$', axis.Regular(25, -4, 4), 'lep_deta'),
    # GraphInfo('lep_dr', '$\Delta R_{\ell\ell}$', axis.Regular(25, 0, 6), 'lep_dr'),
    GraphInfo('lep_mt', '$M_{T}(\ell,\ell)$', axis.Regular(25, 0, 400), 'lep_mt'),
    GraphInfo('lep_cosTheta', r'$cos(\theta_{\ell\ell})$', axis.Regular(25, -1, 1), 'lep_cosTheta'),

    # GraphInfo('jet_mass', '$M_{jj}$ (GeV)', axis.Regular(20, 0, 800), 'jet_mass'),
    # # GraphInfo('jet_dphi', '$\Delta\phi_{jj}$', axis.Regular(25, -np.pi, np.pi), 'jet_dphi'),
    # # GraphInfo('jet_deta', '$\Delta\eta_{jj}$', axis.Regular(25, -4, 4), 'jet_deta'),
    # GraphInfo('jet_dr', '$\Delta R_{jj}$', axis.Regular(25, 0, 6), 'jet_dr'),
    # GraphInfo('jet_mt', '$M_{T}(jj)$', axis.Regular(25, 0, 750), 'jet_mt'),
    # GraphInfo('jet_cosTheta', r'$cos(\theta_{jj})$', axis.Regular(25, -1, 1), 'jet_cosTheta'),

    # GraphInfo('jet1cos', r'$\cos(\theta_{j\ell1})$', axis.Regular(20, -1, 1), 'JetLep1_Cos'),
    # GraphInfo('jet2cos', r'$\cos(\theta_{j\ell2})$', axis.Regular(20, -1, 1), 'JetLep2_Cos'),


    GraphInfo('mt1', '$M_{T}(\ell_{1})$', axis.Regular(20, 0, 500), 'mT_1'),
    GraphInfo('mt2', '$M_{T}(\ell_{2})$', axis.Regular(20, 0, 300), 'mT_2'),
    GraphInfo('mt2_l', '$M_{T}(\ell_{1})$', axis.Regular(20, 1, 150), 'mT2_l'),
    # GraphInfo('cosdphi_1', '$cos(\Delta\phi_{\ell1 Met})$', axis.Regular(25, -1, 1), 'cosdphi_1'),
    # GraphInfo('cosdphi_2', '$cos(\Delta\phi_{\ell2 Met})$', axis.Regular(25, -1, 1), 'cosdphi_2'),
]

jet_pt = [
    GraphInfo('j1Pt', '$p_{T}(j_{1})$ (GeV)', axis.Regular(20, 1, 650), 'j1Pt'),
    GraphInfo('j2Pt', '$p_{T}(j_{2})$ (GeV)', axis.Regular(20, 1, 500), 'j2Pt'),
    GraphInfo('j3Pt', '$p_{T}(j_{3})$ (GeV)', axis.Regular(20, 1, 300), 'j3Pt'),
    GraphInfo('j4Pt', '$p_{T}(j_{4})$ (GeV)', axis.Regular(20, 1, 250), 'j4Pt'),
    GraphInfo('j5Pt', '$p_{T}(j_{5})$ (GeV)', axis.Regular(20, 1, 150), 'j5Pt'),
    GraphInfo('j6Pt', '$p_{T}(j_{6})$ (GeV)', axis.Regular(20, 1, 150), 'j6Pt'),
    # GraphInfo('j7Pt', '$p_{T}(j_{7})$ (GeV)', axis.Regular(20, 1, 150), 'j7Pt'),
    # GraphInfo('j8Pt', '$p_{T}(j_{8})$ (GeV)', axis.Regular(20, 1, 100), 'j8Pt'),

    GraphInfo('b1Pt', '$p_{T}(b_{1})$ (GeV)', axis.Regular(20, 0.1, 500), 'b1Pt'),
    GraphInfo('b2Pt', '$p_{T}(b_{2})$ (GeV)', axis.Regular(20, 0.1, 300), 'b2Pt'),
    GraphInfo('b3Pt', '$p_{T}(b_{3})$ (GeV)', axis.Regular(20, 0.1, 200), 'b3Pt'),
    GraphInfo('b4Pt', '$p_{T}(b_{4})$ (GeV)', axis.Regular(20, 0.1, 100), 'b4Pt'),
]

signal =  [
    GraphInfo('signal', '$Disc_{Signal}$', axis.Regular(15, 0, 1.0), 'Signal'),
    # GraphInfo('signal_fom', '$Disc_{Signal}$', axis.Regular(25, 0, 0.8), 'Signal', output="fom"),
]

plots_ttz =  [
    # GraphInfo('signal', '$Disc_{Signal}$', axis.Regular(15, 0, 1.0), 'Signal'),
    # GraphInfo('mt1', '$M_{T}(\ell_{1})$', axis.Regular(20, 0, 200), 'mT_1'),
    # GraphInfo('mt2', '$M_{T}(\ell_{2})$', axis.Regular(20, 0, 200), 'mT_2'),
    GraphInfo('njets', '$N_{j}$', axis.Regular(12, 0, 12), 'NJets'),
    GraphInfo('nmuons', '$N_{\mu}$', axis.Regular(4, 0, 4), 'NMuons'),
    GraphInfo('nelectrons', '$N_{e}$', axis.Regular(4, 0, 4), 'NElectrons'),
    # GraphInfo('nlooseMus', '$N_{loose\mu}$', axis.Regular(4, 0, 4), 'NlooseMuons'),
    # GraphInfo('nlooseEls', '$N_{loosee}$', axis.Regular(4, 0, 4), 'NlooseElectrons'),
    GraphInfo('nloosebjets', '$N_{looseb}$', axis.Regular(8, 0, 8), 'NlooseBJets'),
    GraphInfo('nmediumbjets', '$N_{mediumb}$', axis.Regular(8, 0, 8), 'NmediumBJets'),
    GraphInfo('ntightbjets', '$N_{tightb}$', axis.Regular(7, 0, 7), 'NtightBJets'),
    GraphInfo('met', '$p_{T}^{miss}$ (GeV)', axis.Regular(16, 0, 400), 'Met'),
    GraphInfo('ht', '$H_{T}$ (GeV)', axis.Regular(16, 0, 1600), 'HT'),
    GraphInfo('ht_b', '$H_{T}(b)$ (GeV)', axis.Regular(15, 0, 600), 'HT_b'),
    GraphInfo('centrality', '$H_{T}/E$', axis.Regular(20, 0, 1), 'centrality'),
    GraphInfo('j1Pt', '$p_{T}(j_{1})$ (GeV)', axis.Regular(15, 0.1, 650), 'j1Pt'),
    GraphInfo('j2Pt', '$p_{T}(j_{2})$ (GeV)', axis.Regular(15, 0.1, 500), 'j2Pt'),
    GraphInfo('j3Pt', '$p_{T}(j_{3})$ (GeV)', axis.Regular(15, 0.1, 300), 'j3Pt'),
    GraphInfo('j4Pt', '$p_{T}(j_{4})$ (GeV)', axis.Regular(10, 0.1, 250), 'j4Pt'),
    # GraphInfo('j5Pt', '$p_{T}(j_{5})$ (GeV)', axis.Regular(10, 0.1, 150), 'j5Pt'),
    # GraphInfo('j6Pt', '$p_{T}(j_{6})$ (GeV)', axis.Regular(10, 0.1, 150), 'j6Pt'),
    # GraphInfo('j7Pt', '$p_{T}(j_{7})$ (GeV)', axis.Regular(10, 0.1, 150), 'j7Pt'),
    # GraphInfo('j8Pt', '$p_{T}(j_{8})$ (GeV)', axis.Regular(10, 0.1, 100), 'j8Pt'),
    GraphInfo('b1Pt', '$p_{T}(b_{1})$ (GeV)', axis.Regular(15, 0.1, 500), 'b1Pt'),
    GraphInfo('b2Pt', '$p_{T}(b_{2})$ (GeV)', axis.Regular(15, 0.1, 300), 'b2Pt'),
    # GraphInfo('b3Pt', '$p_{T}(b_{3})$ (GeV)', axis.Regular(10, 0.1, 200), 'b3Pt'),
    # GraphInfo('b4Pt', '$p_{T}(b_{4})$ (GeV)', axis.Regular(10, 0.1, 100), 'b4Pt'),
    GraphInfo('l1Pt', '$p_{T}(l_{1})$ (GeV)', axis.Regular(15, 0.1, 200), 'l1Pt'),
    GraphInfo('l2Pt', '$p_{T}(l_{2})$ (GeV)', axis.Regular(15, 0, 150), 'l2Pt'),
    GraphInfo('lepMass', '$M_{\ell\ell}$ (GeV)', axis.Regular(15, 0, 750), 'lep_mass'),
    # GraphInfo('ZMass', '$M_{Z}$ (GeV)', axis.Regular(15, 75, 105), 'ZMass'),
    # GraphInfo('met1cos', r'$\cos(\theta_{met\ell1})$', axis.Regular(15, -1, 1), 'cosdphi_1'),
    # GraphInfo('met2cos', r'$\cos(\theta_{met\ell2})$', axis.Regular(15, -1, 1), 'cosdphi_2'),
    # GraphInfo('mT_lin_1', '$M_{Tlin}(\ell_{1})$', axis.Regular(20, 0, 200), 'mT_lin_1'),
    # GraphInfo('mT_lin_2', '$M_{Tlin}(\ell_{2})$', axis.Regular(20, 0, 200), 'mT_lin_2'),

    # GraphInfo('jetMass', '$M_{jj}$ (GeV)', axis.Regular(15, 0, 1500), 'jetMass'),
    # GraphInfo('jetDR', '$\Delta R_{jj}$', axis.Regular(15, 0, 6), 'jetDR'),
    # GraphInfo('lepDR', '$\Delta R_{\ell\ell}$', axis.Regular(15, 0, 6), 'lepDR'),
    # GraphInfo('lepcos', r'$\cos(\theta_{\ell\ell})$', axis.Regular(15, -1, 1), 'LepCos'),
    # GraphInfo('jet1cos', r'$\cos(\theta_{j\ell1})$', axis.Regular(15, -1, 1), 'JetLep1_Cos'),
    # GraphInfo('jet2cos', r'$\cos(\theta_{j\ell2})$', axis.Regular(15, -1, 1), 'JetLep2_Cos'),
]

combine = {
    # 'ttz_cr' : ("CR-ttz", GraphInfo('nbjets', '$N_{b}$', axis.Regular(7, 0, 7), 'NBJets')),
    # 'ttw_cr' : ("Signal", GraphInfo('nbjets', '$N_{b}$', axis.Regular(7, 0, 7), 'NBJets', cuts="Signal<0.8")),
    # 'signal' : ("Signal", GraphInfo('signal', '$Disc_{Signal}$', axis.Regular(20, 0, 1), 'Signal')),
    # 'signal' : ("Signal", GraphInfo('ht_b', '$H_{T}(b)$ (GeV)', axis.Regular(20, 0, 600), 'HT_b')),
    'signal' : ("Signal", GraphInfo('yield', '$H_{T}(b)$ (GeV)', axis.Variable([0, 1000, 2000]), 'HT_b')),
}

dilepton = [
    GraphInfo("HT", r"$H_{T}$", axis.Regular(25, 0, 1500), lambda vg : vg.get_hist('HT')),
    GraphInfo("Met", r"Met", axis.Regular(25, 0, 500), lambda vg : vg.get_hist('Met')),
    GraphInfo('nelectrons', '$N_{e}$', axis.Regular(3, 0, 3), lambda vg: (vg.TightElectron.num(), vg.scale)),
    GraphInfo('centrality', '$H_{T}/E$', axis.Regular(20, 0, 1), lambda vg: vg.get_hist('Centrality')),
    # GraphInfo('l1Pt', '$p_{T}(l_{1})$ (GeV)', axis.Regular(25, 0, 500), lambda vg: vg.TightLepton.get_hist('pt', 0)),
    # GraphInfo('l2Pt', '$p_{T}(l_{2})$ (GeV)', axis.Regular(25, 0, 200), lambda vg: vg.TightLepton.get_hist('pt', 1)),
    GraphInfo("Njet", r"$N_{j}$", axis.Regular(9, 0, 9), lambda vg : (vg.Jets.num(), vg.scale)),
    GraphInfo("Nbjet", r"$N_{b}$", axis.Regular(5, 0, 5), lambda vg: vg.get_hist('NBjets_medium')),
    # GraphInfo("Nbjet", r"$N_{b}$", axis.Regular(5, 0, 5), lambda vg: (vg['NBjets_medium'], vg['wgt_nobtag']*vg.get_sf("Nominal"))),
    # GraphInfo("Nloose", r"$N_{loose}$", axis.Regular(3, 0, 3), lambda vg: (vg['N_loose_el']+vg['N_loose_mu'], vg.scale)),
    # GraphInfo("j1disc", r"$disc(j_{1})$", axis.Regular(25, 0, 1), lambda vg: vg.Jets.get_hist("discriminator", -1)),
    # GraphInfo("flavor", r'flavor', axis.Regular(7, 0, 7), lambda vg: vg.Jets.get_hist("flavor", -1)),
    # GraphInfo("j2disc", r"$disc(j_{2})$", axis.Regular(25, 0, 1), lambda vg: vg.Jets.get_hist("discriminator", 1)),
    # GraphInfo("j3disc", r"$disc(j_{3})$", axis.Regular(25, 0, 1), lambda vg: vg.Jets.get_hist("discriminator", 2)),
    # GraphInfo("j4disc", r"$disc(j_{4})$", axis.Regular(25, 0, 1), lambda vg: vg.Jets.get_hist("discriminator", 3)),
]

multi = [
    GraphInfo("HT", r"$H_{T}$", axis.Regular(25, 0, 1500), lambda vg : vg.get_hist('HT')),
    GraphInfo("Met", r"Met", axis.Regular(25, 0, 500), lambda vg : vg.get_hist('Met')),
    GraphInfo('nelectrons', '$N_{e}$', axis.Regular(4, 0, 4), lambda vg: (vg.TightElectron.num(), vg.scale)),
    GraphInfo('l1Pt', '$p_{T}(l_{1})$ (GeV)', axis.Regular(25, 0, 500), lambda vg: vg.TightLepton.get_hist('pt', 0)),
    GraphInfo('l2Pt', '$p_{T}(l_{2})$ (GeV)', axis.Regular(25, 0, 200), lambda vg: vg.TightLepton.get_hist('pt', 1)),
    GraphInfo('l3Pt', '$p_{T}(l_{3})$ (GeV)', axis.Regular(25, 0, 200), lambda vg: vg.TightLepton.get_hist('pt', 2)),
    GraphInfo("Njet", r"$N_{j}$", axis.Regular(7, 2, 9), lambda vg : (vg.Jets.num(), vg.scale)),
    GraphInfo("Nbjet", r"$N_{b}$", axis.Regular(4, 1, 5), lambda vg: vg.get_hist('NBjets_medium')),
    GraphInfo("Nloose", r"$N_{loose}$", axis.Regular(3, 0, 3), lambda vg: (vg['N_loose_el']+vg['N_loose_mu'], vg.scale)),
]

ttzCR = [
    GraphInfo("HT", r"$H_{T}$", axis.Regular(25, 0, 1500), lambda vg : vg.get_hist('HT')),
    GraphInfo("Met", r"Met", axis.Regular(25, 0, 500), lambda vg : vg.get_hist('Met')),
    GraphInfo('nelectrons', '$N_{e}$', axis.Regular(4, 0, 4), lambda vg: (vg.TightElectron.num(), vg.scale)),
    GraphInfo('l1Pt', '$p_{T}(l_{1})$ (GeV)', axis.Regular(25, 0, 500), lambda vg: vg.TightLepton.get_hist('pt', 0)),
    GraphInfo('l2Pt', '$p_{T}(l_{2})$ (GeV)', axis.Regular(25, 0, 200), lambda vg: vg.TightLepton.get_hist('pt', 1)),
    # GraphInfo('l3Pt', '$p_{T}(l_{3})$ (GeV)', axis.Regular(25, 0, 200), lambda vg: vg.TightLepton.get_hist('pt', 2)),

    # GraphInfo("mu_mass", "", axis.Regular(25, 75, 105), lambda vg : (vg.dimass('TightMuon', 0, "TightMuon", 1), vg.scale)),
    GraphInfo("Njet", r"$N_{j}$", axis.Regular(8, 0, 8), lambda vg : (vg.Jets.num(), vg.scale)),
    # GraphInfo("ZMass", r"$M_{Z}$", axis.Regular(25, 75, 105), lambda vg : vg.get_hist("Zmass")),
    GraphInfo("Nbjet", r"$N_{b}$", axis.Regular(4, 0, 4), lambda vg: vg.get_hist('NBjets_medium')),
    GraphInfo("Nloose", r"$N_{loose}$", axis.Regular(3, 0, 3), lambda vg: (vg['N_loose_el']+vg['N_loose_mu'], vg.scale)),
    # GraphInfo("", r"", axis.Regular(), ),
    # GraphInfo("", r"", axis.Regular(), ),


]
