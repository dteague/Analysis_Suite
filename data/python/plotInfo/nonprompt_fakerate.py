#!/usr/bin/env python3
import numpy as np
import boost_histogram.axis as axis
from analysis_suite.plotting.plotter import GraphInfo

pt_binning = axis.Variable([15, 20, 25, 35, 45, 60, 90]) #[15, 20, 25, 35, 55, 100]
np_ptbins = {
    "Electron": pt_binning,
    "Muon": pt_binning,
}

np_etabins = {
    "Electron": axis.Variable([0.0, 0.8, 1.479, 2.5]), # axis.Variable([0.0, 1.1, 1.479, 1.653, 2.5]),
    "Muon": axis.Variable([0.0, 1.2, 2.1, 2.4]),
}
nonprompt_fake_bins = (np_ptbins['Electron'], np_etabins["Electron"])

op_chan = {"Electron": "Muon", "Muon": "Electron"}

mt = GraphInfo("mt", '$M_{{T}}({})$', axis.Regular(20, 0, 100), lambda vg : vg['AllLepton'].get_hist('mt', 0))

nonprompt = {
    'SideBand' : [
        GraphInfo("tight_low_pt", '$p_{{T}}({}_{{tight}})$', axis.Regular(12, 0, 30), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("tight_high_pt", '$p_{{T}}({}_{{tight}})$', axis.Regular(30, 0, 150), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("tight_low_mt", '$M_{{T}}({}_{{tight}})$', axis.Regular(10, 0, 100), lambda vg : vg['AllLepton'].get_hist('mt', 0)),
        GraphInfo("tight_high_mt", '$M_{{T}}({}_{{tight}})$', axis.Regular(30, 0, 150), lambda vg : vg['AllLepton'].get_hist('mt', 0)),
        GraphInfo("tighteta", '$|\eta({})|$', axis.Regular(25, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),
        GraphInfo("tightphi", '$\phi({})$', axis.Regular(25, -np.pi, np.pi), lambda vg : vg['AllLepton'].get_hist('phi', 0)),
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(30, 0, 150), lambda vg : vg.get_hist("Met")),
        GraphInfo("metphi", '$\phi(MET)$', axis.Regular(20, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
        GraphInfo("iso", '$iso({})$', axis.Regular(20, 0, 0.4), lambda vg : vg['AllLepton'].get_hist('iso', 0)),
        # GraphInfo("mvaTTH", '$mva_{{ttH}}$', axis.Regular(40, -1, 1), lambda vg: vg['TightLepton'].get_hist('mvaTTH', 0)),
    ],
    "Measurement" : [
        GraphInfo("pt", '$p_{{T}}({})$', axis.Regular(20, 0, 100), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("pt_fr", '$p_{{T}}({})$', pt_binning, lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("mt", '$M_{{T}}({})$', axis.Regular(30, 0, 150), lambda vg : vg['AllLepton'].get_hist('mt', 0)),
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(20, 0, 50), lambda vg : vg.get_hist("Met")),
        GraphInfo("ht", '$H_T$', axis.Regular(30, 0, 300), lambda vg : vg.get_hist("HT")),
        GraphInfo("metphi", '$\phi(MET)$', axis.Regular(20, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
        GraphInfo("njets", '$N_j$', axis.Regular(6, 0, 6), lambda vg : (vg.Jets.num(), vg.scale)),
        GraphInfo("iso", '$iso({})$', axis.Regular(20, 0, 0.4), lambda vg : vg['AllLepton'].get_hist('iso', 0)),
        GraphInfo("fr", '', nonprompt_fake_bins, lambda vg : vg['AllLepton'].get_hist2d('pt', 'abseta', 0)),
    ],
    'Closure' : [
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(20, 0, 200), lambda vg : vg.get_hist("Met")),
        GraphInfo("ht", '$H_T$', axis.Regular(15, 0, 600), lambda vg : vg.get_hist("HT")),
        GraphInfo("njets", '$N_j$', axis.Regular(8, 0, 8), lambda vg : (vg.Jets.num(), vg.scale)),
        GraphInfo("nbjets", '$N_b$', axis.Regular(5, 0, 5), lambda vg : (vg.BJets.num(), vg.scale)),

        GraphInfo("pt1", '$p_{{T}}(\ell_{{1}})$', axis.Regular(20, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("pt2", '$p_{{T}}(\ell_{{2}})$', axis.Regular(12, 0, 120), lambda vg : vg['AllLepton'].get_hist('pt', 1)),
        GraphInfo("eta1", '$\eta(\ell_{{1}})$', axis.Regular(10, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),
        GraphInfo("eta2", '$\eta(\ell_{{2}})$', axis.Regular(10, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 1)),
        GraphInfo("pt", '$p_{{T}}(\ell)$', axis.Regular(20, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("eta", '$\eta(\ell)$', axis.Regular(10, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),
    ],

    "DY_closure": [
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(12, 0, 120), lambda vg,chan : vg.get_hist("Met")),
        GraphInfo("ht", '$H_T$', axis.Regular(8, 0, 160), lambda vg,chan : vg.get_hist("HT")),
        GraphInfo("pt", '$p_{{T}}({})$', axis.Regular(12, 0, 120), lambda vg, chan : vg[chan].get_hist('pt', 0)),
        GraphInfo("eta", '$\eta({})$', axis.Regular(12, -2.5, 2.5), lambda vg, chan : vg[chan].get_hist('eta', 0)),
        GraphInfo("mass_z", '$M({0}, {0})$', axis.Regular(20, 70, 110), lambda vg, chan : (vg.mass(op_chan[chan], 0, op_chan[chan], 1), vg.scale)),
        GraphInfo("pt_z1", '$p_{{T}}({}_{{lead}})$', axis.Regular(15, 0, 150), lambda vg, chan : vg[op_chan[chan]].get_hist('pt', 0)),
        GraphInfo("pt_z2", '$p_{{T}}({}_{{sublead}})$', axis.Regular(10, 0, 100), lambda vg, chan : vg[op_chan[chan]].get_hist('pt', 1)),
    ],
}
