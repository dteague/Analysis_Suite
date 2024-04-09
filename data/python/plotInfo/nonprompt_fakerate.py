#!/usr/bin/env python3
import numpy as np
import boost_histogram.axis as axis
from analysis_suite.plotting.plotter import GraphInfo

# pt_binning = axis.Variable([15, 20, 25, 35, 45, 60, 90]) #[15, 20, 25, 35, 55, 100]
# pt_binning = axis.Variable([15, 20, 35, 55, 90])
#[15, 20, 25, 35, 55, 100]
pt_binning = axis.Variable([20, 25, 30, 35, 40, 45, 50, 100])
np_ptbins = {
    "Electron": axis.Variable([20, 25, 40, 50, 100]), # axis.Variable([20, 30, 45, 70]),
    "Muon": axis.Variable([20, 30, 40, 50, 100]), #axis.Variable([20, 35, 70])
}

np_etabins = {
    # "Electron":  axis.Variable([0.0, 2.5]),
    # "Muon": axis.Variable([0.0, 2.4]),
    "Electron": axis.Variable([0.0, 0.8, 1.479, 2.5]),
    "Muon": axis.Variable([0.0, 1.2, 2.1, 2.4]),
}
nonprompt_fake_bins = (np_ptbins['Electron'], np_etabins["Electron"])

op_chan = {"Electron": "Muon", "Muon": "Electron"}

nonprompt = {

    "SideBand_bjet": [
        GraphInfo("mt", '$M_{{T}}({}_{{tight}})$', axis.Regular(16, 0, 200),
                  lambda vg, btag, hlt : (vg['AllLepton']['mt', 0], vg.scale*vg['bjet_scale'][:, btag]*vg[hlt])),
        GraphInfo("pt", '$p_{{T}}({}_{{tight}})$', axis.Regular(30, 0, 150),
                  lambda vg, btag, hlt : (vg['AllLepton']['pt', 0], vg.scale*vg['bjet_scale'][:, btag]*vg[hlt])),
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(30, 0, 150),
                  lambda vg, btag, hlt : (vg["Met"], vg.scale*vg['bjet_scale'][:, btag]*vg[hlt])),
    ],
    "SideBand_nohlt": [
        GraphInfo("mt", '$M_{{T}}({}_{{tight}})$', axis.Regular(16, 0, 200),
                  lambda vg, btag : (vg['AllLepton']['mt', 0], vg.scale*vg['bjet_scale'][:, btag])),
        GraphInfo("pt", '$p_{{T}}({}_{{tight}})$', axis.Regular(30, 0, 150),
                  lambda vg, btag : (vg['AllLepton']['pt', 0], vg.scale*vg['bjet_scale'][:, btag])),
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(30, 0, 150),
                  lambda vg, btag : (vg["Met"], vg.scale*vg['bjet_scale'][:, btag])),

    ],
    "SideBand_none": [
        GraphInfo("mt", '$M_{{T}}({}_{{tight}})$', axis.Regular(25, 0, 250),
                  lambda vg, btag, hlt : (vg['AllLepton']['mt', 0], vg.scale*vg[hlt])),
        GraphInfo("mt_fix", '$M_{{T}}({}_{{tight}})$', axis.Regular(12, 0, 120),
                  lambda vg, btag, hlt : (vg['AllLepton']['mt_fix', 0], vg.scale*vg[hlt])),
        GraphInfo("mass", '$M_{{T}}({}_{{tight}})$', axis.Regular(25, 0, 350),
                  lambda vg, btag, hlt : (vg.mass("AllLepton", 0, 'Jets', 0), vg.scale*vg[hlt])),
        GraphInfo("pt", '$p_{{T}}({}_{{tight}})$', axis.Regular(30, 0, 150),
                  lambda vg, btag, hlt : (vg['AllLepton']['pt', 0], vg.scale*vg[hlt])),
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(30, 0, 150),
                  lambda vg, btag, hlt : (vg["Met"], vg.scale*vg[hlt])),
    ],

    "SideBand": [
        GraphInfo("mt", '$M_{{T}}({}_{{tight}})$', axis.Regular(25, 0, 250),
                  lambda vg : (vg['AllLepton']['mt', 0], vg.scale)),
        GraphInfo("mt_fix", '$M_{{T}}({}_{{tight}})$', axis.Regular(12, 0, 120),
                  lambda vg : (vg['AllLepton']['mt_fix', 0], vg.scale)),
        GraphInfo("pt", '$p_{{T}}({}_{{tight}})$', axis.Regular(30, 0, 150),
                  lambda vg : (vg['AllLepton']['pt', 0], vg.scale)),
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(30, 0, 150),
                  lambda vg : (vg["Met"], vg.scale)),
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
        GraphInfo("j_btag", "", axis.Regular(25, 0, 1), lambda vg : vg["Jets"].get_hist("discriminator", 0)),
        # GraphInfo("l_btag", "", axis.Regular(25, 0, 1), lambda vg : vg["AllLepton"].get_hist("jet_btag", 0)),
        GraphInfo("fr", '', nonprompt_fake_bins, lambda vg : vg['AllLepton'].get_hist2d('pt', 'abseta', 0)),

    ],
    "Measurement_bjet": [
        GraphInfo("pt", '$p_{{T}}({})$', axis.Regular(20, 0, 100), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("pt_fake", '$p_{{T}}({})$', axis.Regular(20, 0, 100), lambda vg : vg['FakeLepton'].get_hist('pt', 0)),
        GraphInfo("pt_tight", '$p_{{T}}({})$', axis.Regular(20, 0, 100), lambda vg : vg['TightLepton'].get_hist('pt', 0)),
        GraphInfo("eta", '$\eta({})$', axis.Regular(25, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),
        GraphInfo("eta_tight", '$\eta({})$', axis.Regular(25, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 0)),
        GraphInfo("rawpt", '$p_{{T}}({})$', axis.Regular(20, 0, 100), lambda vg : vg['AllLepton'].get_hist('rawPt', 0)),
        GraphInfo("rawpt_fake", '$p_{{T}}({})$', axis.Regular(20, 0, 100), lambda vg : vg['FakeLepton'].get_hist('rawPt', 0)),
        GraphInfo("rawpt_tight", '$p_{{T}}({})$', axis.Regular(20, 0, 100), lambda vg : vg['TightLepton'].get_hist('rawPt', 0)),
        GraphInfo("mt", '$M_{{T}}({})$', axis.Regular(30, 0, 150), lambda vg : vg['AllLepton'].get_hist('mt', 0)),
        GraphInfo("mt_fake", '$M_{{T}}({})$', axis.Regular(30, 0, 150), lambda vg : vg['FakeLepton'].get_hist('mt', 0)),
        GraphInfo("mt_tight", '$M_{{T}}({})$', axis.Regular(30, 0, 150), lambda vg : vg['TightLepton'].get_hist('mt', 0)),
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(20, 0, 50), lambda vg : vg.get_hist("Met")),
        GraphInfo("metphi", '$\phi(MET)$', axis.Regular(20, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
        GraphInfo("ptRatio", "ptRatio", axis.Regular(24, 0, 1.2), lambda vg: vg['AllLepton'].get_hist("ptRatio", 0)),
        GraphInfo("ptRatio_fake", "ptRatio", axis.Regular(24, 0, 1.2), lambda vg: vg['FakeLepton'].get_hist("ptRatio", 0)),
        GraphInfo("ptRatio_tight", "ptRatio", axis.Regular(24, 0, 1.2), lambda vg: vg['TightLepton'].get_hist("ptRatio", 0)),
        # GraphInfo("ptRatio2", "ptRatio2", axis.Regular(24, 0, 1.2), lambda vg: vg['FakeLepton'].get_hist("ptRatio2", 0)),
        # GraphInfo("flav_split", "flavor", axis.Regular(25, 0, 1), lambda vg: vg['AllLepton'].get_hist("jet_btag", 0)),
        # GraphInfo("flav_split_tight", "flavor", axis.Regular(25, 0, 1), lambda vg: vg['TightLepton'].get_hist("jet_btag", 0)),
    ],
    "Closure": [
        GraphInfo("pt1", '$p_{{T}}(\ell_{{1}})$', axis.Regular(20, 0, 250), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("pt2", '$p_{{T}}(\ell_{{2}})$', axis.Regular(24, 0, 120), lambda vg : vg['AllLepton'].get_hist('pt', 1)),
        GraphInfo("eta1", '$\eta(\ell_{{1}})$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),
        GraphInfo("eta2", '$\eta(\ell_{{2}})$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 1)),
        GraphInfo("pt", '$p_{{T}}(\ell)$', axis.Regular(20, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("eta", '$\eta(\ell)$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),

        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(20, 0, 300), lambda vg : vg.get_hist("Met")),
        GraphInfo("ht", '$H_T$', axis.Regular(20, 0, 750), lambda vg : vg.get_hist("HT")),
        GraphInfo("ht_old", '$H_T$', axis.Regular(15, 0, 600), lambda vg : vg.get_hist("HT")),
        GraphInfo("metphi", '$\phi(MET)$', axis.Regular(20, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
        GraphInfo("njets", '$N_j$', axis.Regular(6, 0, 6), lambda vg : (vg.Jets.num(), vg.scale)),

        # GraphInfo("flav_split", "flavor", axis.Regular(23, 0, 23), lambda vg: vg['AllLepton'].get_hist("jet_flav", -1)),
        # GraphInfo("flav_split_tight", "flavor", axis.Regular(23, 0, 23), lambda vg: vg['TightLepton'].get_hist("jet_flav", -1)),
    ],
    "FakeClosure": [
        GraphInfo("fake_pt", '$p_{{T}}(\ell)$', axis.Regular(20, 0, 100), lambda vg : vg['FakeLepton'].get_hist('pt', 0)),
        GraphInfo("fake_pt_split", '$p_{{T}}(\ell)$', axis.Variable([20, 25, 35, 50, 100]), lambda vg : vg['FakeLepton'].get_hist('pt', 0)),
        GraphInfo("tight_pt", '$p_{{T}}(\ell)$', axis.Regular(50, 0, 250), lambda vg : vg['TightLepton'].get_hist('pt', 0)),
        GraphInfo("fake_eta", '$p_{{T}}(\ell)$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['FakeLepton'].get_hist('eta', 0)),
        GraphInfo("tight_eta", '$p_{{T}}(\ell)$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 0)),
        GraphInfo("pt", '$p_{{T}}(\ell)$', axis.Regular(20, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("eta", '$\eta(\ell)$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['AllLepton'].get_hist('eta', 0)),
        GraphInfo("fake_rawpt", '$p_{{T}}(\ell)$', axis.Regular(50, 0, 250), lambda vg : vg['FakeLepton'].get_hist('rawPt', 0)),
        GraphInfo("fake_ptRatio", '$p_{{T}}(\ell)$', axis.Regular(50, 0, 1.1), lambda vg : vg['FakeLepton'].get_hist('ptRatio', 0)),
        GraphInfo("tight_rawpt", '$p_{{T}}(\ell)$', axis.Regular(50, 0, 250), lambda vg : vg['TightLepton'].get_hist('rawPt', 0)),
        GraphInfo("tight_ptRatio", '$p_{{T}}(\ell)$', axis.Regular(50, 0, 1.1), lambda vg : vg['TightLepton'].get_hist('ptRatio', 0)),
        # GraphInfo("flav_split", "flavor", axis.Regular(23, 0, 23), lambda vg: vg['AllLepton'].get_hist("jet_flav", -1)),
        # GraphInfo("flav_split_tight", "flavor", axis.Regular(23, 0, 23), lambda vg: vg['TightLepton'].get_hist("jet_flav", -1)),
    ],




    "DY_closure": [
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(12, 0, 120), lambda vg,chan : vg.get_hist("Met")),
        GraphInfo("ht", '$H_T$', axis.Regular(8, 0, 160), lambda vg,chan : vg.get_hist("HT")),
        GraphInfo("mt", '$M_{{T}}({})$', axis.Regular(30, 0, 150), lambda vg,chan : vg[chan].get_hist('mt', 0)),
        GraphInfo("pt", '$p_{{T}}({})$', axis.Regular(12, 0, 120), lambda vg, chan : vg[chan].get_hist('pt', 0)),
        GraphInfo("eta", '$\eta({})$', axis.Regular(12, -2.5, 2.5), lambda vg, chan : vg[chan].get_hist('eta', 0)),
        GraphInfo("mass_z", '$M({0}, {0})$', axis.Regular(20, 70, 110), lambda vg, chan : (vg.mass(op_chan[chan], 0, op_chan[chan], 1), vg.scale)),
        GraphInfo("pt_z1", '$p_{{T}}({}_{{lead}})$', axis.Regular(15, 0, 150), lambda vg, chan : vg[op_chan[chan]].get_hist('pt', 0)),
        GraphInfo("pt_z2", '$p_{{T}}({}_{{sublead}})$', axis.Regular(10, 0, 100), lambda vg, chan : vg[op_chan[chan]].get_hist('pt', 1)),
    ],
}
