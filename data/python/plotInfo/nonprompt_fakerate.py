#!/usr/bin/env python3
import numpy as np
import boost_histogram.axis as axis
from analysis_suite.plotting.plotter import GraphInfo

pt_binning = axis.Variable([15, 20, 30, 55, 70])
# np_ptbins = axis.Variable([15, 20, 25, 35, 50, 70])[15, 20, 25, 35, 50, 70]
np_ptbins = {
    "Electron": pt_binning,  # axis.Variable([15, 20, 25, 35, 45, 60, ]),
    "Muon": pt_binning, #axis.Variable([15, 20, 25, 70]),
}


# np_ptbins = axis.Variable([15, 20, 25, 30, 35, 40, 45, 50, ])
np_etabins = {
    "Electron": axis.Variable([0.0, 1.1, 1.479, 1.653, 2.5]),
    # "Electron": axis.Variable([0.0, 0.8, 1.479, 2.5]),
    "Muon": axis.Variable([0.0, 1.2, 2.1, 2.4]),
}
nonprompt_fake_bins = (np_ptbins['Electron'], np_etabins["Electron"])

ptbins = axis.Variable([10, 15, 25, 35, 50, 70])

nonprompt = {
    'SideBand' : [
        GraphInfo("tightpt", '$p_{{T}}({}_{{tight}})$', axis.Regular(40, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', 0)),
        GraphInfo("tightmt", '$M_{{T}}({}_{{tight}})$', axis.Regular(40, 0, 200), lambda vg : vg['TightLepton'].get_hist('mt', 0)),
        GraphInfo("tighteta", '$|\eta({})|$', axis.Regular(25, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 0)),
        GraphInfo("tightphi", '$|\phi({})|$', axis.Regular(25, -np.pi, np.pi), lambda vg : vg['TightLepton'].get_hist('phi', 0)),
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(40, 0, 200), lambda vg : vg.get_hist("Met")),
        GraphInfo("metphi", '$\phi(MET)$', axis.Regular(40, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
        GraphInfo("mvaTTH", '$mva_{{ttH}}$', axis.Regular(40, -1, 1), lambda vg: vg['TightLepton'].get_hist('mvaTTH', 0)),
    ],
    "Measurement" : [
        GraphInfo("pt", '$p_{{T}}({})$', axis.Regular(40, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        # GraphInfo("pt_fr", '$p_{{T}}({})$', np_ptbins["Electron"], lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        # GraphInfo("mt", '$M_{{T}}({})$', axis.Regular(40, 0, 20), lambda vg : vg['AllLepton'].get_hist('mt', 0)),
        # GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(30, 0, 30), lambda vg : vg.get_hist("Met")),
        # GraphInfo("metphi", '$\phi(MET)$', axis.Regular(20, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
        # GraphInfo("ptrel", '$p_{{T}}^{{rel}}({}, j)$', axis.Regular(40, 0, 20), lambda vg: vg['AllLepton'].get_hist('ptRel', 0)),
        # GraphInfo("ptratio", '$p_{{T}}({})/p_{{T}}(j)$', axis.Regular(48, 0, 1.2), lambda vg: vg['AllLepton'].get_hist('ptRatio', 0)),
        # GraphInfo("mvaTTH", '$mva_{{ttH}}$', axis.Regular(40, -1, 1), lambda vg: vg['AllLepton'].get_hist('mvaTTH', 0)),
        # GraphInfo("iso", "iso", axis.Regular(40, 0, 0.4), lambda vg: vg['AllLepton'].get_hist('iso', 0)),
        GraphInfo("fr", '', nonprompt_fake_bins, lambda vg : vg['AllLepton'].get_hist2d('pt', 'abseta', 0)),
    ],
    'Closure' : [
        GraphInfo("met", '$MET_{{Puppi}}$', axis.Regular(20, 0, 200), lambda vg : vg.get_hist("Met")),
        GraphInfo("ht", '$H_T$', axis.Regular(15, 0, 600), lambda vg : vg.get_hist("HT")),
        GraphInfo("njets", '$N_j$', axis.Regular(8, 0, 8), lambda vg : (vg.Jets.num(), vg.scale)),
        # GraphInfo("pt_tight", '$p_{{T}}(\ell_{{1}})$', axis.Regular(20, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', 0)),
        # GraphInfo("pt_fake", '$p_{{T}}(\ell_{{1}})$', axis.Regular(20, 0, 100), lambda vg : vg['FakeLepton'].get_hist('rawPt', 0)),

        # GraphInfo("eta_tight", '$p_{{T}}(\ell_{{1}})$', axis.Regular(25, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 0)),
        # GraphInfo("eta_fake", '$p_{{T}}(\ell_{{1}})$', axis.Regular(25, -2.5, 2.5), lambda vg : vg['FakeLepton'].get_hist('eta', 0)),
        # GraphInfo("iso_tight", '$p_{{T}}(\ell_{{1}})$', axis.Regular(25, 0, 0.4), lambda vg : vg['TightLepton'].get_hist('iso', 0)),
        # GraphInfo("iso_fake", '$p_{{T}}(\ell_{{1}})$', axis.Regular(25, 0, 0.4), lambda vg : vg['FakeLepton'].get_hist('iso', 0)),
        # GraphInfo("mva_tight", '$p_{{T}}(\ell_{{1}})$', axis.Regular(25, -1, 1), lambda vg : vg['TightLepton'].get_hist('mvaTTH', 0)),
        # GraphInfo("mva_fake", '$p_{{T}}(\ell_{{1}})$', axis.Regular(25, -1, 1), lambda vg : vg['FakeLepton'].get_hist('mvaTTH', 0)),
        # GraphInfo("ptrel_tight", '$p_{{T}}(\ell_{{1}})$', axis.Regular(20, 0, 20), lambda vg : vg['TightLepton'].get_hist('ptRel', 0)),
        # GraphInfo("ptrel_fake", '$p_{{T}}(\ell_{{1}})$', axis.Regular(20, 0, 20), lambda vg : vg['FakeLepton'].get_hist('ptRel', 0)),
        # GraphInfo("ptratio_tight", '$p_{{T}}(\ell_{{1}})$', axis.Regular(24, 0, 1.2), lambda vg : vg['TightLepton'].get_hist('ptRatio', 0)),
        # GraphInfo("ptratio_fake", '$p_{{T}}(\ell_{{1}})$', axis.Regular(24, 0, 1.2), lambda vg : vg['FakeLepton'].get_hist('ptRatio', 0)),

        GraphInfo("pt1", '$p_{{T}}(\ell_{{1}})$', axis.Regular(20, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("pt2", '$p_{{T}}(\ell_{{2}})$', axis.Regular(20, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', 1)),
        GraphInfo("ptrel1", '$p_{{T}}^{{rel}}({}, j)$', axis.Regular(40, 0, 20), lambda vg: vg['AllLepton'].get_hist('ptRel', 0)),
        GraphInfo("ptratio1", '$p_{{T}}({})/p_{{T}}(j)$', axis.Regular(48, 0, 1.2), lambda vg: vg['AllLepton'].get_hist('ptRatio', 0)),
        GraphInfo("mvaTTH1", '$mva_{{ttH}}$', axis.Regular(40, -1, 1), lambda vg: vg['AllLepton'].get_hist('mvaTTH', 0)),
        GraphInfo("iso1", "iso", axis.Regular(40, 0, 0.4), lambda vg: vg['AllLepton'].get_hist('iso', 0)),
        GraphInfo("ptrel2", '$p_{{T}}^{{rel}}({}, j)$', axis.Regular(40, 0, 20), lambda vg: vg['AllLepton'].get_hist('ptRel', 1)),
        GraphInfo("ptratio2", '$p_{{T}}({})/p_{{T}}(j)$', axis.Regular(48, 0, 1.2), lambda vg: vg['AllLepton'].get_hist('ptRatio', 1)),
        GraphInfo("mvaTTH2", '$mva_{{ttH}}$', axis.Regular(40, -1, 1), lambda vg: vg['AllLepton'].get_hist('mvaTTH', 1)),
        GraphInfo("iso2", "iso", axis.Regular(40, 0, 0.4), lambda vg: vg['AllLepton'].get_hist('iso', 1)),

        # GraphInfo("nbjets", '$N_b$', axis.Regular(8, 0, 8), lambda vg : (vg.BJets.num(), vg.scale)),
        # GraphInfo("nloosebjets", '$N_{{bloose}}$', axis.Regular(8, 0, 8), lambda vg : vg.get_hist("N_bloose")),
        # GraphInfo("ntightbjets", '$N_{{btight}}$', axis.Regular(8, 0, 8), lambda vg : vg.get_hist("N_btight")),

        GraphInfo("pt", '$p_{{T}}(\ell)$', axis.Regular(20, 0, 200), lambda vg : vg['FakeLepton'].get_hist('pt', 0)),
        GraphInfo("eta", '$\eta(\ell)$', axis.Regular(20, -2.5, 2.5), lambda vg : vg['FakeLepton'].get_hist('eta', 0)),
    ],
    "Isolation" : [
        GraphInfo("ptrel", '$p_{{T}}^{{rel}}({}, j)$', axis.Regular(40, 0, 20), lambda vg: vg['AllLepton'].get_hist('ptRel', 0)),
        GraphInfo("ptratio", '$p_{{T}}({})/p_{{T}}(j)$', axis.Regular(48, 0, 1.2), lambda vg: vg['AllLepton'].get_hist('ptRatio', 0)),
        GraphInfo("mvaTTH", '$mva_{{ttH}}$', axis.Regular(40, -1, 1), lambda vg: vg['AllLepton'].get_hist('mvaTTH', 0)),
        GraphInfo("iso", "iso", axis.Regular(40, 0, 0.4), lambda vg: vg['AllLepton'].get_hist('iso', 0)),
        GraphInfo("loosemt", '$M_{{T}}({}_{{loose}})$', axis.Regular(20, 0, 200), lambda vg : vg['AllLepton'].get_hist('mt', 0)),
        GraphInfo("pt", '$p_{{T}}({})$', axis.Regular(40, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        GraphInfo("ptcone", '$p_{{T}}^{{cone}}({})$', axis.Regular(40, 0, 200), lambda vg : (vg['AllLepton']['pt', 0]/vg['AllLepton']['ptRatio', 0], vg.scale)),
        GraphInfo("eta", '$\eta({})$', axis.Regular(25, 0, 2.5), lambda vg : vg['AllLepton'].get_hist('abseta', 0)),

        # GraphInfo("loosefr", '', (nonprompt_ptbins, nonprompt_etabins), lambda vg : vg['AllLepton'].get_hist2d('pt', 'abseta', 0)),

        # GraphInfo("pt_stuff", '', (axis.Regular(6, 0.6, 0.9), axis.Regular(10, 0, 10)), lambda vg : vg['AllLepton'].get_hist2d('ptRatio', 'ptRel', 0)),
        # GraphInfo("pt", '$p_{{T}}({})$', axis.Regular(40, 0, 200), lambda vg : vg['AllLepton'].get_hist('pt', 0)),
        # GraphInfo("eta", '$\eta({})$', axis.Regular(25, 0, 2.5), lambda vg : vg['AllLepton'].get_hist('abseta', 0)),

        # GraphInfo("pt", '$p_{{T}}({}_{{tight}})$', nonprompt_ptbins, lambda vg : vg['TightLepton'].get_hist('pt', 0)),
        # GraphInfo("met", '$MET$', axis.Regular(20, 0, 30), lambda vg : vg.get_hist("Met")),
    ],
}
