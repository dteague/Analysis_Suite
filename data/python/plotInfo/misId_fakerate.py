#!/usr/bin/env python3
import boost_histogram.axis as axis
from analysis_suite.plotting.plotter import GraphInfo

cm_ptbins = axis.Variable([15, 50, 75, 100, 150])
cm_etabins = axis.Variable([0.0, 1.1, 1.479, 1.653, 2.5])

def fake_chargeMisId(vg):
    part = vg['TightElectron']
    mask = part['flip', -1]
    return (ak.flatten(part['pt', -1][mask]), ak.flatten(part['abseta', -1][mask])), ak.flatten(part.scale(-1)[mask])

def fake_pt(vg):
    part = vg['TightElectron']
    mask = part['flip', -1]
    return ak.flatten(part['pt', -1][mask]), ak.flatten(part.scale(-1)[mask])

charge_misId = {
    'Measurement' : [
        # GraphInfo("pt_lead", '$p_{{T}}({}_{{lead}})$', axis.Regular(20, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', 0)),
        # GraphInfo("pt_sub", '$p_{{T}}({}_{{sub}})$', axis.Regular(20, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', 1)),
        # GraphInfo("pt_all", '$p_{{T}}(\ell)$', axis.Regular(20, 0, 200), lambda vg : vg['TightLepton'].get_hist('pt', -1)),
        # GraphInfo("eta_lead", '$\eta({}_{{lead}})$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 0)),
        # GraphInfo("eta_sub", '$\eta({}_{{sub}})$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 1)),
        # GraphInfo("eta_all", '$\eta(\ell)$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', -1)),
        # GraphInfo("mass", '$M({})$', axis.Regular(30, 0, 400), lambda vg : (vg.mass("TightLepton", 0, "TightLepton", 1), vg.scale)),
        # GraphInfo("ht", 'HT', axis.Regular(30, 250, 1000), lambda vg : vg.get_hist("HT")),
        # GraphInfo("met", 'MET', axis.Regular(30, 25, 250), lambda vg : vg.get_hist("Met")),
        GraphInfo("all_fr", 'pteta', (cm_ptbins, cm_etabins), lambda vg : vg['TightElectron'].get_hist2d('pt', 'abseta', -1)),
        GraphInfo("flip_fr", 'pteta', (cm_ptbins, cm_etabins), lambda vg : fake_chargeMisId(vg)),
        GraphInfo("pt_allq", '$p_{{T}}(\ell)$', axis.Regular(37, 15, 200), lambda vg : vg['TightElectron'].get_hist('pt', -1)),
        GraphInfo("pt_flipq", '$p_{{T}}(\ell)$', axis.Regular(37, 15, 200), lambda vg : fake_pt(vg)),
    ],
    'Closure' : [
        GraphInfo("pt_lead", '$p_{{T}}({}_{{lead}})$', axis.Regular(30, 0, 150), lambda vg : vg['TightLepton'].get_hist('pt', 0)),
        GraphInfo("pt_sub", '$p_{{T}}({}_{{sub}})$', axis.Regular(30, 0, 150), lambda vg : vg['TightLepton'].get_hist('pt', 1)),
        GraphInfo("pt_all", '$p_{{T}}(\ell)$', axis.Regular(30, 0, 150), lambda vg : vg['TightLepton'].get_hist('pt', -1)),
        GraphInfo("eta_lead", '$\eta({}_{{lead}})$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 0)),
        GraphInfo("eta_sub", '$\eta({}_{{sub}})$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', 1)),
        GraphInfo("eta_all", '$\eta(\ell)$', axis.Regular(26, -2.5, 2.5), lambda vg : vg['TightLepton'].get_hist('eta', -1)),
        GraphInfo("mass", '$M({})$', axis.Regular(20, 70, 115), lambda vg : (vg.mass("TightLepton", 0, "TightLepton", 1), vg.scale)),
        GraphInfo("ht", 'HT', axis.Regular(30, 0, 250), lambda vg : vg.get_hist("HT")),
        GraphInfo("met", 'MET', axis.Regular(25, 0, 50), lambda vg : vg.get_hist("Met")),
        GraphInfo("metphi", '$\phi(MET)$', axis.Regular(20, -np.pi, np.pi), lambda vg : vg.get_hist("Met_phi")),
        # GraphInfo("lhe_ht", 'HT', axis.Variable([70, 100, 200, 400, 600, 800, 1200]), lambda vg : vg.get_hist("LHE_HT")),

    ],
}
