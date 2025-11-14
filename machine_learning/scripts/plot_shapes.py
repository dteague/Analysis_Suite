#!/usr/bin/env python3
import boost_histogram.axis as axis
import numpy as np

from analysis_suite.plotting.plotter import Plotter
from analysis_suite.plotting.hist_getter import HistGetter,GraphInfo
from analysis_suite.commons.constants import all_eras
import analysis_suite.commons.user as user

from analysis_suite.commons.configs import get_ntuple

years = all_eras

workdir = user.workspace_area/'btag_test'/'split_files'
plot_dir = user.workspace_area/'btag_test'/'train_shapes'
plot_dir.mkdir(exist_ok=True)

ntuple = get_ntuple('bdt')

graphs = {
    'njet': GraphInfo(r'$N_{{j}}$', axis.Regular(12, 0, 12), 'NJets'),

    'nbjet_l': GraphInfo(r'$N_{{looseb}}$', axis.Regular(7, 0, 7), 'NlooseBJets'),
    'nbjet_m': GraphInfo(r'$N_{{mediumb}}$', axis.Regular(6, 0, 6), 'NmediumBJets'),
    'nbjet_t': GraphInfo(r'$N_{{tightb}}$', axis.Regular(5, 0, 5), 'NtightBJets'),

    'nmuon_l': GraphInfo(r'$N_{{loose\mu}}$', axis.Regular(2, 0, 2), 'NlooseMuons'),
    'nelec_l': GraphInfo(r'$N_{{loosee}}$', axis.Regular(2, 0, 2), 'NlooseElectrons'),
    'nmuon': GraphInfo(r'$N_{{\mu}}$', axis.Regular(4, 0, 4), 'NMuons'),
    'nelec': GraphInfo(r'$N_{{e}}$', axis.Regular(4, 0, 4), 'NElectrons'),

    'ht': GraphInfo(r'$H_{{T}}$ (GeV)', axis.Regular(25, 0, 1750), 'HT'),
    'ht_b': GraphInfo(r'$H_{{T}}(b)$ (GeV)', axis.Regular(25, 0, 1000), 'HT_b'),
    'met': GraphInfo('MET (GeV)', axis.Regular(25, 0, 500), 'Met'),
    'centrality': GraphInfo('$H_{{T}}/E_{{tot}}$', axis.Regular(25, 0, 1), 'centrality'),

    'j1Pt': GraphInfo('$p_{{T}}(j_{{1}})$', axis.Regular(25, 0, 750), 'j1Pt'),
    'j2Pt': GraphInfo('$p_{{T}}(j_{{2}})$', axis.Regular(25, 0, 400), 'j2Pt'),
    'j3Pt': GraphInfo('$p_{{T}}(j_{{3}})$', axis.Regular(25, 0, 300), 'j3Pt'),
    'j4Pt': GraphInfo('$p_{{T}}(j_{{4}})$', axis.Regular(25, 0, 200), 'j4Pt'),
    'j5Pt': GraphInfo('$p_{{T}}(j_{{5}})$', axis.Regular(25, 0, 150), 'j5Pt'),

    'j1Disc': GraphInfo('$Disc(j_{{1}})$', axis.Regular(25, 0, 1), 'j1Disc'),
    'j2Disc': GraphInfo('$Disc(j_{{2}})$', axis.Regular(25, 0, 1), 'j2Disc'),
    'j3Disc': GraphInfo('$Disc(j_{{3}})$', axis.Regular(25, 0, 1), 'j3Disc'),
    'j4Disc': GraphInfo('$Disc(j_{{4}})$', axis.Regular(25, 0, 1), 'j4Disc'),
    'j5Disc': GraphInfo('$Disc(j_{{5}})$', axis.Regular(25, 0, 1), 'j5Disc'),

    'b1Pt': GraphInfo('$p_{{T}}(b_{{1}})$', axis.Regular(25, 0, 750), 'b1Pt'),
    'b2Pt': GraphInfo('$p_{{T}}(b_{{2}})$', axis.Regular(25, 0, 300), 'b2Pt'),
    'b3Pt': GraphInfo('$p_{{T}}(b_{{3}})$', axis.Regular(25, 0, 200), 'b3Pt'),
    'b4Pt': GraphInfo('$p_{{T}}(b_{{4}})$', axis.Regular(25, 0, 150), 'b4Pt'),

    'b1Disc': GraphInfo('$Disc(b_{{1}})$', axis.Regular(25, 0, 1), 'b1Disc'),
    'b2Disc': GraphInfo('$Disc(b_{{2}})$', axis.Regular(25, 0, 1), 'b2Disc'),
    'b3Disc': GraphInfo('$Disc(b_{{3}})$', axis.Regular(25, 0, 1), 'b3Disc'),
    'b4Disc': GraphInfo('$Disc(b_{{4}})$', axis.Regular(25, 0, 1), 'b4Disc'),

    'l1Pt': GraphInfo(r'$p_{{T}}(\ell_{{1}})$', axis.Regular(25, 0, 500), 'l1Pt'),
    'l2Pt': GraphInfo(r'$p_{{T}}(\ell_{{2}})$', axis.Regular(25, 0, 200), 'l2Pt'),
    'l3Pt': GraphInfo(r'$p_{{T}}(\ell_{{3}})$', axis.Regular(25, 0, 150), 'l3Pt'),

    'lep_mass': GraphInfo(r'$m(\ell_{{1}}, \ell_{{2}})$', axis.Regular(25, 0, 750), 'lep_mass'),
    'lep_dphi': GraphInfo(r'$\Delta\phi(\ell_{{1}}, \ell_{{2}})$', axis.Regular(25, -np.pi, np.pi), 'lep_dphi'),
    'lep_deta': GraphInfo(r'$\Delta\eta(\ell_{{1}}, \ell_{{2}})$', axis.Regular(25, -4, 4), 'lep_deta'),
    'lep_dr': GraphInfo(r'$\Delta R(\ell_{{1}}, \ell_{{2}})$', axis.Regular(25, 0, 5), 'lep_dr'),
    'lep_mt': GraphInfo(r'$m_{{T}}(\ell_{{1}}, \ell_{{2}})$', axis.Regular(25, 0, 500), 'lep_mt'),
    'lep_cosTheta': GraphInfo(r'$cos(\Delta\theta(\ell_{{1}}, \ell_{{2}}))$', axis.Regular(25, -1, 1), 'lep_cosTheta'),

    'jet_mass': GraphInfo('$m(j_{{1}}, j_{{2}})$', axis.Regular(25, 0, 1500), 'jet_mass'),
    'jet_dphi': GraphInfo('$\Delta\phi(j_{{1}}, j_{{2}})$', axis.Regular(25, -np.pi, np.pi), 'jet_dphi'),
    'jet_deta': GraphInfo('$\Delta\eta(j_{{1}}, j_{{2}})$', axis.Regular(25, -4, 4), 'jet_deta'),
    'jet_dr': GraphInfo('$\Delta R(j_{{1}}, j_{{2}})$', axis.Regular(25, 0, 5), 'jet_dr'),
    'jet_mt': GraphInfo('$m_{{T}}(j_{{1}}, j_{{2}})$', axis.Regular(25, 0, 1000), 'jet_mt'),
    'jet_cosTheta': GraphInfo(r'$cos(\Delta\theta(j_{{1}}, j_{{2}}))$', axis.Regular(25, -1, 1), 'jet_cosTheta'),

    'JetLep1_Cos': GraphInfo(r'$cos(\Delta\phi(j_{{1}}, \ell_{{1}}))$', axis.Regular(25, -1, 1), 'JetLep1_Cos'),
    'JetLep2_Cos': GraphInfo(r'$cos(\Delta\phi(j_{{1}}, \ell_{{2}}))$', axis.Regular(25, -1, 1), 'JetLep2_Cos'),

    'bjet_lep1_dr': GraphInfo('$\Delta R(b_{{1}}, \ell_{{1}})$', axis.Regular(25, 0, 5), 'bjet_lep1_dr'),
    'bjet_lep1_mass': GraphInfo('$m(b_{{1}}, \ell_{{1}})$', axis.Regular(25, 0, 1000), 'bjet_lep1_mass'),
    'bjet_lep2_dr': GraphInfo('$\Delta R(b_{{1}}, \ell_{{2}})$', axis.Regular(25, 0, 5), 'bjet_lep2_dr'),
    'bjet_lep2_mass': GraphInfo('$m(b_{{1}}, \ell_{{2}})$', axis.Regular(25, 0, 750), 'bjet_lep2_mass'),

    'mT_1': GraphInfo('$m_{{T}}(\ell_{{1}}, p_{{T}}^{{miss}})$', axis.Regular(25, 0, 500), 'mT_1'),
    'mT_2': GraphInfo('$m_{{T}}(\ell_{{2}}, p_{{T}}^{{miss}})$', axis.Regular(25, 0, 400), 'mT_2'),
    'mT2_l': GraphInfo('$m_{{T2}}(\ell_{{1}},\ell_{{2}}, p_{{T}}^{{miss}})$', axis.Regular(25, 0, 150), 'mT2_l'),
    'cosdphi_1': GraphInfo(r'$cos(\Delta\phi(\ell_{{1}}, p_{{T}}^{{miss}}))$', axis.Regular(25, -1, 1), 'cosdphi_1'),
    'cosdphi_2': GraphInfo(r'$cos(\Delta\phi(\ell_{{2}}, p_{{T}}^{{miss}}))$', axis.Regular(25, -1, 1), 'cosdphi_2'),
}

shapes = {
    'Three Top': ['ttt_nlo'],
    'Four Top': ['4top'],
    'Background': None
}

plotter = Plotter(get_ntuple('bdt'), 'all', sig='4top', bkg='all', outdir=plot_dir)
# for year in years:
year = "all"
plotter.set_year(year)
filename = workdir/year/'test_Nominal_signal.root'
# filename = workdir/'train_Nominal_signal.root'
hist_factory = HistGetter(ntuple, year, filename=filename)
for name, graph in graphs.items():
    print(name)
    hist = hist_factory.get_hist(graph)
    plotter.plot_hist(name, hist, plot_type='shape', shapes=shapes)
