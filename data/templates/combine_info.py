from analysis_suite.plotting.plotter import GraphInfo
import boost_histogram.axis as axis
from .params import cut

dilep_bins = [0.0, 0.1, 0.15, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.94,]
multi_bins = [0.0, 0.3, 0.45, 0.6, 0.94]
dilep_graph = GraphInfo('sig_3top', '', axis.Variable(dilep_bins), '3top_sig')
multi_graph = GraphInfo('sig_3top', '', axis.Variable(multi_bins), '3top_sig')
top4_graph = GraphInfo('sig_4top', '', axis.Regular(1, 0, 1), '4top_sig')
ttz_graph = GraphInfo('HT', '', axis.Regular(20, 250, 750), "HT")

met_graph = GraphInfo('met', 'MET', axis.Regular(15, 50, 400), "Met")
njet_graph = GraphInfo('njet', '$N_{j}$', axis.Regular(7, 2, 9), "NJets")
jetpt_graph = GraphInfo('jetpt', '$p_{T}(j_{1})$', axis.Regular(15, 40, 500), "j1Pt")
nelec_graph = GraphInfo('nElectrons', '$N_{e}$', axis.Regular(3, 0, 3), "NElectrons")

rate_params = ['ttz']
# rate_params = []
regions = {
    "Dilepton": {
        'dir': 'second_train',
        'glob': 'test*',
        'graph': dilep_graph,
        'mask': lambda vg : vg['NMuons']+vg['NElectrons'] == 2,
    },
    "Multi": {
        'dir': 'second_train',
        'glob': 'test*',
        'graph': multi_graph,
        'mask': lambda vg : vg['NMuons']+vg['NElectrons'] > 2,
    },
    "ttzCR": {
        'dir': '',
        'glob': 'processed*ttzCR.root',
        'graph': ttz_graph,
        'mask': None
    },
    "ttttCR": {
        'dir': 'first_train',
        'glob': 'test*',
        'graph': top4_graph,
        'mask': lambda vg : vg['4top_sig'] > cut,
    },
    # 'nelecTest': {
    #     'dir': "",
    #     'glob': 'processed*signal.root',
    #     'graph': nelec_graph,
    #     'mask': None
    # },
    # 'jetTest': {
    #     'dir': "",
    #     'glob': 'processed*signal.root',
    #     'graph': njet_graph,
    #     'mask': None
    # },
    # 'metTest': {
    #     'dir': "",
    #     'glob': 'processed*signal.root',
    #     'graph': met_graph,
    #     'mask': None
    # },
    # "jetptTest": {
    #     'dir': "",
    #     'glob': 'processed*signal.root',
    #     'graph': jetpt_graph,
    #     'mask': None
    # },
}
