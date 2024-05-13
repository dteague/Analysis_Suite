from analysis_suite.plotting.plotter import GraphInfo
import boost_histogram.axis as axis

sig_bins = [0.0, 0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.94, 1.0]
dilep_graph = GraphInfo('sig_3top', '', axis.Variable(sig_bins), '3top_sig')
multi_graph = GraphInfo('sig_3top', '', axis.Variable(sig_bins), '3top_sig')
top4_graph = GraphInfo('sig_4top', '', axis.Regular(1, 0, 1), '4top_sig')
ttz_graph = GraphInfo('HT', '', axis.Regular(20, 250, 750), "HT")

rate_params = ['ttz']

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
    "ttttCR": {
        'dir': 'first_train',
        'glob': 'test*',
        'graph': top4_graph,
        'mask': lambda vg : vg['4top_sig'] > 0.97
    },
    "ttzCR": {
        'dir': '',
        'glob': 'processed*ttzCR.root',
        'graph': ttz_graph,
        'mask': None
    },
}
