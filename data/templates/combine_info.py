from analysis_suite.plotting.hist_getter import GraphInfo
import boost_histogram.axis as axis
from .params import cut

dilep_bins = [0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.9]
multi_bins = [0, 0.24, 0.36, 0.48, 0.6, 0.9]
ht_bins = [250, 275, 300, 350, 400, 450, 550, 700, 850]

dilep_graph = GraphInfo(r'$BDT_{{3top}}$', axis.Variable(dilep_bins), 'BDT')
multi_graph = GraphInfo(r'$BDT_{{3top}}$',axis.Variable(multi_bins), 'BDT')
top4_graph = GraphInfo(r'$BDT_{{4top}}$',axis.Regular(1, 0, 1), 'BDT')
ttz_graph = GraphInfo(r'$H_{{T}}$ (GeV)', axis.Variable(ht_bins), "HT")

rate_params = ['ttz']

regions = {
    "Dilepton": {
        'dir': 'split_files',
        'glob': 'bdt_3top_sig*signal.root',
        'graph': dilep_graph,
        'mask': lambda vg : (vg['NLeps'] == 2),
        'ntuple': ('signal', 'dilep_ntuple'),
    },
    "Multi": {
        'dir': 'split_files',
        'glob': 'bdt_3top_sig*signal.root',
        'graph': multi_graph,
        'mask': lambda vg : (vg['NLeps']>=3),
        'ntuple': ('signal', 'multi_ntuple'),
    },
    "ttzCR": {
        'dir': '',
        'glob': 'processed*ttzCR.root',
        'year_split': True,
        'graph': ttz_graph,
        'mask': None,
        'ntuple': ('ttzCR', 'info'),
    },
    "ttttCR": {
        'dir': 'split_files',
        'glob': 'bdt_4top_sig*signal.root',
        'graph': top4_graph,
        'mask': lambda vg : vg['BDT'] > cut,
        'ntuple': ('signal', 'dilep_ntuple'),
    },
}
