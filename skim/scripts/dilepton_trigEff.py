#!/usr/bin/env python3
from pathlib import Path
import boost_histogram.axis as axis
import numpy as np
from analysis_suite.commons.info import NtupleInfo
from scipy.stats import beta
import gzip

from correctionlib.schemav2 import VERSION, MultiBinning, Category, Correction, CorrectionSet
import correctionlib._core as core

import analysis_suite.commons.user as user
from analysis_suite.plotting.plotter import GraphInfo
import analysis_suite.commons.configs as config
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import hep, plot, plot_colorbar
from analysis_suite.commons.constants import lumi, all_eras
from analysis_suite.commons.info import GroupInfo
import analysis_suite.data.plotInfo.nonprompt_fakerate as pinfo
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.commons.plot_utils import plot, nonratio_plot, ratio_plot, hep

info = NtupleInfo(
    filename = user.hdfs_area / 'workspace/dilep_trigeff/{year}/{workdir}/',
    trees = ['Signal'],
    branches = [
        lambda vg : vg.mergeParticles("AllLeptons", "TightMuon", "TightElectron")
    ],
    color_by_group = {
        'ttbar_lep': 'black',
        'data': 'black',
    }
)


trigger_branch = "Dilepton_trigger"
chans = ['EE', 'EM', 'MM']
mc = 'ttbar_lep'

phi_bins = np.linspace(-np.pi, np.pi, 13)
extra_graphs = {
    "Met":  GraphInfo("met", 'MET (GeV)', axis.Variable([130, 150, 175, 200, 225, 250, 275, 300, 325, 350]),
                      lambda vg: vg.get_hist('Met')),
    "MHT":  GraphInfo("MHT", 'MHT (GeV)', axis.Variable([100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350]),
                      lambda vg: vg.get_hist('MHT')),
}


all_graphs = {
    "EE": {
        "pt2d": GraphInfo("pt2d", '', (axis.Variable([25, 50, 75, 100, 150, 200]), axis.Variable([15, 30, 50, 75, 100, 150, 200])),
                          lambda vg : ((vg['TightElectron']['pt', 0], vg['TightElectron']['pt', 1]), vg.scale)),
        "leadpt":  GraphInfo("leadpt", '$p_{{T}}(e_{{lead}})$ (GeV)', axis.Variable([25, 50, 75, 100, 150, 200]),
                         lambda vg: vg["TightElectron"].get_hist('pt', 0)),
        "subpt":  GraphInfo("subpt", '$p_{{T}}(e_{{sub}})$ (GeV)', axis.Variable([15, 30, 50, 75, 100, 150, 200]),
                         lambda vg: vg["TightElectron"].get_hist('pt', 1)),
        "leadeta":  GraphInfo("leadeta", '$|\eta(e_{{lead}})|$', axis.Variable([-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2, 2.5]),
                         lambda vg: vg["AllLeptons"].get_hist('eta', 0)),
        "subeta":  GraphInfo("subeta", '$|\eta(e_{{sub}})|$', axis.Variable([-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2, 2.5]),
                         lambda vg: vg["AllLeptons"].get_hist('eta', 1)),
        "leadphi":  GraphInfo("leadphi", '$\phi(e_{{lead}})$', axis.Variable(phi_bins),
                         lambda vg: vg["TightElectron"].get_hist('phi', 0)),
        "subphi":  GraphInfo("subphi", '$\phi(e_{{sub}})$', axis.Variable(phi_bins),
                         lambda vg: vg["TightElectron"].get_hist('phi', 1)),
        **extra_graphs
    },
    'EM': {
        "pt2d": GraphInfo("pt2d", '', (axis.Variable([15, 30, 50, 75, 100, 150, 200]), axis.Variable([15, 30, 50, 75, 100, 150, 200])),
                          lambda vg : ((vg['TightElectron']['pt', 0], vg['TightMuon']['pt', 0]), vg.scale)),

        "elpt":  GraphInfo("leadpt", '$p_{{T}}(e)$ (GeV)', axis.Variable([15, 30, 50, 75, 100, 150, 200]),
                         lambda vg: vg["TightElectron"].get_hist('pt', 0)),
        "mupt":  GraphInfo("subpt", '$p_{{T}}(\mu)$ (GeV)', axis.Variable([15, 30, 50, 75, 100, 150, 200]),
                         lambda vg: vg["TightMuon"].get_hist('pt', 0)),
        "eleta":  GraphInfo("el_eta", '$|\eta(e)|$', axis.Variable([-2.5, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2, 2.5]),
                         lambda vg: vg["TightElectron"].get_hist('eta', 0)),
        "mueta":  GraphInfo("mu_eta", '$|\eta(\mu)|$', axis.Variable([-2.4, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2, 2.4]),
                         lambda vg: vg["TightMuon"].get_hist('eta', 0)),
        "elphi":  GraphInfo("el_phi", '$\phi(e)$', axis.Variable(phi_bins),
                         lambda vg: vg["TightElectron"].get_hist('phi', 0)),
        "muphi":  GraphInfo("mu_phi", '$\phi(\mu)$', axis.Variable(phi_bins),
                         lambda vg: vg["TightMuon"].get_hist('phi', 0)),
        **extra_graphs
    },
    "MM": {
        "pt2d": GraphInfo("pt2d", '', (axis.Variable([25, 50, 75, 100, 150, 200]), axis.Variable([15, 30, 50, 75, 100, 150, 200])),
                          lambda vg : ((vg['TightMuon']['pt', 0], vg['TightMuon']['pt', 1]), vg.scale)),
        "leadpt":  GraphInfo("leadpt", '$p_{{T}}(\mu_{{lead}})$ (GeV)', axis.Variable([25, 50, 75, 100, 150, 200]),
                         lambda vg: vg["TightMuon"].get_hist('pt', 0)),
        "subpt":  GraphInfo("subpt", '$p_{{T}}(\mu_{{sub}})$ (GeV)', axis.Variable([15, 30, 50, 75, 100, 150, 200]),
                         lambda vg: vg["TightMuon"].get_hist('pt', 1)),
        "leadeta":  GraphInfo("leadeta", '$|\eta(\mu_{{lead}})|$', axis.Variable([-2.4, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2, 2.4]),
                         lambda vg: vg["TightMuon"].get_hist('eta', 0)),
        "subeta":  GraphInfo("subeta", '$|\eta(\mu_{{sub}})|$', axis.Variable([-2.4, -2.0, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2, 2.4]),
                         lambda vg: vg["TightMuon"].get_hist('eta', 1)),
        "leadphi":  GraphInfo("leadphi", '$\phi(\mu_{{lead}})$', axis.Variable(phi_bins),
                         lambda vg: vg["TightMuon"].get_hist('phi', 0)),
        "subphi":  GraphInfo("subphi", '$\phi(\mu_{{sub}})$', axis.Variable(phi_bins),
                         lambda vg: vg["TightMuon"].get_hist('phi', 1)),
        **extra_graphs
    }
}

def get_totalEff(passhist, allhist):
    alf = (1-0.682689492137)/2
    top = np.sum(passhist.vals)
    bot = np.sum(allhist.vals)
    boterr2 = np.sum(allhist.sumw2)

    p = top*bot/(boterr2+1e-6)+1
    t = (bot-top)*bot/(boterr2+1e-6)+1

    lo, eff, hi = beta.ppf(alf, p, t), beta.mean(p, t), beta.ppf(1 - alf, p, t)

    error2 = ((eff-lo)**2 + (hi-eff)**2)/2
    return eff, np.sqrt(error2)

def getEff(year):

    filename = info.get_filename(year=year)
    print(filename)

    groups = info.get_info().setup_groups([mc, 'data'])

    plotter = Plotter(filename, groups, ntuple=info, year=year)
    plotter.set_groups(bkg=[mc])
    hists = {'all': dict(), 'pass': dict(),}

    for chan in chans:
        graphs = all_graphs[chan]
        print(chan)
        plotter.mask(lambda vg : vg['TightMuon'].num() == chan.count("M"))
        plotter.fill_hists(graphs)
        total_hists = {name: plotter.get_hists(name) for name in graphs.keys()}
        plotter.mask(lambda vg: vg[trigger_branch] == 1, clear=False)
        plotter.fill_hists(graphs)
        pass_hists = {name: plotter.get_hists(name) for name in graphs.keys()}
        hists['pass'][chan] = pass_hists
        hists['all'][chan] = total_hists
        for key, hist in hists['all'][chan].items():
            hist['data'].set_plot_details(['Pass Met Trigger', 'r'])
            hist[mc].set_plot_details(['Pass Met Trigger', 'r'])
        for hist in hists['pass'][chan].values():
            hist['data'].set_plot_details(['Pass Dilepton Trigger', 'b'])
            hist[mc].set_plot_details(['Pass Dilepton Trigger', 'b'])

    return hists

def plot_1d(plotname, h_pass, h_all, ratio, graph, effs=None, ylim=None):
    with ratio_plot(plotname, graph.axis_name, graph.edges(), ratio_bot=0.93, ratio_top=1.07, pad_label='Efficiency',
                    subpad_label='Scale Factor', zero_bot=ylim is None) as ax:
        pad, subpad = ax

        h_pass.color = 'b'
        h_pass.plot_points(pad)
        h_all.color = 'r'
        h_all.plot_points(pad)
        if effs is None:
            eff_val, err_val = get_totalEff(h_pass, h_all)
        else:
            eff_val, err_val = effs
        if ylim is not None:
            pad.set_ylim(*ylim)

        ratio.color = 'b'
        ratio.plot_points(subpad)
        subpad.plot([graph.edges()[0], graph.edges()[-1]], [eff_val, eff_val], color='r')
        subpad.text((graph.edges()[-1]+3*graph.edges()[0])/4, 1.02, f'${round(eff_val, 3)}\pm{round(err_val, 3)}$')
        hep.cms.label(ax=pad, lumi=lumi[year], data=True)
    return eff_val, err_val

def build_pt2d(sf, syst):
    npt1, npt2 = sf.axes.size
    pt1bins = sf.axes[0].edges.tolist()
    pt2bins = sf.axes[1].edges.tolist()
    pt1bins[0] = 15.0
    pt2bins[0] = 15.0
    pt1bins[-1] = "Infinity"
    pt2bins[-1] = "Infinity"

    vals = sf.values().flatten()
    # if syst == "up":
    #     vals = vals + np.sqrt(sf.variances().flatten())
    # elif syst == "down":
    #     vals = vals - np.sqrt(sf.variances().flatten())

    return MultiBinning.parse_obj({
        "nodetype": "multibinning",
        "inputs": ["leading_pt", "subleading_pt"],
        "edges": [pt1bins, pt2bins],
        "content": list(vals),
        "flow": "error",
    })


def build_corr(name, desc, sf):
    # systs = ["nom", "up", "down"]
    return Correction.parse_obj({
        "version": 1,
        "name": name,
        "description": desc,
        "inputs": [
            # {"name": "systematic", "type": "string",
            #  "description": "Central value and shifts (statistical only)"},
            {"name": "leading_pt", "type": "real", "description": "Leading Pt"},
            {"name": "subleading_pt", "type": "real", "description": "Subleading Pt"},
        ],
        "output": {"name": "weight", "type": "real", "description": "Trigger Scale Factor"},
        "data": build_pt2d(sf, 'nom'),
    })



def make_scalefactor(year, ee_scale, em_scale, mm_scale):
    if "2016" in year:
        year += "VFP"
    cset = CorrectionSet.parse_obj({
        "schema_version": VERSION,
        "corrections": [
            build_corr("EE", "MM trigger scale factors", ee_scale),
            build_corr("EM", "MM trigger scale factors", em_scale),
            build_corr("MM", "MM trigger scale factors", mm_scale),
        ],
    })

    outdir = user.analysis_area / 'data/POG/USER/'
    with gzip.open(outdir/f"{year}_UL"/"my_trigger_sf.json.gz", "wt") as fout:
        fout.write(cset.json(exclude_unset=True, indent=4))


workdir = user.workspace_area / 'dilepton_trigger'
workdir.mkdir(exist_ok=True)

latex_name = {"EE": 'e', "MM": '\mu'}

for year in all_eras:
    hists = getEff(year)
    scales = {}
    for chan in chans:
        graph_names = all_graphs[chan]
        chan_dir = workdir / year / chan
        chan_dir.mkdir(exist_ok=True, parents=True)
        for name, graph in graph_names.items():
            allhists = hists['all'][chan][name]
            passhists = hists['pass'][chan][name]
            if '2d' in name:
                mc_eff = Histogram.efficiency(passhists[mc], allhists[mc])
                data_eff = Histogram.efficiency(passhists['data'], allhists['data'])
                scale = data_eff/mc_eff
                scales[chan] = scale.hist
                with plot(chan_dir/f'2d_scale_{chan}.png') as ax:
                    if chan == 'EM':
                        ax.set_xlabel(f"$p_{{T}}({latex_name['EE']})$ [GeV]")
                        ax.set_ylabel(f"$p_{{T}}({latex_name['MM']})$ [GeV]")
                    else:
                        ax.set_xlabel(f"$p_{{T}}(lead {latex_name[chan]})$ [GeV]")
                        ax.set_ylabel(f"$p_{{T}}(sub {latex_name[chan]})$ [GeV]")
                        graph_dim = len(graph.edges())
                        for i in range(0, graph_dim):
                            for j in range(i+2, graph_dim):
                                scale.hist[i,j] = (0, 0)
                    mesh = scale.plot_2d(ax)

                    plot_colorbar(mesh, ax)
                continue
            continue
            mc_eff = Histogram.efficiency(passhists[mc], allhists[mc])
            mc_eff.set_plot_details(["MC Efficiency", 'b'])
            data_eff = Histogram.efficiency(passhists['data'], allhists['data'])
            data_eff.set_plot_details(["Data Efficiency", 'b'])
            scale = data_eff/mc_eff

            mc_effi, mc_err = plot_1d(chan_dir/f'{name}_mcEff_{chan}.png', passhists[mc], allhists[mc], mc_eff, graph)
            data_effi, data_err = plot_1d(chan_dir/f'{name}_dataEff_{chan}.png', passhists['data'], allhists['data'], data_eff, graph)
            tot_scale = data_effi/mc_effi
            scale_err = tot_scale*np.sqrt((data_err/data_effi)**2+(mc_err/mc_effi)**2)
            plot_1d(chan_dir/f'{name}_scale_{chan}.png', data_eff, mc_eff, scale, graph, (tot_scale, scale_err), ylim=(0.8, 1.1))

    make_scalefactor(year, scales['EE'], scales['EM'], scales['MM'])
