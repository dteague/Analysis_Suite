#!/usr/bin/env python3
import argparse
import numpy as np
import awkward as ak
import boost_histogram.axis as axis
from contextlib import contextmanager
import multiprocessing as mp

import warnings
warnings.simplefilter("ignore", UserWarning)

from analysis_suite.plotting.plotter import Plotter
from analysis_suite.plotting.full_getter import GraphInfo, YearGetter
import analysis_suite.commons.configs as config
import analysis_suite.commons.user as user

def plot_2d(filename, vals, graph, year):
    xbins, ybins = graph.edges(0), graph.edges(1)
    with plot(filename) as ax:
        xx = np.tile(xbins, (len(ybins)-1, 1))
        yy = np.tile(ybins, (len(xbins)-1, 1)).T
        mesh = ax.pcolormesh(xx, yy, vals, shading='flat')
        plot_colorbar(mesh, ax)
        ax.set_xscale('log')
        cms_label(ax, year=year)

@contextmanager
def plot_wrapper(hist_factory, plotter, graphs):
    hist_factory.reset_mask()
    yield
    hists = hist_factory.get_hists(graphs, flow=True)
    plotter.plot_hists(hists)

def mask_chan(hist_factory, plotter, chan, subchan, extra):
    plotter.set_extra_text(f'{subchan}_{chan}_{extra}')
    if chan == "MM":
        hist_factory.mask(lambda vg : vg["TightElectron"].num() == 0 )
    elif chan == "EE":
        hist_factory.mask(lambda vg : vg["TightMuon"].num() == 0 )
    else:
        hist_factory.mask(lambda vg : vg["TightElectron"].num()*vg["TightMuon"].num() > 0 )


def plot(year, workdir):
    graphs = {
        'eta_b': GraphInfo('eta', axis.Regular(25, -2.4, 2.4), lambda vg: vg.BJets.get_hist('eta', -1)),
        'pt_b': GraphInfo('pt', axis.Regular(25, 0, 500), lambda vg: vg.BJets.get_hist('pt', -1)),
        'phi_b': GraphInfo('phi', axis.Regular(25, -3.14, 3.14), lambda vg: vg.BJets.get_hist('phi', -1)),
    }

    ntuple = config.get_ntuple('signal', 'multi_ntuple')
    ntuple.remove_group('nonprompt_mc')
    plotter = Plotter(ntuple, year, bkg='all', outdir=workdir/year)

    scaledir = user.workspace_area/'btag_test'
    hist_factory = YearGetter(ntuple, year,
                              scales=['btag_jetlep', 'wz', 'theory_rescale'],
                              workdir=scaledir)

    with plot_wrapper(hist_factory, plotter, graphs):
        hist_factory.mask(lambda vg : vg['N_loose_el']+vg['N_loose_mu'] == 0)
        plotter.set_extra_text(f'all')

    with plot_wrapper(hist_factory, plotter, graphs):
        hist_factory.mask(lambda vg : vg['N_loose_el']+vg['N_loose_mu'] == 0)
        hist_factory.mask(lambda vg: vg['passZVeto'])
        plotter.set_extra_text(f'multi')

    with plot_wrapper(hist_factory, plotter, graphs):
        hist_factory.mask(lambda vg : vg['N_loose_el']+vg['N_loose_mu'] == 0)
        hist_factory.mask(lambda vg: ~vg['passZVeto'])
        plotter.set_extra_text(f'ttz')

if __name__ == "__main__":
    years = ['all']
    # years = ['2016pre', '2016post', '2017', '2018']
    workdir = user.analysis_area/'preselection_test"

    for year in years:
        plot(year, workdir)
