#!/usr/bin/env python3
import numpy as np
import boost_histogram.axis as axis

from analysis_suite.plotting.plotter import GraphInfo
import analysis_suite.commons.configs as config
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import hep, plot, plot_colorbar
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.info import GroupInfo
import analysis_suite.data.plotInfo.nonprompt_fakerate as pinfo
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.commons.user import workspace_area
from scipy import stats

latex_chan = {"Electron": "e", "Muon": "\mu",
              "EE": 'ee', "EM": 'e\mu', 'MM': '\mu\mu'}

def plot_2d(name, hist, chan, **kwargs):
    part = latex_chan[chan]
    with plot(name) as ax:
        mesh = hist.plot_2d(ax, noText=True)
        ax.set_xlabel(f"$\eta({part})$")
        ax.set_ylabel(f"$\phi({part})$")
        plot_colorbar(mesh, ax)
        ax.plot([-2.5, -1.3], [-1.57, -1.57], color='r')
        ax.plot([-2.5, -1.3], [-0.87, -0.87], color='r')
        ax.plot([-1.3, -1.3], [-0.87, -1.57], color='r')
        ax.plot([-2.5, -2.5], [-0.87, -1.57], color='r')

def plot_hem():
    plot_dir = workspace_area/'HEM_test'
    plot_dir.mkdir(exist_ok=True)

    ntuple = config.get_ntuple('signal_loose')
    groups = {'data': ['data']}
    filename = ntuple.get_filename(year='2018')

    eta_bins = axis.Variable(
        [-2.5, -2.1, -1.7, -1.3, -0.9, -0.5, -0.1, 0.1, 0.5, 0.9, 1.3, 1.7, 2.1, 2.5])
    phi_bins = axis.Variable(
        [-np.pi, -2.96705973, -2.61799388, -2.26892803, -1.91986218, -1.57079633,
         -1.22173048, -0.87266463, -0.52359878, -0.17453293,  0.17453293,
         0.52359878,  0.87266463,  1.22173048,  1.57079633,  1.91986218,
         2.26892803,  2.61799388,  2.96705973, np.pi])
    graphs = [
        GraphInfo("occupancy_Muon", '', (eta_bins, phi_bins),
                  lambda vg : vg['TightMuon'].get_hist2d('eta', 'phi', -1)),
        GraphInfo("occupancy_Electron", '', (eta_bins, phi_bins),
                  lambda vg : vg['TightElectron'].get_hist2d('eta', 'phi', -1)),
        GraphInfo("eta_Muon", '', axis.Regular(25, -2.5, 2.5),
                  lambda vg : vg['TightMuon'].get_hist('eta', -1)),
        GraphInfo("eta_Electron", '', axis.Regular(25, -2.5, 2.5),
                  lambda vg : vg['TightElectron'].get_hist('eta', -1)),
    ]

    chans = ['Electron','Muon']
    plotter = Plotter(filename, groups, ntuple=ntuple, year='2018')

    plotter.fill_hists(graphs)
    for chan in chans:
        hists = plotter.get_hists(f'occupancy_{chan}')
        plot_2d(plot_dir/f'occupancy_{chan}.png', hists['data'], chan)
        plotter.plot_stack(f"eta_{chan}", plot_dir/f'eta_{chan}.png', chan=latex_chan[chan])

if __name__ == '__main__':
    plot_hem()
