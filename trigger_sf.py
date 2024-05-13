#!/usr/bin/env python3
import argparse
import numpy as np
from datetime import datetime

from analysis_suite.commons.histogram import Histogram
import analysis_suite.commons.configs as config
from analysis_suite.commons.constants import lumi
from analysis_suite.commons.info import GroupInfo
from analysis_suite.plotting.plotter import Plotter
from analysis_suite.commons.user import workspace_area

import boost_histogram.axis as axis
from analysis_suite.plotting.plotter import GraphInfo

graphs = [
    GraphInfo("pt", '$p_{{T}}({})$', axis.Regular(10, 10, 30), lambda vg : vg['AllLepton'].get_hist('rawPt', 0)),
]

latex_chan = {"Electron": "e", "Muon": "\mu",
              "EE": 'ee', "EM": 'e\mu', 'MM': '\mu\mu'}

name_chan = {"Electron": "el", "Muon": "mu"}

def scale_prescale(vg, hlt, clear=False):
    if clear:
        vg.scale = 1
    else:
        vg.scale = vg[hlt]

def trigger_eff(workdir, ginfo, year, input_dir):
    groups = ginfo.setup_groups(["data", "wjet_ht"])
    ntuple = config.get_ntuple('single_trigger', 'measurement')
    filename = ntuple.get_filename(year=year, workdir=input_dir)

    chans = ['Electron','Muon']


    plotter = Plotter(filename, groups, ntuple=ntuple, year=year)

    plotter.set_groups(bkg='wjet_ht')

    trig_cuts = {"Muon": 20, "Electron": 25}

    for reg in ["low", "high"]:
        for chan in chans:
            latex = latex_chan[chan]
            chan_name = name_chan[chan]

            plotter.mask(lambda vg : vg[chan].num() == 1)
            plotter.fill_hists(graphs, ginfo)
            pt_all = plotter.get_hists('pt')
            print(chan, "all")

            plotter.mask(lambda vg : vg[f"hlt_{reg}_{chan_name}"] > 0, clear=False)
            plotter.scale(lambda vg : scale_prescale(vg, f'hlt_{reg}_{chan_name}'), groups='data')

            plotter.fill_hists(graphs, ginfo)
            pt_hlt = plotter.get_hists('pt')
            print(chan, reg)

            eff_data = Histogram.efficiency(pt_hlt['data'], pt_all['data'])
            eff_mc = Histogram.efficiency(pt_hlt['wjet_ht'], pt_all['wjet_ht'])
            print(np.array2string(eff_data.vals, separator=','))
            print(np.array2string(eff_mc.vals, separator=','))
            print(pt_all['data'].vals)
            print(pt_hlt['data'].vals)


            plotter.scale(lambda vg : scale_prescale(vg, f'hlt_{reg}_{chan_name}', True), groups='data')

if __name__ == "__main__":
    workdir = workspace_area/'trigger_eff'

    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-d', '--workdir', help="directory to run over. If nothing, use date",)
    parser.add_argument('-i', '--input_dir', default=None)
    parser.add_argument('-r', '--run', type=lambda x: [i.strip() for i in x.split(',')],
                        help="Regions to run through (sideband, measurement, closure, dy)")
    args = parser.parse_args()

    workdir /= datetime.now().strftime("%m%d") if args.workdir is None else args.workdir
    # workdir.mkdir(exist_ok=True)

    color_by_group = {
        "data": "black",
        "wjet_ht": "olive",
        "wjets": "olive",
    }
    ginfo = GroupInfo(color_by_group)

    for year in args.years:
        trigger_eff(workdir, ginfo, year, args.input_dir)
