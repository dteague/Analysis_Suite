#!/usr/bin/env python3
import gzip
import argparse
import numpy as np
from correctionlib.schemav2 import VERSION, MultiBinning, Category, Correction, CorrectionSet

import analysis_suite.commons.configs as config
from analysis_suite.commons.histogram import Histogram
import analysis_suite.data.plotInfo.befficiency as pinfo
from analysis_suite.plotting.plotter import Plotter
import analysis_suite.commons.user as user

def get_sf(sf, syst):
    val, err = sf.value, np.sqrt(sf.variance)
    if syst == 'central':
        return val
    else:
        return val + (err if syst=='up' else -err)

def build_pteta(sf, syst):
    npt, neta = sf.axes.size
    ptbins = sf.axes[0].edges.tolist()
    ptbins[-1] = "Infinity"

    return MultiBinning.parse_obj({
        "nodetype": "multibinning",
        "inputs": ["abseta", "pt"],
        "edges": [sf.axes[1].edges.tolist(), ptbins],
        "content": [get_sf(sf.hist[pt, eta], syst) for eta in range(neta) for pt in range(npt)],
        "flow": "error",
    })

def build_scales(sf):
    build_flavor = lambda subsf, syst : Category.parse_obj({
        "nodetype": "category",
        "input": "flavor",
        "content": [ {"key": key, "value": build_pteta(value, syst) } for key, value in subsf.items() ],
    })
    build_wp = lambda subsf, syst : Category.parse_obj({
        "nodetype": "category",
        "input": "working_point",
        "content": [ {"key": key, "value": build_flavor(value, syst) } for key, value in sf.items() ],
    })
    return Category.parse_obj({
        "nodetype": "category",
        "input": "systematic",
        "content": [{"key": syst, "value": build_wp(sf, syst) } for syst in systs],
    })

jet_flavs = {'udsg': 0, 'c': 4, 'b': 5}
wps = {"L": 1, "M": 2, "T": 3}
systs = ['central', 'down', 'up']

def make_efficiencies(year, input_dir):
    weights = {wp: dict() for wp in wps.keys()}

    ntuple = config.get_ntuple('befficiency', 'measurement')
    ginfo = ntuple.get_info()
    graphs = pinfo.ptetas

    plotter = Plotter(ntuple.get_filename(year=year, workdir=input_dir), ginfo.setup_groups(),
                      ntuple=ntuple, year=year)
    plotter.fill_hists(graphs)
    all_groups = ginfo.get_groups()

    for flav_name, flav in jet_flavs.items():
        all_hist = plotter.get_sum(all_groups, f'{flav_name}_all')
        for wp_name, wp_id in wps.items():
            wp_hist = plotter.get_sum(all_groups, f'{flav_name}_{wp_name}')
            weights[wp_name][flav] = Histogram.efficiency(wp_hist, all_hist)
            # print(wp_name, flav)
    #         print((1-weights[wp_name][flav].vals))
    #         print((weights[wp_name][flav].err)/(1-weights[wp_name][flav].vals))
    #         print()
    # return
    cset = CorrectionSet.parse_obj({
        "schema_version": VERSION,
        "corrections": [
            Correction.parse_obj({
                "version": 1,
                "name": "SS",
                "description": "BTagging efficiencies for different flavor jets in 3Top-2lepton-SS signal region",
                "inputs": [
                    {"name": "systematic", "type": "string",
                     "description": "Central value and shifts (statistical only) in efficiency"},
                    {"name": "working_point", "type": "string",
                     "description": "Working point to get efficiency"},
                    {"name": "flavor", "type": "int",
                     "description": "hadron flavor definition: 5=b, 4=c, 0=udsg"},
                    {"name": "abseta", "type": "real", "description": "Jet abseta"},
                    {"name": "pt", "type": "real", "description": "Jet pt"},
                ],
                "output": {"name": "weight", "type": "real", "description": "Efficiency of passing Btag cut"},
                "data": build_scales(weights),
            })
        ],
    })

    outdir = user.analysis_area / 'data/POG/USER/'
    if year == "2016":
        with gzip.open(outdir/"2016preVFP_UL"/"beff.json.gz", "wt") as fout:
            fout.write(cset.json(exclude_unset=True, indent=4))
        with gzip.open(outdir/"2016postVFP_UL"/"beff.json.gz", "wt") as fout:
            fout.write(cset.json(exclude_unset=True, indent=4))
    else:
        with gzip.open(outdir/f"{year}_UL"/"beff.json.gz", "wt") as fout:
            fout.write(cset.json(exclude_unset=True, indent=4))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-i', '--input_dir', default=None)
    args = parser.parse_args()

    for year in args.years:
        make_efficiencies(year, args.input_dir)
