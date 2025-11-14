#!/usr/bin/env python3
import argparse
import boost_histogram.axis as axis
import numpy as np
import warnings
import pickle

warnings.simplefilter("ignore", UserWarning)

from analysis_suite.commons.user import workspace_area
from analysis_suite.commons.constants import all_eras
from analysis_suite.commons.configs import  get_ntuple
from analysis_suite.plotting.hist_getter import GraphInfo, HistGetter

def process(workdir, year):
    ntuple = get_ntuple('btag')
    ginfo = ntuple.get_info()
    graph = {
        'nloose': GraphInfo(r"$N_{{loose}}$", axis.Regular(1, 0, 3), lambda vg: (vg['N_loose_el']+vg['N_loose_mu'], vg.scale)),
    }

    hist_factory = HistGetter(ntuple, year)

    syst_tots = {}
    for syst in hist_factory.systs:
        if syst != 'Nominal' and "LHE_" not in syst and "PS_" not in syst:
            continue
        hist_factory.reset_syst(syst)
        hists = hist_factory.get_hists(graph, fix_negative=True)
        syst_tots[syst] = {}
        for name, hist in hists.items():
            for group, h in hist.items():
                syst_tots[syst][group] = h.vals[0]


    nom = syst_tots.pop("Nominal")
    data = {}
    for syst, groups in syst_tots.items():
        data[syst] = {}
        for group, tot in groups.items():
            data[syst][group] = nom[group]/tot

    # Dump MC scale factors
    scale_file = workdir/f"theory_scales.pkl"
    if scale_file.exists():
        with open(scale_file, "rb") as f:
            scales = pickle.load(f)
    else:
        scales = dict()
    scales[year] = data
    with open(scale_file, "wb") as f:
        pickle.dump(scales, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : all_eras if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-d', '--workdir', required=True,
                        help="directory to run over. If nothing, use date",)
    args = parser.parse_args()

    workdir = workspace_area / args.workdir
    workdir.mkdir(exist_ok=True, parents=True)

    for year in args.years:
        print(year)
        process(workdir, year)
