#!/usr/bin/env python3
import argparse
import uproot
import numpy as np

import analysis_suite.commons.user as user
from analysis_suite.commons.constants import all_eras
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import ratio_plot

def plot_updown(pad, hist):
    pad.hist(x=hist.axis.centers, weights=hist.vals, bins=hist.axis.edges,
             label=hist.name, histtype="step", linewidth=2,
             color=hist.color, edgecolor=hist.darkenColor())

def make_band(f, group, syst, nom_hist, outfile):
    if group not in f[f'{syst}Up'] or group not in f[f'{syst}Down']:
        return
    up = Histogram("", nom_hist.axis, name="Up", color='orangered')
    down = Histogram("", nom_hist.axis, name="Down", color='dodgerblue')
    up += f[f'{syst}Up'][group].to_boost()
    down += f[f'{syst}Down'][group].to_boost()

    up_ratio = Histogram("", nom_hist.axis, color='orangered')
    down_ratio = Histogram("", nom_hist.axis, color='dodgerblue')
    up_ratio += nom_hist/up
    down_ratio += nom_hist/down

    max_val = np.max([up_ratio.vals, down_ratio.vals])-1
    min_val = 1-np.min([up_ratio.vals, down_ratio.vals])
    top = 1+10**np.round(np.log10(max_val))
    bot = 1-10**np.round(np.log10(min_val))

    with ratio_plot(outfile, "BDT", nom_hist.get_xrange(), ratio_top=top, ratio_bot=bot) as ax:
        pad, subpad = ax
        nom_hist.plot_points(pad)
        plot_updown(pad, up)
        plot_updown(pad, down)

        subpad.plot(nom_hist.get_xrange(), [1, 1], color='k', linewidth=2)
        plot_updown(subpad, up_ratio)
        plot_updown(subpad, down_ratio)


def make_all_bands(workdir, year, ntupleName):
    systs = list()
    nominal = dict()
    outdir = workdir/f'band_{year}'
    outdir.mkdir(exist_ok=True)

    for fname in workdir.glob(f'*{year}*root'):
        region = fname.stem[fname.stem.rfind('_')+1:]
        with uproot.open(fname) as f:
            for name, cls in f.classnames().items():
                if "/" in name:
                    continue
                elif "TH1" in cls and "data" not in name:
                    hist = f[name].to_boost()
                    nominal[name[:-2]] = Histogram("", hist.axes[0], name="Nominal")
                    nominal[name[:-2]] += hist
                elif "Up" in name:
                    # if "Muon_tthMVA" not in name: continue
                    systs.append(name[:-4])
            for syst in systs:
                for group, nom_hist in nominal.items():
                    make_band(f, group, syst, nom_hist, outdir/f"{region}_{syst}_{group}.png")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="main", description="")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area / x,
                        help="Working Directory")
    parser.add_argument('-t', '--extra_text', default="")
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016pre", "2016post", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument("--debug", action='store_true')
    args = parser.parse_args()

    combine_dir = args.workdir/'combine'/args.extra_text

    for year in args.years:
        make_all_bands(combine_dir, year, "signal")
