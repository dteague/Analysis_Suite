#!/usr/bin/env python3
import argparse
import uproot
import numpy as np

import analysis_suite.commons.user as user
from analysis_suite.commons.constants import all_eras
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import ratio_plot

def plot_updown(pad, hist, **kwargs):
    pad.hist(x=hist.axis.centers, weights=hist.vals, bins=hist.axis.edges,
             label=hist.name, histtype="step", linewidth=2,
             color=hist.color, edgecolor=hist.darkenColor(), **kwargs)

signal = ['TTTJ', 'TTTW']

def make_band(f, syst, nom_hists, outfile):
    print(outfile.name)
    axis = list(nom_hists.values())[0].axis
    up_color = 'orangered'
    down_color = 'dodgerblue'
    sig = Histogram("", axis, name="Signal", color='black')
    sig_up = Histogram("", axis, color=up_color)
    sig_down = Histogram("", axis, color=down_color)
    bkg = Histogram("", axis, name="Background", color='black')
    bkg_up = Histogram("", axis, name="Up", color=up_color)
    bkg_down = Histogram("", axis, name="Down", color=down_color)

    for group, hist in nom_hists.items():
        if group in signal:
            sig += hist
            if group in f[f'{syst}Up']:
                sig_up += f[f'{syst}Up'][group].to_boost()
                sig_down += f[f'{syst}Down'][group].to_boost()
            else:
                sig_up += hist
                sig_down += hist
        else:
            bkg += hist
            if group in f[f'{syst}Up']:
                bkg_up += f[f'{syst}Up'][group].to_boost()
                bkg_down += f[f'{syst}Down'][group].to_boost()
            else:
                bkg_up += hist
                bkg_down += hist

    scale = 50*int(max(bkg.vals)/(50*max(sig.vals)))
    sig.scale(scale, changeName=True)
    sig_up.scale(scale)
    sig_down.scale(scale)

    sig_up_ratio = sig_up/sig
    sig_down_ratio = sig_down/sig
    bkg_up_ratio = bkg_up/bkg
    bkg_down_ratio = bkg_down/bkg

    sig_up_ratio.color = up_color
    sig_down_ratio.color = down_color
    bkg_up_ratio.color = up_color
    bkg_down_ratio.color = down_color


    max_val = np.max([sig_up_ratio.vals, sig_down_ratio.vals,
                      bkg_up_ratio.vals, bkg_down_ratio.vals])-1
    min_val = 1-np.min([sig_up_ratio.vals, sig_down_ratio.vals,
                      bkg_up_ratio.vals, bkg_down_ratio.vals])
    def rounder(val):
        order = np.floor(np.log10(val))
        rounded = np.round(2*val+10**order, -int(order))/2
        if rounded > 1:
            return 1.
        else:
            return rounded

    top = 1+rounder(max_val)
    bot = 1-rounder(min_val)

    with ratio_plot(outfile, "BDT", sig.get_xrange(), ratio_top=top, ratio_bot=bot) as ax:
        pad, subpad = ax
        sig.plot_points(pad)
        plot_updown(pad, sig_up, linestyle='dashdot')
        plot_updown(pad, sig_down, linestyle='dashdot')
        bkg.plot_points(pad)
        plot_updown(pad, bkg_up)
        plot_updown(pad, bkg_down)

        subpad.plot([axis.edges[0], axis.edges[-1]], [1, 1], color='k', linewidth=2)
        plot_updown(subpad, sig_up_ratio, linestyle='dashdot')
        plot_updown(subpad, sig_down_ratio, linestyle='dashdot')
        plot_updown(subpad, bkg_up_ratio)
        plot_updown(subpad, bkg_down_ratio)

def make_all_bands(workdir, year, ntupleName):
    outdir = workdir/f'band_{year}'
    outdir.mkdir(exist_ok=True)

    for fname in workdir.glob(f'*{year}*root'):
        region = fname.stem[fname.stem.rfind('_')+1:]
        if 'Dilep' not in region:
            continue
        nominal = dict()
        systs = list()
        with uproot.open(fname) as f:
            for name, cls in f.classnames().items():
                if "/" in name:
                    continue
                elif "TH1" in cls and "data" not in name:
                    hist = f[name].to_boost()
                    nominal[name[:-2]] = Histogram("", hist.axes[0], name="Nominal")
                    nominal[name[:-2]] += hist
                elif "Up" in name:
                    systs.append(name[:-4])
            for syst in systs:
                make_band(f, syst, nominal, outdir/f"{region}_{syst}.png")



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
