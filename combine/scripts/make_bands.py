#!/usr/bin/env python3
import argparse
import uproot
import numpy as np

import analysis_suite.commons.user as user
from analysis_suite.commons.constants import all_eras
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import ratio_plot
from analysis_suite.commons.configs import get_inputs
from analysis_suite.data.systs import systematics, get_shape_systs, dummy

def plot_updown(pad, hist, **kwargs):
    pad.hist(x=hist.axis.centers, weights=hist.vals, bins=hist.axis.edges,
             label=hist.name, histtype="step", linewidth=2,
             color=hist.color, edgecolor=hist.darkenColor(), **kwargs)

signal = ['TTTJ', 'TTTW']

def make_band_group(f, workdir, group, syst, nom_hists, outdir, region):
    # Make plot of nominal and up and down variation per group
    hist = nom_hists[group]
    axis = hist.axis
    up_color = 'orangered'
    down_color = 'dodgerblue'
    nom = Histogram("", axis, name="Nominal", color='black')
    up = Histogram("", axis, color=up_color)
    down = Histogram("", axis, color=down_color)

    if group in f[f'{syst}Up']:
        up += f[f'{syst}Up'][group].to_boost()
        down += f[f'{syst}Down'][group].to_boost()
    else:
        up += hist
        down += hist
    nom += hist

    up_ratio = up/nom
    down_ratio = down/nom
    up_ratio.color = up_color
    down_ratio.color = down_color

    max_diff = np.max([up_ratio.vals, down_ratio.vals])-1
    min_diff = 1-np.min([up_ratio.vals, down_ratio.vals])
    padding = 0.2
    diff = np.max([max_diff, min_diff])/(1-padding)

    graph_info = get_inputs(workdir, 'combine_info').regions[region]['graph']
    with ratio_plot(outdir/f"{region}_{syst}_{group}.png", graph_info.axis_name,
                    nom.get_xrange(), ratio_top=1+diff, ratio_bot=1-diff) as ax:
        pad, subpad = ax
        nom.plot_points(pad)
        plot_updown(pad, up)
        plot_updown(pad, down)

        subpad.plot([axis.edges[0], axis.edges[-1]], [1, 1], color='k', linewidth=2)
        plot_updown(subpad, up_ratio)
        plot_updown(subpad, down_ratio)


def make_band(f, workdir, syst, nom_hists, outdir, region, bkg_or_sig='bkg'):
    axis = list(nom_hists.values())[0].axis
    up_color = 'orangered'
    down_color = 'dodgerblue'

    nom = Histogram("", axis, name="Background", color='black')
    up = Histogram("", axis, name=f"{syst} Up", color=up_color)
    down = Histogram("", axis, name=f"{syst} Down", color=down_color)
    up_nosmooth = Histogram("", axis, name="Up", color=up_color)
    down_nosmooth = Histogram("", axis, name="Down", color=down_color)

    isJEC = "JER" in syst or "JEC" in syst
    if isJEC:
        start = syst.find("LOWESS")
        end = start + 6
        if syst[start-4:start-2] == '20':
            start = start - 4
        nosmooth_name = syst[:start]+syst[end:]

    for group, hist in nom_hists.items():
        if (group in signal) == (bkg_or_sig == 'sig'):
            nom += hist
            if group in f[f'{syst}Up']:
                up += f[f'{syst}Up'][group].to_boost()
                down += f[f'{syst}Down'][group].to_boost()
            else:
                up += hist
                down += hist
            if isJEC:
                if group in f[f'{nosmooth_name}Up']:                                        
                    up_nosmooth += f[f'{nosmooth_name}Up'][group].to_boost()
                    down_nosmooth += f[f'{nosmooth_name}Down'][group].to_boost()
                else:
                    up_nosmooth += hist
                    down_nosmooth += hist

    up_ratio = up/nom
    down_ratio = down/nom
    up_ratio.color = up_color
    down_ratio.color = down_color

    padding = 0.2
    diff = np.max([abs(1-up_ratio.vals), abs(1-down_ratio.vals)])/(1-padding)
    vals = np.array([0.5, 0.3, 0.18, 0.1,
                     0.05, 0.03, 0.018, 0.01,
                     0.005, 0.003, 0.0018, 0.001,
                     0.0005, 0.0003, 0.00018, 0.00001,
                     0.00005, 0.00003, 0.000018, 0.00001])
    # vals = np.array([0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001])
    diff2 = vals[np.argmin(np.where(diff>vals, 10, vals-diff))]
    # diff_map = {
    #     "JER2016LOWESS16APV": 0.018,
    #     "JECABSOLUTE2016LOWESS16APV": 0.03,
    #     "JECABSOLUTELOWESS": 0.03,
    #     "JECBBEC12016LOWESS16APV": 0.03,
    #     "JECBBEC1LOWESS": 0.018,
    #     "JECEC22016LOWESS16APV": 0.01,
    #     "JECEC2LOWESS": 0.01,
    #     "JECHF2016LOWESS16APV": 0.01,
    #     "JECHFLOWESS": 0.01,
    #     "JECRELATIVEBALLOWESS": 0.018,
    #     "JECRELATIVESAMPLE2016LOWESS16APV": 0.018,
    #     "JECFLAVORQCDLOWESS": 0.05,
    # }
    # if syst in diff_map:
    #     diff2 = diff_map[syst]
    print(syst, diff2, 1-diff2, 1+diff2)

    graph_info = get_inputs(workdir, 'combine_info').regions[region]['graph']
    with ratio_plot(outdir/f"{region}_{syst}.png", graph_info.axis_name, nom.get_xrange(),
                    ratio_top=1+diff2, ratio_bot=1-diff2) as ax:
        pad, subpad = ax
        nom.plot_points(pad)
        plot_updown(pad, up)
        plot_updown(pad, down)

        subpad.plot([axis.edges[0], axis.edges[-1]], [1, 1], color='k', linewidth=2)
        plot_updown(subpad, up_ratio)
        plot_updown(subpad, down_ratio)

        if isJEC:
            pad.plot([0], [0], linestyle='dotted', color=up_color, label=f'{nosmooth_name} Up', marker=',')
            pad.plot([0], [0], linestyle='dotted', color=down_color, label=f'{nosmooth_name} Down', marker=',')
            plot_updown(subpad, up_nosmooth/nom, linestyle='dotted')
            plot_updown(subpad, down_nosmooth/nom, linestyle='dotted')


def make_all_bands(workdir, combine_dir, year, ntupleName, use_regions):
    outdir = combine_dir/f'band_split_{year}'
    outdir.mkdir(exist_ok=True)

    shape_systs = [syst.get_name(year) for syst in systematics if syst.syst_type == "shape"]

    for fname in combine_dir.glob(f'*{year}*root'):
        region = fname.stem[fname.stem.rfind('_')+1:]
        if region == "card" or (use_regions and region not in use_regions):
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
                elif "Up" in name and name[:-4] in shape_systs:
                    systs.append(name[:-4])

            for syst in systs:
                for name in nominal.keys():
                    make_band_group(f, workdir, name, syst, nominal, outdir, region)
                # make_band(f, workdir, syst, nominal, outdir, region)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="main", description="")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area / x,
                        help="Working Directory")
    parser.add_argument('-t', '--extra_text', default="")
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016pre", "2016post", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument("-r", "--region", type=lambda x : x.split(','), default=None)
    parser.add_argument("--debug", action='store_true')
    args = parser.parse_args()

    combine_dir = args.workdir/'combine'/args.extra_text

    for year in args.years:
        make_all_bands(args.workdir, combine_dir, year, "signal", args.region)
