#!/usr/bin/env python3
import argparse
import uproot
import numpy as np
import multiprocessing as mp

import analysis_suite.commons.user as user
from analysis_suite.commons.constants import all_eras, lumi
from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.plot_utils import ratio_plot
from analysis_suite.commons.configs import get_inputs
from analysis_suite.data.systs import systematics, get_shape_systs, dummy
from analysis_suite.combine.systematics import use_lowess

def plot_updown(pad, hist, **kwargs):
    pad.hist(x=hist.axis.centers, weights=hist.vals, bins=hist.axis.edges,
             label=hist.plot_label, histtype="step", linewidth=2,
             color=hist.color, edgecolor=hist.darkenColor(), **kwargs)

signal = ['TTTJ', 'TTTW']
ratio_range = np.array([0.5, 0.3, 0.18, 0.1,
                    0.05, 0.03, 0.018, 0.01,
                    0.005, 0.003, 0.0018, 0.001,
                    0.0005, 0.0003, 0.00018, 0.00001,
                    0.00005, 0.00003, 0.000018, 0.00001])

def make_band_group(f, workdir, group, syst, nom_hists, outdir, region, year):
    # Make plot of nominal and up and down variation per group
    hist = nom_hists[group]
    axis = hist.axis
    up_color = 'orange'
    down_color = 'blue'
    nom = Histogram(axis, name="Nominal", color='black')
    up = Histogram(axis, color=up_color, name=f'{syst}_Up')
    down = Histogram(axis, color=down_color, name=f'{syst}_Down')
    up_orig = Histogram(axis, color=up_color)
    down_orig = Histogram(axis, color=down_color)

    if group not in f[f'{syst}Up']:
        return
    isJEC = 'LOWESS' in syst
    if isJEC:
        start = syst.find("LOWESS")
        end = start + 6
        if syst[start-4:start-2] == '20':
            start = start - 4
        nosmooth_name = syst[:start]+syst[end:]
        up_orig.plot_label = nosmooth_name+"_Up"
        down_orig.plot_label = nosmooth_name+"_Down"

    if group in f[f'{syst}Up']:
        up += f[f'{syst}Up'][group].to_boost()
        down += f[f'{syst}Down'][group].to_boost()
    if isJEC and group in f[f'{nosmooth_name}Up']:
        up_orig += f[f'{nosmooth_name}Up'][group].to_boost()
        down_orig += f[f'{nosmooth_name}Down'][group].to_boost()

    nom += hist

    up_ratio = up/nom
    up_ratio.values()[up_ratio.values() < 1e-6] = 1.
    down_ratio = down/nom
    down_ratio.values()[down_ratio.values() < 1e-6] = 1.
    up_ratio.color = up_color
    down_ratio.color = down_color

    if up_orig:
        up_orig_ratio = up_orig/nom
        up_orig_ratio.values()[up_orig_ratio.values() < 1e-6] = 1.
        down_orig_ratio = down_orig/nom
        down_orig_ratio.values()[down_orig_ratio.values() < 1e-6] = 1.
        up_orig_ratio.color = up_color
        down_orig_ratio.color = down_color

    if isJEC:
        all_vals = np.array([up_ratio.vals, down_ratio.vals,
                             up_orig_ratio.vals, down_orig_ratio.vals])
    else:
        all_vals = np.array([up_ratio.vals, down_ratio.vals])

    # diff = np.max(np.abs(all_vals-1))
    diff = np.max(np.where(all_vals > 1e-5, np.abs(all_vals-1), 0. ))
    ratio_bins = ratio_range[diff < ratio_range]
    ratio = ratio_range[-1] if len(ratio_bins) == 0 else ratio_bins[-1]

    graph_info = get_inputs(workdir, 'combine_info').regions[region]['graph']
    with ratio_plot(outdir/f"{region}_{syst}_{group}.pdf", graph_info.axis_name,
                    nom.get_xrange(), ratio_top=1+ratio, ratio_bot=1-ratio, lumi=lumi[year]) as ax:
        pad, subpad = ax
        nom.plot_points(pad)
        up.plot_hist(pad)
        down.plot_hist(pad)

        up_ratio.plot_hist(subpad)
        down_ratio.plot_hist(subpad)
        if up_orig:
            alpha = 0.3
            up_orig.plot_hist(pad, linestyle='--', alpha=alpha)
            down_orig.plot_hist(pad, linestyle='--', alpha=alpha)
            up_orig_ratio.plot_hist(subpad, linestyle='--', alpha=alpha)
            down_orig_ratio.plot_hist(subpad, linestyle='--', alpha=alpha)


def make_band(f, workdir, syst, nom_hists, outdir, region, year, bkg_or_sig='bkg'):
    axis = list(nom_hists.values())[0].axis
    up_color = 'orange'
    down_color = 'blue'

    nom = Histogram(axis, name="Background", color='black')
    up = Histogram(axis, color=up_color, name=f'{syst}_Up')
    down = Histogram(axis, color=down_color, name=f'{syst}_Down')
    up_orig = Histogram(axis, color=up_color)
    down_orig = Histogram(axis, color=down_color)

    isJEC = 'LOWESS' in syst or "JER20" in syst
    if isJEC and "JER" not in syst:
        start = syst.find("LOWESS")
        end = start + 6
        if syst[start-4:start-2] == '20':
            start = start - 4
        nosmooth_name = syst[:start]+syst[end:]
        up_orig.plot_label = nosmooth_name+"_Up"
        down_orig.plot_label = nosmooth_name+"_Down"
    if "JER20" in syst:
        nosmooth_name = f"JER{syst[-2:]}"
        up_orig.plot_label = nosmooth_name+"_Up"
        down_orig.plot_label = nosmooth_name+"_Down"
        up.plot_label = "JERLOWESS_up"
        down.plot_label = "JERLOWESS_up"

    for group, hist in nom_hists.items():
        if (group == 'SIG') != (bkg_or_sig == 'sig'):
            continue
        nom += hist
        if group in f[f'{syst}Up']:
            up += f[f'{syst}Up'][group].to_boost()
            down += f[f'{syst}Down'][group].to_boost()
        else:
            up += hist
            down += hist
        if isJEC and group in f[f'{nosmooth_name}Up']:
            up_orig += f[f'{nosmooth_name}Up'][group].to_boost()
            down_orig += f[f'{nosmooth_name}Down'][group].to_boost()
        elif isJEC:
            up_orig += hist
            down_orig += hist

    if up_orig:
        up_orig_ratio = up_orig/nom
        up_orig_ratio.values()[up_orig_ratio.values() < 1e-6] = 1.
        down_orig_ratio = down_orig/nom
        down_orig_ratio.values()[down_orig_ratio.values() < 1e-6] = 1.
        up_orig_ratio.color = up_color
        down_orig_ratio.color = down_color
    up_ratio = up/nom
    up_ratio.values()[up_ratio.values() < 1e-6] = 1.
    down_ratio = down/nom
    down_ratio.values()[down_ratio.values() < 1e-6] = 1.
    up_ratio.color = up_color
    down_ratio.color = down_color

    if isJEC:
        all_vals = np.array([up_ratio.vals, down_ratio.vals,
                             up_orig_ratio.vals, down_orig_ratio.vals])
    else:
        all_vals = np.array([up_ratio.vals, down_ratio.vals])

    if "JEC" in syst and "LOWESS" not in syst:
        return
    if "JER20" not in syst:
        return
    # print(list(up_ratio.vals))
    # print(list(nom.vals))
    # print(syst, np.sqrt(sum((1-up_ratio.vals)**2*nom.vals)/nom.integral()), max(abs(1-up_ratio.vals)))
    # print(syst, up/nom.integral()-down.integral()/nom.integral()))
    # return
    # diff = np.max(np.abs(all_vals-1))
    diff = np.max(np.where(all_vals > 1e-5, np.abs(all_vals-1), 0. ))
    ratio_bins = ratio_range[diff < ratio_range]
    ratio = ratio_range[-1] if len(ratio_bins) == 0 else ratio_bins[-1]

    graph_info = get_inputs(workdir, 'combine_info').regions[region]['graph']
    with ratio_plot(outdir/f"{region}_{syst}.pdf", graph_info.axis_name, nom.get_xrange(),
                    ratio_top=1+ratio, ratio_bot=1-ratio, lumi=lumi[year]) as ax:
        pad, subpad = ax
        nom.plot_points(pad)

        up.plot_hist(pad)
        down.plot_hist(pad)

        up_ratio.plot_hist(subpad)
        down_ratio.plot_hist(subpad)
        if up_orig:
            alpha = 0.3
            up_orig.plot_hist(pad, linestyle='--', alpha=alpha)
            down_orig.plot_hist(pad, linestyle='--', alpha=alpha)
            up_orig_ratio.plot_hist(subpad, linestyle='--', alpha=alpha)
            down_orig_ratio.plot_hist(subpad, linestyle='--', alpha=alpha)


def star_make_bands(filename, workdir, syst, nominal, outdir, region, year):
    with uproot.open(filename) as f:
        members = list(f.keys(recursive=False, cycle=False, filter_classname='TH1*'))
        members.remove('data_obs')
        nominal = {m: Histogram(f[m].to_boost()) for m in members}

        # for name in nominal.keys():
        #     make_band_group(f, workdir, name, syst, nominal, outdir, region, year)
        make_band(f, workdir, syst, nominal, outdir, region, year)



def make_all_bands(workdir, combine_dir, year, ntupleName, region, cores=1):
    outdir = combine_dir/f'band_lowess_{year}'
    outdir.mkdir(exist_ok=True)

    shape_systs = get_shape_systs(year)
    inputs = list()

    filename = combine_dir/f'combine_{year}_{region}.root'
    with uproot.open(filename) as f:
        members = list(f.keys(recursive=False, cycle=False, filter_classname='TH1*'))
        members.remove('data_obs')
        nominal = {m: Histogram(f[m].to_boost()) for m in members}
        systs = [s[:-2] for s in f.keys(recursive=False, cycle=False, filter_name="*Up", filter_classname='TDirectory')]

        lowess_systs = [s.replace('LOWESS', '').replace(year, '') for s in systs if 'LOWESS' in s]
        for syst in systs:
            # keeps = ['CFERR', 'HF', "ISR", "LF", 'MU']
            # keeps = ['JEC']
            # keep = False
            # for k in keeps:
            #     if k in syst:
            #         keep = True
            # if not keep:
            #     continue
            # if syst in lowess_systs :
            #     continue
            inputs.append([filename, workdir, syst, nominal, outdir, region, year])
            if cores != 1:
                continue

            for name in nominal.keys():
                make_band_group(f, workdir, name, syst, nominal, outdir, region, year)
            make_band(f, workdir, syst, nominal, outdir, region, year)

    if cores != 1:
        with mp.Pool(cores) as pool:
            pool.starmap(star_make_bands, inputs)


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
    parser.add_argument("-j", '--cores', default=1, type=int)
    parser.add_argument("--debug", action='store_true')
    args = parser.parse_args()

    combine_dir = args.workdir/'combine'/args.extra_text

    for year in args.years:
        for region in args.region:
            make_all_bands(args.workdir, combine_dir, year, "signal", region, cores=args.cores)
