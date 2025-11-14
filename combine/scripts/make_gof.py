#!/usr/bin/env python3
import subprocess
import argparse
from pathlib import Path
import random
import json
import numpy as np
from contextlib import contextmanager

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import mplhep as hep

plt.style.use([hep.style.CMS])

# For certain CMSSW releases, you'll get a numpy error about subnormals. Shouldn't
# affect us, so just ignore
import warnings
warnings.filterwarnings('ignore')

lumi_dict = {
    "2016pre": 19.52,
    "2016post": 16.81,
    "2016_small": 12.9,
    "2016" : 35.9,
    "2017_e": 36.75,
    "2017" : 41.48,
    "2018" : 59.83,
    "all" : 139.0,
}

# Nice plotting wrapper
@contextmanager
def plot(filename, subplot_info=None, **kwargs):
    if subplot_info is None:
        subplot_info = {}
    fig, ax = plt.subplots(**subplot_info)
    yield ax
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()
    if hasattr(plot, "workdir"):
        filename = f"{plot.workdir}/{filename}"
    fig.savefig(filename, bbox_inches="tight", dpi=300)
    if "extra_format" in kwargs:
        fig.savefig(filename.with_suffix(f'.{kwargs["extra_format"]}'), bbox_inches="tight", dpi=300)
    plt.close(fig)

# Setup mplhep labeling
def cms_label(ax, lumi=None, year=None, hasData=False, label='Internal'):
    if lumi is None and year is None:
        hep.cms.label(ax=ax, label=label, data=hasData)
        return
    if year is not None:
        lumi = lumi_dict[year]
    lumi = round(lumi, 1)
    hep.cms.label(ax=ax, lumi=lumi, label=label, data=hasData)


def produce_gof(output, name, n_runs, n_toys, verbose=False):
    obs_made = Path(f'higgsCombine.{name}_obs.GoodnessOfFit.mH120.123456.root').exists()
    higgs_files = Path('.').glob(f"higgsCombine.{name}*")
    runs = n_runs - len(list(higgs_files)) + obs_made
    option = ""
    # option = '--setParameters mask_Name1_yr2016post_Multi=1'
    # option += ',mask_Name1_yr2016post_ttttCR=1'
    # option += ',mask_Name1_yr2016post_ttzCR=0'
    # option += ',mask_Name1_yr2016post_Dilepton=1'
    # option += ',mask_Name1_yr2016pre_Multi=1'
    # option += ',mask_Name1_yr2016pre_ttttCR=1'
    # option += ',mask_Name1_yr2016pre_ttzCR=0'
    # option += ',mask_Name1_yr2016pre_Dilepton=1'
    # option += ',mask_Name2_Multi=1'
    # option += ',mask_Name2_ttttCR=1'
    # option += ',mask_Name2_ttzCR=0'
    # option += ',mask_Name2_Dilepton=1'
    # option += ',mask_Name3_Multi=1'
    # option += ',mask_Name3_ttttCR=1'
    # option += ',mask_Name3_ttzCR=0'
    # option += ',mask_Name3_Dilepton=1'
    cmd = []
    for i in range(runs):
        seed = random.randint(0, 10000)
        cmd.append(f'combineTool.py -M GoodnessOfFit {ws} -n .{name}_{seed} -s {seed} ' +
                        f'--algo=saturated --rMin -100 --rMax 100 -t {n_toys} --saveToys {option}')
    if not obs_made:
        cmd.append(f'combineTool.py -M GoodnessOfFit {ws} -n .{name}_obs -D data_obs ' +
                    f'--algo=saturated --rMin -100 --rMax 100 --saveToys {option}')

    stdout = subprocess.STDOUT if verbose else subprocess.DEVNULL
    processes = [subprocess.Popen(i, shell=True, stdout=None, stderr=stdout) for i in cmd]
    for p in processes: p.wait()

    subprocess.call(f'combineTool.py -M CollectGoodnessOfFit -o gof_{name}.json --input higgsCombine.{name}* {option}',
                    shell=True, stdout=stdout, stderr=stdout)


def plot_gof(output, name, filetype='pdf'):
    out_file = output/f'gof_{name}'

    year_list = ['2016pre', '2016post', '2017', '2018', '2016', 'all']
    year = 'all'
    for yr in year_list:
        if yr in name:
            year = yr
            break

    nbins = args.nbins
    with open(f'gof_{name}.json') as f:
        data = json.load(f)['120.0']
        pvalue = data['p']
        obs = data['obs'][0]
        toys = data['toy']

        upper, lower = np.percentile(toys, [75, 25])
        iqr = upper - lower
        bin_min = int(lower - 2*iqr - 1)
        if bin_min < 0:
            bin_min = 0
        bin_max = int(upper + 2*iqr + 1)
        if obs > bin_max:
            bin_max = max(toys)

        split_bin = int(nbins*(obs-bin_min)/(bin_max-bin_min))
        width = (bin_max-obs)/(nbins-split_bin)
        bin_min = obs - width*split_bin

    bin_down = np.linspace(bin_min, obs, split_bin+1)
    bin_up = np.linspace(obs, bin_max, nbins-split_bin+1)

    with plot(out_file.with_suffix('.'+filetype)) as pad:
        style = {'linewidth': 1, 'alpha': 0.5, 'edgecolor':'k', 'histtype':'stepfilled'}
        pad.hist(toys, bin_down, **style)
        pad.hist(toys, bin_up, **style)
        pad.plot([obs, obs], [0, pad.get_ylim()[1]], linewidth=3, color='k')
        pad.text(obs, pad.get_ylim()[1], "Observed", va='bottom', ha='center')
        pad.text(0.05, 0.96, f'{len(toys)} Toys', transform=pad.transAxes, ha='left', va='top')
        pad.text(0.05, 0.91, f'p-value = {pvalue:0.3f}', transform=pad.transAxes, ha='left', va='top')
        pad.set_xlabel(r'$-2ln\lambda_{saturated}$')
        pad.set_ylabel(r'A.U.')
        pad.set_ylim(top=pad.get_ylim()[1]*1.3)
        cms_label(pad, year=year, hasData=True, label="Internal")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--toys', type=int, default=150)
    parser.add_argument('-n', '--n_runs', type=int, default=5)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-ft', '--filetype', default='png')
    parser.add_argument('--nbins', type=int, default=50)
    args = parser.parse_args()

    output = Path('output')
    output.mkdir(exist_ok=True)

    # Can probably change this later, but works for my naming convention
    for ws in Path('.').glob("*card.root"):
        name = str(ws)[:str(ws).find("_card")]
        print(name)

        produce_gof(output, name, args.n_runs, args.toys, verbose=args.verbose)
        plot_gof(output, name, filetype=args.filetype)
