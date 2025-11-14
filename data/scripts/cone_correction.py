#!/usr/bin/env python3
import uproot
import numpy as np
import awkward as ak
import argparse
from datetime import datetime

from analysis_suite.commons.constants import lumi
from analysis_suite.commons.info import fileInfo
import analysis_suite.commons.configs as config
from analysis_suite.commons.plot_utils import plot, cms_label
from analysis_suite.commons.user import workspace_area


latex_chan = {"Electron": "e", "Muon": "\mu",
              "EE": 'ee', "EM": 'e\mu', 'MM': '\mu\mu'}

def get_wgt(evt, chan):
    wgt_all = evt["weight"][:,0]
    tight = ak.num(evt[f'Tight{chan}/rawPt'], axis=-1)
    fake = ak.num(evt[f'Tight{chan}/rawPt'], axis=-1)
    return ak.concatenate((wgt_all[tight], wgt_all[fake]))

def get_var(evt, chan, var):
    return ak.concatenate((ak.flatten(evt[f'Tight{chan}/{var}']), ak.flatten(evt[f'Fake{chan}/{var}'])))

def add_hist(old_hist, new_hist):
    if old_hist is None:
        return new_hist
    for chan, hists in new_hist.items():
        for name, hist in hists.items():
            old_hist[chan][name] += hist
    return old_hist

def plot_ptcorr(filename, pt, bins, lumi, chan, fact=None):
    with plot(filename, extra_format='pdf') as ax:
        split = 2./3
        print_spot= split*max(pt) + (1-split)*min(pt)
        width = (bins[1]-bins[0])/2

        ax.hist(bins[:-1], bins, weights=pt, histtype='step')
        # ax.errorbar(bins[:-1]+width, pt, yerr=err, ls="", capsize=5, ecolor='blue')
        ax.set_xlim(-1, 1)
        ax.set_ylim(0.9*min(pt))
        ax.plot([0.65, 0.65], [0, max(pt)*1.5], color='red')

        fact_text = r"$\frac{p_{T}}{p_{T}^{ratio}}$" if fact is None else f"{fact:0.3f}"+r"$\times\frac{p_{T}}{p_{T}^{ratio}}$"
        ax.text(-0.35, print_spot, fact_text)
        ax.text(0.775, print_spot, r"$p_{T}$")
        ax.set_xlabel("$disc_{TTH}$")
        ax.set_ylabel(f"$mean(p_{{T}}({latex_chan[chan]}))$")
        cms_label(ax=ax, lumi=lumi)


parts = ['TightElectron', "FakeElectron"]

branches = ["TightElectron/new_mvaTTH", "TightElectron/rawPt", "TightElectron/ptRatio",
            "FakeElectron/new_mvaTTH", "FakeElectron/rawPt", "FakeElectron/ptRatio",
            "TightMuon/mvaTTH", "TightMuon/rawPt", "TightMuon/ptRatio",
            "FakeMuon/mvaTTH", "FakeMuon/rawPt", "FakeMuon/ptRatio",
            "weight"
            ]
mva_name = {"Electron": "new_mvaTTH", "Muon": "mvaTTH"}

nbins = 57
# mva_bins = np.linspace(-1, 1, nbins+1)
mva_bins = np.linspace(-0.995, 1, nbins+1)
chans = ["Muon", "Electron"]
width = (mva_bins[-1]-mva_bins[0])/nbins
flip_bin = round(0.35/width)


def get_hist(year, sample):
    hists = {
        "Electron" : {
            "raw": np.zeros(nbins),
            "wgt": np.zeros(nbins),
            "jet": np.zeros(nbins),
        },
        "Muon": {
            "raw": np.zeros(nbins),
            "wgt": np.zeros(nbins),
            "jet": np.zeros(nbins),
        }
    }

    ntuple = config.get_ntuple('fake_rate', 'pt_correction')
    filename = next(ntuple.get_filename(year=year).glob("*.root"))
    with uproot.open(filename) as f:
        for sample, sample_evt in f.items():
            sample_evt = f[sample]
            sumw = sample_evt["sumweight"].to_numpy()[0][0]
            xsec = fileInfo.get_xsec(sample)*lumi[year]*1000
            for batch in sample_evt["Measurement"].iterate(branches, step_size=1000000):
                for chan in chans:
                    wgt_all = get_wgt(batch, chan)
                    raw_pt = get_var(batch, chan, "rawPt")
                    ratio = get_var(batch, chan, "ptRatio")
                    jet_pt = raw_pt/ratio
                    mva = get_var(batch, chan, mva_name[chan])

                    for i in range(nbins):
                        mask = (mva >= mva_bins[i]) * (mva < mva_bins[i+1])
                        hists[chan]['wgt'][i] += ak.sum(wgt_all[mask])
                        hists[chan]['raw'][i] += ak.sum(raw_pt[mask]*wgt_all[mask])
                        hists[chan]['jet'][i] += ak.sum(jet_pt[mask]*wgt_all[mask])
    for chan in chans:
        for typ in hists[chan].keys():
            hists[chan][typ] *= xsec/sumw
    return hists


def process_year(year, workdir):
    ntuple = config.get_ntuple('fake_rate', 'pt_correction')
    all_hist = None
    for group, members in groups.items():
        for member in members:
            hist = get_hist(year, member)
            if hist is None:
                continue
            all_hist = add_hist(all_hist, hist)

    corr_fact = {
        "Electron": {
            '2016': 0.850,
            '2017': 0.800,
            '2018': 0.800
        },
        "Muon": {
            '2016': 0.775,
            '2017': 0.725,
            '2018': 0.750
        },
    }

    for chan in chans:
        raw = all_hist[chan]['raw']/all_hist[chan]['wgt']
        jet = all_hist[chan]['jet']/all_hist[chan]['wgt']
        tot = np.concatenate((jet[:-flip_bin], raw[-flip_bin:]))
        raw_corr = tot[-flip_bin]/tot[-flip_bin-1]
        print(chan)
        plot_ptcorr(workdir/f"pt_uncorr_{chan}_{year}.png", tot, mva_bins, lumi[year], chan)
        fact = corr_fact[chan][year]# corr_facts[np.argmin(abs(corr_facts - raw_corr))]
        tot = np.concatenate((fact*jet[:-flip_bin], raw[-flip_bin:]))
        plot_ptcorr(workdir/f"pt_corr_{chan}_{year}.png", tot, mva_bins, lumi[year], chan, fact=fact)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-d', '--workdir', default=datetime.now().strftime("%m%d"),
                        help="directory to run over. If nothing, use date",)
    parser.add_argument('-i', '--input_dir', default=None)
    args = parser.parse_args()

    workdir = workspace_area / args.workdir / "cone_correction"
    workdir.mkdir(exist_ok=True)

    for year in args.years:
        print(year)
        process_year(year, workdir)
