#!/usr/bin/env python3
from pathlib import Path
import awkward as ak
import boost_histogram as bh
import gzip
import argparse
import pickle
import numpy as np
from correctionlib.schemav2 import VERSION, MultiBinning, Category, Correction, CorrectionSet
import correctionlib._core as core
import analysis_suite.commons.user as user

def get_sf(sf, syst):
    val, err = sf.value, sf.variance
    if syst == "up":
        val += err
    elif syst == "up":
        val -= err
    return val/(1-val)

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

def build_corr(name, desc, sf):
    systs = ["nom", "up", "down"]
    return Correction.parse_obj({
        "version": 1,
        "name": name,
        "description": desc,
        "inputs": [
            {"name": "systematic", "type": "string",
             "description": "Central value and shifts (statistical only)"},
            {"name": "abseta", "type": "real", "description": "Jet abseta"},
            {"name": "pt", "type": "real", "description": "Jet pt"},
        ],
        "output": {"name": "weight", "type": "real", "description": "Fake Rate Scale Factor"},
        "data": Category.parse_obj({
            "nodetype": "category",
            "input": "systematic",
            "content": [{"key": syst, "value": build_pteta(sf, syst) } for syst in systs],
        })
    })



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('--nonprompt', required=True)
    parser.add_argument('--charge', required=True)
    args = parser.parse_args()

    basePath = user.analysis_area / "data/POG/USER"

    for year in args.years:
        with open(user.workspace_area/"charge_misId"/f"{args.charge}/charge_misid_rate_scaled_{year}.pickle", "rb") as f:
            charge_fr = pickle.load(f)
        with open(user.workspace_area/"fake_rate"/f"{args.nonprompt}/fr_{year}.pickle", "rb") as f:
            fake_rates = pickle.load(f)
            elec_nonprompt_fr = fake_rates["Electron"]["data_ewk"]
            muon_nonprompt_fr = fake_rates["Muon"]["data_ewk"]

        cset = CorrectionSet.parse_obj({
            "schema_version": VERSION,
            "corrections": [
                build_corr("Nonprompt_electron", "Transfer factor for Nonprompt bkg for electrons", elec_nonprompt_fr),
                build_corr("Nonprompt_muon", "Transfer factor for Nonprompt bkg for muons", muon_nonprompt_fr),
                build_corr("Charge_MisId", "Transfer factor for ChargeMisId in electrons", charge_fr),
            ],
        })

        def write_out(year):
            with gzip.open(basePath / f'{year}_UL' / "fake_rates.json.gz", "wt") as fout:
                fout.write(cset.json(exclude_unset=True, indent=4))
        
        if year == '2016':
            write_out('2016preVFP')
            write_out('2016postVFP')
        write_out(year)
