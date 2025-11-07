#!/usr/bin/env python3
import numpy as np
import pickle

def scale(vg, workdir, group, member, year, syst):
    with open(workdir/f"btag_scales.pkl", "rb") as f:
        scales = pickle.load(f)[year]
    # jet_bin = np.digitize(vg["NJets"], np.arange(1, 10)) - 1
    if group not in scales or not scales[group]:
        return
    if "BJet" not in syst and "Jet_JE" not in syst:
        syst = "Nominal"
    vg.scale = scales[group][syst]
