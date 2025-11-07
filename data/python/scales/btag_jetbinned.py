#!/usr/bin/env python3
import numpy as np
import pickle

def scale(vg, workdir, group, member, year, syst):
    with open(workdir/f"btag_scales_njet.pkl", "rb") as f:
        scales = pickle.load(f)[year]
    if group not in scales or not scales[group]:
        return
    if "BJet" not in syst and "Jet_JE" not in syst:
        syst = "Nominal"
    scale = scales[group][syst].vals
    njets = vg['NJets'] if "Flat" in repr(vg) else vg["Jets"].num()
    jet_bin = np.digitize(njets, np.arange(2,7)) - 1
    vg.scale = scale[jet_bin]
