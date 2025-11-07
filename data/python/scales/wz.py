#!/usr/bin/env python3
import numpy as np
import pickle

def scale(vg, workdir, group, member, year, syst):
    if member != "wzTo3lnu":
        return
    with open(workdir/'wz_scale_factor.pickle', 'rb') as f:
        scales = pickle.load(f)[year].vals
    njets = vg['NJets'] if "Flat" in repr(vg) else vg["Jets"].num()
    jet_bin = np.digitize(njets, np.arange(2,6)) - 1
    vg.scale = scales[jet_bin]
