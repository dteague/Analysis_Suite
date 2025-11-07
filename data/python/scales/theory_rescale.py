#!/usr/bin/env python3
import numpy as np
import pickle

def scale(vg, workdir, group, member, year, syst):
    with open(workdir/f"theory_scales.pkl", "rb") as f:
        scales = pickle.load(f)[year]

    if syst not in scales or group not in scales[syst]:
        return
    vg.scale = scales[syst][group]
