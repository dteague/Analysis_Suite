#!/usr/bin/env python3
import numpy as np
import uproot

@np.vectorize
def likelihood_sig(s, b):
    return np.sqrt(2*(s+b)*np.log(1+s/(b+1e-5))-2*s)

def get_syst_index(filename, systName):
    with uproot.open(filename) as f:
        syst_dir = [key for key in f.keys() if "Systematics" in key][0]
        systNames = [name.member("fName") for name in f[syst_dir]]
    if systName not in systNames:
        return -1
    else:
        return systNames.index(systName)

