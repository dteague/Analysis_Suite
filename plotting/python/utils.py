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
        index_dir = [key for key in f.keys() if "Syst_Index" in key][0]
        syst_index = {systNames[int(name.member("fName"))]: int(name.member("fTitle")) for name in f[index_dir]}
    if systName not in syst_index:
        return None, None
    else:
        return syst_index[systName], systNames.index(systName)
