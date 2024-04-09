#!/usr/bin/env python3
import numpy as np
import uproot

@np.vectorize
def likelihood_sig(s, b):
    return np.sqrt(2*(s+b)*np.log(1+s/(b+1e-5))-2*s)

def get_syst_index(filename, group, systName):
    with uproot.open(filename) as f:
        print(group)
        if group not in f:
            return None, None
        tree = f[group]
        systNames = [name.member("fName") for name in tree["Systematics"]]
        print(np.unique(systNames))
        if systName not in systNames:
            return None, None
        systNum = systNames.index(systName)

        if "Syst_Index" in tree:
            for name in tree["Syst_Index"]:
                if int(name.member("fName")) == systNum:
                    return int(name.member('fTitle')), systNum
            return None, None
        else:
            return systNum, systNum
