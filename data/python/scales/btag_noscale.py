#!/usr/bin/env python3
import numpy as np
import awkward as ak

def scale(vg, workdir, group, member, year, syst):
    vg['wgt_nobtag']
    print(sum(vg['NBjets_medium']*vg.scale)/sum(vg.scale))
    vg._scale = vg.get_sf(vg.syst_name)*ak.to_numpy(vg.arr['wgt_nobtag'])
    print(sum(vg['NBjets_medium']*vg.scale)/sum(vg.scale))
