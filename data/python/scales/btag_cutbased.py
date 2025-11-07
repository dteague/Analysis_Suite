#!/usr/bin/env python3
import numpy as np
import pickle

def scale(vg, workdir, group, member, year, syst):
    vg['wgt_nobtag'], vg['bwgt_cb_m']
    vg._scale = vg.get_sf(vg.syst_name)*vg.arr['wgt_nobtag'] * vg.arr['bwgt_cb_m']
