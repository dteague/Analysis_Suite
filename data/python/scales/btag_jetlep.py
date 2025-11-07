#!/usr/bin/env python3
import numpy as np
import pickle

def scale(vg, workdir, group, member, year, syst):
    with open(workdir/f"btag_scales_njet_lep2.pkl", "rb") as f:
        sfile_2 = pickle.load(f)[year]
    with open(workdir/f"btag_scales_njet_lep3.pkl", "rb") as f:
        sfile_3 = pickle.load(f)[year]
    if group not in sfile_2 or not sfile_2[group]:
        return
    if "BJet" not in syst and "Jet_JE" not in syst:
        syst = "Nominal"
    scale_2 = sfile_2[group][syst].vals
    scale_3 = sfile_3[group][syst].vals
    njets = vg['NJets'] if "Flat" in repr(vg) else vg["Jets"].num()
    # nleps = vg['NMuons']+vg['NElectrons'] if "Flat" in repr(vg) else vg["TightLepton"].num()
    if 'Flat' not in repr(vg):
       nleps = vg["TightLepton"].num()
    elif "NLeps" in vg.branches:
        nleps = vg['NLeps']
    else:
        nleps = vg['NMuons']+vg['NElectrons']
    jet_bin = np.digitize(njets, np.arange(2,7)) - 1

    if np.count_nonzero(nleps == 2):
        vg.scale = (scale_2[jet_bin[nleps==2]], nleps == 2)
    if np.count_nonzero(nleps > 2):
        vg.scale = (scale_3[jet_bin[nleps>2]], nleps > 2)
