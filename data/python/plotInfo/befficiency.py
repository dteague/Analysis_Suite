#!/usr/bin/env python3
import boost_histogram as bh
from analysis_suite.plotting.plotter import GraphInfo
import awkward as ak

def get_pteta(part, flav, wp):
    flav_mask = part.flavor == flav
    wp_mask = part.pass_btag[flav_mask] >= wp
    pt = ak.flatten(part['pt', -1][flav_mask][wp_mask])
    abseta = ak.flatten(part['abseta', -1][flav_mask][wp_mask])
    scale = ak.flatten(part.scale(-1)[flav_mask][wp_mask])

    return ((pt, abseta), scale)

def get_pt(part, flav, wp):
    flav_mask = part.flavor == flav
    wp_mask = part.pass_btag[flav_mask] >= wp
    pt = ak.flatten(part['pt', -1][flav_mask][wp_mask])
    scale = ak.flatten(part.scale(-1)[flav_mask][wp_mask])
    return (pt, scale)

def get_eta(part, flav, wp):
    flav_mask = part.flavor == flav
    wp_mask = part.pass_btag[flav_mask] >= wp
    abseta = ak.flatten(part['abseta', -1][flav_mask][wp_mask])
    scale = ak.flatten(part.scale(-1)[flav_mask][wp_mask])

    return (abseta, scale)


ptbins = bh.axis.Variable([25, 35, 50, 70, 90, 120, 150, 200])
etabins = bh.axis.Regular(5, 0, 2.5)
ptetas = [
    # light Jets
    GraphInfo("udsg_all", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 0, 0)),
    GraphInfo("udsg_L", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 0, 1)),
    GraphInfo("udsg_M", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 0, 2)),
    GraphInfo("udsg_T", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 0, 3)),
    # C-jets
    GraphInfo("c_all", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 4, 0)),
    GraphInfo("c_L", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 4, 1)),
    GraphInfo("c_M", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 4, 2)),
    GraphInfo("c_T", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 4, 3)),
    # B-jets
    GraphInfo("b_all", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 5, 0)),
    GraphInfo("b_L", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 5, 1)),
    GraphInfo("b_M", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 5, 2)),
    GraphInfo("b_T", "", (ptbins, etabins), lambda vg: get_pteta(vg.BJets, 5, 3)),
]
