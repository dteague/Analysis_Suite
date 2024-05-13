#!/usr/bin/env python3
import numpy as np
from copy import copy
import uproot
import pandas as pd
from matplotlib import colors as clr
import boost_histogram as bh

from analysis_suite.commons.constants import all_eras
from analysis_suite.combine.card_maker import Card_Maker
from analysis_suite.combine.hist_writer import HistWriter
from analysis_suite.commons.histogram import Histogram
from analysis_suite.data.systs import systematics, get_shape_systs, dummy
import analysis_suite.commons.user as user
from analysis_suite.machine_learning.XGBoost import XGBoostMaker
from analysis_suite.commons.configs import get_list_systs, get_inputs, get_ntuple_info
from analysis_suite.combine.combine_wrapper import runCombine
from combine.scripts.make_2d_corr import run_single

def get_xgb(workdir, groupDict):
    inputs = get_inputs(workdir)
    ml_runner = XGBoostMaker(inputs.usevars, groupDict, region="signal", systName="Nominal")
    return ml_runner


def invert_groups(dic):
    final_dict = {}
    for group, members in dic.items():
        for member in members:
            final_dict[member] = group
    return final_dict

def parse_file(filename, groups, model):
    comb_ss_hists, sep_ss_hists, comb_ml_hists, sep_ml_hists = {}, {}, {}, {}
    axis = bh.axis.Variable([0.0, 0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 1.0, 1.1])
    bins = axis.edges
    s_ttt = ['tttj', 'tttw']
    systs = dict(get_shape_systs(year))

    with uproot.open(filename) as f:
        for member, group in groups.items():
            if member not in f:
                continue
            df = f[member].arrays(library='pd')
            ss_mask = (df['NMuons'] + df['NElectrons']) == 2
            mva_3top = model.predict(df, output_second).T[1]
            fourtop_mask = model.predict(df, output_first).T[1] < cut
            comb_group = group if group not in s_ttt else 'ttt'

            def fill(hist, mva, weight, mask):
                cr_vals = np.full(np.count_nonzero(~mask), 1.01)
                hist.fill(mva[mask], weight=weight[mask])
                hist.fill(cr_vals, weight=weight[~mask])

            def get_hist(dic, group, syst):
                if syst not in dic:
                    dic[syst] = dict()
                if group not in dic[syst]:
                    dic[syst][group] = Histogram("", axis)
                return dic[syst][group]

            for syst in f['weights'][member]:
                if syst.name == 'index':
                    continue
                systname = systs[syst.name]
                wgt = syst.array()

                fill(get_hist(sep_ss_hists, group, systname), mva_3top[ss_mask], wgt[ss_mask], fourtop_mask[ss_mask])
                fill(get_hist(sep_ml_hists, group, systname), mva_3top[~ss_mask], wgt[~ss_mask], fourtop_mask[~ss_mask])
                fill(get_hist(comb_ss_hists, comb_group, systname), mva_3top[ss_mask], wgt[ss_mask], fourtop_mask[ss_mask])
                fill(get_hist(comb_ml_hists, comb_group, systname), mva_3top[~ss_mask], wgt[~ss_mask], fourtop_mask[~ss_mask])
    return comb_ss_hists, sep_ss_hists, comb_ml_hists, sep_ml_hists

def update_dict(old_dict, update_dict):
    for key, subdict in update_dict.items():
        old_dict[key] = subdict

def setup_hists(workdir, output_first, output_second, cut, year):
    model = get_xgb(workdir, {})
    ginfo = get_ntuple_info("signal", remove='nonprompt_mc')
    groups = invert_groups(ginfo.setup_groups())

    comb_ss, sep_ss, comb_ml, sep_ml = parse_file(output_first/year/"test_Nominal_signal.root", groups, model)
    for f in (workdir/year).glob("test*"):
        if "Nominal" in f.name:
            continue
        tmp1, tmp2, tmp3, tmp4 = parse_file(f, groups, model)
        update_dict(comb_ss, tmp1)
        update_dict(sep_ss, tmp2)
        update_dict(comb_ml, tmp3)
        update_dict(sep_ml, tmp4)

    combine_dir = output_second/'combine'
    combine_dir.mkdir(exist_ok=True)
    with HistWriter(combine_dir / f'signal_{year}_Dilepton.root') as writer:
        for syst, hists in sep_ss.items():
            writer.add_syst(hists, syst=syst, blind=True)
    with HistWriter(combine_dir / f'signal_{year}_Multi.root') as writer:
        for syst, hists in sep_ml.items():
            writer.add_syst(hists, syst=syst, blind=True)
    with HistWriter(combine_dir / f'signal_{year}_Dilepton_sig.root') as writer:
        for syst, hists in comb_ss.items():
            writer.add_syst(hists, syst=syst, blind=True)
    with HistWriter(combine_dir / f'signal_{year}_Multi_sig.root') as writer:
        for syst, hists in comb_ml.items():
            writer.add_syst(hists, syst=syst, blind=True)
    # x10
    with HistWriter(combine_dir / f'signal_{year}_Dilepton_sig_x10.root') as writer:
        for syst, hists in comb_ss.items():
            if 'ttt' in hists:
                hists['ttt'].scale(10)
                hists['4top'].scale(10)
            writer.add_syst(hists, syst=syst, blind=True)
    with HistWriter(combine_dir / f'signal_{year}_Multi_sig_x10.root') as writer:
        for syst, hists in comb_ml.items():
            if 'ttt' in hists:
                hists['ttt'].scale(10)
                hists['4top'].scale(10)
            writer.add_syst(hists, syst=syst, blind=True)
    # x100
    with HistWriter(combine_dir / f'signal_{year}_Dilepton_sig_x100.root') as writer:
        for syst, hists in comb_ss.items():
            if 'ttt' in hists:
                hists['ttt'].scale(10)
                hists['4top'].scale(10)
            writer.add_syst(hists, syst=syst, blind=True)
    with HistWriter(combine_dir / f'signal_{year}_Multi_sig_x100.root') as writer:
        for syst, hists in comb_ml.items():
            if 'ttt' in hists:
                hists['ttt'].scale(10)
                hists['4top'].scale(10)
            writer.add_syst(hists, syst=syst, blind=True)
    # x1000
    with HistWriter(combine_dir / f'signal_{year}_Dilepton_sig_x1000.root') as writer:
        for syst, hists in comb_ss.items():
            if 'ttt' in hists:
                hists['ttt'].scale(10)
                hists['4top'].scale(10)
            writer.add_syst(hists, syst=syst, blind=True)
    with HistWriter(combine_dir / f'signal_{year}_Multi_sig_x1000.root') as writer:
        for syst, hists in comb_ml.items():
            if 'ttt' in hists:
                hists['ttt'].scale(10)
                hists['4top'].scale(10)
            writer.add_syst(hists, syst=syst, blind=True)



def make_card(combinedir, year, region, signals, nosyst=False, sig=False):
    remove = ['nonprompt_mc']
    if region == "Multi":
        remove.append('charge_flip')
    ginfo = get_ntuple_info('signal', remove=remove)
    group_list = [g for g in ginfo.setup_groups().keys() if g not in signals]
    if sig:
        group_list.remove('tttj')
        group_list.remove('tttw')

    with Card_Maker(combinedir, year, region, signals, group_list, 'signal', nosyst=nosyst) as card:
        card.write_preamble()
        if not nosyst:
            card.write_systematics(systematics)
        else:
            card.write_systematics(dummy)
        card.add_stats()
    if nosyst:
        region += "_nosyst"
    return f"{region}=signal_{year}_{region}_card.txt "

def setup_cards(combinedir, year, signals, nosyst=False, sig=False):
    sig_reg = '_sig' if sig else ''
    combine_cmd = "combineCards.py "
    combine_cmd += make_card(output_second/'combine', year, 'Dilepton'+sig_reg, signals, nosyst=nosyst, sig=sig)
    combine_cmd += make_card(output_second/'combine', year, 'Multi'+sig_reg, signals, nosyst=nosyst, sig=sig)
    extra = ('_sig' if sig else "") + ("_nosyst" if nosyst else "")
    final_card = f"final_{year}" + extra + "_card.txt"
    combine_cmd += f'> {final_card}'
    runCombine(combine_cmd)
    return final_card

def sed_file(indir, infile, outfile, intext, outtext):
    with open(indir/infile, 'r') as fin:
        text = fin.read()
    with open(indir/outfile, 'w') as fout:
        fout.write(text.replace(intext, outtext))

cut = 0.97

workdir = user.workspace_area/"October_run"
output_base = user.analysis_area/"corr_test"
output_first = output_base / "4top_first_Dec"
output_second = output_base / f"final_Dec_{cut}"
combine_dir = output_second/'combine'
runCombine.work_dir = combine_dir
years = ['2016pre', '2016post', '2017', '2018']
# years = ['2018']

sep_nosyst = "combineCards.py "
sep_syst = "combineCards.py "
comb_nosyst = "combineCards.py "
comb_syst = "combineCards.py "


for year in years:
    pass
    # setup_hists(workdir, output_first, output_second, cut, year)
    # exit()
    # card = setup_cards(combine_dir, year, ['tttj', 'tttw', '4top'], nosyst=True)
    # sep_nosyst += f'era{year}={card} '
    # run_single(combine_dir, card, year, 'tttj', '4top', no_syst=True, sig_up=100, other_up=5)
    # run_single(combine_dir, card, year, 'tttw', '4top', no_syst=True, sig_up=100, other_up=5)

    # card = setup_cards(combine_dir, year, ['tttj', 'tttw', '4top'], nosyst=False)
    # sep_syst += f'era{year}={card} '
    # run_single(combine_dir, card, year, 'tttj', '4top', no_syst=False, sig_up=100, other_up=5)
    # run_single(combine_dir, card, year, 'tttw', '4top', no_syst=False, sig_up=100, other_up=5)

    # card= setup_cards(combine_dir, year, ['ttt', '4top'], nosyst=True, sig=True)
    # comb_nosyst += f'era{year}={card} '
    # run_single(combine_dir, card, year, 'ttt', '4top', no_syst=True, sig_up=100, other_up=5)

    # card = setup_cards(combine_dir, year, ['ttt', '4top'], nosyst=False, sig=True)
    # comb_syst += f'era{year}={card} '
    # run_single(combine_dir, card, year, 'ttt', '4top', no_syst=False, sig_up=100, other_up=5)

if years == all_eras:
    # sep_nosyst += "> final_all_nosyst_card.txt"
    # runCombine(sep_nosyst)
    # run_single(combine_dir, 'final_all_nosyst_card.txt', 'all', 'tttj', '4top', no_syst=True, sig_up=100, other_up=5)
    # run_single(combine_dir, 'final_all_nosyst_card.txt', 'all', 'tttw', '4top', no_syst=True, sig_up=100, other_up=5)

    # sep_syst += "> final_all_card.txt"
    # runCombine(sep_syst)
    # run_single(combine_dir, 'final_all_card.txt', 'all', 'tttj', '4top', no_syst=False, sig_up=100, other_up=5)
    # run_single(combine_dir, 'final_all_card.txt', 'all', 'tttw', '4top', no_syst=False, sig_up=100, other_up=5)

    # comb_nosyst += "> final_all_sig_nosyst_card.txt"
    # runCombine(comb_nosyst)
    # run_single(combine_dir, 'final_all_sig_nosyst_card.txt', 'all', 'ttt', '4top', no_syst=True, sig_up=100, other_up=5)

    # comb_syst += "> final_all_sig_card.txt"
    # runCombine(comb_syst)
    # run_single(combine_dir, 'final_all_sig_card.txt', 'all', 'ttt', '4top', no_syst=False, sig_up=100, other_up=5)

    # Run 10x
    sed_file(combine_dir, 'final_all_sig_nosyst_card.txt', 'final_all_sig_x10_nosyst_card.txt', '_sig', "_sig_x10")
    run_single(combine_dir, 'final_all_sig_x10_nosyst_card.txt', 'all', 'ttt', '4top', no_syst=True, sig_up=15, other_down=0.3, other_up=1.7, extra='_x10')

    sed_file(combine_dir, 'final_all_sig_card.txt', 'final_all_sig_x10_card.txt', '_sig', "_sig_x10")
    run_single(combine_dir, 'final_all_sig_x10_card.txt', 'all', 'ttt', '4top', no_syst=False, sig_up=15, other_down=0.3, other_up=1.7, extra='_x10')

    # Run 100x
    sed_file(combine_dir, 'final_all_sig_nosyst_card.txt', 'final_all_sig_x100_nosyst_card.txt', '_sig', "_sig_x100")
    run_single(combine_dir, 'final_all_sig_x100_nosyst_card.txt', 'all', 'ttt', '4top', no_syst=True, sig_up=4, other_down=0.8, other_up=1.2, extra='_x100')

    sed_file(combine_dir, 'final_all_sig_card.txt', 'final_all_sig_x100_card.txt', '_sig', "_sig_x100")
    run_single(combine_dir, 'final_all_sig_x100_card.txt', 'all', 'ttt', '4top', no_syst=False, sig_up=4, other_down=0.8, other_up=1.2, extra='_x100')

    # Run 1000x
    sed_file(combine_dir, 'final_all_sig_nosyst_card.txt', 'final_all_sig_x1000_nosyst_card.txt', '_sig', "_sig_x1000")
    run_single(combine_dir, 'final_all_sig_x1000_nosyst_card.txt', 'all', 'ttt', '4top', no_syst=True, sig_down=0, sig_up=2, other_down=0.9, other_up=1.1,extra='_x1000')

    sed_file(combine_dir, 'final_all_sig_card.txt', 'final_all_sig_x1000_card.txt', '_sig', "_sig_x1000")
    run_single(combine_dir, 'final_all_sig_x1000_card.txt', 'all', 'ttt', '4top', no_syst=False,  sig_down=0, sig_up=2, other_down=0.9, other_up=1.1, extra='_x1000')
