#!/usr/bin/env python3
import numpy as np
from copy import copy
import uproot
import pandas as pd
from matplotlib import colors as clr
import boost_histogram as bh

from analysis_suite.commons.constants import lumi
from analysis_suite.commons.plot_utils import plot, hep
from analysis_suite.combine.card_maker import Card_Maker
from analysis_suite.combine.hist_writer import HistWriter
from analysis_suite.commons.histogram import Histogram
from analysis_suite.data.systs import systematics, get_shape_systs, dummy
from analysis_suite.commons.configs import get_list_systs, get_inputs, get_ntuple_info
import analysis_suite.commons.user as user
from analysis_suite.machine_learning.XGBoost import XGBoostMaker

output_base = user.workspace_area/"Jan_ptcut"
s_4top = ['4top']
s_ttt = ['tttj', 'tttw']

groupDict = {
    'Background': np.array(['ttw', 'ttz', 'ttz_m1-10', 'tth', 'ttww', 'ttwz', 'ttzz',
                            'tthh', 'ttzh', 'ttwh', 'wwz', 'wzz', 'www', 'zzz', 'wzg', 'wwg',
                            'vh2nonbb', 'zz4l', 'wpwpjj_ewk', 'ww_doubleScatter', 'wzTo3lnu',
                            'st_twll', 'ggh2zz', 'tzq', 'ttg_hadronic', 'ttg_singleLept',
                            'ttg_dilep', 'zg', 'wg', 'tg']),
    'NotTrained': np.array(['nonprompt']),
    'OnlyTrain': np.array(['ttbar_2l2n', 'ttbar_semilep', 'ttbar_hadronic', 'wjets_ht70-100',
                           'wjets_ht100-200', 'wjets_ht200-400', 'wjets_ht400-600',
                           'wjets_ht600-800', 'wjets_ht800-1200', 'wjets_ht1200-2500',
                           'wjets_ht2500-Inf'])}


def get_xgb(workdir, groupDict):
    # params = {
    #     'colsample_bytree': 0.5283131479037887, 'eta': 0.037489501191992,
    #     'eval_metric': 'logloss', 'gamma': 5.439613578244588, 'max_depth': 5,
    #     'min_child_weight': 4.0, 'n_estimators': 500, 'subsample': 0.7009163217893063,
    # }
    params = {
        'colsample_bytree': 0.5283131479037887, 'eta': 0.037489501191992,
        'eval_metric': 'logloss', 'gamma': 5.439613578244588, 'max_depth': 3,
        'min_child_weight': 1.0, 'n_estimators': 500, 'subsample': 0.9
    }
    inputs = get_inputs(workdir)
    ml_runner = XGBoostMaker(inputs.usevars, groupDict, region="signal", systName="Nominal")
    ml_runner.update_params(params)
    return ml_runner

def get_weights(indir, year='2018'):
    weights = dict()
    with uproot.open(indir/year/"test_Nominal_signal.root") as f:
        wgt_dir = f['weights']
        weights = {samp[:-2]: wgt_dir[samp].arrays(library='pd') for samp in wgt_dir}
    return weights

def get_set(typ, indir, model, name='Signal', year='2018', cut=0.99):
    extra = "_signal" if typ == "test" else ""

    final_df = pd.DataFrame()
    if typ == 'test':
        filename = f"{typ}_Nominal_signal.root"
        file_dir = indir/year
    else:
        filename = f"{typ}_Nominal.root"
        file_dir = indir/'train_files'
    with uproot.open(file_dir/filename) as f:
        for key, cls in f.classnames().items():
            if "/" in key or cls != "TTree":
                continue
            df = f[key].arrays(library='pd')
            classID = np.unique(df['classID'])[0]
            print(key, np.sum(df.split_weight), len(df), np.unique(df.split_weight), np.sum(df.scale_factor), np.sum(df.split_weight)/np.sum(df.scale_factor))
            if "4top" in key:
                if typ != "test": continue
                df['classID'] = 0 if classID == 1 else 1
            elif "tttj" in key or 'tttw' in key:
                df['classID'] = 0 if classID == 1 else 1
            if "4top_sig" not in df.columns:
                pred = model.predict(df, indir).T[1]
                df['4top_sig'] = pred
            # print(key, np.unique(df["sampleName"]), len(df))
            final_df = pd.concat([final_df, df])
    final_df = final_df[final_df["4top_sig"] < cut]
    return final_df


def train(model, output, years=['2018'], signal="Signal", output_name="Signal"):
    model.train()
    model.plot_training_progress()
    for year in years:
        model.apply_model(year, get_auc=True)
        model.roc_curve(year, signal=signal)
        model.plot_overtrain(year, signal)
        model.plot_train_breakdown(year, signal)
        model.plot_train_breakdown(year, signal, use_test=False)
        model.output(output, year, output_name)
    model.output_train(output, output_name)


def run_4top_first(workdir, output, groupDict, years=['2018']):
    groupDict = copy(groupDict)
    # groupDict["Signal"] = s_ttt
    # groupDict['4top'] = s_4top
    groupDict["Background"] = np.append(groupDict["Background"], s_ttt)
    groupDict["Signal"] = s_4top

    model = get_xgb(workdir, groupDict)
    model.set_outdir(output)
    for year in years:
        (output/year).mkdir(exist_ok=True, parents=True)
        model.setup_year(workdir, year)
    train(model, output, output_name="4top_sig", years=years)

def run_3top_second(workdir, indir, output, groupDict, cut, years=['2018']):
    groupDict = copy(groupDict)
    groupDict["Signal"] = s_ttt
    groupDict["NotTrained"] = np.append(groupDict["NotTrained"], s_4top)
    model = get_xgb(workdir, groupDict)
    model.set_outdir(output)
    for year in years:
        (output/year).mkdir(exist_ok=True, parents=True)
        model.test_sets[year] = get_set("test", indir, model, cut=cut)
        model.test_weights[year] = get_weights(indir)
    model.train_set = get_set("train", indir, model, cut=cut)
    model.validation_set = get_set("validation", indir, model, cut=cut)

    train(model, output, output_name="ttt_sig", years=years)

def plot_final_hist(workdir, output_first, output_second, groupDict, cut, year):
    model = get_xgb(workdir, groupDict)
    ginfo = get_ntuple_info("signal", remove='nonprompt_mc')
    groups = ginfo.setup_groups()

    nbins = 31
    axis = bh.axis.Variable([0.0, 0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 1.0, 1.1])
    bins = axis.edges
    out_hists = {key: Histogram("", axis) for key in groups.keys()}
    with uproot.open(output_first/year/"test_Nominal_signal.root") as f:
        for group, members in groups.items():
            for member in members:
                if member not in f:
                    print(f"{member} not found")
                    continue
                df = f[member].arrays(library='pd')
                weight = f['weights'][member]["Nominal"].array()
                mva_4top = model.predict(df, output_first).T[1]
                mva_3top = model.predict(df, output_second).T[1]

                mask = mva_4top < cut

                out_hists[group].fill(mva_3top[mask], weight=weight[mask])
                out_hists[group].fill(np.full(np.count_nonzero(mva_4top >= cut), 1.01), weight=weight[~mask])
                # hist = np.histogram(mva_3top[mask], bins, weights=weight[mask])[0]
                # hist[-1] += np.sum(weight[~mask])
                # out_hists[group] += hist

    combine_dir = output_second/'combine'
    combine_dir.mkdir(exist_ok=True)
    with HistWriter(combine_dir / f'signal_{year}_signal.root') as writer:
        writer.add_syst(out_hists, syst="Nominal", blind=True)

    signals = ['ttt', '4top']
    group_list = [g for g in groups.keys() if g not in signals]
    with Card_Maker(combine_dir, year, 'signal', signals, group_list, 'signal', nosyst=True) as card:
        card.write_preamble()
        card.write_systematics(dummy)
        card.add_stats()
    card = output_second/'combine'/f'signal_{year}_signal_nosyst_card.txt'
    card.rename(card.parent/f"final_{year}_nosyst_card.txt")

    signals = ['ttt']
    group_list = [g for g in groups.keys() if g not in signals]
    with Card_Maker(combine_dir, year, 'signal', ['ttt'], group_list, 'signal', nosyst=True) as card:
        card.write_preamble()
        card.write_systematics(dummy)
        card.add_stats()
    card = output_second/'combine'/f'signal_{year}_signal_nosyst_card.txt'
    card.rename(card.parent/f"final_{year}_nosyst_card_single.txt")


    with plot(output_second/f'final_plot_{year}.png') as ax:
        signal = out_hists.pop('ttt').vals/axis.widths
        tttt = out_hists.pop('4top').vals/axis.widths
        stack = np.array([hist.vals/axis.widths for hist in out_hists.values()])
        colors = [ginfo.get_color(key) for key in out_hists.keys()]
        names = [f'${ginfo.get_legend_name(key)}$' for key in out_hists.keys()]
        ax.hist(x=bins[:-1], bins=bins, weights=signal/np.sum(signal), color='r', label="ttt", linewidth=3, histtype='step')
        ax.hist(x=bins[:-1], bins=bins, weights=tttt/np.sum(tttt), color='orange', label="4top", linewidth=3, histtype='step')
        n, bins, patches = ax.hist(
            weights=stack.T/np.sum(stack), bins=bins, x=np.tile(bins[:-1], (len(stack), 1)).T,
            label=names, histtype='stepfilled', stacked=True,
            color=colors
        )
        # Apply patch to edge colors
        dark = 0.3
        for p, c in zip(patches, colors):
            ec =  [i - dark if i > dark else 0.0 for i in clr.to_rgb(c)]
            if isinstance(p, list):
                p[0].set_ec(ec)
            else:
                p.set(ec=ec)

        ax.set_xlim(0., bins[-1])
        ax.set_xlabel("$disc_{BDT}$")
        ax.set_ylabel("A.U.")
        ax.legend()
        hep.cms.label(ax=ax, lumi=lumi[year], label="Preliminary")

    # model.use_vars.append("4top_sig")

output_first = output_base / "4top_first_Dec"
jcut = 0.93
output_second = output_base / f"final_Dec_{cut}"
years = ['2016pre', '2016post', '2017', '2018']
run_4top_first(workdir, output_first, groupDict, years=years)
run_3top_second(workdir, output_first, output_second, groupDict, cut, years=years)
for year in years:
    plot_final_hist(workdir, output_first, output_second, groupDict, cut, year)
