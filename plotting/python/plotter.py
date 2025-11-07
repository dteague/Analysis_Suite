#!/usr/bin/env python3
import numpy as np
from sklearn.metrics import roc_curve
from pathlib import Path
from copy import copy
# from prettytable import PrettyTable
import warnings

from analysis_suite.commons.histogram import Histogram
from analysis_suite.commons.constants import lumi
import analysis_suite.commons.configs as config
from analysis_suite.commons.plot_utils import plot, nonratio_plot, ratio_plot, cms_label

from .stack import Stack

def get_a_dictval(hists):
    if not hists:
        return None
    else:
        return next(iter(hists.values()))

def signif(x, p):
    if int(x) == x:
        return int(x)
    mags = 10 ** (p - 1 - np.floor(np.log10(x)))
    return np.round(x * mags) / mags


class Plotter:
    color_list = ['#ed2839',  '#0acdff', "#2fbf71",'#f1a340', '#d8b365']

    def __init__(self, ntuple, year, sig="", bkg=None, data='data', out_format='pdf', outdir=".", from_ntuple=True):
        # Setup groups
        self.ntuple = ntuple
        self.ginfo = ntuple.get_info(keep_dd_data=from_ntuple)
        if sig is None:
            self.signame = ""
            self.sigs = []
        elif isinstance(sig, str):
            self.signame = sig
            self.sigs = [sig]
        else:
            self.signame = sig[0]
            self.sigs = sig[1]
        self.data = data
        if bkg is None:
            bkg = []
        elif bkg == "all":
            self.bkg = [g for g in self.ginfo.get_groups() if g not in self.sigs and g != data]
        else:
            self.bkg = bkg
        # Setup output name
        self.ft = out_format
        self.outdir = Path(outdir)
        self.outdir.mkdir(exist_ok=True, parents=True)
        self.extra = str(year)
        self.year = year
        self.lumi = lumi[year]


    def fill_hist_groups(self, hists):
        out_hists = {}
        for member, hist in hists.items():
            group = self.ntuple.get_group_name(member, None, False)
            if group is None:
                continue
            elif group not in out_hists:
                out_hists[group] = hist
            else:
                out_hists[group] += hist
        return out_hists

    def set_hist_groups(self, hists, fix_negative=False):
        out_hists = dict()
        for group, hist in hists.items():
            hist.color = self.ginfo.get_color(group)
            hist.plot_label = self.ginfo.get_legend_name(group)
            hist.axis_name = hist.axis_name
            if fix_negative and min(hist.vals) < 0:
                hist.values()[:] = np.abs(hist.vals)
            out_hists[group] = hist
        if self.signame and self.signame not in out_hists:
            for sig in self.sigs:
                if self.signame not in out_hists:
                    out_hists[self.signame] = out_hists.pop(sig)
                    out_hists[self.signame].plot_label = self.ginfo.get_legend_name(self.signame)
                else:
                    out_hists[self.signame] += out_hists.pop(sig)
        return out_hists

    def set_extra_text(self, extra=None):
        if extra is None:
            self.extra=str(self.year)
        else:
            self.extra = f"{self.year}_{extra}"

    def set_year(self, year):
        self.extra = self.extra.replace(self.year, year)
        self.year = year
        self.lumi = lumi[year]


    def make_stack(self, hists):
        stack = Stack()
        for group in self.bkg:
            stack += self.hists[name][group]
        return stack

    def plot_hist(self, name, hist, **kwargs):
        self.plot_hists({name: hist}, **kwargs)

    def plot_hists(self, hist_dict, plot_type="stack", **kwargs):
        for name, hists in hist_dict.items():
            group_hists = self.set_hist_groups(hists)
            if plot_type == "shape":
                self.plot_shape(group_hists, name, **kwargs)
            elif plot_type == "stack_shape":
                self.plot_stack_shape(group_hists, name, **kwargs)
            elif plot_type == 'stack':
                self.plot_stack(group_hists, name, **kwargs)
            elif plot_type == 'signals':
                self.plot_signals(group_hists, name, **kwargs)
            else:
                raise AttributeError(f"Not correct Output type: {plot_type}")

    def add_total_systs(self, stack, nom_hists, hist_systs):
        tot_syst = 0
        for syst_dict in hist_systs:
            for group, syst_shift in syst_dict.items():
                if group in self.sigs:
                    continue
                if isinstance(syst_shift, float):
                    stack.add_syst(syst_shift*nom_hists[group])
                    # print(group, syst_shift, (syst_shift*nom_hists[group]).vals)
                    tot_syst += ((syst_shift*nom_hists[group]).integral()/stack.integral())**2
                else:
                    stack.add_syst(syst_shift)
                    # print(syst_shift.vals)
                    tot_syst += (syst_shift.integral()/stack.integral())**2
        print(tot_syst, np.sqrt(tot_syst))

    def add_group_systs(self, hist, hist_group, hist_systs):
        for syst_dict in hist_systs:
            for group, syst_shift in syst_dict.items():
                if hist_group != group:
                    continue
                if isinstance(syst_shift, float):
                    hist.add_syst(syst_shift*hist)
                else:
                    hist.add_syst(syst_shift)

    def plot_stack(self, hists, name, **kwargs):
        a_hist = get_a_dictval(hists)
        if a_hist is None:
            print(f'Problem plotting {name}: skipping')
            return
        axis = a_hist.axes
        axis_name = a_hist.axis_name
        if "Regular" in repr(axis[0]):
            yaxis_label = f'Events / {signif(axis[0].widths[0], 2):g}'
            yaxis_label += " GeV" if "GeV" in axis_name else " units"
            normed = False
        else:
            yaxis_label = f'Events / bin'
            normed = True

        if 'normed' in kwargs:
            normed = kwargs['normed']
        data = hists.pop(self.data, Histogram())
        signal = hists.pop(self.signame, Histogram())
        make_ratio = bool(data)

        stack = Stack(*axis)
        ratio = Histogram(*axis, color="black")
        band = Histogram(*axis, color="black")
        for bkg in self.bkg:
            if bkg in hists:
                stack += hists[bkg]

        if signal:
            width = signal.axis.widths
            if normed:
                sig_scale = np.max(stack.vals/width)/np.max(signal.vals/width)
            else:
                sig_scale = np.max(stack.vals)/np.max(signal.vals)
            rounder = 25 if sig_scale >= 25 else 1
            sig_scale = int(sig_scale//rounder*rounder)
            signal.scale(sig_scale, for_plot=True)

        plot_syst = 'syst_hists' in kwargs or 'syst_err' in kwargs
        if 'syst_hists' in kwargs and not isinstance(kwargs['syst_hists'], bool):
            self.add_total_systs(stack, hists, kwargs['syst_hists'])
        elif plot_syst:
            stack.add_syst()
            stack.variances()[:] = kwargs['syst_err']

        # Lower plot
        plot_func = ratio_plot if make_ratio else nonratio_plot
        if make_ratio:
            band = stack.get_self_ratio()
            ratio += data/stack

        file_name = self.outdir/f"{name}_{self.extra}.{self.ft}"
        with plot_func(file_name, axis_name, a_hist.get_xrange(),
                       pad_label=yaxis_label, lumi=self.lumi, data=bool(data),
                       **kwargs) as ax:
            if make_ratio:
                pad, subpad = ax
            else:
                pad, subpad = ax, None

            if kwargs.get('log', False):
                pad.set_yscale('log')
            if kwargs.get('text', False):
                pad.text(0.05, 0.95, kwargs.get('text'), transform=pad.transAxes, ha='left', va='top', size=30)
            stack.plot_stack(pad, normed=normed)
            signal.plot_hist(pad, normed=normed, linestyle='--', color='k')
            data.plot_points(pad, data=True, normed=normed)
            # Ratio
            ratio.plot_points(subpad, data=True)
            # Error bands
            if plot_syst:
                stack.plot_error_band(pad, color='k', hatch=r'///', systs=True, normed=normed)
                band.plot_error_band(subpad, color='k', hatch=r'///', systs=True)
            else:
                stack.plot_error_band(pad, normed=normed)
                band.plot_error_band(subpad)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # cms_label(pad, lumi=self.lumi, hasData=data)


    def plot_stack_shape(self, hists, outname, plot_sigs=None, **kwargs):
        if plot_sigs is None:
            raise Exception("Need list of signals to plot")
        a_hist = get_a_dictval(hists)
        axis = a_hist.axes

        signals = []
        stack = Stack(*axis)
        for group, hist in hists.items():
            if group in plot_sigs:
                signal = hist
                signals.append(signal/(signal.integral()))
            else:
                stack += hist
        stack.normalize()

        file_name = f"{self.outdir}/{outname}_{self.extra}.{self.ft}"
        with nonratio_plot(file_name, a_hist.axis_name, a_hist.get_xrange(),
                           pad_label="A.U.", **kwargs) as pad:
            stack.plot_stack(pad)
            # stack.plot_error_band(pad)
            for signal in signals:
                signal.plot_hist(pad)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cms_label(pad, lumi=self.lumi)

    def plot_shape(self, hists, outname, shapes=None, **kwargs):
        a_hist = get_a_dictval(hists)
        axis = a_hist.axes

        if shapes is None:
            shapes_hists = hists
        else:
            shape_hists = {group: Histogram(*axis, name=group, color=Plotter.color_list[i])
                           for i, group in enumerate(shapes)}
            groups = list(hists.keys())
            other_group = None
            for shape_group, names in shapes.items():
                if names is None:
                    other_group = shape_group
                    continue
                for name in names:
                    shape_hists[shape_group] += hists[name]
                    groups.pop(groups.index(name))
            for group in groups:
                shape_hists[other_group] += hists[group]

        file_name = f"{self.outdir}/{outname}_shape_{self.extra}.{self.ft}"
        with nonratio_plot(file_name, a_hist.axis_name, a_hist.get_xrange(), **kwargs) as pad:
            for group, hist in shape_hists.items():
                hist.plot_shape(pad)
            cms_label(pad, lumi=self.lumi, label='')

    def plot_signals(self, hists, outname, **kwargs):
        a_hist = get_a_dictval(hists)
        axis = a_hist.axes

        file_name = f"{self.outdir}/{outname}_sigs_{self.extra}.{self.ft}"
        with nonratio_plot(file_name, a_hist.axis_name, a_hist.get_xrange(), **kwargs) as pad:
            for hist in hists.values():
                hist.plot_hist(pad)
            cms_label(pad, lumi=self.lumi)
