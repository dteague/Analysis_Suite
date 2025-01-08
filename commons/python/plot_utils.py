#!/usr/bin/env python3
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import logging
from mpl_toolkits.axes_grid1 import make_axes_locatable
from contextlib import contextmanager
import warnings
import mplhep as hep

plt.style.use([hep.style.CMS, hep.style.firamath])

@contextmanager
def ratio_plot(filename, xlabel, binning, **kwargs):
    plot_inputs = {"nrows": 2, "ncols": 1, "sharex": True, 'figsize': (11,11),
                   "gridspec_kw": {"hspace": 0.1, "height_ratios": [3,1]}}
    fig, ax = plt.subplots(**plot_inputs)
    yield ax
    setup_ticks(*ax)
    axisSetup(ax[0], ax[1], xlabel=xlabel, binning=binning, **kwargs)
    ax[0].legend(ncols=2, fontsize='x-small')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()
    if hasattr(plot, "workdir"):
        filename = f"{plot.workdir}/{filename}"
    fig.savefig(filename, bbox_inches="tight", dpi=300)
    # if "extra_format" in kwargs:
        # fig.savefig(filename.with_suffix(f'.{kwargs["extra_format"]}'), bbox_inches="tight", dpi=300)
    plt.close(fig)

@contextmanager
def nonratio_plot(filename, xlabel, binning, **kwargs):
    fig, ax = plt.subplots()
    yield ax
    setup_ticks(ax)
    axisSetup(ax, xlabel=xlabel, binning=binning, **kwargs)
    if kwargs.get('legend', True):
        ax.legend()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()
    if hasattr(plot, "workdir"):
        filename = f"{plot.workdir}/{filename}"
    fig.savefig(filename, bbox_inches="tight", dpi=300)
    if "extra_format" in kwargs:
        fig.savefig(filename.with_suffix(f'.{kwargs["extra_format"]}'), bbox_inches="tight", dpi=300)

    plt.close(fig)

@contextmanager
def plot(filename, subplot_info=None, **kwargs):
    if subplot_info is None:
        subplot_info = {}
    fig, ax = plt.subplots(**subplot_info)
    yield ax
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()
    if hasattr(plot, "workdir"):
        filename = f"{plot.workdir}/{filename}"
    fig.savefig(filename, bbox_inches="tight", dpi=300)
    if "extra_format" in kwargs:
        fig.savefig(filename.with_suffix(f'.{kwargs["extra_format"]}'), bbox_inches="tight", dpi=300)
    plt.close(fig)

def plot_colorbar(cf, ax, barpercent=5):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=f"{barpercent}%", pad=0.05)
    plt.gcf().colorbar(cf, ax=cax, cax=cax)

def color_options(color):
    cvec = clr.to_rgb(color)
    dark = 0.3
    return {"color": color, "edgecolor": [i - dark if i > dark else 0.0 for i in cvec]}

def setup_ticks(pad, subpad=None):
    if subpad is not None:
        ticks(subpad)
        subpad.tick_params(direction="in")
    ticks(pad)


def ticks(pad):
    pad.minorticks_on()
    pad.tick_params(direction="in", length=9, top=True, right=True)
    pad.tick_params(direction="in", length=4, which='minor', top=True,
                    right=True)

def axisSetup(pad, subpad=None, xlabel="", binning=None, ratio_top=2.0, ratio_bot=0.0, zero_bot=True, normed=False,
              **kwargs):
    if subpad is None:
        xpad = pad
    else:
        xpad = subpad
        subpad.set_ylabel(kwargs.get("subpad_label", "Data/MC"))
        subpad.set_ylim(top=ratio_top, bottom=ratio_bot)

    if xlabel:
        xpad.set_xlabel(xlabel, loc="right")
    if binning is not None:
        xpad.set_xlim(binning[0], binning[-1])
    if zero_bot:
        pad.set_ylim(bottom=0.)

    xpad.ticklabel_format(useOffset=False)
    if 'pad_label' in kwargs:
        pad.set_ylabel(kwargs['pad_label'], loc='top')
    elif normed:
        pad.set_ylabel('Events/bins', loc='top')
    else:
        pad.set_ylabel('Events', loc='top')


def cms_label(ax, lumi=None, year=None, hasData=False):
    if lumi is None and year is None:
        hep.cms.label(ax=ax, label="Work in Progress", data=hasData)
        return
    if year is not None:
        from .constants import lumi
        lumi = lumi[year]
    hep.cms.label(ax=ax, lumi=lumi, label="Work in Progress", data=hasData)
