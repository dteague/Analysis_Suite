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
    ax[0].legend()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()
    if hasattr(plot, "workdir"):
        filename = f"{plot.workdir}/{filename}"
    fig.savefig(filename, bbox_inches="tight")
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
    fig.savefig(filename, bbox_inches="tight")
    plt.close(fig)

@contextmanager
def plot(filename, *args, **kwargs):
    fig, ax = plt.subplots(*args, **kwargs)
    yield ax
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()
    if hasattr(plot, "workdir"):
        filename = f"{plot.workdir}/{filename}"
    fig.savefig(filename, bbox_inches="tight")
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
    xpad = pad if subpad is None else subpad
    if xlabel:
        xpad.set_xlabel(xlabel)
    if binning is not None:
        xpad.set_xlim(binning[0], binning[-1])
    if zero_bot:
        pad.set_ylim(bottom=0.)

    if 'pad_label' in kwargs:
        pad.set_ylabel(kwargs['pad_label'])
    elif normed:
        pad.set_ylabel('Events/bins')
    else:
        pad.set_ylabel('Events')
    right_align_label(pad.get_yaxis(), isYaxis=True)
    right_align_label(xpad.get_xaxis())

    if subpad is not None:
        subpad.set_ylabel(kwargs.get("subpad_label", "Data/MC"))
        subpad.set_ylim(top=ratio_top, bottom=ratio_bot)


def right_align_label(axis, isYaxis=False):
    label = axis.get_label()
    if isYaxis:
        label.set_y(1.0)
    else:
        label.set_x(1.0)
    label.set_horizontalalignment('right')
    axis.set_label(label)
