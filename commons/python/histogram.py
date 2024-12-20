import numpy as np
from copy import copy
from matplotlib import colors as clr
import boost_histogram as bh
from boost_histogram.accumulators import WeightedSum as bh_weights

import warnings
warnings.simplefilter("ignore", UserWarning)
from scipy.stats import beta

class Histogram:
    def __init__(self, *args, **kwargs):
        if len(args) == 0:
            args = (bh.axis.Regular(1, 0, 1),)
        self.hist = bh.Histogram(*args, storage=bh.storage.Weight())
        self.breakdown = dict()
        self.color = kwargs.get('color', 'k')
        self.name = kwargs.get('name', "")
        self.draw_scale = 1

    def __add__(self, right):
        return Histogram(*self.meta).set(self.hist + right.hist)

    def __sub__(self, right):
        return self + (-1)*right

    def __mul__(self, right):
        return Histogram(*self.meta).set(right*self.hist)

    def __rmul__(self, right):
        return right*self

    def __imul__(self, right):
        self.hist *= right
        return self

    def __iadd__(self, right):
        if isinstance(right, Histogram):
            self._set_hist(right.hist)
            self.breakdown.update({mem: bh_weights()
                                   for mem in right.breakdown.keys()
                                   if mem not in self.breakdown})
            for mem, info in right.breakdown.items():
                self.breakdown[mem] += info
        elif isinstance(right, bh.Histogram):
            self._set_hist(right)
        return self

    def _set_hist(self, hist):
        if not self:
            self.hist = copy(hist)
        else:
            self.hist += hist

    def __truediv__(self, denom):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ratio = np.nan_to_num(self.vals/denom.vals)
            error2 = ratio**2*(self.err_ratio + denom.err_ratio)

        return_obj = Histogram(*self.hist.axes)
        return_obj.hist.values()[:] = ratio
        return_obj.hist.variances()[:] = error2
        return_obj.set_metadata(self)
        return return_obj

    def __bool__(self):
        return not self.hist.empty()

    def __getstate__(self):
        return {
            "hist": self.hist,
            "color": self.color,
            "name": self.name,
            "breakdown": self.breakdown,
            "draw_sc": self.draw_sc
        }

    def __setstate__(self, state):
        self.hist = state["hist"]
        self.color = state["color"]
        self.name = state["name"]
        self.breakdown = state["breakdown"]
        self.draw_sc = state["draw_sc"]

    def __getattr__(self, attr):
        if attr == 'axis':
            return self.hist.axes[0]
        elif attr == 'axes':
            return self.hist.axes
        elif attr == 'vals':
            return self.hist.view().value
        elif attr == 'err':
            return np.sqrt(self.hist.view().variance)
        elif attr == 'sumw2':
            return self.hist.view().variance
        elif attr == 'err_ratio':
            return np.nan_to_num(self.sumw2/(self.vals**2+1e-6))
        else:
            raise Exception(f"{attr} was not found!")

    def set(self, hist):
        if isinstance(hist, Histogram):
            self.hist = hist.hist
        elif isinstance(hist, bh.Histogram):
            self.hist = hist
        else:
            raise Exception("Problem")

    def meta(self):
        return {
            "color": self.color,
            "name": self.name,
        }

    @staticmethod
    def efficiency(top, bot, asymm=False):
        alf = (1-0.682689492137)/2
        aa = top.vals*bot.vals/(bot.sumw2+1e-6)+1
        bb = (bot.vals-top.vals)*bot.vals/(bot.sumw2+1e-6)+1

        lo = np.array([beta.ppf(alf, p, t) for p, t in zip(aa, bb)])
        hi = np.array([beta.ppf(1 - alf, p, t) for p, t in zip(aa, bb)])
        if asymm:
            eff = lo
            error2 = ((hi-lo)/2)**2
        else:
            eff = np.array([beta.mean(p, t) for p, t in zip(aa, bb)])
            error2 = ((eff-lo)**2 + (hi-eff)**2)/2

        return_obj = Histogram(*top.hist.axes)
        return_obj.hist.values()[:] = eff
        return_obj.hist.variances()[:] = error2
        return return_obj

    def move_overflow(self):
        if len(self.axes) == 1:
            self.hist[-1] += self.hist[bh.overflow]
            self.hist[bh.overflow] = (0, 0)
        else:
            last_y = -1 if self.axes.size[1] > 1 else 0
            last_x = -1 if self.axes.size[0] > 1 else 0
            for i in range(self.axes.size[0]):
                self.hist[i, last_y] += self.hist[i, bh.overflow]
                self.hist[i, bh.overflow] = (0, 0)
            for i  in range(self.axes.size[1]):
                self.hist[last_x, i] += self.hist[bh.overflow, i]
                self.hist[bh.overflow, i] = (0, 0)
            self.hist[last_x, last_y] += self.hist[bh.overflow, bh.overflow]
            self.hist[bh.overflow, bh.overflow] = (0, 0)

    def project(self, ax):
        new_hist = Histogram(self.hist.axes[0])
        new_hist.hist = self.hist.project(ax)
        new_hist.set_metadata(self)
        return new_hist

    def fill(self, *vals, weight, flow=True, member=None):
        self.hist.fill(*vals, weight=weight)
        if member is not None:
            self.breakdown[member] = bh_weights().fill(weight)
        if flow and sum(self.hist.shape) > 1:
            self.move_overflow()

    def get_name(self):
        if self.draw_scale == 1:
            return self.name
        else:
            str_scale = str(scale) if isinstance(scale, int) else f'{scale:0.2f}'
            return f"{self.name} x {str_scale}"

    def set_plot_details(self, name, color):
        name = f'${name}$' if '\\' in name else name
        self.name = name
        self.color = color
            # name = group_info.get_legend_name(self.group)
            # self.name = f'${name}$' if '\\' in name else name
            # self.color = group_info.get_color(self.group)

    def darkenColor(self, color=None):
        if color is None:
            color = self.color
        cvec = clr.to_rgb(color)
        dark = 0.3
        return [i - dark if i > dark else 0.0 for i in cvec]

    def get_xrange(self):
        return [self.axis.edges[0], self.axis.edges[-1]]

    def scale(self, scale, for_plot=False):
        if forPlot:
            self.draw_scale *= scale
        else:
            self.hist *= scale
            for mem, info in self.breakdown.items():
                self.breakdown[mem] *= scale

    def integral(self, flow=True):
        return self.hist.sum(flow=flow).value

    def plot_points(self, pad, normed=False, **kwargs):
        if not self or pad is None:
            return
        mask = self.vals > 0.
        if normed:
            vals = self.vals/self.axis.widths
            err = self.err/self.axis.widths
        else:
            vals = self.vals
            err = self.err
        pad.errorbar(x=self.axis.centers[mask], xerr=self.axis.widths[mask]/2,
                     y=self.draw_sc*vals[mask], ecolor=self.color,
                     yerr=self.draw_sc*err[mask], fmt='o',
                     color=self.color, barsabove=True, label=self.get_name(),
                     markersize=4, **kwargs)

    def plot_hist(self, pad, normed=False, **kwargs):
        if not self or pad is None:
            return
        mask = self.vals > 0.
        if normed:
            vals = self.vals/self.axis.widths
            err = self.err/self.axis.widths
        else:
            vals = self.vals
            err = self.err
        pad.hist(x=self.axis.centers[mask], weights=self.draw_sc*self.vals[mask], bins=self.axis.edges,
                 label=self.get_name(), histtype="step", linewidth=3, color=self.color, **kwargs)

    def plot_2d(self, pad, noText=False, **kwargs):
        import matplotlib.patheffects as path_effects
        if not self or pad is None:
            return

        xx = np.tile(self.axes[0].edges, (len(self.axes[1])+1, 1))
        yy = np.tile(self.axes[1].edges, (len(self.axes[0])+1, 1)).T
        vals = self.vals.T
        vals_masked = np.ma.masked_where(vals < 1e-5, vals)
        color_plot = pad.pcolormesh(xx, yy, vals_masked, shading='flat', **kwargs)
        if noText:
            return color_plot

        xstart, xend = self.get_xrange()
        min_size = (xend-xstart)/9
        min_ysize = (self.axes[1].edges[-1]-self.axes[1].edges[0])/14

        for j, y in enumerate(self.axes[1].centers):
            offset = False
            for i, x in enumerate(self.axes[0].centers):
                ha = 'center'
                if offset:
                    offset = False
                elif self.axis.widths[i-1] < min_size and self.axis.widths[i] < min_size:
                    offset = True

                ytot = y - offset*min_ysize
                if self.vals[i,j] < 1e-5:
                    continue
                val_str = f'{self.vals[i,j]:.3f}\n$\pm${self.err[i,j]:.3f}'
                text = pad.text(x, ytot, val_str, fontsize='x-small', ha=ha, va='center')
        return color_plot


    def plot_band(self, pad, normed=False, asymm=False, **kwargs):
        if not self or pad is None:
            return
        if normed:
            vals = self.vals/self.axis.widths
            err = self.err/self.axis.widths
        else:
            vals = self.vals
            err = self.err
        bottom = vals - err
        pad.hist(weights=2*err, x=self.axis.centers,
                 bins=self.axis.edges, bottom=bottom,
                 histtype='stepfilled', color=self.color,
                 align='mid', stacked=True, hatch='//',
                 alpha=0.4, label=self.get_name(), **kwargs)

    def plot_shape(self, pad, normed=False, **kwargs):
        if not self or pad is None:
            return
        pad.hist(x=self.axis.centers, weights=self.draw_sc*self.vals, bins=self.axis.edges,
                 label=self.get_name(), histtype="stepfilled", linewidth=1.5, density=True, alpha=0.3,
                 hatch="///", color=self.color, edgecolor=self.darkenColor())

    def get_int_err(self, sqrt_err=False, roundDigit=2):
        tot = self.hist.sum()
        err = np.sqrt(tot.variance) if sqrt_err else tot.variance
        return np.round(np.array([tot.value, err]), roundDigit)
