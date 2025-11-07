import numpy as np
from copy import copy
from matplotlib import colors as clr
import boost_histogram as bh
from boost_histogram.accumulators import WeightedSum as bh_weights

colors = {
    'blue': '#3f90da',
    'b': '#3f90da',
    'yellow': '#ffa90e',
    'y': '#ffa90e',
    'red': '#bd1f01',
    'r': '#bd1f01',
    'grey': '#94a4a2',
    'gray': '#94a4a2',
    'purple': '#832db6',
    'brown': '#a96b59',
    'orange': '#e76300',
    'tan': '#b9ac70',
    'darkgrey': '#717581',
    'skyblue': '#92dadd',
    'black': '#000000',
    'k': '#000000',
}

import warnings
warnings.simplefilter("ignore", UserWarning)
from scipy.stats import beta

class Histogram(bh.Histogram):
    def __init__(self, *args, **kwargs):
        if len(args) == 0:
            axes = (bh.axis.Regular(1, 0, 1),)
        elif isinstance(args[0], bh.Histogram):
            axes = args[0].axes
        elif not isinstance(args[0], bh.axis.Axis):
            axes = (bh.axis.Regular(*args),)
        else:
            axes = args
        super().__init__(*axes, storage=bh.storage.Weight())

        if len(args) > 0 and isinstance(args[0], bh.Histogram):
            self.values()[:] = args[0].values()
            self.variances()[:] = args[0].variances()

        self.breakdown = dict()
        color = kwargs.get('color', 'k')
        self.color = color if color[0] == '#' else colors[color]
        self.plot_label = kwargs.get('name', "")
        self.axis_name = kwargs.get('axis_name', "")
        self.draw_scale = 1
        self.syst_error = None

    def __bool__(self):
        return not self.empty()

    def __sub__(self, right):
        return super().__add__(-1*right)

    def __truediv__(self, right):
        if not isinstance(right, Histogram):
            return super().__mul__(1/right)
        # error propagation for independent a, b:
        # c = a / b: var(c) = (var(a)/a^2 + var(b)/b^2) c^2
        output = self.copy(deep=False)
        variance = self.variances(flow=True)/(self.values(flow=True)**2+1e-10)
        variance += right.variances(flow=True)/(right.values(flow=True)**2+1e-10)
        output.values(flow=True)[:] /= (right.values(flow=True)+1e-10)
        output.variances(flow=True)[:] = variance*output.values(flow=True)**2
        return output

    def __getattr__(self, attr):
        if attr == 'axis':
            return self.axes[0]
        elif attr == 'vals':
            return self.move_overflow(copy(self.values(flow=True)))
        elif attr == 'sumw2':
            return self.move_overflow(copy(self.variances(flow=True)))
        elif attr == 'err':
            return np.sqrt(self.move_overflow(copy(self.variances(flow=True))))
        elif attr == 'syst_err2':
            return self.move_overflow(copy(self.syst_error))
        elif attr == 'total_err':
            stats = self.sumw2
            syst = self.syst_err2
            return np.sqrt(stats+syst)
        else:
            return getattr(super(), attr)
            # raise Exception(f'{attr} was not found')

    def get_self_ratio(self):
        output = Histogram(self.axis)
        output.values()[:] = np.where(abs(self.vals) > 1e-5, 1., 0.)
        output.variances()[:] = self.sumw2/(self.vals**2 + 1e-10)
        if self.syst_error is not None:
            output.syst_error = np.zeros(self.axis.extent)
            output.syst_error[1:-1] += self.syst_err2/(self.vals**2 + 1e-10)
        return output

    def add_syst(self, hist=None):
        if self.syst_error is None:
            self.syst_error = np.zeros(self.axis.extent)
        if hist is not None:
            self.syst_error += hist.values(flow=True)**2

    def set_data(self, val, sumw2=None):
        self.values()[:] = val
        if sumw2 is not None:
            self.variances()[:] = sumw2


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

        return_obj = Histogram(*top.axes)
        return_obj.values()[:] = eff
        return_obj.variances()[:] = error2
        return return_obj

    def move_overflow(self, vals):
        if len(self.axes) == 1:
            vals[-2] += vals[-1]
            vals[1] += vals[0]
            return vals[1:-1]
        else:
            vals[1] += vals[0]
            vals[-2] += vals[-1]
            vals[1:-1,1] += vals[1:-1,0]
            vals[1:-1,-2] += vals[1:-1,-1]
            return vals[1:-1,1:-1]

    def fill(self, *args, **kwargs):
        keep_overflow = kwargs.pop('flow', False)
        member = kwargs.pop('member', False)
        if member:
            wgt = np.array(kwargs['weight'])
            self.breakdown[member] = [bh_weights(sum(wgt), sum(wgt**2)), len(wgt)]
        super().fill(*args, **kwargs)
        if not keep_overflow:
            self.values(flow=True)[:] = np.pad(self.vals, 1)
            self.variances(flow=True)[:] = np.pad(self.sumw2, 1)
        else:
            self.values(flow=True)[:] = np.pad(self.values(flow=False)[:], 1)
            self.variances(flow=True)[:] = np.pad(self.variances()[:], 1)

    def get_name(self):
        name = f'${self.plot_label}$' if '\\' in self.plot_label else self.plot_label
        if self.draw_scale == 1:
            return name
        else:
            scale = self.draw_scale
            str_scale = str(scale) if isinstance(scale, int) else f'{scale:0.2f}'
            return f"{name} x {str_scale}"

    def set_plot_details(self, name, color):
        self.plot_label = name
        self.color = colors[color]

    def darkenColor(self, color=None):
        return colors['k']
        # if color is None:
        #     color = self.color
        # cvec = clr.to_rgb(color)
        # dark = 0.3
        # return [i - dark if i > dark else 0.0 for i in cvec]

    def get_xrange(self):
        return [self.axis.edges[0], self.axis.edges[-1]]

    def scale(self, scale, for_plot=False):
        if for_plot:
            self.draw_scale = scale
        else:
            self *= scale

    def integral(self, flow=True):
        return self.sum(flow=flow).value

    def plot_points(self, pad, normed=False, use_color=False, **kwargs):
        if not self or pad is None:
            return
        mask = self.vals > 0.
        width = self.axis.widths if normed else 1.
        vals = self.vals/width
        err = self.err/width
        color = self.color if use_color else 'k'
        xerr = None if kwargs.pop("data", False) else self.axis.widths[mask]/2
        a = pad.errorbar(x=self.axis.centers[mask], xerr=xerr,
                     y=self.draw_scale*vals[mask],
                     yerr=self.draw_scale*err[mask], fmt='o',
                     color=color, label=self.get_name(),
                     markersize=4, **kwargs)

    def plot_hist(self, pad, normed=False, **kwargs):
        if not self or pad is None:
            return
        mask = self.vals > 0.
        width = self.axis.widths if normed else 1.
        vals = self.vals/width
        err = self.err/width
        color = colors[kwargs.pop('color')] if 'color' in kwargs else self.color
        pad.hist(x=self.axis.centers[mask], weights=self.draw_scale*vals[mask], bins=self.axis.edges,
                 label=self.get_name(), histtype="step", linewidth=2, color=color, **kwargs)

    def plot_2d(self, pad, noText=False, **kwargs):
        if not self or pad is None:
            return

        xx = np.tile(self.axes[0].edges, (len(self.axes[1])+1, 1))
        yy = np.tile(self.axes[1].edges, (len(self.axes[0])+1, 1)).T
        vals = self.vals.T
        vals_masked = np.ma.masked_where(vals < 1e-5, vals)
        color_plot = pad.pcolormesh(xx, yy, vals_masked, shading='flat', rasterized=True,
                                    **kwargs)

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
                val_str = f'{self.vals[i,j]:.3f}\n' + fr'$\pm${self.err[i,j]:.3f}'
                text = pad.text(x, ytot, val_str, fontsize='x-small', ha=ha, va='center')
        return color_plot


    def plot_error_band(self, pad, color='r', hatch='+', systs=False, normed=False, **kwargs):
        if not self or pad is None:
            return
        width = self.axis.widths if normed else 1.
        if systs:
            err = self.total_err/width
            name = r"$Stats \oplus Systs$"
        else:
            err = self.err/width
            name = "Stats"
        bottom = self.vals/width - err
        pad.hist(weights=2*err, x=self.axis.centers, bins=self.axis.edges,
                 bottom=bottom, histtype='stepfilled', color=color, align='mid',
                 stacked=True, hatch=hatch, alpha=0.3, label=name, **kwargs)

    def plot_shape(self, pad, normed=False, **kwargs):
        if not self or pad is None:
            return
        pad.hist(x=self.axis.centers, weights=self.draw_scale*self.vals*self.axis.widths, bins=self.axis.edges,
                 label=self.get_name(), histtype="stepfilled", linewidth=1.5, density=True, alpha=0.3,
                 hatch="///", color=self.color, edgecolor=self.darkenColor())

    def get_int_err(self, sqrt_err=False, roundDigit=2):
        tot = self.sum(flow=True)
        err = np.sqrt(tot.variance) if sqrt_err else tot.variance
        return np.round(np.array([tot.value, err]), roundDigit)
