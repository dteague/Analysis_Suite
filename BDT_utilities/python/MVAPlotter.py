"""
.. module:: MVAPlotter
   :synopsis: Takes MvaMaker output and creates graphs
.. moduleauthor:: Dylan Teague
"""
import sys
from os.path import isfile
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
import xgboost as xgb
import uproot4 as uproot
# import shap


class MVAPlotter(object):
    """Wrapper for a bunch of plotting and outputing functions

    Args:
      groups(list): List of group names
      do_show(bool): Bool for if matplotlib should show its graphs
      save_dir(string): Directory name to save all plots
      lumi(float): Luminosity (in ipb)
      color_dict(dict: dict): Dictionary for setting line color for groups
      work_set(pandas.DataFrame): DataFrame with event/variable data
    """
    def __init__(self, workdir, groups, lumi=140000, is_train=False):
        self.groups = list(groups.keys())
        self.do_show = False
        self.save_dir = workdir
        self.lumi = lumi/1000
        colors = ['#f15854', '#faa43a', '#60bd68', '#5da5da']
        #  ['#CC0000', '#99FF00', '#FFCC00', '#3333FF']
        #  ['#462066', '#FFB85F', '#FF7A5A', '#00AAA0']

        self.color_dict = {grp: col for grp, col in zip(self.groups, colors)}
        self.color_dict["All Backgrounds"] = colors[len(self.groups)]
        self.work_set = pd.DataFrame()
        work_file, other_file = ("train.root", "test.root") \
            if is_train else  ("test.root", "train.root")
        
        with uproot.open(f'"{workdir}/{work_file}"') as f:
            for group, samples in groups.items():
                for sample in [s for s in samples if s in f]:
                    temp = f[sample].arrays(library="pd")
                    temp["classID"] = self.groups.index(group)
                    self.work_set = pd.concat((self.work_set, temp))

        other = 0
        with uproot.open(f'{workdir}/{other_file}') as f:
            column = self.work_set.columns[0]
            for sample in sum(list(groups.values()), []):
                if sample not in f:
                    continue
                other += len(f[sample].arrays([column]))

        print(self.work_set.columns)
        print(other)
        
        self._ratio = (1. + 1.*other/len(self.work_set))
        self.work_set.scale_factor *= self._ratio
        self.work_set["final_factor"] = self.work_set.scale_factor*lumi
        

    def __len__(self):
        return len(self.groups)

    def set_show(self, do_show):
        """**Set do_show to a value**

        Args:
          do_show(bool): If Matplotlib show show plots
        """
        self.do_show = do_show

    def add_variable(self, name, arr):
        """**Add variable to MVAPlotter object**

        If one needs to make a unique, composite variable, this allows
        the user to add it to the DataFrame

        Args:
          name(string): Name of variable for the DataFrame
          arr(numpy.ndarray): Array of variable values for the DataFrame
        """

        if len(arr) != len(self.work_set):
            print("bad!")
            sys.exit()
        self.work_set.insert(1, name, arr)

    def apply_cut_max_bdt(self, sample, is_max=True):
        """**Cut on which BDT is maximum**

        Cut on file by only keeping the events where the BDT with the
        maximum value matches the sample provided. There is a `is_max`
        flag to allow the inverse cut

        Args:
          sample(string): Sample BDT that must be maximized
          is_max(bool, optional): whether to keep only maxBDT or all BUT, defaults to True
        """
        bdt = list()
        for group in self.groups:
            if group == "Signal":
                bdt.append(self.work_set[f'BDT.{group}'])
            else:
                bdt.append(1-self.work_set[f'BDT.{group}'])

        max_bdt = np.argmax(bdt, axis=0)
        if is_max:
            self.work_set = self.work_set[max_bdt == self.groups.index(sample)]
        else:
            self.work_set = self.work_set[max_bdt != self.groups.index(sample)]

    def apply_cut(self, cut):
        """**Cut DataFrame using Root Style cut string**

        Args:
          cut(string): Cut string
        """
        if cut.find("<") != -1:
            tmp = cut.split("<")
            self.work_set = self.work_set[self.work_set[tmp[0]] < float(tmp[1])]
        elif cut.find(">") != -1:
            tmp = cut.split(">")
            self.work_set = self.work_set[self.work_set[tmp[0]] > float(tmp[1])]
        elif cut.find("==") != -1:
            tmp = cut.split("==")
            self.work_set = self.work_set[self.work_set[tmp[0]] == float(tmp[1])]
        else:
            print("Problem!")
            sys.exit()

    def get_sample(self, groups=None):
        """**Get DataFrame matching the groups**

        .. note:: groups set to None (the default) means the whole
                  DataFrame is returned

        Args:
          groups(list, optional): list of groups to include, defaults to None

        Returns:
          pandas.DataFrame: DataFrame including only groups specified (None means all
          groups)
        """
        # Not list, just single instance
        if isinstance(groups, str):
            return self.work_set[self.work_set["classID"] == self.groups.index(groups)]
        elif not isinstance(groups, (list, np.ndarray)):
            return self.work_set

        # if is a list
        final_set = pd.DataFrame(columns=self.work_set.columns)
        for group in groups:
            tmp_set = self.work_set[self.work_set["classID"] == self.groups.index(group)]
            final_set = pd.concat((final_set, tmp_set))
        return final_set

    def get_variable(self, var, groups=None):
        """**Get numpy array of variable**

        .. note:: This function uses `get_sample` so groups=None means all groups included

        Args:
          var(string): Variable returned
          groups(list, optional): list of groups to include, defaults to None

        Returns:
          numpy.ndarray: Array of the variable (under the group constraint)

        """
        return self.get_sample(groups)[var]

    def get_hist(self, var, bins, groups=None):
        """**Get Histogram of variable**

        Get numpy histogram of a variable

        Args:
          var(string): Variable to put in the histogram
          bins(numpy.ndarray): numpy array of bins
          groups(list, optional): list of groups to include, defaults to None (means all)

        Returns:
          tuple - ndarray of histogram, ndarray of bins: Numpy histogram of the variable

        """
        final_set = self.get_sample(groups)
        return np.histogram(final_set[var], bins=bins, weights=final_set.final_factor)[0]

    def get_hist_err2(self, var, bins, groups=None):
        """**Get Error squared histogram for variable**

        Get numpy histogram of variable squared (poisson err if sqrt)

        Args:
          var(string): Variable to put in the histogram
          bins(numpy.ndarray): numpy array of bins
          groups(list, optional): list of groups to include, defaults to None (means all)

        Returns:
          tuple - ndarray of histogram, ndarray of bins: Numpy histogram of the variable

        """
        final_set = self.get_sample(groups)
        return np.histogram(final_set[var], bins=bins,
                            weights=final_set.final_factor**2)[0]

    def get_hist_2d(self, var1, var2, bins, groups=None):
        """**Get numpy 2D histogram for 2 variables**

        Args:
          var1(string): Variable 1 for the x-axis
          var2(string): Variable 2 for the y-axis
          bins(numpy.ndarray): numpy array of bins
          groups(list, optional): list of groups to include, defaults to None (means all)

        Returns:
          tuple - ndarray of histogram, ndarray of bins: Numpy histogram of the variable

        """
        final_set = self.get_sample(groups)
        return np.histogram2d(x=final_set[var1], y=final_set[var2], bins=bins,
                              weights=final_set.final_factor)[0]

    def get_hist_err2_2d(self, var1, var2, bins, groups=None):
        """**Get numpy 2d histogram for 2 variables squared**

        Poisson error can be found by sqrting this histogram

        Args:
          var1(string): Variable 1 for the x-axis
          var2(string): Variable 2 for the y-axis
          bins(numpy.ndarray): numpy array of bins
          groups(list, optional): list of groups to include, defaults to None (means all)

        Returns:
          tuple - ndarray of histogram, ndarray of bins: Numpy histogram of the variable

        """
        final_set = self.get_sample(groups)
        return np.histogram2d(x=final_set[var1], y=final_set[var2], bins=bins,
                              weights=final_set.final_factor**2)[0]

    def plot_func(self, sig, bkg, var, bins, extra_name="", scale=True):
        """**Plot arbitrary variable**

        Args:
          sig(string): Primary group to be plotted
          bkg(string): Other groups to be plotted against sig
          var(string): Variable graphed
          bins(numpy.ndarray): numpy array of bins
          extra_name(string, optional): name to append to filename, defaults to ""
          scale(bool, optional): whether to scale sig up to bkg, defaults to True
        """
        sig_hist = self.get_hist(var, bins, sig)
        bkg_hist = self.get_hist(var, bins, bkg)
        #scale = self._find_scale(max(sig_hist), max(bkg_hist)) if scale else 1.
        scale = sum(bkg_hist)/sum(sig_hist)
        bkg_name = "All Backgrounds" if len(bkg) > 1 else bkg[0]
        # sig_name = sig if scale == 1 else "{} x {}".format(sig, scale)
        sig_name = sig
        if extra_name:
            extra_name = f'_{extra_name}'

        # Make plot
        fig, ax = plt.subplots()
        ax.hist(x=bins[:-1], weights=sig_hist*scale, bins=bins, label=sig_name,
                histtype="step", linewidth=1.5, color=self.color_dict[sig])
        ax.hist(x=bins[:-1], weights=bkg_hist, bins=bins, label=bkg_name,
                histtype="step", linewidth=1.5, color=self.color_dict[bkg_name])
        ax.legend(loc='upper left')
        ax.set_xlabel(var)
        ax.set_ylabel("Events/bin")
        ax.set_title(f'Lumi = {self.lumi} ifb')
        fig.tight_layout()
        plt.savefig(f'{self.save_dir}/{var}{extra_name}.png')
        if self.do_show:
            plt.show()
        plt.close()

    def plot_func_2d(self, samples, var1, var2, bins1, bins2, name,
                     lines=None):
        """**plot 2 arbitrary variable**

        .. note:: This can't plot 2 groups at once!

        Args:
          samples: type samples:
          var1(string): Variable 1 for the x-axis
          var2(string): Variable 2 for the y-axis
          bins1(numpy.ndarray): numpy array of bins for the x-axis
          bins2(numpy.ndarray): numpy array of bins for the y-axis
          name(string): name used for the plot
          lines(list): Coordinate for point where lines will be draw, defaults to None
        """
        grp = self.get_sample(samples)
        fig, ax = plt.subplots()
        hist2d = ax.hist2d(grp[var1], grp[var2], [bins1, bins2],
                           weights=grp.final_factor, cmap=plt.cm.jet)
        if lines is not None:
            ax.plot([lines[0], lines[0]], [bins2[0],bins2[-1]], 'r-')
            ax.plot([bins1[0],bins1[-1]],  [lines[1], lines[1]], 'r-')

        fig.colorbar(hist2d[-1])
        ax.set_xlabel(var1)
        ax.set_ylabel(var2)
        ax.set_title(f'Lumi = {self.lumi} ifb')
        fig.tight_layout()
        plt.savefig(f'{self.save_dir}/2D_{name}.png')
        if self.do_show:
            plt.show()
        plt.close()

    def get_fom(self, sig, bkg, var, bins, sb_denom=True, reverse=False):
        """**Return FoM histogram**

        Args:
          sig(string): Primary group to be plotted
          bkg(string): Other groups to be plotted against sig
          var(string): Variable used to calculate Figure of Merit
          bins(numpy.ndarray): numpy array of bins
          sb_denom(bool, optional): If False, uses S/sqrt(B) instead of S/sqrt(S+B), defaults to True
          reverse(bool, optional): if True, calculates FoM using LESS than not greater than, defaults to False

        Returns:
          list: histogram of the Figure of merit (only weights)

        """
        drt = 1 if not reverse else -1
        sig_hist = self.get_hist(var, bins, sig)[::drt]
        bkg_hist = self.get_hist(var, bins, bkg)[::drt]
        n_sig = [np.sum(sig_hist[i:]) for i in range(len(bins))]
        if sb_denom:
            n_tot = [np.sum(bkg_hist[i:])+np.sum(sig_hist[i:])
                     for i in range(len(bins))]
        else:
            n_tot = [np.sum(bkg_hist[i:])+np.sum(sig_hist[i:])
                     for i in range(len(bins))]
        return [s/math.sqrt(t) if t > 0 else 0 for s, t in zip(n_sig, n_tot)][::drt]

    def plot_fom(self, sig, bkg, var, bins, extra_name="", sb_denom=True,
                 reverse=False):
        """**Plots Figure of Merits for variable scan**

        Plot figure of merit by creating a "signal region" by only
        keeping events with a value of var greater than what value is
        being scanned. One can use a cut that keeps events LESS than
        the cut for var using the reverse flag.

        Args:
          sig(string): Primary group to be plotted
          bkg(string): Other groups to be plotted against sig
          var(string): Variable graphed
          bins(numpy.ndarray): numpy array of bins
          extra_name(string, optional): name to append to filename, defaults to ""
          sb_denom(bool, optional): If False, uses S/sqrt(B) instead of S/sqrt(S+B), defaults to True
          reverse(bool, optional): if True, calculates FoM using LESS than not greater than, defaults to False

        Returns:
          list: List of (FOM, bin) for the max FOM
        """
        if extra_name:
            extra_name = f'_{extra_name}'
        fom = self.get_fom(sig, bkg, var, bins, sb_denom, reverse)
        fom_maxbin = bins[fom.index(max(fom))]

        fig, ax = plt.subplots()

        fom_plot = ax.plot(bins, fom, label= f'$S/\sqrt{S+B}={max(fom):.3f}$\n cut={fom_maxbin:.2f}',
                           color='#fd6a02')
        ax.plot(np.linspace(bins[0], bins[-1], 5), [max(fom)]*5,
                linestyle=':', color=fom_plot[-1].get_color())
        ax.set_xlabel("BDT value", horizontalalignment='right', x=1.0)
        ax.set_ylabel("A.U.", horizontalalignment='right', y=1.0)

        ax2 = ax.twinx()
        sig_hist = self.get_hist(var, bins, sig)
        bkg_hist = self.get_hist(var, bins, bkg)
        bkg_name = "All Backgrounds" if len(bkg) > 1 else bkg[0]
        sig_plot = ax2.hist(x=bins[:-1], weights=sig_hist, bins=bins, histtype="step",
                 linewidth=1.5, color=self.color_dict[sig], label=sig,
                 density=True)
        bkg_plot = ax2.hist(x=bins[:-1], weights=bkg_hist, bins=bins, histtype="step",
                            linewidth=1.5, color=self.color_dict[bkg_name], density=True,
                            label=bkg_name)

        if reverse:
            ax.set_title("Reversed Cumulative Direction")

        lines, labels = ax2.get_legend_handles_labels()
        lines2, labels2 = ax.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='center left')
        fig.tight_layout()
        plt.savefig("%s/StoB%s.png" % (self.save_dir, extra_name))
        if self.do_show:
            plt.show()
        plt.close()
        return (fom, fom_maxbin)

    def plot_fom_2d(self, sig, var1, var2, bins1, bins2, extra_name=""):
        """**Plot Figure of Merit using a 2D scan**

        Figure of merit is calculated by cutting on var1 and var2,
        using the events with values greater than the cuts to define
        the "signal region".
        Returns info on max bin and max StoB

        Args:
          sig(string): Group to plot in the graph
          var1(string): Variable 1 for the x-axis
          var2(string): Variable 2 for the y-axis
          bins1(numpy.ndarray): numpy array of bins for the x-axis
          bins2(numpy.ndarray): numpy array of bins for the y-axis
          extra_name(string, optional): name to append to filename, defaults to ""

        Returns:
          list: List of (FOM, x value, y value) for the max FOM

        """
        if extra_name:
            extra_name = f'_{extra_name}'
        grp_id = self.groups.index(sig)

        zvals = []
        xbins = list(bins1[:-1])*(len(bins2)-1)
        ybins = np.array([[i]*(len(bins1)-1) for i in bins2[:-1]]).flatten()
        max_fom = (0, -1, -1)
        for valy in bins2[:-1]:
            # df_cut2 = self.work_set[self.work_set[var2] >= valy]
            df_cut2 = self.work_set[self.work_set[var2] <= valy] # only 4top
            for valx in bins1[:-1]:
                final = df_cut2[df_cut2[var1] >= valx]
                s = np.sum(final[final["classID"] == grp_id].final_factor)
                b = np.sum(final[final["classID"] != grp_id].final_factor)
                fom = s/math.sqrt(s+b) if b+s > 0 and s > 0 else 0
                if fom > max_fom[0]:
                    max_fom = (fom, valx, valy)
                zvals.append(fom)
        plt.hist2d(xbins, ybins, [bins1, bins2], weights=zvals,
                   cmap=plt.cm.jet)
        plt.plot([max_fom[1], max_fom[1]], [0, 1], 'r-')
        plt.plot([0, 1], [max_fom[2], max_fom[2]], 'r-')
        plt.colorbar()
        plt.xlabel(var1)
        plt.ylabel(var2)
        plt.savefig(f'{self.save_dir}/stob2D{extra_name}.png')
        plt.close()

        return max_fom

    def approx_likelihood(self, sig, bkg, var, bins):
        """**Get Basic max Likelihood significance**

        Args:
          sig(string): Primary group to be plotted
          bkg(string): Other groups to be plotted against sig
          var(string): Variable used to calculate
          bins(numpy.ndarray): numpy array of bins

        Returns:
          float: The Maximum Log Likelihood approximation (without errors)
        """
        sig_hist = self.get_hist(var, bins, sig)
        bkg_hist = self.get_hist(var, bins, bkg)
        term1, term2 = 0, 0
        for sig_val, bkg_val in zip(sig_hist, bkg_hist):
            if bkg_val <= 0 or sig_val <= 0:
                continue
            term1 += (sig_val+bkg_val)*math.log(1+sig_val/bkg_val)
            term2 += sig_val
        return math.sqrt(2*(term1 - term2))

    def make_roc(self, sig, bkg, var, extra_name=""):
        """Make and Save ROC curve for variable

        Args:
          sig(string): Primary group to be plotted
          bkg(string): Other groups to be plotted against sig
          var(string): Variable used to plot ROC curve (name of group)
          extra_name(string, optional): name to append to filename, defaults to ""
        """
        
        final_set = pd.concat((self.get_sample(sig), self.get_sample(bkg)))
        pred = final_set[var].array
        if extra_name:
            extra_name = f'_{extra_name}'


        truth = [1 if i == self.groups.index(sig) else 0
                 for i in final_set["classID"].array]
        fpr, tpr, _ = roc_curve(truth, pred)
        auc = roc_auc_score(truth, pred)

        # plot
        fig, ax = plt.subplots()
        ax.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
        ax.plot(np.linspace(0, 1, 5), np.linspace(0, 1, 5), linestyle=':')
        ax.legend()
        ax.set_xlabel("False Positive Rate", horizontalalignment='right', x=1.0)
        ax.set_ylabel("True Positive Rate", horizontalalignment='right', y=1.0)
        fig.tight_layout()
        plt.savefig(f'{self.save_dir}/roc_curve.BDT.{var}{extra_name}.png')
        if self.do_show:
            plt.show()
        plt.close()

    def write_out_root(self, filename, chan="SS"):
        """**Write DataFrame out to ROOT file compatible with VVPlotter**

        Args:
          filename(string): Name the ROOT file is saved as
          input_tree(string): input_tree filename used for getting sumweight graphs
          chan(string, optional): Folder plots are put into, defaults to "SS"
        """
        import uproot

        rm_groups = ["groupName", "classID", "finalWeight"]
        write_set = self.work_set.drop(columns=rm_groups)
        type_by_name = {name: "float32" if name[0] != "N" else "int"
                        for name in write_set.keys()}
        with uproot.recreate(f'{self.save_dir}/{filename}.root',
                             compression=uproot.ZLIB(4)) as outfile:
            for group in np.unique(self.work_set.groupName):
                tree_name = f'{group}_{chan}'
                try:
                    outfile[tree_name] = uproot.newtree(type_by_name)
                    outfile[tree_name].extend(
                        write_set[self.work_set.groupName == group].to_dict("list"))
                except:
                    continue

    def write_out(self, outdir):
        full_dir = f'{self.save_dir}/{outdir}'
        if not os.path.isdir(full_dir):
            os.mkdir(full_dir)
        for group in np.unique(self.work_set.groupName):
            outname = f'{full_dir}/{group}.parquet'
            self.work_set[self.work_set.groupName == group].to_parquet(outname, compression="gzip")

    def print_info(self, var, subgroups):
        """**Print out basic statistics information for a variable**

        The information for each sample in a group includes the mean,
        standard deviation, kurtosis, total weighted sum, total raw
        events.

        Args:
          var(string): Variable used to print statistics
          subgroups(list): list of all the samples for the information to be printed about
        """
        from scipy.stats import kurtosis
        info = list()
        for grp_name in subgroups:
            final_set = self.work_set[self.work_set.groupName == grp_name]
            variable = final_set[var]
            info.append([np.mean(variable), np.std(variable),
                         kurtosis(variable), np.sum(final_set.final_factor),
                         len(final_set), grp_name])

        info = sorted(info, reverse=True)
        print("| name      | events | raw events | mean+-std | kurtosis |")
        print("-"*50)
        for arr in info:
            print(f'|{arr[5]:10} | {arr[3]:.2f}| {arr[4]} | {arr[0]:.2f}+-{arr[1]:.2f} | {arr[2]:0.2f} |')
        print("-"*50)

    def plot_all_shapes(self, var, bins, extra_name=""):
        """**Plot all groups (normalized to 1) to compare shapes**

        Args:
          var(string): Variable graphed
          bins(numpy.ndarray): numpy array of bins
          extra_name(name: name to append to filename, defaults to "", optional): name to append to filename, defaults to ""
        """
        if extra_name:
            extra_name = f'_{extra_name}'

        fig, ax = plt.subplots()
        for group in self.groups:
            ax.hist(x=bins[:-1], weights=self.get_hist(var, bins, group),
                    bins=bins, label=group, histtype="step", linewidth=1.5,
                    density=True)
        ax.legend()
        ax.set_xlabel(var)
        ax.set_ylabel("A.U.")
        ax.set_title(f'Lumi = {self.lumi} ifb')
        fig.tight_layout()
        plt.savefig(f'{self.save_dir}/{var}{extra_name}.png')
        plt.close()

    def _find_scale(self, s, b):
        """**Find scale factor for signal to match rough size of background**

        Increases the scale for signal in a roughly logorithmic
        pattern (1, 5, 10 , 50, 100, ...) to find where the signal and
        background are roughly the same. This is used to make sure the
        plots are scaled well enough to see signal and background

        Args:
          s(float): Number of signal events
          b(float): Number of background events

        Returns:
          int: Scale to be used
        """
        scale = 1
        prev_s = 1
        while b//(scale*s) != 0:
            prev_s = scale
            if int(math.log10(scale)) == math.log10(scale):
                scale *= 5
            else:
                scale *= 2
        return prev_s
