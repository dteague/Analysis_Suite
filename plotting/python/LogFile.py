from prettytable import PrettyTable
import prettyprinter as pp
import math
import numpy as np
from copy import copy

from analysis_suite.commons.configs import asymptotic_sig, sig_fig

pp.install_extras()


BKG = 0
SIGNAL = 1
DATA = 2
TOTAL = 3

class LogFile:
    """Wrapper for Logfile for a plot

    Attributes
    ----------
    plotTable : PrettyTable
        PrettyTable holding all the histogram information
    output_name : string
        Name/path of logfile to be created
    analysis : string
        Analysis running over
    selection : string
        Selection of analysis running over
    lumi : float
        Luminosity of this run (in ipb, but converted to ifb in class)
    hists : list of lists
        List of lists with Integral and Error^2 indexed by the constants 
        BKG, SIGNAL, DATA, TOTAL
    callTime : string
        String of the time this script was called (any format)
    command : string
        Command used to start this script
    name : string
        Name of the histogram
    
    """
    callTime = ""
    command = ""
    def __init__(self, name, lumi, graph_obj=None):
        self.plotTable = PrettyTable(["Plot Group", "Weighted Events", "Error"])
        self.breakTable = PrettyTable(["Plot Group", "Sample", "Weighted Events", "Error"])
        self.output_name = f'{name}_info.log'
        self.lumi = lumi
        self.hists = [np.array([0., 0.]) for i in range(4)] 
        self.graph_obj = graph_obj
        self.name = name

    def add_mc(self, stacker):
        """Add background data to this class

        Parameters
        ----------
        drawOrder : list of tuples (string, GenericHist)
            List of all background hists with their names
        """
        for hist in stacker.stack:
            integral, error = hist.get_int_err(True)
            self.plotTable.add_row([hist.name, integral, error])
            self.hists[BKG] += hist.get_int_err()
        self.hists[TOTAL] += self.hists[BKG]

    def add_signal(self, signal):
        """Add signal data to this class

        Parameters
        ----------
        signal : GenericHist
            Histogram of signal information
        groupName : string
            Name used to label the singal
        """
        self.hists[SIGNAL] += signal.get_int_err()
        integral, error = signal.get_int_err(True)
        self.plotTable.add_row([signal.name, integral, error])
        self.hists[TOTAL] += self.hists[SIGNAL]

    @staticmethod
    def add_metainfo(callTime, command):
        """Set specific metadata for output file

        Parameters
        ----------
        callTime : string
            Time script was call (not formated, must be done before here)
        command : string
            Full commandline string use for this run

        """
        LogFile.callTime = callTime
        LogFile.command = command

    def get_sqrt_err(self, idx):
        """Grab the Integral and error on that information

        Parameters
        ----------
        idx : int
            int pointed to self.hists list
        """
        hist = copy(self.hists[idx])
        hist[1] = np.sqrt(hist[1])
        return hist

    def write_out(self, path, isLatex=False):
        """Write out all current information to the objects output file

        Parameters
        ----------
        isLatex : bool, optional
            Whether table should be written out in latex or org style table
        """
        with open(f'{path}/{self.output_name}', 'w') as out:
            out.write("<html><pre><code>\n")
            out.write('-' * 80 + '\n')
            out.write(f'Script called at {LogFile.callTime} \n')
            out.write(f'The command was: {LogFile.command} \n')
            out.write(f'The name of this Histogram is: {self.name} \n')
            out.write('-' * 80 + '\n')
            # out.write(f'Selection: {self.analysis,}/{self.selection}\n')
            out.write(f'Luminosity: {self.lumi:0.2f} fb^{{-1}}\n')
            if isLatex:
                out.write('\n' + self.plotTable.get_latex_string() + '\n'*2)
            else:
                out.write('\n' + self.plotTable.get_string() + '\n'*2)

            if self.hists[TOTAL].any():
                tot, err = self.get_sqrt_err(TOTAL)
                out.write(f"Total sum of Monte Carlo: {sig_fig(tot)} +/- {sig_fig(err)} \n")
            if self.hists[SIGNAL].any():
                tot, err = self.get_sqrt_err(BKG)
                out.write(f"Total sum of background Monte Carlo: {sig_fig(tot)} +/- {sig_fig(err)} \n")
                tot, err = self.get_sig_bkg_ratio()
                out.write(f"Ratio S/(S+B): {sig_fig(tot)} +/- {sig_fig(err)} \n")
                tot, err = self.get_likelihood()
                out.write(f"Ratio S/sqrt(B): {sig_fig(tot)} +/- {sig_fig(err)} \n")
            if self.hists[DATA].any():
                out.write(f'Number of events in data {self.hists[DATA][0]} \n')

            if True:
                out.write('\n')
                out.write(self.breakTable.get_latex_string() if isLatex else self.breakTable.get_string())
                out.write('\n'*2)

            if self.graph_obj is not None:
                out.write("\n")
                pp.pprint(self.graph_obj, stream=out)

            out.write("</code></pre></html>\n")


    def get_sig_bkg_ratio(self):
        """Get S/B with its error

        Returns
        -------
        tuple
            tuple S/B and its error
        """
        sig, sigErr = self.hists[SIGNAL]
        bkg, bkgErr = self.hists[BKG]
        s_b = sig/bkg
        s_b_err = s_b * math.sqrt(sigErr/sig**2 + bkgErr/ bkg**2)
        return (s_b, s_b_err)

    def get_likelihood(self):
        """Get Figure of merit

        Returns
        -------
        tuple
            tuple Figure of Merit (S/sqrt(S+B)) and its error
        """
        sig, sigErr = self.hists[SIGNAL]
        bkg, bkgErr = self.hists[BKG]
        s_sqrtb = asymptotic_sig(sig, bkg)
        s_sqrtb_err = s_sqrtb * math.sqrt(sigErr / sig**2 + 0.25*bkgErr/bkg**2)
        return (s_sqrtb, s_sqrtb_err)


    def add_breakdown(self, group, break_dict, roundDigit=2):
        sorted_keys = sorted(break_dict.keys(), key=lambda key: break_dict[key].value, reverse=True)
        for sample in sorted_keys:
            events = round(break_dict[sample].value, roundDigit)
            err = round(math.sqrt(break_dict[sample].variance), roundDigit)
            self.breakTable.add_row([group, sample, events, err])
            group = "" # So the group only shows up once in the table
