#!/usr/bin/env python3
import numpy as np
import prettytable

class Card_Maker:
    def __init__(self, path, year, cr, signals, plot_groups, variable, nosyst=False):
        if 'data' in plot_groups:
            plot_groups.remove('data')
        self.variable = variable
        self.path = path
        self.year = year
        self.cr = cr
        self.channels = np.array([cr])
        self.signals = signals
        self.plot_groups = np.array(plot_groups)
        self.nosyst = nosyst

    def __enter__(self):
        region = self.cr + ("_nosyst" if self.nosyst else "")
        self.f = open(f"{self.path}/{self.variable}_{self.year}_{region}_card.txt", 'w')
        return self

    def __exit__(self, type, value, traceback):
        self.f.close()

    def _tab(self, arr, sep=" "):
        return sep.join(arr.astype(str))

    def write(self, line):
        self.f.write(line)
        self.f.write("\n")

    def end_section(self):
        self.write("-"*50)

    def write_preamble(self):
        # Specify numbers of groups
        self.write(f"imax *  number of channels")
        self.write(f"jmax *  number of backgrounds plus signals minus 1")
        self.write(f"kmax *  number of nuisance parameters (sources of systematical uncertainties)")
        self.end_section()

        # Specify shape locations
        self.write(f"shapes * * {self.variable}_{self.year}_$CHANNEL.root $PROCESS $SYSTEMATIC/$PROCESS")
        self.end_section()

        # Specify channels and number of events
        self.write("bin " + " ".join(self.channels))
        self.write("observation " + self._tab(np.full(len(self.channels), -1)))
        self.end_section()

        table = prettytable.PrettyTable()
        table.header=False
        table.border=False
        table.left_padding_width=0

        # Specify channel and plot names with MC counts
        table.add_column('', ['bin', "process", "process", "rate"])
        for chan in self.channels:
            for p in np.arange(0, -len(self.signals), -1):
                table.add_column("", [chan, self.signals[abs(p)], p, -1])
            for p in np.arange(len(self.plot_groups)):
                table.add_column("", [chan, self.plot_groups[p], p+1, -1])
        self.write(table.get_string(align='l'))
        self.end_section()


    def write_systematics(self, syst_list):
        all_groups = np.concatenate((self.signals, self.plot_groups))
        syst_list = list(filter(lambda x: x.good_syst(all_groups, self.year), syst_list))

        table = prettytable.PrettyTable()
        table.header=False
        table.border=False
        table.left_padding_width=0

        # Specify systematics
        for syst in syst_list:
            syst_row = syst.output(all_groups, self.year)
            if syst_row is not None:
                table.add_row(syst_row)
        table.align='l'
        self.write(table.get_string(align='l'))
        self.write("syst_error group = " + " ".join([syst.get_name(self.year) for syst in syst_list]))

    def add_rateParam(self, group):
        self.write(f'rate_{group} rateParam * {group} 1.0')

    def add_stats(self, autoStats=0):
        self.write(f"* autoMCStats {autoStats}")
