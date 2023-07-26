#!/usr/bin/env python3
import numpy as np
import prettytable

class Card_Maker:
    def __init__(self, path, year, cr, plot_groups, variable):
        if 'data' in plot_groups:
            plot_groups.remove('data')
        self.variable = variable
        self.path = path
        self.year = year
        self.cr = cr
        self.channels = np.array([cr])
        self.plot_groups = np.array(plot_groups)
        self.nChans = len(self.channels)
        self.nGroups = len(self.plot_groups)

    def __enter__(self):
        self.f = open(f"{self.path}/{self.variable}_{self.year}_{self.cr}_card.txt", 'w')
        return self

    def __exit__(self, type, value, traceback):
        self.f.close()

    def tab_list(self, inlist):
        return '\t'+'\t'.join(inlist.astype(str))

    def write(self, line):
        self.f.write(line)
        self.f.write("\n")

    def end_section(self):
        self.write("-"*50)

    def write_systematics(self, syst_list):
        # Specify numbers of groups
        self.write(f"imax {self.nChans}  number of channels")
        self.write(f"jmax {self.nGroups - 1}  number of backgrounds plus signals minus 1")
        self.write(f"kmax {len(syst_list)} number of nuisance parameters (sources of systematical uncertainties)")
        self.end_section()

        # Specify shape locations
        self.write(f"shapes * * {self.variable}_{self.year}_$CHANNEL.root $PROCESS $SYSTEMATIC/$PROCESS")
        self.end_section()

        # Specify channels and number of events
        self.write("bin" + self.tab_list(self.channels))
        self.write("observation" + self.tab_list(-1*np.ones(self.nChans)))
        self.end_section()

        table = prettytable.PrettyTable()
        table.header=False
        table.border=False
        table.left_padding_width=0

        # Specify channel and plot names with MC counts
        table.add_column('', ['bin', "process", "process", "rate"])
        for chan in self.channels:
            for p in np.arange(self.nGroups):
                table.add_column("", [chan, self.plot_groups[p], p, -1])
        self.write(table.get_string(align='l'))
        self.end_section()

        # Specify systematics
        table.clear()
        for syst in syst_list:
            table.add_row(syst.output(self.plot_groups, [self.year]))
        table.align='l'
        self.write(table.get_string(align='l'))
        self.write("syst_error group = " + " ".join([syst.name for syst in syst_list]))
        self.write("* autoMCStats 0")

    def add_rateParam(self, group):
        self.write(f'rate_{group} rateParam * {group} 1.0')
