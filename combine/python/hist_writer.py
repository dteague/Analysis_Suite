#!/usr/bin/env python3
import uproot
from analysis_suite.commons.histogram import Histogram

class HistWriter:
    def __init__(self, outfile):
        self.outfile = outfile

    def __enter__(self):
        self.f = uproot.recreate(self.outfile)
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def add_syst(self, hists, ginfo, syst="Nominal", blind=True):
        syst = syst.replace("_up", "Up").replace("_down", "Down")
        for group, hist in hists.items():
            outgroup = ginfo.get_combine_name(group)
            if group == 'data':
                if syst == "Nominal" and not blind:
                    self.f['data_obs'] = hist.hist
            else:
                outname = outgroup if syst == "Nominal" else f'{syst}/{outgroup}'
                self.f[outname] = hist.hist

        if syst == "Nominal" and blind:
            data_obs = self.get_totmc(hists).hist
            data_obs.variances()[:] = data_obs.values()
            self.f['data_obs'] = data_obs

    def get_totmc(self, hists):
        first = next(iter(hists.values()))
        all_hist = Histogram("", first.axis)
        for group, hist in hists.items():
            if group == 'data':
                continue
            all_hist += hist
        return all_hist

    def close(self):
        self.f.close()


