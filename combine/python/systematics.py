#!/usr/bin/env python3
from dataclasses import dataclass, field
from analysis_suite.commons.constants import all_eras

year_convert = {
    "2016pre": "16APV",
    "2016post": "16",
    "2017": "17",
    "2018": "18",
}



@dataclass
class Systematic:
    name: str
    syst_type: str = 'shape'
    corr : bool = False
    dan_name: str = None
    systs: dict = field(default_factory=dict)
    default_groups: list = None

    def __post_init__(self):
        self.systs = {# year: dict() for year in all_eras
                      }
        self.dan_name = self.name

    def output(self, group_list, year, chan='all'):
        if chan in self.systs[year]:
             syst = self.systs[year][chan]
        elif 'all' in self.systs[year]:
             syst = self.systs[year]['all']
        else:
            return None

        line = [self.get_name(year), self.syst_type]
        for group in group_list:
            if any([g in self.systs[year]['groups'] for g in [group, 'all']]):
                line.append(syst)
            else:
                line.append('-')
        return line

    def good_syst(self, group_list, year):
        if isinstance(group_list, str):
            group_list = [group_list]
        if year not in self.systs:
            return False
        systs = self.systs[year]['groups']
        for group in group_list:
            if group in systs or "all" in systs:
                return True
        return False

    def add(self, syst=1, groups=None, year=all_eras, chan='all'):
        groups = groups if groups is not None else Systematic.default_groups
        year = [year] if isinstance(year, str) else year
        chans = [chan] if isinstance(chan, str) else chan

        for yr in year:
            if yr not in self.systs:
                self.systs[yr] = {'groups': groups}
            for chan in chans:
                self.systs[yr][chan] = syst
        return self

    def get_name(self, year):
        name = self.dan_name
        if "JEC" in name or "JER" in name:
            if not self.corr:
                name += year[:3]
            name += "LOWESS"
        if not self.corr:
            name += year_convert[year]
        return name

    def dname(self, name):
        self.dan_name = name
        return self
