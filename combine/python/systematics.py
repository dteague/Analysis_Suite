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
    syst_type: str
    corr : bool = False
    dan_name: str = None
    systs: dict = field(default_factory=list)

    def __post_init__(self):
        self.systs = {year: dict() for year in all_eras}
        self.dan_name = self.name

    def output(self, group_list, year):
        name = self.get_name(year)
        systs = self.systs[year]
        line = [name, self.syst_type]
        for group in group_list:
            if group in systs:
                line.append(systs[group])
            elif "all" in systs:
                line.append(systs["all"])
            else:
                line.append('-')
        return line

    def good_syst(self, group_list, year):
        systs = self.systs[year]
        for group in group_list:
            if group in systs or "all" in systs:
                return True
        return False

    def add(self, syst, groups="all", year=all_eras):
        syst_dict = {group: syst for group in groups} if isinstance(groups, list) else {groups: syst}
        year = [year] if isinstance(year, str) else year
        for yr in year:
            self.systs[yr].update(syst_dict)

        return self

    def get_name(self, year):
        name = self.dan_name
        if not self.corr:
            name += year_convert[year]
        return name

    def dname(self, name):
        self.dan_name = name
        return self
