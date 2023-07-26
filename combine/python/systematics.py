#!/usr/bin/env python3
from dataclasses import dataclass, field

@dataclass
class Systematic:
    name: str
    syst_type: str
    systs: dict = field(default_factory=list)

    def __post_init__(self):
        self.systs = {"2016pre": dict(), "2016post": dict(), "2017": dict(), "2018": dict()}

    def output(self, group_list, years):
        line = [self.name, self.syst_type]
        for year in years:
            for group in group_list:
                if group in self.systs[year]:
                    line.append(self.systs[year][group])
                elif "all" in self.systs[year]:
                    line.append(str(self.systs[year]["all"]))
                else:
                    line.append('-')
        return line

    def add(self, syst, groups="all", year=["2016pre", "2016post", "2017", "2018"]):
        if isinstance(year, list):
            for yr in year:
                if isinstance(groups, list):
                    self.systs[yr].update({group: syst for group in groups})
                else:
                    self.systs[yr][groups] = syst
        else:
            year = str(year)
            if isinstance(groups, list):
                self.systs[year].update({group: syst for group in groups})
            else:
                self.systs[year][groups] = syst

        return self
