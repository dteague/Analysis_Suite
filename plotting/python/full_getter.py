#!/usr/bin/env python3
from .hist_getter import HistGetter, GraphInfo
from analysis_suite.commons.constants import all_eras

class YearGetter:
    def __init__(self, ntuple_info, years, **kwargs):
        self.factories = []
        if years == "all":
            years = all_eras
        elif isinstance(years, str):
            years = [years]
        for year in years:
            self.factories.append(HistGetter(ntuple_info, year, **kwargs))

    def run_all_years(self, key, *args, **kwargs):
        output = None
        for factory in self.factories:
            year_output = getattr(factory, key)(*args, **kwargs)
            if year_output is None:
                continue
            elif output is None:
                output = year_output
            else:
                def recursive_add(left, right):
                    if isinstance(left, dict):
                        for k, v in left.items():
                            if k not in right:
                                return
                            recursive_add(v, right[k])
                        for k in [k for k in right.keys() if k not in left]:
                            left[k] = right[k]
                    else:
                        left += right
                recursive_add(output, year_output)
        if output is not None:
            return output

    def __getattr__(self, key):
        return lambda *args, **kwargs: self.run_all_years(key, *args, **kwargs)
