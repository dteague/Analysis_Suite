from analysis_suite.commons.histogram import Histogram

import numpy as np

class Stack(Histogram):
    def __init__(self, bin_info, stack_by_int=False):
        super().__init__("", bin_info)
        self.stack = list()
        self.options = {"stacked": True, "histtype": "stepfilled"}
        self.stack_by_int = stack_by_int

    def __iadd__(self, right):
        if self.stack_by_int:
            idx = self._get_index(right.integral())
            self.stack.insert(idx, right)
        else:
            self.stack.append(right)
        return super().__iadd__(right)

    def _get_index(self, integral):
        if not self.stack:
            return 0
        else:
            return np.argmax(np.array([s.integral() for s in self.stack]) < integral)

    def __getitem__(self, item):
        for hist in self.stack:
            if hist.group == item:
                return hist

    def recalculate_stack(self):
        vals = np.zeros(len(self.axis))
        sumw2 = np.zeros(len(self.axis))
        for hist in self.stack:
            vals += hist.vals
            sumw2 += hist.sumw2
        self.hist.view().value = vals
        self.hist.view().variance = sumw2


    def plot_stack(self, pad, **kwargs):
        if not self.stack:
            return
        n, bins, patches = pad.hist(
            weights=np.array([h.vals for h in self.stack]).T, bins=self.axis.edges,
            x=np.tile(self.axis.centers, (len(self.stack), 1)).T,
            color=[h.color for h in self.stack],
            label=[h.name for h in self.stack],
            **self.options, **kwargs
        )

        # Apply patch to edge colors
        edgecolors = [self.darkenColor(h.color) for h in self.stack]
        for p, ec in zip(patches, edgecolors):
            if isinstance(p, list):
                p[0].set_ec(ec)
            else:
                p.set(ec=ec)
