#!/usr/bin/env python3

def plot_stack(self, name, outfile, chan=None, **kwargs):
    graph = self.graphs[name]
    signal = self.hists[name].get(self.sig, Histogram(""))
    data = Histogram("")
    if self.data in self.hists[name]:
        data = copy(self.hists[name].get(self.data, Histogram("")))
        data.color = 'k'
    stack = self.make_stack(name)
    plotter = ratio_plot if data else nonratio_plot

    axis_name = graph.axis_name if chan is None else graph.axis_name.format(chan)
    with plotter(self.workdir/outfile, axis_name, stack.get_xrange(), **kwargs) as ax:
        ratio = Histogram("Ratio", graph.bin_tuple, color="black")
        band = Histogram("Ratio", graph.bin_tuple, color="plum")
        error = Histogram("Stat Errors", graph.bin_tuple, color="plum")

        if data:
            ratio += data/stack
        band += stack/stack
        for hist in stack.stack:
            error += hist

        if data:
            pad, subpad = ax
        else:
            pad, subpad = ax, None

        #upper pad
        if self.scale_signal:
            # sig_scale = 250
            # rounder = 250
            # sig_scale = np.max(stack.vals)/np.max(signal.vals)
            # sig_scale = int(sig_scale//rounder*rounder)
            sig_scale = 750
            signal.scale(sig_scale, changeName=True, forPlot=True)

        stack.plot_stack(pad)
        signal.plot_points(pad)
        data.plot_points(pad)
        error.plot_band(pad)

        if (region := kwargs.get('region', False)) and data:
            (subpad if subpad else pad).text(graph.edges()[0], -0.7, region.format(chan))

        # ratio pad
        ratio.plot_points(subpad)
        band.plot_band(subpad)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hep.cms.label(ax=pad, lumi=self.lumi, data=data)
