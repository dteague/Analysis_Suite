#!/usr/bin/env python3
import sys
import datetime
import subprocess
import logging

import analysis_suite.commons.configs as config
import analysis_suite.commons.constants as constants
from analysis_suite.plotting.html_utils import writeHTML, get_plot_area, make_plot_paths
import analysis_suite.commons.user as user
import analysis_suite.data.plotInfo.plotInfo as plots_module

from .plotter import Plotter
from .LogFile import LogFile

def setup(cli_args):
    callTime = str(datetime.datetime.now())
    command = ' '.join(sys.argv)
    LogFile.add_metainfo(callTime, command)

    basePath = get_plot_area(cli_args.region, cli_args.name, cli_args.workdir)
    make_plot_paths(basePath)

    ginfo = config.get_ntuple_info(cli_args.region)
    plots = getattr(plots_module, cli_args.plots)
    if cli_args.hists != ['all']:
        plots = [graph for graph in plots if graph.name in cli_args.hists]

    argList = list()
    allSysts = ["Nominal"]
    for year in cli_args.years:
        for syst in allSysts:
            if cli_args.type == 'ntuple':
                filename = ntuple.get_filename(year=year, workdir=input_dir)
                filename = user.hdfs_area / f'workspace/{cli_args.region}/{year}'
            else:
                filename = cli_args.workdir / year / f'{cli_args.type}_{syst}_{cli_args.region}.root'
            outpath = basePath / year
            if syst != "Nominal":
                outpath = outpath / syst
            make_plot_paths(outpath)
            for plot in plots:
                argList.append((filename, outpath, plot, cli_args.signal, ginfo, year, syst, cli_args))
    return argList



def run(infile, outpath, graph, signalName, ginfo, year, syst, args):
    logging.info(f'Processing {graph.name} for year {year} and systematic {syst}')
    kwargs = {'ratio_bot': args.ratio_range[0], 'ratio_top': args.ratio_range[1]}

    inputs = config.get_inputs(args.workdir)

    groups = ginfo.setup_groups()
    bkg = list(filter(lambda x: x not in [signalName, 'data'] ,groups.keys()))

    blind = args.region == "Signal" and not args.unblind
    plotter = Plotter(infile, groups, year=year, ginfo=ginfo, scale=True)

    plotter.set_groups(sig=signalName, bkg=bkg, data='data')
    for run_group, scale_func in inputs.scale.items():
        if run_group in groups:
            plotter.scale(lambda vg: scale_func(vg, year), groups=run_group)
    plotter.fill_hists(graph)

    graph_name = f'{graph.name}'
    plotter.plot_from_graph(graph, outpath/'plots'/f'{graph.name}.png', **kwargs)
    # graph_name = f'{graph.name}_shape'
    # plotter.plot_shape(graph.name, outpath/'plots'/f'{graph_name}.png', **kwargs)

    subprocess.call('convert {0}.png {0}.pdf'.format(outpath/'plots'/graph_name), shell=True)

    # setup log file
    logger = LogFile(graph_name, constants.lumi[year], graph)
    for group, hist in plotter.get_hists(graph.name).items():
        logger.add_breakdown(group, hist.breakdown)
    logger.add_mc(plotter.make_stack(graph.name))
    if signalName:
        logger.add_signal(plotter.get_hists(graph.name, signalName))
    logger.write_out(outpath/'logs')


def cleanup(cli_args):
    basePath = get_plot_area(cli_args.region, cli_args.name, cli_args.workdir)
    # combined page
    writeHTML(basePath, cli_args.name, constants.years)
    for year in cli_args.years:
        systs = []
        # systs = [i[1] for i in get_files(cli_args, year) if i[1] != "Nominal"]
        yearAnalysis = f'{cli_args.name}/{year}'
        writeHTML(basePath / year, yearAnalysis, systs)
        for syst in systs:
            writeHTML(basePath / year / syst, f'{yearAnalysis}/{syst}')
    logging.critical(user.website+str(basePath.relative_to(user.www_area)))
