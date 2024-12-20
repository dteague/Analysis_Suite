#!/usr/bin/env python3
import sys
import datetime
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

    ntuple = config.get_ntuple(cli_args.region, obj=cli_args.region_type)

    argList = list()
    allSysts = ["Nominal"]
    for year in cli_args.years:
        for syst in allSysts:
            if cli_args.type == 'ntuple':
                filename = ntuple.get_filename(year=year)
            else:
                filename = cli_args.workdir / year / f'{cli_args.type}_{syst}_{cli_args.region}.root'
            outpath = basePath / year
            if syst != "Nominal":
                outpath = outpath / syst
            make_plot_paths(outpath)
            if cli_args.type == 'ntuple':
                argList.append((filename, outpath, None, cli_args.signal, year, syst, cli_args))
            else:
                plots = getattr(plots_module, cli_args.plots)
                if cli_args.hists != ['all']:
                    plots = [graph for graph in plots if graph.name in cli_args.hists]
                for plot in plots:
                    argList.append((filename, outpath, [plot], cli_args.signal, year, syst, cli_args))
    return argList



def run(infile, outpath, graphs, signalName, year, syst, args):
    # logging.info(f'Processing {graphs.name} for year {year} and systematic {syst}')
    kwargs = {
        'ratio_bot': args.ratio_range[0],
        'ratio_top': args.ratio_range[1],
        'extra_format': 'pdf',
    }

    scales = config.get_inputs(args.workdir, "scales")

    if args.type == "ntuple":
        ntuple = config.get_ntuple(args.region, obj=args.region_type)
        ginfo = ntuple.get_info(keep_dd_data=True, remove=['nonprompt_mc'])
        graphs = getattr(plots_module, args.plots)
        if args.hists != ['all']:
            graphs = [graph for graph in graphss if graph.name in args.hists]
    else:
        ginfo = config.get_ntuple_info(args.region, remove=['nonprompt_mc'])
        ntuple = None

    groups = ginfo.setup_groups()
    bkg = list(filter(lambda x: x not in [signalName, 'data'] ,groups.keys()))

    plotter = Plotter(infile, groups, ntuple=ntuple, year=year, ginfo=ginfo, scale=True)
    plotter.set_groups(sig=signalName, bkg=bkg, data='data')
    for scaler in scales.scale_list:
        scaler(plotter, year, syst)

    plotter.fill_hists(graphs)

    for graph in graphs:
        plotter.plot_from_graph(graph, outpath/'plots'/f'{graph.name}.png', **kwargs)

        # setup log file
        logger = LogFile(graph.name, constants.lumi[year], graph)
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
