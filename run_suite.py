#!/usr/bin/env python3
import multiprocessing as mp
import warnings
import logging
import sys
import argparse
import pkgutil
import shutil

import analysis_suite.commons.user as user

warnings.filterwarnings('ignore')

def get_cli():
    parser = argparse.ArgumentParser(prog="main", description="Central script for running tools in the Analysis suite")
    parser.add_argument('tool', type=str, help="Tool to run",
                    choices=['mva', 'plot', 'flatten', 'combine'])
    ##################
    # Common options #
    ##################
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area / x,
                        help="Working Directory")
    parser.add_argument("-j", type=int, default=1, help="Number of cores")
    parser.add_argument("--log", type=str, default="ERROR",
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help="Set debug status (currently not used)")
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016pre", "2016post", "2017", "2018"] if x == "all" \
                                   else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument("-s", "--systs", default="Nominal",
                        type=lambda x : [i.strip() for i in x.split(',')],
                        help="Systematics to be used")
    parser.add_argument("--unblind", action="store_true")

    if len(sys.argv) == 1:
        pass
    elif sys.argv[1] == 'flatten':
        ntupleInfo = [ f.name for f in pkgutil.iter_modules(path=[user.analysis_area/"ntuple_info"]) if not f.ispkg]
        parser.add_argument('-n', '--ntuple', required=True, choices= ntupleInfo,
                            help="Ntuple info class used for make root files")
        parser.add_argument('-i', '--inputs')
    elif sys.argv[1] == "mva":
        parser.add_argument('-t', '--train', action="store_true")
        parser.add_argument('-m', '--model', required=True, choices=['DNN', 'TMVA', 'XGBoost', "CutBased"],
                            help="Run the training")
        parser.add_argument("--save", action='store_true', help="Save training set")
        parser.add_argument("--plot", action='store_true')
        parser.add_argument('-r', '--region', default="signal",)
        parser.add_argument('--rerun', action='store_true')
    elif sys.argv[1] == "plot":
        parser.add_argument('-p', '--plots', default="plots")
        parser.add_argument('-n', '--name', default=None,
                            help='Name of directory used for storing the plots')
        parser.add_argument("--hists", default="all",
                            type=lambda x : [i.strip() for i in x.split(',')],
                            help="Pick specific histogram to plot")
        parser.add_argument("--drawStyle", type=str, default='stack',
                            help='Way to draw graph',
                            choices=['stack', 'compare', 'sigratio'])
        parser.add_argument("--ratio_range", nargs=2, default=[0.4, 1.6], type=float,
                            help="Ratio min ratio max (default 0.5 1.5)")
        parser.add_argument('-t', '--type', default='processed', choices=['processed', 'test', 'train', 'ntuple'],
                            help='Type of file in the analysis change')
        parser.add_argument('-r', '--region', default='Signal',
                            help='Region that this histogram will be plotted from')
    else:
        pass

    # Combos
    if len(sys.argv) > 1 and (sys.argv[1] == "plot" or sys.argv[1] == "combine"):
        parser.add_argument("-sig", "--signal", type=str, default='',
                            help="Name of the group to be made into the Signal")

    return parser.parse_args()

def setup(workdir):
    workdir.mkdir(exist_ok=True)
    inputs = workdir/'inputs.py'
    if not inputs.exists():
        shutil.copy(user.analysis_area/'data/python/inputs.py', inputs)

if __name__ == "__main__":
    cli_args = get_cli()
    setup(cli_args.workdir)
    logging.basicConfig(level=getattr(logging, cli_args.log, None))

    ##############
    # Setup jobs #
    ##############

    if cli_args.tool == "mva":
        from analysis_suite.machine_learning import job_main
    elif cli_args.tool == "plot":
        from analysis_suite.plotting import job_main
    elif cli_args.tool == "flatten":
        from analysis_suite.flatten import job_main
    elif cli_args.tool == "combine":
        from analysis_suite.combine import job_main
    else:
        exit()


    argList = job_main.setup(cli_args)
    func = job_main.run

    #############
    # Start Job #
    #############
    if cli_args.j == 1:
        [func(*al) for al in argList]
    else:
        with mp.Pool(cli_args.j) as pool:
            pool.starmap(func, argList)

    job_main.cleanup(cli_args)
