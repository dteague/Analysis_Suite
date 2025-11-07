#!/usr/bin/env python3
import argparse
import shutil
import subprocess
from analysis_suite.commons.setup_functions import get_analyses
import analysis_suite.commons.user as user

def setup_runfile():
    runfile_dir = user.analysis_area / "runfiles"
    if not runfile_dir.exists():
        runfile_dir.mkdir()
    return runfile_dir

def setup_args(runfile_dir):
    runfile_options = set()
    for f in runfile_dir.glob("*.dat"):
        name = f.stem
        runfile_options.add(name[:name.rfind("_", 0, name.rfind("_"))])

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", required=True,
                        choices = list(runfile_options),
                        help="output filename")
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                                else [i.strip() for i in x.split(',')],
                        help="Year to use")
    parser.add_argument('-t', '--types', default="mc,data",
                        type=lambda x : [i.strip() for i in x.split(',')],)
    parser.add_argument('-r', '--rerun', action='store_true')
    parser.add_argument('-a', '--analysis', required=True, choices=get_analyses())
    parser.add_argument('-d', '--dryrun', action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    runfile_dir = setup_runfile()
    args = setup_args(runfile_dir)

    with open(user.analysis_area/'data'/'.analyze_info', 'w') as f:
        f.write(args.analysis)

    farmout_requirements = '(TARGET.OpSysAndVer =?= "AlmaLinux9" || TARGET.OpSysAndVer =?= "CentOS9")'

    for typ in args.types:
        for year in args.years:
            analysis_dir = f"{args.analysis}_{year}_{args.filename}_{typ}"
            runfile = (runfile_dir / f'{args.filename}_{typ}_{year}.dat').resolve()
            if not runfile.exists():
                print(f"{runfile} not found!")
                continue

            if not args.rerun:
                shutil.rmtree(user.submit_area/f"{analysis_dir}-analyze", ignore_errors=True)
                shutil.rmtree(user.scratch_area/f"{analysis_dir}-analyze", ignore_errors=True)
                shutil.rmtree(user.hdfs_area/analysis_dir, ignore_errors=True)

            farmout_call = f"farmoutAnalysisJobs {analysis_dir}"
            farmout_call += f' --input-file-list={runfile}'
            farmout_call += f' --infer-cmssw-path'
            farmout_call += f' --fwklite analyze.py'
            farmout_call += f' --input-basenames-not-unique'
            farmout_call += rf" --base-requirements='{farmout_requirements}'"
            farmout_call += f' --use-singularity=rhel9'
            farmout_call += f' --output-dir={user.eos_area/analysis_dir}'
            farmout_call += " --resubmit-failed-jobs" if args.rerun else ""
            if args.dryrun:
                print(farmout_call)
            else:
                subprocess.call(farmout_call, shell=True)
