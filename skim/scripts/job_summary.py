#!/usr/bin/env python3
import argparse
from prettytable import PrettyTable

import analysis_suite.commons.user as user

def dirs(directory):
    return [d.stem for d in directory.iterdir() if d.is_dir()]

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', action='store_true')
parser.add_argument("-y", "--years", help="Year to check")
parser.add_argument("-a", "--analysis", help="Analysis used")
parser.add_argument("-f", "--filename", help="runfile to check")
args = parser.parse_args()

hdfs_dirs = set(dirs(user.hdfs_area))
scratch_dirs = set([d.replace("-analyze", "") for d in dirs(user.submit_area)])
run_jobs = sorted(hdfs_dirs.intersection(scratch_dirs))

if args.data:
    run_jobs = filter(lambda x: "data" in x, run_jobs)
if args.analysis is not None:
    run_jobs = filter(lambda x: args.analysis in x, run_jobs)
if args.filename is not None:
    run_jobs = filter(lambda x: args.filename in x, run_jobs)
if args.years is not None:
    run_jobs = filter(lambda x: args.years in x, run_jobs)

table = PrettyTable(["Submitted", "Completed", "Diff", "Job Name"])
table.align = "r"
table.align["Job Name"] = "l"

for job in run_jobs:
    hdfs_num = len(list((user.hdfs_area / job).glob("*root")))
    nfs_num = len(list((user.submit_area / (job+"-analyze")).glob("analyze-*")))
    table.add_row([nfs_num, hdfs_num, nfs_num-hdfs_num, job])


print(table.get_string())
