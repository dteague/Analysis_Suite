#!/usr/bin/env python3
import argparse
import subprocess
import analysis_suite.commons.user as user
from pathlib import Path
import shutil
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", required=True,
                    help="output filename")
parser.add_argument("-y", "--years", required=True,
                    type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                               else [i.strip() for i in x.split(',')],
                    help="Year to use")
parser.add_argument('-a', '--analysis', required=True,
                    help='Analysis used in running jobs')
parser.add_argument('-t', '--types', default="mc,data",
                    type=lambda x : [i.strip() for i in x.split(',')],
                    help="Types of files to hadd, namely MC, Data, or both (will be in separate files)")
parser.add_argument('-d', '--workdir', default=datetime.now().strftime("%y%m%d"),
                    help="Output file to store/label the root files (default is date)")
parser.add_argument('-r', '--rerun', action="store_true",
                    help="Add rerun files to list of hadd'ed files")
args = parser.parse_args()

base_dir = user.hdfs_workspace / args.filename
base_dir.mkdir(exist_ok=True)

for typ in args.types:
    for year in args.years:
        output_dir = base_dir / year / args.workdir
        output_dir.mkdir(parents=True, exist_ok=True)
        input_dir = user.hdfs_area / f"{args.analysis}_{year}_{args.filename}_{typ}"
        if not input_dir.exists():
            continue
        files = f"{input_dir}/*root "
        if args.rerun and Path(f"{input_dir}_rerun").exists():
            files += f"{input_dir}_rerun/*root"
        subprocess.call(f"hadd -f -v 1 {args.filename}_{typ}_{year}.root {files}", shell=True)
        shutil.move(f'{args.filename}_{typ}_{year}.root', output_dir)
