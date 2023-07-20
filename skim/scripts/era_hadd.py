#!/usr/bin/env python3
import argparse
import subprocess
import analysis_suite.commons.user as user

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", required=True,
                    help="output filename")
parser.add_argument("-y", "--years", required=True,
                    type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                               else [i.strip() for i in x.split(',')],
                    help="Year to use")
parser.add_argument('--era', required=True)
parser.add_argument('-a', '--analysis', help="Analysis name")
parser.add_argument("-o", "--output", required=True,
                     help="Output file name (root not needed)")
args = parser.parse_args()

for year in args.years:
    output_dir = user.hdfs_area/f"{args.analysis}_{year}_{args.filename}"
    info_dir = user.submit_area/f"{args.analysis}_{year}_{args.filename}-analyze"

    hadd_files = ""
    for file_dir in info_dir.glob("analyze-*"):
        input_file = file_dir / f'{file_dir.name}.inputs'
        with open(input_file) as f:
            if args.era in f.read() and (output_dir/f"{file_dir.name}.root").exists():
                hadd_files += str(output_dir/f"{file_dir.name}.root") + " "


    subprocess.call(f"hadd -f -v 1 {args.output}_{year}.root {hadd_files}", shell=True)
