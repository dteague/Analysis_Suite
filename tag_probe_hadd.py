#!/usr/bin/env python3
import argparse
import subprocess
import shutil
import analysis_suite.commons.user as user

analysis = 'TagProbe'
mc_key = "DYJetsToLL_M-50"
mc_name = "DY_madgraph"
data = {
    "2016pre": [f"Run2016{sub}" for sub in 'BCDEF'],
    "2016post": [f"Run2016{sub}" for sub in 'FGH'],
    "2017": [f"Run2017{sub}" for sub in 'BCDEF'],
    '2018': [f"Run2018{sub}" for sub in 'ABCD'],
}
selection = 'dy_ll'
hdfs_out = '/store/user/dteague/TagAndProbe'

parser = argparse.ArgumentParser()
parser.add_argument("-y", "--years", required=True,
                    type=lambda x : ["2016", "2017", "2018"] if x == "all" \
                               else [i.strip() for i in x.split(',')],
                    help="Year to use")
args = parser.parse_args()

def get_list(name, groups):
    output_dir = user.hdfs_area/name
    info_dir = user.scratch_area/f"{name}-analyze"

    output = {group: [] for group in groups}
    for file_dir in info_dir.glob("analyze-*"):
        input_file = file_dir / f'{file_dir.name}.inputs'
        with open(input_file) as f:
            filename = f.read()
        for group in groups:
            rootfile = output_dir/f"{file_dir.name}.root"
            if group in filename and rootfile.exists():
                output[group].append(str(rootfile))
                break
    return output

def hdfs_mkdir(directory):
    subprocess.call(f"hdfs dfs -mkdir -p {hdfs_out}/{directory}", shell=True)

def hdfs_cp(files, out):
    subprocess.call(f'hadd -f -v 0 -j 10 tnp.root {" ".join(files)}', shell=True)
    shutil.move('tnp.root', user.hdfs_area/"TagAndProbe"/out/'tnp.root')


for year in args.years:
    name = f"{analysis}_{year}_{selection}"
    data_files = get_list(f'{name}_data', data[year])
    mc_files = get_list(f'{name}_mc', [mc_key])

    for subera in data[year]:
        hdfs_mkdir(f'{year}/{subera}')
        hdfs_cp(data_files[subera], f'{year}/{subera}')

    hdfs_mkdir(f'{year}/{mc_name}')
    hdfs_cp(mc_files[mc_key], f'{year}/{mc_name}')
