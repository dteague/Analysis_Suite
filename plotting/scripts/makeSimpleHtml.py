#!/usr/bin/env python3
import argparse
from pathlib import Path
from shutil import copy

import analysis_suite.commons.user as user #
from analysis_suite.plotting.html_utils import writeHTML


def writeOut(input_dir, name):
    output_dir = user.www_area / input_dir
    plot_dir = output_dir / "plots"

    output_dir.mkdir(exist_ok=True, parents=True)
    plot_dir.mkdir(exist_ok=True, parents=True)

    plots = user.analysis_area / input_dir
    for plot in plots.glob("*png"):
        copy(plot, plot_dir)

    writeHTML(output_dir, name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path_to_files', type=str, required=True)
    parser.add_argument('-n', '--name', type=str, required=True)
    args = parser.parse_args()

    sub_dirs = [str(f) for f in (user.analysis_area/args.path_to_files).glob("*") if f.is_dir()]
    if sub_dirs:
        for f in sub_dirs:
            writeOut(f[f.index("workspace"):], args.name)
    else:
        writeOut(args.path_to_files, args.name)

