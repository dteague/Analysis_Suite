#!/usr/bin/env python3
import uproot
import numpy as np
from pathlib import Path
import argparse

def scale_daniel_file(infile, outfile, scale, scale_groups=None):
    output = {}
    print(infile)
    with uproot.open(infile) as f:
        dirname = f.keys(recursive=False)[0]
        f = f[dirname]
        groups = [g for g in f.keys(recursive=False, cycle=False) if "_" not in g]
        groups.append('data_obs')
        systs = [g.split("_")[-1] for g in f.keys(recursive=False, cycle=False) if "_" in g]
        systs.remove('obs')
        systs = np.unique(systs)
        for group in groups:
            if scale_groups is None or group in scale_groups:
                output[group] = scale*f[group].to_boost()
            else:
                output[group] = f[group].to_boost()

        for syst in systs:
            for group in groups:
                hist = f.get(f'{group}_{syst}', None)
                if hist is None:
                    continue
                if scale_groups is None or group in scale_groups:
                    output[f'{group}_{syst}'] = scale*hist.to_boost()
                else:
                    output[f'{group}_{syst}'] = hist.to_boost()

    with uproot.recreate(outfile) as f:
        for key, val in output.items():
            f[key] = val

def scale_dylan_file(infile, outfile, scale, scale_groups=None):
    output = {}
    print(infile)
    with uproot.open(infile) as f:
        groups = [k for k, cls in f.classnames(recursive=False, cycle=False).items() if "TH1" in cls]
        systs = [k for k, cls in f.classnames(recursive=False, cycle=False).items() if "TDirectory" in cls]
        for group in groups:
            if scale_groups is None or group in scale_groups:
                output[group] = scale*f[group].to_boost()
            else:
                output[group] = f[group].to_boost()

        for syst in systs:
            output[syst] = {}
            for group, hist in f[syst].items(cycle=False):
                if scale_groups is None or group in scale_groups:
                    output[syst][group] = scale*hist.to_boost()
                else:
                    output[syst][group] = hist.to_boost()

    def write(f, value, key=""):
        if isinstance(value, dict):
            for k, val in value.items():
                write(f, val, f'{key}/{k}')
        else:
            f[key] = value

    with uproot.recreate(outfile) as f:
        write(f, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--datacard')
    parser.add_argument("-o", '--outdir', default=None)
    parser.add_argument("-s", '--scale', required=True)
    args = parser.parse_args()

    if args.outdir is None:
        outdir = Path(f"corr_{args.scale}")
    else:
        outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True, parents=True)

    files = []
    newcard = []
    with open(args.datacard, 'r') as f:
        for line in f.readlines():
            if "shapes" not in line:
                newcard.append(line)
                continue
            sline = line.split()
            infile = Path(sline[3])
            outfile = infile.name
            if "input" in outfile:
                outfile = outfile.replace("input", sline[2].split('_')[-1])
                line = line.replace(sline[4], sline[4].split('/')[-1])
                line = line.replace(sline[5], sline[5].split('/')[-1])
            line = line.replace(str(infile), outfile)
            newcard.append(line)

            files.append((infile, outfile))

    with open(outdir/args.datacard, 'w') as f:
        for line in newcard:
            f.write(line)

    for infile, outfile in files:
        if "combine" in outfile:
            scale_dylan_file(infile, outdir/outfile, float(args.scale))
        else:
            scale_daniel_file(infile, outdir/outfile, float(args.scale))
