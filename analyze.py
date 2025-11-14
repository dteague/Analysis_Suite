#!/usr/bin/env python3
import os
import argparse
from pathlib import Path
from contextlib import contextmanager
import yaml

from analysis_suite.commons.info import fileInfo
import analysis_suite.commons.setup_functions as setup
import analysis_suite.commons.user as user

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")

@contextmanager
def rOpen(filename, option=""):
    rootfile = ROOT.TFile(filename, option)
    yield rootfile
    try:
        rootfile.Close()
    except:
        print("problem closing file, will ignore, but double check!")


def get_shape_systs():
    from analysis_suite.data.systs import systematics
    systs = [syst.name for syst in systematics if syst.syst_type == "shape"]
    return list(dict.fromkeys(systs))


def setInputs(inputs):
    root_inputs = ROOT.TList()
    for key, data in inputs.items():
        if isinstance(data, (str, int, float, bool)):
            root_inputs.Add(ROOT.TNamed(key, str(data)))
            continue
        subList = ROOT.TList()
        subList.SetName(key)
        if isinstance(data, dict):
            for dataName, dataVal in data.items():
                subList.Add(ROOT.TNamed(dataName, str(dataVal)))
        elif isinstance(data, list):
            for subitem in data:
                subList.Add(ROOT.TNamed(subitem, subitem))
        root_inputs.Add(subList)
    return root_inputs


def getSumW(infiles):
    output = ROOT.TH1F('sumweight', 'sumweight', 14, 0, 14)
    LHESCALE, PDF, ALPHAZ = 1, 10, 12
    fChain = ROOT.TChain()
    for fname in files:
        fChain.Add(f"{fname}/Runs")
    for entry in fChain:
        sumW = entry.genEventSumw if hasattr(entry, 'genEventSumw') else -1
        output.Fill(0, sumW)
        if hasattr(entry, "LHEScaleSumw"):
            for i, scale in enumerate(entry.LHEScaleSumw):
                output.Fill(LHESCALE+i, scale*sumW)
        if hasattr(entry, "LHEPdfSumw"):
            if len(entry.LHEPdfSumw) >= 101:
                pdf = sorted([entry.LHEPdfSumw[i] for i in range(101)])
                err = (pdf[85] - pdf[15])/2
                output.Fill(PDF, (pdf[50]-err)*sumW)
                output.Fill(PDF+1, (pdf[50]+err)*sumW)
            # Alpha Z sumweight
            if len(entry.LHEPdfSumw) == 103:
                output.Fill(ALPHAZ, entry.LHEPdfSumw[101]*sumW)
                output.Fill(ALPHAZ+1, entry.LHEPdfSumw[102]*sumW)
            else:
                output.Fill(ALPHAZ, sumW)
                output.Fill(ALPHAZ+1, sumW)
    return output

def run_jetmet(infiles, inputs):
    meta = inputs['MetaData']
    keep_drop_file = f"kd_{'data' if meta['isData'] else 'mc'}.txt"
    year_dict = {
        '2016preVFP':  'UL2016_preVFP',
        '2016postVFP': 'UL2016',
        '2017':        'UL2017',
        '2018':        'UL2018',
    }
    jmeCorrections = createJMECorrector(
        (not meta['isData']), year_dict[meta['Year']], meta['DataRun'],
        "Merged", "AK4PFchs", noGroom=False
    )
    nentries = inputs['NEvents'] if inputs['NEvents'] > 0 else None
    p = PostProcessor(".", infiles, modules=[jmeCorrections()], provenance=True, maxEntries=nentries,
                      # cut="nMuon+nElectron>=2&&Sum$(Jet_pt)>200",
                      branchsel=user.analysis_area/'data'/keep_drop_file)
    p.run()

def run_ntuple(analysis, infiles, outfile, inputs):
    rInputs = setInputs(inputs)
    if int(args.verbose) > 0:
        print(yaml.dump(inputs, indent=4, default_flow_style=False))

    # Run Selection
    fChain = ROOT.TChain()
    for fname in infiles:
        fChain.Add(f"{fname}/Events")

    selector = getattr(ROOT, analysis)()
    with rOpen(outfile, "RECREATE") as rOutput:
        selector.SetInputList(rInputs)
        selector.setOutputFile(rOutput)
        fChain.Process(selector, "")
        # Output
        anaFolder = selector.getOutdir()
        anaFolder.WriteObject(getSumW(files), 'sumweight')
        for tree in [tree.tree for tree in selector.getTrees()]:
            anaFolder.WriteObject(tree, tree.GetTitle())
        for i in selector.GetOutputList():
            anaFolder.WriteObject(i, i.GetName())


if __name__ == "__main__":
    if (analysis_choices := setup.get_analyses()) is not None:
        analysis_choices.remove("BaseSelector")

    parser = argparse.ArgumentParser(prog="main")
    parser.add_argument("-i", "--infile", default ="No Input File")
    parser.add_argument("-o", "--outfile", default="output.root",)
    parser.add_argument("--local", action='store_true',
                        help="Add if file is local (because current code gets"
                        "file metadata from file name)")
    parser.add_argument("-v", "--verbose", default=-1,
                        help="Current levels are 1 for progress bar, 4 for shortened run (10 000 events)"
                        "and 9 for max output (only on 3 events)")
    parser.add_argument("-a", "--analysis", default=None, choices=analysis_choices)
    parser.add_argument("-j", "--cores", default = 1, type=int,
                        help="Number of cores to run over")
    parser.add_argument('-ns', '--no_syst', action='store_true',
                        help="Run with no systematics")
    parser.add_argument("-n", "--number_events", type=int, default=-1)
    args = parser.parse_args()

    inputfile = args.infile if (env := os.getenv("INPUT")) is None else env
    outputfile = args.outfile if (env := os.getenv("OUTPUT")) is None else env

    if "root" == inputfile[inputfile.rindex(".")+1:]:
        files = [inputfile]
    else:
        with open(inputfile) as f:
            files = [line.strip() for line in f]

    details = setup.get_details(files[0], outputfile, args.analysis)
    analysis = details['Analysis']
    get_shapes = (analysis in ["ThreeTop", 'BScale']) and not args.no_syst

    # Setup inputs
    inputs = dict()
    inputs["MetaData"] = details
    inputs["Verbosity"] = args.verbose
    inputs["NEvents"] = args.number_events
    inputs["Systematics"] = get_shape_systs() if get_shapes else []

    # if analysis == "ThreeTop":
    #     run_jetmet(files, inputs)
        # files = [Path(f.replace(".root", "_Skim.root")).name for f in files]
    run_ntuple(analysis, files, outputfile, inputs)
    # # Remove skimmed file if done
    # if analysis == "ThreeTop":
    #     for f in files:
    #         Path(f).unlink()
