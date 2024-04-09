#!/usr/bin/env python3
import os
import argparse
from contextlib import contextmanager
import yaml

from analysis_suite.commons.info import fileInfo
import analysis_suite.commons.setup_functions as setup
import analysis_suite.commons.user as user

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")

@contextmanager
def rOpen(filename, option=""):
    rootfile = ROOT.TFile(filename, option)
    yield rootfile
    rootfile.Close()

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

def get_info_local(filename):
    year="2016"
    return {"year": year, "sampleName": "data", "selection": "test", "dataset": "DoubleMuon",
            "run": ""}

def get_info_general(filename):
    sampleName = filename.split('/')

    yearDict = {"UL16NanoAODAPVv": "2016preVFP",
                "UL16NanoAODv": "2016postVFP",
                "UL17": "2017",
                "UL18": "2018",
                "Run2016" : "2016postVFP",
                "Run2017" : "2017",
                "Run2018" : "2018",
                }
    data_regions = {
        "DoubleMuon" : "DoubleMuon",
        "SingleMuon" : "SingleMuon",

        "MuonEG" : "MuonEG",

        "DoubleEG": "DoubleEG",
        "EGamma": "DoubleEG",
        "SingleElectron": "SingleElectron",
    }

    year = None
    for yearName, yearkey in yearDict.items():
        if yearName in filename:
            if yearkey == "2016postVFP" and "HIPM" in filename:
                year = "2016preVFP"
            else:
                year = yearkey
            break

    dataset = "None"
    for dset, name in data_regions.items():
        if dset in sampleName:
            dataset = name
            break
    run = ""
    if "Run201" in filename:
        i = filename.index("Run201")+7
        run = filename[i:i+1]
    isUL = "UL"  in filename
    return {"year": year, "selection": "From_DAS", 'run': run,
            "sampleName": sampleName, "dataset": dataset}

def run_multi(start, evts, files, inputs, selector):
    fChain = ROOT.TChain()
    for fname in files:
        fChain.Add(f"{fname}/Events")

    selector = getattr(ROOT, inputs["MetaData"]["Analysis"])()
    rInputs = setInputs(inputs)

    with rOpen(f"tmp_{start}.root", "RECREATE") as rOutput:
        selector.SetInputList(rInputs)
        selector.setOutputFile(rOutput)
        fChain.Process(selector, "", evts, start)

        anaFolder = selector.getOutdir()
        for tree in [tree.tree for tree in selector.getTrees()]:
            anaFolder.WriteObject(tree, tree.GetTitle())
        for i in selector.GetOutputList():
            anaFolder.WriteObject(i, i.GetName())
    return

if __name__ == "__main__":
    analysis_dict = {}
    if (analysis_choices := setup.get_analyses()) is not None:
        analysis_choices.remove("BaseSelector")
        analysis_dict['choices'] = analysis_choices

    parser = argparse.ArgumentParser(prog="main")
    parser.add_argument("-i", "--infile", default ="No Input File")
    parser.add_argument("-o", "--outfile", default="output.root",)
    parser.add_argument("--local", action='store_true',
                        help="Add if file is local (because current code gets"
                        "file metadata from file name)")
    parser.add_argument("-v", "--verbose", default=-1,
                        help="Current levels are 1 for progress bar, 4 for shortened run (10 000 events)"
                        "and 9 for max output (only on 3 events)")
    parser.add_argument("-a", "--analysis", **analysis_dict)
    parser.add_argument("-j", "--cores", default = 1, type=int,
                        help="Number of cores to run over")
    parser.add_argument('-ns', '--no_syst', action='store_true',
                        help="Run with no systematics")
    parser.add_argument("-n", "--number_events", default=-1)
    args = parser.parse_args()

    inputfile = args.infile if (env := os.getenv("INPUT")) is None else env
    outputfile = args.outfile if (env := os.getenv("OUTPUT")) is None else env

    if "root" == inputfile[inputfile.rindex(".")+1:]:
        files = [inputfile]
    else:
        with open(inputfile) as f:
            files = [line.strip() for line in f]

    testfile = files[0]
    if "root://" in testfile and "user" not in testfile:
        details = get_info_general(testfile)
    else:
        details = get_info_local(testfile)

    groupName = fileInfo.get_group(details["sampleName"])
    if args.analysis:
        analysis = args.analysis
    else:
        with open(user.analysis_area/"data/.analyze_info") as f:
            analysis = f.readline().strip()
    get_shapes = analysis == "ThreeTop" and not args.no_syst

    # Setup inputs
    inputs = dict()
    inputs["MetaData"] = {
        "DAS_Name": '/'.join(details["sampleName"]),
        "Group": groupName,
        'Analysis': analysis,
        'Dataset': details['dataset'],
        'Selection': details["selection"],
        'Year': details["year"],
        'Xsec': 1,
        'isData': True,
        'DataRun': details['run']
    }
    if not args.local:
        inputs['MetaData'].update({'Xsec': fileInfo.get_xsec(groupName), 'isData': fileInfo.is_data(groupName)})
    inputs["Verbosity"] = args.verbose
    inputs["NEvents"] = args.number_events
    inputs["Systematics"] = get_shape_systs() if get_shapes else []

    # Possibly need to fix for fakefactor stuff
    rInputs = setInputs(inputs)
    if int(args.verbose) > 0:
        print(yaml.dump(inputs, indent=4, default_flow_style=False))

    # Run Selection
    fChain = ROOT.TChain()
    for fname in files:
        fChain.Add(f"{fname}/Events")

    import time
    now = time.time()
    if args.cores == 1:
        selector = getattr(ROOT, analysis)()
        with rOpen(outputfile, "RECREATE") as rOutput:
            selector.SetInputList(rInputs)
            selector.setOutputFile(rOutput)
            fChain.Process(selector, "")
            ## Output
            anaFolder = selector.getOutdir()
            anaFolder.WriteObject(getSumW(files), 'sumweight')
            for tree in [tree.tree for tree in selector.getTrees()]:
                anaFolder.WriteObject(tree, tree.GetTitle())
            for i in selector.GetOutputList():
                anaFolder.WriteObject(i, i.GetName())
