#!/usr/bin/env python3
import os
import argparse
from analysis_suite.commons import FileInfo
import analysis_suite.commons.configs as configs
import tarfile
from pathlib import Path

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")

def setInputs(inputs):
    root_inputs = ROOT.TList()
    for key, data in inputs.items():
        subList = ROOT.TList()
        subList.SetName(key)
        if isinstance(data, dict):
            for dataName, dataVal in data.items():
                subList.Add(ROOT.TNamed(dataName, str(dataVal)))
        elif isinstance(data, list):
            for subitem in data:
                subList.Add(ROOT.TNamed(subitem, subitem))
        elif isinstance(data, (str, int, float, bool)):
            root_inputs.Add(ROOT.TNamed(key, str(data)))
            continue
        root_inputs.Add(subList)
    return root_inputs


def getSumW(infiles):
    runChain = ROOT.TChain()
    sumweight = ROOT.TH1F("sumweight", "sumweight", 1, 0, 1)
    for fname in infiles:
        runChain.Add(f"{fname}/metaInfo/metaInfo")

    if runChain.GetBranchStatus("summedWeights"):
        runChain.Draw("0>>sumweight",  "summedWeights")
    else:
        pass
    return sumweight

# Use for file in hdfs area
def get_info_local(filename):
    sampleName = filename.split('/')
    return {"year": "2017", "sampleName": sampleName}

def get_info_general(filename):
    sampleName = filename.split('/')

    yearDict = {"UL16NanoAODAPVv": "2016preVFP",
                "UL16NanoAODv": "2016postVFP",
                "UL17": "2017",
                "UL18": "2018",
                "Run2016" : "2016",
                "Run2017" : "2017",
                "Run2018" : "2018",
                }
    year = None
    for yearName in yearDict.keys():
        if yearName in filename:
            year = yearDict[yearName]
            break

    return {"year": year, "selection": "From_DAS",
            "sampleName": sampleName}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="main")
    parser.add_argument("-i", "--infile", default = "blah.in")
    parser.add_argument("-o", "--outfile", default="output.root")
    parser.add_argument("--test", action='store_true')
    parser.add_argument("-v", "--verbose", default=-1)
    parser.add_argument('-a', '--analysis', default='')
    args = parser.parse_args()

    inputfile = args.infile if (env := os.getenv("INPUT")) is None else env
    outputfile = args.outfile if (env := os.getenv("OUTPUT")) is None else env

    if "root" == inputfile[inputfile.rindex(".")+1:]:
        files = [inputfile]
    else:
        with open(inputfile) as f:
            files = [line.strip() for line in f]

    print(files)
    testfile = files[0]
    if "root://" in testfile and "user" not in testfile:
        details = get_info_general(testfile)
    else:
        details = get_info_local(testfile)
    info = FileInfo(**details)
    # groupName = info.get_group(details["sampleName"])
    groupName = "Test"
    datadir = Path(os.getenv("CMSSW_BASE"))/"src"/"analysis_suite"/"data"
    if args.analysis:
        analysis = args.analysis
    else:
        with open(datadir /".analyze_info") as f:
            analysis = f.readline().strip()
    print(groupName)

    # Setup inputs
    inputs = dict()
    inputs["MetaData"] = {
        "DAS_Name": '/'.join(details["sampleName"]),
        "Group": groupName,
        'Analysis': analysis,
        # 'Selection': details["selection"],
        'Year': details["year"],
    }
    if args.test:
        inputs['MetaData'].update({'Xsec': 1, 'isData': False})
    else:
        inputs['MetaData'].update({'Xsec': info.get_xsec(groupName), 'isData': info.is_data(groupName)})
    inputs["Verbosity"] = args.verbose
    inputs["Systematics"] = configs.get_shape_systs()
    # Possibly need to fix for fakefactor stuff
    rInputs = setInputs(inputs)

    # Run Selection
    fChain = ROOT.TChain()
    for fname in files:
        fChain.Add(f"{fname}/eemm/ntuple")

    selector = getattr(ROOT, analysis)()

    with configs.rOpen(outputfile, "RECREATE") as rOutput:
        selector.SetInputList(rInputs)
        selector.setOutputFile(rOutput)
        fChain.Process(selector, "")

        ## Output
        anaFolder = selector.getOutdir()
        anaFolder.WriteObject(getSumW(files), "sumweight")
        for tree in [tree.tree for tree in selector.getTrees()]:
            anaFolder.WriteObject(tree, tree.GetTitle())
        for i in selector.GetOutputList():
            anaFolder.WriteObject(i, i.GetName())
