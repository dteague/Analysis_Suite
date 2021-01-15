#!/usr/bin/env python3
import os
from analysis_suite.commons import FileInfo
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine( "gErrorIgnoreLevel = 1001;")


def setInputs(inputs, list_THX=[]):
    root_inputs = ROOT.TList()
    for th1 in list_THX:
        root_inputs.Add(th1)
    for key, value in inputs.items():
        root_inputs.Add(ROOT.TNamed(key, str(value)))
    return root_inputs


inputfile = os.getenv("INPUT")
outputfile = os.getenv("OUTPUT")

if os.getenv("INPUT") is None:
    inputfile = "blah.in"
if os.getenv("OUTPUT") is None:
    outputfile = "output.root"


files = list()
with open(inputfile) as f:
    for line in f:
        files.append(line.strip())

print(files)
sampleName = files[0].split('/')
analysisName = sampleName[sampleName.index("user")+2]

analysis, year, selection = analysisName.split("_")
info = FileInfo(analysis=analysis, year=year, selection=selection)

inputs = dict()
groupName = info.get_group(sampleName)
inputs["Group"] = groupName
inputs["xsec"] = info.get_xsec(groupName)
inputs['Analysis'] = analysis
inputs['Selection'] = selection
inputs["Year"] = year

fChain = ROOT.TChain()
for fname in files:
    fChain.Add(f"{fname}/Events")

selector = getattr(ROOT, "ThreeTop")()
rInputs = setInputs(inputs)
selector.SetInputList(rInputs)
fChain.Process(selector, "")

print("here")

## Output
rOutput = ROOT.TFile(outputfile, "RECREATE")
rOutput.cd()
anaFolder = ROOT.gDirectory.mkdir(groupName)
anaFolder.cd()
for i in selector.GetOutputList():
    i.Write()
    
runChain = ROOT.TChain()
sumweight = ROOT.TH1F("sumweight", "sumweight", 1, 0, 1)
for fname in files:
    runChain.Add(f"{fname}/Runs")
runChain.Draw("0>>sumweight",  "genEventSumw")
anaFolder.cd()
sumweight.Write()

rOutput.Close()
