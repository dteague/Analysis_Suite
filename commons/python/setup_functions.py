#!/usr/bin/env python3
from xml.dom.minidom import parse
from analysis_suite.commons.info import fileInfo

from .user import analysis_area

def get_analyses():
    xml_filename = (analysis_area/'skim/src/classes_def.xml').resolve()
    if xml_filename.exists():
        xml_classes = parse(str(xml_filename))
        return [c.getAttribute("name") for c in xml_classes.getElementsByTagName('class')]
    else: # For on cluster where xml doesn't exists
        return None

def get_info_local(filename):
    return {
        "Year": '2016',
        "DAS_Name": "data",
        "Selection": "test",
        "Dataset": "DoubleMuon",
        "DataRun": "",
    }

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

    groupName = fileInfo.get_group(sampleName)

    return {
        "Year": year,
        "Selection": "From_DAS",
        'DataRun': run,
        "DAS_Name": filename,
        "Dataset": dataset,
        'Xsec': fileInfo.get_xsec(groupName),
        'isData': fileInfo.is_data(groupName),
        "Group": groupName,
    }

def get_details(inputfile, outputfile, analysis=None):
    if "root://" in inputfile and "user" not in inputfile:
        details = get_info_general(inputfile)
    else:
        details = get_info_local(inputfile)

    if analysis is None:
        with open(analysis_area/"data/.analyze_info") as f:
            analysis = f.readline().strip()
    details['Analysis'] = analysis

    return details
