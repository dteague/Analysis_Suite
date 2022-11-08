#!/usr/bin/env python3
import json
from pathlib import Path
from shutil import copy

import analysis_suite.commons.user as user

def writeHTML(path, name, subdir=None):
    """**Creates a useable HTML page of plots**

    Copies the necessary files (located in the html directory)

    Parameters
    ----------
    path : string
        Where the files should be saved
    name : string
        Name of the analysis used for titling
    channels : list of strings, optional
        List of channels to create the differnt sub-webpages
    """
    if subdir is None:
        subdir = []

    for filename in (user.analysis_area / "commons/html").glob("*"):
        copy(filename, path)

    info = { "Title" : name,
             "Subdirectory" : subdir}
    with open(f'{path}/extraInfo.json', 'w') as outjson:
        outjson.write(json.dumps(info))

def get_plot_area(name, path=None):
    www_path = user.www_area/name
    if path:
        www_path /= path.stem
    www_path /= time.strftime("%Y_%m_%d")
    return www_path

def make_plot_paths(path):
    (path/"plots").mkdir(exist_ok=True, parents=True)
    (path/"logs").mkdir(exist_ok=True, parents=True)
