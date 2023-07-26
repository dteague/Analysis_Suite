#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

measurement = NtupleInfo(
    filename = user.hdfs_workspace / 'befficiency/{year}/{workdir}/',
    trees = ['Signal'],
    color_by_group = {
        "ttt": "crimson",
        'nonprompt_mc': 'blue',
        "xg": "indigo",
        "ttw": "olivedrab",
        "tth": "goldenrod",
        "ttz": "steelblue",
        "ttXY": "teal",
        "rare": "deeppink",
        "tttt": "tomato",
    },
    cut = [
        lambda vg: vg["passZVeto_fake"] == 0,
    ],
)
