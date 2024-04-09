#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

# All of the trees to be processed
trees = ["Test"]

# Must have the name info
info = NtupleInfo(
    # Input files. Can be a single file or directory
    filename = user.hdfs_area / 'workspace/signal_region/{year}/{workdir}/',
    trees = trees,
    # Any cut applied on the tree using the VarGetter class
    cut=[
        lambda vg : vg['variable'] > 100,
    ],
    # Groups used with the color used for plotting
    color_by_group = {
        'group1': 'red',
        'group2': 'orange',
        'group3': 'yellow',
        'group4': 'green',
        'group5': 'blue',
        'group6': 'indigo',
        'group7': 'violet',
    }
)

# Mainly used for fake rates, so data needs to be renamed as a different bkg
info.set_groups_trees('Test', ['group1', 'group2'])
