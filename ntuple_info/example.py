#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

# All of the trees to be processed
trees = ["Test"]

# Must have the name info
info = NtupleInfo(
    # Input files. Can be a single file or directory
    filename = user.hdfs_area / 'workspace/signal_region/{year}',
    trees = trees,
    # Any cut applied on the tree using the VarGetter class
    cut=lambda vg : vg['variable'] > 100
)

# Mainly used for fake rates, so data needs to be renamed as a different bkg
info.add_change('Test', {'ttbar': 'TT'})
