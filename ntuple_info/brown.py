#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user
from copy import deepcopy
import awkward as ak



info = NtupleInfo(
    filename = '',
    trees = [],
    chan="",
    color_by_group = {
        'data': 'k',
        '4top': 'blue',
        'tth': 'darkgrey',
        'top': 'skyblue',
        'abcdnn': 'orange',
        'ttt_nlo': 'k',
    }
)
