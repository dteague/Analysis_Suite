#!/usr/bin/env python3
from analysis_suite.commons.info import NtupleInfo
import analysis_suite.commons.user as user

measurement = NtupleInfo(
    filename = user.hdfs_workspace / 'befficiency/{year}/{workdir}/',
    trees = ['Signal'],
)
