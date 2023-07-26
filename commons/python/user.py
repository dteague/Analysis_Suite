#!/usr/bin/env python3
from pathlib import Path
import socket

lxplus_name = 'dteague'
user_name = Path().owner()

eos_area = Path(f'/store/user/{lxplus_name}')
hdfs_area = Path(f'/hdfs/store/user/{lxplus_name}')
hdfs_workspace = hdfs_area / 'workspace'

submit_area = Path(f'/nfs_scratch/{user_name}')
scratch_area = Path(f'/scratch/{user_name}')
combine_area = submit_area / "CMSSW_10_2_13/src"

analysis_area = Path(f'{__file__}').resolve().parents[2]
workspace_area = analysis_area / 'workspace'

xrd_tag = "root://cms-xrd-global.cern.ch/"

hostname = socket.gethostname()
if 'hep.wisc.edu' in hostname:
    www_area = Path.home()/'public_html'
    website = f'https://hep.wisc.edu/~{user_name}/'
elif 'uwlogin' in hostname or 'lxplus' in hostname:
    www_area = Path(f'/eos/home-{lxplus_name[0]}/{lxplus_name}/www')
    website = f'https://{lxplus_name}.web.cern.ch/{lxplus_name}/'
else:
    www_area = ""
    website = ""

def setup():
    workspace_area.mkdir(exist_ok=True)
    hdfs_workspace.mkdir(exist_ok=True)
    if not combine_area.exists():
        print("To setup combine, please Run the following commands")
        print(f"""
            export SCRAM_ARCH=slc7_amd64_gcc700
            cd {submit_area}
            cmsrel CMSSW_10_2_13
            cd CMSSW_10_2_13/src
            cmsenv
            git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
            git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
            pushd HiggsAnalysis/CombinedLimit
            git fetch origin
            git checkout v8.2.0
            popd
            scramv1 b clean; scramv1 b # always make a clean build
        """)

if __name__ == "__main__":
    setup()
