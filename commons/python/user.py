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

if __name__ == "__main__":
    setup()
