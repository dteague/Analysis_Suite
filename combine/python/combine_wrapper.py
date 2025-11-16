#!/usr/bin/env python3
import subprocess
import warnings
from analysis_suite.commons.histogram import Histogram

warnings.simplefilter("ignore", UserWarning)
from statsmodels.nonparametric.smoothers_lowess import lowess

def runCombine(command, output=True, error=subprocess.STDOUT, workdir=None):
    cwd = getattr(runCombine, 'work_dir', ".")
    if workdir is not None:
        cwd = workdir
    was_output = False
    was_error = False
    with subprocess.Popen([command], shell=True, stderr=error,
                          cwd=cwd, stdout=subprocess.PIPE) as process:
        # process.wait()
        was_error = process.returncode

        if output or was_error:
            was_output = bool(process.stdout.peek())
            for line in process.stdout:
                print(line.decode('utf8'), end="")
            if was_output:
                print("\n")
        else:
            for line in process.stdout:
                pass
    if was_error:
        raise Exception("Error in combine code!!")


def smooth_hist(nom, up, down, frac=0.67, it=5, symm=False):
    centers = nom.axis.centers
    if symm:
        up_ratio = 1+(up.vals-down.vals)/(2*nom.vals+1e-5)
        down_ratio = 1+(down.vals-up.vals)/(2*nom.vals+1e-5)
    else:
        up_ratio = (up.vals+1e-5)/(nom.vals+1e-5)
        down_ratio = (down.vals+1e-5)/(nom.vals+1e-5)

    if len(centers) > 2:
        up_ratio_lowess = lowess(up_ratio, centers, frac=frac, it=it).T[1]
        down_ratio_lowess = lowess(down_ratio, centers, frac=frac, it=it).T[1]
    else:
        up_ratio_lowess = up_ratio
        down_ratio_lowess = down_ratio

    up_lowess = Histogram(nom.axis)
    up_lowess.set_data(up_ratio_lowess*nom.vals, up.variances())
    down_lowess = Histogram(nom.axis)
    down_lowess.set_data(down_ratio_lowess*nom.vals, down.variances())

    if symm:
        up.set_data(up_ratio*nom.vals)
        down.set_data(down_ratio*nom.vals)
    return up_lowess, down_lowess
