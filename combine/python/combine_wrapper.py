#!/usr/bin/env python3
import subprocess

def runCombine(command, output=True, error=subprocess.STDOUT):
    cwd = getattr(runCombine, 'work_dir', ".")
    was_output = False
    was_error = False
    with subprocess.Popen([command], shell=True, stderr=error,
                          cwd=cwd, stdout=subprocess.PIPE) as process:
        process.wait()
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
