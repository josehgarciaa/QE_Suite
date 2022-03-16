

import subprocess


def run_pw(inpfile=None, logfile=None, shell=True, qe_path=""):
    try:
        with open(inpfile) as f:
            pass
    except:
        raise FileNotFoundError

    if logfile is None:
        logfile = inpfile+".log"

    proc= subprocess.run([qe_path+"pw.x -inp "+inpfile+" |tee "+logfile ], shell=True);


