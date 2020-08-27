#!/usr/bin/env python

"""
Generate qsub script

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (09.03.2019)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import subprocess
import datetime

def create_qsub_script(ofname, ncpus, mem, wall_time, project, email_address):

    f = open(ofname, "w")

    f.write("#!/bin/bash\n")
    f.write("\n")
    f.write("#PBS -l wd\n")
    f.write("#PBS -l ncpus=%d\n" % (ncpus))
    f.write("#PBS -l mem=%s\n" % (mem))
    f.write("#PBS -l walltime=%s\n" % (wall_time))
    f.write("#PBS -q normal\n")
    f.write("#PBS -P %s\n" % (project))
    f.write("#PBS -j oe\n")
    f.write("#PBS -M %s\n" % (email_address))
    f.write("#PBS -l storage=gdata/w35+gdata/wd9\n")
    f.write("\n")
    f.write("\n")
    f.write("\n")
    f.write("\n")
    f.write("source activate sci\n")
    f.write("module add netcdf/4.7.1\n")
    f.write("python ./run_site_comparison.py --qsub\n")
    f.write("\n")

    f.close()

    os.chmod(ofname, 0o755)


if __name__ == "__main__":

    #------------- Change stuff ------------- #
    project = "w35"
    ofname = "benchmark_cable_qsub.sh"
    ncpus = 2
    mem = "32GB"
    wall_time = "00:30:00"
    email_address = "mdekauwe@gmail.com"
    # ------------------------------------------- #

    create_qsub_script(ofname, ncpus, mem, wall_time, project, email_address)
