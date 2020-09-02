#!/usr/bin/env python

"""
Get CABLE, build the executables, setup the run directory, run CABLE ... this is
a wrapper around the actual scripts (see scripts directory)

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (09.03.2019)"
__email__ = "mdekauwe@gmail.com"

import os
import shutil
import sys
import glob
import datetime
import subprocess
import numpy as np
from optparse import OptionParser


from user_options import *

sys.path.append("scripts")
from get_cable import GetCable
from build_cable import BuildCable
from run_cable_site import RunCable


parser = OptionParser()
parser.add_option("--qsub", action="store_true", default=False,
                  help="Run qsub script?")
parser.add_option("-s", "--skipsrc", action="store_true", default=False,
                  help="Rebuild src?")

(options, args) = parser.parse_args()

if options.qsub == False and options.skipsrc == False:

    #
    ## Get CABLE ...
    #
    G = GetCable(src_dir=src_dir, user=user)
    G.main(repo_name=repos[0], trunk=trunk) # Default is True

    # Run on a users branch, not integration
    if repos[1] != "integration":
        get_user_branch = True
    else:
        get_user_branch = False

    if share_branch:
        get_user_branch = False

    G.main(repo_name=repos[1], trunk=False, user_branch=get_user_branch,
           share_branch=share_branch) # integration branch

    #
    ## Build CABLE ...
    #
    B = BuildCable(src_dir=src_dir, NCDIR=NCDIR, NCMOD=NCMOD, FC=FC,
                   CFLAGS=CFLAGS, LD=LD, LDFLAGS=LDFLAGS)
    #B.main(repo_name=repos[0])

    if share_branch:
        print(os.path.basename(repos[1]))
        B.main(repo_name=os.path.basename(repos[1]))
    else:
        B.main(repo_name=repos[1])


#

#
## Run CABLE for each science config, for each repo
#

if not os.path.exists(run_dir):
    os.makedirs(run_dir)
os.chdir(run_dir)

cable_aux = os.path.join("../", aux_dir)
for repo_id, repo in enumerate(repos):
    cable_src = os.path.join(os.path.join("../", src_dir), repo)
    R = RunCable(met_dir=met_dir, log_dir=log_dir,
                 output_dir=output_dir, restart_dir=restart_dir,
                 aux_dir=cable_aux, namelist_dir=namelist_dir,
                 met_subset=met_subset, cable_src=cable_src, mpi=mpi,
                 num_cores=num_cores)
    for sci_id, sci_config in enumerate(sci_configs):
        R.main(sci_config, repo_id, sci_id)



os.chdir(cwd)
