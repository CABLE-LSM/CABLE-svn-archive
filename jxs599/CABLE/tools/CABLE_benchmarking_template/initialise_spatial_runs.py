#!/usr/bin/env python

"""
Checkout CABLE repositories, build executables and generate a qsub script to
wrap cable benchmarking scripts when used on raijin.

We can't use the run_comparison script as raijin nodes have no internet
connection.

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (12.03.2019)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import datetime
import subprocess
from optparse import OptionParser

from user_options import *

sys.path.append("scripts")
from get_cable import GetCable
from build_cable import BuildCable
from generate_qsub_script import create_qsub_script

parser = OptionParser()
parser.add_option("-s", "--skipbuild", action="store_true", default=False,
                  help="Rebuild src?")
parser.add_option("-g", "--skipget", action="store_true", default=False,
                  help="Get src?")

(options, args) = parser.parse_args()


if options.skipget == False:
    #
    ## Get CABLE ...
    #
    G = GetCable(src_dir=src_dir, user=user)
    G.main(repo_name=repos[0], trunk=trunk) # Default is True
    G.main(repo_name=repos[1], trunk=False) # integration branch

elif options.skipbuild == False:

    #
    ## Build CABLE ...
    #
    B = BuildCable(src_dir=src_dir, NCDIR=NCDIR, NCMOD=NCMOD, FC=FC,
                   CFLAGS=CFLAGS, LD=LD, LDFLAGS=LDFLAGS)
    B.main(repo_name=repos[0])
    B.main(repo_name=repos[1])



#------------- Change stuff ------------- #
tmp_ancillary_dir = "global_files" # GSWP3 grid/mask file, temporarily

met_dir = "/g/data/wd9/MetForcing/Global/GSWP3_2017/"
start_yr = 1901
end_yr = 1901
walltime = "0:30:00"
qsub_fname = "qsub_wrapper_script_simulation.sh"
#------------- Change stuff ------------- #

cable_aux = os.path.join("../", aux_dir)
for repo_id, repo in enumerate(repos):
    cable_src = os.path.join(os.path.join("../", src_dir), repo)
    for sci_id, sci_config in enumerate(sci_configs):

        R = RunCable(met_dir=met_dir, log_dir=log_dir, output_dir=output_dir,
                     restart_dir=restart_dir, aux_dir=aux_dir, spin_up=spin_up,
                     cable_src=cable_src, qsub_fname=qsub_fname,
                     met_data=met_data, nml_fname=nml_fname, walltime=walltime,
                     tmp_ancillary_dir=tmp_ancillary_dir)
        R.initialise_stuff()
        R.setup_nml_file()
