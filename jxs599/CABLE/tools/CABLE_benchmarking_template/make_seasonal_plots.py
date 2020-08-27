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

import numpy as np
from user_options import *

sys.path.append("scripts")
from benchmark_seasonal_plot import main as seas_plot

#
## Make seasonal plots ...
#
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

ofdir = os.path.join(run_dir, output_dir)
all_files = glob.glob(os.path.join(ofdir, "*.nc"))
sites = np.unique([os.path.basename(f).split(".")[0].split("_")[0] \
                   for f in all_files])
for site in sites:
    print(site)
    for sci_id, sci_config in enumerate(sci_configs):

        old_fname = glob.glob("%s/%s_*_R%d_S%d_out.nc" % \
                        (ofdir, site, 0, sci_id))[0]
        new_fname = glob.glob("%s/%s_*_R%d_S%d_out.nc" % \
                        (ofdir, site, 1, sci_id))[0]
        plot_fname = os.path.join(plot_dir, "%s_S%d.png" % (site, sci_id))
        seas_plot(old_fname, new_fname, plot_fname)
