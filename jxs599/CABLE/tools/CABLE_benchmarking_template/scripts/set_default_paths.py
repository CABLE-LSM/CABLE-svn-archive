#!/usr/bin/env python

"""
Set default paths on various machines...

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (06.06.2020)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import subprocess
import datetime

def set_paths(nodename):

    if "Mac" in nodename or "imac" in nodename:
        NCDIR = '/opt/local/lib/'
        NCMOD = '/opt/local/include/'
        FC = 'gfortran'
        FCMPI = 'gfortran-mp-9'
        CFLAGS = '-O2'
        LD = "'-lnetcdf -lnetcdff'"
        LDFLAGS = "'-L/opt/local/lib -O2'"

        #
        ## Met paths ...
        #
        #met_dir = "/Users/mdekauwe/research/CABLE_runs/met_data/plumber_met"
        met_dir = "/Users/mdekauwe/research/plumber_test"

    elif "unsw" in nodename:
        cmd = "module load netcdf-c/4.4.1.1-intel"
        cmd = "module load netcdf-f/4.4.4-intel"
        error = subprocess.call(cmd, shell=True)
        if error == 1:
            raise("Error loading netcdf libs")

        #NCDIR = '/share/apps/netcdf/intel/4.1.3/lib'
        #NCMOD = '/share/apps/netcdf/intel/4.1.3/include'
        NCDIR = '/share/apps/netcdf-f/intel/4.4.4/lib'
        NCMOD = '/share/apps/netcdf-f/intel/4.4.4/include'

        FC = 'ifort'
        FCMPI = 'mpif90'
        CFLAGS = '-O2'
        LD = "'-lnetcdf -lnetcdff'"
        LDFLAGS = "'-L/opt/local/lib -O2'"

        #
        ## Met paths ...
        #
        #met_dir = ("/srv/ccrc/data45/z3509830/CABLE_runs/Inputs/"
        #           "PLUMBER_sites/met")
        met_dir = ("/srv/ccrc/data04/z3509830/Fluxnet_data/"
                   "All_flux_sites_processed_PLUMBER2/"
                   "Post-processed_PLUMBER2_outputs/Nc_files/Met")
    else:
        sys.path.append("/opt/Modules/v4.3.0/init/")
        import python as mod
        #exec(open('/opt/Modules/v4.3.0/init/python.py').read())
        ver = "4.7.1"
        mod.module('unload', 'netcdf')
        mod.module('unload', 'openmpi')
        mod.module('load', 'netcdf/%s' % (ver))
        mod.module('load', 'intel-compiler/2019.3.199')
        mod.module('load', 'intel-mpi/2019.6.166')

        NCDIR = '/apps/netcdf/%s/lib' % (ver)
        NCMOD = '/apps/netcdf/%s/include' % (ver)
        FCMPI = 'mpif90'
        FC = 'ifort'
        CFLAGS = '-O2'
        LD = "'-lnetcdf -lnetcdff'"
        LDFLAGS = "'-L/opt/local/lib -O2'"
        #
        ## Met paths ...
        #
        #met_dir = ("/g/data1/w35/Shared_data/Observations/Fluxnet_data/"
        #           "FLUXNET2015/Processed_data/Missing_10%_Gapfill_20%/Daily")
        met_dir = ("/g/data/w35/Shared_data/Observations/Fluxnet_data/"
                   "Post-processed_PLUMBER2_outputs/Nc_files/Met")

    return (met_dir, NCDIR, NCMOD, FC, FCMPI, CFLAGS, LD, LDFLAGS)
