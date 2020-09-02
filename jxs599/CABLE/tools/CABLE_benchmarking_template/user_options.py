import datetime
import os
import sys
import shutil

now = datetime.datetime.now()
date = now.strftime("%d_%m_%Y")
cwd = os.getcwd()
(sysname, nodename, release, version, machine) = os.uname()

sys.path.append("scripts")
from set_default_paths import set_paths


#------------- User set stuff ------------- #

#
## Qsub stuff ... ignore this block if not running a qsub script
#
project = "p66"
qsub_fname = "benchmark_cable_qsub.sh"
ncpus = 2
mem = "32GB"
wall_time = "01:30:00"
email_address = "srb001@csiro.au"

#
## Repositories to test, default is head of the trunk against personal repo.
## But if trunk is false, repo1 could be anything
#
user = "jxs599"
trunk = True
repo1 = "Trunk_%s" % (date)
#repo1 = "Trunk"
share_branch = True
repo2 = "test_jxs599"
repos = [repo1, repo2]


#
## user directories ...
#
src_dir = "src"
aux_dir = "src/CABLE-AUX"
run_dir = "runs"
log_dir = "logs"
plot_dir = "plots"
output_dir = "outputs"
restart_dir = "restart_files"
namelist_dir = "namelists"

if not os.path.exists(src_dir):
    os.makedirs(src_dir)

#
## Needs different paths for NCI, storm ... this is set for my mac
## comment out the below and set your own, see scripts/set_default_paths.py
#
(met_dir, NCDIR, NCMOD, FC, FCMPI, CFLAGS, LD, LDFLAGS) = set_paths(nodename)

#
## Met files ...
#
#met_subset = ['FI-Hyy_1996-2014_FLUXNET2015_Met.nc',\
#              'AU-Tum_2002-2017_OzFlux_Met.nc']
#met_subset = ['TumbaFluxnet.1.4_met.nc']

# Till fixed
#met_dir = "/g/data/w35/mgk576/research/CABLE_runs/met/Ozflux"
met_subset = ['AU-Tum_2002-2017_OzFlux_Met.nc','AU-How_2003-2017_OzFlux_Met.nc']
#met_subset = [] # if empty...run all the files in the met_dir

met_subset = [ 'IT-Amp_2003-2006_LaThuile_Met.nc' ] 
# ~Plumber sites from Anna's 2016 paper
met_subset = [ 'US-Blo_2000-2006_FLUXNET2015_Met.nc', \
 'HU-Bug_2003-2006_LaThuile_Met.nc', 'PT-Esp_2002-2004_LaThuile_Met.nc' ]
# 'US-FPe_2000-2006_LaThuile_Met.nc', 'US-Ha1_1992-2012_FLUXNET2015_Met.nc', \
# 'FR-Hes_1997-2006_LaThuile_Met.nc', 'AU-How_2003-2017_OzFlux_Met.nc',      \
# 'US-Ho1_1996-2004_LaThuile_Met.nc', 'FI-Hyy_1996-2014_FLUXNET2015_Met.nc', \
# 'BW-Ma1_2000-2000_LaThuile_Met.nc',    'ID-Pag_2002-2003_LaThuile_Met.nc', \
# 'US-Syv_2002-2008_FLUXNET2015_Met.nc', 'AU-Tum_2002-2017_OzFlux_Met.nc',   \
# 'ZA-Kru_2000-2002_FLUXNET2015_Met.nc', 'NL-Loo_1997-2013_FLUXNET2015_Met.nc', \
# 'US-UMB_2000-2014_FLUXNET2015_Met.nc' ]

#
## science configs
#
sci1 = {
        "cable_user%GS_SWITCH": "'medlyn'",
}

sci2 = {
        "cable_user%GS_SWITCH": "'leuning'",
}

sci3 = {
        "cable_user%FWSOIL_SWITCH": "'Haverd2013'",
}

sci4 = {
        "cable_user%FWSOIL_SWITCH": "'standard'",
}

sci5 = {
        "cable_user%GS_SWITCH": "'medlyn'",
        "cable_user%FWSOIL_SWITCH": "'Haverd2013'",
}

sci6 = {
        "cable_user%GS_SWITCH": "'leuning'",
        "cable_user%FWSOIL_SWITCH": "'Haverd2013'",
}


sci7 = {
        "cable_user%GS_SWITCH": "'medlyn'",
        "cable_user%FWSOIL_SWITCH": "'standard'",
}

sci8 = {
        "cable_user%GS_SWITCH": "'leuning'",
        "cable_user%FWSOIL_SWITCH": "'standard'",
}


#sci_configs = [sci1, sci2, sci3, sci4, sci5, sci6, sci7, sci8]
sci_configs = [sci2]

#
## MPI stuff
#
mpi = True
num_cores = ncpus # set to a number, if None it will use all cores...!

# ----------------------------------------------------------------------- #
