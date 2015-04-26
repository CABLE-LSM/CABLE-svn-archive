#!/bin/bash

cd /srv/ccrc/data45/z3509830/CABLE_runs/Inputs/

DIRS=`find * -maxdepth 1 -type d -path "PLUMBER_sites_slope_*"`


for IN_DIR in $DIRS
do


    in_file="./Site_files/"${IN_DIR}".txt"



cat > $in_file << EOF
# Single sites over which CABLE is run are defined here. Lines NOT starting with 
# a "#" are interpreted by the run script (run.ksh) as meaningful and used as 
# varibles internally (this includes blank/space lines). 
#
# The first meaningful line is the short name of the site used only to name the 
# sub-directory containg output. The second and third meaningful lines are 
# passed to CABLE as arguments. These are the relative paths to the met. data 
# and the pool size.   
#
#
#Sites
Amplero
./../../Inputs/${IN_DIR}/AmpleroFluxnet.1.4_met.nc

Blodgett
./../../Inputs/${IN_DIR}/BlodgettFluxnet.1.4_met.nc

Bugac
./../../Inputs/${IN_DIR}/BugacFluxnet.1.4_met.nc

ElSaler
./../../Inputs/${IN_DIR}/ElSalerFluxnet.1.4_met.nc

ElSaler2
./../../Inputs/${IN_DIR}/ElSaler2Fluxnet.1.4_met.nc

Espirra
./../../Inputs/${IN_DIR}/EspirraFluxnet.1.4_met.nc

FortPeck
./../../Inputs/${IN_DIR}/FortPeckFluxnet.1.4_met.nc

Harvard
./../../Inputs/${IN_DIR}/HarvardFluxnet.1.4_met.nc

Hesse
./../../Inputs/${IN_DIR}/HesseFluxnet.1.4_met.nc

Howard
./../../Inputs/${IN_DIR}/HowardFluxnet.1.4_met.nc

Howlandm
./../../Inputs/${IN_DIR}/HowlandmFluxnet.1.4_met.nc

Hyytiala
./../../Inputs/${IN_DIR}/HyytialaFluxnet.1.4_met.nc

Kruger
./../../Inputs/${IN_DIR}/KrugerFluxnet.1.4_met.nc

Loobos
./../../Inputs/${IN_DIR}/LoobosFluxnet.1.4_met.nc

Merbleue
./../../Inputs/${IN_DIR}/MerbleueFluxnet.1.4_met.nc

Mopane
./../../Inputs/${IN_DIR}/MopaneFluxnet.1.4_met.nc

Palang
./../../Inputs/${IN_DIR}/PalangFluxnet.1.4_met.nc

Sylvania
./../../Inputs/${IN_DIR}/SylvaniaFluxnet.1.4_met.nc

Tumba
./../../Inputs/${IN_DIR}/TumbaFluxnet.1.4_met.nc

UniMich
./../../Inputs/${IN_DIR}/UniMichFluxnet.1.4_met.nc

#
#
#
# sitename - etc 
#Hyytiala
#
# met. file - etc
#data/sample_met/Hyytiala1997-2005_csiro.nc
#
# CNP pools size - etc 
#data/surface_data/poolcnpInHyytiala.csv
EOF

done
