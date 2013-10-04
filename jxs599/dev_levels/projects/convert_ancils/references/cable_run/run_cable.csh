#!/bin/csh
#############################################################################
###  after a bit of book-keeping this csh script will compile CABLE       ###
###  from the "../src" directory and deposit the executable binary in     ###
###  this directory. CABLE is then run over each of the sites specified   ### 
###  in "main.nml" and output data moved into a created directory (out)   ###
###  labelled by the name of the site. finally, an "R" script to produce  ###
###  plots of flux data is called.                                        ### 
#############################################################################

######################################################################
### possibly helpful book-keeping to allow for immediate execution ### 
### and avoid accidently over-writing data - up to a point !!      ###
######################################################################
#if ( -d out.9 ) then
#   echo "time to organize your out/ directory"
#   exit
#endif
#foreach i ( 8 7 6 5 4 3 2 1 )
#   @ j = $i + 1 
#   if ( -d out.$i ) then
#      mv out.$i out.$j
#   endif
#end
#if ( -d out ) then
#   mv out out.1
#endif


######################################################################
###  set up directory for this run - make output directory and     ###
### delete old cable binary if it exists.                          ###
######################################################################
#mkdir out
#if ( -f cable ) then
#   rm -f cable
#endif

######################################################################
###  manouvre, compile cable, manouvre some more                   ###
######################################################################
#cd ../src/
#if ( -f cable ) then
#   rm -f cable
#endif
#make
#mv cable ../run/
#make clean
#cd ../run/src/
cd src/

######################################################################
### call CABLE.R in batch mode to avoid going into R first, and    ###
### then clean up this directory (these files have already been    ###
### dealt with in R-script )                                       ###
######################################################################

R CMD BATCH --slave CABLE.R
cd ../
rm -f out_cable.nc
rm -f restart_out.nc
rm -f log_cable.txt
rm -f fort.66
gprof cable > gprof.1.
cd src/

R CMD BATCH --slave CABLE.R
cd ../
rm -f out_cable.nc
rm -f restart_out.nc
rm -f log_cable.txt
rm -f fort.66
gprof cable > gprof.2.
cd src/

R CMD BATCH --slave CABLE.R
cd ../
rm -f out_cable.nc
rm -f restart_out.nc
rm -f log_cable.txt
rm -f fort.66
gprof cable > gprof.3.
cd src/







