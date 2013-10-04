#!/bin/ksh 

#############################################################################
###  as well as a bit of book-keeping this script will compile CABLE      ###
###  from source code  and deposit the executable binary in               ###
###  this directory. CABLE is then run over each of the sites specified   ### 
###  in "sites_main.nml" and output data moved into a created directory   ###
###  labelled by the name of the site. An "R" script to produce           ###
###  plots of flux data is called.                                        ### 
#############################################################################

#PBS -l walltime=2700
#PBS -l vmem=3072MB
#PBS -l ncpus=1
#PBS -l jobfs=2GB
#PBS -q express 
#PBS -p 10
#PBS -j oe
#PBS -N SS 

######################################################################
######################################################################
###  main script starts here. main script (run.ksh) called         ###
###  from command line and calls above funtions                    ###
######################################################################
######################################################################

#turning this on for debugging script can be useful
#set -x

#this doesnt work
#if we are on vayu (potentially cherx etc) or a qsub-ed job add module(s)
if [[ `uname -n | cut -c 1-4` == 'vayu'  ]] || [ -n "$PBS_JOBID" ] ; then
   module add R
fi

#if qsub-ed create flag file to work around user-process limit 
if [ -n "$PBS_JOBID" ]; then
   cd $PBS_O_WORKDIR
fi

#source file containing all the functions called from here
. ./functions.ksh

#so long as not qsub-ed AND script is submitted with no args 
if [ -z "$PBS_JOBID" ] && [[ $# == 0 ]]; then
   banner_welcome
fi

######################################################################
###  set up directory for this run - make output directory and     ###
### possibly helpful book-keeping to allow for immediate execution ### 
### and avoid accidently over-writing data - up to a point !!      ###
######################################################################

   if [[ $1 != 'plot' ]]; then
      book_keeping
      if [[ -d out ]]; then
         print "Adding to existing out directory" 
      else          
         mkdir out
     fi      
   fi

######################################################################
###  manouvre, compile cable, manouvre some more                   ###
###  NB. qsub-ed jobs assume an executable already exists          ###
######################################################################
if [[ ! -f cable ]]; then
   force_build=true
fi           

if [[ -z "$PBS_JOBID" ]] || [[ forcebuild==true ]]; then
   if [[ $1 == 'build' ]] || [[ $1 == 'all' ]]; then
      print "\n*** BUILDING CABLE ***\n"
      build_cable
   fi
fi

########################################################################
##### call CABLE.R in batch mode to avoid going into R first, and    ###
##### then clean up this directory (these files have already been    ###
##### dealt with in R-script )                                       ###
########################################################################
if [[ $1 == 'plot' ]] || [[ $2 == 'plot' ]] || [[ -n "$PBS_JOBID" ]] || [[ $1 == 'all' ]]; then
   run_cable
fi

#########################################################################
###### call plot.R in batch mode to avoid going into R first          ###
###### pdfs will be left in out/*sitename* directory                  ###
###### NB. "plot"s data from out/. therefore if plotting only out/    ###
###### must exist                                                     ###
#########################################################################

#if [[ $1 == 'plot' ]] || [[ $2 == 'plot' ]] || [[ -n "$PBS_JOBID" ]] || [[ $1 == 'all' ]]; then
#   plot_cable
#   #qsubed job hangs on this, removal of .wait at end of run allows qsub of next site
#   if [ -n "$PBS_JOBID" ]; then
#      rm -f ../../.wait
#   fi    
#fi

