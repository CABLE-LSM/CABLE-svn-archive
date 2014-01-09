#!/bin/csh 

#PBS -l walltime=2700
#PBS -l vmem=3072MB
#PBS -l ncpus=1
#PBS -l jobfs=2GB
#PBS -q express 
#PBS -p 10
#PBS -j oe
#PBS -N qs_CABLE 


#############################################################################
###  after a bit of book-keeping this csh script will compile CABLE       ###
###  from the "../src" directory and deposit the executable binary in     ###
###  this directory. CABLE is then run over each of the sites specified   ### 
###  in "main.nml" and output data moved into a created directory (out)   ###
###  labelled by the name of the site. finally, an "R" script to produce  ###
###  plots of flux data is called.                                        ### 
#############################################################################


if ( ($?PBS_JOBID) ) then
   cd $PBS_O_WORKDIR
endif    

source ~/.login

if ( ! ($?PBS_JOBID) ) then
   echo "This script is tested on NCI machines and requires R module."
   echo " " 
   echo "Please edit this script directly to modify the default behaviour described below."
   echo " " 
   echo "As well as a bit of book-keeping this csh script will:" 
   echo      " " 
   echo      "   1. compile CABLE from the ../src directory"       
   echo      "   2. run CABLE over each of the sites specified in main.nml"  
   echo      "   3. move output data into a new directory (out/*sitename*)"   
   echo      "   4. finally, flux data is plotted and left in out/*sitename*. "
   echo      " "
   echo      "     NB. plotted data compares the run version of CABLE (out_cable.nc),"  
   echo      "         the previous version of CABLE (old_cable.nc), and observations" 
   echo      " "
   echo      "Hit enter to proceed, enter Q to quit, H for more info "
   echo      " "
   
   set response = $<
   if ($response == 'Q') then
      echo "Adios"
      exit
   else if ($response == 'H') then
      echo "Adios - this is not really populated yet, but ..."
      echo " if not supplied by the user in this directory, old_cable.nc by default"
      echo " will be the same as observed data." 
      exit
   endif
endif

   ######################################################################
   ###  set up directory for this run - make output directory and     ###
   ######################################################################

if ( ! ($?PBS_JOBID) ) then
   if ( $1 != 'plot' ) then

      ### possibly helpful book-keeping to allow for immediate execution ### 
      ### and avoid accidently over-writing data - up to a point !!      ###
      if ( -d out.9 ) then
         echo "time to organize your out/ directory"
         exit
      endif
      foreach i ( 8 7 6 5 4 3 2 1 )
         @ j = $i + 1 
         if ( -d out.$i ) then
            mv out.$i out.$j
         endif
      end
      if ( -d out ) then
         mv out out.1
      endif
   
   endif
endif

if ( $1 != 'plot' ) then
   mkdir out
endif
   

   ######################################################################
   ###  manouvre, compile cable, manouvre some more                   ###
   ######################################################################

if ( ! ($?PBS_JOBID) ) then

   if ( $1 != 'plot' ) then
  
      if ( $1 == 'run') then

         echo 'cable executable already exists. Just run it.'        

      else   

         if ( -f cable ) then
            rm -f cable
         endif
      
         set mypwd = `pwd`
         
         pushd ../build
         
            build.ksh 
         
            if ( -f cable ) then
               echo 'CABLE executable built successfully and  will be copied to directory:'
               echo $mypwd 
            else
               echo 'ERROR: build unsuccessful'     
               exit
            endif
         
            /bin/cp cable $mypwd 
            
         popd
      endif
   
   endif
   
endif
      
######################################################################
### call CABLE.R in batch mode to avoid going into R first, and    ###
### then clean up this directory (these files have already been    ###
### dealt with in R-script )                                       ###
######################################################################

if ( $1 != 'plot' && $1 != 'make' ) then

   #remove any trace of previous runs
   rm -f fort.66 *nc 
   rm -f *.txt *00.bin *00.dat
   cd src/
   
   R CMD BATCH --slave CABLE.R
   
   cp CABLE.Rout ../out

   cd ../
   if ( -f out_cable.nc ) then
      echo 'CABLE appears to have run successfully.'
   else
      echo 'ERROR: run unsuccessful'     
   endif
   
   rm -f out_cable.nc restart_out.nc log_cable.txt fort.66

  ./cable_ncdf.pl ebal_tot_cncheck00    
  ./cable_ncdf.pl ebal_tot00
   ncl cable_map.ncl
  
  ./cable_ncdf.pl wbal_tot00
   ncl cable_map_w.ncl

endif

#######################################################################
#### call plot.R in batch mode to avoid going into R first          ###
#### pdfs will be left in out/*sitename* directory                  ###
#######################################################################

if ( $1 == 'plot' || $2 == 'plot') then
   cd src/
   R CMD BATCH --slave plot.R
   cp plot.Rout ../out
endif

if ( ($?PBS_JOBID) ) then
   cd src/
   R CMD BATCH --slave plot.R
   cp plot.Rout ../out
endif    

