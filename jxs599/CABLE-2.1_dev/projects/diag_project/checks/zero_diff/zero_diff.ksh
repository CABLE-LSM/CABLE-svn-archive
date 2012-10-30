#!/bin/ksh

# NAME
#     zero_diff.ksh -- check that the CABLE model is bitwise identical between 
#                      two revisions where this is indeed inteneded
#
# SYNOPSIS
#     zero_diff.ksh -r [-n ncpus]
#     zero_diff.ksh -d [-n ncpus]

# DESCRIPTION
#     Executing the script zero_diff.ksh following an approriately          
#     configured run, will compare binary output produed by the run TO      
#     binary data from a standard run. Should the two versions be bitwise 
#     identical, the printed result of zero_diff.ksh will be ZERO. Else it 
#     will be greater than zero.
#
#     CABLE reads a namelist file (cable.nml) at runtime to configure 
#     the model. Setting the variable cable_user%RUN_DIAG_LEVEL to "zero" will 
#     result in a CALL to subroutine cable_diag() being executed. This  
#     subroutine is passed (amongst other things) the data to be written to file
#     and the node #(N-1) from which it is written (starting from 00). For a
#     single node (serial) run this is 00. The data is written to a file(s) 
#     called FLUXES00.bin [to FLUXES(N-1).bin]. An additional file is output,
#     FLUXES.dat, which contains information regarding FLUXES.bin necessary
#     to interpret/read this binary file in post-processing.
#
#     The data recorded to make the "zero difference" assessment is written to
#     FLUXES**.bin as follows. At each timestep, on each active tile, we sum the
#     sensible and latent heat fluxes. e.g. For a single-site with one active 
#     tile, and running for 1 year @ 30 minute intervals, we have 17520 numbers,
#     which are the sum of the fluxes. (For a typical ACCESS run for 1 day, 
#     ~576,000 numbers.) These numbers/files should be recorded for a known 
#     stable version. By default these should be renamed to std.bin & std.dat.
#     The files should be reproduced using a build from the code you wish to 
#     verify as "neutral".   
#
#     Following your RUN, the output file(s) [FLUXES**.bin] are left in your 
#     execution directory. Manouvere zero_diff.ksh, these flux files you wish to
#     compare, and the base files (std.bin & std.dat) into any same directory.  
#     Just sym-linking the files will do. Execute zero_diff.ksh. 
#
#     zero_diff.ksh is just a wrapper for the fortran program diff_main. This 
#     program has been built on vayu, and is executed using the full path name
#     in zero_diff.ksh. You might be lucky and the binary proves portable to 
#     another system. Otherwise contact cable_help@nf.nci.org.au to get the
#     source code.
#
#    -r  executes fortran binaries built with standard, single REAL precision, 
#        real*4. This should be used with FLUXES.bin and std.bin data produced
#        using the same precision
# 
#    -d  executes fortran binaries built with DOUBLE precision, real*8. 
#        This should be used with FLUXES.bin and std.bin data produced
#        using the same precision
# 
#    -n  NCPUS Where a parallel run has resulted in NCPUS output files from 
#        FLUXES00.bin to FLUXES(NCPUS-1).bin. The specified NCPUS directs the 
#        script to concatenate NCPUS files into a single file for processing. 
# 

# these are the teo data files to compare
set -A basename std FLUXES00	
CABLE_tools='/projects/access/CABLE-AUX/tools'
spec_tools='zero_diff'

if [[ $1 == '-r' ]]; then

   # if a run across multiple processors was conducted, this call 
   # concatenates all those files 
   if [[ $2 == '-n' ]]; then
      $CABLE_tools/$spec_tools/diag_cat_nodes_r4 ${basename[1]} $3 
   else
      set basename[1]="FLUXES"	
   fi  
   
   ### execute f90 program to do all the work desired
   $CABLE_tools/$spec_tools/diff_main_r4 ${basename[0]} ${basename[1]}

elif [[ $1 == '-d' ]]; then

   if [[ $2 == '-n' ]]; then
      $CABLE_tools/$spec_tools/diag_cat_nodes_r8 ${basename[1]} $3 
   else
      set basename[1]="FLUXES"	
   fi  
   
   $CABLE_tools/$spec_tools/diff_main_r8 ${basename[0]} ${basename[1]}

else
   print 'Supported KINDS are real*4 (-r) and real*8 (-d)'        
fi
































