#!/bin/csh -x

#=============Job Setup==================

setenv DIR ~ste69f/access/xaiyl         # Directory - where files are 
setenv RUNID xaiyl                      # Runid or Jobid from UMUI
setenv YR 5                             # Number of Years - run or otherwise
setenv REINIT 3                         # 1/3 - For processing when output files are 1/3 months in length
#setenv RESUB 3                          # Not used
#setenv Pmonth pa                        # Ext. letter for Files containing Monthly Means
#setenv Pdaily pe                        # Ext. letter for Files containing Daily Means & Min/Max

#=============Processing=================

setenv CNV2NC n                         # convert UM ppfiles to NetCDF before running umplot
#setenv PROC n                           # Not used

#=============Temp Varname===============

setenv tname temp_3                     # temp_`number` - variable name for tscrn

#=============Model Setup================

setenv MODEL m                          # c/m - cable or moses
setenv RES 48                           # 48/96 - n48 or n96 resolution
setenv MASK 1                           # 1/2 - Which land-sea mask ?
setenv CAL 360                          # 360/365/366
setenv CASA n                           # y/n - did you run with CASA-CNP ?

#=============UMPLOT Setup===============

setenv SPLIT n                          # y/n - ONLY if YR is a multiple of 5 and  >= 10
setenv TAY n                            # y/n - For comparision plots eg Taylor Plot, please set LDIR is TAY y
setenv LDIR ~kow014/access/xaank_run305                      
                                        # default: $HOME
                                        # ONLY needs to be set to something other than $HOME if $TAY is y E.g.
                                        # MOSES: ~ste69f/access/xaiyl_m48_20yrs .OR. ~ste69f/access/xagpb_n96_m_20yrs

#=============Format Setup===============

setenv JPEG n                           # y/n
#setenv CHFMT n                          # Not used yet
#setenv FMT pdf                          # pdf/eps/jpg - Not used yet

exit

