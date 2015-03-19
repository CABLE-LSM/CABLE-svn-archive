#!/usr/bin/python

#import python modules
import os
import sys

#import local, application specific modules
from TestSuite_common import Tprint
from TestSuite_dirs import root, root_app, src, bin, run, binSerial, binParallel 
from TestSuite_dirs import trunk, UM, offline, UMrun, UMsrc, SVNURL, cfgs
from TestSuite_dirs import TestUM, TestSerial, TestMPI 

def TestSuite_runner( cfg ):

###############################################################################

   # Open logfile and set echo to screen w Tprint()
   logfname = "runlog"
   mlogfile = open(root + "/" + logfname, "a") 
   Tprint(mlogfile, "\nStart.\n")
   
   os.chdir( root )
   Tprint(mlogfile, "Run Applications.\n")

   for i in range( len( cfg.path ) ):  

      #For serial runs      
      if(str(cfg.mode[i]) == '1'):

         if TestSerial is True: 
            Tprint(mlogfile, str( "Run Serial Model at " + cfg.name[i] ) )
            rundir = str( run + "/" + cfg.path[i] )
            os.system("/bin/cp " + binSerial + "/cable " + rundir ) 
#jhan: remove hard-wiringand put in config
            os.system("/bin/cp " + binSerial + "/.tr* " + rundir )
            os.system("/bin/cp " + cfgs + cfg.path[i] + "/run_cable " + rundir ) 
            os.system("/bin/cp " + cfgs + cfg.path[i] + "/cable.nml " + rundir ) 
            # GoTo rundir and execute
            os.chdir( rundir )
#jhan: remove hard-wiringand put in config
            os.system("ln -s .trunk_sumbal_TumbaFluxnet.1.3_met.nc .trunk_sumbal" ) 
            os.system("qsub run_cable" ) 
            Tprint(mlogfile, str( "Run Serial Model: qsubbed" ) )
         
      #For parallel runs      
      if(str(cfg.mode[i]) == '2'):
         
         if TestMPI is True: 
            Tprint(mlogfile, str( "Run Parallel Model at " + cfg.name[i] ) )
            # cp executable and namelist to rundir
            rundir = str( run + "/" + cfg.path[i] )
            os.system("/bin/cp " + binParallel + "/cable-mpi " + rundir ) 
            os.system("/bin/cp " + binParallel + "/.tr* " + rundir ) 
            os.system("/bin/cp " + cfgs + cfg.path[i] + "/run_cable-mpi " + rundir ) 
            os.system("/bin/cp " + cfgs + cfg.path[i] + "/cable.nml " + rundir ) 
        
            # GoTo rundir and execute
            os.chdir( rundir )
            os.system("ln -s /g/data1/wd9/MetForcing/Global/GSWP2/ gswp" ) 
            os.system("ln -s .trunk_sumbal_GSWP2 .trunk_sumbal" )
            os.system("qsub run_cable-mpi" ) 
            Tprint(mlogfile, str( "Run Parallel Model: qsubbed" ) )
 
# jhan: running UM build/ run manually from UMUI job vaheb
      #For UM runs      
      # UM executable is already in the right place         
      #if(str(cfg.mode[i]) == '3'):
      #   if TestUM is True: 
      #      
      #      Tprint(mlogfile, "Run UM Model" )
      #      os.chdir( "/home/599/jxs599/umui_runs/" + UMrun )
      #   
      #      # Run UM 
      #      os.system("qsub umuisubmit_run > " + root + "/um_runlog")
               



