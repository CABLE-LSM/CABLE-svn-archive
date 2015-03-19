#!/usr/bin/python 
__author__ = 'Jhan Srbinovsky'

# Build & Run various CABLE applications described in config file 

# usage: ./TestSuite_main -f <config file>

#import python modules
import os
import sys
import subprocess

#import local, application specific modules
from TestSuite_common import Tprint, CleanSlate, TestSuite_structure
from TestSuite_input import Locate_cfg, Configfile_interpreter
from TestSuite_build import TestSuite_builder
from TestSuite_run import TestSuite_runner 
from TestSuite_dirs import root, TestUM, TestSerial, TestMPI

###############################################################################

def main(argv):
   
   print "\n\nStart with a clean slate\n\n"
   CleanSlate()

   # Open logfile and set echo to screen w Tprint()
   logfname = "mainlog"
   mlogfile = open(root + "/" + logfname, "a") 
   Tprint(mlogfile, "\nStart.\n")
   
   # Locate config file (ifile), log file (ofile)
   # Defaults to using default.cfg as input, default.out as output (currently no output) 
   ifile = [] 
   ofile = []
   Locate_cfg( argv, ifile, ofile )
   Tprint(mlogfile, "Config file:")
   Tprint(mlogfile, str(ifile))

   # class to hold all config file input
   class Configs(object):
      def __init__(self):
         # Declare flags as mutable:
         self.mode = []
         self.name = []
         self.path = []
   cfg = []
   cfg = Configs() 
   
   # Read config file and fill class instance cfg
   Configfile_interpreter( ifile, cfg )
   
   Tprint( mlogfile, "\nApplications in config file to be evaluated:\n" ) 
   for i in range( len( cfg.name ) ):   
      Tprint(mlogfile,cfg.name[i] )
   
   ##print "\nSet up structure for building/running all Apps described in " + \
   ##"config file. These will be cleaned up following execution of the TestSuite"
   Tprint( mlogfile, "\nSet up structure for building/running all Apps described in " + \
   "config file.\n" ) 
   TestSuite_structure( cfg, mlogfile )
   
   #! Build Applications
   Tprint( mlogfile, "\nBuild Apps described in config file.\n" ) 
   TestSuite_builder( cfg )
   
   #! Run Applications
   Tprint( mlogfile, "\nRun Apps described in config file.\n" ) 
   TestSuite_runner( cfg )

   mlogfile.close()
   #! we may simply be able to use run.ksh here
   #! worry about qsub later. first run on a single node   

###############################################################################

   # Evaluate Success of Applications and write to ofile [default.out] 
   
   
   
   
#Atm_Step: Timestep                      2
   
   





# comment out if interactive
#################################################################################
if __name__ == "__main__":
   main(sys.argv[1:])

################################################################################


