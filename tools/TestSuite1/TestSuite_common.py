#!/usr/bin/python
import os
from TestSuite_dirs import root, root_app, UMrun, UMsrc, uHome, project,user 
from TestSuite_dirs import src, bin, binSerial, binParallel, run 

###############################################################################

def Tprint( filename, message ):
   print message
   filename.write(message + "\n") 

###############################################################################

def CleanSlate():  
   os.system("/bin/rm -f " + root + "/*log*" ) 
   os.system("/bin/rm -f " + root + "/*pyc" ) 
   os.system("/bin/rm -fr " + root_app ) 
   #os.system("/bin/rm -fr " + uHome +"/umui_runs/" + UMrun ) 
   #shortUser = "/short/" + project + "/" + user 
   #os.system("/bin/rm -fr " + shortUser + "/" + UMsrc ) 
   #os.system("/bin/rm -fr " + shortUser + "/UM_ROUTDIR/" + user + "/" + UMsrc ) 
   
################################################################################

def TestSuite_structure( cfg, logfile ):

   Tprint( logfile, "mkdir structure\n" )
   Tprint( logfile, root_app )
   Tprint( logfile, src )
   Tprint( logfile, bin )
   Tprint( logfile, binSerial )
   Tprint( logfile, binParallel )
   
   os.system("/bin/mkdir -p " + root_app ) 
   os.system("/bin/mkdir -p " + src ) 
   os.system("/bin/mkdir -p " + bin ) 
   os.system("/bin/mkdir -p " + binSerial ) 
   os.system("/bin/mkdir -p " + binParallel ) 
   os.system("/bin/mkdir -p " + run ) 
   
   for i in range( len( cfg.path ) ):   
      os.system("/bin/mkdir -p " + run + "/" + cfg.path[i] ) 
      Tprint( logfile, str( run + "/" +  cfg.path[i] ) )
               




