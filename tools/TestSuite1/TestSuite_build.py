#!/usr/bin/python

#import python modules
import os
import sys

#import local, application specific modules
from TestSuite_common import Tprint
from TestSuite_dirs import root, root_app, src, bin, run, binSerial, binParallel 
from TestSuite_dirs import trunk, UM, offline, UMrun, UMsrc, SVNURL, cfgs
from TestSuite_dirs import UMcfgs, TestUM, TestSerial, TestMPI, uHome, project,user 

def TestSuite_builder( cfg ):

   # Open logfile and set echo to screen w Tprint()
   logfname = "buildlog"
   mlogfile = open(root + "/" + logfname, "a") 
   Tprint(mlogfile, "\nStart.\n")
   
   # checkout trunk (or URL to test - remove hardwiring)
   os.chdir(src)
   
   #use this on raijin
   Tprint( mlogfile, str( "checking out " + SVNURL ) ) 
   os.system("/usr/bin/svn co " + SVNURL ) 
   
   Tprint(mlogfile, "Building Models: \n")

   # overwrite build scripts checked out 
   #os.system("/bin/cp " + cfgs + "/build.ksh " + offline ) 
   #os.system("/bin/cp " + cfgs + "/build_mpi.ksh " + offline ) 
   
   # offline
   os.chdir(offline)
   
   # serial version
   if TestSerial is True: 
      Tprint(mlogfile, "Building Serial Model:")
      os.system("./build.ksh >> " + root + "/SerialBuildlog" )
      FileExists = os.path.isfile('cable') 
      if FileExists is True:
         Tprint(mlogfile, "Building Serial Model: DONE")
      else:
         sys.exit()
      Tprint(mlogfile, "Copying Serial model to the bin/ directory and cleaning up\n" )
      os.system("/bin/cp cable " + binSerial ) 
#jhan: remove hard-wiringand put in config
      os.system("/bin/cp .trunk_sumbal_TumbaFluxnet.1.3_met.nc " + binSerial ) 
      os.system("/bin/rm -fr .tmp" ) 

   # parallel version
   if TestMPI is True:
      Tprint(mlogfile, "Building Parallel Model:")
      os.system("./build_mpi.ksh >> " + root + "/ParallelBuildlog"  )
      FileExists = os.path.isfile('cable-mpi') 
      if FileExists is True:
         Tprint(mlogfile, "Building Parallel Model: DONE")
      else:
         sys.exit()
      Tprint(mlogfile, "Copying Parallel model to the bin/ directory and cleaning up\n" )
      os.system("/bin/cp cable-mpi " + binParallel ) 
      os.system("/bin/cp .trunk_sumbal_GSWP2 " + binParallel) 
      os.system("/bin/rm -fr .tmp" ) 

      # UM 
   if TestUM is True: 
      os.chdir(UM)
      Tprint(mlogfile, "Building libcable for UM Model:")

      os.system("./build.ksh >> " + root + "/libcableBuildlog" ) 
      os.system("rm -f " + uHome + "/short/" + UMsrc + "/new_sumbal " )
      os.system("rm -f " + uHome + "/short/" + UMsrc + "/bin/" + UMsrc + ".exe" )
      os.system("rm -f " + uHome + "/UM_ROUTDIR/" + UMsrc + "/ummodel/bin/" + UMsrc + ".exe" )
      Tprint( mlogfile, "Then move to accessdev for manual UM build....\n") 
   
# jhan: running UM build/ run manually from UMUI job vaheb
      ## cp UM runscripts to execute from
      #os.system("/bin/cp -r " + UMcfgs + "/" + UMrun + " " + uHome + "/umui_runs/ " ) 
      #
      ## for now we are using ppsrc in place at ~/UM_ROUTDIR
      ## cp UM Extracted src directory 
      ##os.system("/bin/cp -r -p " + UMcfgs+"/"+UMsrc+".tar.gz" + " /short/"+project+"/"+user+"/UM_ROUTDIR/"+ user+"/ " ) 
   #jh#an: run this manually for now
      ##os.system("/bin/rm " /short/"+project+"/"+user+"/UM_ROUTDIR/"+ user + "/" + UMsrc +"/ummodel/bin/") 
      #sys.exit() 
      ##os.chdir( "/short/"+project+"/"+user+"/UM_ROUTDIR/"+ user )
      ##os.system("/bin/tar xvfz " +UMsrc+".tar.gz" ) 
      #os.chdir( uHome + "/umui_runs/" + UMrun )
      #
      ## BUild UM 
      #os.system("./umuisubmit_compile > " + root + "/um_buildlog")
   
   mlogfile.close()
################################################################################




