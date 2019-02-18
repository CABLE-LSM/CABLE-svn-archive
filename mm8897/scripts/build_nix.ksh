# This is a much shortened version of cables build.ksh
# which is called by nix to compile cable.

# It can also be called manmually after entering the nix-shell 
# which is the normal procedure to work with changing cable sources.
# 
clean_build()
{
      print '\ncleaning up\n'
      print '\n\tPress Enter too continue buiding, Control-C to abort now.\n'
      read dummy 
      rm -fr .tmp
}



build_status()
{
   if [[ -f .tmp/cable ]]; then
   	mv .tmp/cable nix_cable
   	print '\nBUILD OK\n'
   else
      print '\nOooops. Something went wrong\n'        
      print '\nKnown build issues:\n'        
      print '\nSome systems require additional library. \n'        
      print '\nEdit Makefile_offline; add -lnetcdff to LD = ...\n'        
   fi
   exit
}
      
build_build()
{
   if [[ ${dosvn} -eq 1 ]] ; then
       # write file for consumption by Fortran code
       # get SVN revision number 
       CABLE_REV=`svn info | grep Revis |cut -c 11-18`
       if [[ $CABLE_REV = "" ]]; then
	   echo "this is not an svn checkout"
	   CABLE_REV=0
	   echo "setting CABLE revision number to " $CABLE_REV 
       fi         
       print $CABLE_REV > ~/.cable_rev
       # get SVN status 
       CABLE_STAT=`svn status`
       print $CABLE_STAT >> ~/.cable_rev
   fi
 
   if [[ ! -d .tmp ]]; then
      mkdir .tmp
   fi
   
   if [[ -f nix_cable ]]; then
      print '\ncable executable exists. copying to dated backup file\n' 
      mv nix_cable nix_cable.`date +%d.%m.%y`
   fi
   
   # directories contain source code
   CORE="../core/biogeophys"
   DRV="."
   CASA="../core/biogeochem"
   
   /bin/cp -p $CORE/*90 ./.tmp
   /bin/cp -p $DRV/*90 ./.tmp
   /bin/cp -p $CASA/*90 ./.tmp
   
   print "\n\n\tPlease note: CASA-CNP files are included in build only for " 
   print "\ttechnical reasons. Implementation is not officially available with" 
   print "\tthe release of CABLE 2.0\n"
    
   /bin/cp -p Makefile_offline  ./.tmp
   
   cd .tmp/
   make -f Makefile_offline
 
}

###########################################
## build_nix.ksh - MAIN SCRIPT STARTS HERE   ##
###########################################

if [[ $1 = 'clean' ]]; then
   clean_build
fi
export FC=gfortran
#export CFLAGS='-g -O0 -x f95-cpp-input -Wall'
export CFLAGS='-x f95-cpp-input'
export LD='-lnetcdff'
export LDFLAGS="-L ${NCDIR} -O2"
build_build
cd ../
build_status
