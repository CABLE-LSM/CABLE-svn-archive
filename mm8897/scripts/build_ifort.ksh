#!/bin/ksh


## test georgia with gcc for zlib hdf5 and netcdf-c and ifort for netcdf-fortran
clean_build()
{
      print '\ncleaning up\n'
      print '\n\tPress Enter too continue buiding, Control-C to abort now.\n'
      read dummy 
      rm -fr .tmpifort
}

build_status()
{
   if [[ -f .tmpifort/cable ]]; then
   	mv .tmpifort/cable ifort_cable
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
 
   if [[ ! -d .tmpifort ]]; then
      mkdir .tmpifort
   fi
   
   if [[ -f ifort_cable ]]; then
      print '\ncable executable exists. copying to dated backup file\n' 
      mv ifort_cable ifort_cable.`date +%d.%m.%y`
   fi
   
   # directories contain source code
   CORE="../core/biogeophys"
   DRV="."
   CASA="../core/biogeochem"
   
   /bin/cp -p $CORE/*90 ./.tmpifort
   /bin/cp -p $DRV/*90 ./.tmpifort
   /bin/cp -p $CASA/*90 ./.tmpifort
   
   print "\n\n\tPlease note: CASA-CNP files are included in build only for " 
   print "\ttechnical reasons. Implementation is not officially available with" 
   print "\tthe release of CABLE 2.0\n"
    
   /bin/cp -p Makefile_offline  ./.tmpifort
   
   cd .tmpifort/
   make -f Makefile_offline
 
}

###########################################
## build.ksh - MAIN SCRIPT STARTS HERE   ##
###########################################
source /home/mm/intel/bin/compilervars.sh intel64
if [[ $1 = 'clean' ]]; then
   clean_build
fi
export NCDIR='/home/mm/local/lib'
export NCMOD='/home/mm/local/include'
export LD_LIBRARY_PATH='${NCDIR}:{LD_LIBRARY_PATH}' 
export FC=ifort
export CFLAGS='-O0 -g -fp-model precise -check all'
#export CFLAGS='-fp-model precise'
export LD='-L${NCDIR} -lnetcdff -L${NCDIR} -lnetcdf'
export LDFLAGS='-L/home/mm/local/lib -O2'
build_build
cd ../
build_status
