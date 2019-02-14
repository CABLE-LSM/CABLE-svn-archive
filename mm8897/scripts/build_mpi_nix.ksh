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
   if [[ -f .mpitmp/cable-mpi ]]; then
   	mv .mpitmp/cable-mpi .
   	print '\nBUILD OK\n'
   else
      print '\nOooops. Something went wrong\n'        
      print '\nKnow build issues:\n'        
      print '\nSome systems require additional library. \n'        
      print '\nEdit Makefile_offline; add -lnetcdff to LD = ...\n'        
   fi
   exit
}


      
build_build()
{

   # write file for consumption by Fortran code
   # get SVN revision number 
   CABLE_REV=`svn info | grep Revis |cut -c 11-18`
   if [[ $CABLE_REV="" ]]; then
      echo "this is not an svn checkout"
      CABLE_REV=0
      echo "setting CABLE revision number to " $CABLE_REV 
   fi         
   print $CABLE_REV > ~/.cable_rev
   # get SVN status 
   CABLE_STAT=`svn status`
   print $CABLE_STAT >> ~/.cable_rev
 
   if [[ ! -d .mpitmp ]]; then
      mkdir .mpitmp
   fi
   
   if [[ -f cable-mpi ]]; then
      print '\ncable-mpi executable exists. copying to a dated backup file\n' 
      mv cable-mpi cable-mpi.`date +%d.%m.%y`
   fi
   
   CORE="../core/biogeophys"
   DRV="."
   CASA="../core/biogeochem"
   
   /bin/cp -p $CORE/*90 ./.mpitmp
   /bin/cp -p $DRV/*90 ./.mpitmp
   /bin/cp -p $CASA/*90 ./.mpitmp
   
   print "\n\n\tPlease note: CASA-CNP files are included in build only for " 
   print "\ttechnical reasons. Implementation is not officially available with" 
   print "\tthe release of CABLE 2.0\n"
    
   /bin/cp -p Makefile_mpi  ./.mpitmp
   
  cd .mpitmp/

   make -f Makefile_mpi
}

###########################################
## build_nix.ksh - MAIN SCRIPT STARTS HERE   ##
###########################################

if [[ $1 = 'clean' ]]; then
   clean_build
fi
if [[ ! -d ~/CABLE-AUX ]];then
   set_up_CABLE_AUX
else
   print "\n~/CABLE-AUX is at least present.\n"
fi


   
#export FC=gfortran
export FC=mpif90
export CFLAGS='-O2 -x f95-cpp-input'
export LD='-lnetcdff'
export LDFLAGS="-L ${NCDIR} -O2"
build_build
cd ../
build_status
