#!/bin/ksh

known_hosts()
{
   set -A kh vayu 
}

## vayu.nci.org.au
host_vayu()
{
   NCDF_ROOT=/apps/netcdf/3.6.3
   export NCDIR=$NCDF_ROOT'/lib/Intel'
   export NCMOD=$NCDF_ROOT'/include/Intel'
   export FC=ifort
   export CFLAGS='-O2 -fp-model precise -ftz -fpe0'
   if [[ $1 = 'debug' ]]; then      
      export CFLAGS='-O0 -traceback -g -fp-model precise -ftz -fpe0' 
   fi
   build_build
   cd ../
   build_status
}



## unknown machine, user entering options stdout 
host_read()
{
   print "\n\tWhat is the root path of your NetCDF library" \
         "and .mod file. "
   print "\tRemember these have to be created by the same " \
         "Fortran compiler you" 
   print "\twant to use to build CABLE. e.g./usr/local/intel"
   read NCDF_ROOT
   
   print "\n\tWhat is the path, relative to this root, of " \
         "your NetCDF library." 
   print "\te.g. lib"
   read NCDF_DIR
   export NCDIR=$NCDF_ROOT/$NCDF_DIR
   
   print "\n\tWhat is the path, relative to this root, of " \
         "your NetCDF .mod file."
   print "\te.g. include"
   read NCDF_MOD
   export NCMOD=$NCDF_ROOT/$NCDF_MOD

   print "\n\tWhat is the Fortran compiler you wish to use."
   print "\te.g. ifort, gfortran"
   read FC 
   export FC

   print "\n\tWhat are the approriate compiler options"
   print "\te.g.(ifort) -O2 -fp-model precise "
   read CFLAGS 
   export CFLAGS
}


host_write()
{
   print '#!/bin/ksh' > junk
   print '' >> junk
   print 'known_hosts()' >> junk
   print '{' >> junk
   print '   set -A kh' ${kh[*]} $HOST_MACH >> junk
   print '}' >> junk
   print '' >> junk
   print '' >> junk
   print '## '$HOST_COMM >> junk
   print 'host_'$HOST_MACH'()' >> junk
   print '{' >> junk
   print '   export NCDIR='$NCDF_ROOT'/'$NCDF_DIR >> junk
   print '   export NCMOD='$NCDF_ROOT'/'$NCDF_MOD >> junk
   print '   export FC='$FC >> junk
   print '   export CFLAGS='$CFLAGS >> junk
   print '   build_build' >> junk
   print '   cd ../' >> junk
   print '   build_status' >> junk
   print '}' >> junk
   print '' >> junk
   print '' >> junk
}


not_recognized()
{  
   print "\n\n\tThis is not a recognized host for which we " \
         "know the location of the" 
   print "\tnetcdf distribution and correct compiler switches."

   print "\n\tPlease enter these details as prompted, and the " \
         "script will be " 
   print "\tupdated accordingly. " 
   print "\n\tIf this is a common machine for CABLE users, " \
         "please email"
   print "\n\t\t cable_help@nci.org.au "  
   print "\n\talong with your new build.ksh so that we can " \
         "update the script "
   print "\tfor all users. "
   print "\n\tTo enter compile options for this build press " \
         "enter, otherwise " 
   print "\tControl-C to abort script."           
   
   host_read

   print "\n\tPlease supply a comment include the new build " \
         "script." 
   print "\n\tGenerally the host URL e.g. vayu.nci.org.au "
   read HOST_COMM
   
   build_build
}


do_i_no_u()
{
   integer kmax=${#kh[*]}
   integer k=0
   typeset -f subr
   
   while [[ $k -lt $kmax ]]; do
      if [[ $HOST_MACH = ${kh[$k]} ]];then
         print 'Host recognized'
         subr=host_${kh[$k]}
         $subr
      fi        
      (( k = k + 1 ))
   done 
}


build_status()
{
   if [[ -f .tmp/cable ]]; then
   	mv .tmp/cable .
   	print '\nBUILD OK\n'
   else
      print '\nOooops. Something went wrong\n'        
   fi
}


      
i_do_now()
{
      cd ../
      host_write
      tail -n +7 build.ksh > build.ksh.tmp
      cat junk build.ksh.tmp > build.ksh.new
      mv build.ksh.new build.ksh
      chmod u+x build.ksh 
      rm -f build.ksh.tmp build.ksh.new junk 
      build_status
}


build_build()
{
   print 'building'
   if [[ $1 = 'clean' ]]; then
      rm -fr .tmp
      exit
   fi
   
   if [[ ! -d .tmp ]]; then
      mkdir .tmp
   fi
   
   if [[ -f cable ]]; then
      print '\ncable executable exists. copying to cable.bu\n' 
      mv cable cable.bu
   fi
   
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
## build.ksh - MAIN SCRIPT STARTS HERE   ##
###########################################

known_hosts

HOST_MACH=`uname -n | cut -c 1-4`

do_i_no_u

not_recognized

i_do_now

