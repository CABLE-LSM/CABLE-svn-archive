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
}



## unknown machine, user entering options stdout 
host_read()
{
   print "\nWhat is the root path of your NetCDF library and .mod file. Remember these have to be created by the same Fortran compiler you want to use to build CABLE. e.g./usr/local/intel"
   read NCDF_ROOT
   
   print "\nWhat is the path, relative to this root, of your NetCDF library. e.g. lib"
   read NCDF_DIR
   export NCDIR=$NCDF_ROOT/$NCDF_DIR
   
   print "\nWhat is the path, relative to this root, of your NetCDF .mod file. e.g. include"
   read NCDF_MOD
   export NCMOD=$NCDF_ROOT/$NCDF_MOD

   print "\nWhat is the Fortran compiler you wish to use. e.g. ifort, gfortran"
   read FC 
   export FC
   print "\nWhat are the approriate Fortran compiler to use. e.g.(ifort) -O2 -fp-model precise "
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
   print 'export NCDIR='$NCDF_ROOT'/'$NCDF_DIR >> junk
   print 'export NCMOD='$NCDF_ROOT'/'$NCDF_MOD >> junk
   print 'export FC='$FC >> junk
   print 'export CFLAGS='$CFLAGS >> junk
   print '}' >> junk
   print '' >> junk
   print '' >> junk
}


not_recognized()
{  
  print "this is not a recognized host for which we know the location of netcdf distribution and correct compiler switches"
  print "\nPlease enter these details as prompted. If this is a common location for CABLE users then Please email cable_help@nci.org.au so we can update the script. Otherwise you can modify the script yourself to avoid this message re-occuring. Just search build.ksh for instances of /'HOST_MACH example/' "

   print "To enter compile options for this build, press enter, otherwise Control-C to abort script"           
   
   host_read

   print "\nIf CABLE builds OK, should i enter this configuraion into the build script permanently? Y or N "
   read response_buildconfig 
   
   if [[ $response_buildconfig = 'Y' ]]; then
      write_buildconfig='true'  
      print "\nPlease supply a comment include the new build script. General the host URL e.g. vayu.nci.org.au "
      read HOST_COMM
   
   else
      write_buildconfig='false'  
   fi

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
         subr=host_${kh[0]}
         $subr
      fi        
      (( k = k + 1 ))
   done 
}


build_status()
{
   if [[ -f cable ]]; then
   	print '\nBUILD OK\n'
   	mv cable ../
   fi
}


      
i_do_now()
{
   if [[ -f cable ]]; then
   	print '\nBUILD OK\n'
   	mv cable ../
      
      if [[ $write_buildconfig = 'true'  ]]; then
         cd ../
         host_write
         tail -n +7 build.ksh > build.ksh.tmp
         cat junk build.ksh.tmp > build.ksh.new
         mv build.ksh.new build.ksh 
         rm -f build.ksh.tmp build.ksh.new 
      fi   
   fi
}


build_build()
{
   print 'building'
#   if [[ $1 = 'clean' ]]; then
#      rm -fr .tmp
#      exit
#   fi
#   
#   if [[ ! -d .tmp ]]; then
#      mkdir .tmp
#   fi
#   
#   if [[ -f cable ]]; then
#      print '\ncable executable exists. copying to cable.bu\n' 
#      mv cable cable.bu
#   fi
#   
#   CORE="../core/biogeophys"
#   DRV="."
#   CASA="../core/biogeochem"
#   
#   /bin/cp -p $CORE/*90 ./.tmp
#   /bin/cp -p $DRV/*90 ./.tmp
#   /bin/cp -p $CASA/*90 ./.tmp
#   
#   print '\nPlease note: CASA-CNP files are included in build only for technical reasons. Implementation is not officially available with the release of CABLE 2.0\n' 
#   /bin/cp -p Makefile_offline  ./.tmp
#   
   cd .tmp/
#   
#   make -f Makefile_offline
}

###########################################
## build.ksh - MAIN SCRIPT STARTS HERE   ##
###########################################

known_hosts

HOST_MACH=`uname -n | cut -c 1-4`

do_i_no_u

not_recognized
touch cable 
i_do_now

