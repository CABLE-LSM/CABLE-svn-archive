#!/bin/ksh

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
}


## shine-cl.nexus.csiro.au 
host_shin()
{
   NETCDF_ROOT=/usr/local/intel
   export NCDIR=$NCDF_ROOT/lib
   export NCMOD=$NCDF_ROOT/include
   export FC='ifort'
   export CFLAGS='-O2 -fp-model precise -ftz -fpe0'
   if [[ $1 = 'debug' ]]; then      
      export CFLAGS='-O0 -traceback -g -fp-model precise -ftz -fpe0' 
   fi
}


## NOTE: architecture on burnet and cherax doesn't like -ftz -fpe0 
## cherax.hpsc.csiro.au 
host_cher()
{
   NETCDF_ROOT=/usr/local/intel
   export NCDIR=$NCDF_ROOT/lib
   export NCMOD=$NCDF_ROOT/include
   export FC='ifort'
   export CFLAGS='-O2 -fp-model precise' 
}


## NOTE: architecture on burnet and cherax doesn't like -ftz -fpe0 
## burnet.hpsc.csiro.au 
host_burn()
{
   NETCDF_ROOT=/usr/local/intel
   export NCDIR=$NCDF_ROOT/lib
   export NCMOD=$NCDF_ROOT/include
   export FC='ifort'
   export CFLAGS='-O2 -fp-model precise' 
}
## unknown machine, user entering options stdout 
host_read()
{
   print "\nWhat is the root path of your NetCDF library and .mod file. Remember these have to be created by the same Fortran compiler you want to use to build CABLE. e.g./usr/local/intel"
   read NETCDF_ROOT
   print "\nWhat is the path, relative to this root, of your NetCDF library. e.g. lib"
   read NCDIR
   export NCDIR
   print "\nWhat is the path, relative to this root, of your NetCDF .mod file. e.g. include"
   read NCMOD
   export NCMOD
   print "\nWhat is the Fortran compiler you wish to use. e.g. ifort, gfortran"
   read FC 
   export FC
   print "\nWhat are the approriate Fortran compiler to use. e.g.(ifort) -O2 -fp-model precise "
   read CFLAGS 
   export CFLAGS
}



##HOST_MACH example
host_xxxx()
{
   NCDF_ROOT=
   export NCDIR=$NCDF_ROOT/
   export NCMOD=$NCDF_ROOT/
   export FC= 
   export CFLAGS= 
   if [[ $1 = 'debug' ]]; then      
      export CFLAGS='' 
   fi
}


recognized_host()
{
   if [[ $HOST_MACH = 'vayu' || \
         $HOST_MACH = 'cher' || \
         $HOST_MACH = 'shin' || \
         ##HOST_MACH example
         #$HOST_MACH = 'xxxx' || \
         $HOST_MACH = 'burn'    \
            ]]; then
   
      if [[ $HOST_MACH = 'vayu' ]]; then
         host_vayu $1
      fi       
      if [[ $HOST_MACH = 'shin' ]]; then
         print 'if cond 1'
         host_shin $1
      fi       
      if [[ $HOST_MACH = 'cher' ]]; then
         host_cher $1
      fi       
      if [[ $HOST_MACH = 'burn' ]]; then
         host_burn $1
      fi       
      ##HOST_MACH example
      #if [[ $HOST_MACH = 'xxxx' ]]; then
      #   host_xxxx $1
      #fi       
   
   else
     
     print "this is not a recognized host for which we know the location of netcdf distribution and correct compiler switches"
     print "\nPlease enter these details as prompted. If this is a common location for CABLE users then Please email cable_help@nci.org.au so we can update the script. Otherwise you can modify the script yourself to avoid this message re-occuring. Just search build.ksh for instances of /'HOST_MACH example/' "
   
   print "To enter compile options for this build, press enter, otherwise Control-C to abort script"           
      if [[ $HOST_MACH = 'read' ]]; then
         host_read
      fi       
   fi
}

## build.ksh - MAIN SCRIPT STARTS HERE

HOST_MACH=`uname -n | cut -c 1-4`

##HOST_MACH example. run the above command in  shell, or hard-wire $HOST_MACH
#HOST_MACH='xxxx'

recognized_host $1

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

print '\nPlease note: CASA-CNP files are included in build only for technical reasons. Implementation is not officially available with the release of CABLE 2.0\n' 
/bin/cp -p Makefile_offline  ./.tmp

cd .tmp/

make -f Makefile_offline
   
if [[ -f cable ]]; then
	print '\nBUILD OK\n'
	mv cable ../
fi

