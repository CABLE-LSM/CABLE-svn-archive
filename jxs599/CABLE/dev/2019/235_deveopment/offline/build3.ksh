#!/bin/ksh

#build script for CABLE
#Simplified from previous build script. Only set up for raijin 
#Usage: 
#> ./build3
#   OR
#> ./build3 mpi
#builds MPI version of CABLE

## gadi.nci.org.au
host_gadi()
{
   . /etc/bashrc
   module purge
   module add intel-compiler/2019.5.281
   module add netcdf/4.6.3

   export FC='ifort'
   export NCDIR=$NETCDF_ROOT'/lib/Intel'
   export NCMOD=$NETCDF_ROOT'/include/Intel'
   export CFLAGS='-O2 -fp-model precise'
   if [[ $1 = 'debug' ]]; then
      export CFLAGS='-O0 -traceback -g -fp-model precise -ftz -fpe0'
   fi
   export LDFLAGS='-L'$NCDIR' -O0'
   export LD='-lnetcdf -lnetcdff'
   build_build
   cd ../
   build_status
}





## raijin.nci.org.au
host_raij()
{
   module unload intel-fc netcdf
   if [[ $1 = 'mpi' ]]; then
   		module del openmpi
      module add intel-mpi/5.1.3.210 intel-fc/17.0.1.132 openmpi/1.10.2
      export FC='mpif90'
	 else
   		module add intel-fc/16.0.1.150
      export FC=$F90
      export FC='ifort'
   fi
	 module add netcdf/4.3.3.1

   #export NCDIR=$NETCDF_ROOT'/lib/Intel'
   export NCDIR=$NETCDF_ROOT'/lib/'
#   echo "$NETCDF_ROOT"
#   exit
   export NCMOD=$NETCDF_ROOT'/include/'
   #export NCMOD=$NETCDF_ROOT'/include/Intel'
   #export CFLAGS='-O0 -fp-model precise'
   #export CFLAGS='-O0 -g'
   #if [[ $1 = 'debug' ]]; then
       export CFLAGS='-O0 -traceback -g -fp-model precise -ftz -fpe0 -check all,noarg_temp_created'
   #fi
   export LDFLAGS='-L'$NCDIR' '
   export LD='-lnetcdf -lnetcdff'
   
   if [[ $1 = 'mpi' ]]; then
    build_build mpi
   else
    build_build
   fi 
   cd ../
   build_status
}

clean_build()
{
      rm -fr .tmp
      print '\ncleaning up\n'
      print '\n\tPress Enter too continue buiding, Control-C to abort now.\n'
      read dummy
}


build_build()
{
   if [[ ! -d .tmp ]]; then
      mkdir .tmp
   fi

   if [[ -f cable ]]; then
      print '\ncable executable exists. copying to dated backup file\n'
      mv cable cable.`date +%d.%m.%y`
   fi

   # directories contain source code
   ALB="../science/albedo"
   RAD="../science/radiation"
   CAN="../science/canopy"
   CNP="../science/casa-cnp"
   GWH="../science/gw_hydro"
   MIS="../science/misc"
   ROU="../science/roughness"
   SOI="../science/soilsnow"
   OFF="../offline"
   UTI="../util"
   DIA="../util/diag"
   PAR="../params"
   SLI="../science/sli"
   POP="../science/pop"
   /bin/cp -p $ALB/*90 ./.tmp
   /bin/cp -p $CAN/*90 ./.tmp
   /bin/cp -p $CNP/*90 ./.tmp
   /bin/cp -p $GWH/*90 ./.tmp
   /bin/cp -p $MIS/*90 ./.tmp
   /bin/cp -p $RAD/*90 ./.tmp
   /bin/cp -p $ROU/*90 ./.tmp
   /bin/cp -p $SOI/*90 ./.tmp
   /bin/cp -p $SLI/*90 ./.tmp
   /bin/cp -p $POP/*90 ./.tmp
   /bin/cp -p $OFF/*90 ./.tmp

/bin/cp -p $UTI/*90 ./.tmp
/bin/cp -p $DIA/*90 ./.tmp
/bin/cp -p $PAR/*90 ./.tmp

/bin/cp -p Makefile3_offline  ./.tmp
/bin/cp -p Makefile3_mpi  ./.tmp

cd .tmp/
if [[ $1 = 'mpi' ]]; then
   make -f Makefile3_mpi
else   
	 make -f Makefile3_offline
fi

}

build_status()
{
   if [[ -f .tmp/cable ]]; then
   	mv .tmp/cable .
   	print '\nBUILD OK\n'
   elif [[ -f .tmp/cable-mpi ]]; then
   	mv .tmp/cable-mpi .
   	print '\nBUILD OK\n'
   else
      print '\nOooops. Something went wrong\n'
      print '\nKnown build issues:\n'
      print '\nSome systems require additional library. \n'
      print '\nEdit Makefile_offline; add -lnetcdff to LD = ...\n'
   fi
   exit
}
###########################################
## build.ksh - MAIN SCRIPT STARTS HERE   ##
###########################################

if [[ $1 = 'clean' ]]; then
   clean_build
fi

if [[ $1 = 'mpi' ]]; then
  echo "Building cable_mpi" 
	host_gadi mpi
	#host_raij mpi
else
	host_gadi 
	#host_raij 
fi
