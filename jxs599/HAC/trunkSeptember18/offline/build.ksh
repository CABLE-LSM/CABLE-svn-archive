#!/bin/ksh

## raijin.nci.org.au
host_raij()
{
   module add intel-fc
   module add netcdf

   export NCDIR=$NETCDF_ROOT'/lib/Intel'
   export NCMOD=$NETCDF_ROOT'/include/Intel'
   export FC=$F90
   export CFLAGS='-O0 -fp-model precise'
   if [[ $1 = 'debug' ]]; then
      export CFLAGS='-O0 -traceback -g -fp-model precise -ftz -fpe0'
   fi
   export LDFLAGS='-L'$NCDIR' -O0'
   export LD='-lnetcdf -lnetcdff'
   build_build
   cd ../
}


clean_build()
{
      print '\ncleaning up\n'
      print '\n\tPress Enter too continue buiding, Control-C to abort now.\n'
      read dummy
      rm -fr .tmp
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
CBM="../science/cbl_model_driver_offline.F90"
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
#FRO="../util/FromHAC"
PAR="../params"
SLI="../science/sli"
POP="../science/pop"
/bin/cp -p $ALB/*90 ./.tmp
/bin/cp -p $CAN/*90 ./.tmp
/bin/cp -p $CNP/*90 ./.tmp
/bin/cp -p $CBM     ./.tmp
/bin/cp -p $GWH/*90 ./.tmp
/bin/cp -p $MIS/*90 ./.tmp
/bin/cp -p $RAD/*90 ./.tmp
/bin/cp -p $RAD/*inc ./.tmp
/bin/cp -p $ROU/*90 ./.tmp
/bin/cp -p $SOI/*90 ./.tmp
#/bin/cp -p $SUR/*90 ./.tmp
/bin/cp -p $SLI/*90 ./.tmp
/bin/cp -p $POP/*90 ./.tmp
/bin/cp -p $OFF/*90 ./.tmp

/bin/cp -p $UTI/*90 ./.tmp
/bin/cp -p $DIA/*90 ./.tmp
#/bin/cp -p $FRO/*90 ./.tmp
/bin/cp -p $PAR/*90 ./.tmp

/bin/cp -p Makefile_offline  ./.tmp

cd .tmp/

make -f Makefile_offline

mv cable ../
}

###########################################
## build.ksh - MAIN SCRIPT STARTS HERE   ##
###########################################

if [[ $1 = 'clean' ]]; then
   clean_build
fi

module purge

host_raij
