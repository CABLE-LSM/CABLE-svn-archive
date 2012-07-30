#!/bin/ksh

. /apps/modules/Modules/default/init/ksh

if [[ ! -d mpitmp ]]; then
   mkdir mpitmp
fi

if [[ -f cable-mpi ]]; then
   print '\ncable-mpi executable exists. copying to cable-mpi.bu\n' 
   mv cable-mpi cable-mpi.bu
fi

CORE="../core/biogeophys"
DRV="."
CASA="../core/biogeochem"

/bin/cp -p $CORE/*90 ./mpitmp
/bin/cp -p $DRV/*90 ./mpitmp
/bin/cp -p $CASA/*90 ./mpitmp

print '\nPlease note: CASA-CNP files are included in build only for technical reasons. Implementation is not officially available with the release of CABLE 2.0\n' 
/bin/cp -p Makefile_mpi  ./mpitmp

cd mpitmp/

module add netcdf/3.6.3 openmpi
make -f Makefile_mpi
   
if [[ -f cable-mpi ]]; then
   cp cable-mpi ../
	print '\nBUILD OK\n'
	print '\nThe following executable has been copied to your ~/CABLE_AUX/run/ directory.\n'
	echo `ls -lh cable-mpi`
	mv cable-mpi ~/CABLE-AUX/run/
fi

