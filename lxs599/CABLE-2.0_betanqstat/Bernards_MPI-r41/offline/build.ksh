#!/bin/ksh

. /apps/modules/Modules/default/init/ksh

if [[ ! -d tmp ]]; then
   mkdir tmp
fi

if [[ -f cable ]]; then
   print '\ncable executable exists. copying to cable.bu\n' 
   mv cable cable.bu
fi

CORE="../core/biogeophys"
DRV="."
CASA="../core/biogeochem"

/bin/cp -p $CORE/*90 ./tmp
/bin/cp -p $DRV/*90 ./tmp
/bin/cp -p $CASA/*90 ./tmp

print '\nPlease note: CASA-CNP files are included in build only for technical reasons. Implementation is not officially available with the release of CABLE 2.0\n' 
/bin/cp -p Makefile_offline  ./tmp

cd tmp/

module add netcdf/3.6.3 
make -f Makefile_offline
   
if [[ -f cable ]]; then
   cp cable ../
	print '\nBUILD OK\n'
	print '\nThe following executable has been copied to your ~/CABLE_AUX/run/ directory.\n'
	echo `ls -lh cable`
	mv cable ~/CABLE-AUX/run/
fi

