#!/bin/ksh

if [[ ! -d tmp ]]; then
   mkdir tmp
fi

if [[ -f cable ]]; then
   print '\ncable executable exists. copying to cable.bu\n' 
   mv cable cable.bu
fi

CORE="../core"
DRV="."
CASA="./CASA"

/bin/cp -p $CORE/*90 ./tmp
/bin/cp -p $DRV/*90 ./tmp
/bin/cp -p $CASA/*90 ./tmp

print '\nPlease note: CASA-CNP files are included in build only for technical reasons. Implementation is not officially available with the release of CABLE 2.0\n' 
/bin/cp -p Make_CABLE-offline  ./tmp

cd tmp/

make -f Make_CABLE-offline
   
if [[ -f cable ]]; then
   mv cable ../
fi

