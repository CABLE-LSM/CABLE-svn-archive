#!/bin/ksh

if [[ ! -d tmp ]]; then
   mkdir tmp
fi

COM="../../core/src/common/"
SCI="../../core/src/science/"
DRV="../../offline/src/"
CNP="../../core/pkg/CASA/"

if [[ `uname -s` == "Darwin" ]]; then
   /bin/cp -p $COM/*90 ./tmp
   /bin/cp -p $DRV/*90 ./tmp
   /bin/cp -p $SCI/*90 ./tmp
   /bin/cp -p $CNP/*90 ./tmp
   /bin/cp -p Make_CABLE-offline  ./tmp
else
   /bin/cp --preserve $COM/*90 ./tmp
   /bin/cp --preserve $DRV/*90 ./tmp
   /bin/cp --preserve $SCI/*90 ./tmp
   /bin/cp --preserve $CNP/*90 ./tmp
   /bin/cp Make_CABLE-offline --preserve ./tmp
fi

cd tmp/
   make -f Make_CABLE-offline
   
if [[ -f cable ]]; then
   mv cable ../
fi

