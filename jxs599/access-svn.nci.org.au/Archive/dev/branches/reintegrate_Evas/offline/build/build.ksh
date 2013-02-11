#!/bin/ksh

if [[ ! -d tmp ]]; then
   mkdir tmp
fi

COM="../../core/src/common/"
SCI="../../core/src/science/"
DRV="../../offline/src/"
CNP="../../core/pkg/CASA/"

/bin/cp --preserve $COM/*90 ./tmp
/bin/cp --preserve $DRV/*90 ./tmp
/bin/cp --preserve $SCI/*90 ./tmp
/bin/cp --preserve $CNP/*90 ./tmp

/bin/cp Make_CABLE-offline --preserve ./tmp

cd tmp/
   make -f Make_CABLE-offline
   
if [[ -f cable ]]; then
   mv cable ../
fi

