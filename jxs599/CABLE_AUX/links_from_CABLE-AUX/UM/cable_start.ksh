#!/bin/ksh

. ~/.kshrc

if [[ ! -d $CABLE_AUX/CABLE-AUX ]]; then
   mkdir $CABLE_AUX/CABLE-AUX
   cp -r /projects/access/CABLE-AUX/UM $CABLE_AUX/CABLE-AUX
   cp -r /projects/access/CABLE-AUX/core $CABLE_AUX/CABLE-AUX
fi

if [[ ! -d /short/p66/`whoami`/$RUNID ]]; then
   mkdir /short/p66/`whoami`/$RUNID
fi
if [[ ! -f /short/p66/`whoami`/$RUNID/cable.nml ]]; then
   cp $CABLE_AUX/CABLE-AUX/UM/cable.nml /short/p66/`whoami`/$RUNID
fi


