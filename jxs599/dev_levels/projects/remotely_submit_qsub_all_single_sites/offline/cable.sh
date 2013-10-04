#!/bin/csh  

#PBS -q express
#PBS -l vmem=1000mb
#PBS -l walltime=2000
cd $CABLE_SHINE/bCABLE-1.9/offline/build/
make -f Make_CABLE-offline
mv cable $CABLE_SHINE/bCABLE-1.9/offline/test_offline/
cd $CABLE_SHINE/bCABLE-1.9/offline/test_offline/
./run_qs.csh
