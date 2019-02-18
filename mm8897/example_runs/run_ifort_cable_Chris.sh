#!/bin/bash
version=$1
rm cable.nml
ln -s "Chris_cable.nml" "cable.nml"
source /home/mm/intel/bin/compilervars.sh -arch intel64
NCDIR=/home/mm/local/lib
NCMOD=/home/mm/local/include
LD_LIBRARY_PATH=${NCDIR}:${LD_LIBRARY_PATH} 
../cable2.0-trunk/offline/ifort_cable

