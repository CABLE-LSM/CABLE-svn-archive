#!/bin/bash

# Run LIS
#PBS -m ae
#PBS -P w35
#PBS -l walltime=172800
#PBS -l mem=9000MB
#PBS -l ncpus=28
#PBS -j oe
#PBS -q normalbw
#PBS -l wd
#PBS -l other=gdata1

module load intel-mpi/5.1.3.210
module load netcdf

ulimit -s unlimited



#set data dirs
if [ ! -h surface_data ]; then
   ln -s /g/data1/w35/mrd561//CABLE/CABLE_AUX-dev/offline surface_data
fi

if [ ! -h gswp ]; then
   ln -s /g/data1/wd9/MetForcing/Global/GSWP3/Corr.20C2.05 gswp
fi

year=1901
while [ "$year" -le 2010 ]; do

   sh ./create_cable-nml.sh -y $year -e false -g true -a monthly
   mpirun -n 28 ./cable-mpi
   year=$(( $year + 1 ))

done

