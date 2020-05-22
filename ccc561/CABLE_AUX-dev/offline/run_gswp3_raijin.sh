#!/bin/bash

# Run LIS
#PBS -m ae
#PBS -P w35
#PBS -l walltime=102800
#PBS -l mem=50GB
#PBS -l ncpus=28
#PBS -j oe
#PBS -q normalbw
#PBS -l wd
#PBS -l other=gdata1

module load intel-mpi/5.1.3.210
module load netcdf

ulimit -s unlimited

year=1901
while [ "$year" -le 2010 ]; do

   ./create_cable-nml.sh -y $year -e false -g true -a monthly
   mpirun -n 28 ./cable-mpi

   year=$(( $year + 1 ))

done

