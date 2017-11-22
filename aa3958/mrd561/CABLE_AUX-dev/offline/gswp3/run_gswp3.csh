#!/bin/csh

# Run LIS
#PBS -m ae
#PBS -P w35
#PBS -l walltime=45000
#PBS -l mem=7500MB
#PBS -l ncpus=16
#PBS -j oe
#PBS -q express
#PBS -l wd
#PBS -l other=gdata1

module load intel-mpi/4.1.1.036
module load netcdf/4.2.1.1

set year = 1901
while ( $year <= 2010 )
   sh ./create_cable-nml.sh -y $year -e false -g true -a monthly
   mpirun -n 16 ./cable-mpi
   @ year++
end

