#!/bin/bash
#SBATCH --time=29:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=120G
#SBATCH --exclusive
#SBATCH --qos=express

# --time= - about 16 minutes per year/spin of simulation on pearcey with np=20

# remove whatever versions of these modules are in the environment
module del netcdf openmpi intel-cc intel-fc
# load the safest default versions
module add intel-cc intel-fc netcdf openmpi
# for reference see what's in there
module list

cd /scratch1/gol124/support/CABLE/run
ulimit -s unlimited
# to enable coredumps
ulimit -c unlimited
ulimit -a

cp /datastore/wan028/pearcey-data/run-CABLE-2.0_mpi_bp2971/nml-trendy-runspin/cable_CNP_1901.nml cable.nml
cp /datastore/wan028/pearcey-data/run-CABLE-2.0_mpi_bp2971/trendy/poolcnp_out1860.csv poolcnp_in.csv
cp /datastore/wan028/run_cable2.0_mpi/trendy/S1/c/varyingLAI/restart_out_1901_S1.nc restart_in.nc

i=1 

while [ $i -le 1 ]
do
    yr=1901
    echo "Starting year $yr spin $i"
    time mpirun -v -np 20 ./cable-mpi
    cp -p restart_out.nc restart_in.nc
    cp -p restart_out.nc  output/restart_spin${i}_${yr}_cn.nc
    cp -p poolcnp_out.csv poolcnp_in.csv
    mv poolcnp_out.csv   output/poolcnp_out${i}_${yr}_cn.csv
    mv cnpflux${yr}.csv  output/cnpflux${i}_${yr}_cn.csv
    mv out_cable.nc      output/outcable${i}_${yr}_cn.nc
    i=`expr $i + 1`
done

