#!/bin/bash
# jtk561 - script to compile the code
module load intel-fc
module load intel-cc
module load netcdf/4.2.1.1
rm -f *.o
rm -f *.mod
rm -f run_prog
FFLAGS="-fpe0 -O0 -g -traceback -fp-model precise"
ifort -I${NETCDF}/include -L${NETCDF}/lib -lnetcdff -o run_prog $FFLAGS TypeDef.f90 GlobalDefs.f90 Constants.f90 DateFunctions.f90 Utils.f90 Pointer.f90 SubDiurnalMetModule.f90 GetSubDiurnalMetModule.f90 fileio_module.f90 run_prog.f90
