#!/bin/bash

export NCLIB='/share/apps/netcdf/intel/4.1.3/lib'
export NCDIR='/share/apps/netcdf/intel/4.1.3/lib'
export NCMOD='/share/apps/netcdf/intel/4.1.3/include'
export FC=ifort
export CFLAGS='-O0 -fp-model precise -fpe0 -g -traceback -nostand -check all,nobounds,noarg_temp_created -debug all -shared -fpic ' 
export LDFLAGS='-L/share/apps/intel/Composer/lib/intel64 -L/share/apps/netcdf/intel/4.1.3/lib '
export LD='-lnetcdf -lnetcdff'


mkdir .mytmp

tmp_build=$(pwd)/.mytmp


cp *90 .mytmp/
cp ../core/biogeophys/*90 .mytmp/.
cp ../core/biogeochem/*90 .mytmp/.

cp Makefile_share .mytmp/.

cd .mytmp

make -f Makefile_share

ifort -O0 -fp-model precise -fpe0 -g -traceback -nostand -check all,nobounds,noarg_temp_created -debug all -shared -fpic $LDFLAGS -L${NCLIB} -o cable_python_driver.so cable_python_driver.F90 *.so -lnetcdf -lnetcdff




echo '  now please run :   export LD_LIBRARY_PATH=$(pwd)/.mytmp:$LD_LIBRARY_PATH '


