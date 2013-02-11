#!/bin/ksh

if [[ ! -d tmp ]]; then
   mkdir tmp
fi

libpath="~/CABLE-UM/libcable.a"

if [[ -f $libpath ]]; then
   print '\nCABLE library exists at $libpath. Copying to $libpath.bu\n' 
   mv $libath $libpath.bu 
fi

CORE="../core"
DRV="."

/bin/cp -p $CORE/*90 ./tmp
/bin/cp -p $DRV/*90 ./tmp

#print '\nPlease note: CASA-CNP files are included in build only for technical reasons. Implementation is not officially available with the release of CABLE 2.0\n' 
/bin/cp -p Make_CABLE-UM ./tmp

cd tmp/

make -f Make_CABLE-UM
   
if [[ -f cable_explicit_driver.o ]]; then
   print '\nCompile appears successful. Now build library\n'
else
   print'\nCompile failed\n'
   exit
fi

## make library from CABLE object files

/usr/bin/ar r libcable.a cable_explicit_driver.o cable_implicit_driver.o   \
   cable_rad_driver.o cable_hyd_driver.o cable_common.o  \
   cable_define_types.o  \
   cable_soilsnow.o cable_air.o cable_albedo.o cable_radiation.o  \
   cable_roughness.o cable_carbon.o cable_canopy.o cable_cbm.o    \
   cable_um_tech.o cable_um_init_subrs.o cable_um_init.o 

if [[ -f libcable.a ]]; then
   print '\nLibrary build successful. Copying libcable.a to ~/CABLE-UM.\n'
else
   print '\nBuild failed\n'
   exit
fi

if [[ ! -d ~/CABLE-UM ]]; then
   mkdir ~/CABLE-UM
fi

/bin/cp -p libcable.a ~/CABLE-UM

if [[ -f ~/CABLE-UM/libcable.a ]]; then
   print "\nYour timestamped library shoulb thiis one:\n"
   echo `ls -alt ~/CABLE-UM/libcable.a`
   print '\nDONE\n'
   exit
else
   print '\nSomething went wrong!\n'
   exit
fi

