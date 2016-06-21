#!/bin/csh -x

if ( $HOSTNAME == ruby ) then
 set FLDSUBSET = ~dix043/src/python/um/um_fields_subset.py
else
 set FLDSUBSET = ~mrd599/src/python/um/um_fields_subset.py
 set CONV2NC   = ~lxs599/umutils/conv2nc.tcl
endif

set a=a
set rid=$RUNID
set ext=$1 #pf
if ( $CPL == y ) then
 set flist=`ls $rid.$ext-??????????`
else
 set flist=`ls $rid$a.$ext?????`
 if (${#flist} == 0) then
  set flist=`ls $rid$a.$ext???????`
  #set flist=`ls $rid$a.$ext?????`
 endif
endif

foreach file ( $flist )
 if ( ! -e $file.nc ) then
  #if ( $ext == pe ) then
  # python $FLDSUBSET -i $file -o $file.sub -v 3236
  # echo "Converting " $file " to netCDF: " $file.nc
  # $CONV2NC -i $file.sub -o $file.nc
  # rm $file.sub
  #else
   echo "Converting " $file " to netCDF: " $file.nc
   $CONV2NC -i $file -o $file.nc
  #endif
 endif
end

exit

