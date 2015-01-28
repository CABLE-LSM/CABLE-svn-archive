#!/bin/csh -x

if ( $HOSTNAME == cherax ) then
 set FLDSUBSET = ~dix043/src/python/um/um_fields_subset.py
else
 set FLDSUBSET = ~mrd599/src/python/um/um_fields_subset.py
endif

set a=a
set rid=$RUNID
set ext=$1 #pf
if ( $CPL == y ) then
 set flist=`ls $rid.$ext-??????????`
else
 set flist=`ls $rid$a.$ext?????`
endif

foreach file ( $flist )
 if ( ! -e $file.nc ) then
  if ( $ext == pe ) then
   python $FLDSUBSET -i $file -o $file.sub -v 3236
   echo "Converting " $file " to netCDF: " $file.nc
   $CONV2NC -i $file.sub -o $file.nc
   rm $file.sub
  else
   echo "Converting " $file " to netCDF: " $file.nc
   $CONV2NC -i $file -o $file.nc
  endif
 endif
end

exit

