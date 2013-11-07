#!/bin/csh

set a=a
set rid=$RUNID
set ext=$1 #pf
set flist=`ls $rid$a.$ext?????`

foreach file ( $flist )
 if ( ! -e $file.nc ) then
  echo "Converting " $file " to netCDF: " $file.nc
  $CONV2NC -i $file -o $file.nc
 endif
end

exit

