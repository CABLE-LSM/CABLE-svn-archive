#!/bin/csh 
# Lauren 23 Apr 13

set a=a
set rid=$RUNID
set pfrom=( $Pdaily )
set pto=( pb )

set plist=`ls $rid$a.$pfrom*.nc` 

foreach pfile ( $plist )

  # rename .$pfrom to .$pto
  set newname=`echo $pfile | sed -e s/\.$pfrom/\.$pto/`
  #echo $pfile $newname

if ( ! -e $newname ) then
  echo "Extracting " $newname " from " $pfile
 if ($RUNID == saaca) then
  ncks --mk_rec_dmn t_1 $pfile test.nc
  cdo selname,temp,temp_1 test.nc $newname
  rm test.nc
 else
 if ($RUNID == saaba || $RUNID == saadb) then
  cdo selname,temp_12,temp_13,temp_14 $pfile $newname
  ncrename -v temp_12,temp   $newname
  ncrename -v temp_13,temp_1 $newname
 else
  cdo selname,temp,temp_1 $pfile $newname
 endif
 endif

else

echo "File " $newname "Already Exists!"

endif
end

exit

