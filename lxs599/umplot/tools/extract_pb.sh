#!/bin/csh 
# Lauren 23 Apr 13

set a=a
set rid=$RUNID
set pfrom=( $Ptemp1 ) # pe
set pto=( $Ptemps )   # pb or pt

if ( $pfrom != $pto ) then

set plist=`ls $rid$a.$pfrom*.nc` 

if ( ${#plist} > 0 ) then
foreach pfile ( $plist )

  # rename .$pfrom to .$pto
  set newname=`echo $pfile | sed -e s/\.$pfrom/\.$pto/`
  #echo $pfile $newname

if ( ! -e $newname ) then

 echo "Extracting " $newname " from " $pfile

 #if ($RUNID == saaca) then
 # ncks --mk_rec_dmn t_1 $pfile test.nc
 # cdo selname,temp,temp_1 test.nc $newname
 # rm test.nc
 #else
 # if ($RUNID == saaba || $RUNID == saadb) then
 #  cdo selname,temp_12,temp_13,temp_14 $pfile $newname
 #  ncrename -v temp_12,temp   $newname
 #  ncrename -v temp_13,temp_1 $newname
 # else
 #  #cdo selname,temp,temp_1 $pfile $newname

   cdo selname,$txname,$tiname $pfile $newname
   ncrename -v $txname,tmax $newname
   ncrename -v $tiname,tmin $newname

 #  if ( $txname == "temp_1" ) then
 #   ncrename -v $txname,temp_3 $newname
 #   ncrename -v $tiname,temp_1 $newname
 #   ncrename -v temp_3,temp   $newname
 #  else
 #   if ($txname != "temp") then
 #    ncrename -v $txname,temp   $newname
 #   endif
 #   if ($tiname != "temp_1") then
 #    ncrename -v $tiname,temp_1 $newname
 #   endif
 #  endif #txname == "temp_1"
 # endif #rid=saaba
 #endif #rid=saaca

else  # ! newname
 echo "File " $newname "Already Exists!"
endif # ! newname

end   # foreach
# echo "No files"
endif # if plist

else

 echo "(1)     File Extensions are the Same"
 echo "(1)     pfrom = " $pfrom " and pto = " $pto
 set pbhead=`ls $rid$a.$pto?????.nc | head -1`
 set nvar=`cdo npar $pbhead`
 if ( $nvar <= 2 ) then
  echo "(1)     Script will rename variables in " $pto
  echo "(1)     If this is Incorrect, Try setting Ptemps to pt"
  ncrename -v $txname,tmax $newname
  ncrename -v $tiname,tmin $newname
 else
  echo "(1)     temps NOT extract from" $Ptemp1
  echo "(1)     Try setting Ptemps to pt"
 endif

endif

exit

