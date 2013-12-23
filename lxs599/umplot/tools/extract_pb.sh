#!/bin/csh 
# Lauren 23 Apr 13

set a=a
set rid=$RUNID
set pfrom=( $Ptemp1 ) # pe
set pto=( $Ptemps )   # pb or pt

set plist=`ls $rid$a.$pfrom*.nc`

if ( $pfrom != $pto ) then

 if ( ${#plist} > 0 ) then
  foreach pfile ( $plist )
   # rename .$pfrom to .$pto
   set newname=`echo $pfile | sed -e s/\.$pfrom/\.$pto/`
   if ( ! -e $newname ) then
    echo "Extracting " $newname " from " $pfile
    # ncks --mk_rec_dmn t_1 $pfile test.nc
    cdo selname,$txname,$tiname $pfile $newname
    ncrename -v $txname,tmax $newname
    ncrename -v $tiname,tmin $newname
   else  # ! newname
    echo "(1)     File " $newname "Already Exists!"
   endif # ! newname
  end   # foreach
   # echo "No files"
 endif # if plist

else # pfrom pto

 echo "(1)     File Extensions are the Same"
 echo "(1)     pfrom = " $pfrom " and pto = " $pto
 if ( ${#plist} > 0 ) then
  set pbhead=`ls $rid$a.$pto?????.nc | head -1`
  set nvar=`cdo npar $pbhead`
  if ( $nvar <= 2 ) then
   echo "(1)     Script will rename variables in " $pto
   echo "(1)     If this is Incorrect, Try setting Ptemps to pt"
   foreach pfile ( $plist )
    ncrename -v $txname,tmax $pfile
    ncrename -v $tiname,tmin $pfile
   end
  else # nvar
   echo "(1)     temps NOT extract from" $Ptemp1
   echo "(1)     Try setting Ptemps to pt"
  endif # nvar
 endif # plist

endif # pfrom pto

exit

