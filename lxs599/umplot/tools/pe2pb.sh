#!/bin/csh
# Lauren 3 Oct 12

set a=a
set rid=$RUNID
set pfrom=( $Ptemp1 ) #pe
set pto=( $Ptemps )   #pb

set plist=`ls $rid$a.$pfrom*.nc`

foreach pfile ( $plist )
  # rename .$pfrom to .$pto
  set newname=`echo $pfile | sed -e s/\.$pfrom/\.$pto/`
  if ( ! -e $newname ) then
   echo "Moving " $pfile " to " $newname
   mv $pfile $newname
  else
   echo "File " $newname "Already Exists!"
  endif
end

exit

