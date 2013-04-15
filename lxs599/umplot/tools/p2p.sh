#!/bin/csh -x
# Lauren 3 Oct 12

set a=a
set rid=$RUNID
set pfrom=( pb )
set pto=( pm )

set plist=`ls $rid$a.$pfrom*.nc` # pm (monthly values)

foreach pfile ( $plist )
  # rename .$pfrom to .$pto
  set newname=`echo $pfile | sed -e s/\.$pfrom/\.$pto/`
  #echo $pfile $newname
  mv $pfile $newname
end

exit

