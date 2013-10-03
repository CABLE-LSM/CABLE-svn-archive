#!/bin/csh
# Lauren 3 Oct 12

set a=a
set rid=$RUNID
set pfrom=( pf )
set pto=( pe )

set plist=`ls $rid$a.$pfrom*.nc`

foreach pfile ( $plist )
  # rename .$pfrom to .$pto
  set newname=`echo $pfile | sed -e s/\.$pfrom/\.$pto/`
  echo "Moving " $pfile " to " $newname
  mv $pfile $newname
end

exit

