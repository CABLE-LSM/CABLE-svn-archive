#!/bin/csh
# Lauren 3 Oct 12

set a=a
set rid=$RUNID
set pfrom=( $1 ) #pa #$Pmonth
set pto=( $2 )   #pm

if ( $pfrom != $pto ) then

set plist=`ls $rid$a.$pfrom*.nc`

if ( ${#plist} > 0 ) then
 foreach pfile ( $plist )
  # rename .$pfrom to .$pto
  set newname=`echo $pfile | sed -e s/\.$pfrom/\.$pto/`
  if (! -e $newname) then
   echo "Moving " $pfile " to " $newname
   mv $pfile $newname
  else
   echo "File " $newname "Already Exists!"
  endif
 end
endif

else
 echo "(1)     Error: File Extensions are the Same"
 echo "(1)     pfrom = " $pfrom " and pto = " $pto
endif

exit

