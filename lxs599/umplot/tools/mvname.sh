#!/bin/csh 

set a=a
set mvf=$1   #saabb # RUNID 1
set mvt=$2   #saaba # RUNID 2
set ext=pf   # if only one file extension

set filelist=`ls $mvf$a.p*`
#set filelist=`ls $mvf$a.$ext*`
#set filelist=`ls $mvf$a.pfk[012345678]*`

foreach pfile ( $filelist )
 set newname=`echo $pfile | sed -e s/$mvf/$mvt/`
 if ( ! -e $newname ) then
  echo "Moving " $pfile " to " $newname
  mv $pfile $newname
 endif
end

exit

