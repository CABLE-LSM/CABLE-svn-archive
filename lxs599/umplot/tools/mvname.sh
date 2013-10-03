#!/bin/csh 

set a=a
set mvf=saabb # RUNID 1
set mvt=saaba # RUNID 2
set ext=pf    # if only one file extension

set filelist=`ls $mvf$a.p*`
#set filelist=`ls $mvf$a.$ext*`
#set filelist=`ls $mvf$a.pfk[012345678]*`

foreach pfile ( $filelist )
 set newname=`echo $pfile | sed -e s/$mvf/$mvt/`
 echo "Moving " $pfile " to " $newname
 mv $pfile $newname
end

exit

