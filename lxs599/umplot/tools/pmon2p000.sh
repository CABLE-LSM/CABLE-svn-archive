#!/bin/csh

set a=a
set rid=$RUNID

# Just change the EXT letter i.e. pb,pc,pm etc
set pext = ( e )
#set pext = ( b c f g )

foreach ext ( $pext )

set i=1
set plist=`ls $rid$a.p$ext*.nc`

# change months to numbers in .p
set newlist=`echo $plist | sed -e 's/jan/010/g'`
set newlist=`echo $newlist | sed -e 's/feb/020/g'`
set newlist=`echo $newlist | sed -e 's/mar/030/g'`
set newlist=`echo $newlist | sed -e 's/apr/040/g'`
set newlist=`echo $newlist | sed -e 's/may/050/g'`
set newlist=`echo $newlist | sed -e 's/jun/060/g'`
set newlist=`echo $newlist | sed -e 's/jul/070/g'`
set newlist=`echo $newlist | sed -e 's/aug/080/g'`
set newlist=`echo $newlist | sed -e 's/sep/090/g'`
set newlist=`echo $newlist | sed -e 's/oct/100/g'`
set newlist=`echo $newlist | sed -e 's/nov/110/g'`
set newlist=`echo $newlist | sed -e 's/dec/120/g'`

foreach pfile ( $plist )
  echo $pfile $newlist[$i]
  mv $pfile $newlist[$i]
  @ i++
end

end

exit

