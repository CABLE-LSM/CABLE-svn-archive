#!/bin/csh

set a=a
set rid=$RUNID

# Just change the EXT letter i.e. pb,pc,pm etc
set pext = ( e )

foreach ext ( $pext )

set i=1
set plist=`ls $rid$a.p$ext*.nc`

# change numbers to months in .p
set newlist=`echo $plist | sed -e 's/010/jan/g'`
set newlist=`echo $newlist | sed -e 's/020/feb/g'`
set newlist=`echo $newlist | sed -e 's/030/mar/g'`
set newlist=`echo $newlist | sed -e 's/040/apr/g'`
set newlist=`echo $newlist | sed -e 's/050/may/g'`
set newlist=`echo $newlist | sed -e 's/060/jun/g'`
set newlist=`echo $newlist | sed -e 's/070/jul/g'`
set newlist=`echo $newlist | sed -e 's/080/aug/g'`
set newlist=`echo $newlist | sed -e 's/090/sep/g'`
set newlist=`echo $newlist | sed -e 's/100/oct/g'`
set newlist=`echo $newlist | sed -e 's/110/nov/g'`
set newlist=`echo $newlist | sed -e 's/120/dec/g'`

foreach pfile ( $plist )
  #echo $pfile $newlist[$i]
  mv $pfile $newlist[$i]
  @ i++
end

end

exit

