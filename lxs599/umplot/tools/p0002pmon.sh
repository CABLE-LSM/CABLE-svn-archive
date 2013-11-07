#!/bin/csh

set a=a
set rid=$RUNID

# Just change the EXT letter i.e. pb,pc,pm etc
set pext = ( pe pb )

foreach ext ( $pext )

set i=1
set plist=`ls $rid$a.$ext*.nc`

# change numbers to months in .p
#set newlist=`echo $plist | sed -e 's/010/jan/g'`
#set newlist=`echo $newlist | sed -e 's/020/feb/g'`
#set newlist=`echo $newlist | sed -e 's/030/mar/g'`
#set newlist=`echo $newlist | sed -e 's/040/apr/g'`
#set newlist=`echo $newlist | sed -e 's/050/may/g'`
#set newlist=`echo $newlist | sed -e 's/060/jun/g'`
#set newlist=`echo $newlist | sed -e 's/070/jul/g'`
#set newlist=`echo $newlist | sed -e 's/080/aug/g'`
#set newlist=`echo $newlist | sed -e 's/090/sep/g'`
#set newlist=`echo $newlist | sed -e 's/100/oct/g'`
#set newlist=`echo $newlist | sed -e 's/110/nov/g'`
#set newlist=`echo $newlist | sed -e 's/120/dec/g'`

set newlist=`echo $plist | sed -e 's/001/jan/g'`
set newlist=`echo $newlist | sed -e 's/002/feb/g'`
set newlist=`echo $newlist | sed -e 's/003/mar/g'`
set newlist=`echo $newlist | sed -e 's/004/apr/g'`
set newlist=`echo $newlist | sed -e 's/005/may/g'`
set newlist=`echo $newlist | sed -e 's/006/jun/g'`
set newlist=`echo $newlist | sed -e 's/007/jul/g'`
set newlist=`echo $newlist | sed -e 's/008/aug/g'`
set newlist=`echo $newlist | sed -e 's/009/sep/g'`
set newlist=`echo $newlist | sed -e 's/010/oct/g'`
set newlist=`echo $newlist | sed -e 's/011/nov/g'`
set newlist=`echo $newlist | sed -e 's/012/dec/g'`

foreach pfile ( $plist )
  if ( ! -e $newlist[$i] ) then
   echo "Moving " $pfile " to " $newlist[$i]
   mv $pfile $newlist[$i]
  endif
  @ i++
end

end

exit

