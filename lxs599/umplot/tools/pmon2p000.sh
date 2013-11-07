#!/bin/csh

set a=a
set rid=$RUNID

# Just change the EXT letter i.e. pb,pc,pm etc
#set pext = ( pb pc pe pf pg )
set pext = ( p )

foreach ext ( $pext )

set i=1
set plist=`ls $rid$a.$ext*.nc`

# change months to numbers in .p
#set newlist=`echo $plist | sed -e 's/jan/010/g'`
#set newlist=`echo $newlist | sed -e 's/feb/020/g'`
#set newlist=`echo $newlist | sed -e 's/mar/030/g'`
#set newlist=`echo $newlist | sed -e 's/apr/040/g'`
#set newlist=`echo $newlist | sed -e 's/may/050/g'`
#set newlist=`echo $newlist | sed -e 's/jun/060/g'`
#set newlist=`echo $newlist | sed -e 's/jul/070/g'`
#set newlist=`echo $newlist | sed -e 's/aug/080/g'`
#set newlist=`echo $newlist | sed -e 's/sep/090/g'`
#set newlist=`echo $newlist | sed -e 's/oct/100/g'`
#set newlist=`echo $newlist | sed -e 's/nov/110/g'`
#set newlist=`echo $newlist | sed -e 's/dec/120/g'`

set newlist=`echo $plist | sed -e 's/jan/001/g'`
set newlist=`echo $newlist | sed -e 's/feb/002/g'`
set newlist=`echo $newlist | sed -e 's/mar/003/g'`
set newlist=`echo $newlist | sed -e 's/apr/004/g'`
set newlist=`echo $newlist | sed -e 's/may/005/g'`
set newlist=`echo $newlist | sed -e 's/jun/006/g'`
set newlist=`echo $newlist | sed -e 's/jul/007/g'`
set newlist=`echo $newlist | sed -e 's/aug/008/g'`
set newlist=`echo $newlist | sed -e 's/sep/009/g'`
set newlist=`echo $newlist | sed -e 's/oct/010/g'`
set newlist=`echo $newlist | sed -e 's/nov/011/g'`
set newlist=`echo $newlist | sed -e 's/dec/012/g'`

foreach pfile ( $plist )
  if ( ! -e $newlist[$i] ) then
   echo "Moving " $pfile " to " $newlist[$i]
   mv $pfile $newlist[$i]
  else
   echo "File" $newlist[$i] "Already Exists!"
  endif
  @ i++
end

end

exit

