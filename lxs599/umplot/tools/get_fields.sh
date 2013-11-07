#!/bin/csh -x

set a=a
set ext = ( $1 ) #( pm )
set fld = ( $2 ) # precip,$tname,temp,sh,lh,solar,longwave...
#precip,$tname,temp,sh,lh,solar,longwave,ilr,field203,field208,field1526,field1527,field1894,field1895,field1893,field322,snowdepth,temp_1,field202_1,field141,field1532,field1533,field1386,sm_1,sm,field1388,field1389,field1390,field1523,soiltemp

# ./pmon2p000.sh # ./pm2num.sh
#set pfile=`ls $DIR/$RUNID$a.pa*0.nc`

set pfile=`ls $DIR/$RUNID$a.$ext?????.nc`

foreach file ( $pfile )

# cdo showname $file.nc
 cdo selname,$fld $file.nc $file.sub.nc

# cdo selname,precip,$tname,temp,sh,lh,solar,longwave,ilr,field203,field208,field1526,field1527,field1894,field1895,field1893,field322,snowdepth,temp_1,field202_1,field141,field1532,field1533,field1386,sm_1,sm,field1388,field1389,field1390,field1523,soiltemp $file.nc $file.sub.nc

end

exit

