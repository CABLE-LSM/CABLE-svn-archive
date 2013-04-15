#!/bin/csh -x

set a=a
# ./pm2num.sh
set pmfile=`ls $DIR/$RUNID$a.pm*0.nc`
#set pafile=`ls $DIR/$RUNID$a.pa*0`

foreach mfile ( $pmfile )

# cdo showname $mfile.nc
cdo selname,precip,$tname,temp,sh,lh,solar,longwave,ilr,field203,field208,field1526,field1527,field1894,field1895,field1893,field322,snowdepth,temp_1,field202_1,field141,field1532,field1533,field1386,sm_1,sm,field1388,field1389,field1390,field1523,soiltemp $mfile.nc $mfile.sub.nc

end

exit

