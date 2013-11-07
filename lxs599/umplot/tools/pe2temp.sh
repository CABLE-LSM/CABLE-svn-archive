#!/bin/csh -x

# set YR before running
set rid=$RUNID # uaawa
set a=a
set ext=( pe )

set filelist=`ls $rid$a.$ext*.nc`

foreach file ( $filelist )
 ncks -v $txname,$tiname $file $file.sub
 #ncks -v temp_12,temp_13 $file $file.sub
end

set yr    = (h8 h9 i0 i1 i2 i3)
#set yr    = (h8 h9 i0 i1 i2 i3 i4 i5 i6 i7 i8 i9 j0 j1 j2 j3 j4 j5 j6 j7 j8 j9 k0)
set mon   = (jan feb mar apr may jun jul aug sep oct nov dec)
set count = (01 02 03 04 05 06 07 08 09 10 11 12)
set i = 1

foreach yri ( $yr )
 foreach moni ( $mon )
  if ( -e $rid$a.$ext${yri}${moni}.nc.sub ) then
   mv $rid$a.$ext${yri}${moni}.nc.sub $rid$a.pb${yri}${count[$i]}.nc
   cdo monmean $rid$a.pb${yri}${count[$i]}.nc pb${yri}${count[$i]}.nc
  endif
@ i = $i + 1
 end
@ i = 1
end

@ numYr = ($YR * 12)
set pblist=`ls $rid$a.pb*.nc | head -$numYr`
cdo copy $pblist Tempseries_${YR}yrs.nc
cdo yseasavg Tempseries_${YR}yrs.nc Tseasonal_means_${YR}yrs.nc
cdo ymonavg Tempseries_${YR}yrs.nc Tmonthly_means_${YR}yrs.nc

exit

