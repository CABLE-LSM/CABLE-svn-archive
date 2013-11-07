#!/bin/csh
# Lauren 6 Sept 12 - note: might rearrange - rename first then see if need convert to nc.

#module load python
#module load netcdf
module load nco

set i=1
set a=a
set rid=$RUNID

#set mext=$Pmonth #pa # monthly files    e.g. pm
#set dext=$Pdaily #pe # daily files      e.g. pa
#set 1ext=$Ptemp1 #pe # from: max/min T files e.g. pb
#set bext=$Ptemps #pb #   to: max/min T files e.g. pb
#set text=$Ptimes #pf # timeseries files e.g. pe
#set jext=$PcasaC #pj # carbon files e.g. pj

set palist=`ls $rid$a.pa?????` # pm (monthly values)
set pblist=`ls $rid$a.pc?????` # pb (daily tmin, tmax)
set pelist=`ls $rid$a.pe?????` # pa (daily values)
set pflist=`ls $rid$a.pf?????` # pe (timeseries)
set pglist=`ls $rid$a.pg?????` # pg (timeseries)
set pjlist=`ls $rid$a.pj?????` # pj (monthly values)

## change months to numbers in .pa
#set newlist=`echo $palist | sed -e 's/jan/010/g'`
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

foreach afile ( $palist )
  if (! -e $afile.nc) then
#  if (! -e $newlist[$i].nc) then
     ~ste69f/umutils/conv2nc.tcl -i $afile -o $afile.nc
#     ~ste69f/umutils/conv2nc.tcl -i $afile -o $newlist[$i].nc
  endif
  # rename .pa to .pm
  set newname=`echo $afile.nc | sed -e 's/\.pa/\.pm/'`
#  set newname=`echo $newlist[$i].nc | sed -e 's/\.pa/\.pm/'`
  mv $afile.nc $newname
#  mv $newlist[$i].nc $newname
  # swap fields; field322 ad field322_1 (sea ice and land albedo)
#  ncrename -v field322_1,field322_2 $newname
#  ncrename -v field322,field322_1 $newname
#  ncrename -v field322_2,field322 $newname
#  @ i++
end

# change months to numbers in .pb
set newlist=`echo $pblist | sed -e 's/jan/010/g'`
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

#set i=1

foreach bfile ( $pblist )
  if (! -e $newlist[$i].nc) then
     ~ste69f/umutils/conv2nc.tcl -i $bfile -o $newlist[$i].nc
  endif
  # rename .pc to .pb  
  set newname=`echo $newlist[$i].nc | sed -e 's/\.pc/\.pb/'`
  mv $newlist[$i].nc $newname
  # swap fields; field temp and field temp_1 (min and max temp)
#  ncrename -v temp_1,temp_2 $newlist[$i].nc
#  ncrename -v temp,temp_1 $newlist[$i].nc
#  ncrename -v temp_2,temp $newlist[$i].nc
  @ i++
end

# change months to numbers in .pe file (only required for casa spinup)
set newlist=`echo $pelist | sed -e 's/jan/010/g'`
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

set i=1

foreach efile ( $pelist )
  if (! -e $newlist[$i].nc) then
     ~ste69f/umutils/conv2nc.tcl -i $efile -o $newlist[$i].nc
  endif
  # rename .pe to .pa
  set newname=`echo $newlist[$i].nc | sed -e 's/\.pe/\.pa/'`
  mv $newlist[$i].nc $newname
  @ i++
end

# change months to numbers in .pf file
set newlist=`echo $pflist | sed -e 's/jan/010/g'`
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

set i=1

foreach ffile ( $pflist )
  if (! -e $newlist[$i].nc) then
     ~ste69f/umutils/conv2nc.tcl -i $ffile -o $newlist[$i].nc
     #python ~dix043/src/python/um/um_timeseries.py -i $ffile -o $newlist[$i].nc
  endif
  # rename .pf to .pe
  set newname=`echo $newlist[$i].nc | sed -e 's/\.pf/\.pe/'`
  mv $newlist[$i].nc $newname
  @ i++
end

# change months to numbers in .pg file
set newlist=`echo $pglist | sed -e 's/jan/010/g'`
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

set i=1

foreach gfile ( $pglist )
  if (! -e $newlist[$i].nc) then
     ~ste69f/umutils/conv2nc.tcl -i $gfile -o $newlist[$i].nc
     #python ~dix043/src/python/um/um_timeseries.py -i $gfile -o $newlist[$i].nc
  endif
  @ i++
end

# change months to numbers in .pj
set newlist=`echo $pjlist | sed -e 's/jan/010/g'`
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

set i=1

foreach jfile ( $pjlist )
  if (! -e $newlist[$i].nc) then
     ~ste69f/umutils/conv2nc.tcl -i $jfile -o $newlist[$i].nc
  endif
  @ i++
end


exit

