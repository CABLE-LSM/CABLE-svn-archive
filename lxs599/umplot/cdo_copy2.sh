#!/bin/csh

set a=a
set dir=$DIR
set rid=$RUNID
set decade=(h i j k)
set year=(0 1 2 3 4 5 6 7 8 9)

dmget $DIR/$rid$a.pm*.nc

#set count_all=`ls h[0123456789]*.nc i[0123456789]*.nc j[0123456789]*.nc k[0123456789]*.nc | wc -l`

#if ($count_all != $YR) then
 foreach dec ( $decade )
  foreach yr ( $year )
   if ( -e $dir/$rid$a.pm$dec$yr\jan.nc || -e $dir/$rid$a.pm$dec$yr\dec.nc ) then
    set count=`ls $dir/$rid$a.pm$dec$yr*.nc | wc -l`
   else
    set count=0
   endif
   if ($count != 12) then
    if ( -e $dir/$rid$a.pm$dec$yr\jan.nc ) then 
     if ($count == 1) then
      cp $dir/$rid$a.pm$dec$yr\jan.nc $dec$yr.nc
     else if ($count == 2) then
      cdo copy $dir/$rid$a.pm$dec$yr\jan.nc $dir/$rid$a.pm$dec$yr\feb.nc $dec$yr.nc
     else if ($count == 3) then
      cdo copy $dir/$rid$a.pm$dec$yr\jan.nc $dir/$rid$a.pm$dec$yr\feb.nc $dir/$rid$a.pm$dec$yr\mar.nc $dec$yr.nc
     else if ($count == 4) then
      cdo copy $dir/$rid$a.pm$dec$yr\jan.nc $dir/$rid$a.pm$dec$yr\feb.nc $dir/$rid$a.pm$dec$yr\mar.nc $dir/$rid$a.pm$dec$yr\apr.nc $dec$yr.nc
     else if ($count == 5) then
      cdo copy $dir/$rid$a.pm$dec$yr\jan.nc $dir/$rid$a.pm$dec$yr\feb.nc $dir/$rid$a.pm$dec$yr\mar.nc $dir/$rid$a.pm$dec$yr\apr.nc $dir/$rid$a.pm$dec$yr\may.nc $dec$yr.nc
     else if ($count == 6) then
      cdo copy $dir/$rid$a.pm$dec$yr\jan.nc $dir/$rid$a.pm$dec$yr\feb.nc $dir/$rid$a.pm$dec$yr\mar.nc $dir/$rid$a.pm$dec$yr\apr.nc $dir/$rid$a.pm$dec$yr\may.nc $dir/$rid$a.pm$dec$yr\jun.nc $dec$yr.nc
     else if ($count == 7) then
      cdo copy $dir/$rid$a.pm$dec$yr\jan.nc $dir/$rid$a.pm$dec$yr\feb.nc $dir/$rid$a.pm$dec$yr\mar.nc $dir/$rid$a.pm$dec$yr\apr.nc $dir/$rid$a.pm$dec$yr\may.nc $dir/$rid$a.pm$dec$yr\jun.nc $dir/$rid$a.pm$dec$yr\jul.nc $dec$yr.nc
     else if ($count == 8) then
      cdo copy $dir/$rid$a.pm$dec$yr\jan.nc $dir/$rid$a.pm$dec$yr\feb.nc $dir/$rid$a.pm$dec$yr\mar.nc $dir/$rid$a.pm$dec$yr\apr.nc $dir/$rid$a.pm$dec$yr\may.nc $dir/$rid$a.pm$dec$yr\jun.nc $dir/$rid$a.pm$dec$yr\jul.nc $dir/$rid$a.pm$dec$yr\aug.nc $dec$yr.nc
     else if ($count == 9) then
      cdo copy $dir/$rid$a.pm$dec$yr\jan.nc $dir/$rid$a.pm$dec$yr\feb.nc $dir/$rid$a.pm$dec$yr\mar.nc $dir/$rid$a.pm$dec$yr\apr.nc $dir/$rid$a.pm$dec$yr\may.nc $dir/$rid$a.pm$dec$yr\jun.nc $dir/$rid$a.pm$dec$yr\jul.nc $dir/$rid$a.pm$dec$yr\aug.nc $dir/$rid$a.pm$dec$yr\sep.nc $dec$yr.nc
     else if ($count == 10) then
      cdo copy $dir/$rid$a.pm$dec$yr\jan.nc $dir/$rid$a.pm$dec$yr\feb.nc $dir/$rid$a.pm$dec$yr\mar.nc $dir/$rid$a.pm$dec$yr\apr.nc $dir/$rid$a.pm$dec$yr\may.nc $dir/$rid$a.pm$dec$yr\jun.nc $dir/$rid$a.pm$dec$yr\jul.nc $dir/$rid$a.pm$dec$yr\aug.nc $dir/$rid$a.pm$dec$yr\sep.nc $dir/$rid$a.pm$dec$yr\oct.nc $dec$yr.nc
     else if ($count == 11) then    
      cdo copy $dir/$rid$a.pm$dec$yr\jan.nc $dir/$rid$a.pm$dec$yr\feb.nc $dir/$rid$a.pm$dec$yr\mar.nc $dir/$rid$a.pm$dec$yr\apr.nc $dir/$rid$a.pm$dec$yr\may.nc $dir/$rid$a.pm$dec$yr\jun.nc $dir/$rid$a.pm$dec$yr\jul.nc $dir/$rid$a.pm$dec$yr\aug.nc $dir/$rid$a.pm$dec$yr\sep.nc $dir/$rid$a.pm$dec$yr\oct.nc $dir/$rid$a.pm$dec$yr\nov.nc $dec$yr.nc
     endif
    else if( -e $dir/$rid$a.pm$dec$yr\dec.nc ) then
     if ($count == 1) then
      cp $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 2) then
      cdo copy $dir/$rid$a.pm$dec$yr\nov.nc $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 3) then
      cdo copy $dir/$rid$a.pm$dec$yr\oct.nc $dir/$rid$a.pm$dec$yr\nov.nc $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 4) then
      cdo copy $dir/$rid$a.pm$dec$yr\sep.nc $dir/$rid$a.pm$dec$yr\oct.nc $dir/$rid$a.pm$dec$yr\nov.nc $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 5) then
      cdo copy $dir/$rid$a.pm$dec$yr\aug.nc $dir/$rid$a.pm$dec$yr\sep.nc $dir/$rid$a.pm$dec$yr\oct.nc $dir/$rid$a.pm$dec$yr\nov.nc $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 6) then
      cdo copy $dir/$rid$a.pm$dec$yr\jul.nc $dir/$rid$a.pm$dec$yr\aug.nc $dir/$rid$a.pm$dec$yr\sep.nc $dir/$rid$a.pm$dec$yr\oct.nc $dir/$rid$a.pm$dec$yr\nov.nc $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 7) then
      cdo copy $dir/$rid$a.pm$dec$yr\jun.nc $dir/$rid$a.pm$dec$yr\jul.nc $dir/$rid$a.pm$dec$yr\aug.nc $dir/$rid$a.pm$dec$yr\sep.nc $dir/$rid$a.pm$dec$yr\oct.nc $dir/$rid$a.pm$dec$yr\nov.nc $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 8) then
      cdo copy $dir/$rid$a.pm$dec$yr\may.nc $dir/$rid$a.pm$dec$yr\jun.nc $dir/$rid$a.pm$dec$yr\jul.nc $dir/$rid$a.pm$dec$yr\aug.nc $dir/$rid$a.pm$dec$yr\sep.nc $dir/$rid$a.pm$dec$yr\oct.nc $dir/$rid$a.pm$dec$yr\nov.nc $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 9) then
      cdo copy $dir/$rid$a.pm$dec$yr\apr.nc $dir/$rid$a.pm$dec$yr\may.nc $dir/$rid$a.pm$dec$yr\jun.nc $dir/$rid$a.pm$dec$yr\jul.nc $dir/$rid$a.pm$dec$yr\aug.nc $dir/$rid$a.pm$dec$yr\sep.nc $dir/$rid$a.pm$dec$yr\oct.nc $dir/$rid$a.pm$dec$yr\nov.nc $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 10) then
      cdo copy $dir/$rid$a.pm$dec$yr\mar.nc $dir/$rid$a.pm$dec$yr\apr.nc $dir/$rid$a.pm$dec$yr\may.nc $dir/$rid$a.pm$dec$yr\jun.nc $dir/$rid$a.pm$dec$yr\jul.nc $dir/$rid$a.pm$dec$yr\aug.nc $dir/$rid$a.pm$dec$yr\sep.nc $dir/$rid$a.pm$dec$yr\oct.nc $dir/$rid$a.pm$dec$yr\nov.nc $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 11) then
      cdo copy $dir/$rid$a.pm$dec$yr\feb.nc $dir/$rid$a.pm$dec$yr\mar.nc $dir/$rid$a.pm$dec$yr\apr.nc $dir/$rid$a.pm$dec$yr\may.nc $dir/$rid$a.pm$dec$yr\jun.nc $dir/$rid$a.pm$dec$yr\jul.nc $dir/$rid$a.pm$dec$yr\aug.nc $dir/$rid$a.pm$dec$yr\sep.nc $dir/$rid$a.pm$dec$yr\oct.nc $dir/$rid$a.pm$dec$yr\nov.nc $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
     endif
    endif
   else 
    if ( -e $dir/$rid$a.pm$dec$yr\jan.nc && -e $dir/$rid$a.pm$dec$yr\dec.nc ) then
     cdo -b 64 copy $dir/$rid$a.pm$dec$yr\jan.nc $dir/$rid$a.pm$dec$yr\feb.nc $dir/$rid$a.pm$dec$yr\mar.nc $dir/$rid$a.pm$dec$yr\apr.nc $dir/$rid$a.pm$dec$yr\may.nc $dir/$rid$a.pm$dec$yr\jun.nc $dir/$rid$a.pm$dec$yr\jul.nc $dir/$rid$a.pm$dec$yr\aug.nc $dir/$rid$a.pm$dec$yr\sep.nc $dir/$rid$a.pm$dec$yr\oct.nc $dir/$rid$a.pm$dec$yr\nov.nc $dir/$rid$a.pm$dec$yr\dec.nc $dec$yr.nc
    endif
   endif
  end
 end
#endif

@ numyr = $YR + 1
set count_all=`ls h[0123456789]*.nc i[0123456789]*.nc j[0123456789]*.nc k[0123456789]*.nc | wc -l`

if ($SPLIT == n) then
 if ($count_all == $numyr) then

  set short = `ls h[0123456789]*.nc i[0123456789]*.nc j[0123456789]*.nc k[0123456789]*.nc | head -${numyr}`
  cdo copy $short Mmonthly_means_${YR}yrs.nc
  cdo yearavg Mmonthly_means_${YR}yrs.nc yearly_means_${YR}yrs.nc

 else

  set short = `ls h[0123456789]*.nc i[0123456789]*.nc j[0123456789]*.nc k[0123456789]*.nc | head -${YR}`  
  cdo copy $short Mmonthly_means_${YR}yrs.nc  
  cdo yearavg Mmonthly_means_${YR}yrs.nc yearly_means_${YR}yrs.nc

 endif

rm h[0123456789]*.nc i[0123456789]*.nc j[0123456789]*.nc k[0123456789]*.nc

endif

exit
