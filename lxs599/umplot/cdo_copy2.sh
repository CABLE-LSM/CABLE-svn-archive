#!/bin/csh

unalias ls

set a=a
set dir=$DIR
set rid=$RUNID
set ext=$Pmonth
set decade=(h i j k)
set year=(0 1 2 3 4 5 6 7 8 9)

rm h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc

#dmget $DIR/$rid$a.$ext?????.nc

#set count_all=`ls h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc | wc -l`

#if ($count_all != $YR) then
 foreach dec ( $decade )
  foreach yr ( $year )
   if ( -e $dir/$rid$a.$ext$dec$yr\jan.nc || -e $dir/$rid$a.$ext$dec$yr\dec.nc ) then
    set count=`ls $dir/$rid$a.$ext$dec$yr???.nc | wc -l`
   else
    set count=0
   endif
   if ($count != 12) then
    setenv fullstrt F
    if ( -e $dir/$rid$a.$ext$dec$yr\jan.nc ) then 
     if ($count == 1) then
      cp $dir/$rid$a.$ext$dec$yr\jan.nc $dec$yr.nc
     else if ($count == 2) then
      cdo copy $dir/$rid$a.$ext$dec$yr\jan.nc $dir/$rid$a.$ext$dec$yr\feb.nc $dec$yr.nc
     else if ($count == 3) then
      cdo copy $dir/$rid$a.$ext$dec$yr\jan.nc $dir/$rid$a.$ext$dec$yr\feb.nc $dir/$rid$a.$ext$dec$yr\mar.nc $dec$yr.nc
     else if ($count == 4) then
      cdo copy $dir/$rid$a.$ext$dec$yr\jan.nc $dir/$rid$a.$ext$dec$yr\feb.nc $dir/$rid$a.$ext$dec$yr\mar.nc $dir/$rid$a.$ext$dec$yr\apr.nc $dec$yr.nc
     else if ($count == 5) then
      cdo copy $dir/$rid$a.$ext$dec$yr\jan.nc $dir/$rid$a.$ext$dec$yr\feb.nc $dir/$rid$a.$ext$dec$yr\mar.nc $dir/$rid$a.$ext$dec$yr\apr.nc $dir/$rid$a.$ext$dec$yr\may.nc $dec$yr.nc
     else if ($count == 6) then
      cdo copy $dir/$rid$a.$ext$dec$yr\jan.nc $dir/$rid$a.$ext$dec$yr\feb.nc $dir/$rid$a.$ext$dec$yr\mar.nc $dir/$rid$a.$ext$dec$yr\apr.nc $dir/$rid$a.$ext$dec$yr\may.nc $dir/$rid$a.$ext$dec$yr\jun.nc $dec$yr.nc
     else if ($count == 7) then
      cdo copy $dir/$rid$a.$ext$dec$yr\jan.nc $dir/$rid$a.$ext$dec$yr\feb.nc $dir/$rid$a.$ext$dec$yr\mar.nc $dir/$rid$a.$ext$dec$yr\apr.nc $dir/$rid$a.$ext$dec$yr\may.nc $dir/$rid$a.$ext$dec$yr\jun.nc $dir/$rid$a.$ext$dec$yr\jul.nc $dec$yr.nc
     else if ($count == 8) then
      cdo copy $dir/$rid$a.$ext$dec$yr\jan.nc $dir/$rid$a.$ext$dec$yr\feb.nc $dir/$rid$a.$ext$dec$yr\mar.nc $dir/$rid$a.$ext$dec$yr\apr.nc $dir/$rid$a.$ext$dec$yr\may.nc $dir/$rid$a.$ext$dec$yr\jun.nc $dir/$rid$a.$ext$dec$yr\jul.nc $dir/$rid$a.$ext$dec$yr\aug.nc $dec$yr.nc
     else if ($count == 9) then
      cdo copy $dir/$rid$a.$ext$dec$yr\jan.nc $dir/$rid$a.$ext$dec$yr\feb.nc $dir/$rid$a.$ext$dec$yr\mar.nc $dir/$rid$a.$ext$dec$yr\apr.nc $dir/$rid$a.$ext$dec$yr\may.nc $dir/$rid$a.$ext$dec$yr\jun.nc $dir/$rid$a.$ext$dec$yr\jul.nc $dir/$rid$a.$ext$dec$yr\aug.nc $dir/$rid$a.$ext$dec$yr\sep.nc $dec$yr.nc
     else if ($count == 10) then
      cdo copy $dir/$rid$a.$ext$dec$yr\jan.nc $dir/$rid$a.$ext$dec$yr\feb.nc $dir/$rid$a.$ext$dec$yr\mar.nc $dir/$rid$a.$ext$dec$yr\apr.nc $dir/$rid$a.$ext$dec$yr\may.nc $dir/$rid$a.$ext$dec$yr\jun.nc $dir/$rid$a.$ext$dec$yr\jul.nc $dir/$rid$a.$ext$dec$yr\aug.nc $dir/$rid$a.$ext$dec$yr\sep.nc $dir/$rid$a.$ext$dec$yr\oct.nc $dec$yr.nc
     else if ($count == 11) then    
      cdo copy $dir/$rid$a.$ext$dec$yr\jan.nc $dir/$rid$a.$ext$dec$yr\feb.nc $dir/$rid$a.$ext$dec$yr\mar.nc $dir/$rid$a.$ext$dec$yr\apr.nc $dir/$rid$a.$ext$dec$yr\may.nc $dir/$rid$a.$ext$dec$yr\jun.nc $dir/$rid$a.$ext$dec$yr\jul.nc $dir/$rid$a.$ext$dec$yr\aug.nc $dir/$rid$a.$ext$dec$yr\sep.nc $dir/$rid$a.$ext$dec$yr\oct.nc $dir/$rid$a.$ext$dec$yr\nov.nc $dec$yr.nc
     endif
    else if( -e $dir/$rid$a.$ext$dec$yr\dec.nc ) then
     if ($count == 1) then
      cp $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 2) then
      cdo copy $dir/$rid$a.$ext$dec$yr\nov.nc $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 3) then
      cdo copy $dir/$rid$a.$ext$dec$yr\oct.nc $dir/$rid$a.$ext$dec$yr\nov.nc $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 4) then
      cdo copy $dir/$rid$a.$ext$dec$yr\sep.nc $dir/$rid$a.$ext$dec$yr\oct.nc $dir/$rid$a.$ext$dec$yr\nov.nc $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 5) then
      cdo copy $dir/$rid$a.$ext$dec$yr\aug.nc $dir/$rid$a.$ext$dec$yr\sep.nc $dir/$rid$a.$ext$dec$yr\oct.nc $dir/$rid$a.$ext$dec$yr\nov.nc $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 6) then
      cdo copy $dir/$rid$a.$ext$dec$yr\jul.nc $dir/$rid$a.$ext$dec$yr\aug.nc $dir/$rid$a.$ext$dec$yr\sep.nc $dir/$rid$a.$ext$dec$yr\oct.nc $dir/$rid$a.$ext$dec$yr\nov.nc $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 7) then
      cdo copy $dir/$rid$a.$ext$dec$yr\jun.nc $dir/$rid$a.$ext$dec$yr\jul.nc $dir/$rid$a.$ext$dec$yr\aug.nc $dir/$rid$a.$ext$dec$yr\sep.nc $dir/$rid$a.$ext$dec$yr\oct.nc $dir/$rid$a.$ext$dec$yr\nov.nc $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 8) then
      cdo copy $dir/$rid$a.$ext$dec$yr\may.nc $dir/$rid$a.$ext$dec$yr\jun.nc $dir/$rid$a.$ext$dec$yr\jul.nc $dir/$rid$a.$ext$dec$yr\aug.nc $dir/$rid$a.$ext$dec$yr\sep.nc $dir/$rid$a.$ext$dec$yr\oct.nc $dir/$rid$a.$ext$dec$yr\nov.nc $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 9) then
      cdo copy $dir/$rid$a.$ext$dec$yr\apr.nc $dir/$rid$a.$ext$dec$yr\may.nc $dir/$rid$a.$ext$dec$yr\jun.nc $dir/$rid$a.$ext$dec$yr\jul.nc $dir/$rid$a.$ext$dec$yr\aug.nc $dir/$rid$a.$ext$dec$yr\sep.nc $dir/$rid$a.$ext$dec$yr\oct.nc $dir/$rid$a.$ext$dec$yr\nov.nc $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 10) then
      cdo copy $dir/$rid$a.$ext$dec$yr\mar.nc $dir/$rid$a.$ext$dec$yr\apr.nc $dir/$rid$a.$ext$dec$yr\may.nc $dir/$rid$a.$ext$dec$yr\jun.nc $dir/$rid$a.$ext$dec$yr\jul.nc $dir/$rid$a.$ext$dec$yr\aug.nc $dir/$rid$a.$ext$dec$yr\sep.nc $dir/$rid$a.$ext$dec$yr\oct.nc $dir/$rid$a.$ext$dec$yr\nov.nc $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
     else if ($count == 11) then
      cdo copy $dir/$rid$a.$ext$dec$yr\feb.nc $dir/$rid$a.$ext$dec$yr\mar.nc $dir/$rid$a.$ext$dec$yr\apr.nc $dir/$rid$a.$ext$dec$yr\may.nc $dir/$rid$a.$ext$dec$yr\jun.nc $dir/$rid$a.$ext$dec$yr\jul.nc $dir/$rid$a.$ext$dec$yr\aug.nc $dir/$rid$a.$ext$dec$yr\sep.nc $dir/$rid$a.$ext$dec$yr\oct.nc $dir/$rid$a.$ext$dec$yr\nov.nc $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
     endif
    endif
   else 
    if ( -e $dir/$rid$a.$ext$dec$yr\jan.nc && -e $dir/$rid$a.$ext$dec$yr\dec.nc ) then
     setenv fullstrt T
     cdo -b 64 copy $dir/$rid$a.$ext$dec$yr\jan.nc $dir/$rid$a.$ext$dec$yr\feb.nc $dir/$rid$a.$ext$dec$yr\mar.nc $dir/$rid$a.$ext$dec$yr\apr.nc $dir/$rid$a.$ext$dec$yr\may.nc $dir/$rid$a.$ext$dec$yr\jun.nc $dir/$rid$a.$ext$dec$yr\jul.nc $dir/$rid$a.$ext$dec$yr\aug.nc $dir/$rid$a.$ext$dec$yr\sep.nc $dir/$rid$a.$ext$dec$yr\oct.nc $dir/$rid$a.$ext$dec$yr\nov.nc $dir/$rid$a.$ext$dec$yr\dec.nc $dec$yr.nc
    endif
   endif
  end
 end
#endif

@ numyr = $YR + 1
set count_all=`ls h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc | wc -l`

if ($SPLIT == n) then
 if ($count_all == $numyr && $fullstrt == F) then

  set short = `ls h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc | head -${numyr}`
  cdo copy $short Mmonthly_means_${YR}yrs.nc
  cdo yearmean Mmonthly_means_${YR}yrs.nc yearly_means_${YR}yrs.nc

 else

  set short = `ls h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc | head -${YR}`  
  cdo copy $short Mmonthly_means_${YR}yrs.nc  
  cdo yearmean Mmonthly_means_${YR}yrs.nc yearly_means_${YR}yrs.nc

 endif

rm h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc

endif

exit
