#!/bin/csh 

unalias ls

set a=a
set dir=$DIR
set rid=$RUNID
set ext=$Pmonth
set decade=(h i j k)
set year=(0 1 2 3 4 5 6 7 8 9)

rm h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc

dmget $DIR/$rid$a.$ext?????.nc

foreach dec ( $decade )
  foreach yr ( $year )

    if ( -e $dir/$rid$a.$ext$dec$yr\001.nc || -e $dir/$rid$a.$ext$dec$yr\012.nc ) then
     set count=`ls $dir/$rid$a.$ext$dec$yr???.nc | wc -l`
    else
     set count=0
    endif

    if ($count != 12 && $count != 0) then
     setenv fullstrt F
     cdo mergetime $dir/$rid$a.$ext$dec$yr???.nc $dec$yr.nc
    else 
     #echo $count
     setenv fullstrt T
     if ( -e $dir/$rid$a.$ext$dec$yr\001.nc && -e $dir/$rid$a.$ext$dec$yr\012.nc ) then
     # setenv fullstrt T
      cdo -b 64 mergetime $dir/$rid$a.$ext$dec$yr???.nc $dec$yr.nc
     #else
     # echo "(1)     Error in cdo_merge for Mmonthly file"
     endif

    endif

  end
end

@ numyr = $YR + 1
set count_all=`ls h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc | wc -l`

if ($SPLIT == n) then
 if ($count_all == $numyr && $fullstrt == F) then

  set short = `ls h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc | head -${numyr}`
  cdo mergetime $short Mmonthly_means_${YR}yrs.nc
  cdo yearmean Mmonthly_means_${YR}yrs.nc yearly_means_${YR}yrs.nc

 else

  set short = `ls h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc | head -${YR}`  
  cdo mergetime $short Mmonthly_means_${YR}yrs.nc  
  #cdo mergetime h?.nc Mmon_dec_70s.nc  
  #cdo mergetime i?.nc Mmon_dec_80s.nc  
  #cdo mergetime j?.nc Mmon_dec_90s.nc  
  #cdo mergetime k?.nc Mmon_dec_00s.nc  
  cdo yearmean Mmonthly_means_${YR}yrs.nc yearly_means_${YR}yrs.nc

 endif

rm h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc

endif

exit

