#!/bin/csh 

unalias ls

set a=a
set dir=$DIR
set rid=$RUNID
set ext=$Pmonth
if ($VN >= 85) then
 set decade=(000 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030)
 #set decade=(185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210)
else
 if ($rid == vaoyl) then
  set decade=(1 2 3 4 5 6 7 8 9)
 else
  set decade=(h i j k)
 endif
endif
set year=(0 1 2 3 4 5 6 7 8 9)

rm [0123456789][0123456789].nc 
rm h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc
rm [012][0123456789][0123456789][0123456789].nc

#dmget $DIR/$rid$a.$ext?????.nc

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

  end # yr
end   # dec

@ numyr = $YR + 1
if ($VN >= 85) then
 set count_all=`ls [012][0123456789][0123456789][0123456789].nc | wc -l`
else
 set count_all=`ls h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc | wc -l`
endif

if ($SPLIT == n) then
 if ($count_all == $numyr && $fullstrt == F) then

  if ($VN >= 85) then
    set short = `ls [012][0123456789][0123456789][0123456789].nc | head -${numyr}`
  else
  if ($rid == vaoyl) then
    set short = `ls [0123456789][0123456789].nc | head -${YR}`
  else
    set short = `ls h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc | head -${numyr}`
  endif
  endif
  cdo -b 64 mergetime $short Mmonthly_means_${YR}yrs.nc
  cdo yearmonmean Mmonthly_means_${YR}yrs.nc yearly_means_${YR}yrs.nc

 else

  if ($VN >= 85) then
    set short = `ls [012][0123456789][0123456789][0123456789].nc | head -${YR}`  
  else
  if ($rid == vaoyl) then
    set short = `ls [0123456789][0123456789].nc | head -${YR}`  
  else
  set short = `ls h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc | head -${YR}`  
  endif
  endif
  cdo mergetime $short Mmonthly_means_${YR}yrs.nc  
  cdo yearmonmean Mmonthly_means_${YR}yrs.nc yearly_means_${YR}yrs.nc

 endif

rm [0123456789][0123456789].nc
rm h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc
rm [012][0123456789][0123456789][0123456789].nc

endif

exit

