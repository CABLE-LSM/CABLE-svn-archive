#!/bin/csh 

unalias ls # or add /bin ?

set a=a
set dir=$DIR
set rid=$RUNID
set ext=$Pmonth
set flist=`ls $dir/$rid.$ext-??????????.nc`
set fone = $flist[1]
@ year=`echo $fone | tail -c14 | head -c4`  # year 0001
set fnum=`ls $dir/$rid.$ext-??????????.nc | wc -l`
@ nof = ${fnum} / 12

rm y????.nc

dmget $DIR/$rid.$ext-??????????.nc

while ( $year <= $nof )

    if ( -e $dir/$rid.$ext-*0$year\001???.nc || -e $dir/$rid.$ext-*0$year\012???.nc ) then
     set count=`ls $dir/$rid.$ext-*0$year??????.nc | wc -l`
    else
     set count=0
    endif

    if ($count != 12 && $count != 0) then
     setenv fullstrt F
     cdo mergetime $dir/$rid.$ext-*0$year??????.nc tmp.nc
    else 
     setenv fullstrt T
     if ( -e $dir/$rid.$ext-*0$year\001???.nc && -e $dir/$rid.$ext-*0$year\012???.nc ) then
      cdo -b 64 mergetime $dir/$rid.$ext-*0$year??????.nc tmp.nc
     endif
    endif

    if ( $year < 10 ) then
     mv tmp.nc y000$year.nc
    else
     if ( $year < 100 ) then
      mv tmp.nc y00$year.nc
     else
      if ( $year < 1000 ) then
       mv tmp.nc y0$year.nc
      else
       mv tmp.nc y$year.nc
      endif
     endif
    endif

 @ year++
end

@ numyr = $YR + 1
set count_all=`ls y????.nc | wc -l`

if ($SPLIT == n) then
 if ($count_all == $numyr && $fullstrt == F) then

  set short = `ls y????.nc | head -${numyr}`
  cdo mergetime $short Mmonthly_means_${YR}yrs.nc
  cdo yearmean Mmonthly_means_${YR}yrs.nc yearly_means_${YR}yrs.nc

 else

  set short = `ls y????.nc | head -${YR}`  
  cdo mergetime $short Mmonthly_means_${YR}yrs.nc  
  cdo yearmean Mmonthly_means_${YR}yrs.nc yearly_means_${YR}yrs.nc

 endif

rm y????.nc

endif

exit

