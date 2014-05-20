#!/bin/csh
# ==========================
# Lauren Stevens 15 Oct 2013
# ==========================

set a    = a
set rid  = $RUNID
set dir  = $DIR
set xlist= ( pc pb pf pa pe pj pi )
#set xlist=( pc $Pmonth $Ptemp1 $Ptemps $Ptimes $PcasaC $Pmonth $Pdaily )

set yr2  = ( 01 23 45 67 89 )
set yr5  = ( 01234 56789 )
set dec1 = ( h i j k )
set i = 1

cd $dir

foreach ext ( $xlist )

  dmget $rid$a.$ext*

  if ( $ext == $Ptimes ) then
   set plist=`ls $rid$a.$ext?????_*swlw.nc`
   if ( ${#plist} != 0 ) then
    if ( ! -e $rid\_$ext\nc_swlw.tar ) then
     tar cvf - $rid$a.$ext?????_*swlw.nc > $rid\_$ext\nc_swlw.tar
    endif
   endif
  endif

foreach dec ( $dec1 )

 if ( $ext != $Pdaily && $ext != pi ) then

  foreach yr ( $yr5 )
   set flist=`ls $rid$a.$ext${dec}[${yr}]???`
   set nlist=`ls $rid$a.$ext${dec}[${yr}]???.nc`
   if ( ${#flist} != 0 ) then
    if ( ! -e $rid\_$ext\_${dec}${i}.tar ) then
     tar cvf - $rid$a.$ext${dec}[${yr}]??? > $rid\_$ext\_${dec}${i}.tar
    endif
   endif
   if ( ${#nlist} != 0 ) then
    if ( ! -e $rid\_$ext\_${dec}$i\nc.tar ) then
     tar cvf - $rid$a.$ext${dec}[${yr}]???.nc > $rid\_$ext\_${dec}$i\nc.tar
    endif
   endif
   @ i++
  end
  @ i = 1

 else # if == pdaily,pi

  #echo "Daily Files"
  foreach yr ( $yr2 )
   set flist=`ls $rid$a.$ext${dec}[${yr}]???`
   set nlist=`ls $rid$a.$ext${dec}[${yr}]???.nc`
   if ( ${#flist} != 0 ) then
    if ( ! -e $rid\_$ext\_${dec}${yr}.tar ) then
     tar cvf - $rid$a.$ext${dec}[${yr}]??? > $rid\_$ext\_${dec}${yr}.tar
    endif
   endif
   if ( ${#nlist} != 0 ) then
    if ( ! -e $rid\_$ext\_${dec}${yr}\nc.tar ) then
     tar cvf - $rid$a.$ext${dec}[${yr}]???.nc > $rid\_$ext\_${dec}${yr}\nc.tar
    endif
   endif
  end

 endif # if != pdaily
end # foreach dec
end # foreach ext

#if ($SPLIT == y) then
# foreach blck ( $block )
#  cd $dir/block$blck\_5yrs
#  tar cvf - *5yr*ps > plots_blk$blck\_5yrs.tar
#  if ($CHFMT == y) then
#   tar cvf - *5yr*jpg > plots_jpg_blk$blck\_5yrs.tar
#  endif
#  tar cvf - *5yr*nc > nc_blk$blck\_5yrs.tar
# end
# cd $dir
#else
# cd $dir/plots
# tar cvf - *${YR}yr*ps > plots_${YR}yrs.tar
# if ($CHFMT == y) then
#  cd $dir/plots_jpg
#  tar cvf - *${YR}yr*jpg > plots_jpg_${YR}yrs.tar
# endif
# cd $dir
# tar cvf - *${YR}yr*nc > umplot_nc_${YR}yrs.tar
#endif

#ls [AmsTy]*${YR}yr*.nc Me*${YR}yr*.nc
#ls [ATy]*${YR}yr*.nc Me*${YR}yr*.nc mm*${YR}yr*.nc
#ls [Ay]*${YR}yr*.nc Me*${YR}yr*.nc mm*${YR}yr*.nc T[ei]*${YR}yr*.nc
#ls A*${YR}yr*.nc Me*${YR}yr*.nc mm*${YR}yr*.nc T[ei]*${YR}yr*.nc

exit

