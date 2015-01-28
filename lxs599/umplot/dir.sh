#!/bin/csh

if (! -d plots) then
 echo "(0)     Mkdir plots/" 
 mkdir plots
 #mkdir plots_${YR}yrs
 #echo ""
else
 echo "(0)     Directory plots/ Already Exists"
 echo "(0)"
endif

if (${CHFMT} == y ) then
 echo "(0)     Run is producing Jpeg files"
 if (! -d plots_jpg) then
  echo "(0)     Mkdir plots_jpg/" 
  mkdir plots_jpg
  #echo ""
 else
  echo "(0)     Directory plots_jpg/ Already Exists"
  #echo ""
 endif
endif

exit
