#!/bin/csh

if (! -d plots) then
echo "(0)     Mkdir plots/" 
mkdir plots
else
echo "(0)     Directory plots/ Already Exists"
endif

if (${CHFMT} == y ) then
echo "(0)     Run is producing Jpeg files"
if (! -d plots_jpg) then
echo "(0)     Mkdir plots_jpg/" 
mkdir plots_jpg
else
echo "(0)     Directory plots_jpg/ Already Exists"
endif
endif

exit
