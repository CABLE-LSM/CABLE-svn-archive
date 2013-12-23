#!/bin/csh -x

set file=$1
set rotate=$2

ps2epsi file.ps file.eps

if ( $rotate == "n" ) then
 convert -density 500 file.eps file.jpg
else
 convert -density 500 -rotate 270 file.eps file.jpg
endif

exit

