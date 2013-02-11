#!/bin/csh 

echo "This script htmlizes your fortran code for easy browsing. "
echo " "
echo "All you need to do is give me the absolute path of your CABLE checkout."

set src_path = $<
set re_path = $src_path"/"

set htmlin = "CABLE_html"

./cable_prepare.pl $re_path

echo '$dir_html="'$htmlin'/";' >  test

cat cable.f2h_start test cable.f2h_end > cable.f2h
rm -f test 
./f90tohtml cable.f2h 

