#!/bin/csh
rm -f cat_Nnodes
make
make clean
mv cat_Nnodes ../../bin
echo ""
echo ""
echo "assuming successful binary has been moved to ../../bin"
echo ""
echo "Please hit enter to continue"
echo ""
set test = "$<"
tput reset
