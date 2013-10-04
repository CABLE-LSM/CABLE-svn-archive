#!/bin/csh

#move to the "new" directory and get the name of the dump, given the new code tested 
cd new
set newdump = `ls ../../$1/$1a.da*`
echo $newdump

#open this new dump, convert the desired field to name = "t"
xconv $newdump 
mv t new.nc

#run ferrret script to compute diff and record as ps file
ferret -script diff.jnl


