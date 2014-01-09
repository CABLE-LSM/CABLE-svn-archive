#!/bin/csh

cd new
if -f wabqla.daj0120 then
   rm -f wabqla.daj0120
endif   
ls -al ../../wabql/wabqla.daj0120
ln -s ../../wabql/wabqla.daj0120
xconv wabqla.daj0120
mv t new.nc
ferret -script diff.jnl

