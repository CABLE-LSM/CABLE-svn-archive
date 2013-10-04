#!/bin/ksh 

RUN_CABLE clean
rm -f Sum*
rm -f ../build/cable 
##on home2 n(ETWORK)fs does not like cp --preserve
#rm -fr ../build/tmp
RUN_CABLE build
RUN_CABLE run 
cable_diff.pl


