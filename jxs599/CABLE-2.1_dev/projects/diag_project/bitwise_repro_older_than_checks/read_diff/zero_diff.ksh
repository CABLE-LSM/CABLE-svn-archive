#!/bin/ksh

basename="std"	
basename2="FLUXES00"	

#### write a little file for consumption by f90 code, interpreting command line args
print $basename > input.dat
print $basename2 >> input.dat

### execute f90 program to do all the work desired
./diff_main





































