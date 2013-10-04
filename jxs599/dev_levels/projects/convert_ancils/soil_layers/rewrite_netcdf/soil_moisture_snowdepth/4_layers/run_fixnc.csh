#!/bin/csh
#############################################################################
###  after a bit of book-keeping this csh script will compile CABLE       ###
###  from the "../src" directory and deposit the executable binary in     ###
###  this directory. CABLE is then run over each of the sites specified   ### 
###  in "main.nml" and output data moved into a created directory (out)   ###
###  labelled by the name of the site. finally, an "R" script to produce  ###
###  plots of flux data is called.                                        ### 
#############################################################################

######################################################################
### call plot.R in batch mode to avoid going into R first          ###
######################################################################
R CMD BATCH --slave fixnc.R


