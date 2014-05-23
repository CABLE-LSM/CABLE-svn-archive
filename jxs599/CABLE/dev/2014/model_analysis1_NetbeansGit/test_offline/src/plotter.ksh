#!/bin/ksh
   cd src/
   print '\n*** PLOTTING CABLE FLUX DATA ETC ***\n'
   print '\nThis may take some time.\nIf desirable turn off unneccessary plots in plot_main.nml or atleast choose PNG files\n'
   R CMD BATCH --slave plot.R
   print '\n*** FINISHED PLOTTING CABLE DATA ***\n'



