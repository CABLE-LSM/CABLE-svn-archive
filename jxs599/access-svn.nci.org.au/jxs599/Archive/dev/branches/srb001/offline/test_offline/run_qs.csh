#!/bin/csh 

set qs = 'qsub_'
set qs_filename_base = 'qsmain_'
set qs_filename_suffix = '.nml'

cd ../
   if ( -e core ) then
      rm -f core
   endif   
   ln -s ../core
cd test_offline/

R CMD BATCH --slave Qsub.R

foreach i ( `cat qsj.j` )

   mkdir $qs$i
   cp * $qs$i
   cp -r src/ $qs$i
   cd $qs$i
   ln -s ../data
   mv $qs_filename_base$i$qs_filename_suffix main.nml
   cat plot_main.nml >> main.nml
   /opt/pbs/bin/qsub -S /bin/csh run.csh 
   cd ../

end      






