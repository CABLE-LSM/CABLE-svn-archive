#!/bin/ksh

# RUN_CABLE qsub
qsub_avail()
{
   qstat -u $uid | grep $uid > .qu
   nqu=`wc -l .qu | cut -c $nch`
   
   if [[ $nqu < $nproc ]]; then 
     
      mkdir $qs$i
      mkdir $qs$i/src
      cp RUN_CABLE $qs$i
      cp src/main.txt src/sourced.txt $qs$i/src
      cp src/*R $qs$i/src/
      cp src/*ksh $qs$i/src/  
      mv src/$qs_filename_base$i$qs_filename_suffix $qs$i
      
      cd $qs$i
      ln -s ../data
      mv $qs_filename_base$i$qs_filename_suffix sites_main.txt
      
      /opt/pbs/bin/qsub -S /bin/ksh src/QRUN.ksh
      
      cd ../ 
   else
      sleep 1m
      qsub_avail 
   fi
}


quit_option()   
{
   print "\nIf you wish to quit now, just hit enter\c" 
}  

quit_response()
{      
      print "\nAdios"
      exit
}


help_option()
{
   print "\nIf you wish to see more information on command line arguments to RUN_CABLE , enter "'"H"'" " 
} 




