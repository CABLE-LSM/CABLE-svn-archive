#!/bin/ksh

build_cable()
{      
   if [[ -f cable ]]; then
      rm -f cable
   fi  
   
   cd $1 
      
      if [[ -f cable ]]; then
         rm -f cable
      fi  

      ./build.ksh 
   
      if [[ -f cable ]]; then
         print '\n*** CABLE BUILD SUCCESSFULL ***\n\nExecutable  will be copied to directory:\n'
         print '\t'$mypwd"\n" 
         /bin/cp cable $mypwd 
      else
         print '\n*** ERROR: BUILD FAIELED ***\n'
         exit
      fi 
   
   cd $mypwd      
}
 


run_cable() 
{
   #remove any trace of previous runs
   rm -f fort.66 *nc 
   rm -f *.txt *00.bin *00.dat
   cd $1/src
   
   print '\n*** RUNNING CABLE ***\n'
   
   R CMD BATCH --slave CABLE.R
   print '\n*** CABLE RUN FINISHED ***\n'
   
   cp CABLE.Rout ../out

   cd ../
   
   if [[ -f out_cable.nc ]]; then
      print '\n*** CABLE RUN (appears) SUCCESSFULL ***\n'
   else
      print '\n*** ERROR: RUN FAILED ***\n'     
   fi  
   
   rm -f out_cable.nc restart_out.nc log_cable.txt fort.66
}



plot_cable()
{
   cd src/
   print '\n*** PLOTTING CABLE FLUX DATA ETC ***\n'
   print '\nThis may take some time.\nIf desirable turn off unneccessary plots in plot_main.nml or atleast choose PNG files\n'
   R CMD BATCH --slave plot.R
   print '\n*** FINISHED PLOTTING CABLE DATA ***\n'
   cp plot.Rout ../out
}



qsub_avail()
{
   qstat -u $uid | grep $uid > .qu
   nqu=`wc -l .qu | cut -c $nch`
   
   if [[ $nqu < $nproc ]]; then 
     
      mkdir $qs$i
      mkdir $qs$i/src
      cp RUN_CABLE $qs$i
      cp main.nml sourced.nml $qs$i
      cp src/*R $qs$i/src/
      cp src/*ksh $qs$i/src/  
      mv src/$qs_filename_base$i$qs_filename_suffix $qs$i
      
      cd $qs$i
      ln -s ../data
      mv $qs_filename_base$i$qs_filename_suffix sites_main.nml
      
      /opt/pbs/bin/qsub -S /bin/ksh RUN_CABLE
      
      cd ../ 
   else
      sleep 1m
      qsub_avail 
   fi
}



#--- writes welcome banner to screen outlining basic operation of RUN_CABLE
banner_welcome()
{
   clear

   banner_welcome_txt
   config_run

} 
 
banner_welc_sub()
{
   print "\n\t1. Building a CABLE executable from source code."       
   print "\n\t2. Running CABLE over each of the sites specified in"  
   print "\tsites_main.nml. Either consecutively in your current shell or"  
   print "\tconcurrently as qsub jobs on for example vayu."  
   print "\n\t3. Plotting flux data computed by CABLE for these sites."
}



banner_welcome_txt()
{   
   print "\nThis script is currenntly capable of:" 

   banner_welc_sub
   
   print "\nIn addition to being able to execute these tasks independently,"
   print "RUN_CABLE can also be configured to perform combinations of these" 
   print "tasks. RUN_CABLE can be configured via arguments at the" 
   print "command-line, or via this script.\n" 
   print "Hit enter to  display the the various incantations of RUN_CABLE.\c"
   read duumy
}


  
config_run()
{

   clear
   config_text
  
   read response
   
   if [[ $response == '1' ]]; then
      RUN_CABLE build 
   fi
  
   if [[ $response == '2' ]]; then
      RUN_CABLE build run
   fi
  
   if [[ $response == '3' ]]; then
      RUN_CABLE run 
   fi
  
   if [[ $response == '4' ]]; then
      RUN_CABLE plot
   fi
  
   if [[ $response == '5' ]]; then
      RUN_CABLE run plot 
   fi
  
   if [[ $response == '6' ]]; then
      RUN_CABLE all 
   fi

   if [[ $response == '7' ]]; then
      RUN_CABLE qsub 
   fi

   if [[ $response == '8' ]]; then
      RUN_CABLE clean
      print "\nsqueaky"
      quit_response 
   fi

   if [[ $response == '9' ]]; then
      RUN_CABLE more_help
      print "\n NA - yet!"
      quit_response 
   fi

   if [[ $response == '' ]]; then
      quit_response 
   fi

   if [[ $response == "P" ]]; then
      config_text > junk
      quit_response 
   fi
}



config_text()
{
   print "\nGenerically: RUN_CABLE [OPTION1] [OPTION2]\n "
   
   print "\nRUN_CABLE without any additional arguments will bring you"
   print "back here"
   print "\n1. RUN_CABLE build "
   print "\t\tBuild CABLE only. Useful in code development, and also to "  
   print "\t\tcreate an executable for a qsub job."
   
   print "\n2. RUN_CABLE build run"
   print "\t\tBuild CABLE and then run it."
   
   print "\n3. RUN_CABLE run "
   print "\t\tCABLE executable already exists, just run it."
   print "\t\t\tUseful if you are switching data sets, or have a pre-built"  
   print "\t\tspecial version of CABLE"
   
   print "\n4. RUN_CABLE plot "
   print "\t\tPlot existing CABLE data (data must be in out/ directory)."
   
   print "\n5. RUN_CABLE run plot "
   print "\t\tRun existing CABLE executable and plot this new data"  
   print "\t\t(will be left in a newly created out/ directory)."
   
   print "\n6. RUN_CABLE all "
   print "\t\tDo everything. Build CABLE, then run it, then plot it."
  
   print "\n7. RUN_CABLE qsub "
   print "\t\tLike "RUN_CABLE all" but on multiple processors. One site per"
   print "\t\tprocessor. Note: Only forces build if executable is not " 
   print "\t\tavailable. Submits as many sites as it can to the queue in" 
   print "\t\tseperate run directories qsub_X (X=1,2,3,....)[NOTE:most NCI" 
   print "\t\tusers have a default limit of 8 jobs at a time]"
   
   print "\n8. RUN_CABLE clean"
   print "\t\tGive me back my nice clean directory."
   
   print "\n9. RUN_CABLE help [OPTION]"
   print "\t\tGet further help on a particular option."
   
   print "\nIf you really really like paper and want to print, you can print" 
   print "this page to a file by Entering "'"P"'" "

   print "\nEnter option number to proceed from here OR Enter to quit"
}

 
help()
{     
   config_run
}
         
book_keeping()
{      
   if [[ -d out.9 ]]; then
      print "\n\ntime to organize your out/ directory"
      exit
   fi 
   
   i=8; ((j=i+1)); k=0;
   while [ $k -lt 8 ]
      do 
         if [[ -d out.$i ]]; then
            mv out.$i out.$j
         fi 
         ((j = j - 1)) 
         ((i = i - 1)) 
         ((k = k + 1)) 
      done

   if [[ -d out ]]; then
      mv out out.1
   fi
}

clean()
{
   rm -fr out* cable src/qsj.j src/*out qs* bu/  
   rm -f build  ../core nohup.out *out *csv *txt .qu  
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




