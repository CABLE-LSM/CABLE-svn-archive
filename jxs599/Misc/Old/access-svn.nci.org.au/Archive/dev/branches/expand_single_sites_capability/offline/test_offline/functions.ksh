#!/bin/ksh

build_cable()
{      
   if [[ -f cable ]]; then
      rm -f cable
   fi  

   mypwd=`pwd`
   cd ../build
      
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
   cd src/
   
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
   #qsubed job hangs on this, removal of .wait at end of run allows qsub of next site
   if [ -n "$PBS_JOBID" ]; then
      rm -f ../../.wait
   fi    
}




#--- writes welcome banner to screen outlining basic operation of run.ksh
banner_welcome()
{
   banner_welcome_txt
   banner_options
   quit_option   

   typeset -u response 
   read response 
   
   if [[ $response == 'C' ]]; then
      run_tut
   elif [[ $response == 'H' ]]; then
      cheat_sheet 
   elif [[ $response == 'CC' ]]; then
      config_run
   else
      quit_response
   fi 
} 
 
banner_welc_sub()
{
   print "\n\t1. Building a CABLE executable from source code"       
   print "\n\t2. Running CABLE over each of the sites specified in sites_main.nml"  
   print "\n\t3. Plotting flux data computed by CABLE for these sites"
}



banner_welcome_txt()
{   
   print "\nThis script is currenntly capable of:" 

   banner_welc_sub

   print "\n./run.ksh without any additional arguments will print this banner and offer the following choices:" 
}



banner_options()
{
   print "\nIf you wish to continue with this interactive configuration, which is kind of a "'"CABLE tutorial"'", enter "'"C"'"" 
   print "\nAlternatively, you might find it easier to configure the behaviour of run.ksh via command line arguments," 
   print "bypassing this interactive process." 
   help_option 
   config_option 
}


 
run_tut()
{
   clear     
   print "\nAs we said before, this script is basically capable of:" 
   
   banner_welc_sub

   print "\nIt can do all of these things one after the other, or you can choose to do a combination thereof. "
   print "Additionally this can be achieved from the command line as well, and once familiar with the various "
   print "incantations of ./run.ksh you will probably not need to come back to this "'"tutorial"'" "
   
   help_option
   config_option
   print "\nNote however that you can com back to this configuration going through help pages first anyway.\n"   
   quit_option   
   
   typeset -u response 
   read response 
 
   if [[ $response == 'H' ]]; then
      run_help
   elif [[ $response == 'CC' ]]; then
      config_run
   else
      quit_response 
   fi 

}

config_option()
{
   print "\nIf you already know what you want to do you can probably go ahead and configure run.ksh here by entering "'"CC"'"."
}        
 
  
config_run()
{
   print "\nWhich option should i proceed with:\n"
   print "\n1. ./run.ksh build "
   print "2. ./run.ksh build run"
   print "3. ./run.ksh run "
   print "4. ./run.ksh plot "
   print "5. ./run.ksh run plot "
   print "6. ./run.ksh all "
   
   read response
   
   if [[ $response == '1' ]]; then
      ./run.ksh build 
   fi
  
   if [[ $response == '2' ]]; then
      ./run.ksh build run
   fi
  
   if [[ $response == '3' ]]; then
      ./run.ksh run 
   fi
  
   if [[ $response == '4' ]]; then
      ./run.ksh plot
   fi
  
   if [[ $response == '5' ]]; then
      ./run.ksh run plot 
   fi
  
   if [[ $response == '6' ]]; then
      ./run.ksh all 
   fi

}

run_help()
{
   clear      
   #jhan:later can use a2ps to make this printable
   print "\nFor a quick and nasty cheat sheet, enter "'"C"'" "
   print "\nTo learn more about the build process, enter "'"B"'" "
   print "\nTo learn more about the run process, enter "'"R"'" "
   print "\nTo learn more about the plot process, enter "'"P"'" \n"
   
   config_option   
   quit_option   
   
   typeset -u response 
   read response 
 
   if [[ $response == 'C' ]]; then
      cheat_sheet
   elif [[ $response == 'B' ]]; then
      help_build
   elif [[ $response == 'R' ]]; then
      help_run
   elif [[ $response == 'P' ]]; then
      help_plot
   elif [[ $response == 'CC' ]]; then
      config_run
   else
      quit_response
   fi 
 
}  

roundnaround()
{
   clear
   help_option
   config_option
   quit_option

   typeset -u response 
   read response 
 
   if [[ $response == 'H' ]]; then
      run_help
   elif [[ $response == 'CC' ]]; then
      config_run
   else
      quit_response
   fi
}


cheat_sheet()
{
   #jhan: can add more later - including qsub etc
   clear
   print "\nNAME: run.ksh - optionally runs/builds/plots CABLE offline"  
   print "\nSYNOPSIS: run.ksh [OPTION1] [OPTION2]\n "
   print "\nrun.ksh accepts at most two arguments, all of which are pretty obvious\n"
   print "\n./run.ksh build "
   print "\t\tBuild CABLE only. Useful in code development, and also to create an executable for a qsub job."
   print "\n./run.ksh build run"
   print "\t\tBuild CABLE and then run it."
   print "\n./run.ksh run "
   print "\t\tCABLE executable already exists, just run it."
   print "\t\t\tUseful if you are switching data sets, or have a pre-built special version of CABLE"
   print "\n./run.ksh plot "
   print "\t\tPlot existing CABLE data (must be in out/)."
   print "\n./run.ksh run plot "
   print "\t\tRun existing CABLE executable and plot this new data (will naturally be in out/)."
   print "\n./run.ksh all "
   print "\t\tDo everything. Build CABLE, then run it, then plot it."
   print "\nHit enter to continue"
   read dummy 
  
   roundnaround

}   

help_build()
{
   print "\nunder construction help build" 
   print "\n\nHit enter to continue"
   read dummy 
   roundnaround
}
   
help_run()
{
   print "\nunder construction help run" 
   print "\n\nHit enter to continue"
   read dummy 
   roundnaround
}
   
help_plot()
{
   print "\nunder construction help plot" 
   print "\n\nHit enter to continue"
   read dummy 
   roundnaround
}
   
#   print      "   4. Plotting, flux data is plotted and left in out/*sitename*. "
#   print      "   3. move output data into a new directory (out/*sitename*)"   
#   print      "     NB. plotted data compares the run version of CABLE (out_cable.nc),"  
#   print      "         the previous version of CABLE (old_cable.nc), and observations" 
#   print      " "
#   print      "Hit enter to proceed, enter Q to quit, H for more info "
#   print      " "

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

quit_option()   
{
   print "\nIf you wish to quit now, just hit enter" 
}  

quit_response()
{      
      print "\nAdios"
      exit
}


help_option()
{
   print "\nIf you wish to see more information on command line arguments to run.ksh , enter "'"H"'" " 
} 




