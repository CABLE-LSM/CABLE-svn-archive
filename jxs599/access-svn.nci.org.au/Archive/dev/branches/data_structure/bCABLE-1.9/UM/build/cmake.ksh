#!/bin/ksh


# build CABLE
bld_libcable()
{
   ./bld_libcable-UM
}

## build UM with new CABLE library
bld_UM()
{
   cd ~/UM_ROUTDIR/jxs599/laaar/
      rm -f ummodel/bin/laaar.exe   
      ./bld.csh
      /bin/cp -f ummodel/bin/laaar.exe ~/UM_ROUTDIR/jxs599/laaaf/ummodel/bin/laaaf.exe
   cd ../
   
   clear
   ls -alt ~/CABLE-UM/libcable.a
   echo ''
   
   ls -alt ~/UM_ROUTDIR/jxs599/laaaf/ummodel/bin/laaaf.exe
   echo ''
}
   ## prepare and run UM manually
run_UM()
{
   cd ~/umui_runs
   umuirunID="laaaf-117125231" #16 nodes
   #umuirunID="laaaf-319143631" #1 node - express queue - mem limitations
   #umuirunID="laaaf-321120720" #1 node - normal queue
   
   ## rm existing jobscripts
   if [[ -d $umuirunID ]]; then
      rm -fr $umuirunID/
   fi 
   
   ## copy jobscripts
   /bin/cp -r E$umuirunID/ $umuirunID
   cd $umuirunID
    
   ## rm existing testing output 
#   if [[ -f ~/UM_ROUTDIR/laaaf/canopy_flux00.bin ]]; then
      rm -fr ~/UM_ROUTDIR/laaaf/ts*
      rm -fr ~/UM_ROUTDIR/laaaf/lat*
      rm -fr ~/UM_ROUTDIR/laaaf/lon*
#   fi

   ## run UM manually
   ./FCM_MAIN_SCR
}





if [[ $1 == 'make' ]]; then
   make -f Make_CABLE-UM

elif [[ $1 == 'build' ]]; then

   if [[ $2 == 'cable' ]]; then
      bld_libcable
   elif [[ $2 == 'um' ]]; then
      bld_UM
   else 
      bld_libcable
      bld_UM
   fi

elif [[ $1 == 'run' ]]; then

   if [[ $2 == 'bld_um' ]]; then
      bld_UM
      run_UM
   else
      run_UM
   fi

elif [[ $1 == 'help' ]]; then
   echo './cmake           -  by default, with no arguments given, cmake ' 
   echo '                               1. makes cable; '
   echo '                               2. builds cable library;' 
   echo '                               3. builds UM with new cable library;'
   echo '                               4. if setup - runs UM($expid) '
   echo ''
   echo './cmake help      -  prints this doc'
   echo ''
   echo './cmake make      -  compile CABLE without linking. useful for syntax checking'
   echo ''
   echo './cmake build     -  by default, with no further arguments given, cmake build' 
   echo '                               1. makes cable; '
   echo '                               2. builds cable library;' 
   echo '                               3. builds UM with new cable library;'
   echo ''
   echo './cmake build cable     -  ' 
   echo '                               1. makes cable; '
   echo '                               2. builds cable library;' 
   echo '                               3. builds UM with new cable library;'
   echo ''
   echo './cmake build um        -  assumes you have the cable library already built and in place ' 
   echo '                           and proceeds to  build UM with this cable library.'
   echo ''
   echo './cmake run          -  by default, with no further arguments given, cmake run' 
   echo '                               runs the UM executable as specified in configuration at '
   echo '                               the top of this file. ' 
   echo ''
   echo './cmake run bld_um    -  assumes you have the cable library already built and in place ' 
   echo '                           and proceeds to  build UM with this cable library and the run.'


else   
  
   bld_libcable
   bld_UM
   run_UM

fi



