#!/bin/ksh 

###Users MAY change this section

# Most NCI users are restricted to 8 qsub processes at a time. $nproc limits the
# number of single-sites qsubbed at a time. the remaining sites listed in 
# sites_main.nml are qsubed automatically as an earlier job finishes, meaning 
# that this script has a maximum of $nproc jobs running at a time.
# if this script is qsubed via "run.ksh qsub" as reccomended, then this number is $nproc+1

nproc="8"

qs="qsub_"
qs_filename_base='qsmain_'
qs_filename_suffix='.nml'
###END - Users MAY change this section

#if [ -n "$PBS_JOBID" ]; then
#   cd $PBS_O_WORKDIR
#fi
. src/functions.ksh

cd ../
   if [[ -e core ]]; then
      rm -f core
   fi 
   ln -s ../core
cd test_offline/

if [[ -e build ]]; then
   rm -f build 
fi 
ln -s ../build

if [[ -d $qs"1" ]]; then

   if [[ -d bu ]]; then
      print'\nDirectory bu/ already exists, implying you have qsub-ed CABLE 
      over singlesites at least twice before. Decide what you want to keep, 
      move it out of the way, rm qs*, and  execute ./run.ksh again'
      exit
   else
      print '\nqsub_* directories already exist, implying you have qsub-ed CABLE over singlesites before. These will be moved to a back up directory
         bu/ for you, but just this time. Next time you will have to move
         it yourself'
      mkdir bu
      mv qs* bu
  fi 
fi   

cd src/
module add R
R CMD BATCH --slave Qsub.R
cd ../

if [[ $nproc -lt 10 ]]; then
   nch=1
elif  [[ $nproc -gt 9 ]] && [[ $nproc -lt 100 ]]; then
   nch=2
else
   print "\nWOW! Thousands or more processes hey?\n"
   exit
fi

uid=`whoami`


for i in `cat src/qsj.j`; do 

   qsub_avail
    
done

  



