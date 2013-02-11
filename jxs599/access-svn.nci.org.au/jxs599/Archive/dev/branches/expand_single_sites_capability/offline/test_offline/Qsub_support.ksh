#!/bin/ksh 
#set -x

#PBS -l walltime=2700
#PBS -l vmem=3072MB
#PBS -l ncpus=1
#PBS -l jobfs=2GB
#PBS -q express 
#PBS -p 10
#PBS -j oe
#PBS -N wrapper 


###Users MAY change this section

# Most NCI users are restricted to 8 qsub processes at a time. $nproc limits the
# number of single-sites qsubbed at a time. the remaining sites listed in 
# sites_main.nml are qsubed automatically as an earlier job finishes, meaning 
# that this script has a maximum of $nproc jobs running at a time.
# if this script is qsubed via "run.ksh qsub" as reccomended, then this number is $nproc+1
nproc='6' 

qs="qsub_"
qs_filename_base='qsmain_'
qs_filename_suffix='.nml'
###END - Users MAY change this section

if [ -n "$PBS_JOBID" ]; then
   cd $PBS_O_WORKDIR
fi

cd ../
   if [[ -e core ]]; then
      rm -f core
   fi 
   ln -s ../core
cd test_offline/


if [[ -d $qs$i ]]; then
   #jhan 
   if [[ -d bu ]]; then
      echo 'directory bu/ already exists. Decide what you want to keep, 
         move it out of the way, rm qs*, and  execute ./run_qs.csh again'
      exit
   else
      echo $qs$i, ' already exists. This will be moved to a back up directory
         bu/ for you, but just this time. Next time you will have to move
         it yourself'
      mkdir bu
      mv qs* bu
  fi 
fi   


module add R
R CMD BATCH --slave Qsub.R

   
for i in `cat qsj.j`; do 

   if [[ -f cable ]] || [[ $i == 1 ]]; then
      echo 'Using the cable executable present'     
   else   
      echo 'Before qsubing we need to build a cable executable'     
   fi
   
    
   touch .wait
   if [[ $i < $nproc ]] || [[ $i == $nproc ]] ; then
      
      mkdir $qs$i
      cp * $qs$i
      cp -r src/ $qs$i
      cd $qs$i
      ln -s ../data
      mv $qs_filename_base$i$qs_filename_suffix main.nml
      cat plot_main.nml >> main.nml
      
      /opt/pbs/bin/qsub -S /bin/ksh run.ksh       
      
      cd ../
     
      if [[ $i < $nproc ]]; then
         qspace=true
      else   
         qspace=false
      fi
         
   else
      qspace=false
   fi

  
   while ( ! ($qspace) ) 
   do
      if [[ -f .wait ]]; then
         qspace=false
      else
         qspace=true
      fi
   done  
   
   if [[ $i > $nproc ]]; then
      mkdir $qs$i
      cp * $qs$i
      cp -r src/ $qs$i
      cd $qs$i
      ln -s ../data
      mv $qs_filename_base$i$qs_filename_suffix main.nml
      cat plot_main.nml >> main.nml
      
      /opt/pbs/bin/qsub -S /bin/ksh run.ksh       
      
      cd ../
   fi 

     
done





