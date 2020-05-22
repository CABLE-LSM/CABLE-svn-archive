#!/bin/bash

rundir=$1
runtype=$2

if [ -z "$runtype" ]; then

   echo "runtype (sli,ssgw,default) not set"
   echo "Assuming you want ssgw, it is the best!"
   runtype=ssgw

fi

mkdir $rundir

cp run_model.sh ${rundir}/run_${rundir}.sh
cp ${runtype}_create_cable-nml.sh ${rundir}/.

cd $rundir

ln -s ${runtype}_create_cable-nml.sh create_cable-nml.sh

if [ -f README ]; then

   echo "README made by ${USER} to test CMIP6-GM2 on $(date '+%d-%m-%Y')" >> README

fi

summary="Ran  $rundir which is $runtype on $(date '+%d-%m-%Y  %H:%M') "

echo $summary >> README



