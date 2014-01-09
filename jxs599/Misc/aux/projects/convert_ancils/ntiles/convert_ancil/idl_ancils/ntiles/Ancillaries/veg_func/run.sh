#!/bin/csh

set k = 'vegfunc'
set Ntiles = 9 #irrelevant here
set opft = 5 
set npft = 13 

set exitFLAG = 1

set runFLAG = 'fullrun'
set filein = 'qrparm.veg.func_seas'

#set runFLAG = 'readonly'
#set filein = 'vegfunc_13'

echo ''
echo 'this run is ' $runFLAG
echo ''
mkdir ancils 
mkdir $k 
cp Template/* $k/
if $filein == 'vegfunc_13' then
   cp vegfunc_13 $k
endif   
cd $k
      
set fileout = $k'_'$npft

if $filein == 'qrparm.veg.func_seas' then
   ln -s ~/data/$filein
endif
      
echo 'nmain, "'"$filein"'", "'"$fileout"'", '$Ntiles', fixhdrin,intcin,realcin,levdepcin,rowdepcin,coldepcin,fieldcin,extracin, temphistin,cfi1in,cfi2in,cfi3in,lookupin,datain,fixhdrout,intcout,realcout,levdepcout,rowdepcout,coldepcout,fieldcout,extracout, temphistout,cfi1out,cfi2out,cfi3out,lookupout,data, '$npft', '$opft', "'"$runFLAG"'"' > ancil_batch
if ( $exitFLAG > 0) then
   echo 'exit' >> ancil_batch
endif
      
source ./idl_setup
setenv IDL_STARTUP ancil_batch
idl
unsetenv IDL_STARTUP


if ( $exitFLAG > 0) then
   cp $fileout ../
   cp  ancil.txt ../
   cp tow* ..
   cp ans* ..
   cp dis* ..
   cd ../
   rm -fr $k
   rm -fr ancils
endif














