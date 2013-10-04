#!/bin/csh

set k = 'tilefrac'
set Ntiles = 17 
set exitFLAG = 1

mkdir ancils 
mkdir $k 
cp Template/* $k/
cd $k
      
set filein = 'qrparm.veg.frac_igbp'
set fileout = $k'_'$Ntiles
     
ln -s ~/data/$filein
      
echo 'nmain,"'"$filein"'", "'"$fileout"'", '$Ntiles', fixhdrin,intcin,realcin,levdepcin,rowdepcin,coldepcin,fieldcin,extracin, temphistin,cfi1in,cfi2in,cfi3in,lookupin,datain,fixhdrout,intcout,realcout,levdepcout,rowdepcout,coldepcout,fieldcout,extracout, temphistout,cfi1out,cfi2out,cfi3out,lookupout,data' > ancil_batch
if ( $exitFLAG > 0) then
   echo 'exit' >> ancil_batch
endif
      
source ./idl_setup
setenv IDL_STARTUP ancil_batch
idl
unsetenv IDL_STARTUP


if ( $exitFLAG > 0) then
   cp $fileout ../
   cp run.log ../
   cd ../
   rm -fr $k
   rm -fr ancils
   rm -f run.log
endif

