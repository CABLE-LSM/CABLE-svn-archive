#!/bin/csh

set Ntiles = 17 

mkdir ancils 
foreach k ('depth' 'mass' )
      mkdir $k
      mkdir $k
      cp Template/* $k
      cd $k
      if  ( $k == 'depth') then
         set filein = 'qrparm.snowd1'
         set fileout = 'snowdepth_'$Ntiles
      endif

      if  ( $k == 'mass') then
         set filein = 'qrparm.snowm1'
         set fileout = 'snowmass_'$Ntiles
      endif

      ln -s ~/data/$filein
      
      echo 'nmain,"'"$filein"'", "'"$fileout"'", '$Ntiles', fixhdrin,intcin,realcin,levdepcin,rowdepcin,coldepcin,fieldcin,extracin, temphistin,cfi1in,cfi2in,cfi3in,lookupin,datain,fixhdrout,intcout,realcout,levdepcout,rowdepcout,coldepcout,fieldcout,extracout, temphistout,cfi1out,cfi2out,cfi3out,lookupout,data' > ancil_batch
      echo 'exit' >> ancil_batch
      
      source ./idl_setup
      setenv IDL_STARTUP ancil_batch
      idl
      unsetenv IDL_STARTUP
      cp $fileout ../ancils
      cd ../
      rm -fr $k
end


foreach k ('tsoil' 'sthf' 'smcl' )
   mkdir $k
   foreach i (1 2 3 4)
      mkdir $k/$i
      cp Template/* $k/$i
      cd $k/$i
      
      set filein = 'qrparm.'$k$i'_tile'
      set fileout = $k$i'_'$Ntiles
     
      ln -s ~/data/$filein
      
      echo 'nmain,"'"$filein"'", "'"$fileout"'", '$Ntiles', fixhdrin,intcin,realcin,levdepcin,rowdepcin,coldepcin,fieldcin,extracin, temphistin,cfi1in,cfi2in,cfi3in,lookupin,datain,fixhdrout,intcout,realcout,levdepcout,rowdepcout,coldepcout,fieldcout,extracout, temphistout,cfi1out,cfi2out,cfi3out,lookupout,data' > ancil_batch
      echo 'exit' >> ancil_batch
      
      source ./idl_setup
      setenv IDL_STARTUP ancil_batch
      idl
      unsetenv IDL_STARTUP
      cp $fileout ../../ancils
      cd ../../
   end
   rm -fr $k
end


