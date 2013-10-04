
;pro nmain'qrparm.veg.frac_igbp', 'vegfrac.nc', 'vegfrac9', 9, fixhdrin,intcin,realcin,levdepcin,rowdepcin,coldepcin,fieldcin,extracin, temphistin,cfi1in,cfi2in,cfi3in,lookupin,datain,fixhdrout,intcout,realcout,levdepcout,rowdepcout,coldepcout,fieldcout,extracout, temphistout,cfi1out,cfi2out,cfi3out,lookupout,data 
;---nmain srgs = existing ancillary, netcdf modified 
;---data, out ancil, ntiles

pro nmain, filein, fileout, ntiles,                       $ 
            fixhdrin,intcin,realcin,levdepcin,rowdepcin,coldepcin,fieldcin,extracin, $
            temphistin,cfi1in,cfi2in,cfi3in,lookupin,datain, $
            fixhdrout,intcout,realcout,levdepcout,rowdepcout,coldepcout,fieldcout,extracout, $
            temphistout,cfi1out,cfi2out,cfi3out,lookupout,data, npft, opft, runflag 

   nlat = 73
   nlon = 96 
   nfields_updated = 2
   ntimes = 12
   nfields_init = 1
   nfield_pft = 2

   print, 'at least made it to start of Ntiles_main'

   ntiles_anc, filein, fileout, ntiles, nlat, nlon, $
            fixhdrin,intcin,realcin,levdepcin,rowdepcin,coldepcin,fieldcin,extracin, $
            temphistin,cfi1in,cfi2in,cfi3in,lookupin,datain, $
            fixhdrout,intcout,realcout,levdepcout,rowdepcout,coldepcout,fieldcout,extracout, $
            temphistout,cfi1out,cfi2out,cfi3out,lookupout,data, npft,nfields_updated, ntimes,$
            nfields_init, nfield_pft, opft, runflag


   print, 'made it to the last line'
;help
;   exit
end
