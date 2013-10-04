
;---subr. to modify old ancil. fix_llokup in separate file

pro fix_hdr, fixhdrin, fixhdrout, ntiles, nlat, nlon, npft,nfields_updated, ntimes,  nfields_init
   fixhdrout=fixhdrin
   ;see tech. paper pof3.pdf - this corresponds to the 2nd dimnsion in the lookup table -> # fields
   ;therefore for vars with ntiles initialized only once,   fixhdrout(151)= long(ntiles)
   ;for veg_func ancil this is (5*12) + (5*12) + 12, npfts*times updated (12 * over a year)
   fixhdrout(151)= long64( nfields_updated*(npft * ntimes) + ( nfields_init* ntimes) )
   fixhdrout(152)= long64( nfields_updated*(npft * ntimes) + ( nfields_init* ntimes) )
   ;dimension of data in full,   fixhdrout(160)= long(ntiles) * long(nlat) * long(nlon)
   ;for veg_func ancil this is ( (5*12) + (5*12) + 12 ) * 7008(=nlat*nlon)
   fixhdrout(160)= long64( fixhdrout(151) * nlat * nlon )
   ;fixhdrout(161)= fixhdrout(160) 
end

;==============================================================================

pro fix_intc, intcin, intcout, ntiles, npft, nfield_pft
   intcout=intcin
   ;# diff. fields in dump, for ntiles intcout(14)= long(ntiles)
   ;here 1 = the remainig field wich is independent of field type
   ;intcout(14)= long64( ( nfield_pft * npft ) + 1 )
   intcout(14)= long64( 0 )
end

;==============================================================================

pro   fix_others, realcin,levdepcin,rowdepcin,coldepcin, $
            fieldcin,extracin,temphistin,cfi1in,cfi2in,cfi3in, $
            realcout, levdepcout, rowdepcout, coldepcout, $
             fieldcout,extracout,temphistout,cfi1out,cfi2out, cfi3out
   realcout = realcin
   levdepcout = levdepcin
   rowdepcout =  rowdepcin
   coldepcout =  coldepcin
   fieldcout =  fieldcin
   extracout =  extracin
   temphistout =  temphistin
   cfi1out =  cfi1in
   cfi2out =  cfi2in
   cfi3out =  cfi3in
end

;==============================================================================

pro    ntiles_anc, filein, fileout, ntiles, nlat, nlon, $
            fixhdrin,intcin,realcin,levdepcin,rowdepcin,coldepcin,fieldcin,extracin, $
            temphistin,cfi1in,cfi2in,cfi3in,lookupin,datain, $
            fixhdrout,intcout,realcout,levdepcout,rowdepcout,coldepcout,fieldcout,extracout, $
            temphistout,cfi1out,cfi2out,cfi3out,lookupout,data, npft,nfields_updated, ntimes,$
            nfields_init, nfield_pft, opft, runflag
   ;---read in std ancil we  want to reproduce
   ;---fixhdrin avail. in this scope via arg. dec., mem. alloc. in rdancil
   
   rdancil, filein, fixhdrin,intcin,realcin,levdepcin,rowdepcin,coldepcin, $
         fieldcin,extracin,temphistin,cfi1in,cfi2in,cfi3in,lookupin,datain, runflag, $
         /swap

   if (runflag ne 'readonly') then begin
     ;---make changes to header
   
      fix_hdr, fixhdrin, fixhdrout, ntiles, nlat, nlon, npft,nfields_updated, ntimes,  nfields_init
   
      fix_intc, intcin, intcout, ntiles, npft, nfield_pft
   
      fix_lookup, lookupin, lookupout, ntiles, nlat, nlon, npft,nfields_updated, ntimes,  nfields_init
   
      fix_others, realcin,levdepcin,rowdepcin,coldepcin, $
               fieldcin,extracin,temphistin,cfi1in,cfi2in,cfi3in, $
               realcout, levdepcout, rowdepcout, coldepcout, $
                fieldcout,extracout,temphistout,cfi1out,cfi2out, cfi3out
      
      ;---re-write data for ntiles
      nrec=(size(lookupin))[2]
      rewritedata, datain, data, nrec, ntimes, nlat, nlon, npft, opft
      
      ;---print ancil stuff to files
      prancil, fixhdrin, fixhdrout, intcin, intcout, lookupin, lookupout, $
              realcin,levdepcin,rowdepcin,coldepcin, $
              fieldcin,extracin,temphistin,cfi1in,cfi2in,cfi3in, $
              realcout, levdepcout, rowdepcout, coldepcout, $
               fieldcout,extracout,temphistout,cfi1out,cfi2out, cfi3out, datain, data

      ;---write ancil
      ;-----------------------------------------------------
      print,'Writing out ancillary file: ',fileout
   
      wrancil, fileout, fixhdrout,intcout,realcout,levdepcout,rowdepcout, $
             coldepcout, fieldcout,extracout,temphistout,cfi1out,cfi2out,$
             cfi3out,lookupout, data, npft, ntimes, /swap;, /makewellformed
   endif
      
end


