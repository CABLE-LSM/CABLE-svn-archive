
;---subr. to modify old ancil. fix_llokup in separate file

pro fix_hdr, fixhdrin, fixhdrout, ntiles, nlat, nlon
   fixhdrout=fixhdrin
   fixhdrout(151)= long(ntiles) 
   fixhdrout(160)= long(ntiles) * long(nlat) * long(nlon) 
end

;==============================================================================

pro fix_intc, intcin, intcout, ntiles
   intcout=intcin
;   intcout(7)= long(1) ;atmos. levels i don't know why i thought i need to do this
   intcout(14)= long(ntiles) 
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
            temphistout,cfi1out,cfi2out,cfi3out,lookupout,data, npft 

   ;---read in std ancil we  want to reproduce
   ;---fixhdrin avail. in this scope via arg. dec., mem. alloc. in rdancil

   rdancil, filein, fixhdrin,intcin,realcin,levdepcin,rowdepcin,coldepcin, $
         fieldcin,extracin,temphistin,cfi1in,cfi2in,cfi3in,lookupin,datain,$
         /swap
      
   ;---make changes to header

   fix_hdr, fixhdrin, fixhdrout, ntiles, nlat, nlon

   fix_intc, intcin, intcout, ntiles

   fix_lookup, lookupin, lookupout, ntiles, nlat, nlon

  fix_others, realcin,levdepcin,rowdepcin,coldepcin, $
            fieldcin,extracin,temphistin,cfi1in,cfi2in,cfi3in, $
           realcout, levdepcout, rowdepcout, coldepcout, $
             fieldcout,extracout,temphistout,cfi1out,cfi2out, cfi3out
   
   ;---re-write data for ntiles
   nrec=(size(lookupin))[2]
   rewritedata, datain, data, nrec, ntiles, nlat, nlon, npft
;
;   ;---print ancil stuff to files
  
;   prancil, fixhdrin, fixhdrout, intcin, intcout, lookupin, lookupout, $
;            realcin,levdepcin,rowdepcin,coldepcin, $
;            fieldcin,extracin,temphistin,cfi1in,cfi2in,cfi3in, $
;            realcout, levdepcout, rowdepcout, coldepcout, $
;             fieldcout,extracout,temphistout,cfi1out,cfi2out, cfi3out, datain, data
;   
;;help
;;---write ancil
;;-----------------------------------------------------
;   wrancil, fileout, fixhdrin,intcin,realcin,levdepcin,rowdepcin,coldepcin, $
;         fieldcin,extracin,temphistin,cfi1in,cfi2in,cfi3in,lookupin,data,$
;         /swap
;print,'Writing out ancillary file: ',fileout

   wrancil, fileout, fixhdrout,intcout,realcout,levdepcout,rowdepcout, $
           coldepcout, fieldcout,extracout,temphistout,cfi1out,cfi2out,$
           cfi3out,lookupout, data, /swap;, /makewellformed
end


