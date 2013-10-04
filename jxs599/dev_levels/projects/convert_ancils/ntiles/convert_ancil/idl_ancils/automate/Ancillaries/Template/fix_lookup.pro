
pro fix_lookup, lookupin, lookupout, ntiles, nlat, nlon
   lookupout = lon64arr(64,ntiles)
   ;lookupout = lookupin
   for i=0, 63 do begin ;std # elements in dim=1 for lookup
      for j=0, ntiles-1 do begin
         case i of
            14: begin
               ;lookupout[i,j]=long(nlat)*long(nlon)
               lookupout[i,j]=nlat*nlon
            end   
            17: begin
               lookupout[i,j]=long(nlat)
            end   
            18: begin
               lookupout[i,j]=long(nlon)
            end   
           ;22: 1391 etc check this code in each case 
            28: begin
               if j eq 0 then begin
                  lookupout[i,j]=long(2048)   ;why 2048
               endif else begin
                  lookupout[i,j]= lookupout[i,j-1] + long(4096) ;is 8192(=disklen) always the case
               endelse   
            end   
            ;this sets disklen,th. data1 alloc. -> data assign.
            29: begin
               lookupout[i,j]=long(4096)
            end   
            39: begin
               fbase = long(nlat) * long(nlon)
               if j eq 0 then begin
                  lookupout[i,j]=long(1)
               endif else begin
                  lookupout[i,j]= lookupout[i,j-1] + fbase 
               endelse   
            end   
            42: begin
               lookupout[i,j]=j+1
            end   
            else: begin
               lookupout[i,j]=lookupin[i,0]
            end   
         endcase
      endfor
   endfor
end
;==============================================================================

