
pro fix_lookup, lookupin, lookupout, ntiles, nlat, nlon, npft,nfields_updated, ntimes,  nfields_init, nfield_pft
   ;=132 for npft=5 original
   ftemp=  (nfields_updated*npft * ntimes) + ( nfields_init* ntimes) 
;   print, 'fix_lookup:ftemp ', ftemp
   lookupout = lon64arr(64,ftemp)
   stj= intarr(ftemp)
   ;see tech. paper pof3.pdf - this corresponds to the 2nd dimnsion in the lookup table -> # fields
   ;therefore for vars with ntiles initialized only once,   fixhdrout(151)= long(ntiles)
   ;for veg_func ancil this is (5*12) + (5*12) + 12, npft*times updated (12 * over a year)
   ;different fields in init. update
   ndiff =  (nfields_updated*npft) +  nfields_init 
   ;lookupout = lookupin
   laicode=long64(217)
   htcode=long64(218)
   condcode=long64(213)
   recctr=0
   ;--- 
   for i=0, ntimes-1 do begin
      for j=0, ndiff-1 do begin
         if j lt npft then begin 
            stj[recctr] = 1
         endif   
         if ( j ge npft AND j lt (2*npft) ) then begin
            stj[recctr] = 2
         endif   
         if j ge (2*npft)then begin 
            stj[recctr] = 3
         endif   
         recctr = recctr+1
      endfor
   endfor

   for i=0, 63 do begin ;std # elements in dim=1 for lookup
      for j=0, ftemp-1 do begin
   
         case i of

            1: begin
               fcondold=0
               for k=1, ntimes do begin
                  fcond = k * ndiff
                  if j lt fcond AND j ge fcondold then begin
                     lookupout[i,j]=long64(k)
                     fcondold=fcond
                  endif 
               endfor
            end  
             
            7: begin
               fcondold=0
               for k=1, ntimes do begin
                  fcond = k * ndiff
                  if j lt fcond AND j ge fcondold then begin
                     lookupout[i,j]=long64(k+1)
                     if k ge ntimes then begin
                        lookupout[i,j]=long64(1)
                     endif           
                     fcondold=fcond
                  endif 
               endfor
            end   

            14: begin
               lookupout[i,j]=long64(nlat*nlon)
            end   
            17: begin
               lookupout[i,j]=long64(nlat)
            end   
            18: begin
               lookupout[i,j]=long64(nlon)
            end   
            22: begin
               case stj[j] of
                  1: lookupout[i,j]=long64(1392) 
                  2: lookupout[i,j]=long64(1393) 
                  3: lookupout[i,j]=long64(1384) 
               endcase
            end   
           ;22: 1391 etc check this code in each case 
            28: begin
               if j eq 0 then begin
                  lookupout[i,j]=long64(21013) ;+ long64(8192)  ;why 2048
               endif else begin
                  lookupout[i,j]= lookupout[i,j-1] + long64(8192) ;is 8192(=disklen) always the case
               endelse   
            end   
            ;this sets disklen,th. data1 alloc. -> data assign.
            29: begin
               lookupout[i,j]=long64(8192)
            end   
            39: begin
               fbase =long64(nlat *nlon)
               if j eq 0 then begin
                  lookupout[i,j]=long64(1)
               endif else begin
                  lookupout[i,j]= lookupout[i,j-1] + fbase 
               endelse   
            end   
            ;stashcode of var, un this case (217(*5), 218(*5), 213) (*12)  
            41: begin
               case stj[j] of
                  1: lookupout[i,j]= laicode
                  2: lookupout[i,j]= htcode
                  3: lookupout[i,j]= condcode
               endcase
            end   
            ;in this ancil case we can set this to 0
            42: begin
               lookupout[i,j]=long64(0)
            end   
            else: begin
                    ;we can do this here as all=8192
               lookupout[i,j]=lookupin[i,0]
            end   
         endcase
      endfor
   endfor
end
;==============================================================================















