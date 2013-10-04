pro rewritedata, data, redata,nrec, ntiles, nlat, nlon,npft

   length = double(nlat * nlon)
   ch = replicate(0,ntiles)
   datatile =  dblarr(nrec,length)
   redatatile = dblarr(ntiles,length)
   temp = double(ntiles * length)
   redata= dblarr(temp)
   oldntiles=nrec

   for i=0L,nrec-1 do begin
      for j=0L,length-1 do begin
         jj = (i*length) + j
        datatile[i,j]=data[jj]
      endfor    
   endfor  
  

   for i=0L,ntiles-1 do begin
      for j=0L,length-1 do begin
         ;---first of new tiles is equal to npfts
         if(i lt npft) then begin  
            redatatile[i,j]=datatile[i,j]
         endif

         ;---extra tiles set to zero
         if( i ge npft AND i lt (npft + ntiles - oldntiles) ) then begin 
            redatatile[i,j] = 0.0 
         endif
         
         ;---rewrite ancil. w last matching orig. 1st 5tiles repeated twice - BUT
         ;---half the frac. (net result should be same) 
         ;-----------------------------------------
         ;---1st 5 tiles repeated in postions 1:5 
         ;if(i <= 5) Lnewfield[,,i] = 0.5 * rdata$field[,,i] 
         ;---1st 5 tiles repeated in postions 6:10
         ;if(i > 5 & i <= 10 ) Lnewfield[,,i] = 0.5 * rdata$field[,,i-5] 
         ;---tiles in postions 11:13 set to 0 (i.e. they have no effect)
         ;if(i > 10 & i <= 13 )  Lnewfield[,,i] = 0.0 
                                                                            
         ;---tiles13:17 set to original 5:9(non-veg tiles)
         if(i ge (npft + ntiles - oldntiles) ) then begin 
            redatatile[i,j] = datatile[i-(ntiles - oldntiles),j] 
            if(ch[i]<1) then begin
               print,'datatile(13:16) = ',redatatile[i,j]
               ch[i]=2
            endif
         endif

      endfor    
   endfor  


   for i=0L,ntiles-1 do begin
      for j=0L,length-1 do begin
         jj = (i*length) + j
        redata[jj]=redatatile[i,j]
      endfor    
   endfor  
  
  
  
    
end
















