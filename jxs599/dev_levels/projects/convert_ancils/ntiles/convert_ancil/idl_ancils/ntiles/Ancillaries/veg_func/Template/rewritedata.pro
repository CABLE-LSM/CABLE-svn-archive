pro rewritedata, data, redata,nrec, ntimes, nlat, nlon,npft,opft
   ;here we want to separate into npft per timestep
   fpft = long64 ( opft + opft + 1) 
   newfpft = long64 ( npft + npft + 1 ) 
   temp =  long64( ( npft *ntimes ) +  ( npft * ntimes )  +  ntimes ) 
   length = long64( nlat * nlon )
   ftemp= long64( temp*length ) 

   redata= dblarr(ftemp)
   datatile =  dblarr(ntimes,fpft,length)
   redatatile = dblarr(ntimes,newfpft,length)

   for i=0LL,ntimes-1LL do begin
      for j=0LL,fpft-1LL do begin
         for k=0LL,length-1LL do begin
            kk = (i*fpft*length) + (j*length) + k
            datatile[i,j,k]=data[kk]
         endfor    
      endfor    
   endfor  

   for i=0LL,ntimes-1LL do begin
      for j=0LL,newfpft-1LL do begin
         for k=0LL,length-1LL do begin
            ;---first of new tiles is equal to npfts
            if(j lt opft) then begin  
               redatatile[i,j,k]=datatile[i,j,k]
            endif
            if(j ge opft AND j lt (2*opft) ) then begin  
               redatatile[i,j,k]=datatile[i,j-opft,k]
            endif
            if(j ge (2*opft) AND j lt (npft) ) then begin  
               redatatile[i,j,k]=0.
            endif
            if(j ge (npft) AND j lt (npft+opft) ) then begin  
               jdash= j - (npft-opft)
               redatatile[i,j,k] = datatile[i,jdash,k]
            endif
            if(j ge (npft+opft) AND j lt (npft+(2*opft)) ) then begin  
               jdash=j - (npft-opft)
               redatatile[i,j,k] = datatile[i,jdash-opft,k]
            endif
            if(j ge (npft+(2*opft)) AND j lt (2*npft) ) then begin  
               redatatile[i,j,k] = 0.
            endif
            if(j ge 2*npft ) then begin  
               redatatile[i,j,k] = datatile[i,fpft-1,k]
            endif
         endfor   
      endfor  
   endfor  

   for i=0LL,ntimes-1LL do begin
      for j=0LL,newfpft-1LL do begin
         for k=0LL,length-1LL do begin
            kk = (i*newfpft*length) + (j*length) + k
            redata[kk]=redatatile[i,j,k]
         endfor    
      endfor    
   endfor  
   
end
















