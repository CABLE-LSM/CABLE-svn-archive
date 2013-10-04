

pro prit, vname,vold, vnew 

   filepr = 'ancil.txt'
   openw,1,filepr,/append
   fileprd = 'ancil.txt.diff'
   openw,2,fileprd,/append

;   openw,3,'in.txt',/append
;   openw,4,'out.txt',/append
;   printf,3,vold
;   printf,4,vnew
;   printf,3,''
;   printf,3,''
;   printf,4,''
;   printf,4,''
;   close,3
;   close,4
   ;---lookup is 2D array [64,ntiles]
   if vname eq 'lookup'  then begin
      tnew=vnew[1,*]
      told=vold[1,*]
   endif else begin
      tnew=vnew
      told=vold
   endelse
   nold= 0LL
   mnew= 0LL
   nold = n_elements(told)
   mnew = n_elements(tnew)
   n=nold
   m=mnew
   if n ne m then begin 
      print, 'size of ', vname, '  differs !' 
      print, 'size of old: ', nold
      print, 'size of new: ', mnew
      print, 'proceed assuming size of old: although misses data'
      print, 'restructure old'
   endif


if vname ne 'data'  then begin
    
   newold = vnew

   if vname eq 'lookup'  then begin
      for i=0,63 do begin
         for j=0, m-1 do begin
            if j lt n  then begin
               newold[i,j] = vold[i,j]
            endif else begin
               newold[i,j] = 1394931 
            endelse  
         endfor
      endfor
   endif else begin
      for i=0L,m-1 do begin
            if i lt n then begin
                newold[i] = vold[i]
            endif else begin
               newold[i] = 1394931 
            endelse  
      endfor
   endelse
   
   if vname eq 'lookup'  then begin
      printf,1, vname, ', old, new', n, '* 64'
      printf,1,'-----------------------------------' 
      printf,2, vname, ', old, new', n, '* 64'
      printf,2,'-----------------------------------' 
      for i=0,63 do begin
         for j=0, m-1 do begin
            printf,1, newold[i,j], vnew[i,j]
            if newold[i,j] ne vnew[i,j] then begin
               printf,2, i, j, newold[i,j], vnew[i,j]
            endif   
         endfor
      endfor
   endif else begin
      printf,1, vname, ', old, new', n
      printf,1,'-----------------------------------' 
      printf,2, vname, ', old, new', n
      printf,2,'-----------------------------------' 

      for i=0L,m-1 do begin
         printf,1,i 
         printf,1, newold[i], vnew[i]
            if newold[i] ne vnew[i] then begin
               printf,2, i, newold[i], vnew[i]
            endif   
      endfor

   endelse

   print,''   
   printf,1,''   
   printf,1,''  
   printf,2,''  
   printf,2,''  
   close, 1
   close, 2
endif

end

;==============================================================================
pro prancil, fixhdrin, fixhdrout, intcin, intcout, lookupin, lookupout, $
            realcin,levdepcin,rowdepcin,coldepcin, $
            fieldcin,extracin,temphistin,cfi1in,cfi2in,cfi3in, $
            realcout, levdepcout, rowdepcout, coldepcout, $
             fieldcout,extracout,temphistout,cfi1out,cfi2out,cfi3out, datain,data
      prit, 'fixhdr', fixhdrin, fixhdrout
      prit, 'intc', intcin, intcout
      prit, 'realc', realcin, realcout
      prit, 'levdepc', levdepcin, levdepcout
      prit, 'rowdepc', rowdepcin, rowdepcout
      prit, 'coldepc', coldepcin, coldepcout
      prit, 'fieldc', fieldcin, fieldcout
      prit, 'extrac', extracin, extracout
      prit, 'temphist', temphistin, temphistout
      prit, 'cfi1', cfi1in, cfi1out
      prit, 'cfi2', cfi2in, cfi2out
      prit, 'cfi3', cfi3in, cfi3out
      prit, 'lookup', lookupin, lookupout
      prit, 'data', datain, data 

end

;==============================================================================


