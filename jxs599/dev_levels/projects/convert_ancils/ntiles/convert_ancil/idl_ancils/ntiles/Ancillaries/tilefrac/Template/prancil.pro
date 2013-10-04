

pro prit, vname,vold, vnew 

   filepr = 'ancil.txt'
   openw,1,filepr,/append
   fileprd = 'ancil.txt.diff'
   openw,2,fileprd,/append

   ;---lookup is 2D array [64,ntiles]
   if vname eq 'lookup'  then begin
      tnew=vold[1,*]
      told=vnew[1,*]
   endif else begin
      tnew=vold
      told=vnew
   endelse

   n = n_elements(told)
   m = n_elements(tnew)

   if n ne m then begin 
      print, 'size of ', vname, 'differs !' 
      print, 'size of vold: ', n
      print, 'size of new: ', m
      print, 'proceed assuming size of old'
   endif 
   
   if vname eq 'lookup'  then begin
      printf,1, vname, ', old, new', n, '* 64'
      printf,1,'-----------------------------------' 
      printf,2, vname, ', old, new', n, '* 64'
      printf,2,'-----------------------------------' 
      for i=0,63 do begin
         for j=0, n-1 do begin
            printf,1, vold[i,j], vnew[i,j]
            if vold[i,j] ne vnew[i,j] then begin
               printf,2, i, j, vold[i,j], vnew[i,j]
            endif   
         endfor
      endfor
   endif else begin
      printf,1, vname, ', old, new', n
      printf,1,'-----------------------------------' 
      printf,2, vname, ', old, new', n
      printf,2,'-----------------------------------' 

      for i=0L,n-1 do begin
         printf,1,i 
         printf,1, vold[i], vnew[i]
            if vold[i] ne vnew[i] then begin
               printf,2, i, vold[i], vnew[i]
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


