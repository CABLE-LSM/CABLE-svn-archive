
   ;rdancil_readit(unit,fixhdr,bytes,swap,int0 ,150,2)
FUNCTION rdancil_readit,unit,fixhdr,bytes,swap,proto,index,dims
   
;   print,'fixhdr', fixhdr[index-1]

   if fixhdr[index-1] le 0 then return,proto*0
   size=1LL
   for i=1LL,dims do size=size*fixhdr[index-1+i]
   ;answer=reform(replicate(proto,size), fixhdr[index:index+dims-1])
   rdanswer= replicate(proto,size)
   if ( fixhdr[index+dims-1] ne fixhdr[index] ) then begin
      answer= lon64arr(fixhdr[index:index+dims-1])
   endif else begin
      answer= replicate(proto,size)
   endelse    
;   print,'[index:index+dims-1]', index,index+dims-1
;   print,'fixhdr[index:index+dims-1]', fixhdr[index:index+dims-1]
;   print,'size', size
;   print,'size(answer)', size(answer)
;   fpos = (fixhdr[index-1]-1)*bytes
;   print, 'fpos  ', fpos
   point_lun,unit,(fixhdr[index-1]-1)*bytes

   readu,unit,rdanswer
   k=0L
   for j=0, fixhdr[index+dims-1] -1 do begin
      if ( fixhdr[index+dims-1] ne fixhdr[index] ) then begin
         for i=0, fixhdr[index] -1 do begin
            ;print, i, j, k
            answer[i,j] = rdanswer[k]
            k=k+1
         endfor
      endif  else begin  
         answer = rdanswer
      endelse    
   endfor

   if swap then answer=swap_endian(answer)

   if size gt 10000 then begin 
      openw,11,'answer.txt'
      for j=0, fixhdr[index+dims-1] -1 do begin
         for i=0, fixhdr[index] -1 do begin
            printf,11, i,j, answer[i,j]
         endfor
      endfor
      close,11
   endif

   return,answer
end

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------


PRO rdancil,filename, $
            fixhdr,intc,realc,levdepc,rowdepc,coldepc,fieldc,extrac, $
            temphist,cfi1,cfi2,cfi3,lookup,data,runflag, $
            swap=swap,bits=bits

   if n_elements(bits) eq 0 then bits=64
   if n_elements(swap) eq 0 then swap=0

   case bits of
      64: begin
         real0=0D0
         int0=0LL
      end
      32: begin
         real0=0.
         int0=0L
      end
   endcase
   bytes=bits/8

   openr,unit,filename,/get_lun

   fixhdr=replicate(int0,256)
   readu,unit,fixhdr

   if swap then fixhdr=swap_endian(fixhdr)
   intc    =rdancil_readit(unit,fixhdr,bytes,swap,int0 ,100,1)
   realc   =rdancil_readit(unit,fixhdr,bytes,swap,real0,105,1)
   levdepc =rdancil_readit(unit,fixhdr,bytes,swap,real0,110,2)
   rowdepc =rdancil_readit(unit,fixhdr,bytes,swap,int0 ,115,2)
   coldepc =rdancil_readit(unit,fixhdr,bytes,swap,int0 ,120,2)
   fieldc  =rdancil_readit(unit,fixhdr,bytes,swap,int0 ,125,2)
   extrac  =rdancil_readit(unit,fixhdr,bytes,swap,int0 ,130,1)
   temphist=rdancil_readit(unit,fixhdr,bytes,swap,int0 ,135,1)
   cfi1    =rdancil_readit(unit,fixhdr,bytes,swap,int0 ,140,1)
   cfi2    =rdancil_readit(unit,fixhdr,bytes,swap,int0 ,142,1)
   cfi3    =rdancil_readit(unit,fixhdr,bytes,swap,int0 ,144,1)
   lookup  =rdancil_readit(unit,fixhdr,bytes,swap,int0 ,150,2)

   ndata=0L
   nrec=(size(lookup))[2]
   for i=0L,nrec-1L do begin
       ndata=ndata+(lookup[14,i]>0)
   ;    print, 'ndata  ', ndata
   endfor

   if (runflag ne 'readonly') then begin
      data=dblarr(ndata)
      ;jhan{
      ;mydlen = 4096L
      ;mylen = 7008L
      ;ndata = mylen * 9 
      ;;}
      ;start=2048L
      
      for i=0,nrec-1 do begin
        start=lookup[28,i]
         if (start ne -99L) then begin
           disklen=lookup[29,i]    
           length=lookup[14,i]
            if disklen ne 0L and start ne 0 then begin 
               point_lun,unit,start*bytes    
            endif
            if disklen eq 0L then disklen=length
        
            data1=dblarr(disklen)
            readu,unit,data1
            if swap then data1=swap_endian(data1)
            ;for j=0L,mylen-1 do begin
            for j=0L,length-1 do begin
               jj = (i*length) + j
               data[jj]=double(data1[j])
            endfor        
         endif
      endfor
   endif   
   close,unit
   free_lun,unit

end

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------




















