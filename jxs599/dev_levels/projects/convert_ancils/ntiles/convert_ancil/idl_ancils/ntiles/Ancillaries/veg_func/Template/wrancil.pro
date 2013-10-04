
PRO wrancil_writeit, array, unit, fixhdr, bytes,swap,proto,index,dims
   s=size(array)
   i1=index-1
   i2=index+dims-1
   opoint=fixhdr[i1:i2]

   if s[0] eq 0 then begin
      point=[-32768,replicate(1,dims)]
      changehdr=(opoint[0] ne 0 and opoint[0] ne -32768)
   endif else begin
      if s[0] ne dims then message,'WRANCIL: array dimensions array wrong'
    
      point_lun,-unit,pos
      point=[pos/bytes+1,s[1:dims]]

      size=s[s[0]+2]
      towrite=replicate(proto,size)

      allatonce = 1

      if (allatonce gt 0) then begin
         towrite[*]=array
      endif else begin
         if (s[0] lt 2) then begin
            towrite[*]=array
         endif        
         k=0L 
         if (s[0] eq 2) then begin
            for i=0, s[1]-1 do begin
               for j=0, s[2]-1 do begin
                  towrite[k] = array[i,j]     
                  k=k+1
               endfor        
            endfor        
         endif        
      endelse        
   
     if size gt 10000 then begin
        openw,11,'towrite.txt'
           k=0L 
           for i=0, s[1]-1 do begin
              for j=0, s[2]-1 do begin
                 printf,11, i,j, towrite[k] 
                 k=k+1
              endfor        
           endfor
        close,11
     endif 
   
     if swap then towrite=swap_endian(towrite)
     
      writeu,unit,towrite
    
      changehdr=(max(abs(point-opoint)) ne 0)
   endelse

   if changehdr then begin
      print,'Note: size/pointer info in fixed header changed at element: ',index
      print,'Old: ',opoint
      print,'New: ',point
      fixhdr[i1:i2]=point
   endif

end
;-----------------------------------------------------------------------

PRO wrancil,filename, $
            fixhdr,intc,realc,levdepc,rowdepc,coldepc,fieldc,extrac, $
            temphist,cfi1,cfi2,cfi3,lookup,data, npft, ntimes, $
            swap=swap,bits=bits, $
            makewellformed=makewellformed,sector=sector

if n_elements(makewellformed) eq 0 then makewellformed=0
if n_elements(sector) eq 0 then sector=2048
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

openw,unit,filename,/get_lun

writeu,unit,fixhdr ; this write is just to get file offset right.

wrancil_writeit, intc    ,unit,fixhdr,bytes,swap,int0 ,100,1
wrancil_writeit, realc   ,unit,fixhdr,bytes,swap,real0,105,1
wrancil_writeit, levdepc ,unit,fixhdr,bytes,swap,real0,110,2
wrancil_writeit, rowdepc ,unit,fixhdr,bytes,swap,int0 ,115,2
wrancil_writeit, coldepc ,unit,fixhdr,bytes,swap,int0 ,120,2
wrancil_writeit, fieldc  ,unit,fixhdr,bytes,swap,int0 ,125,2
wrancil_writeit, extrac  ,unit,fixhdr,bytes,swap,int0 ,130,1
wrancil_writeit, temphist,unit,fixhdr,bytes,swap,int0 ,135,1
wrancil_writeit, cfi1    ,unit,fixhdr,bytes,swap,int0 ,140,1
wrancil_writeit, cfi2    ,unit,fixhdr,bytes,swap,int0 ,142,1
wrancil_writeit, cfi3    ,unit,fixhdr,bytes,swap,int0 ,144,1

if makewellformed then point_lun,-unit,startoflookup

point_lun,-unit,jhpos
print, 'jhpos pre write lookup',jhpos
wrancil_writeit, lookup  ,unit,fixhdr,bytes,swap,int0 ,150,2
;  print,'lu14  ', lookup[14,155:170]
point_lun,-unit,jhpos
print, 'jhpos post write lookup',jhpos

print,'jhpos: start data ',lookup[28,0]*bytes

length=lookup[14,0]
;data1=dblarr(length)
index=0L

newnrec = (npft + npft + 1) * ntimes 
openw,13,'disklen.txt'
printf,13, 'start disklen lookup[28] lookup]29] '
for i=0L,(newnrec-1) do begin
   start=lookup[28,i]
   ;print, 'data: start', start
    if (start ne -99) then begin
        disklen=lookup[29,i]
        length=lookup[14,i]

        wellformed=(disklen ne 0 and start ne 0)
        
        point_lun,-unit,pos
        start=pos/bytes
        disklen=ceil(length/float(sector))*sector
        printf,13, i,start, disklen, lookup[28,i], lookup[29,i] 
        
        lookup[29,i]=disklen
        lookup[28,i]=start

        if makewellformed then begin
            point_lun,-unit,pos
            start=pos/bytes
            disklen=ceil(length/float(sector))*sector
            lookup[29,i]=disklen
            lookup[28,i]=start
            print, 'make wellformed',pos
        endif else if wellformed then begin
            point_lun,unit,start*bytes
            print, 'wellformed',start*bytes
        endif else begin
            disklen=length
            print, 'NOT wellformed',disklen
        endelse
         ;point_lun,-unit,jhpos
         ;print, 'jhpos pre write data1',jhpos

        data1=replicate(real0,disklen)
        for j=0L,length-1 do begin
            jj = (i*length) + j
            data1[j]=data[jj]
         endfor        

        if swap then data1=swap_endian(data1)
        writeu,unit,data1    
         ;point_lun,-unit,jhpos
         ;print, 'jhpos post write data1',jhpos
         ;print, ''
    endif
endfor

if makewellformed then begin
    point_lun,unit,startoflookup
    wrancil_writeit, lookup  ,unit,fixhdr,bytes,swap,int0 ,150,2
endif

point_lun,unit,0

fhdr=replicate(int0,256)
fhdr[*]=fixhdr

if swap then fhdr=swap_endian(fhdr)
writeu,unit,fhdr ; now write fixed header including changes
close,unit
close,13

end


