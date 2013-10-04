FUNCTION rdancil_readit,unit,fixhdr,bytes,swap,proto,index,dims
if fixhdr[index-1] le 0 then return,proto*0
size=1L
for i=1,dims do size=size*fixhdr[index-1+i]
answer=reform(replicate(proto,size),fixhdr[index:index+dims-1])
point_lun,unit,(fixhdr[index-1]-1)*bytes
readu,unit,answer
if swap then answer=swap_endian(answer)
return,answer
end

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------

PRO rdancil,filename, $
            fixhdr,intc,realc,levdepc,rowdepc,coldepc,fieldc,extrac, $
            temphist,cfi1,cfi2,cfi3,lookup,data, $
            swap=swap,bits=bits

if n_elements(bits) eq 0 then bits=64
if n_elements(swap) eq 0 then swap=0

case bits of
    64: begin
        real0=0D0
        fl0=0.
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
endfor
;jhan{
;mydlen = 4096L
;mylen = 7008L
;ndata = mylen * 9 
length=lookup[14,0]
data=dblarr(ndata)
data1=fltarr(length)
;}

index=0L
start=2048L

for i=0,nrec-1 do begin
  start=lookup[28,i]
;jhan{
     ;disklen = mydlen
     ;length = mylen
     disklen=lookup[29,i]    
     length=lookup[14,i]
    ;print,'rdancil:i, lookup[29/14/28,i]', i+1, disklen, length, start
;}
    if (start ne -99L) then begin
        if disklen ne 0L and start ne 0 $
          then point_lun,unit,start*bytes    
        if disklen eq 0L then disklen=length
        readu,unit,data1
        if swap then data1=swap_endian(data1)
       ;for j=0L,mylen-1 do begin
       for j=0L,length-1 do begin
            jj = (i*length) + j
            data[jj]=double(data1[j])
         endfor        
        index=index+length
    endif
endfor

close,unit
free_lun,unit

end

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------




















