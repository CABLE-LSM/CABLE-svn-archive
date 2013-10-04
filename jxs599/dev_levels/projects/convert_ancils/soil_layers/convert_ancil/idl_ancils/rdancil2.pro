FUNCTION rdancil2_readit,unit,fixhdr,bytes,swap,proto,index,dims
if fixhdr[index-1] le 0 then return,proto*0
size=1L
for i=1,dims do size=size*fixhdr[index-1+i]
answer=reform(replicate(proto,size),fixhdr[index:index+dims-1])
point_lun,unit,(fixhdr[index-1]-1)*bytes
readu,unit,answer
if swap then answer=swap_endian(answer)
return,answer
end


PRO rdancil2,filename, $
             fixhdr,intc,realc,levdepc,rowdepc,coldepc,fieldc,extrac, $
             temphist,cfi1,cfi2,cfi3,lookup,data, $
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

intc    =rdancil2_readit(unit,fixhdr,bytes,swap,int0 ,100,1)
realc   =rdancil2_readit(unit,fixhdr,bytes,swap,real0,105,1)
levdepc =rdancil2_readit(unit,fixhdr,bytes,swap,real0,110,2)
rowdepc =rdancil2_readit(unit,fixhdr,bytes,swap,int0 ,115,2)
coldepc =rdancil2_readit(unit,fixhdr,bytes,swap,int0 ,120,2)
fieldc  =rdancil2_readit(unit,fixhdr,bytes,swap,int0 ,125,2)
extrac  =rdancil2_readit(unit,fixhdr,bytes,swap,int0 ,130,1)
temphist=rdancil2_readit(unit,fixhdr,bytes,swap,int0 ,135,1)
cfi1    =rdancil2_readit(unit,fixhdr,bytes,swap,int0 ,140,1)
cfi2    =rdancil2_readit(unit,fixhdr,bytes,swap,int0 ,142,1)
cfi3    =rdancil2_readit(unit,fixhdr,bytes,swap,int0 ,144,1)
lookup  =rdancil2_readit(unit,fixhdr,bytes,swap,int0 ,150,2)

nrec=(size(lookup))[2]
data=ptrarr(nrec)

for i=0L,nrec-1L do begin
    start=lookup[28,i]
    if (start ne -99) then begin
        disklen=lookup[29,i]    
        length=lookup[14,i]
        pack=(lookup[20,i] mod 10 ne 0)
        floats=((not pack) and abs(lookup[38,i]) eq 1)
        
        if pack then length=disklen
        if disklen ne 0 and start ne 0 $
          then point_lun,unit,start*bytes    
        if disklen eq 0 then disklen=length

        if floats then proto=real0 else proto=int0
        data1=replicate(proto,disklen)
        readu,unit,data1
        data[i]=ptr_new(data1[0:length-1])
    endif
endfor
if swap then data=swap_endian(data)

close,unit
free_lun,unit
end
