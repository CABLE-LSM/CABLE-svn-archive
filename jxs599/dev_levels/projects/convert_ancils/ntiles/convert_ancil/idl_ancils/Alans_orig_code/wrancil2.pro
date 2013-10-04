PRO wrancil2_writeit,array,unit,fixhdr,bytes,swap,proto,index,dims

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
    towrite[*]=array
    
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

PRO wrancil2,filename, $
             fixhdr,intc,realc,levdepc,rowdepc,coldepc,fieldc,extrac, $
             temphist,cfi1,cfi2,cfi3,lookup,data, $
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

wrancil2_writeit, intc    ,unit,fixhdr,bytes,swap,int0 ,100,1
wrancil2_writeit, realc   ,unit,fixhdr,bytes,swap,real0,105,1
wrancil2_writeit, levdepc ,unit,fixhdr,bytes,swap,int0 ,110,2
wrancil2_writeit, rowdepc ,unit,fixhdr,bytes,swap,int0 ,115,2
wrancil2_writeit, coldepc ,unit,fixhdr,bytes,swap,int0 ,120,2
wrancil2_writeit, fieldc  ,unit,fixhdr,bytes,swap,int0 ,125,2
wrancil2_writeit, extrac  ,unit,fixhdr,bytes,swap,int0 ,130,1
wrancil2_writeit, temphist,unit,fixhdr,bytes,swap,int0 ,135,1
wrancil2_writeit, cfi1    ,unit,fixhdr,bytes,swap,int0 ,140,1
wrancil2_writeit, cfi2    ,unit,fixhdr,bytes,swap,int0 ,142,1
wrancil2_writeit, cfi3    ,unit,fixhdr,bytes,swap,int0 ,144,1

if makewellformed then point_lun,-unit,startoflookup
wrancil2_writeit, lookup  ,unit,fixhdr,bytes,swap,int0 ,150,2

; now point to start of data
point_lun,-unit,pos

if makewellformed then begin
    pos=ceil(pos/float(sector))*sector    
    point_lun,unit,pos
endif
fixhdr[159]=1+pos/bytes

totdisklen=0L
address=1L
for i=0L,(size(lookup))[2]-1 do begin
    start=lookup[28,i]
    if (start ne -99) then begin
        disklen=lookup[29,i]
        length=lookup[14,i]
        pack=(lookup[20,i] mod 10 ne 0)
        if pack then length=disklen
        
        wellformed=(disklen ne 0 and start ne 0)
        
        if makewellformed then begin
            point_lun,-unit,pos
            start=pos/bytes
            disklen=ceil(length/float(sector))*sector
            lookup[29,i]=disklen
            lookup[28,i]=start
            lookup[39,i]=address
        endif else if wellformed then begin
            point_lun,unit,start*bytes
        endif else begin
            disklen=length
        endelse
        
        if ((not pack) and abs(lookup[38,i]) eq 1) $
          then proto=real0 else proto=int0
        ;; proto=(*data[i])[0]*0
        
        totdisklen=totdisklen+disklen
        address=address+length
        
        data1=replicate(proto,disklen)
        data1[0:length-1]=*data[i]
        if swap then data1=swap_endian(data1)

        writeu,unit,data1    
    endif
endfor

if makewellformed then begin
    point_lun,unit,startoflookup
    wrancil2_writeit, lookup  ,unit,fixhdr,bytes,swap,int0 ,150,2
endif

point_lun,unit,0

fhdr=replicate(int0,256)
fhdr[*]=fixhdr
fhdr[160]=totdisklen

if swap then fhdr=swap_endian(fhdr)
writeu,unit,fhdr ; now write fixed header including changes
close,unit
end
