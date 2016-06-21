load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

ln48 = addfile("/home/cmar/ste69f/umplot/data/landfrac_n48.nc","r")
lv2  = addfile("/home/cmar/ste69f/umplot/data/landfrac_N96.nc","r")
lv1  = addfile("/home/cmar/ste69f/umplot/data/landfrac_ACCESS_N96.nc","r")

lsm1 = ln48->lsm
lsm2 = lv2->lsm
lsm3 = lv1->lsm

lsm_1   = where(lsm1.eq.1,0,lsm1)
lsm_1p1 = where(lsm_1.eq.0,-999,lsm_1)
lsm_1p1@_FillValue = -999

lsm_2   = where(lsm2.eq.1,0,lsm2)
lsm_2p1 = where(lsm_2.eq.0,-999,lsm_2)
lsm_2p1@_FillValue = -999

lsm_3   = where(lsm3.eq.1,0,lsm3)
lsm_3p1 = where(lsm_3.eq.0,-999,lsm_3)
lsm_3p1@_FillValue = -999

;print(lsm1&latitude)
;printVarSummary(lsm1)
;printMinMax(lsm1,True)
;(0)     Land fraction in grid box: min=0   max=1

lats = (/-90,-66,-23,23,66/)
late = (/-66,-23,23,66,90/)

do i =0,4,1
 lsm_4 = lsm1(:,:,{lats(i):late(i)},:)
 ;printVarSummary(lsm_4)
 print(" ")
 print(""+lats(i)+","+late(i)+": "+sum(lsm_4))
 ;print(""+lats(i)+","+late(i)+": "+num(lsm_4.gt.0))
 print(int2flt(num(lsm_4.gt.0))/int2flt(num(lsm1.gt.0)))
 delete(lsm_4)
end do

; N48 ================
print(" ")
print("===== N48 =====")
print("Number of Grid Points     : "+num(lsm1.ge.0))
 ;(0)     7008
print("Number of Land Points     : "+num(lsm1.gt.0))
 ;(0)     2950
print("Number of Full Land Points: "+num(lsm1.eq.1))
 ;(0)     1902
perc1 = ((int2flt(num(lsm1.gt.0))-int2flt(num(lsm1.eq.1)))/int2flt(num(lsm1.gt.0)))*100
print("Number of Coastal Points  : "+(num(lsm1.gt.0)-num(lsm1.eq.1))+" ("+perc1+"%)")
 ;(0)     1048 (35.5254%)
;printMinMax(lsm_1,True)
 ;(0)     min=0   max=0.996825
print(" ")
print("Min & Max Land Fractions (excl. 0 & 1):")
printMinMax(lsm_1p1,True)
 ;(0)     min=0.00295823   max=0.996825
print(" ")
;print(sum(lsm1))
 ;(0)     2361.75
perc_1 = sum(lsm_1)/(int2flt(num(lsm1.gt.0))-int2flt(num(lsm1.eq.1)))
print("Sum of Coastal Fracs: "+sum(lsm_1)+", Avg = "+perc_1+" (Divided by #Cpts)")

; N96 v2 =============
print(" ")
print("===== N96 v2 =====")
print("Number of Grid Points     : "+num(lsm2.ge.0))
 ;(0)     27840
print("Number of Land Points     : "+num(lsm2.gt.0))
 ;(0)     10865
print("Number of Full Land Points: "+num(lsm2.eq.1))
 ;(0)     8431
perc2 = ((int2flt(num(lsm2.gt.0))-int2flt(num(lsm2.eq.1)))/int2flt(num(lsm2.gt.0)))*100
print("Number of Coastal Points  : "+(num(lsm2.gt.0)-num(lsm2.eq.1))+" ("+perc2+"%)")
 ;(0)     2434 (22.4022%)
;printMinMax(lsm_2,True)
 ;(0)     min=0   max=0.996825
print(" ")
print("Min & Max Land Fractions (excl. 0 & 1):")
printMinMax(lsm_2p1,True)
 ;(0)     min=0.00295823   max=0.996825
print(" ")
;print(sum(lsm2))
 ;(0)
perc_2 = sum(lsm_2)/(int2flt(num(lsm2.gt.0))-int2flt(num(lsm2.eq.1)))
print("Sum of Coastal Fracs: "+sum(lsm_2)+", Avg = "+perc_2+" (Divided by #Cpts)")

; N96 v1 =============
print(" ")
print("===== N96 v1 =====")
print("Number of Grid Points     : "+num(lsm3.ge.0))
 ;(0)     27840
print("Number of Land Points     : "+num(lsm3.gt.0))
 ;(0)     11022
print("Number of Full Land Points: "+num(lsm3.eq.1))
 ;(0)     8478
perc3 = ((int2flt(num(lsm3.gt.0))-int2flt(num(lsm3.eq.1)))/int2flt(num(lsm3.gt.0)))*100
print("Number of Coastal Points  : "+(num(lsm3.gt.0)-num(lsm3.eq.1))+" ("+perc3+"%)")
 ;(0)     2544 (23.0811%)
;printMinMax(lsm_3,True)
 ;(0)     min=0   max=0.996825
print(" ")
print("Min & Max Land Fractions (excl. 0 & 1):")
printMinMax(lsm_3p1,True)
 ;(0)     min=0.00295823   max=0.996825
print(" ")
;print(sum(lsm3))
 ;(0)
perc_3 = sum(lsm_3)/(int2flt(num(lsm3.gt.0))-int2flt(num(lsm3.eq.1)))
print("Sum of Coastal Fracs: "+sum(lsm_3)+", Avg = "+perc_3+" (Divided by #Cpts)")

print(" ")

end

