load "/home/cmar/ste69f/umplot/wgt_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

fin   = addfile("/home/cmar/ste69f/access/xakfz_5yrs_50kSnow/seasonal_means_5yrs.nc","r")
ain   = addfile("/home/cmar/ste69f/umplot/data/areacella.nc","r")
lin   = addfile("/home/cmar/ste69f/umplot/data/landfrac_v2.nc","r")
prec1 = fin->precip
area  = ain->areacella
lfr   = lin->lsm
larea = area*lfr(0,0,:,:)
prec  = dim_avg_n_Wrap(prec1,0)
prec  = prec*3600*24
lat   = prec&latitude
lat&latitude = prec&latitude

print("")
print(" Global Area: "+sum(area) +" m2")
print("   Land Area: "+sum(larea)+" m2")
print("")

;; Test --------------------------
;
;wgt1  = get_wgt(prec)
;wgtl1 = wgt1*lfr(0,0,:,:)
;pr_g1 = wgt_areaave2(prec,wgt1 ,0)
;pr_l1 = wgt_areaave2(prec,wgtl1,0)
;
;wgt2  = area
;wgtl2 = larea
;pr_g2 = wgt_areaave2(prec,wgt2 ,0)
;pr_l2 = wgt_areaave2(prec,wgtl2,0)
;
;print("Globe w wgt : "+pr_g1+" mm/day")
;print("Globe w area: "+pr_g2+" mm/day")
;print("Globe diff  : "+(pr_g1-pr_g2)+" mm/day")
;print("")
;print(" Land w wgt : "+pr_l1+" mm/day")
;print(" Land w area: "+pr_l2+" mm/day")
;print(" Land diff  : "+(pr_l1-pr_l2)+" mm/day")
;print("")

;; Test --------------------------

; Latitude Bands -----------------

pi_c  = 3.14159265359
re    = 6.37122e06
sarea = 4*pi_c*re*re
print("SurfArea Earth: "+(/sarea/)+" m2")

; tropics
fd     = lat({24})
frad   = (fd*pi_c)/180.
sin_fr = sin(frad)
hgtt   = re*sin(frad)
segmt  = 2*pi_c*re*hgtt
;sfrt   = segmt/sarea ;hgtt/(2*re)
print("  Half Tropics: "+(/segmt/)+" m2")
;print(""+(/sfrt/)+"")

; poles
fp     = lat({66})
prad   = (fp*pi_c)/180.
sin_fp = sin(prad)
hgtp   = re*(1-sin(prad))
segmp  = 2*pi_c*re*hgtp
;sfrp   = segmp/sarea ;hgtp/(2*re)
print("      One Pole: "+(/segmp/)+" m2")
;print(""+(/sfrp/)+"")

;midlats
hgtm = re - hgtt - hgtp ; or re(sin_fp-sin_fr)
segmm= 2*pi_c*re*hgtm
;sfrm = segmm/sarea ;hgtm/(2*re)
print("    One MidLat: "+(/segmm/)+" m2")
;print(""+(/sfrm/)+"")

tota = 2*segmm + 2*segmp + 2*segmt
print("    Total Area: "+(/tota/)+" m2")
print("")

end

