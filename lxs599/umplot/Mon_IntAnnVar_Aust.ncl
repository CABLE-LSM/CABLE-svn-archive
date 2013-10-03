; -----------------------------------------------------------------------------
; Lauren Stevens 2 June 2011
; Interannual Variability
; NCL - XY Plots
; -----------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Set-Up ----------------------------------------------------------------------

;procedure IntVar_Mon(fileDir:string,OtherDir:string,numYr:string,Mstp:string)
;function  IntVar()

begin

cabDir = systemfunc("echo $CABLE")
mosDir = systemfunc("echo $MOSES")
numYr  = systemfunc("echo $YR")
resol  = systemfunc("echo $RES")
lcable = systemfunc("echo $MODEL")
vMask  = systemfunc("echo $MASK")
split  = systemfunc("echo $SPLIT")
tname  =  systemfunc("echo $tname")
if (split.eq."y") then
blck  = True
blk   = systemfunc("echo $BLOCK")
block = stringtoint(blk)
else
blck  = False
end if

;dir  = systemfunc("echo $DIR")
;rid  = systemfunc("echo $RUNID")
;ldir = systemfunc("echo $LDIR")
;lmod = systemfunc("echo $MODEL")
;if (lmod.eq."c") then
; cable = True
;else
; cable = False
;end if

ObsDir = "/home/cmar/ste69f/umplot/obs/"
LsmDir = "/home/cmar/ste69f/umplot/data/"

moses_file = addfile(mosDir+"/Mmonthly_means_"+numYr+"yrs.nc","r")
cable_file = addfile(cabDir+"/Mmonthly_means_"+numYr+"yrs.nc","r")
tscrn_obs  = addfile(ObsDir+"ERAi_monavg_t2m.nc","r")
pr_obs     = addfile(ObsDir+"ERA_INT_pr_8908.nc","r")
era_mask   = addfile(ObsDir+"ERAi_LSM_invertlat.nc","r")

if (resol.eq.48) then
lsm_file   = addfile(LsmDir+"LSM_UMn48.nc","r")
else
if (vMask.eq.2) then
lsm_file   = addfile(LsmDir+"mask_accessV2_n96.nc","r")
else
lsm_file   = addfile(LsmDir+"mask_accessV1_n96.nc","r")
end if
end if

; Variables ------------------------------------------------------------------

if (isfilevar(moses_file,"tscrn")) then
tscrnM1         = moses_file ->tscrn
else
if (isfilevar(moses_file,"tas")) then
tscrnM1         = moses_file ->tas
else
tscrnM1         = moses_file ->$tname$
end if
end if
if (isfilevar(cable_file,"tscrn")) then
tscrnC1         = cable_file ->tscrn
else
if (isfilevar(cable_file,"tas")) then
tscrnC1         = cable_file ->tas
else
tscrnC1         = cable_file ->$tname$
end if
end if
era_t2m1        = tscrn_obs  ->t2m
if (isfilevar(moses_file,"field77")) then
prM1            = moses_file ->field77
else
prM1            = moses_file ->precip
end if
if (isfilevar(cable_file,"field77")) then
prC1            = cable_file ->field77
else
prC1            = cable_file ->precip
end if
era_pr1         = pr_obs     ->pr
lsm             = lsm_file   ->lsm
emsk1           = era_mask   ->lsm
Mtime           = moses_file ->time
Ctime           = cable_file ->time

Myear           = stringtoint(str_get_field(str_get_field(Mtime@units,3," "),1,"-"))
Cyear           = stringtoint(str_get_field(str_get_field(Ctime@units,3," "),1,"-"))
Mmon            = stringtoint(str_get_field(str_get_field(Mtime@units,3," "),2,"-"))
Cmon            = stringtoint(str_get_field(str_get_field(Mtime@units,3," "),2,"-"))

if (blck) then
if (block.eq.2) then
Myear = Myear + 5
Cyear = Cyear + 5
end if 
if (block.eq.3) then
Myear = Myear + 10
Cyear = Cyear + 10
end if
if (block.eq.4) then
Myear = Myear + 15
Cyear = Cyear + 15
end if
end if

if (isfilevar(moses_file,"field77")) then
prM1            = prM1*48
else
prM1            = prM1*24*60*60
end if
if (isfilevar(cable_file,"field77")) then
prC1            = prC1*48
else
prC1            = prC1*24*60*60
end if
era_pr1         = era_pr1*24*60*60
era_pr1!1       = "latitude"
era_pr1!2       = "longitude"

emsk            = emsk1*emsk1@scale_factor + emsk1@add_offset
emsk!0          = "time"
emsk!1          = "latitude"
emsk!2          = "longitude"
emsk&latitude   = emsk1&latitude
emsk&longitude  = emsk1&longitude
emsk&time       = emsk1&time

nmon = stringtoint(numYr)*12
nyrs = stringtoint(numYr)
tpr = dimsizes(era_pr1)
ttp = dimsizes(era_t2m1)

tscrnM          = tscrnM1(0:nmon-1,:,:,:)
tscrnC          = tscrnC1(0:nmon-1,:,:,:)
if ( nmon .gt. ttp(0) ) then
 era_t2m         = era_t2m1(0:ttp(0)-1,:,:)
 era_pr          = era_pr1(0:tpr(0)-1,:,:)
else
 era_t2m         = era_t2m1(0:nmon-1,:,:)
 era_pr          = era_pr1(0:nmon-1,:,:)
end if
prM             = prM1(0:nmon-1,:,:,:)
prC             = prC1(0:nmon-1,:,:,:)

LandTm          = tscrnM 
LandTc          = tscrnC
LandTe          = era_t2m
LandPm          = prM
LandPc          = prC
LandPe          = era_pr

; Land-Only: Mask Ocean ------------------------------------------------------

do k=0,nmon-1,1
;land_only=mask(tscrnM,lsmask,1)
LandTm(k:k,:,:,:) = where(lsm.gt.0,tscrnM(k:k,:,:,:),-9999)
LandTc(k:k,:,:,:) = where(lsm.gt.0,tscrnC(k:k,:,:,:),-9999)
;LandTe(k:k,:,:)   = where(emsk.eq.1,era_t2m(k:k,:,:),-9999)
LandPm(k:k,:,:,:) = where(lsm.gt.0,prM(k:k,:,:,:),-9999)
LandPc(k:k,:,:,:) = where(lsm.gt.0,prC(k:k,:,:,:),-9999)
;LandPe(k:k,:,:)   = where(emsk.eq.1,era_pr(k:k,:,:),-9999)
end do

do k=0,ttp(0)-1,1
LandTe(k:k,:,:)   = where(emsk.eq.1,era_t2m(k:k,:,:),-9999)
end do
do k=0,tpr(0)-1,1
LandPe(k:k,:,:)   = where(emsk.eq.1,era_pr(k:k,:,:),-9999)
end do

LandTm@_FillValue = -9999
LandTc@_FillValue = -9999
LandTe@_FillValue = -9999
LandPm@_FillValue = -9999
LandPc@_FillValue = -9999
LandPe@_FillValue = -9999

; Region --------------------------------------------------------------------

; if (.not.region) then
;   Area = "Global"
;   latS = -90
;   latN =  90
;   lonL =   0
;   lonR = 360
;  else
;   Area = "Australia"
   latS = -45
   latN = -10
   lonL = 110
   lonR = 160
;  end if

;  res@mpMinLatF  = latS
;  res@mpMaxLatF  = latN
;  res@mpMinLonF  = lonL
;  res@mpMaxLonF  = lonR
;  if(region) then
;  res@gsnAddCyclic = False
;  end if

latitude          = tscrnM&latitude
latitude!0        = "latitude"
latitude&latitude = tscrnM&latitude
lat_era           = era_t2m&latitude
lat_era!0         = "latitude"
lat_era&latitude  = era_t2m&latitude

lat               = latitude({latS:latN})
late              = lat_era({latS:latN})

; Weighted Averages ---------------------------------------------------------

wgt_tscrnM      = wgt_areaave(LandTm(:,:,{latS:latN},{lonL:lonR}),cos(lat*0.0174532),1.0,0)
;wgt_tscrnM      = wgt_areaave(LandTm,cos(tscrnM&latitude*0.0174532),1.0,0)
wgt_tscrnM!0    = "time"
wgt_tscrnM!1    = "ht"
wgt_tscrnM&time = tscrnM&time
wgt_tscrnM&ht   = tscrnM&ht
wgt_tscrnC      = wgt_areaave(LandTc(:,:,{latS:latN},{lonL:lonR}),cos(lat*0.0174532),1.0,0)
;wgt_tscrnC      = wgt_areaave(LandTc,cos(tscrnC&latitude*0.0174532),1.0,0)
wgt_tscrnC!0    = "time"
wgt_tscrnC!1    = "ht"
wgt_tscrnC&time = tscrnC&time
wgt_tscrnC&ht   = tscrnC&ht
wgt_erat2m      = wgt_areaave(LandTe(:,{latS:latN},{lonL:lonR}),cos(late*0.0174532),1.0,0)
;wgt_erat2m      = wgt_areaave(LandTe,cos(era_t2m&latitude*0.0174532),1.0,0)
wgt_erat2m!0    = "time"
wgt_erat2m&time = era_t2m&time

wgt_prM         = wgt_areaave(LandPm(:,:,{latS:latN},{lonL:lonR}),cos(lat*0.0174532),1.0,0)
;wgt_prM         = wgt_areaave(LandPm,cos(prM&latitude*0.0174532),1.0,0)
wgt_prM!0       = "time"
wgt_prM!1       = "ht"
wgt_prM&time    = prM&time
wgt_prM&ht      = prM&surface
wgt_prC         = wgt_areaave(LandPc(:,:,{latS:latN},{lonL:lonR}),cos(lat*0.0174532),1.0,0)
;wgt_prC         = wgt_areaave(LandPc,cos(prC&latitude*0.0174532),1.0,0)
wgt_prC!0       = "time"
wgt_prC!1       = "ht"
wgt_prC&time    = prC&time
wgt_prC&ht      = prC&surface
wgt_prE         = wgt_areaave(LandPe(:,{latS:latN},{lonL:lonR}),cos(late*0.0174532),1.0,0)
;wgt_prE         = wgt_areaave(LandPe,cos(era_pr&latitude*0.0174532),1.0,0)
wgt_prE!0       = "time"
wgt_prE&time    = era_pr&time

; Plot ----------------------------------------------------------------------

strtyr = 1978

Mstart = (Myear - strtyr) *12 + (Mmon - 1)
Mend   = Mstart + nyrs*12 -1
Cstart = (Cyear - strtyr) *12 + (Cmon - 1)
Cend   = Cstart + nyrs*12 -1
Estart = (1989 - strtyr) *12
if (nmon .gt. tpr(0)) then
Eend   = Estart + 20*12 -1
else
Eend   = Estart + nyrs*12 -1
end if

data                 = new((/3,396/),double) ; data starts at 1978-2010
data(0,:)            = -9999
data(1,:)            = -9999 
data(2,:)            = -9999
data(0,Mstart:Mend)  = wgt_tscrnM(:,0)
data(1,Cstart:Cend)  = wgt_tscrnC(:,0)
data(2,Estart:Eend)  = wgt_erat2m(0:239)

data2                = new((/3,396/),double)
data2(0,:)           = -9999
data2(1,:)           = -9999
data2(2,:)           = -9999
data2(0,Mstart:Mend) = wgt_prM(:,0)
data2(1,Cstart:Cend) = wgt_prC(:,0)
data2(2,Estart:Eend) = wgt_prE(0:239) 

yyyymm            = yyyymm_time(1978,2010,"integer")
yyyyfrac          = yyyymm_to_yyyyfrac(yyyymm,0)

wks               = gsn_open_wks("ps","aust_IntvarMon_"+numYr+"yrs")

;Ystart = min((/min((/Myear,1989/)),Cyear/))
if (lcable.eq."c") then
Ystart = Cyear
end if
if (lcable.eq."m") then
Ystart = Myear
end if

res                            = True
res@gsnDraw                    = False
res@gsnFrame                   = False
res@trXMinF                    = Ystart                              ;1985 ;1989
res@trXMaxF                    = Ystart + stringtoint(numYr)         ;2005 ;1998
;res@xyDashPattern             = 0                                   ; make all curves solid
res@xyLineColors               = (/"blue","red"/)                    ;,"green"/)
res@xyLineThicknesses          = (/2.5,2.5,2.5/)
res@vpHeightF                  = 0.4
res@vpWidthF                   = 0.8
res@tmLabelAutoStride          = True
;res@tmXBLabelConstantSpacingF = 1.0
res@tmXBTickSpacingF           = 5.0
res@tiXAxisFontHeightF         = 0.020
res@tmYLLabelFontHeightF       = 0.015
res@tmXBLabelFontHeightF       = 0.015

print("Plotting Inter-Annual Variability (Monthly Means) for Australia")
print("")

plot              = new(2,graphic)

res@tiYAxisString = "Temperature (deg K)"
res@tiXAxisString = "Year"
plot(0)           = gsn_csm_xy(wks,yyyyfrac,data,res)
res@tiYAxisString = "Precipitation (mm/day)"
res@tiXAxisString = "Year"
plot(1)           = gsn_csm_xy(wks,yyyyfrac,data2,res)

; Legend ---------------------------------------------------------------------

lgres                    = True
lgres@lgLineColors       = (/"blue","red","black"/)    ; colors for legend lines
;lgres@lgDashIndex       = 0.
lgres@vpWidthF           = 0.2                         ; width of legend (NDC)
lgres@vpHeightF          = 0.1                         ; height of legend (NDC)
lgres@lgPerimOn          = False                       ; turn off perimeter
;lgres@lgDashIndexes     = (/0,12/)                    ; dash indexes
lgres@lgLineLabelStrings = (/"",""/)                   ; no labels in lines

lgd                      = gsn_create_legend(wks,3,(/"Moses","Cable","Obs"/),lgres)
addres                   = True
addres@amParallelPosF    = 0.35
addres@amOrthogonalPosF  = 0.78
add                      = gsn_add_annotation(plot(1),lgd,addres)

; Panel Plot -----------------------------------------------------------------

pres                      = True
pres@gsnMaximize          = True
pres@txString             = "Australian Interannual Variability - Monthly"
pres@gsnPanelBottom       = 0.18
;pres@txFontHeightF       = 0.03
;pres@pmLegendDisplayMode = "Always"
gsn_panel(wks,plot,(/2,1/),pres)

end

