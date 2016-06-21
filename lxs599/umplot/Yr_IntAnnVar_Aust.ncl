; -----------------------------------------------------------------------------
; Lauren Stevens 2 June 2011
; Interannual Variability
; NCL - XY Plots
; -----------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Set-Up ----------------------------------------------------------------------

;procedure IntVar(fileDir:string,OtherDir:string,numYr:string,Mstp:string)
;function  IntVar()

begin

hdir = systemfunc("echo $hdir")
usr = systemfunc("echo $USERID")
mosDir = systemfunc("echo $MOSES")
cabDir = systemfunc("echo $CABLE")
numYr  = systemfunc("echo $YR")
resol  = systemfunc("echo $RES")
lcable = systemfunc("echo $MODEL")
vMask  = systemfunc("echo $MASK")
split  = systemfunc("echo $SPLIT")
tname =  systemfunc("echo $tname")
tstep =  systemfunc("echo $TSTEP")
if (split.eq."y") then
blck  = True
blk   = systemfunc("echo $BLOCK")
block = stringtoint(blk)
else
blck  = False
end if

;dir  = systemfunc("echo $DIR")
;rid  = systemfunc("echo $RUNID")
;bdir = systemfunc("echo $BDIR")
;lmod = systemfunc("echo $MODEL")
;if (lmod.eq."c") then
; cable = True
;else
; cable = False
;end if

ObsDir = hdir+"/"+usr+"/umplot/obs/"
LsmDir = hdir+"/"+usr+"/umplot/data/"

if (isfilepresent(cabDir+"/yearly_means_"+numYr+"yrs.nc") .and. \
    isfilepresent(mosDir+"/yearly_means_"+numYr+"yrs.nc")) then

;system("dmget "+mosDir+"/yearly_means_"+numYr+"yrs.nc")
moses_year = addfile(mosDir+"/yearly_means_"+numYr+"yrs.nc","r")         ; 20 years  ; 1980-2000
;system("dmget "+cabDir+"/yearly_means_"+numYr+"yrs.nc")
cable_year = addfile(cabDir+"/yearly_means_"+numYr+"yrs.nc","r")         ; 20 years  ; 1985-2005
t2mobs_yr  = addfile(ObsDir+"ERAi_yearavg_t2m.nc","r")                   ; 21.3 years; 1989-2010
pr_obs_yr  = addfile(ObsDir+"ERAINT_year_pr8908.nc","r")                 ; 20 years  ; 1989-2008
era_mask   = addfile(ObsDir+"ERAi_LSM_invertlat.nc","r")

if (resol.eq.48) then
lsm_file   = addfile(LsmDir+"LSM_UMn"+resol+".nc","r")
end if
if (resol.eq.96) then
if (vMask.eq.2) then
lsm_file   = addfile(LsmDir+"mask_accessV2_n"+resol+".nc","r")
else
lsm_file   = addfile(LsmDir+"mask_accessV1_n"+resol+".nc","r")
end if
else
lsm_file   = addfile(LsmDir+"mask_access_N"+resol+".nc","r")
end if

; Variables ------------------------------------------------------------------

if (isfilevar(moses_year,"tscrn")) then
 tMyr1           = moses_year ->tscrn
else
 if (isfilevar(moses_year,"tas")) then
  tMyr1           = moses_year ->tas
 else
  tMyr1           = moses_year ->$tname$
 end if
end if
if (isfilevar(cable_year,"tscrn")) then
 tCyr1           = cable_year ->tscrn
else
 if (isfilevar(cable_year,"tas")) then
  tCyr1           = cable_year ->tas
 else
  tCyr1           = cable_year ->$tname$
 end if
end if
tEyr1           = t2mobs_yr  ->t2m
if (isfilevar(moses_year,"field77")) then
pMyr1           = moses_year ->field77
else
if (isfilevar(moses_year,"precip")) then
pMyr1           = moses_year ->precip
else
pMyr1           = moses_year ->pr
end if
end if
if (isfilevar(cable_year,"field77")) then
pCyr1           = cable_year ->field77
else
if (isfilevar(cable_year,"precip")) then
pCyr1           = cable_year ->precip
else
pCyr1           = cable_year ->pr
end if
end if
pEyr1           = pr_obs_yr  ->pr
lsm             = lsm_file   ->lsm
emsk1           = era_mask   ->lsm
if (isfilevar(cable_year,"time")) then
Mtime           = moses_year ->time
Ctime           = cable_year ->time
else
Mtime           = moses_year ->t
Ctime           = cable_year ->t
end if
;Mtim1 = systemfunc("cdo showyear "+mosDir+"/yearly_means_"+numYr+"yrs.nc")
;Ctim1 = systemfunc("cdo showyear "+cabDir+"/yearly_means_"+numYr+"yrs.nc")
;Mtim2 = stringtoint(str_split(Mtim1," "))
;Ctim2 = stringtoint(str_split(Ctim1," "))
;;Mtim2@units = ""
;;Ctim2@units = ""

Myear           = stringtoint(str_get_field(str_get_field(Mtime@units,3," "),1,"-"))
Cyear           = stringtoint(str_get_field(str_get_field(Ctime@units,3," "),1,"-"))
Mmon            = stringtoint(str_get_field(str_get_field(Mtime@units,3," "),2,"-"))
Cmon            = stringtoint(str_get_field(str_get_field(Ctime@units,3," "),2,"-"))

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

if (isfilevar(moses_year,"field77")) then
pMyr1           = pMyr1*tstep
else
pMyr1           = pMyr1*24*60*60
end if
if (isfilevar(cable_year,"field77")) then
pCyr1           = pCyr1*tstep
else
pCyr1           = pCyr1*24*60*60
end if
pEyr1           = pEyr1*24*60*60
pEyr1!1         = "latitude"
pEyr1!2         = "longitude"

emsk            = emsk1*emsk1@scale_factor + emsk1@add_offset
emsk!0          = "time"
emsk!1          = "latitude"
emsk!2          = "longitude"
emsk&latitude   = emsk1&latitude
emsk&longitude  = emsk1&longitude
emsk&time       = emsk1&time

nyrs = stringtoint(numYr)
ttp = dimsizes(tEyr1)
tpr = dimsizes(pEyr1)

if (Mmon.ne.1.and.split.eq."n") then 
 tMyr           = tMyr1(0:nyrs,:,:,:)
 pMyr           = pMyr1(0:nyrs,:,:,:)
else
 tMyr           = tMyr1(0:nyrs-1,:,:,:)
 pMyr           = pMyr1(0:nyrs-1,:,:,:)
end if
if (Cmon.ne.1.and.split.eq."n") then
 tCyr           = tCyr1(0:nyrs,:,:,:)
 pCyr           = pCyr1(0:nyrs,:,:,:)
else
 tCyr           = tCyr1(0:nyrs-1,:,:,:)
 pCyr           = pCyr1(0:nyrs-1,:,:,:)
end if
if (nyrs .gt. ttp(0)) then
tEyr            = tEyr1(0:ttp(0)-1,:,:)
pEyr            = pEyr1(0:tpr(0)-1,:,:)
else
tEyr            = tEyr1(0:nyrs-1,:,:)
pEyr            = pEyr1(0:nyrs-1,:,:)
end if

LandTm          = tMyr 
LandTc          = tCyr
LandTe          = tEyr
LandPm          = pMyr
LandPc          = pCyr
LandPe          = pEyr

; Land-Only: Mask Ocean ------------------------------------------------------

do k=0,nyrs-1,1
LandTm(k:k,:,:,:) = where(lsm.gt.0,tMyr(k:k,:,:,:),-9999)
LandTc(k:k,:,:,:) = where(lsm.gt.0,tCyr(k:k,:,:,:),-9999)
;LandTe(k:k,:,:)   = where(emsk.eq.1,tEyr(k:k,:,:),-9999)
LandPm(k:k,:,:,:) = where(lsm.gt.0,pMyr(k:k,:,:,:),-9999)
LandPc(k:k,:,:,:) = where(lsm.gt.0,pCyr(k:k,:,:,:),-9999)
;LandPe(k:k,:,:)   = where(emsk.eq.1,pEyr(k:k,:,:),-9999)
end do

if (nyrs .gt. ttp(0)) then
 do k=0,ttp(0)-1,1
  LandTe(k:k,:,:)   = where(emsk.eq.1,tEyr(k:k,:,:),-9999)
 end do
 do k=0,tpr(0)-1,1
  LandPe(k:k,:,:)   = where(emsk.eq.1,pEyr(k:k,:,:),-9999)
 end do
else
 do k=0,nyrs-1,1
  LandTe(k:k,:,:)   = where(emsk.eq.1,tEyr(k:k,:,:),-9999)
  LandPe(k:k,:,:)   = where(emsk.eq.1,pEyr(k:k,:,:),-9999)
 end do
end if

LandTm@_FillValue = -9999
LandTc@_FillValue = -9999
LandTe@_FillValue = -9999
LandPm@_FillValue = -9999
LandPc@_FillValue = -9999
LandPe@_FillValue = -9999

; Region --------------------------------------------------------------------

; Area = "Austrlia"
 latS = -45
 latN = -10
 lonL = 110
 lonR = 160

latitude          = tMyr&latitude
latitude!0        = "latitude"
latitude&latitude = tMyr&latitude
lat               = latitude({latS:latN})
lat_era           = tEyr&latitude
lat_era!0         = "latitude" 
lat_era&latitude  = tEyr&latitude
late              = lat_era({latS:latN})

;print(cos(tMyr&latitude*0.0174532))
;print(cos(lat*0.0174532))

; Weighted Averages ---------------------------------------------------------

wgt_tscrnM      = wgt_areaave(LandTm(:,:,{latS:latN},{lonL:lonR}),cos(lat*0.0174532),1.0,0)
;wgt_tscrnM      = wgt_areaave(LandTm,cos(tMyr&latitude*0.0174532),1.0,0)
wgt_tscrnM!0    = "time"
wgt_tscrnM!1    = "ht"
if (isfilevar(moses_year,"time")) then
wgt_tscrnM&time = tMyr&time
else
wgt_tscrnM&time = tMyr&t
end if
wgt_tscrnM&ht   = tMyr&ht
wgt_tscrnC      = wgt_areaave(LandTc(:,:,{latS:latN},{lonL:lonR}),cos(lat*0.0174532),1.0,0)
;wgt_tscrnC      = wgt_areaave(LandTc,cos(tCyr&latitude*0.0174532),1.0,0)
wgt_tscrnC!0    = "time"
wgt_tscrnC!1    = "ht"
if (isfilevar(cable_year,"time")) then
wgt_tscrnC&time = tCyr&time
else
wgt_tscrnC&time = tCyr&t
end if
wgt_tscrnC&ht   = tCyr&ht
wgt_erat2m      = wgt_areaave(LandTe(:,{latS:latN},{lonL:lonR}),cos(late*0.0174532),1.0,0)
;wgt_erat2m      = wgt_areaave(LandTe,cos(tEyr&latitude*0.0174532),1.0,0)
wgt_erat2m!0    = "time"
wgt_erat2m&time = tEyr&time

wgt_prM         = wgt_areaave(LandPm(:,:,{latS:latN},{lonL:lonR}),cos(lat*0.0174532),1.0,0)
;wgt_prM         = wgt_areaave(LandPm,cos(pMyr&latitude*0.0174532),1.0,0)
wgt_prM!0       = "time"
wgt_prM!1       = "ht"
if (isfilevar(moses_year,"time")) then
wgt_prM&time    = pMyr&time
else
wgt_prM&time    = pMyr&t
end if
wgt_prM&ht      = pMyr&surface
wgt_prC         = wgt_areaave(LandPc(:,:,{latS:latN},{lonL:lonR}),cos(lat*0.0174532),1.0,0)
;wgt_prC         = wgt_areaave(LandPc,cos(pCyr&latitude*0.0174532),1.0,0)
wgt_prC!0       = "time"
wgt_prC!1       = "ht"
if (isfilevar(cable_year,"time")) then
wgt_prC&time    = pCyr&time
else
wgt_prC&time    = pCyr&t
end if
wgt_prC&ht      = pCyr&surface
wgt_prE         = wgt_areaave(LandPe(:,{latS:latN},{lonL:lonR}),cos(late*0.0174532),1.0,0)
;wgt_prE         = wgt_areaave(LandPe,cos(pEyr&latitude*0.0174532),1.0,0)
wgt_prE!0       = "time"
wgt_prE&time    = pEyr&time

; Plot ----------------------------------------------------------------------

strtyr = 1978

if (Myear .ge. strtyr) then
Mstart = (Myear - strtyr)
else
Mstart = Myear
end if
if (Mmon.ne.1.and.split.eq."n") then
Mend   = Mstart + nyrs
else
Mend   = Mstart + nyrs -1
end if
if (Cyear .ge. strtyr) then
Cstart = (Cyear - strtyr)
else
Cstart = Cyear
end if
if (Cmon.ne.1.and.split.eq."n") then
Cend   = Cstart + nyrs
else
Cend   = Cstart + nyrs -1
end if
Estart = (1989 - strtyr)
if (nyrs .gt. ttp(0)) then
Eend   = Estart + 20 -1
else
Eend   = Estart + nyrs -1
end if

data                 = new((/3,33/),double) ; data starts at 1978-2010
data(0,:)            = -9999
data(1,:)            = -9999 
data(2,:)            = -9999
data(0,Mstart:Mend)  = wgt_tscrnM(:,0)
data(1,Cstart:Cend)  = wgt_tscrnC(:,0)
if (nyrs .gt. ttp(0)) then
data(2,Estart:Eend)  = wgt_erat2m(0:19)
else
data(2,Estart:Eend)  = wgt_erat2m(0:nyrs-1)
end if
;data&time            = ispan(1978,2010,1)

data2                = new((/3,33/),double)
data2(0,:)           = -9999
data2(1,:)           = -9999
data2(2,:)           = -9999
data2(0,Mstart:Mend) = wgt_prM(:,0)
data2(1,Cstart:Cend) = wgt_prC(:,0)
if (nyrs .gt. tpr(0)) then
data2(2,Estart:Eend) = wgt_prE(0:19) 
else
data2(2,Estart:Eend) = wgt_prE(0:nyrs-1) 
end if
 
yyyymm            = yyyymm_time(1978,2010,"integer")
yyyyfrac          = yyyymm_to_yyyyfrac(yyyymm,0)
if (Cyear .ge. strtyr) then
yyyy              = ispan(1978,2010,1)
else
yyyy              = ispan(Cyear,Cyear+32,1)
end if

wks               = gsn_open_wks("ps","aust_IntvarYear_"+numYr+"yrs")

if (lcable.eq."c") then
Ystart = Cyear
end if
if (lcable.eq."m") then
Ystart = Myear
end if
;Ystart = min((/min((/Myear,1989/)),Cyear/))

res                            = True
res@gsnDraw                    = False
res@gsnFrame                   = False
res@trXMinF                    = Ystart                 ;1985 ;1989
res@trXMaxF                    = Ystart + stringtoint(numYr)         ;2005 ;1998
;res@xyDashPattern             = 0                      ; make all curves solid
res@xyLineColors               = (/"blue","red"/)       ;,"green"/)
res@xyLineThicknesses          = (/2.5,2.5,2.5/)
res@vpHeightF                  = 0.4
res@vpWidthF                   = 0.8
res@tmLabelAutoStride          = True
;res@tmXBLabelConstantSpacingF = 1.0
res@tmXBTickSpacingF           = 5.0
res@tiXAxisFontHeightF         = 0.020
res@tmYLLabelFontHeightF       = 0.015
res@tmXBLabelFontHeightF       = 0.015

print("")
print("Plotting Inter-Annual Variability (Yearly Means) for Australia")
print("")

plot              = new(2,graphic)

res@tiYAxisString = "Temperature (deg K)"
res@tiXAxisString = "Year"
plot(0)           = gsn_csm_xy(wks,yyyy,data,res)
res@tiYAxisString = "Precipitation (mm/day)"
res@tiXAxisString = "Year"
plot(1)           = gsn_csm_xy(wks,yyyy,data2,res)

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
pres@txString             = "Australian Interannual Variability - Yearly"
pres@gsnPanelBottom       = 0.18
;pres@txFontHeightF       = 0.03
;pres@pmLegendDisplayMode = "Always"
gsn_panel(wks,plot,(/2,1/),pres)

else

print("")
print("Yearly Files not Present")
print("")

end if

end

