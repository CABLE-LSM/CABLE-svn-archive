; -----------------------------------------------------------------
; Lauren Stevens 2 June 2011
; Interannual Variability
; NCL - XY Plots
; -----------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
usr    = systemfunc("echo $USERID")
loadscript("~"+usr+"/umplot/plot_bar.ncl")

; Set-Up ----------------------------------------------------------

;procedure BarChart(pltDir:string,:string,Yr:string,pltType:string,Mstp:string)

begin

;usr    = systemfunc("echo $USERID")
;filDir = systemfunc("echo $DIR")
numYr  = systemfunc("echo $YR")
resol  = systemfunc("echo $RES")
mosDir = systemfunc("echo $MOSES")
cabDir = systemfunc("echo $CABLE")
vMask  = systemfunc("echo $MASK")
tname =  systemfunc("echo $tname")

;rid  = systemfunc("echo $RUNID")
;OtherDir = systemfunc("echo $BDIR")
;lmod = systemfunc("echo $MODEL")
;if (lmod.eq."c") then
; cable = True
;else
; cable = False
;end if

ObsDir = "~"+usr+"/umplot/obs/"
LsmDir = "~"+usr+"/umplot/data/"

moses_file = addfile(mosDir+"/seasonal_means_"+numYr+"yrs.nc","r")
cable_file = addfile(cabDir+"/seasonal_means_"+numYr+"yrs.nc","r")

;pr_obs     = addfile(ObsDir+"ERA_INT_pr_8908.nc","r")
;tscrn_obs  = addfile(ObsDir+"ERAi_monavg_t2m.nc","r")
;tscrn_obs  = addfile(ObsDir+"trderi5tp2mmn1.nc","r")
tscrn_obs  = addfile(ObsDir+"ERAi_tp2m.nc","r")
;lsm_file   = addfile(LsmDir+"LSM_UMn"+resol+".nc","r")
;orog_file  = addfile(LsmDir+"orography_UMn"+resol+".nc","r")
era_mask   = addfile(ObsDir+"ERAi_LSM_invertlat.nc","r")

if (resol.eq.48) then
 lsm_file   = addfile(LsmDir+"LSM_UMn48.nc","r")
 orog_file  = addfile(LsmDir+"orography_UMn48.nc","r")
 lf_file   = addfile(LsmDir+"landfrac_n48.nc","r")
else
 if (vMask.eq.2) then
  lsm_file   = addfile(LsmDir+"mask_accessV2_n96.nc","r")
  orog_file  = addfile(LsmDir+"qrparm.orog_hg2.nc","r")
  ;orog_file  = addfile(LsmDir+"qrparm.orog_accessV2.nc","r")
  lf_file   = addfile(LsmDir+"landfrac_v2.nc","r")
 else
  lsm_file   = addfile(LsmDir+"mask_accessV1_n96.nc","r")
  orog_file  = addfile(LsmDir+"qrparm.orog_v1.nc","r")
  lf_file   = addfile(LsmDir+"landfrac_ACCESS_N96.nc","r")
 end if
end if
area_file = addfile(LsmDir+"areacella.nc","r")

; Variables -------------------------------------------------------

;prM1              = moses_file ->field77
;prC1              = cable_file ->field77
;era_pr1           = pr_obs     ->pr
if (isfilevar(moses_file,"tscrn")) then
tscrnM1           = moses_file ->tscrn ; tas
else
if (isfilevar(moses_file,"tas")) then
tscrnM1           = moses_file ->tas
else
tscrnM1           = moses_file ->$tname$ ; tas
end if
end if
if (isfilevar(cable_file,"tscrn")) then
tscrnC1           = cable_file ->tscrn ; tas
else
if (isfilevar(cable_file,"tas")) then
tscrnC1           = cable_file ->tas
else
tscrnC1           = cable_file ->$tname$ ; tas
end if
end if
era_t2m1          = tscrn_obs  ->a ;t2m
lsm               = lsm_file   ->lsm
lfr               = lf_file   ->lsm
;lsm               = lf_file   ->lsm
orog              = orog_file  ->ht
emsk1             = era_mask   ->lsm
larea             = area_file  ->areacella

;Mtime             = moses_file ->t;ime
;Ctime             = cable_file ->t;ime
;Myear             = stringtoint(str_get_field(str_get_field(Mtime@units,3," "),1,"-"))
;Cyear             = stringtoint(str_get_field(str_get_field(Ctime@units,3," "),1,"-"))

emsk              = emsk1*emsk1@scale_factor + emsk1@add_offset
emsk!0            = "time"
emsk!1            = "latitude"
emsk!2            = "longitude"
emsk&latitude     = emsk1&latitude
emsk&longitude    = emsk1&longitude
emsk&time         = emsk1&time
;;emsk_rg           = linint2(emsk&longitude,emsk&latitude,emsk,True,tscrnM1&longitude,tscrnM1&latitude,0)

if (isdim(tscrnM1,"latitude")) then
latitude          = tscrnM1&latitude
latitude!0        = "latitude"
latitude&latitude = tscrnM1&latitude
else
latitude          = tscrnM1&lat
latitude!0        = "latitude"
latitude&latitude = tscrnM1&lat
end if
lat_era           = era_t2m1&latitude
lat_era!0         = "latitude"
lat_era&latitude  = era_t2m1&latitude

; ANN-------------------------------------------------------------

tscrnM          = dim_avg_n_Wrap(tscrnM1,0)
tscrnC          = dim_avg_n_Wrap(tscrnC1,0)
era_t2m         = dim_avg_n_Wrap(era_t2m1,0)
era_t2m         = era_t2m +273.16
era_t2m1        = era_t2m1+273.16
era_t2m1@units  = "deg K"
era_t2m@units   = "deg K"
era_t2m!0       = "latitude"
era_t2m!1       = "longitude"
era_t2m&latitude  = era_t2m1&latitude
era_t2m&longitude = era_t2m1&longitude
;copy_VarMeta(era_t2m,era_t2m)

; Land-Only: Mask Ocean -------------------------------------------

LandTm          = tscrnM1 
LandTc          = tscrnC1
LandTe          = era_t2m1
;OrogTm          = tscrnM
;OrogTc          = tscrnC
;OrogTe          = era_t2m

dimc = dimsizes(dimsizes(tscrnC1))
dimm = dimsizes(dimsizes(tscrnM1))
do i= 0,3,1
if (dimm.eq.4) then
LandTm(i,:,:,:) = where(lsm(0,:,:,:).gt.0,tscrnM1(i,:,:,:),-9999)
else
LandTm(i,:,:) = where(lsm(0,0,:,:).gt.0,tscrnM1(i,:,:),-9999)
end if
LandTc(i,:,:,:) = where(lsm(0,:,:,:).gt.0,tscrnC1(i,:,:,:),-9999)
LandTe(i,:,:)   = where(emsk(0,:,:).eq.1,era_t2m1(i,:,:),-9999)
end do
;LandTm(:,:,:) = where(lsm(0,:,:,:).gt.0,tscrnM(:,:,:),-9999)
;LandTc(:,:,:) = where(lsm(0,:,:,:).gt.0,tscrnC(:,:,:),-9999)
;LandTe(:,:)   = where(emsk(0,:,:).eq.1,era_t2m(:,:),-9999)
LandTm@_FillValue = -9999
LandTc@_FillValue = -9999
LandTe@_FillValue = -9999

;OrogTm(:,:,:) = where(orog(0,:,:,:).lt.1500,tscrnM(:,:,:),-9999)
;OrogTc(:,:,:) = where(orog(0,:,:,:).lt.1500,tscrnC(:,:,:),-9999)
;OrogTe(:,:)   = where((0,:,:).lt.1500,era_t2m(:,:),-9999)          ; *****
;OrogTm@_FillValue = -9999
;OrogTc@_FillValue = -9999
;OrogTe@_FillValue = -9999

LandTmA = dim_avg_n_Wrap(LandTm,0)
LandTcA = dim_avg_n_Wrap(LandTc,0)
LandTeA = dim_avg_n_Wrap(LandTe,0)

; Weighted Averages -----------------------------------------------

mlat      = new((/2,6/),float)
mlat(:,0) = (/-90, 90/)    ; Global
mlat(:,1) = (/ 66, 90/)    ; High Northern Latitudes (hNlat)
mlat(:,2) = (/ 23, 66/)    ; Mid Northern Latitudes  (mNlat)
mlat(:,3) = (/-23, 23/)    ; Tropics                 (Trop)
mlat(:,4) = (/-66,-23/)    ; Mid Southern Latitudes  (mSlat)
mlat(:,5) = (/-90,-66/)    ; High Southern Latitudes (hSlat)
dlat      = dimsizes(mlat)
;print(mlat)

rad = 4.0*atan(1.0) / 180.  ; pi/180. = 0.01745329 ; radi<-->deg

gdata            = new((/21/),double)
gdata            = -9999
gdata@_FillValue = -9999
; gdata #7 will be for b1500m (orog mask), currently set to -9999 (missing)

land_data            = new((/21/),double)
land_data            = -9999
land_data@_FillValue = -9999
; land_data #7 will be for b1500m (orog mask), currently set to -9999 (missing)

re    = 6.37122e06
con   = re * rad
clat  = cos(tscrnC&latitude*rad)

dNames = getvardims(tscrnC)
dims = dimsizes(tscrnC)
nsiz = dimsizes(dims)
nlat = dims(nsiz-2)
mlon = dims(nsiz-1)

dlon = (tscrnC&$dNames(nsiz-1)$(2)-tscrnC&$dNames(nsiz-1)$(1))
dlati = (tscrnC&latitude(2)-tscrnC&latitude(1))
dx   = con * dlon * clat
dy   = con * dlati
dydx = dy * dx

;wgt  = new((/nlat,mlon/),double)
;wgt!0 = "lat"
;wgt!1 = "lon"
;wgt&lat = tscrnC&latitude
;wgt&lon = tscrnC&$dNames(nsiz-1)$ ;longitude
;wgt  = conform(wgt,dydx,0)
wgt   = larea
wgto  = wgt*(1-lfr(0,0,:,:))
wgtl  = wgt*lfr(0,0,:,:)
wgtlm = wgt*(lsm(0,0,:,:)) ;lsm != lfr
wgtom = wgt*(1-lsm(0,0,:,:))
wgtlm!0 = "lat"
wgtlm!1 = "lon"
wgtlm&lat = tscrnC&latitude
wgtlm&lon = tscrnC&$dNames(nsiz-1)$ ;longitude

k    = 0
do i = 0,dlat(1)-1,1
latS       = mlat(0,i)
latN       = mlat(1,i)
lat        = latitude({latS:latN})
late       = lat_era({latS:latN})
if (dimm.eq.4) then
wgt_M      = wgt_areaave2(tscrnM(:,{latS:latN},:),wgt({latS:latN},:),0)
else
wgt_M      = wgt_areaave2(tscrnM({latS:latN},:),wgt({latS:latN},:),0)
end if
;wgt_M      = wgt_areaave(tscrnM(:,{latS:latN},:),cos(lat*rad),1.0,0)
if (dimm.eq.4) then
wgt_M!0    = "ht"
wgt_M@ht   = tscrnM1&ht
end if
wgt_C      = wgt_areaave2(tscrnC(:,{latS:latN},:),wgt({latS:latN},:),0)
;wgt_C      = wgt_areaave(tscrnC(:,{latS:latN},:),cos(lat*rad),1.0,0)
wgt_C!0    = "ht"
wgt_C@ht   = tscrnC1&ht
wgt_E      = wgt_areaave(era_t2m({latS:latN},:),cos(late*rad),1.0,0)
wgt_E!0    = "ht"
wgt_E@ht   = tscrnC1&ht
if (dimm.eq.4) then
wgt_lM     = wgt_areaave2(LandTmA(:,{latS:latN},:),wgtlm({latS:latN},:),0)
else
wgt_lM     = wgt_areaave2(LandTmA({latS:latN},:),wgtlm({latS:latN},:),0)
end if
;wgt_lM     = wgt_areaave(LandTmA(:,{latS:latN},:),cos(lat*rad),1.0,0)
if (dimm.eq.4) then
wgt_lM!0   = "ht"
wgt_lM@ht  = tscrnM1&ht
end if
wgt_lC     = wgt_areaave2(LandTcA(:,{latS:latN},:),wgtlm({latS:latN},:),0)
;wgt_lC     = wgt_areaave(LandTcA(:,{latS:latN},:),cos(lat*rad),1.0,0)
wgt_lC!0   = "ht"
wgt_lC@ht  = tscrnC1&ht
wgt_lE     = wgt_areaave(LandTeA({latS:latN},:),cos(late*rad),1.0,0)
wgt_lE!0   = "ht"
wgt_lE@ht  = tscrnC1&ht
gdata(k)       = wgt_M  -273.16
gdata(k+1)     = wgt_C  -273.16
gdata(k+2)     = wgt_E  -273.16
land_data(k)   = wgt_lM -273.16
land_data(k+1) = wgt_lC -273.16
land_data(k+2) = wgt_lE -273.16
k              = k + 3
delete(latS)
delete(latN)
delete(lat)
delete(late)
delete(wgt_M)
delete(wgt_C)
delete(wgt_E)
delete(wgt_lM)
delete(wgt_lC)
delete(wgt_lE)
end do
;gdata(18) = land_data(0)
;gdata(19) = land_data(1)
;gdata(20) = land_data(2)

land_bias_ann            = new((/14/),double)
land_bias_ann            = -9999
land_bias_ann@_FillValue = -9999
bias_ann                 = new((/14/),double)
bias_ann                 = -9999
bias_ann@_FillValue      = -9999

m    = 0
do j = 0,dimsizes(land_bias_ann)-1,2
land_bias_ann(j)   = land_data(m)  -land_data(m+2)
land_bias_ann(j+1) = land_data(m+1)-land_data(m+2)
bias_ann(j)        = gdata(m)  -gdata(m+2)
bias_ann(j+1)      = gdata(m+1)-gdata(m+2)
m    = m + 3
end do

print("")
print("Plotting Annual Barcharts")
print("")

PlotBar("Global Annual Screen Temperature Means"    ,gdata        ,"bar_tscrn_gann_"+numYr+"yrs"    ,"ps",False)
PlotBar("Global Annual Screen Temperature Biases"   ,bias_ann     ,"bar_tscrn_gbias_"+numYr+"yrs"   ,"ps",False)
PlotBar("Annual Screen Temperature Means for Land" ,land_data    ,"bar_tscrn_land_"+numYr+"yrs"    ,"ps",False)
PlotBar("Annual Screen Temperature Biases for Land",land_bias_ann,"bar_tscrn_landbias_"+numYr+"yrs","ps",False)

; seasonal ========================================================

sdata            = new((/4,21/),double)
sdata            = -9999
sdata@_FillValue = -9999
; gdata #7 will be for b1500m (orog mask), surrently set to -9999 (missing)

land_sdata            = new((/4,21/),double)
land_sdata            = -9999
land_sdata@_FillValue = -9999
; land_sdata #7 will be for b1500m (orog mask), surrently set to -9999 (missing)

do j = 0,3,1          ; four seasons
k    = 0
do i = 0,dlat(1)-1,1
latS       = mlat(0,i)
latN       = mlat(1,i)
lat        = latitude({latS:latN})
late       = lat_era({latS:latN})
if (dimm.eq.4) then
wgt_M      = wgt_areaave2(tscrnM1(j,:,{latS:latN},:),wgt({latS:latN},:),0)
else
wgt_M      = wgt_areaave2(tscrnM1(j,{latS:latN},:),wgt({latS:latN},:),0)
end if
;wgt_M      = wgt_areaave(tscrnM1(j,:,{latS:latN},:),cos(lat*rad),1.0,0)
if (dimm.eq.4) then
wgt_M!0    = "ht"
wgt_M@ht   = tscrnM1&ht
end if
wgt_C      = wgt_areaave2(tscrnC1(j,:,{latS:latN},:),wgt({latS:latN},:),0)
;wgt_C      = wgt_areaave(tscrnC1(j,:,{latS:latN},:),cos(lat*rad),1.0,0)
wgt_C!0    = "ht"
wgt_C@ht   = tscrnC1&ht
wgt_E      = wgt_areaave(era_t2m1(j,{latS:latN},:),cos(late*rad),1.0,0)
wgt_E!0    = "ht"
wgt_E@ht   = tscrnC1&ht
if (dimm.eq.4) then
wgt_lM     = wgt_areaave2(LandTm(j,:,{latS:latN},:),wgtlm({latS:latN},:),0)
else
wgt_lM     = wgt_areaave2(LandTm(j,{latS:latN},:),wgtlm({latS:latN},:),0)
end if
;wgt_lM     = wgt_areaave(LandTm(j,:,{latS:latN},:),cos(lat*rad),1.0,0)
if (dimm.eq.4) then
wgt_lM!0   = "ht"
wgt_lM@ht  = tscrnM1&ht
end if
wgt_lC     = wgt_areaave2(LandTc(j,:,{latS:latN},:),wgtlm({latS:latN},:),0)
;wgt_lC     = wgt_areaave(LandTc(j,:,{latS:latN},:),cos(lat*rad),1.0,0)
wgt_lC!0   = "ht"
wgt_lC@ht  = tscrnC1&ht
wgt_lE     = wgt_areaave(LandTe(j,{latS:latN},:),cos(late*rad),1.0,0)
wgt_lE!0   = "ht"
wgt_lE@ht  = tscrnC1&ht
sdata(j,k)        = wgt_M   -273.16
sdata(j,k+1)      = wgt_C   -273.16
sdata(j,k+2)      = wgt_E   -273.16
land_sdata(j,k)   = wgt_lM  -273.16
land_sdata(j,k+1) = wgt_lC  -273.16
land_sdata(j,k+2) = wgt_lE  -273.16
k    = k + 3
delete(latS)
delete(latN)
delete(lat)
delete(late)
delete(wgt_M)
delete(wgt_C)
delete(wgt_E)
delete(wgt_lM)
delete(wgt_lC)
delete(wgt_lE)
end do
end do

land_bias_sea            = new((/4,14/),double)
land_bias_sea            = -9999
land_bias_sea@_FillValue = -9999
bias_sea                 = new((/4,14/),double)
bias_sea                 = -9999
bias_sea@_FillValue      = -9999
nLen                     = dimsizes(land_bias_sea)

m = 0
do p = 0,nLen(1)-1,2
land_bias_sea(:,p)   = land_sdata(:,m)  -land_sdata(:,m+2)
land_bias_sea(:,p+1) = land_sdata(:,m+1)-land_sdata(:,m+2)
bias_sea(:,p)        = sdata(:,m)  -sdata(:,m+2)
bias_sea(:,p+1)      = sdata(:,m+1)-sdata(:,m+2)
m = m+3
end do

print("Plotting Seasonal BarCharts")
print("")

PlotBar("Global Seasonal Screen Temperature Means"   ,sdata             ,"bar_tscrn_gsea_"+numYr+"yrs"       ,"ps",True)
PlotBar("Global Seasonal Screen Temperature Biases"  ,bias_sea          ,"bar_tscrn_gseabias_"+numYr+"yrs"   ,"ps",True)
PlotBar("Seasonal Screen Temperature Means for Land" ,land_sdata        ,"bar_tscrn_sealand_"+numYr+"yrs"    ,"ps",True)
PlotBar("Seasonal Screen Temperature Biases for Land",land_bias_sea     ,"bar_tscrn_sealandbias_"+numYr+"yrs","ps",True)
;PlotBar("DJF Screen Temperature Biases for Land"     ,land_bias_sea(0,:),"bar_tscrn_djflandbias_"+numYr+"yrs" ,"ps",False)
;PlotBar("MAM Screen Temperature Biases for Land"     ,land_bias_sea(1,:),"bar_tscrn_mamlandbias_"+numYr+"yrs" ,"ps",False)
;PlotBar("JJA Screen Temperature Biases for Land"     ,land_bias_sea(2,:),"bar_tscrn_jjalandbias_"+numYr+"yrs" ,"ps",False)
;PlotBar("SON Screen Temperature Biases for Land"     ,land_bias_sea(3,:),"bar_tscrn_sonlandbias_"+numYr+"yrs" ,"ps",False)

; ================================================================

end

