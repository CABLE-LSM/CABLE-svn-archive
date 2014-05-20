; Lauren Stevens 02 Sept 2011
; NCL script to produce Energy and Water Balances.
; Observations: ERA Interim. and see Further down.
; Varaibles are : 
; Precip, P
; Evap, E
; Surface Runoff, R
; Drainage, D
; Screen Temp, T
; Net Rad, Rn
; SW and LW components
; Laten and Sensible Heat, LH & SH
; Water balance  = wbal = P-E-(R+D)
; Energy balance = ebal = Rn-S-L

; ENVIRO VARS --------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
usr = systemfunc("echo $USERID")
loadscript("~"+usr+"/umplot/regrid.ncl")
loadscript("~"+usr+"/umplot/tables.ncl")

;procedure GMcTables()

begin

;usr = systemfunc("echo $USERID")
numYr  = systemfunc("echo $YR")
filDir = systemfunc("echo $DIR")
ModDir = systemfunc("echo $PWD")
;ModDir = systemfunc("echo $DIRW")
resol  = systemfunc("echo $RES")
vMask  = systemfunc("echo $MASK")
;Runid  = systemfunc("echo $RUNID")
;cabdir = systemfunc("echo $CABLE")
;mosdir = systemfunc("echo $MOSES")
tname =  systemfunc("echo $tname")
txname =  systemfunc("echo $txname")
tiname =  systemfunc("echo $tiname")

; FILES --------------------------------------------------------------

;moses_file = addfile(mosdir+"/seasonal_means_"+numYr+"yrs.nc","r") 
;cable_file = addfile(cabdir+"/seasonal_means_"+numYr+"yrs.nc","r")

ObsDir = "~"+usr+"/umplot/obs/"
PltDir = "~"+usr+"/umplot/data/"

sea_file = addfile(ModDir+"/seasonal_means_"+numYr+"yrs.nc","r")
tsea_file = addfile(ModDir+"/Tseasonal_means_"+numYr+"yrs.nc","r")
if (resol.eq.48) then 
lsm_file = addfile(PltDir+"LSM_UMn"+resol+".nc","r")
frac_file = addfile(PltDir+"landfrac_n48.nc","r")
else
if (vMask.eq.2) then
lsm_file = addfile(PltDir+"mask_accessV2_n"+resol+".nc","r")
frac_file = addfile(PltDir+"landfrac_v2.nc","r")
else
lsm_file = addfile(PltDir+"mask_accessV1_n"+resol+".nc","r")
frac_file = addfile(PltDir+"landfrac_ACCESS_N96.nc","r")
end if
end if

area_file = addfile(PltDir+"areacella.nc","r")
;couple = "n" ;systemfunc("echo $CPL")

; SET-UP -------------------------------------------------------------

if (isfilevar(sea_file,"tscrn")) then
 tscrn  = sea_file ->tscrn     ;tas
else
 tscrn  = sea_file ->$tname$   ;tas
end if
if (isfilevar(tsea_file,"tmax")) then
 tmax   = tsea_file ->tmax
 tmin   = tsea_file ->tmin
else
 tmax   = tsea_file ->$txname$ ;temp     ;tasmax
 tmin   = tsea_file ->$tiname$ ;temp_1   ;tasmin
end if
sh     = sea_file ->sh        ;hfss
;rnet1  = sea_file ->field3333 ;field202_1
precip = sea_file ->precip    ;pr
;if (isfilevar(sea_file,"evspsblsoi")) then
if (isfilevar(sea_file,"field1526")) then
 evap1 = sea_file ->field1526 ;evspsblsoi
 evap2 = sea_file ->field1527 ;evspsblveg
end if
if (isfilevar(sea_file,"field1528")) then
 evap3 = sea_file ->field1528 ;evspsblsubl
end if
lh     = sea_file->lh         ;hfls
lsm    = lsm_file  ->lsm
lfrac  = frac_file ->lsm
larea  = area_file ->areacella
;lsm    = frac_file ->lsm
srunoff= sea_file->field1532  ; mrros
drunoff= sea_file->field1533  ; mrro;mrross-is it processed?

rsd1   = sea_file->solar
rld1   = sea_file->longwave
rsds   = sea_file->field203
rlds   = sea_file->ilr
rsus   = (-1)*(rsd1-rsds)
rlus   = (-1)*(rld1-rlds)
copy_VarMeta(rsd1,rsus)
copy_VarMeta(rld1,rlus)

;rsds   = sea_file->rsds
;rsus   = sea_file->rsus
;rlds   = sea_file->rlds
;rlus   = sea_file->rlus

;if (isfilevar(moses_file,"tscrn")) then
;tscrn_m  = moses_file ->tscrn      ;tas
;else
;tscrn_m  = moses_file ->$tname$    ;tas
;end if
;tmax_m   = moses_file ->temp       ;tasmax
;tmin_m   = moses_file ->temp_1     ;tasmin
;sh_m     = moses_file ->sh         ;hfss
;;rnet1_m  = moses_file ->field3333 ;field202_1
;precip_m = moses_file ->precip     ; pr
;;if (isfilevar(moses_file,"evspsblsoi")) then
;if (isfilevar(moses_file,"field1526")) then
; evap1_m = moses_file ->field1526  ;evspsblsoi
; evap2_m = moses_file ->field1527  ;evspsblveg
;end if
;lh_m     = moses_file ->lh         ;hfls
;srunoff_m= moses_file ->field1532  ;mrros
;drunoff_m= moses_file ->field1533  ;mrro;mrross

;rsd1_m  = moses_file ->solar
;rsds_m = moses_file ->field203
;rld1_m  = moses_file ->longwave
;rlds_m = moses_file ->ilr

;rsus_m = (-1)*(rsd1_m-rsds_m)
;rlus_m = (-1)*(rld1_m-rlds_m)

;;rsds_m = moses_file ->rsds
;;rsus_m = moses_file ->rsus
;;rlds_m = moses_file ->rlds
;;rlus_m = moses_file ->rlus

rnet = rsd1 + rld1
;rnet = rsds - rsus + rlds - rlus
;rnet_m = rsd1_m + rld1_m
;;rnet_m = rsds_m - rsus_m + rlds_m - rlus_m  
copy_VarMeta(rsd1,rnet)
;copy_VarMeta(rsd1_m,rnet_m)

trunoff = srunoff ;drunoff
trunoff = trunoff + drunoff
;drunoff = trunoff - srunoff
lhf     = lh

;trunoff_m = srunoff_m ;drunoff_m
;trunoff_m = trunoff_m + drunoff_m
;;drunoff_m = trunoff_m - srunoff_m
;lhf_m     = lh_m

rlam  = 2501000             ;
cls   = 1.1335              ;
rad   = 4.0*atan(1.0)/180.

lat          = precip&latitude
lat!0        = "latitude"
lat&latitude = precip&latitude

precip       = precip*24*60*60
precip@units = "mm/day"
;precip_m       = precip_m*24*60*60
;precip_m@units = "mm/day"

if (isvar("couple")) then
if (couple.eq."y") then
 evap         = (lh/rlam)*24*60*60
 evap@units   = "mm/day"
; evap_m         = (lh_m/rlam)*24*60*60
; evap_m@units   = "mm/day"
else
 evap         = sea_file ->evspsbl
 evap         = (evap)*24*60*60
 ;evap         = (evap1+evap2)*24*60*60
 evap@units   = "mm/day"
; evap_m         = moses_file ->evspsbl
; evap_m         = (evap_m)*24*60*60
; ;evap_m         = (evap1_m+evap2_m)*24*60*60
; evap_m@units   = "mm/day"
end if
else
if (isvar("evap3")) then
 evap         = (evap1+evap2+evap3)*24*60*60
 evap@units   = "mm/day"
else
; evap         = sea_file ->evspsbl
; evap         = (evap)*24*60*60
 evap         = (evap1+evap2)*24*60*60
 evap@units   = "mm/day"
; evap_m         = (evap1_m+evap2_m)*24*60*60
; evap_m@units   = "mm/day"
end if
end if
copy_VarMeta(evap1,evap)

lh            = (lh/rlam)*24*60*60
lh@units      = "mm/day"
lh_snw        = (lh/cls)
copy_VarMeta(lh,lh_snw)
lh_snw@units  = "mm/day"
srunoff       = srunoff*24*60*60
srunoff@units = "mm/day"
drunoff       = drunoff*24*60*60
drunoff@units = "mm/day"
trunoff       = trunoff*24*60*60
trunoff@units = "mm/day"

;lh_m            = (lh_m/rlam)*24*60*60
;lh_m@units      = "mm/day"
;lh_snw_m        = (lh_m/cls)
;copy_VarMeta(lh_m,lh_snw_m)
;lh_snw_m@units  = "mm/day"
;srunoff_m       = srunoff_m*24*60*60
;srunoff_m@units = "mm/day"
;drunoff_m       = drunoff_m*24*60*60
;drunoff_m@units = "mm/day"
;trunoff_m       = trunoff_m*24*60*60
;trunoff_m@units = "mm/day"

prec         = dim_avg_n_Wrap(precip,0)
radnet       = dim_avg_n_Wrap(rnet,0)
sensh        = dim_avg_n_Wrap(sh,0)
lathf        = dim_avg_n_Wrap(lhf,0)
tscr         = dim_avg_n_Wrap(tscrn,0)
tmx          = dim_avg_n_Wrap(tmax,0)
tmn          = dim_avg_n_Wrap(tmin,0)
evp          = dim_avg_n_Wrap(evap,0)
srf          = dim_avg_n_Wrap(srunoff,0)
dr           = dim_avg_n_Wrap(drunoff,0)
trf          = dim_avg_n_Wrap(trunoff,0)
lath         = dim_avg_n_Wrap(lh,0)
lath_snw     = dim_avg_n_Wrap(lh_snw,0)
rsd          = dim_avg_n_Wrap(rsds,0)
rsu          = dim_avg_n_Wrap(rsus,0)
rld          = dim_avg_n_Wrap(rlds,0)
rlu          = dim_avg_n_Wrap(rlus,0)

;prec_m       = dim_avg_n_Wrap(precip_m,0)
;radnet_m     = dim_avg_n_Wrap(rnet_m,0)
;sensh_m      = dim_avg_n_Wrap(sh_m,0)
;lathf_m      = dim_avg_n_Wrap(lhf_m,0)
;tscr_m       = dim_avg_n_Wrap(tscrn_m,0)
;tmx_m        = dim_avg_n_Wrap(tmax_m,0)
;tmn_m        = dim_avg_n_Wrap(tmin_m,0)
;evp_m        = dim_avg_n_Wrap(evap_m,0)
;srf_m        = dim_avg_n_Wrap(srunoff_m,0)
;dr_m         = dim_avg_n_Wrap(drunoff_m,0)
;trf_m        = dim_avg_n_Wrap(trunoff_m,0)
;lath_m       = dim_avg_n_Wrap(lh_m,0)
;lath_snw_m   = dim_avg_n_Wrap(lh_snw_m,0)
;rsd_m        = dim_avg_n_Wrap(rsds_m,0)
;rsu_m        = dim_avg_n_Wrap(rsus_m,0)
;rld_m        = dim_avg_n_Wrap(rlds_m,0)
;rlu_m        = dim_avg_n_Wrap(rlus_m,0)

; balances ----------------------
;wbal         = precip-evap-(srunoff+drunoff)
;wbal         = prec-lath-(srf+dr)
wbal         = prec-evp-(srf+dr)
copy_VarMeta(prec,wbal)
;wbal_m       = prec_m-evp_m-(srf_m+dr_m)
obal         = prec-evp
copy_VarMeta(prec,obal)
;obal_m       = prec_m-evp_m
ebal         = radnet-sensh-lathf 
ebal2        = rsd-rsu+rld-rlu-sensh-lathf 
;ebal_m       = radnet_m-sensh_m-lathf_m 
; -------------------------------

dNames = getvardims(prec)
dims = dimsizes(prec)
nsiz = dimsizes(dims)
nlat = dims(nsiz-2)
mlon = dims(nsiz-1)
;lat  = prec&$dNames(nsiz-2)$
;lon  = prec&$dNames(nsiz-1)$

re    = 6.37122e06
con   = re * rad
clat  = cos(prec&latitude*rad)

dlon = (prec&longitude(2)-prec&longitude(1))
dlat = (prec&latitude(2)-prec&latitude(1))
dx   = con * dlon * clat
dy   = con * dlat
dydx = dy * dx

;wgt   = new((/nlat,mlon/),double)
;wgt   = conform(wgt,dydx,0)
wgt = larea
wgto  = wgt*(1-lfrac(0,0,:,:))
wgtl  = wgt*lfrac(0,0,:,:)
wgtlm = wgt*(lsm(0,0,:,:))
wgtom = wgt*(1-lsm(0,0,:,:))
wgt!0 = "lat"
wgt!1 = "lon"
wgt&lat = precip&latitude
wgt&lon = precip&longitude
wgtl!0 = "lat"
wgtl!1 = "lon"
wgtl&lat = precip&latitude
wgtl&lon = precip&longitude
wgtlm!0 = "lat"
wgtlm!1 = "lon"
wgtlm&lat = precip&latitude
wgtlm&lon = precip&longitude
wgtl!0 = "lat"
wgtl!1 = "lon"
wgtl&lat = precip&latitude
wgtl&lon = precip&longitude
wgtlm!0 = "lat"
wgtlm!1 = "lon"
wgtlm&lat = precip&latitude
wgtlm&lon = precip&longitude
wgto!0 = "lat"
wgto!1 = "lon"
wgto&lat = precip&latitude
wgto&lon = precip&longitude
wgtom!0 = "lat"
wgtom!1 = "lon"
wgtom&lat = precip&latitude
wgtom&lon = precip&longitude

; ANN WBAL -----------------------------------------------------------

; have a external chv ?
mlat      = new((/2,6/),float)
mlat(:,0) = (/-90, 90/)    ; chv=0, Global
mlat(:,1) = (/ 66, 90/)    ; chv=1, High Northern Latitudes (hNlat)
mlat(:,2) = (/ 23, 66/)    ; chv=2, Mid Northern Latitudes  (mNlat)
mlat(:,3) = (/-23, 23/)    ; chv=3, Tropics                 (Trop)
mlat(:,4) = (/-66,-23/)    ; chv=4, Mid Southern Latitudes  (mSlat)
mlat(:,5) = (/-90,-66/)    ; chv=5, High Southern Latitudes (hSlat)
Dlat      = dimsizes(mlat)
wlat = (/"glob","highN","midN","trops","midS","highS"/)

;do k = 0,Dlat-1,1

latS       = mlat(0,chv)
latN       = mlat(1,chv)
latr       = lat({latS:latN})

;region -> Land ---------------------------

pr_l           = wgt_areaave2(prec(:,{latS:latN},:),wgtl({latS:latN},:),0)
;pr_l_m         = wgt_areaave2(prec_m({latS:latN},:),wgtl({latS:latN},:),0)
sh_l           = wgt_areaave2(sensh(:,{latS:latN},:),wgtl({latS:latN},:),0)
;sh_l_m         = wgt_areaave2(sensh_m({latS:latN},:),wgtl({latS:latN},:),0)
rn_l           = wgt_areaave2(radnet(:,{latS:latN},:),wgtl({latS:latN},:),0)
;rn_l_m         = wgt_areaave2(radnet_m({latS:latN},:),wgtl({latS:latN},:),0)

ld_l           = wgt_areaave2(rld(:,{latS:latN},:),wgtl({latS:latN},:),0)
;ld_l_m         = wgt_areaave2(rld_m({latS:latN},:),wgtl({latS:latN},:),0)
lu_l           = wgt_areaave2(rlu(:,{latS:latN},:),wgtl({latS:latN},:),0)
;lu_l_m         = wgt_areaave2(rlu_m({latS:latN},:),wgtl({latS:latN},:),0)
sd_l           = wgt_areaave2(rsd(:,{latS:latN},:),wgtl({latS:latN},:),0)
;sd_l_m         = wgt_areaave2(rsd_m({latS:latN},:),wgtl({latS:latN},:),0)
su_l           = wgt_areaave2(rsu(:,{latS:latN},:),wgtl({latS:latN},:),0)
;su_l_m         = wgt_areaave2(rsu_m({latS:latN},:),wgtl({latS:latN},:),0)

ts_l           = wgt_areaave2(tscr(:,{latS:latN},:),wgtl({latS:latN},:),0)
ts_l           = ts_l-273.16
;ts_l_m         = wgt_areaave2(tscr_m({latS:latN},:),wgtl({latS:latN},:),0)
;ts_l_m         = ts_l_m-273.16
tx_l           = wgt_areaave2(tmx(:,{latS:latN},:),wgtl({latS:latN},:),0)
tx_l           = tx_l-273.16
;tx_l_m         = wgt_areaave2(tmx_m({latS:latN},:),wgtl({latS:latN},:),0)
;tx_l_m         = tx_l_m-273.16
tn_l           = wgt_areaave2(tmn(:,{latS:latN},:),wgtl({latS:latN},:),0)
tn_l           = tn_l-273.16
;tn_l_m         = wgt_areaave2(tmn_m({latS:latN},:),wgtl({latS:latN},:),0)
;tn_l_m         = tn_l_m-273.16

lhf_l          = wgt_areaave2(lathf(:,{latS:latN},:),wgtl({latS:latN},:),0)
;lhf_l_m        = wgt_areaave2(lathf_m({latS:latN},:),wgtl({latS:latN},:),0)
ev_l           = wgt_areaave2(evp(:,{latS:latN},:),wgtl({latS:latN},:),0)
;ev_l           = wgt_areaave2(evp(:,{latS:latN},:),wgtlm({latS:latN},:),0)
;ev_l_m         = wgt_areaave2(evp_m({latS:latN},:),wgtl({latS:latN},:),0)
;;ev_l_m         = wgt_areaave2(evp_m({latS:latN},:),wgtlm({latS:latN},:),0)
lh_l            = wgt_areaave2(lath(:,{latS:latN},:),wgtl({latS:latN},:),0)
lh_l2           = wgt_areaave2(lath_snw(:,{latS:latN},:),wgtl({latS:latN},:),0)
;lh_l_m          = wgt_areaave2(lath_m({latS:latN},:),wgtl({latS:latN},:),0)
;lh_l2_m         = wgt_areaave2(lath_snw_m({latS:latN},:),wgtl({latS:latN},:),0)

rf_l           = wgt_areaave2(srf(:,{latS:latN},:),wgtlm({latS:latN},:),0)
;rf_l_m         = wgt_areaave2(srf_m({latS:latN},:),wgtlm({latS:latN},:),0)
dr_l           = wgt_areaave2(dr(:,{latS:latN},:),wgtlm({latS:latN},:),0)
;dr_l_m         = wgt_areaave2(dr_m({latS:latN},:),wgtlm({latS:latN},:),0)
trf_l          = wgt_areaave2(trf(:,{latS:latN},:),wgtlm({latS:latN},:),0)
;trf_l_m         = wgt_areaave2(trf_m({latS:latN},:),wgtlm({latS:latN},:),0)

tot_l          = wgt_areaave2(wbal(:,{latS:latN},:),wgtl({latS:latN},:),0)
;tot_l          = wgt_areaave2(wbal(:,{latS:latN},:),wgtlm({latS:latN},:),0)
;print((/tot_l,tot_l2/))

resi_l   = pr_l-ev_l-rf_l-dr_l
resi_l2  = pr_l-lh_l-rf_l-dr_l 
resi_l3  = pr_l-lh_l2-rf_l-dr_l

;print("Land: "+(/tot_l,resi_l,resi_l2,resi_l3/))

;resi_lm  = pr_l_m-ev_l_m-rf_l_m-dr_l_m
;resi_lm2 = pr_l_m-lh_l_m-rf_l_m-dr_l_m
;resi_lm3 = pr_l_m-lh_l2_m-rf_l_m-dr_l_m

;region -> Ocean ---------------------------

pr_o           = wgt_areaave2(prec(:,{latS:latN},:),wgto({latS:latN},:),0)
;pr_o_m         = wgt_areaave2(prec_m({latS:latN},:),wgto({latS:latN},:),0)
sh_o           = wgt_areaave2(sensh(:,{latS:latN},:),wgto({latS:latN},:),0)
;sh_o_m         = wgt_areaave2(sensh_m({latS:latN},:),wgto({latS:latN},:),0)
rn_o           = wgt_areaave2(radnet(:,{latS:latN},:),wgto({latS:latN},:),0)
;rn_o_m         = wgt_areaave2(radnet_m({latS:latN},:),wgto({latS:latN},:),0)

ld_o           = wgt_areaave2(rld(:,{latS:latN},:),wgto({latS:latN},:),0)
;ld_o_m         = wgt_areaave2(rld_m({latS:latN},:),wgto({latS:latN},:),0)
lu_o           = wgt_areaave2(rlu(:,{latS:latN},:),wgto({latS:latN},:),0)
;lu_o_m         = wgt_areaave2(rlu_m({latS:latN},:),wgto({latS:latN},:),0)
sd_o           = wgt_areaave2(rsd(:,{latS:latN},:),wgto({latS:latN},:),0)
;sd_o_m         = wgt_areaave2(rsd_m({latS:latN},:),wgto({latS:latN},:),0)
su_o           = wgt_areaave2(rsu(:,{latS:latN},:),wgto({latS:latN},:),0)
;su_o_m         = wgt_areaave2(rsu_m({latS:latN},:),wgto({latS:latN},:),0)

ts_o           = wgt_areaave2(tscr(:,{latS:latN},:),wgto({latS:latN},:),0)
ts_o           = ts_o-273.16
;ts_o_m         = wgt_areaave2(tscr_m({latS:latN},:),wgto({latS:latN},:),0)
;ts_o_m         = ts_o_m-273.16
tx_o           = wgt_areaave2(tmx(:,{latS:latN},:),wgto({latS:latN},:),0)
tx_o           = tx_o-273.16
;tx_o_m         = wgt_areaave2(tmx_m({latS:latN},:),wgto({latS:latN},:),0)
;tx_o_m         = tx_o_m-273.16
tn_o           = wgt_areaave2(tmn(:,{latS:latN},:),wgto({latS:latN},:),0)
tn_o           = tn_o-273.16
;tn_o_m         = wgt_areaave2(tmn_m({latS:latN},:),wgto({latS:latN},:),0)
;tn_o_m         = tn_o_m-273.16

lhf_o          = wgt_areaave2(lathf(:,{latS:latN},:),wgto({latS:latN},:),0)
;lhf_o_m        = wgt_areaave2(lathf_m({latS:latN},:),wgto({latS:latN},:),0)
ev_o           = wgt_areaave2(evp(:,{latS:latN},:),wgto({latS:latN},:),0)
;ev_o           = wgt_areaave2(evp(:,{latS:latN},:),wgtom({latS:latN},:),0)
;ev_o_m         = wgt_areaave2(evp_m({latS:latN},:),wgto({latS:latN},:),0)
;;ev_o_m         = wgt_areaave2(evp_m({latS:latN},:),wgtom({latS:latN},:),0)
lh_o           = wgt_areaave2(lath(:,{latS:latN},:),wgto({latS:latN},:),0)
lh_o2          = wgt_areaave2(lath_snw(:,{latS:latN},:),wgto({latS:latN},:),0)
;lh_o_m         = wgt_areaave2(lath_m({latS:latN},:),wgto({latS:latN},:),0)
;lh_o2_m        = wgt_areaave2(lath_snw_m({latS:latN},:),wgto({latS:latN},:),0)

rf_o           = wgt_areaave2(srf(:,{latS:latN},:),wgtom({latS:latN},:),0)
;rf_o_m         = wgt_areaave2(srf_m({latS:latN},:),wgtom({latS:latN},:),0)
dr_o           = wgt_areaave2(dr(:,{latS:latN},:),wgtom({latS:latN},:),0)
;dr_o_m         = wgt_areaave2(dr_m({latS:latN},:),wgtom({latS:latN},:),0)
trf_o          = wgt_areaave2(trf(:,{latS:latN},:),wgtom({latS:latN},:),0)
;trf_o_m         = wgt_areaave2(trf_m({latS:latN},:),wgtom({latS:latN},:),0)
; -9999 or _FillVal as no values of runoff over ocean

tot_o          = wgt_areaave2(obal(:,{latS:latN},:),wgto({latS:latN},:),0)
;tot_o          = wgt_areaave2(obal(:,{latS:latN},:),wgtom({latS:latN},:),0)

resi_o   = pr_o - ev_o
resi_o2  = pr_o - lh_o
resi_o3  = pr_o - lh_o2

;print("Ocean: "+(/tot_o,resi_o,resi_o2,resi_o3/))

;resi_om  = pr_o_m - ev_o_m
;resi_o2m = pr_o_m - lh_o_m
;resi_o3m = pr_o_m - lh_o2_m

;;region -> Global ---------------------------------------

pr_g   = wgt_areaave2(prec(:,{latS:latN},:),wgt({latS:latN},:),0)
;pr_g_m = wgt_areaave2(prec_m({latS:latN},:),wgt({latS:latN},:),0)
sh_g   = wgt_areaave2(sensh(:,{latS:latN},:),wgt({latS:latN},:),0)
;sh_g_m = wgt_areaave2(sensh_m({latS:latN},:),wgt({latS:latN},:),0)
rn_g   = wgt_areaave2(radnet(:,{latS:latN},:),wgt({latS:latN},:),0)
;rn_g_m = wgt_areaave2(radnet_m({latS:latN},:),wgt({latS:latN},:),0)

ld_g   = wgt_areaave2(rld(:,{latS:latN},:),wgt({latS:latN},:),0)
;ld_g_m = wgt_areaave2(rld_m({latS:latN},:),wgt({latS:latN},:),0)
lu_g   = wgt_areaave2(rlu(:,{latS:latN},:),wgt({latS:latN},:),0)
;lu_g_m = wgt_areaave2(rlu_m({latS:latN},:),wgt({latS:latN},:),0)
sd_g   = wgt_areaave2(rsd(:,{latS:latN},:),wgt({latS:latN},:),0)
;sd_g_m = wgt_areaave2(rsd_m({latS:latN},:),wgt({latS:latN},:),0)
su_g   = wgt_areaave2(rsu(:,{latS:latN},:),wgt({latS:latN},:),0)
;su_g_m = wgt_areaave2(rsu_m({latS:latN},:),wgt({latS:latN},:),0)

ts_g   = wgt_areaave2(tscr(:,{latS:latN},:),wgt({latS:latN},:),0)
ts_g   = ts_g-273.16
;ts_g_m = wgt_areaave2(tscr_m({latS:latN},:),wgt({latS:latN},:),0)
;ts_g_m = ts_g_m-273.16
tx_g   = wgt_areaave2(tmx(:,{latS:latN},:),wgt({latS:latN},:),0)
tx_g   = tx_g-273.16
;tx_g_m = wgt_areaave2(tmx_m({latS:latN},:),wgt({latS:latN},:),0)
;tx_g_m = tx_g_m-273.16
tn_g   = wgt_areaave2(tmn(:,{latS:latN},:),wgt({latS:latN},:),0)
tn_g   = tn_g-273.16
;tn_g_m = wgt_areaave2(tmn_m({latS:latN},:),wgt({latS:latN},:),0)
;tn_g_m = tn_g_m-273.16

lhf_g  = wgt_areaave2(lathf(:,{latS:latN},:),wgt({latS:latN},:),0)
;lhf_g_m= wgt_areaave2(lathf_m({latS:latN},:),wgt({latS:latN},:),0)

;ev_g   = wgt_areaave(evp(:,{latS:latN},:),cos(lat({latS:latN},:)*rad),1.0,0)
;ev_g_m  = wgt_areaave(evp_m({latS:latN},:),cos(lat({latS:latN},:)*rad),1.0,0)
ev_g    = wgt_areaave2(evp(:,{latS:latN},:),wgt({latS:latN},:),0)
;ev_g_m  = wgt_areaave2(evp_m({latS:latN},:),wgt({latS:latN},:),0)

lh_g   = wgt_areaave2(lath(:,{latS:latN},:),wgt({latS:latN},:),0)
lh_g2  = wgt_areaave2(lath_snw(:,{latS:latN},:),wgt({latS:latN},:),0)
;lh_g_m  = wgt_areaave2(lath_m({latS:latN},:),wgt({latS:latN},:),0)
;lh_g2_m = wgt_areaave2(lath_snw_m({latS:latN},:),wgt({latS:latN},:),0)

rf_g   = wgt_areaave2(srf(:,{latS:latN},:),wgt({latS:latN},:),0)
;rf_g_m = wgt_areaave2(srf_m({latS:latN},:),wgt({latS:latN},:),0)
dr_g   = wgt_areaave2(dr(:,{latS:latN},:),wgt({latS:latN},:),0)
;dr_g_m = wgt_areaave2(dr_m({latS:latN},:),wgt({latS:latN},:),0)
trf_g   = wgt_areaave2(trf(:,{latS:latN},:),wgt({latS:latN},:),0)
;trf_g_m = wgt_areaave2(trf_m({latS:latN},:),wgt({latS:latN},:),0)

tot_g  = wgt_areaave(wbal(:,{latS:latN},:),cos(lat({latS:latN})*rad),1.0,0)
;tot_g  = wgt_areaave2(wbal(:,{latS:latN},:),wgt({latS:latN},:),0)
;tot_g_m = wgt_areaave2(wbal({latS:latN},:),wgt({latS:latN},:),0)

resi_g  = pr_g-ev_g-rf_g-dr_g       ; NOTE: -rf_g-dr_g  == -rf_l-dr_l 
resi_g2 = pr_g-lh_g-rf_g-dr_g       ; But not for coupled runs bc of ocean
resi_g3 = pr_g-lh_g2-rf_g-dr_g
;resi_g  = pr_g-ev_g-rf_l-dr_l       ; NOTE: -rf_g-dr_g  == -rf_l-dr_l
;resi_g2 = pr_g-lh_g-rf_l-dr_l
;resi_g3 = pr_g-lh_g2-rf_l-dr_l

;print("Global: "+(/tot_g,resi_g,resi_g2,resi_g3/))

;resi_gm  = pr_g_m-ev_g_m-rf_g_m-dr_g_m       ; NOTE: -rf_g-dr_g  == -rf_l-dr_l
;resi_g2m = pr_g_m-lh_g_m-rf_g_m-dr_g_m      ; But not for coupled runs bc of ocean
;resi_g3m = pr_g_m-lh_g2_m-rf_g_m-dr_g_m

;end do

; Tables -----------------------------------------------------------------------

; Obs. Est. from Budyko 1978 (bu), Henning 1989 (hn), and NOAA (Smith et. al. 2008, no)
not_l = 8.5
not_o = 16.1
not_g = 13.9
bus_l = 63 *0.484583  ; Ly/d --> Wm-2
bus_o = 25 *0.484583
bus_g = 36 *0.484583
bul_l = 74 *0.484583
bul_o = 225*0.484583
bul_g = 181*0.484583
bur_l = 137*0.484583
bur_o = 249*0.484583
bur_g = 216*0.484583
hns_l = 77 *0.484583  ; Ly/d --> Wm-2
hns_o = 29 *0.484583
hns_g = 43 *0.484583
hnl_l = 71 *0.484583
hnl_o = 186*0.484583
hnl_g = 152*0.484583
hnr_l = 149*0.484583
hnr_o = 213*0.484583
hnr_g = 195*0.484583

eb_l = rn_l-sh_l-lhf_l
eb_o = rn_o-sh_o-lhf_o
eb_g = rn_g-sh_g-lhf_g

eb_l2= sd_l-su_l-sh_l-lhf_l-lu_l+ld_l
eb_o2= sd_o-su_o-sh_o-lhf_o-lu_o+ld_o
eb_g2= sd_g-su_g-sh_g-lhf_g-lu_g+ld_g

;print("ebal : "+(/eb_l ,eb_o ,eb_g/))
;print("ebal2: "+(/eb_l2,eb_o2,eb_g2/))

land_vals   = (/sd_l,su_l,sh_l,lhf_l,lu_l,ld_l,rn_l,eb_l/)
ocean_vals  = (/sd_o,su_o,sh_o,lhf_o,lu_o,ld_o,rn_o,eb_o/)
global_vals = (/sd_g,su_g,sh_g,lhf_g,lu_g,ld_g,rn_g,eb_g/)
land_vals@_FillValue   = -999
ocean_vals@_FillValue  = -999
global_vals@_FillValue = -999

ebal_tab(land_vals,ocean_vals,global_vals,numYr,wlat(chv))
if (chv.ge.0 .and. chv.le.5) then
system("mv ebal_table" + "_" + numYr + "yrs.eps ebal_"+wlat(chv)+"_" + numYr + "yrs.eps")
end if

delete(land_vals)
delete(ocean_vals)
delete(global_vals)

;land_vals   = (/-999,ts_l ,not_l,-999,sh_l  ,bus_l,hns_l, \
;                -999,lhf_l,bul_l,hnl_l ,-999,rn_l ,bur_l,hnr_l/)
;ocean_vals  = (/-999,ts_o ,not_o,-999,sh_o  ,bus_o,hns_o, \
;                -999,lhf_o,bul_o,hnl_o ,-999,rn_o ,bur_o,hnr_o/)
;global_vals = (/-999,ts_g ,not_g,-999,sh_g  ,bus_g,hns_g, \
;                -999,lhf_g,bul_g,hnl_g ,-999,rn_g ,bur_g,hnr_g/)
;land_vals@_FillValue   = -999
;ocean_vals@_FillValue  = -999
;global_vals@_FillValue = -999
;
;ebal_compr(land_vals,ocean_vals,global_vals,numYr)
;
;delete(land_vals)
;delete(ocean_vals)
;delete(global_vals)

;=========================
; Obs. Est. from Baumgartner & Reichel 1975 (br) and Legates & Willmott 1990 (lw)
brp_l = 2.03
lwp_l = 2.05
brp_o = 2.92
brp_g = 2.66
bre_l = 1.31
bre_o = 3.22
bre_g = 2.66
brr_l = 0.73

land_vals   = (/pr_l,ev_l,rf_l,dr_l,resi_l,lh_l,resi_l2/);,lh_l2,resi_l3/)
ocean_vals  = (/pr_o,ev_o,rf_o,dr_o,resi_o,lh_o,resi_o2/);,lh_o2,resi_o3/)
global_vals = (/pr_g,ev_g,rf_g,dr_g,resi_g,lh_g,resi_g2/);,lh_g2,resi_g3/)
land_vals@_FillValue   = -999
ocean_vals@_FillValue  = -999
global_vals@_FillValue = -999

wbal_tab(land_vals,ocean_vals,global_vals,numYr,wlat(chv))
if (chv.ge.0 .and. chv.le.5) then
system("mv wBal_table" + "_" + numYr + "yrs.eps wBal_"+wlat(chv)+"_" + numYr + "yrs.eps")
end if


delete(land_vals)
delete(ocean_vals)
delete(global_vals)

;land_vals   = (/-999,pr_l,brp_l,lwp_l,-999,ev_l,bre_l,-999,rf_l,-999,dr_l,brr_l/)
;ocean_vals  = (/-999,pr_o,brp_o,-999 ,-999,lh_o,bre_o,-999,-999,-999,-999,-999/)
;global_vals = (/-999,pr_g,brp_g,-999 ,-999,lh_g,bre_g,-999,-999,-999,-999,-999/)
;land_vals@_FillValue   = -999
;ocean_vals@_FillValue  = -999
;global_vals@_FillValue = -999
;
;wbal_compr(land_vals,ocean_vals,global_vals,numYr)
;
;delete(land_vals)
;delete(ocean_vals)
;delete(global_vals)

;=========================
; Obs. Estimates from NOAA (Smith et. al. 2008 ,no)
;not_l = 8.5
;not_o = 16.1
;not_g = 13.9

land_vals   = (/ts_l,tx_l,tn_l/)
ocean_vals  = (/ts_o,tx_o,tn_o/)
global_vals = (/ts_g,tx_g,tn_g/)
land_vals@_FillValue   = -999
ocean_vals@_FillValue  = -999
global_vals@_FillValue = -999

temp_tab(land_vals,ocean_vals,global_vals,numYr,wlat(chv))
if (chv.ge.0 .and. chv.le.5) then
system("mv temp_table" + "_" + numYr + "yrs.eps temp_"+wlat(chv)+"_" + numYr + "yrs.eps")
end if


delete(land_vals)
delete(ocean_vals)
delete(global_vals)

;land_vals   = (/-999,ts_l,not_l,-999,tx_l,-999,tn_l/)
;ocean_vals  = (/-999,ts_o,not_o,-999,tx_o,-999,tn_o/)
;global_vals = (/-999,ts_g,not_g,-999,tx_g,-999,tn_g/)
;land_vals@_FillValue   = -999
;ocean_vals@_FillValue  = -999
;global_vals@_FillValue = -999
;
;temp_compr(land_vals,ocean_vals,global_vals,numYr)
;
;delete(land_vals)
;delete(ocean_vals)
;delete(global_vals)

;=========================
; Obs. Estimates from Stephens et. al. (s)
ssd_l = -999
ssd_o = -999
ssd_g = 188
ssu_l = -999
ssu_o = -999
ssu_g = 23
ssh_l = -999
ssh_o = -999
ssh_g = 24
slu_l = -999
slu_o = -999
slu_g = 398
slh_l = -999
slh_o = -999
slh_g = 88
sld_l = -999
sld_o = -999
sld_g = 345.6

;land_vals   = (/-999 ,sd_l ,ssd_l,-999 ,su_l ,ssu_l, \
;                -999 ,sh_l ,ssh_l,-999,lhf_l,slh_l, \
;                -999 ,lu_l ,slu_l,-999 ,ld_l ,sld_l/)
;ocean_vals  = (/-999 ,sd_o ,ssd_o,-999 ,su_o ,ssu_o, \
;                -999 ,sh_o ,ssh_o,-999,lhf_o,slh_o, \
;                -999 ,lu_o ,slu_o,-999 ,ld_o ,sld_o/)
;global_vals = (/-999 ,sd_g ,ssd_g,-999 ,su_g ,ssu_g, \
;                -999 ,sh_g ,ssh_g,-999,lhf_g,slh_g, \
;                -999 ,lu_g ,slu_g,-999 ,ld_g ,sld_g/)
;land_vals@_FillValue   = -999
;ocean_vals@_FillValue  = -999
;global_vals@_FillValue = -999
;
;ecomps_compr(land_vals,ocean_vals,global_vals,numYr)

;=========================

end

