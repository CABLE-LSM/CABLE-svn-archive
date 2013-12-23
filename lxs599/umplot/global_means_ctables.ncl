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
loadscript("~"+usr+"/umplot/wgt_func.ncl")

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
cabdir = systemfunc("echo $CABLE")
mosdir = systemfunc("echo $MOSES")
tname =  systemfunc("echo $tname")
txname =  systemfunc("echo $txname")
tiname =  systemfunc("echo $tiname")

; FILES --------------------------------------------------------------

moses_file = addfile(mosdir+"/seasonal_means_"+numYr+"yrs.nc","r") 
toses_file = addfile(mosdir+"/Tseasonal_means_"+numYr+"yrs.nc","r") 
if (isfilepresent(mosdir+"/clt_yseasmean.nc")) then
coses_file = addfile(mosdir+"/clt_yseasmean.nc","r") 
end if
cable_file = addfile(cabdir+"/seasonal_means_"+numYr+"yrs.nc","r")
table_file = addfile(cabdir+"/Tseasonal_means_"+numYr+"yrs.nc","r")

ObsDir = "~"+usr+"/umplot/obs/"
PltDir = "~"+usr+"/umplot/data/"

erai_file = addfile(ObsDir+"ERAi_tp2m.nc","r")
erat_file = addfile(ObsDir+"ERAi_tmax_tmin.nc","r")
elsm_file = addfile(ObsDir+"ERAi_LSM_invertlat.nc","r")

;sea_file = addfile(ModDir+"/seasonal_means_"+numYr+"yrs.nc","r")
if (resol.eq.48) then 
lsm_file = addfile(PltDir+"LSM_UMn"+resol+".nc","r")
frac_file = addfile(PltDir+"landfrac_n48.nc","r")
else
if (vMask.eq.2) then
 lsm_file = addfile(PltDir+"mask_accessV2_n"+resol+".nc","r")
 frac_file = addfile(PltDir+"landfrac_v2.nc","r")
 area_file = addfile(PltDir+"areacella_fx_ACCESS1-0_historical_r0i0p0.nc","r")
else
lsm_file = addfile(PltDir+"mask_accessV1_n"+resol+".nc","r")
frac_file = addfile(PltDir+"landfrac_ACCESS_N96.nc","r")
end if
end if

;couple = "n" ;systemfunc("echo $CPL")

; SET-UP -------------------------------------------------------------

if (isfilevar(cable_file,"tscrn")) then
tscrn  = cable_file ->tscrn     ;tas
else
tscrn  = cable_file ->$tname$   ;tas
end if
if (isfilevar(table_file,"tmax")) then
 tmax   = table_file ->tmax      ;tasmax
 tmin   = table_file ->tmin      ;tasmin
else
 tmax   = table_file ->$txname$ ;temp      ;tasmax
 tmin   = table_file ->$tiname$ ;temp_1    ;tasmin
end if
etmp = erai_file->a
emax = erat_file->mx2t
emin = erat_file->mn2t

sh     = cable_file ->sh        ;hfss
;rnet1  = cable_file ->field3333 ;field202_1
precip = cable_file ->precip    ;pr
;if (isfilevar(cable_file,"evspsblsoi")) then
if (isfilevar(cable_file,"field1526")) then
 evap1 = cable_file ->field1526 ;evspsblsoi
 evap2 = cable_file ->field1527 ;evspsblveg
end if
if (isfilevar(cable_file,"field1528")) then
 evap3 = cable_file ->field1528 ;evspsblsubl
end if
lh     = cable_file->lh         ;hfls
elsm   = elsm_file ->lsm
lsm    = lsm_file  ->lsm
lfrac  = frac_file ->lsm
;lsm    = frac_file ->lsm
area   = area_file->areacella
srunoff= cable_file->field1532  ; mrros
drunoff= cable_file->field1533  ; mrro;mrross-is it processed?
clt    = cable_file->field30_1
clt    = clt*100

rsd1   = cable_file->solar
rld1   = cable_file->longwave
rsds   = cable_file->field203
rlds   = cable_file->ilr
rsus   = (-1)*(rsd1-rsds)
rlus   = (-1)*(rld1-rlds)

;rsds   = cable_file->rsds
;rsus   = cable_file->rsus
;rlds   = cable_file->rlds
;rlus   = cable_file->rlus

if (isfilevar(moses_file,"tscrn")) then
tscrn_m  = moses_file ->tscrn      ;tas
else
if (isfilevar(moses_file,"tas")) then
tscrn_m  = moses_file ->tas
else
tscrn_m  = moses_file ->$tname$
end if
end if
if (isfilevar(moses_file,"tasmax")) then
tmax_m   = moses_file ->tasmax
else
 if (isfilevar(toses_file,"tmax")) then
  tmax_m   = toses_file ->tmax
 else
  tmax_m   = toses_file ->$txname$ ;temp
 end if
end if
if (isfilevar(moses_file,"tasmin")) then
tmin_m   = moses_file ->tasmin
else
 if (isfilevar(toses_file,"tmin")) then
  tmin_m   = toses_file ->tmin
 else
  tmin_m   = toses_file ->$tiname$ ;temp_1
 end if
end if
if (isfilevar(moses_file,"hfss")) then
sh_m     = moses_file ->hfss
else
sh_m     = moses_file ->sh         ;hfss
end if
;rnet1_m  = moses_file ->field3333 ;field202_1
if (isfilevar(moses_file,"pr")) then
precip_m = moses_file ->pr
else
precip_m = moses_file ->precip     ; pr
end if
if (isfilevar(moses_file,"evspsblsoi")) then
 evap1_m = moses_file ->evspsblsoi
 evap2_m = moses_file ->evspsblveg
end if
if (isfilevar(moses_file,"field1526")) then
 evap1_m = moses_file ->field1526  ;evspsblsoi
 evap2_m = moses_file ->field1527  ;evspsblveg
end if
if (isfilevar(moses_file,"field1528")) then
 evap3_m = moses_file ->field1528  ;evspsblsubl
end if
if (isfilevar(moses_file,"hfls")) then
lh_m     = moses_file ->hfls
else
lh_m     = moses_file ->lh         ;hfls
end if
diml = dimsizes(lh_m)
rank = dimsizes(getvardims(lh_m))
if (isfilevar(moses_file,"mrro")) then
 srunoff_m = moses_file ->mrros
 drunoff_m = moses_file ->mrro ;mrross
 trunoff_m = drunoff_m
 drunoff_m = trunoff_m - srunoff_m
else
 if (isfilevar(moses_file,"field1532")) then
  srunoff_m = moses_file ->field1532  ;mrros
  drunoff_m = moses_file ->field1533  ;mrro;mrross
  trunoff_m = srunoff_m
  trunoff_m = trunoff_m + drunoff_m
 else
  if (rank.eq.4) then
   srunoff_m = new((/diml(0),diml(1),diml(2),diml(3)/),float)
   drunoff_m = new((/diml(0),diml(1),diml(2),diml(3)/),float)
   trunoff_m = new((/diml(0),diml(1),diml(2),diml(3)/),float)
  else
   srunoff_m = new((/diml(0),diml(1),diml(2)/),float)
   drunoff_m = new((/diml(0),diml(1),diml(2)/),float)
   trunoff_m = new((/diml(0),diml(1),diml(2)/),float)
  end if
  srunoff_m = 0.
  drunoff_m = 0.
  trunoff_m = 0.
 end if
end if
if (isfilepresent(mosdir+"/clt_yseasmean.nc")) then
if (isfilevar(coses_file,"clt")) then
clt_m = coses_file->clt
end if
else
if (isfilevar(moses_file,"field30_1")) then
clt_m = moses_file->field30_1
clt_m = clt_m*100
else
if (rank.eq.4) then
clt_m = new((/diml(0),diml(1),diml(2),diml(3)/),float)
else
clt_m = new((/diml(0),diml(1),diml(2)/),float)
end if
clt_m = 0.
end if
end if

if (isfilevar(moses_file,"rsds")) then
 rsds_m = moses_file ->rsds
 rsus_m = moses_file ->rsus
 rlds_m = moses_file ->rlds
 rlus_m = moses_file ->rlus
 rnet_m = rsds_m - rsus_m + rlds_m - rlus_m
else
 rsd1_m = moses_file ->solar
 rsds_m = moses_file ->field203
 rld1_m = moses_file ->longwave
 rlds_m = moses_file ->ilr
 rsus_m = (-1)*(rsd1_m-rsds_m)
 rlus_m = (-1)*(rld1_m-rlds_m)
 rnet_m = rsd1_m + rld1_m
end if

rnet   = rsd1 + rld1
;rnet   = rsds - rsus + rlds - rlus

trunoff = srunoff ;drunoff
trunoff = trunoff + drunoff
;drunoff = trunoff - srunoff

lhf     = lh
lhf_m     = lh_m

rlam  = 2501000             ;
cls   = 1.1335              ;
rad   = 4.0*atan(1.0)/180.

lat          = precip&latitude
lat!0        = "latitude"
lat&latitude = precip&latitude

precip       = precip*24*60*60
precip@units = "mm/day"
precip_m       = precip_m*24*60*60
precip_m@units = "mm/day"

if (isvar("couple")) then
if (couple.eq."y") then
 evap         = (lh/rlam)*24*60*60
 evap@units   = "mm/day"
 evap_m         = (lh_m/rlam)*24*60*60
 evap_m@units   = "mm/day"
else
 evap         = cable_file ->evspsbl
 evap         = (evap)*24*60*60
 ;evap         = (evap1+evap2)*24*60*60
 evap@units   = "mm/day"
 evap_m         = moses_file ->evspsbl
 evap_m         = (evap_m)*24*60*60
 ;evap_m         = (evap1_m+evap2_m)*24*60*60
 evap_m@units   = "mm/day"
end if
else
if (isvar("evap3")) then
 evap         = (evap1+evap2+evap3)*24*60*60
 evap@units   = "mm/day"
 evap_m         = (evap1_m+evap2_m+evap3_m)*24*60*60
 evap_m@units   = "mm/day"
else
; evap         = cable_file ->evspsbl
; evap         = (evap)*24*60*60
 evap         = (evap1+evap2)*24*60*60
 evap@units   = "mm/day"
 evap_m         = (evap1_m+evap2_m)*24*60*60
 evap_m@units   = "mm/day"
end if
end if

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

lh_m            = (lh_m/rlam)*24*60*60
lh_m@units      = "mm/day"
lh_snw_m        = (lh_m/cls)
copy_VarMeta(lh_m,lh_snw_m)
lh_snw_m@units  = "mm/day"
srunoff_m       = srunoff_m*24*60*60
srunoff_m@units = "mm/day"
drunoff_m       = drunoff_m*24*60*60
drunoff_m@units = "mm/day"
trunoff_m       = trunoff_m*24*60*60
trunoff_m@units = "mm/day"

prec         = dim_avg_n_Wrap(precip,0)
radnet       = dim_avg_n_Wrap(rnet,0)
sensh        = dim_avg_n_Wrap(sh,0)
lathf        = dim_avg_n_Wrap(lhf,0)
tscr         = dim_avg_n_Wrap(tscrn,0)
tmx          = dim_avg_n_Wrap(tmax,0)
tmn          = dim_avg_n_Wrap(tmin,0)
etm          = dim_avg_n_Wrap(etmp,0)
emx          = dim_avg_n_Wrap(emax,0)
emn          = dim_avg_n_Wrap(emin,0)
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
cld          = dim_avg_n_Wrap(clt,0)

prec_m       = dim_avg_n_Wrap(precip_m,0)
radnet_m     = dim_avg_n_Wrap(rnet_m,0)
sensh_m      = dim_avg_n_Wrap(sh_m,0)
lathf_m      = dim_avg_n_Wrap(lhf_m,0)
tscr_m       = dim_avg_n_Wrap(tscrn_m,0)
tmx_m        = dim_avg_n_Wrap(tmax_m,0)
tmn_m        = dim_avg_n_Wrap(tmin_m,0)
evp_m        = dim_avg_n_Wrap(evap_m,0)
srf_m        = dim_avg_n_Wrap(srunoff_m,0)
dr_m         = dim_avg_n_Wrap(drunoff_m,0)
trf_m        = dim_avg_n_Wrap(trunoff_m,0)
lath_m       = dim_avg_n_Wrap(lh_m,0)
lath_snw_m   = dim_avg_n_Wrap(lh_snw_m,0)
rsd_m        = dim_avg_n_Wrap(rsds_m,0)
rsu_m        = dim_avg_n_Wrap(rsus_m,0)
rld_m        = dim_avg_n_Wrap(rlds_m,0)
rlu_m        = dim_avg_n_Wrap(rlus_m,0)
cld_m        = dim_avg_n_Wrap(clt_m,0)

;wbal         = precip-evap-(srunoff+drunoff)
wbal         = prec-evp-(srf+dr)
wbal_m       = prec_m-evp_m-(srf_m+dr_m)
obal         = prec-evp
obal_m       = prec_m-evp_m
ebal         = radnet-sensh-lathf 
ebal_m       = radnet_m-sensh_m-lathf_m 

dNames = getvardims(prec)
dims   = dimsizes(prec)
nsiz   = dimsizes(dims)
nlat   = dims(nsiz-2)
mlon   = dims(nsiz-1)
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

wgt   = new((/nlat,mlon/),double)
wgt   = conform(wgt,dydx,0)
;wgt = area
wgto  = wgt*(1-lfrac(0,0,:,:))
wgtl  = wgt*lfrac(0,0,:,:)
wgtlm = wgt*(lsm(0,0,:,:))
wgtom = wgt*(1-lsm(0,0,:,:))

wgte  = get_wgt(emx)
;wgtle = wgte*elsm(0,:,:)
;wgtlo = wgte*(1-elsm(0,:,:))
wgtle= where(elsm(0,:,:).gt.0,wgte*elsm(0,:,:),-9999)
wgtlo= where(elsm(0,:,:).le.0,wgte*(1 - elsm(0,:,:)),-9999)


; ANN WBAL -----------------------------------------------------------
;region -> Land ---------------------------

pr_l           = wgt_areaave2(prec    ,wgtl,0)
sh_l           = wgt_areaave2(sensh   ,wgtl,0)
rn_l           = wgt_areaave2(radnet  ,wgtl,0)
pr_l_m         = wgt_areaave2(prec_m  ,wgtl,0)
sh_l_m         = wgt_areaave2(sensh_m ,wgtl,0)
rn_l_m         = wgt_areaave2(radnet_m,wgtl,0)

ld_l           = wgt_areaave2(rld   ,wgtl,0)
lu_l           = wgt_areaave2(rlu   ,wgtl,0)
sd_l           = wgt_areaave2(rsd   ,wgtl,0)
su_l           = wgt_areaave2(rsu   ,wgtl,0)
cl_l           = wgt_areaave2(cld   ,wgtl,0)
ld_l_m         = wgt_areaave2(rld_m ,wgtl,0)
lu_l_m         = wgt_areaave2(rlu_m ,wgtl,0)
sd_l_m         = wgt_areaave2(rsd_m ,wgtl,0)
su_l_m         = wgt_areaave2(rsu_m ,wgtl,0)
cl_l_m         = wgt_areaave2(cld_m ,wgtl,0)

ts_l           = wgt_areaave2(tscr  ,wgtl,0)
ts_l           = ts_l-273.16
tx_l           = wgt_areaave2(tmx   ,wgtl,0)
tx_l           = tx_l-273.16
tn_l           = wgt_areaave2(tmn   ,wgtl,0)
tn_l           = tn_l-273.16
ts_l_m         = wgt_areaave2(tscr_m,wgtl,0)
ts_l_m         = ts_l_m-273.16
tx_l_m         = wgt_areaave2(tmx_m ,wgtl,0)
tx_l_m         = tx_l_m-273.16
tn_l_m         = wgt_areaave2(tmn_m ,wgtl,0)
tn_l_m         = tn_l_m-273.16

et_l           = wgt_areaave2(etm   ,wgtle,0)
;et_l           = et_l-273.16
ex_l           = wgt_areaave2(emx   ,wgtle,0)
ex_l           = ex_l-273.16
en_l           = wgt_areaave2(emn   ,wgtle,0)
en_l           = en_l-273.16

;if () then
;
;end if
lhf_l          = wgt_areaave2(lathf     ,wgtl,0)
ev_l           = wgt_areaave2(evp       ,wgtl,0) ;wgtlm
lh_l           = wgt_areaave2(lath      ,wgtl,0)
lh_l2          = wgt_areaave2(lath_snw  ,wgtl,0)
lhf_l_m        = wgt_areaave2(lathf_m   ,wgtl,0)
ev_l_m         = wgt_areaave2(evp_m     ,wgtl,0) ;wgtlm
lh_l_m         = wgt_areaave2(lath_m    ,wgtl,0)
lh_l2_m        = wgt_areaave2(lath_snw_m,wgtl,0)

rf_l           = wgt_areaave2(srf  ,wgtl,0) ;wgtlm
dr_l           = wgt_areaave2(dr   ,wgtl,0) ;wgtlm
trf_l          = wgt_areaave2(trf  ,wgtl,0) ;wgtlm
rf_l_m         = wgt_areaave2(srf_m,wgtl,0) ;wgtlm
dr_l_m         = wgt_areaave2(dr_m ,wgtl,0) ;wgtlm
trf_l_m        = wgt_areaave2(trf_m,wgtl,0) ;wgtlm

;Mtot            = where(lsm(0,0,:,:).gt.0,wbal,-9999)
Mtot            = where(lsm(:,0,:,:).gt.0,wbal,-9999)
Mtot@_FillValue = -9999
tot_l           = wgt_areaave(Mtot,cos(lat*rad),1.0,0)
;tot_l2          = wgt_areaave2(wbal,wgtl,0)
;tot_l3          = wgt_areaave2(wbal,wgtlm,0)

resi_l  = pr_l-ev_l-rf_l-dr_l
resi_l2 = pr_l-lh_l-rf_l-dr_l 
resi_l3 = pr_l-lh_l2-rf_l-dr_l

resi_lm  = pr_l_m-ev_l_m-rf_l_m-dr_l_m
resi_lm2 = pr_l_m-lh_l_m-rf_l_m-dr_l_m
resi_lm3 = pr_l_m-lh_l2_m-rf_l_m-dr_l_m

;region -> Ocean ---------------------------

pr_o           = wgt_areaave2(prec    ,wgto,0)
sh_o           = wgt_areaave2(sensh   ,wgto,0)
rn_o           = wgt_areaave2(radnet  ,wgto,0)
pr_o_m         = wgt_areaave2(prec_m  ,wgto,0)
sh_o_m         = wgt_areaave2(sensh_m ,wgto,0)
rn_o_m         = wgt_areaave2(radnet_m,wgto,0)

ld_o           = wgt_areaave2(rld  ,wgto,0)
lu_o           = wgt_areaave2(rlu  ,wgto,0)
sd_o           = wgt_areaave2(rsd  ,wgto,0)
su_o           = wgt_areaave2(rsu  ,wgto,0)
cl_o           = wgt_areaave2(cld  ,wgto,0)
ld_o_m         = wgt_areaave2(rld_m,wgto,0)
lu_o_m         = wgt_areaave2(rlu_m,wgto,0)
sd_o_m         = wgt_areaave2(rsd_m,wgto,0)
su_o_m         = wgt_areaave2(rsu_m,wgto,0)
cl_o_m         = wgt_areaave2(cld_m,wgto,0)

ts_o           = wgt_areaave2(tscr  ,wgto,0)
ts_o           = ts_o-273.16
tx_o           = wgt_areaave2(tmx   ,wgto,0)
tx_o           = tx_o-273.16
tn_o           = wgt_areaave2(tmn   ,wgto,0)
tn_o           = tn_o-273.16
ts_o_m         = wgt_areaave2(tscr_m,wgto,0)
ts_o_m         = ts_o_m-273.16
tx_o_m         = wgt_areaave2(tmx_m ,wgto,0)
tx_o_m         = tx_o_m-273.16
tn_o_m         = wgt_areaave2(tmn_m ,wgto,0)
tn_o_m         = tn_o_m-273.16

et_o           = wgt_areaave2(etm,wgtlo,0)
;et_o           = et_o-273.16
ex_o           = wgt_areaave2(emx,wgtlo,0)
ex_o           = ex_o-273.16
en_o           = wgt_areaave2(emn,wgtlo,0)
en_o           = en_o-273.16

;if () then
;
;end if
lhf_o          = wgt_areaave2(lathf     ,wgto,0)
ev_o           = wgt_areaave2(evp       ,wgto,0) ; wgtom
lh_o           = wgt_areaave2(lath      ,wgto,0)
lh_o2          = wgt_areaave2(lath_snw  ,wgto,0)
lhf_o_m        = wgt_areaave2(lathf_m   ,wgto,0)
ev_o_m         = wgt_areaave2(evp_m     ,wgto,0) ; wgtom
lh_o_m         = wgt_areaave2(lath_m    ,wgto,0)
lh_o2_m        = wgt_areaave2(lath_snw_m,wgto,0)

rf_o           = wgt_areaave2(srf  ,wgto,0) ; wgtom
dr_o           = wgt_areaave2(dr   ,wgto,0) ; wgtom
trf_o          = wgt_areaave2(trf  ,wgto,0) ;wgtom
rf_o_m         = wgt_areaave2(srf_m,wgto,0) ; wgtom
dr_o_m         = wgt_areaave2(dr_m ,wgto,0) ; wgtom
trf_o_m        = wgt_areaave2(trf_m,wgto,0) ;wgtom

;;OtRf           = where(lsm(0,0,:,:).eq.0,trf,-9999)
;OtRf           = where(lsm(:,0,:,:).eq.0,trf,-9999)
;OtRf@_FillValue= -9999
;trf_o          = wgt_areaave(OtRf,cos(lat*rad),1.0,0)

;;Otot           = where(lsm(0,0,:,:).eq.0,obal,-9999)
;Otot           = where(lsm(:,0,:,:).eq.0,obal,-9999)
;Otot@_FillValue= -9999
;tot_o          = wgt_areaave(Otot,cos(lat*rad),1.0,0)
tot_o          = wgt_areaave2(obal,wgto,0)
tot_o_m        = wgt_areaave2(obal_m,wgto,0)

resi_o   = pr_o - ev_o
resi_o2  = pr_o - lh_o
resi_o3  = pr_o - lh_o2

resi_om  = pr_o_m - ev_o_m
resi_o2m = pr_o_m - lh_o_m
resi_o3m = pr_o_m - lh_o2_m

;;region -> Global ---------------------------------------

pr_g        = wgt_areaave2(prec    ,wgt,0)
sh_g        = wgt_areaave2(sensh   ,wgt,0)
rn_g        = wgt_areaave2(radnet  ,wgt,0)
pr_g_m      = wgt_areaave2(prec_m  ,wgt,0)
sh_g_m      = wgt_areaave2(sensh_m ,wgt,0)
rn_g_m      = wgt_areaave2(radnet_m,wgt,0)

ld_g        = wgt_areaave2(rld  ,wgt,0)
lu_g        = wgt_areaave2(rlu  ,wgt,0)
sd_g        = wgt_areaave2(rsd  ,wgt,0)
su_g        = wgt_areaave2(rsu  ,wgt,0)
cl_g        = wgt_areaave2(cld  ,wgt,0)
ld_g_m      = wgt_areaave2(rld_m,wgt,0)
lu_g_m      = wgt_areaave2(rlu_m,wgt,0)
sd_g_m      = wgt_areaave2(rsd_m,wgt,0)
su_g_m      = wgt_areaave2(rsu_m,wgt,0)
cl_g_m      = wgt_areaave2(cld_m,wgt,0)

ts_g        = wgt_areaave2(tscr  ,wgt,0)
ts_g        = ts_g-273.16
tx_g        = wgt_areaave2(tmx   ,wgt,0)
tx_g        = tx_g-273.16
tn_g        = wgt_areaave2(tmn   ,wgt,0)
tn_g        = tn_g-273.16
ts_g_m      = wgt_areaave2(tscr_m,wgt,0)
ts_g_m      = ts_g_m-273.16
tx_g_m      = wgt_areaave2(tmx_m ,wgt,0)
tx_g_m      = tx_g_m-273.16
tn_g_m      = wgt_areaave2(tmn_m ,wgt,0)
tn_g_m      = tn_g_m-273.16

et_g        = wgt_areaave2(etm, wgte,0)
;et_g        = et_g-273.16
ex_g        = wgt_areaave2(emx, wgte,0)
ex_g        = ex_g-273.16
en_g        = wgt_areaave2(emn, wgte,0)
en_g        = en_g-273.16

;if () then
;
;end if
lhf_g       = wgt_areaave2(lathf     ,wgt,0)
ev_g        = wgt_areaave( evp       ,cos(lat*rad),1.0,0)
lh_g        = wgt_areaave2(lath      ,wgt,0)
lh_g2       = wgt_areaave2(lath_snw  ,wgt,0)
lhf_g_m     = wgt_areaave2(lathf_m   ,wgt,0)
ev_g_m      = wgt_areaave( evp_m     ,cos(lat*rad),1.0,0)
lh_g_m      = wgt_areaave2(lath_m    ,wgt,0)
lh_g2_m     = wgt_areaave2(lath_snw_m,wgt,0)

rf_g        = wgt_areaave( srf  ,cos(lat*rad),1.0,0)
dr_g        = wgt_areaave( dr   ,cos(lat*rad),1.0,0)
trf_g       = wgt_areaave( trf  ,cos(lat*rad),1.0,0)
rf_g_m      = wgt_areaave2(srf_m,wgt,0)
dr_g_m      = wgt_areaave2(dr_m ,wgt,0)
trf_g_m     = wgt_areaave( trf_m,cos(lat*rad),1.0,0)

tot_g  = wgt_areaave(wbal  ,cos(lat*rad),1.0,0)
tot_g_m= wgt_areaave(wbal_m,cos(lat*rad),1.0,0)

resi_g  = pr_g-ev_g-rf_g-dr_g       ; NOTE: -rf_g-dr_g  == -rf_l-dr_l 
resi_g2 = pr_g-lh_g-rf_g-dr_g       ; But not for coupled runs bc of ocean
resi_g3 = pr_g-lh_g2-rf_g-dr_g
;resi_g  = pr_g-ev_g-rf_l-dr_l       ; NOTE: -rf_g-dr_g  == -rf_l-dr_l
;resi_g2 = pr_g-lh_g-rf_l-dr_l
;resi_g3 = pr_g-lh_g2-rf_l-dr_l

resi_gm  = pr_g_m-ev_g_m-rf_g_m-dr_g_m       ; NOTE: -rf_g-dr_g  == -rf_l-dr_l
resi_g2m = pr_g_m-lh_g_m-rf_g_m-dr_g_m      ; But not for coupled runs bc of ocean
resi_g3m = pr_g_m-lh_g2_m-rf_g_m-dr_g_m

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

land_vals   = (/ts_l_m ,ts_l ,not_l,sh_l_m,sh_l  ,bus_l,hns_l, \
                lhf_l_m,lhf_l,bul_l,hnl_l ,rn_l_m,rn_l ,bur_l,hnr_l/)
ocean_vals  = (/ts_o_m ,ts_o ,not_o,sh_o_m,sh_o  ,bus_o,hns_o, \
                lhf_o_m,lhf_o,bul_o,hnl_o ,rn_o_m,rn_o ,bur_o,hnr_o/) 
global_vals = (/ts_g_m ,ts_g ,not_g,sh_g_m,sh_g  ,bus_g,hns_g, \
                lhf_g_m,lhf_g,bul_g,hnl_g ,rn_g_m,rn_g ,bur_g,hnr_g/) 
land_vals@_FillValue   = -999
ocean_vals@_FillValue  = -999
global_vals@_FillValue = -999

;print((/land_vals/))
ebal_compr(land_vals,ocean_vals,global_vals,numYr,"glob","1")

delete(land_vals)
delete(ocean_vals)
delete(global_vals)

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

land_vals   = (/pr_l_m,pr_l,brp_l ,lwp_l,lh_l_m,lh_l,bre_l, \
                rf_l_m,rf_l,dr_l_m,dr_l ,brr_l/)
ocean_vals  = (/pr_o_m,pr_o,brp_o , -999,lh_o_m,lh_o,bre_o, \
                -999  ,-999, -999 , -999, -999/)
global_vals = (/pr_g_m,pr_g,brp_g , -999,lh_g_m,lh_g,bre_g, \
                -999  ,-999, -999 , -999, -999/)
land_vals@_FillValue   = -999
ocean_vals@_FillValue  = -999
global_vals@_FillValue = -999

;print((/land_vals/))
wbal_compr(land_vals,ocean_vals,global_vals,numYr,"glob","1")

delete(land_vals)
delete(ocean_vals)
delete(global_vals)

;=========================
; Obs. Estimates from NOAA (Smith et. al. 2008 ,no)
; from above
;not_l = 8.5
;not_o = 16.1
;not_g = 13.9

land_vals   = (/ts_l_m,ts_l,not_l,et_l,tx_l_m,tx_l,ex_l,tn_l_m,tn_l,en_l/)
ocean_vals  = (/ts_o_m,ts_o,not_o,et_o,tx_o_m,tx_o,ex_o,tn_o_m,tn_o,en_o/)
global_vals = (/ts_g_m,ts_g,not_g,et_g,tx_g_m,tx_g,ex_g,tn_g_m,tn_g,en_g/)
land_vals@_FillValue   = -999
ocean_vals@_FillValue  = -999
global_vals@_FillValue = -999

;print((/land_vals/))
temp_compr(land_vals,ocean_vals,global_vals,numYr,"glob","1")

delete(land_vals)
delete(ocean_vals)
delete(global_vals)

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

land_vals   = (/sd_l_m ,sd_l ,ssd_l,su_l_m ,su_l ,ssu_l, \
                sh_l_m ,sh_l ,ssh_l,lhf_l_m,lhf_l,slh_l, \
                lu_l_m ,lu_l ,slu_l,ld_l_m ,ld_l ,sld_l/)
ocean_vals  = (/sd_o_m ,sd_o ,ssd_o,su_o_m ,su_o ,ssu_o, \
                sh_o_m ,sh_o ,ssh_o,lhf_o_m,lhf_o,slh_o, \
                lu_o_m ,lu_o ,slu_o,ld_o_m ,ld_o ,sld_o/)
global_vals = (/sd_g_m ,sd_g ,ssd_g,su_g_m ,su_g ,ssu_g, \
                sh_g_m ,sh_g ,ssh_g,lhf_g_m,lhf_g,slh_g, \
                lu_g_m ,lu_g ,slu_g,ld_g_m ,ld_g ,sld_g/)
land_vals@_FillValue   = -999
ocean_vals@_FillValue  = -999
global_vals@_FillValue = -999

;print((/land_vals/))
ecomps_compr(land_vals,ocean_vals,global_vals,numYr,"glob","1")

;=========================

end

