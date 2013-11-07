; Lauren Stevens 02 Sept 2011
; NCL script to produce Water Balances.
; Observations: ERA Interim. 
; Varaibles are : 
; Precip P
; Evap E
; Surface Runoff R
; Drainage D
; Water balance = wbal = P-E-R-D

; ENVIRO VARS -----------------------------------------------------


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;procedure WTable()

begin

numYr  = systemfunc("echo $YR")
filDir = systemfunc("echo $DIR")
ModDir = systemfunc("echo $PWD")
;ModDir = systemfunc("echo $DIRW")
resol  = systemfunc("echo $RES")
vMask  = systemfunc("echo $MASK")
;Runid  = systemfunc("echo $RUNID")
;cable  = systemfunc("echo $CABLE")
;moses  = systemfunc("echo $MOSES")

; FILES -----------------------------------------------------------
;moses_file = addfile(moses+"/season_`YR`yrs.nc","r") 
;cable_file = addfile(cable+"/season_`YR`yrs.nc","r")

ObsDir = "/home/cmar/ste69f/umplot/obs/"
DatDir = "/home/cmar/ste69f/umplot/data/"

sea_file = addfile(ModDir+"/seasonal_means_"+numYr+"yrs.nc","r")
if (isfilepresent(ModDir+"/Rseasonal_means_"+numYr+"yrs.nc")) then
 rsea_file = addfile(ModDir+"/Rseasonal_means_"+numYr+"yrs.nc","r")
end if
;era_file = addfile(ObsDir+"trderi5tp2mmn1.nc","r")
if (resol.eq.48) then 
 lsm_file = addfile(DatDir+"LSM_UMn"+resol+".nc","r")
 ;oro_file = addfile(DatDir+"orography_UMn"+resol+".nc","r")
 frac_file = addfile(DatDir+"landfrac_n48.nc","r")
else
 if (vMask.eq.2) then
  lsm_file = addfile(DatDir+"mask_accessV2_n"+resol+".nc","r")
  frac_file = addfile(DatDir+"landfrac_v2.nc","r")
  area_file = addfile(DatDir+"areacella_fx_ACCESS1-0_historical_r0i0p0.nc","r")
 else
  lsm_file = addfile(DatDir+"mask_accessV1_n"+resol+".nc","r")
  frac_file = addfile(DatDir+"landfrac_ACCESS_N96.nc","r")
 end if
end if

;couple = systemfunc("echo $CPL")

; SET-UP ----------------------------------------------------------

if (isfilevar(sea_file,"sh")) then
 sh      = sea_file ->sh ;hfss
else
 sh      = sea_file ->hfss
end if
rnet1   = sea_file ->field202_1
precip  = sea_file ->precip
if (isfilevar(sea_file,"field1526")) then
 evap1   = sea_file ->field1526
 evap2   = sea_file ->field1527
end if
if (isfilevar(sea_file,"field1528")) then
 evap3   = sea_file ->field1528
end if
lh      = sea_file ->lh ;hfls
if (isfilevar(sea_file,"field184")) then
 mflx   = sea_file->field184
end if
lsm     = lsm_file ->lsm
lfrac   = frac_file->lsm
;lsm     = frac_file->lsm
area    = area_file->areacella
if(isfilepresent(ModDir+"/Rseasonal_means_"+numYr+"yrs.nc")) then
 srunoff = rsea_file ->field1532
 drunoff = rsea_file ->field1533
else
 srunoff = sea_file ->field1532
 drunoff = sea_file ->field1533
end if

rsd = sea_file->solar
;rsu = sea_file->field203
rld = sea_file->longwave
;rlu = sea_file->ilr

;rsds = rsu
;rsus = (-1)*(rsd-rsu)
;rlds = rlu
;rlus = (-1)*(rld-rlu)

rnet = rsd+rld

trunoff = srunoff
trunoff = trunoff + drunoff
lhf     = lh

if (isvar("mflx")) then
 mflx = mflx*24*60*60
 mflx@units = "mm/day"
end if

rlam  = 2501000             ;
cls   = 1.1335              ;
rad   = 4.0*atan(1.0)/180.
re    = 6.37122e06
con   = re * rad

lat          = precip&latitude
lat!0        = "latitude"
lat&latitude = precip&latitude

clat  = cos( lat * rad )

precip       = precip*24*60*60
precip@units = "mm/day"

if (isvar("couple")) then
 if (couple.eq."y") then
  evap         = (lh/rlam)*24*60*60
  evap@units   = "mm/day"
 else
  if (isvar("evap3")) then
   evap         = (evap1+evap2+evap3)*24*60*60
   evap@units   = "mm/day"
  else
   evap         = (evap1+evap2)*24*60*60
   evap@units   = "mm/day"
  end if
 end if
else
if (isvar("evap3")) then
 evap         = (evap1+evap2)*24*60*60
 ;evap         = (evap1+evap2+evap3)*24*60*60
 evap@units   = "mm/day"
else
 evap         = (evap1+evap2)*24*60*60
 evap@units   = "mm/day"
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

radnet       = dim_avg_n_Wrap(rnet,0)
sensh        = dim_avg_n_Wrap(sh,0)
lathf        = dim_avg_n_Wrap(lhf,0)
prec         = dim_avg_n_Wrap(precip,0)
evp          = dim_avg_n_Wrap(evap,0)
srf          = dim_avg_n_Wrap(srunoff,0)
dr           = dim_avg_n_Wrap(drunoff,0)
trf          = dim_avg_n_Wrap(trunoff,0)
lath         = dim_avg_n_Wrap(lh,0)
lath_snw     = dim_avg_n_Wrap(lh_snw,0)

;wbal         = precip-evap-(srunoff+drunoff)
wbal         = prec-evp-(srf+dr) ; +riverout ?
obal         = prec-evp ; +riverout

dName= getvardims(prec)
dims = dimsizes(prec)
nsiz = dimsizes(dims)
nlat = dims(nsiz-2)
mlon = dims(nsiz-1)

dlon = (prec&$dName(nsiz-1)$(2)-prec&$dName(nsiz-1)$(1))
dlat = (lat(2)-lat(1))
dx   = con * dlon * clat
dy   = con * dlat
dydx = dy * dx

wgt  = new((/nlat,mlon/),typeof(prec))
wgt  = conform(wgt,dydx,0)
;wgt = area
wgto = wgt*(1-lfrac(0,0,:,:))
wgtl = wgt*lfrac(0,0,:,:)
;wgtom= wgt*(1-lsm(0,0,:,:))
;wgtlm= wgt*lsm(0,0,:,:)

; ANN WBAL --------------------------------------------------------
;region -> Global (Land Only) ---------------------------

pr_l    = wgt_areaave2(prec    ,wgtl,0)
sh_l    = wgt_areaave2(sensh   ,wgtl,0)
rn_l    = wgt_areaave2(radnet  ,wgtl,0)
lhf_l   = wgt_areaave2(lathf   ,wgtl,0)
ev_l    = wgt_areaave2(evp     ,wgtl,0) ;wgtlm
lh_l    = wgt_areaave2(lath    ,wgtl,0)
lh_l2   = wgt_areaave2(lath_snw,wgtl,0)
rf_l    = wgt_areaave2(srf     ,wgtl,0) ;wgtlm
dr_l    = wgt_areaave2(dr      ,wgtl,0) ;wgtlm
trf_l   = wgt_areaave2(trf     ,wgtl,0) ;wgtlm
tot_l   = wgt_areaave2(wbal    ,wgtl,0)

resi_l  = pr_l-ev_l-rf_l-dr_l
resi_l2 = pr_l-lh_l-rf_l-dr_l 
resi_l3 = pr_l-lh_l2-rf_l-dr_l

;region -> Global (Ocean Only) ---------------------------

pr_o    = wgt_areaave2(prec    ,wgto,0)
sh_o    = wgt_areaave2(sensh   ,wgto,0)
rn_o    = wgt_areaave2(radnet  ,wgto,0)
lhf_o   = wgt_areaave2(lathf   ,wgto,0)
ev_o    = wgt_areaave2(evp     ,wgto,0) ;wgtom
lh_o    = wgt_areaave2(lath    ,wgto,0)
lh_o2   = wgt_areaave2(lath_snw,wgto,0)
rf_o    = wgt_areaave2(srf     ,wgto,0) ;wgtom
dr_o    = wgt_areaave2(dr      ,wgto,0) ;wgtom
trf_o   = wgt_areaave2(trf     ,wgto,0) ;wgtom
tot_o   = wgt_areaave2(obal    ,wgto,0)

resi_o  = pr_o - ev_o
resi_o2 = pr_o - lh_o
resi_o3 = pr_o - lh_o2

;region -> Global ---------------------------------------

pr_g    = wgt_areaave2(prec    ,wgt,0)
sh_g    = wgt_areaave2(sensh   ,wgt,0)
rn_g    = wgt_areaave2(radnet  ,wgt,0)
lhf_g   = wgt_areaave2(lathf   ,wgt,0)
ev_g    = wgt_areaave2(evp     ,wgt,0)
lh_g    = wgt_areaave2(lath    ,wgt,0)
lh_g2   = wgt_areaave2(lath_snw,wgt,0)
rf_g    = wgt_areaave2(srf     ,wgt,0)
dr_g    = wgt_areaave2(dr      ,wgt,0)
trf_g   = wgt_areaave2(trf     ,wgt,0)
tot_g   = wgt_areaave2(wbal    ,wgt,0)

resi_g  = pr_g-ev_g-rf_g-dr_g       ; NOTE: -rf_g-dr_g  == -rf_l-dr_l 
resi_g2 = pr_g-lh_g-rf_g-dr_g       ; But not for coupled runs bc of ocean
resi_g3 = pr_g-lh_g2-rf_g-dr_g
;resi_g  = pr_g-ev_g-rf_l-dr_l       ; NOTE: -rf_g-dr_g  == -rf_l-dr_l
;resi_g2 = pr_g-lh_g-rf_l-dr_l
;resi_g3 = pr_g-lh_g2-rf_l-dr_l

; Table -----------------------------------------------------------------------

land_vals   = (/pr_l,ev_l,rf_l,dr_l,resi_l,lh_l,resi_l2,lh_l2,resi_l3,sh_l,lhf_l,rn_l/)
ocean_vals  = (/pr_o,ev_o,rf_o,dr_o,resi_o,lh_o,resi_o2,lh_o2,resi_o3,sh_o,lhf_o,rn_o/)
global_vals = (/pr_g,ev_g,rf_g,dr_g,resi_g,lh_g,resi_g2,lh_g2,resi_g3,sh_g,lhf_g,rn_g/)
;global_vals = (/pr_g,ev_o,rf_l,dr_l,resi_o,lh_g,resi_g2,lh_g2,resi_g3/)

land  = "" + decimalPlaces(land_vals,3,True)
ocean = "" + decimalPlaces(ocean_vals,3,True)
globe = "" + decimalPlaces(global_vals,3,True)
;print(land)
;print(ocean)
;print(globe)

; Header-------------------------
ncr1  = (/1,12/) ;9
x1    = (/0.085,0.995/) ;.104
y1    = (/0.800,0.855/)
;text1 = (/"Precip","Ev+Es","SurfRun","Drainage","Res 1", \
text1 = (/"Precip","Ev+s+i","SurfRun","Drainage","Res 1", \
          "Evp(LH)","Res 2", \
          "E(LHsnw)","Res 3",\
          "SH","LH","Rnet"/)

; First Column-------------------
ncr2  = (/3,1/)
x2    = (/0.005,0.085/)
y2    = (/0.415,0.800/)
text2 = (/"Land","Ocean","Global"/)

; Main Table Body----------------
ncr3        = (/3,12/) ;9
x3          = (/0.085,0.995/)
y3          = (/0.415,0.800/)
text_values = (/land,ocean,globe/)

pltName = "wbal_table" + "_" + numYr + "yrs"
wks     = gsn_open_wks("ps",pltName)

res1               = True
res1@txFontHeightF = 0.012
res1@gsFillColor   = (/"deepskyblue2","lightgoldenrod1","slateblue2","lightseagreen","tomato","lightgoldenrod1","tomato","lightgoldenrod1","tomato","lightskyblue","lightskyblue","lightskyblue"/)
res2               = True
res2@txFontHeightF = 0.012
res2@gsFillColor   = (/"forestgreen","deepskyblue2","darkorchid2"/)
res3               = True
;res3@gsnMaximize   = True
res3@txFontHeightF = 0.012
res3@gsFillColor   = (/ (/"white","white","white","white","white","white","white","gray","gray","white","white","white"/), \
                        (/"white","gray","gray","gray","gray","white","white","gray","gray","white","white","white"/),     \
                        (/"white","gray","gray","gray","gray","white","gray","gray","gray","white","white","white"/)   /)
                        ;(/"white","gray","gray","gray","gray","white","white","white","white"/)   /)
                        ;(/"white","gray","white","white","gray","white","white","white","white"/)   /)

print("Creating Water Balance Table")
print("")
gsn_table(wks,ncr1,x1,y1,text1,res1)
gsn_table(wks,ncr2,x2,y2,text2,res2)
gsn_table(wks,ncr3,x3,y3,text_values,res3)

system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("rm "+pltName+".ps")
;system("convert -density 500 -rotate 270 "+pltName+".ps "+pltName+".jpg")

end

