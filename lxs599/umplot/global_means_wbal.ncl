; Lauren Stevens 02 Sept 2011
; NCL script to produce Water Balances.
; Observations: ERA Interim. 
; Varaibles are : 
; Precip P
; Evap E
; Surface Runoff R
; Drainage D
; Water balance = wbal = P-E-R-D

; ENVIRO VARS ---------------------------------------------------------------------


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;procedure Table()

begin

numYr  = systemfunc("echo $YR")
filDir = systemfunc("echo $DIR")
ModDir = systemfunc("echo $PWD")
resol  = systemfunc("echo $RES")
vMask  = systemfunc("echo $MASK")
;Runid  = systemfunc("echo $RUNID")
;cable  = systemfunc("echo $CABLE")
;moses  = systemfunc("echo $MOSES")

; FILES ---------------------------------------------------------------------------
;moses_file = addfile(moses+"/season_`YR`yrs.nc","r") 
;cable_file = addfile(cable+"/season_`YR`yrs.nc","r")

ObsDir = "/home/cmar/ste69f/umplot/obs/"
PltDir = "/home/cmar/ste69f/umplot/data/"

sea_file = addfile(ModDir+"/seasonal_means_"+numYr+"yrs.nc","r")
if (isfilepresent(ModDir+"/Rseasonal_means_"+numYr+"yrs.nc")) then
rsea_file = addfile(ModDir+"/Rseasonal_means_"+numYr+"yrs.nc","r")
end if
;era_file = addfile(ObsDir+"trderi5tp2mmn1.nc","r")
if (resol.eq.48) then 
lsm_file = addfile(PltDir+"LSM_UMn"+resol+".nc","r")
;oro_file = addfile(PltDir+"orography_UMn"+resol+".nc","r")
else
if (vMask.eq.2) then
lsm_file = addfile(PltDir+"mask_accessV2_n"+resol+".nc","r")
else
lsm_file = addfile(PltDir+"mask_accessV1_n"+resol+".nc","r")
end if
end if

couple = systemfunc("echo $CPL")

; SET-UP --------------------------------------------------------------------------

precip  = sea_file ->precip
if (isfilevar(sea_file,"field1526")) then
evap1   = sea_file ->field1526
evap2   = sea_file ->field1527
end if
lh      = sea_file ->lh
lsm     = lsm_file ->lsm
if(isfilepresent(ModDir+"/Rseasonal_means_"+numYr+"yrs.nc")) then
srunoff = rsea_file ->field1532
drunoff = rsea_file ->field1533
else
srunoff = sea_file ->field1532
drunoff = sea_file ->field1533
end if

trunoff = srunoff
trunoff = trunoff + drunoff

rlam  = 2500000             ;
cls   = 1.1335              ;
rad   = 4.0*atan(1.0)/180.

lat          = precip&latitude
lat!0        = "latitude"
lat&latitude = precip&latitude

precip       = precip*24*60*60
precip@units = "mm/day"

if (isvar("couple")) then
if (couple.eq."y") then
evap         = (lh/2500000)*24*60*60
evap@units   = "mm/day"
else
evap         = (evap1+evap2)*24*60*60
evap@units   = "mm/day"
end if
else
evap         = (evap1+evap2)*24*60*60
evap@units   = "mm/day"
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

prec         = dim_avg_n_Wrap(precip,0)
evp          = dim_avg_n_Wrap(evap,0)
srf          = dim_avg_n_Wrap(srunoff,0)
dr           = dim_avg_n_Wrap(drunoff,0)
trf          = dim_avg_n_Wrap(trunoff,0)
lath         = dim_avg_n_Wrap(lh,0)
lath_snw     = dim_avg_n_Wrap(lh_snw,0)

;wbal         = precip-evap-(srunoff+drunoff)
wbal         = prec-evp-(srf+dr)
obal         = prec-evp

; ANN WBAL -------------------------------------------------------------------------
;region -> Global (Land Only) ---------------------------

Mpr            = where(lsm(:,0,:,:).gt.0,prec,-9999)
Mpr@_FillValue = -9999
pr_l           = wgt_areaave(Mpr,cos(lat*rad),1.0,0)

Mev             = where(lsm(:,0,:,:).gt.0,evp,-9999)
Mev@_FillValue  = -9999
ev_l            = wgt_areaave(Mev,cos(lat*rad),1.0,0)

Mlh             = where(lsm(:,0,:,:).gt.0,lath,-9999)
Mlh2            = where(lsm(:,0,:,:).gt.0,lath_snw,-9999)
Mlh@_FillValue  = -9999
Mlh2@_FillValue = -9999
lh_l            = wgt_areaave(Mlh,cos(lat*rad),1.0,0)
lh_l2           = wgt_areaave(Mlh2,cos(lat*rad),1.0,0)

MRf            = where(lsm(:,0,:,:).gt.0,srf,-9999)
MRf@_FillValue = -9999
rf_l           = wgt_areaave(MRf,cos(lat*rad),1.0,0)

MDr            = where(lsm(:,0,:,:).gt.0,dr,-9999)
MDr@_FillValue = -9999
dr_l           = wgt_areaave(MDr,cos(lat*rad),1.0,0)

MtRf            = where(lsm(:,0,:,:).gt.0,trf,-9999)
MtRf@_FillValue = -9999
trf_l           = wgt_areaave(MtRf,cos(lat*rad),1.0,0)

Mtot            = where(lsm(:,0,:,:).gt.0,wbal,-9999)
Mtot@_FillValue = -9999
tot_l           = wgt_areaave(Mtot,cos(lat*rad),1.0,0)

resi_l  = pr_l-ev_l-rf_l-dr_l
resi_l2 = pr_l-lh_l-rf_l-dr_l 
resi_l3 = pr_l-lh_l2-rf_l-dr_l

;region -> Global (Ocean Only) ---------------------------

Opr  = where(lsm(:,0,:,:).eq.0,prec,-9999)
Opr@_FillValue  = -9999
pr_o = wgt_areaave(Opr,cos(lat*rad),1.0,0)

Oev   = where(lsm(:,0,:,:).eq.0,evp,-9999)
Oev@_FillValue  = -9999
ev_o  = wgt_areaave(Oev,cos(lat*rad),1.0,0)

Olh   = where(lsm(:,0,:,:).eq.0,lath,-9999)
Olh2  = where(lsm(:,0,:,:).eq.0,lath_snw,-9999)
Olh@_FillValue  = -9999
Olh2@_FillValue  = -9999
lh_o  = wgt_areaave(Olh,cos(lat*rad),1.0,0)
lh_o2 = wgt_areaave(Olh2,cos(lat*rad),1.0,0)

ORf  = where(lsm(:,0,:,:).eq.0,srf,-9999)
ORf@_FillValue  = -9999
rf_o = wgt_areaave(ORf,cos(lat*rad),1.0,0)

ODr  = where(lsm(:,0,:,:).eq.0,dr,-9999)
ODr@_FillValue  = -9999
dr_o = wgt_areaave(ODr,cos(lat*rad),1.0,0)

OtRf  = where(lsm(:,0,:,:).eq.0,trf,-9999)
OtRf@_FillValue  = -9999
trf_o = wgt_areaave(OtRf,cos(lat*rad),1.0,0)

Otot  = where(lsm(:,0,:,:).eq.0,obal,-9999)
Otot@_FillValue  = -9999
tot_o = wgt_areaave(Otot,cos(lat*rad),1.0,0)

resi_o  = pr_o - ev_o
resi_o2 = pr_o - lh_o
resi_o3 = pr_o - lh_o2

;;region -> Global ---------------------------------------

Gpr  = prec
pr_g = wgt_areaave(Gpr,cos(lat*rad),1.0,0)

Gev   = evp
ev_g  = wgt_areaave(Gev,cos(lat*rad),1.0,0)

Glh   = lath
Glh2  = lath_snw
lh_g  = wgt_areaave(Glh,cos(lat*rad),1.0,0)
lh_g2 = wgt_areaave(Glh2,cos(lat*rad),1.0,0)

GRf  = srf 
rf_g = wgt_areaave(GRf,cos(lat*rad),1.0,0)

GDr  = dr
dr_g = wgt_areaave(GDr,cos(lat*rad),1.0,0)

GtRf  = trf
trf_g = wgt_areaave(GtRf,cos(lat*rad),1.0,0)

Gtot  = wbal
tot_g = wgt_areaave(Gtot,cos(lat*rad),1.0,0)

resi_g  = pr_g-ev_g-rf_g-dr_g       ; NOTE: -rf_g-dr_g  == -rf_l-dr_l 
resi_g2 = pr_g-lh_g-rf_g-dr_g       ; But not for coupled runs bc of ocean
resi_g3 = pr_g-lh_g2-rf_g-dr_g
;resi_g  = pr_g-ev_g-rf_l-dr_l       ; NOTE: -rf_g-dr_g  == -rf_l-dr_l
;resi_g2 = pr_g-lh_g-rf_l-dr_l
;resi_g3 = pr_g-lh_g2-rf_l-dr_l

; Table -----------------------------------------------------------------------

land_vals   = (/pr_l,ev_l,rf_l,dr_l,resi_l,lh_l,resi_l2,lh_l2,resi_l3/)
ocean_vals  = (/pr_o,ev_o,rf_o,dr_o,resi_o,lh_o,resi_o2,lh_o2,resi_o3/)
global_vals = (/pr_g,ev_g,rf_g,dr_g,resi_g,lh_g,resi_g2,lh_g2,resi_g3/)
;global_vals = (/pr_g,ev_o,rf_l,dr_l,resi_o,lh_g,resi_g2,lh_g2,resi_g3/)

land  = "" + land_vals
ocean = "" + ocean_vals
globe = "" + global_vals
;print(land)
;print(ocean)
;print(globe)

; Header-------------------------
ncr1  = (/1,9/)
x1    = (/0.104,0.995/)
y1    = (/0.800,0.855/)
text1 = (/"Precip","Evap","Surf Runoff","Drainage","Res 1", \
          "Evap(LH)","Res 2", \
          "Evp(LHsnw)","Res 3"/)

; First Column-------------------
ncr2  = (/3,1/)
x2    = (/0.005,0.104/)
y2    = (/0.415,0.800/)
text2 = (/"Land","Ocean","Global"/)

; Main Table Body----------------
ncr3        = (/3,9/)
x3          = (/0.104,0.995/)
y3          = (/0.415,0.800/)
text_values = (/land,ocean,globe/)

pltName = "wbal_table" + "_" + numYr + "yrs"
wks     = gsn_open_wks("ps",pltName)

res1               = True
res1@txFontHeightF = 0.012
res1@gsFillColor   = (/"deepskyblue2","lightgoldenrod1","slateblue2","lightseagreen","tomato","lightgoldenrod1","tomato","lightgoldenrod1","tomato"/)
res2               = True
res2@txFontHeightF = 0.012
res2@gsFillColor   = (/"forestgreen","deepskyblue2","darkorchid2"/)
res3               = True
;res3@gsnMaximize   = True
res3@txFontHeightF = 0.012
res3@gsFillColor   = (/ (/"white","white","white","white","white","white","white","gray","gray"/), \
                        (/"white","gray","gray","gray","gray","white","white","gray","gray"/),     \
                        (/"white","gray","gray","gray","gray","white","gray","gray","gray"/)   /)
                        ;(/"white","gray","gray","gray","gray","white","white","white","white"/)   /)
                        ;(/"white","gray","white","white","gray","white","white","white","white"/)   /)

print("Creating Water Balance Table")
gsn_table(wks,ncr1,x1,y1,text1,res1)
gsn_table(wks,ncr2,x2,y2,text2,res2)
gsn_table(wks,ncr3,x3,y3,text_values,res3)

system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("rm "+pltName+".ps")
;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")

end

