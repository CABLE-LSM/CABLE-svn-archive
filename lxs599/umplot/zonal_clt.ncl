; Lauren Stevens 4 Apr 2013
;----------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
hdir = systemfunc("echo $hdir")
usr = systemfunc("echo $USERID")
loadscript(hdir+"/"+usr+"/umplot/regrid.ncl")
loadscript(hdir+"/"+usr+"/umplot/wgt_func.ncl")
;----------------------------------------------------------------

begin

numYr = systemfunc("echo $YR")
mDir  = systemfunc("echo $TOP")
cDir  = systemfunc("echo $MID")
tDir  = hdir+"/"+usr+"/access/access13/amip" ;systemfunc("echo $MID")
cname = systemfunc("echo $CNM")
resol = systemfunc("echo $RES")

;mland = False
;mland = True

if ( cname .eq. "A" ) then
 rtype = "amip"
 a14 = False
end if
if ( cname.eq. "H" ) then
 rtype = "hist"
 a14 = True
end if
if ( cname.eq. "PI" ) then
 rtype = "pi"
 a14 = True
end if

PltDir    = hdir+"/"+usr+"/umplot/data/"
ObsDir    = hdir+"/"+usr+"/umplot/obs/"
raijinDir = "/short/p66/lxs599/from_dcc/lxs599/umplot/analysis/access13/amip/clouds/"

fmos  = addfile(mDir+"/monthly_means_"+numYr+"yrs.nc","r")
fcab  = addfile(cDir+"/monthly_means_"+numYr+"yrs.nc","r")
fcab3 = addfile(raijinDir+"/clt_yseasmean.nc","r")
;fcab3 = addfile(tDir+"/clt_yseasmean.nc","r")
fobs  = addfile(ObsDir+"/ERAi_clouds_yseasmean.nc","r")
if (resol .eq. 96) then
lsm_file  = addfile(PltDir+"mask_accessV2_n"+resol+".nc","r")
frac_file = addfile(PltDir+"landfrac_N"+resol+".nc","r")
else
lsm_file  = addfile(PltDir+"mask_access_N"+resol+".nc","r")
frac_file = addfile(PltDir+"landfrac_N"+resol+".nc","r")
end if
area_file = addfile(PltDir+"areacella_N"+resol+".nc","r")
era_mask  = addfile(ObsDir+"ERAi_LSM_invertlat.nc","r")

if (isfilevar(fcab,"field30_1")) then
 clt1  = fcab     ->field30_1
else
 clt1  = fcab     ->clt
end if
if (isfilevar(fmos,"field30_1")) then
 clt0  = fmos     ->field30_1 ;clt
else
 clt0  = fmos     ->clt
end if
if (isfilevar(fcab3,"field30_1")) then
 clt3  = fcab3    ->field30_1
else
 clt3  = fcab3    ->clt
end if
clto  = fobs     ->tcc
lsm   = lsm_file ->lsm
lfrac = frac_file->lsm
emsk1 = era_mask ->lsm
garea = area_file->areacella

;clt0 = clt0/100.
clt3 = clt3/100.

;printVarSummary(clt1)
;printVarSummary(clt0)
;printVarSummary(clt3)
;printVarSummary(clto)

emsk           = emsk1*emsk1@scale_factor + emsk1@add_offset
emsk!0         = "time"
emsk!1         = "latitude"
emsk!2         = "longitude"
emsk&latitude  = emsk1&latitude
emsk&longitude = emsk1&longitude
emsk&time      = emsk1&time

zclt1 = dim_avg_n_Wrap(clt1(:,0,:,:),0)
if (dimsizes(dimsizes(clt0)).eq.4) then
 zclt0 = dim_avg_n_Wrap(clt0(:,0,:,:),0)
else
 zclt0 = dim_avg_n_Wrap(clt0(:,:,:),0)
end if
zclt3 = dim_avg_n_Wrap(clt3,0)
zclto = dim_avg_n_Wrap(clto,0)
zclto!0 = "latitude"
zclto!1 = "longitude"
zclt1!0 = "latitude"
zclt1!1 = "longitude"
zclt3!0 = "latitude"
zclt3!1 = "longitude"

wgt     = get_wgt(clt1)
wgte    = get_wgt(clto)
wgto    = wgt*(1-lfrac(0,0,:,:))
wgtl    = wgt*lfrac(0,0,:,:)
wgtlm   = wgt*(lsm(0,0,:,:))
wgtom   = wgt*(1-lsm(0,0,:,:))
wgtel   = wgte*emsk(0,:,:)
wgteo   = wgte*(1-emsk(0,:,:))

;rad     = 4.0*atan(1.0) / 180.
;lat     = clt1&latitude
;lat!0   = "lat"
;lat&lat = clt1&latitude
;late     = clto&latitude
;late!0   = "lat"
;late&lat = clto&latitude

; Zonal Mean ====================

zstdd  = new((/3,145/),typeof(clt0))
zstdd  = 0.0
zstdd(0,:) = zonalAve((zclt0(:,:)))
zstdd(1,:) = zonalAve((zclt1(:,:)))
zstdd(2,:) = (/zonalAve((zclt3(:,:)))/)
zsdtod     = zonalAve(dble2flt(zclto(:,:)))
;print((/zstdd(0,:)/))
;print((/zstdd(1,:)/))
;print((/zstdd(2,:)/))
;print((/dim_avg_n_Wrap(dble2flt(zclt1(:,:)),1)/))

wks3  = gsn_open_wks("ps","zonal_clt_"+numYr+"yrs")
res2  = True
;res2@ = (/0.,.5,1.,1.5/)
res2@trYMinF = 0.0
res2@trYMaxF = 1.2
;res2@xyLineColors           = (/"black","blue","red"/)
res2@xyLineColors           = (/"blue","red","green4"/)
res2@xyLineThicknesses      = (/2,2,2,2/)
res2@xyDashPattern          = 0.
res2@pmLegendDisplayMode    = "Always"            ; turn on legend
res2@pmLegendSide           = "Bottom"               ; Change location of
res2@pmLegendParallelPosF   = .85                 ; move units right
res2@pmLegendOrthogonalPosF = -0.3                ; move units down
res2@pmLegendWidthF         = 0.09                ; Change width and
res2@pmLegendHeightF        = 0.11                ; height of legend.
res2@lgPerimOn              = False               ; turn off box around
res2@lgLabelFontHeightF     = .018                ; label font height
res2@xyExplicitLegendLabels = (/"A1.0","A1.1","A1.3"/)
;res2@xyExplicitLegendLabels = (/"A1.3","A1.1","A1.0"/)
;res2@xyExplicitLegendLabels = (/"ERAi","A1.3","A1.1","A1.0"/)

yzon  = clt1&latitude
yzoo  = clto&latitude

res3= True
;res2@ = (/0.,.5,1.,1.5/)
res3@trYMinF = 0.0
res3@trYMaxF = 1.2
res3@xyLineThicknesses   = 2
;res3@xyLineColors           = "green"
res3@pmLegendDisplayMode    = "Always"            ; turn on legend
res3@pmLegendSide           = "Bottom"               ; Change location of
res3@pmLegendParallelPosF   = .85                 ; move units right
res3@pmLegendOrthogonalPosF = -0.32               ; move units down
res3@pmLegendWidthF         = 0.09                ; Change width and
res3@pmLegendHeightF        = 0.11                ; height of legend.
res3@lgPerimOn              = False               ; turn off box around
res3@lgLabelFontHeightF     = .018                ; label font height
res3@xyExplicitLegendLabels = "ERAi"
res3@tiYAxisString = "" 
res3@tmYRLabelsOn = False
res3@tmXTLabelsOn = False

  print("")
  print("Plotting AMOJ Zonal Cloud Fraction")
  print("")

plot3 = gsn_csm_x2y2(wks3,yzon,yzoo,zstdd(:,:),zsdtod,res2,res3)
;res2@pmLegendSide           = "Top"               ; Change location of
;res3@pmLegendSide           = "Top"               ; Change location of
;res3@pmLegendOrthogonalPosF = -0.22               ; move units down

end

