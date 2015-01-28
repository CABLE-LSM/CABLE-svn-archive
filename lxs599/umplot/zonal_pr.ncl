; Lauren Stevens 4 Apr 2013
;----------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
usr = systemfunc("echo $USERID")
loadscript("~"+usr+"/umplot/regrid.ncl")
loadscript("~"+usr+"/umplot/wgt_func.ncl")
;----------------------------------------------------------------

begin

numYr = systemfunc("echo $YR") 
;wDir  = systemfunc("echo $PWD")
mDir  = systemfunc("echo $MOSES")
cDir  = systemfunc("echo $CABLE")
tDir  = "~ste69f/access/access13/amip" ;systemfunc("echo $CABLE3")
;tDir  = "~"+usr+"/access/access13/amip" ;systemfunc("echo $CABLE3")
cname = systemfunc("echo $CNM")

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

PltDir    = "~"+usr+"/umplot/data/"
ObsDir    = "~"+usr+"/umplot/obs/"

fmos  = addfile(mDir+"/monthly_means_"+numYr+"yrs.nc","r")
fcab  = addfile(cDir+"/monthly_means_"+numYr+"yrs.nc","r")
fcab3 = addfile(tDir+"/monthly_means_27yrs.nc","r")
;fcab3 = addfile("/home/cmar/ste69f/access/access13/"+rtype+"/yearly_means_27yrs.nc","r")
fobs  = addfile(ObsDir+"ERAi_precip.nc","r")
lsm_file  = addfile(PltDir+"mask_accessV2_n96.nc","r")
frac_file = addfile(PltDir+"landfrac_v2.nc","r")
area_file = addfile(PltDir+"areacella_fx_ACCESS1-0_historical_r0i0p0.nc","r")
era_mask  = addfile(ObsDir+"ERAi_LSM_invertlat.nc","r")

if (isfilevar(fcab,"precip")) then
 pr1  = fcab    ->precip
else
 if (isfilevar(fcab,"tot_precip")) then
  pr1  = fcab    ->tot_precip
 else
  pr1  = fcab    ->pr ; field5226
 end if
end if
if (isfilevar(fmos,"precip")) then
 pr0  = fmos     ->precip
else
 if (isfilevar(fmos,"tot_precip")) then
  pr0  = fmos     ->tot_precip
 else
  pr0  = fmos     ->pr ; field5226
 end if
end if
if (isfilevar(fcab3,"precip")) then
 pr3  = fcab3    ->precip
else
 if (isfilevar(fcab3,"tot_precip")) then
  pr3  = fcab3    ->tot_precip
 else
  pr3  = fcab3    ->pr
 end if
end if
pro  = fobs     ->a
lsm   = lsm_file ->lsm
lfrac = frac_file->lsm
emsk1 = era_mask ->lsm
garea = area_file->areacella

pr1 = pr1*24*60*60
pr0 = pr0*24*60*60
pr3 = pr3*24*60*60

;printVarSummary(pr1)
;printVarSummary(pr0)
;printVarSummary(pr3)
;printVarSummary(pro)

emsk           = emsk1*emsk1@scale_factor + emsk1@add_offset
emsk!0         = "time"
emsk!1         = "latitude"
emsk!2         = "longitude"
emsk&latitude  = emsk1&latitude
emsk&longitude = emsk1&longitude
emsk&time      = emsk1&time

zpr1 = dim_avg_n_Wrap(pr1(:,0,:,:),0)
if (dimsizes(dimsizes(pr0)).eq.4) then
 zpr0 = dim_avg_n_Wrap(pr0(:,0,:,:),0)
else
 zpr0 = dim_avg_n_Wrap(pr0(:,:,:),0)
end if
zpr3 = dim_avg_n_Wrap(pr3,0)
zpro = dim_avg_n_Wrap(pro,0)
;printVarSummary(zpro)
zpro!0 = "latitude"
zpro!1 = "longitude"
zpr1!0 = "latitude"
zpr1!1 = "longitude"
zpr3!0 = "latitude"
zpr3!1 = "longitude"

wgt     = get_wgt(pr1)
wgte    = get_wgt(pro)
wgto    = wgt*(1-lfrac(0,0,:,:))
wgtl    = wgt*lfrac(0,0,:,:)
wgtlm   = wgt*(lsm(0,0,:,:))
wgtom   = wgt*(1-lsm(0,0,:,:))
wgtel   = wgte*emsk(0,:,:)
wgteo   = wgte*(1-emsk(0,:,:))

;lat     = pr1&latitude
;lat!0   = "lat"
;lat&lat = pr1&latitude
;late     = pro&latitude
;late!0   = "lat"
;late&lat = pro&latitude

; Zonal Mean ====================

zmean  = new((/3,145/),typeof(pr0))
zmean  = 0.0
zmean(0,:) = zonalAve((zpr0(:,:)))
zmean(1,:) = zonalAve((zpr1(:,:)))
zmean(2,:) = (/zonalAve((zpr3(:,:)))/)
zobmn      = zonalAve((zpro(:,:)))
;print((/zmean(0,:)/))
;print((/zmean(1,:)/))
;print((/zmean(2,:)/))
;print((/dim_avg_n_Wrap(dble2flt(zpr1(:,:)),1)/))

wks3  = gsn_open_wks("ps","zonal_pr_"+numYr+"yrs")
res2  = True
;res2@ = (/0.,.5,1.,1.5/)
res2@trYMinF = 0.0
res2@trYMaxF = 8.
;res2@xyLineColors           = (/"black","blue","red"/)
res2@xyLineColors           = (/"blue","red","green4"/)
res2@xyLineThicknesses      = (/2,2,2,2/)
res2@xyDashPattern          = 0.
res2@pmLegendDisplayMode    = "Always"            ; turn on legend
res2@pmLegendSide           = "Top"               ; Change location of
res2@pmLegendParallelPosF   = .85                 ; move units right
res2@pmLegendOrthogonalPosF = -0.3                ; move units down
res2@pmLegendWidthF         = 0.09                ; Change width and
res2@pmLegendHeightF        = 0.11                ; height of legend.
res2@lgPerimOn              = False               ; turn off box around
res2@lgLabelFontHeightF     = .018                ; label font height
res2@xyExplicitLegendLabels = (/"A1.0","A1.1","A1.3"/)
;res2@xyExplicitLegendLabels = (/"A1.3","A1.1","A1.0"/)
;res2@xyExplicitLegendLabels = (/"ERAi","A1.3","A1.1","A1.0"/)

yzon  = pr1&latitude
yzoo  = pro&latitude

res3= True
;res2@ = (/0.,.5,1.,1.5/)
res3@trYMinF = 0.0
res3@trYMaxF = 8.
res3@xyLineThicknesses   = 2
;res3@xyLineColors           = "green"
res3@pmLegendDisplayMode    = "Always"            ; turn on legend
res3@pmLegendSide           = "Top"               ; Change location of
res3@pmLegendParallelPosF   = .85                 ; move units right
res3@pmLegendOrthogonalPosF = -0.20               ; move units down
res3@pmLegendWidthF         = 0.09                ; Change width and
res3@pmLegendHeightF        = 0.11                ; height of legend.
res3@lgPerimOn              = False               ; turn off box around
res3@lgLabelFontHeightF     = .018                ; label font height
res3@xyExplicitLegendLabels = "ERAi"
res3@tiYAxisString = "" 
res3@tmYRLabelsOn = False
res3@tmXTLabelsOn = False

  print("")
  print("Plotting AMOJ Zonal Precip")
  print("")
plot3 = gsn_csm_x2y2(wks3,yzon,yzoo,zmean(:,:),zobmn,res2,res3)
;res2@pmLegendSide           = "Top"               ; Change location of
;res3@pmLegendSide           = "Top"               ; Change location of
;res3@pmLegendOrthogonalPosF = -0.22               ; move units down

end

