;==============================================================================
; Lauren Stevens 19 July 2011
;==============================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;usr      = systemfunc("echo $USERID")
;loadscript("~"+usr+"/umplot/plot_mmdc.ncl")
;==============================================================================

begin

usr      = systemfunc("echo $USERID")
numYr    = systemfunc("echo $YR")
modl     = systemfunc("echo $MODEL")
filDir   = systemfunc("echo $PWD")
;filDir   = systemfunc("echo $DIRW")
jpg      = systemfunc("echo $CHFMT")
if(jpg.eq."y") then
 jpeg = True
else
 jpeg = False
end if

tsDir    = "~"+usr+"/umplot/obs/"
home     = "~"+usr+"/"
flxobs   = "/home/cmar/ste69f/obsdata/fluxnet/PALS/flxnt_mon/"

;obs ---
 flux_ma = addfile(tsDir  + "MA96_lu_u0_mm_flx.nc"    ,"r")
 met_ma  = addfile(tsDir  + "MA96_mm_met.nc"          ,"r")
 flux_tb = addfile(tsDir  + "AnnCycles_tumbarumba.nc" ,"r")
 hyy_new = addfile(flxobs + "mmdc_ts_Hyytiala.nc"     ,"r")
 lo_new  = addfile(flxobs + "mmdc_ts_Loobos.nc"       ,"r")
 nb_new  = addfile(flxobs + "mmdc_ts_Boreas.nc"       ,"r")
 bv_new  = addfile(flxobs + "mmdc_ts_Bondv.nc"        ,"r")
 th_new  = addfile(flxobs + "mmdc_ts_Tharandt.nc"     ,"r")
 hv_new  = addfile(flxobs + "mmdc_ts_Harvard.nc"      ,"r")
 tb_new  = addfile(flxobs + "mmdc_ts_Tumba.nc"        ,"r")
; model ---
 mosfil  = addfile(filDir+"/AnnualCycle_"+numYr+"yrs.nc" ,"r")

;OBS===========================================================================

; order: bv hv hy lo lw nb th vi wb ma tb hay
;        0   1  2  3  4  5  6 7  8  9  10 11
fils   = systemfunc("ls "+tsDir+"*-*_lu_u0_mm_flx.nc")
f      = addfiles(fils,"r")
ListSetType(f,"join")
sens   = f[:]->H  ;(0:11)
lat    = f[:]->LH ;(0:11)

fils2  = systemfunc("ls "+tsDir+"*-*_mm_met.nc")
f2     = addfiles(fils2,"r")
ListSetType(f2,"join")
radnet = f2[:]->Rn  ;(0:11)
;temp2m = f2[:]->Ta  ;(0:11)
;preci  = f2[:]->PPT ;(0:11)

sens_ma = flux_ma ->H
lat_ma  = flux_ma ->LH
rad_ma  = met_ma  ->Rn
;temp_ma = met_ma  ->Ta
;prec_ma = met_ma  ->PPT
sens_tb = flux_tb ->qh_8yr_mean    ;Qh    ;H
lat_tb  = flux_tb ->qle_8yr_mean   ;Qle   ;LH
rad_tb  = flux_tb ->rnet_8yr_mean  ;Rnet  ;Rn
;temp_tb = flux_tb ->tair_8yr_mean  ;Tair  ;Ta
;prec_tb = flux_tb ->rainf_8yr_mean ;Rainf ;PPT
;sens_tb = tb_new ->Qh
;lat_tb  = tb_new ->Qle
;rad_tb  = tb_new ->Rnet
;temp_tb = tb_new ->Tair
;prec_tb = tb_new ->Rainf
sens_hy = hyy_new ->Qh
lat_hy  = hyy_new ->Qle
rad_hy  = hyy_new ->Rnet
;temp_hy = hyy_new ->Tair
;prec_hy = hyy_new ->Rainf
sens_lo = lo_new  ->Qh
lat_lo  = lo_new  ->Qle
rad_lo  = lo_new  ->Rnet
;temp_lo = lo_new  ->Tair
;prec_lo = lo_new  ->Rainf
sens_nb = nb_new  ->Qh
lat_nb  = nb_new  ->Qle
rad_nb  = nb_new  ->Rnet
;temp_nb = nb_new  ->Tair
;prec_nb = nb_new  ->Rainf
sens_th = th_new  ->Qh
lat_th  = th_new  ->Qle
;rad_th  = th_new  ->Rnet ; missing
;temp_th = th_new  ->Tair
;prec_th = th_new  ->Rainf
sens_hv = hv_new  ->Qh
lat_hv  = hv_new  ->Qle
rad_hv  = hv_new  ->Rnet
;temp_hv = hv_new  ->Tair
;prec_hv = hv_new  ->Rainf
sens_bv = bv_new  ->Qh
lat_bv  = bv_new  ->Qle
rad_bv  = bv_new  ->Rnet
;temp_bv = bv_new  ->Tair
;prec_bv = bv_new  ->Rainf

sens_tb!0 = "axis_1"
copy_VarCoords(sens_tb,lat_tb)
copy_VarCoords(sens_tb,rad_tb)
sens_hy!0 = "axis_1"
sens_hy!1 = "axis_2"
sens_hy!2 = "axis_3"
copy_VarCoords(sens_hy,lat_hy)
copy_VarCoords(sens_hy,rad_hy)
copy_VarCoords(sens_hy,sens_lo)
copy_VarCoords(sens_hy, lat_lo)
copy_VarCoords(sens_hy, rad_lo)
copy_VarCoords(sens_hy,sens_hv)
copy_VarCoords(sens_hy, lat_hv)
copy_VarCoords(sens_hy, rad_hv)
copy_VarCoords(sens_hy,sens_bv)
copy_VarCoords(sens_hy, lat_bv)
copy_VarCoords(sens_hy, rad_bv)
copy_VarCoords(sens_hy,sens_nb)
copy_VarCoords(sens_hy, lat_nb)
copy_VarCoords(sens_hy, rad_nb)
copy_VarCoords(sens_hy,sens_th)
copy_VarCoords(sens_hy, lat_th)

;OBS===========================================================================

sensht            = new((/11,12/),double)
latht             = new((/11,12/),double)
rad               = new((/11,12/),double)
;temp              = new((/11,12/),double)
;precip            = new((/11,12/),double)
sensht            = -9999
latht             = -9999
rad               = -9999
;temp              = -9999
;precip            = -9999
sensht@_FillValue = -9999
latht@_FillValue  = -9999
rad@_FillValue    = -9999
;temp@_FillValue   = -9999
;precip@_FillValue = -9999
sensht!0 = "axis_0"
 latht!0 = "axis_0"
   rad!0 = "axis_0"
sensht!1 = "axis_1"
 latht!1 = "axis_1"
   rad!1 = "axis_1"

do k=0,8,1
 sensht(k,:)       = sens(k,0:11)*1000000/(31*24*60*60)   ;H
 latht(k,:)        = lat(k,0:11)*1000000/(31*24*60*60)    ;LH
 rad(k,:)          = radnet(k,0:11)*1000000/(31*24*60*60) ;Rn
 ;temp(k,:)         = temp2m(k,0:11)                       ;Ta
 ;precip(k,:)       = preci(k,0:11)                        ;PPT
end do

sensht(9,:)       = sens_ma(0:11)*1000000/(31*24*60*60)   ;H
latht(9,:)        = lat_ma(0:11)*1000000/(31*24*60*60)    ;LH
rad(9,:)          = rad_ma(0:11)*1000000/(31*24*60*60)    ;Rn
;temp(9,:)         = temp_ma(0:11)                         ;Ta
;precip(9,:)       = prec_ma(0:11)                         ;PPT

sensht(10,:)      = sens_tb(0:11)                         ;H
latht(10,:)       = lat_tb(0:11)                          ;LH
rad(10,:)         = rad_tb(0:11)                          ;Rn
;temp(10,:)        = temp_tb(0:11)                         ;Ta
;precip(10,:)      = prec_tb(0:11)                         ;PPT

sensht(2,:)       = dim_avg_n_Wrap(sens_hy(:,:,0),1) ;H
latht(2,:)        = dim_avg_n_Wrap(lat_hy(:,:,0),1)  ;LH
rad(2,:)          = dim_avg_n_Wrap(rad_hy(:,:,0),1)  ;Rn
;temp(2,:)         = dim_avg_n_Wrap(temp_hy(:,:,0),1) ;Ta
;precip(2,:)       = dim_avg_n_Wrap(prec_hy(:,:,0),1) ;PPT
sensht(3,:)       = dim_avg_n_Wrap(sens_lo(:,:,0),1) ;H
latht(3,:)        = dim_avg_n_Wrap(lat_lo(:,:,0),1)  ;LH
rad(3,:)          = dim_avg_n_Wrap(rad_lo(:,:,0),1)  ;Rn
;temp(3,:)         = dim_avg_n_Wrap(temp_lo(:,:,0),1) ;Ta
;precip(3,:)       = dim_avg_n_Wrap(prec_lo(:,:,0),1) ;PPT
sensht(5,:)       = dim_avg_n_Wrap(sens_nb(:,:,0),1) ;H
latht(5,:)        = dim_avg_n_Wrap(lat_nb(:,:,0),1)  ;LH
rad(5,:)          = dim_avg_n_Wrap(rad_nb(:,:,0),1)  ;Rn
;temp(5,:)         = dim_avg_n_Wrap(temp_nb(:,:,0),1) ;Ta
;precip(5,:)       = dim_avg_n_Wrap(prec_nb(:,:,0),1) ;PPT
sensht(0,:)       = dim_avg_n_Wrap(sens_bv(:,:,0),1) ;H
latht(0,:)        = dim_avg_n_Wrap(lat_bv(:,:,0),1)  ;LH
rad(0,:)          = dim_avg_n_Wrap(rad_bv(:,:,0),1)  ;Rn
;temp(0,:)         = dim_avg_n_Wrap(temp_bv(:,:,0),1) ;Ta
;precip(0,:)       = dim_avg_n_Wrap(prec_bv(:,:,0),1) ;PPT
sensht(1,:)       = dim_avg_n_Wrap(sens_hv(:,:,0),1) ;H
latht(1,:)        = dim_avg_n_Wrap(lat_hv(:,:,0),1)  ;LH
rad(1,:)          = dim_avg_n_Wrap(rad_hv(:,:,0),1)  ;Rn
;temp(1,:)         = dim_avg_n_Wrap(temp_hv(:,:,0),1) ;Ta
;precip(1,:)       = dim_avg_n_Wrap(prec_hv(:,:,0),1) ;PPT
sensht(6,:)       = dim_avg_n_Wrap(sens_th(:,:,0),1) ;H
latht(6,:)        = dim_avg_n_Wrap(lat_th(:,:,0),1)  ;LH
;rad(6,:)          = dim_avg_n_Wrap(rad_th(:,:,0),1)  ;Rn
;temp(6,:)         = dim_avg_n_Wrap(temp_th(:,:,0),1) ;Ta
;precip(6,:)       = dim_avg_n_Wrap(prec_th(:,:,0),1) ;PPT

;;Lest 7/1/14 - For AMOJ papers
;sensht(9,:)       = -9999;sens_ma(0:11)*1000000/(31*24*60*60)   ;H
;latht(9,:)        = -9999;lat_ma(0:11)*1000000/(31*24*60*60)    ;LH
;rad(9,:)          = -9999;rad_ma(0:11)*1000000/(31*24*60*60)    ;Rn
;temp(9,:)         = -9999;temp_ma(0:11)                         ;Ta
;precip(9,:)       = -9999;prec_ma(0:11)                         ;PPT
;sensht(10,:)      = -9999;sens_tb(0:11)                         ;H
;latht(10,:)       = -9999;lat_tb(0:11)                          ;LH
;rad(10,:)         = -9999;rad_tb(0:11)                          ;Rn
;temp(10,:)        = -9999;temp_tb(0:11)                         ;Ta
;precip(10,:)      = -9999;prec_tb(0:11)                         ;PPT
;sensht(7,:)       = -9999 ; vie
;latht(7,:)        = -9999 ; vie
;rad(7,:)          = -9999 ; vie
;temp(7,:)         = -9999 ; vie
;precip(7,:)       = -9999 ; vie
;;Lest 7/1/14 - For AMOJ papers

;MODEL==========================================

if (isfilevar(mosfil,"hfss")) then
 msh     = mosfil ->hfss      ; sh
 mlh     = mosfil ->hfls      ; lh
 mrad    = mosfil ->field3333 ; field202
 ;mta     = mosfil ->tas       ; temp
 ;mpr     = mosfil ->field5226 ; tot_precip
else
 msh     = mosfil ->sh
 mlh     = mosfil ->lh
 mrad    = mosfil ->field202
 ;mta     = mosfil ->temp
 ;mpr     = mosfil ->tot_precip
end if

if (isfilevar(mosfil,"solar")) then
 mrad_1 = mosfil ->solar    ;field3333
 mrad_2 = mosfil ->longwave ;field3333
 mrad1  = mrad_1+mrad_2
 copy_VarMeta(mrad_1,mrad1)
end if

if (isfilevar(mosfil,"hfss")) then
dNames  = getfilevardims(mosfil,"hfss") ;sh
dNlh    = getfilevardims(mosfil,"hfls") ;lh
dNrad   = getfilevardims(mosfil,"field3333")
else
dNames  = getfilevardims(mosfil,"sh")
dNlh    = getfilevardims(mosfil,"lh")
dNrad   = getfilevardims(mosfil,"field202")
end if

; order in file: 31 sites currently
;man,hay,dwaters,bats,india,hyy,afN,afS,tumb,thar,bond,ding,lwas,walk,loo,harv,bore,viel
;ant1,ant2,ant3,rus1,alk,afr,grl,ant,sau1,sau2,sau3,rus2,sau4

mrad_ro = new((/11,12/),double)
if (isvar("mrad1") .and. sum(mrad).eq.0.) then ; if (modl.eq. "m") then
 mrad_ro(0,:)  = mrad1($dNames(1)$|10,$dNames(0)$|:) ; bv
 mrad_ro(1,:)  = mrad1($dNames(1)$|15,$dNames(0)$|:) ; hv
 mrad_ro(2,:)  = mrad1($dNames(1)$|5 ,$dNames(0)$|:) ; hy
 mrad_ro(3,:)  = mrad1($dNames(1)$|14,$dNames(0)$|:) ; lo
 mrad_ro(4,:)  = mrad1($dNames(1)$|12,$dNames(0)$|:) ; lw
 mrad_ro(5,:)  = mrad1($dNames(1)$|16,$dNames(0)$|:) ; nb
 mrad_ro(6,:)  = mrad1($dNames(1)$|9 ,$dNames(0)$|:) ; th
 mrad_ro(7,:)  = mrad1($dNames(1)$|17,$dNames(0)$|:) ; vi
 mrad_ro(8,:)  = mrad1($dNames(1)$|13,$dNames(0)$|:) ; wb
 mrad_ro(9,:)  = mrad1($dNames(1)$|0 ,$dNames(0)$|:) ; ma
 mrad_ro(10,:) = mrad1($dNames(1)$|8 ,$dNames(0)$|:) ; tb
else
 mrad_ro(0,:)  = mrad($dNames(1)$|10,$dNames(0)$|:) ; bv
 mrad_ro(1,:)  = mrad($dNames(1)$|15,$dNames(0)$|:) ; hv
 mrad_ro(2,:)  = mrad($dNames(1)$|5 ,$dNames(0)$|:) ; hy
 mrad_ro(3,:)  = mrad($dNames(1)$|14,$dNames(0)$|:) ; lo
 mrad_ro(4,:)  = mrad($dNames(1)$|12,$dNames(0)$|:) ; lw
 mrad_ro(5,:)  = mrad($dNames(1)$|16,$dNames(0)$|:) ; nb
 mrad_ro(6,:)  = mrad($dNames(1)$|9 ,$dNames(0)$|:) ; th
 mrad_ro(7,:)  = mrad($dNames(1)$|17,$dNames(0)$|:) ; vi
 mrad_ro(8,:)  = mrad($dNames(1)$|13,$dNames(0)$|:) ; wb
 mrad_ro(9,:)  = mrad($dNames(1)$|0 ,$dNames(0)$|:) ; ma
 mrad_ro(10,:) = mrad($dNames(1)$|8 ,$dNames(0)$|:) ; tb
end if

msh_ro = new((/11,12/),double)
msh_ro(0,:)  = msh($dNames(1)$|10,$dNames(0)$|:) ; bv
msh_ro(1,:)  = msh($dNames(1)$|15,$dNames(0)$|:) ; hv
msh_ro(2,:)  = msh($dNames(1)$|5 ,$dNames(0)$|:) ; hy
msh_ro(3,:)  = msh($dNames(1)$|14,$dNames(0)$|:) ; lo
msh_ro(4,:)  = msh($dNames(1)$|12,$dNames(0)$|:) ; lw
msh_ro(5,:)  = msh($dNames(1)$|16,$dNames(0)$|:) ; nb
msh_ro(6,:)  = msh($dNames(1)$|9 ,$dNames(0)$|:) ; th
msh_ro(7,:)  = msh($dNames(1)$|17,$dNames(0)$|:) ; vi
msh_ro(8,:)  = msh($dNames(1)$|13,$dNames(0)$|:) ; wb
msh_ro(9,:)  = msh($dNames(1)$|0 ,$dNames(0)$|:) ; ma
msh_ro(10,:) = msh($dNames(1)$|8 ,$dNames(0)$|:) ; tb

mlh_ro = new((/11,12/),double)
mlh_ro(0,:)  = mlh($dNames(1)$|10,$dNames(0)$|:) ; bv
mlh_ro(1,:)  = mlh($dNames(1)$|15,$dNames(0)$|:) ; hv
mlh_ro(2,:)  = mlh($dNames(1)$|5 ,$dNames(0)$|:) ; hy
mlh_ro(3,:)  = mlh($dNames(1)$|14,$dNames(0)$|:) ; lo
mlh_ro(4,:)  = mlh($dNames(1)$|12,$dNames(0)$|:) ; lw
mlh_ro(5,:)  = mlh($dNames(1)$|16,$dNames(0)$|:) ; nb
mlh_ro(6,:)  = mlh($dNames(1)$|9 ,$dNames(0)$|:) ; th
mlh_ro(7,:)  = mlh($dNames(1)$|17,$dNames(0)$|:) ; vi
mlh_ro(8,:)  = mlh($dNames(1)$|13,$dNames(0)$|:) ; wb
mlh_ro(9,:)  = mlh($dNames(1)$|0 ,$dNames(0)$|:) ; ma
mlh_ro(10,:) = mlh($dNames(1)$|8 ,$dNames(0)$|:) ; tb

;==============================================================================

sh_ANN            = new((/11,12/),double)
;sh_ANN            = msh_ro ;(axis_4|0:10,axis_1|:)
sh_ANN            = msh_ro - sensht
lh_ANN            = new((/11,12/),double)
;lh_ANN            = mlh_ro ;(axis_4|0:10,axis_1|:)
lh_ANN            = mlh_ro - latht
rad_ANN           = new((/11,12/),double)
;rad_ANN           = mrad_ro ;(axis_4|0:10,axis_1|:)
rad_ANN           = mrad_ro - rad 

;----------------------------------------------------------------------------

;y = fspan(0,23.5,48)
y = ispan(1,12,1)

;Plot-------------------------------------------------------------------------

print("")
print("Plotting MM-AnnualC Biases for All Sites")

pltName = "annualcycle_"+numYr+"yrs"
wks     = gsn_open_wks("ps","annualcycle_"+numYr+"yrs")
res                    = True
res@gsnFrame           = False
res@gsnDraw            = False
res@xyLineColors       = (/"black","deepskyblue1","red","green","mediumorchid1","cyan","magenta",\
                           "firebrick","orange","SlateBlue","navy","SkyBlue"/)
res@xyDashPatterns     = (/0,0,0,0,0,0,0,0,0,0,1,2/)
res@xyLineThicknesses  = (/2,2,2,2,2,2,2,2,2,2,2,2/)
res@gsnYRefLine        = 0.
;res@gsnLeftString      = ""
;res@tiYAxisString      = "W/m2"

res@tiYAxisFontHeightF = 0.02
res@tiYAxisOffsetXF    = 0.005
res@tiYAxisString      = "W/m2"
res@tiXAxisString      = "months"
if (modl.eq."m") then
 res@gsnLeftString     = "Annual - A1.0"
else
 res@gsnLeftString     = "Annual - A1.1"
end if
res@gsnLeftStringFontHeightF    =  0.02
res@gsnLeftStringOrthogonalPosF = -0.10
res@gsnLeftStringParallelPosF   =  0.05

mxf    = (/0.,0.,0./)
;mxf!0  = "axis_0"
;mxf&axis_0 = (/1,2,3/)
mxf(0) = dble2flt(max(rad_ANN))
mxf(1) = dble2flt(max(sh_ANN))
mxf(2) = dble2flt(max(lh_ANN))
mnf    = (/0.,0.,0./)
;mnf!0  = "axis_0"
;mnf&axis_0 = (/1,2,3/)
mnf(0) = dble2flt(min(rad_ANN))
mnf(1) = dble2flt(min(sh_ANN))
mnf(2) = dble2flt(min(lh_ANN))
res@trYMinF                = floor(min(mnf))- 5 ;-90
res@trYMaxF                = ceil(max(mxf)) + 5 ;120

;-----------------------------------------------------------------------------

;res@gsnMaximize = True
plot = new(3,graphic)
plot(0) = gsn_csm_xy(wks,y,rad_ANN,res)
plot(1) = gsn_csm_xy(wks,y,sh_ANN,res)
 res@pmLegendDisplayMode    = "Always"            ; turn on legend
 res@pmLegendSide           = "Bottom"            ; Change location of
 res@pmLegendParallelPosF   = .15                 ; move units right
 res@pmLegendOrthogonalPosF = -0.55               ; move units down
 res@pmLegendWidthF         = 0.16                ; Change width and
 res@pmLegendHeightF        = 0.22                ; height of legend.
 res@lgPerimOn              = False               ; turn off box around
 res@lgLabelFontHeightF     = .015                ; label font height
 res@xyExplicitLegendLabels = (/"bond","harv","hyyt","loob","lwas","nbor","thar","vie","wbra","man","tumb","hay"/)
plot(2) = gsn_csm_xy(wks,y,lh_ANN,res)

; bv hv hy lo lw nb th vi wb ma tb
; Legend ---------------------------------------------------------------------

;colors = (/"black","blue","red","green","purple","cyan","magenta","firebrick","orange","SlateBlue","navy","SkyBlue"/)
;labels = (/"bondville","harvardforest","hyytiala","loobos","littlewashita","nsaboreas","tharandt","vielsalm","walkerbranch","manaus","tumbarumba"/)

;lgres                    = True
;lgres@vpWidthF           = 0.2                         ; width of legend (NDC)
;lgres@vpHeightF          = 0.5                         ; height of legend (NDC)
;;lgres@lgLineColors       = (/"black","blue","red","green","purple","cyan","magenta","firebrick","orange","SlateBlue","navy","SkyBlue"/)    ; colors for legend lines
;;lgres@lgDashIndex       = 0.
;lgres@lgPerimOn          = False                       ; turn off perimeter
;;lgres@lgDashIndexes     = (/0,0,0,0,0,0,0,0,0,0,1,2/)  ; dash indexes
;;lgres@lgLineLabelStrings = (/"",""/)                   ; no labels in lines

;xpos = (/0.05,0.05,0.05,0.25,0.25,0.25,0.5,0.5,0.5,0.75,0.75/)
;ypos = (/0.05,.1,.2,0.05,.1,.2,0.05,.1,.2,0.05,.1/)

;lgd                      = gsn_create_legend(wks,11,lst,lgres)

;do i=0,10,1
;lgres@lgLineColors = colors(i)
;gsn_legend_ndc(wks,1,labels(i),xpos(i),ypos(i),lgres)
;;gsn_legend_ndc(wks,11,labels,xpos,ypos,lgres)
;end do

;addres                   = True
;addres@amParallelPosF    = 0.05
;addres@amOrthogonalPosF  = 1.05
;add                      = gsn_add_annotation(plot(1),lgd,addres)

;-----------------------------------------------------------------------------

pres = True
pres@gsnMaximize = True
pres@gsnPanelFigureStrings = (/"Rnet","SH","LH"/)
pres@gsnPanelFigureStringsBackgroundFillColor = "Transparent"
;pres@gsnFrame = False
;pres@gsnPanelBottom = 0.13

gsn_panel(wks,plot,(/1,3/),pres)

;-----------------------------------------------------------------------------

;legend = create "Legend" legendClass wks
;"vpXF"                     : 0.26
;"vpYF"                     : 0.37
;"vpWidthF"                 : 0.5                     ; width
;"vpHeightF"                : 0.3                     ; height
;"lgPerimOn"                : False                   ; no perimeter
;"lgItemCount"              : 11                      ; how many
;"lgLineLabelStrings"       : (/"bondville","harvardforest","hyytiala","loobos","littlewashita","nsaboreas","tharandt","vielsalm","walkerbranch","manaus","tumbarumba"/)  ; labels
;"lgLabelsOn"               : False                   ; no default labels
;"lgLineLabelFontHeightF"   : 0.015                   ; font height
;"lgDashIndexes"            : (/0,0,0,0,0,0,0,0,0,0,1,1/)            ; line paters
;"lgLineColors"             : (/"black","blue","red","green","purple","cyan","magenta","firebrick","orange","SlateBlue","navy","SkyBlue"/)
;"lgMonoLineLabelFontColor" : True                    ; one label color
;"lgItemOrder"              : (/10,9,8,7,6,5,4,3,2,1,0/)
;;"lgItemPlacement"          : "ExplicitPlacement"
;;"lgItemPositions "         : (/0.1,0.2,0.3,0.1,0.2,0.3,0.1,0.2,0.3,0.1,0.2/)
;end create

colors = (/"black","deepskyblue1","red","green","mediumorchid1","cyan","magenta","firebrick","orange","SlateBlue","navy","forestgreen"/)
;labels = (/"bondville","harvardforest","hyytiala","loobos","littlewashita","nsaboreas","tharandt","vielsalm","walkerbranch","manaus","tumbarumba"/)
labels = (/"Tumbarumba","Manaus","WalkerBranch","Vielsalm","Tharandt","NsaBoreas","LittleWashita","Loobos","Hyytiala","HarvardForest","Bondville"/)

lgres                    = True
lgres@vpWidthF           = 0.6                         ; width of legend (NDC)
lgres@vpHeightF          = 0.7                         ; height of legend (NDC)
;;lgres@lgDashIndex       = 0.
;lgres@lgPerimOn          = False                       ; turn off perimeter
lgres@lgDashIndexes     = (/0,0,0,0,0,0,0,0,0,0,1,1/)  ; dash indexes
lgres@lgLineThicknesses  = (/2,2,2,2,2,2,2,2,2,2,2/)
lgres@lgLineLabelStrings = (/"","","","","","","","","","",""/)                   ; no labels in lines
;lgres@lgLineLabelStrings =  labels
lgres@lgItemOrder        = (/10,9,8,7,6,5,4,3,2,1,0/)
lgres@lgOrientation      = "Horizontal"
lgres@lgLabelAngleF      = 270 
;lgres@gsnMaximize         = True
;lgres@gsnPaperOrientation = "landscape"
;lgres@wkOrientation = "seascape"
lgres@lgBoxMinorExtentF      = 0.2

xpos = 0.30  ;(/0.05,0.05,0.05,0.25,0.25,0.25,0.5,0.5,0.5,0.75,0.75/)
ypos = 0.95 ;(/0.05,.1,.2,0.05,.1,.2,0.05,.1,.2,0.05,.1/)

;lgd                      = gsn_create_legend(wks,11,lst,lgres)

;do i=0,10,1
lgres@lgLineColors = colors;(i)
;lgres@lgLineLabelFontColor = colors
;lgres@lgMonoLineLabelFontColor = True
;gsn_legend_ndc(wks,1,labels(i),xpos(i),ypos(i),lgres)
gsn_legend_ndc(wks,11,labels,xpos,ypos,lgres)
;end do

;draw(legend)
frame(wks)

print(" ")

;-----------------------------------------------------------------------------

if (jpeg) then
 system("ps2epsi "+pltName+".ps "+pltName+".esp")
 system("convert -density 500 "+pltName+".eps "+pltName+".jpg")
 ;system("convert -density 1000 "+pltName+".ps "+pltName+".jpg")
 ;system("convert -density 1000 -geometry 1200x1600 "+pltName+".ps "+pltName+".jpg")
end if

end

