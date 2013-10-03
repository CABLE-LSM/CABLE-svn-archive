;==============================================================================
; Lauren Stevens 19 July 2011
;==============================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "/home/cmar/ste69f/umplot/plot_mmdc.ncl"
;==============================================================================

;procedure mmdc_flux(LongName:string,mosDir:string,cabDir:string)

begin

;LongName = "hyytiala"
numYr    = systemfunc("echo $YR")
filDir   = systemfunc("echo $DIRW")
jpg = systemfunc("echo $JPEG")
if(jpg.eq."y") then
jpeg = True
else
jpeg = False
end if

tsDir    = "/home/cmar/ste69f/umplot/obs/"
home     = "/home/cmar/ste69f/"

;system("dmget " + tsDir + LongName +"/mmdc*_roll.nc" )

;obs
flux_ma  = addfile(tsDir + "MA96_lu_u0_mm_flx.nc"        ,"r")
met_ma   = addfile(tsDir + "MA96_mm_met.nc"              ,"r")
flux_tb  = addfile(tsDir + "AnnCycles_tumbarumba.nc" ,"r")
; model
mosfil   = addfile(filDir+"/AnnualCycle_"+numYr+"yrs.nc"        ,"r")
;igbp     = addfile(home  + "igbp/"+ LongName +"/MMDC_"+ LongName +"_roll.nc","r")

;OBS===========================================================================

fils   = systemfunc("ls "+tsDir+"*-*_lu_u0_mm_flx.nc")
f      = addfiles(fils,"r")
ListSetType(f,"join")
sens   = f[:]->H;(0:11)
lat    = f[:]->LH;(0:11)

fils2  = systemfunc("ls "+tsDir+"*-*_mm_met.nc")
f2     = addfiles(fils2,"r")
ListSetType(f2,"join")
radnet = f2[:]->Rn;(0:11)
temp2m = f2[:]->Ta;(0:11)
preci  = f2[:]->PPT;(0:11)

sens_ma = flux_ma ->H
lat_ma  = flux_ma ->LH
rad_ma  = met_ma  ->Rn
temp_ma = met_ma  ->Ta
prec_ma = met_ma  ->PPT
sens_tb = flux_tb ->qh_8yr_mean    ;Qh    ;H
lat_tb  = flux_tb ->qle_8yr_mean   ;Qle   ;LH
rad_tb  = flux_tb ->rnet_8yr_mean  ;Rnet  ;Rn
temp_tb = flux_tb ->tair_8yr_mean  ;Tair  ;Ta
prec_tb = flux_tb ->rainf_8yr_mean ;Rainf ;PPT

;OBS===========================================================================

sensht            = new((/11,12/),double)
latht             = new((/11,12/),double)
rad               = new((/11,12/),double)
temp              = new((/11,12/),double)
precip            = new((/11,12/),double)
sensht            = -9999
latht             = -9999
rad               = -9999
temp              = -9999
precip            = -9999
sensht@_FillValue = -9999
latht@_FillValue  = -9999
rad@_FillValue    = -9999
temp@_FillValue   = -9999
precip@_FillValue = -9999

do k=0,8,1
sensht(k,:)       = sens(k,0:11)*1000000/(31*24*60*60)   ;H
latht(k,:)        = lat(k,0:11)*1000000/(31*24*60*60)    ;LH
rad(k,:)          = radnet(k,0:11)*1000000/(31*24*60*60) ;Rn
temp(k,:)         = temp2m(k,0:11) ;Ta
precip(k,:)       = preci(k,0:11)  ;PPT
end do
sensht(9,:)       = sens_ma(0:11)*1000000/(31*24*60*60) ;H
latht(9,:)        = lat_ma(0:11)*1000000/(31*24*60*60)  ;LH
rad(9,:)          = rad_ma(0:11)*1000000/(31*24*60*60)  ;Rn
temp(9,:)         = temp_ma(0:11) ;Ta
precip(9,:)       = prec_ma(0:11) ;PPT
sensht(10,:)      = sens_tb(0:11) ;H
latht(10,:)       = lat_tb(0:11)  ;LH
rad(10,:)         = rad_tb(0:11)  ;Rn
temp(10,:)        = temp_tb(0:11) ;Ta
precip(10,:)      = prec_tb(0:11) ;PPT

;---------------------
;qh      = igbp  ->Qh
;qle     = igbp  ->Qle
;tair    = igbp  ->Tair
;rainf   = igbp  ->Rainf
;----------------------

;MODEL==========================================

;csh     = cabfil ->hfss
msh     = mosfil ->hfss
;clh     = cabfil ->hfls
mlh     = mosfil ->hfls
;crad    = cabfil ->field3333
mrad    = mosfil ->field3333
;cta     = cabfil ->tas
mta     = mosfil ->tas
;cpr     = cabfil ->field5226
mpr     = mosfil ->field5226

dNames  = getfilevardims(mosfil,"hfss")
dNlh    = getfilevardims(mosfil,"hfls")
dNrad   = getfilevardims(mosfil,"field3333")
;print(dNames)

mrad_ro = new((/11,12/),double)
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

;sens_all           = (sensht1+sensht2+sensht3)/3
;sens_2             = (sensht2+sensht3)/2
;copy_VarMeta(sensht1,sens_all)
;copy_VarMeta(sensht1,sens_2)
;lh_all             = (latht1+latht2+latht3)/3
;lh_2               = (latht2+latht3)/2
;copy_VarMeta(latht1,lh_all)
;copy_VarMeta(latht1,lh_2)

sh_ANN            = new((/11,12/),double)
sh_ANN            = msh_ro ;(axis_4|0:10,axis_1|:)
sh_ANN            = sh_ANN - sensht
lh_ANN            = new((/11,12/),double)
lh_ANN            = mlh_ro ;(axis_4|0:10,axis_1|:)
lh_ANN            = lh_ANN - latht
rad_ANN           = new((/11,12/),double)
rad_ANN           = mrad_ro ;(axis_4|0:10,axis_1|:)
rad_ANN           = rad_ANN - rad 

;rad_2              = rad3 ;(rad2+rad3)/2
;copy_VarMeta(rad1,rad_all)
;copy_VarMeta(rad1,rad_2)
;temp_all           = (temp1+temp2+temp3)/3 +273.16
;temp_2             = (temp2+temp3)/2 +273.16
;copy_VarMeta(temp1,temp_all)
;copy_VarMeta(temp1,temp_2)
;pr_all             = (precip1+precip2+precip3)/3 *48
;pr_2               = (precip2+precip3)/2 *48
;copy_VarMeta(precip1,pr_all)
;copy_VarMeta(precip1,pr_2)

;qh_ta    = dim_avg_n_Wrap(qh,0)
;sens     = dim_avg_n_Wrap(sens_all,0)
;sens2    = dim_avg_n_Wrap(sens_2,0)
;qh_ta!0  = "month"
;qh_ta!1  = "timestep"
;sens!0   = "month"
;sens!1   = "timestep"
;sens2!0  = "month"
;sens2!1  = "timestep"
;msh!0    = "month"
;msh!1    = "timestep"
;csh!0    = "month"
;csh!1    = "timestep"

;qle_ta   = dim_avg_n_Wrap(qle,0)
;lats     = dim_avg_n_Wrap(lh_all,0)
;lats2    = dim_avg_n_Wrap(lh_2,0)
;qle_ta!0 = "month"
;qle_ta!1 = "timestep"
;lats!0   = "month"
;lats!1   = "timestep"
;lats2!0  = "month"
;lats2!1  = "timestep"
;mlh!0    = "month"
;mlh!1    = "timestep"
;clh!0    = "month"
;clh!1    = "timestep"

;tair_ta   = dim_avg_n_Wrap(tair,0)
;temp      = dim_avg_n_Wrap(temp_all,0)
;tem2      = dim_avg_n_Wrap(temp_2,0)
;tair_ta!0 = "month"
;tair_ta!1 = "timestep"
;temp!0    = "month"
;temp!1    = "timestep"
;tem2!0    = "month"
;tem2!1    = "timestep"
;mta!0     = "month"
;mta!1     = "timestep"
;cta!0     = "month"
;cta!1     = "timestep"

;rainf     = rainf*24*60*60
;mpr       = mpr *48
;cpr       = cpr *48
;rf_ta     = dim_avg_n_Wrap(rainf,0)
;prec      = dim_avg_n_Wrap(pr_all,0)
;prec2     = dim_avg_n_Wrap(pr_2,0)
;rf_ta!0   = "month"
;rf_ta!1   = "timestep"
;prec!0    = "month"
;prec!1    = "timestep"
;prec2!0   = "month"
;prec2!1   = "timestep"
;mpr!0     = "month"
;mpr!1     = "timestep"
;cpr!0     = "month"
;cpr!1     = "timestep"

;rad_ta    = dim_avg_n_Wrap(,0)
;rnet      = dim_avg_n_Wrap(rad_all,0)
;rnet2     = dim_avg_n_Wrap(rad_2,0)
;rad_ta!0  = "month"
;rad_ta!1  = "timestep"
;rnet!0    = "month"
;rnet!1    = "timestep"
;rnet2!0   = "month"
;rnet2!1   = "timestep"
;mrad!0    = "month"
;mrad!1    = "timestep"
;crad!0    = "month"
;crad!1    = "timestep"

;----------------------------------------------------------------------------

;y = fspan(0,23.5,48)
y = ispan(1,12,1)

; July
;data      = new((/4,48/),double)
;data(0,:) = qh_ta(6,:)            ; igbp
;data(1,:) = msh(6,:)              ; moses
;data(2,:) = csh(6,:)              ; cable
;data(3,:) = sens(6,:)             ; fluxnet

;data1      = new((/4,48/),double)
;data1(0,:) = qle_ta(6,:)           ; igbp
;data1(1,:) = mlh(6,:)              ; moses
;data1(2,:) = clh(6,:)              ; cable
;data1(3,:) = lats(6,:)             ; fluxnet

;data7      = new((/4,48/),double)
;data7(0,:) = tair_ta(6,:)          ; igbp
;data7(1,:) = mta(6,:)              ; moses
;data7(2,:) = cta(6,:)              ; cable
;data7(3,:) = temp(6,:)             ; fluxnet

;data8      = new((/4,48/),double)
;data8(0,:) = -9999 ;rad_ta(6,:)           ; igbp
;data8(1,:) = mrad(6,:)             ; moses
;data8(2,:) = crad(6,:)             ; cable
;data8(3,:) = rnet2(6,:)             ; fluxnet

;data9      = new((/4,48/),double)
;data9(0,:) = rf_ta(6,:)            ; igbp
;data9(1,:) = mpr(6,:)              ; moses
;data9(2,:) = cpr(6,:)              ; cable
;data9(3,:) = prec(6,:)             ; fluxnet
;
; Jan
;data2      = new((/4,48/),double)
;data2(0,:) = qh_ta(0,:)            ; igbp
;data2(1,:) = msh(0,:)              ; moses
;data2(2,:) = csh(0,:)              ; cable
;data2(3,:) = sens2(0,:)            ; fluxnet
;
;data3      = new((/4,48/),double)
;data3(0,:) = qle_ta(0,:)           ; igbp
;data3(1,:) = mlh(0,:)              ; moses
;data3(2,:) = clh(0,:)              ; cable
;data3(3,:) = lats2(0,:)            ; fluxnet
;
;data4      = new((/4,48/),double)
;data4(0,:) = tair_ta(0,:)          ; igbp
;data4(1,:) = mta(0,:)              ; moses
;data4(2,:) = cta(0,:)              ; cable
;data4(3,:) = tem2(0,:)             ; fluxnet
;
;data5      = new((/4,48/),double)
;data5(0,:) = -9999 ;rad_ta(0,:)           ; igbp
;data5(1,:) = mrad(0,:)             ; moses
;data5(2,:) = crad(0,:)             ; cable
;data5(3,:) = rnet2(0,:)            ; fluxnet
;
;data6      = new((/4,48/),double)
;data6(0,:) = rf_ta(0,:)            ; igbp
;data6(1,:) = mpr(0,:)              ; moses
;data6(2,:) = cpr(0,:)              ; cable
;data6(3,:) = prec2(0,:)            ; fluxnet
;
;Plot-------------------------------------------------------------------------

print("Plotting MM-AnnualC Biases for All Sites")

;ts_flux(LongName,data2,data3,data,data1)
;ts_rad(LongName,data5,data4,data8,data7)
;ts_prec(LongName,data6,data9)

pltName = "annualc_"+numYr+"yrs"
wks  = gsn_open_wks("ps","annualc_"+numYr+"yrs")
res = True
res@gsnFrame = False
res@gsnDraw  = False
res@xyLineColors           = (/"black","deepskyblue1","red","green","mediumorchid1","cyan","magenta","firebrick","orange","SlateBlue","navy","SkyBlue"/)
res@xyDashPatterns         = (/0,0,0,0,0,0,0,0,0,0,1,2/)
res@xyLineThicknesses      = (/2,2,2,2,2,2,2,2,2,2,2,2/)
res@gsnYRefLine            = 0.
;res@gsnLeftString          = ""
;res@tiYAxisString          = "W/m2"

res@tiYAxisFontHeightF = 0.02
res@tiYAxisOffsetXF    = 0.005
res@tiYAxisString      = "W/m2"
res@tiXAxisString      = "months"
;if (rid.eq."saaca") then
;res@gsnLeftString      = "Annual - A1.0"
;else
res@gsnLeftString      = "Annual" ; - A1.1"
;end if
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
;print((/floor(min(mnf)),ceil(max(mxf))/))
res@trYMinF                = floor(min(mnf));-90
res@trYMaxF                = ceil(max(mxf)) ;120

;res@pmLegendDisplayMode    = "Always"            ; turn on legend
;res@pmLegendSide           = "Top"            ; Change location of
;res@pmLegendParallelPosF   = .25                 ; move units right
;res@pmLegendOrthogonalPosF = -0.2                ; move units down
;res@pmLegendWidthF         = 0.09                ; Change width and
;res@pmLegendHeightF        = 0.21                ; height of legend.
;res@lgPerimOn              = False               ; turn off box around
;res@lgLabelFontHeightF     = .008                ; label font height
;res@xyExplicitLegendLabels = (/"bondville","harvardforest","hyytiala","littlewashita","loobos","nsaboreas","tharandt","vielsalm","walkerbranch","manaus","tumbarumba"/)         ; create explicit labels

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
 res@xyExplicitLegendLabels = (/"bond","harv","hyyt","loob","lwas","nbor","thar","","wbra","","","hay"/)
plot(2) = gsn_csm_xy(wks,y,lh_ANN,res)

; bv hv hy lw lo nb th vi wb ma tb
; Legend ---------------------------------------------------------------------

;colors = (/"black","blue","red","green","purple","cyan","magenta","firebrick","orange","SlateBlue","navy","SkyBlue"/)
;labels = (/"bondville","harvardforest","hyytiala","littlewashita","loobos","nsaboreas","tharandt","vielsalm","walkerbranch","manaus","tumbarumba"/)

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

;legend = create "Legend" legendClass wks
;"vpXF"                     : 0.26
;"vpYF"                     : 0.37
;"vpWidthF"                 : 0.5                     ; width
;"vpHeightF"                : 0.3                     ; height
;"lgPerimOn"                : False                   ; no perimeter
;"lgItemCount"              : 11                      ; how many
;"lgLineLabelStrings"       : (/"bondville","harvardforest","hyytiala","littlewashita","loobos","nsaboreas","tharandt","vielsalm","walkerbranch","manaus","tumbarumba"/)  ; labels
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
;labels = (/"bondville","harvardforest","hyytiala","littlewashita","loobos","nsaboreas","tharandt","vielsalm","walkerbranch","manaus","tumbarumba"/)
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

 if (jpeg) then
 system("convert -density 1000 "+pltName+".ps "+pltName+".jpg")
 ;system("convert -density 1000 -geometry 1200x1600 "+pltName+".ps "+pltName+".jpg")
 end if

end

