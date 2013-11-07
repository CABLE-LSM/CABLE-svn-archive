;==============================================================================
; Lauren Stevens 8 July 2011
;==============================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/home/cmar/ste69f/umplot/plot_mmdc.ncl"
;==============================================================================

;procedure hay_ts()

begin

FlxName  = "hay"
LongName = "hay"
numYr    = systemfunc("echo $YR")
mossf    = systemfunc("echo $moshay")
cabsf    = systemfunc("echo $cabhay")
jpg = systemfunc("echo $CHFMT")
if(jpg.eq."y") then
jpeg = True
else
jpeg = False
end if

tsDir  = "/home/cmar/ste69f/umplot/obs/"
home   = "/home/cmar/ste69f/"

system("dmget " + tsDir + "hayobsjuly.nc")
flux1  = addfile(tsDir + "hayobsjuly.nc" ,"r")
;flux2  = addfile(tsDir + LongName +"/.nc" ,"r")
mosfil = addfile(mossf,"r")
cabfil = addfile(cabsf,"r")
;igbp   = addfile(home  + "MMDC_hay_roll.nc","r")

;==============================================================================

sensht1 = flux1 ->h
latht1  = flux1 ->lh
;rad1    = flux1 ->Rn
temp1   = flux1 ->ta
;precip1 = flux1 ->PPT
;---------------------
;qh      = igbp  ->Qh
;qle     = igbp  ->Qle
;tair    = igbp  ->Tair
;rainf   = igbp  ->Rainf
;----------------------
;csh     = cabfil ->hfss
;msh     = mosfil ->hfss
;clh     = cabfil ->hfls
;mlh     = mosfil ->hfls
;;crad    = cabfil ->field3333
;;mrad    = mosfil ->field3333
;cta     = cabfil ->tas
;mta     = mosfil ->tas
;;cpr     = cabfil ->field5226
;;mpr     = mosfil ->field5226

if (isfilevar(cabfil,"hfss")) then
csh     = cabfil ->hfss
else
csh     = cabfil ->sh
end if
if (isfilevar(mosfil,"hfss")) then
msh     = mosfil ->hfss
else
msh     = mosfil ->sh
end if
if (isfilevar(cabfil,"hfls")) then
clh     = cabfil ->hfls
else
clh     = cabfil ->lh
end if
if (isfilevar(mosfil,"hfls")) then
mlh     = mosfil ->hfls
else
mlh     = mosfil ->lh
end if
;if (isfilevar(cabfil,"field3333")) then
;crad    = cabfil ->field3333
;else
;crad    = cabfil ->field202
;end if
if (isfilevar(mosfil,"field3333")) then
mrad    = mosfil ->field3333
else
 mrad    = mosfil ->field202
 mswr    = mosfil ->solar
 mlwr    = mosfil ->longwave

 dNames  = getvardims(mrad)
 dName2  = getvardims(mswr)
 dims    = dimsizes(mrad)

 mrn1    = mswr+mlwr
 copy_VarMeta(mswr,mrn1)
 mrn1&$dName2(1)$ = ispan(0,47,6)
 mrn2    = linint1_n_Wrap(mrn1&$dName2(1)$,mrn1,True,mrad&$dNames(1)$,0,1)

 ; Roll 4ts for Hyytiala
 mrnt    = new((/dims(0),dims(1)/),typeof(mswr))
 mrnt(:,0:6)  = mrn2(:,41:47)
 mrnt(:,7:47) = mrn2(:,0:40)
end if
if (isfilevar(cabfil,"tas")) then
cta     = cabfil ->tas
else
cta     = cabfil ->temp
end if
if (isfilevar(mosfil,"tas")) then
mta     = mosfil ->tas
else
mta     = mosfil ->temp
end if
if (isfilevar(cabfil,"field5226")) then
cpr     = cabfil ->field5226
else
cpr     = cabfil ->tot_precip
end if
;if (isfilevar(mosfil,"field5226")) then
;mpr     = mosfil ->field5226
;else
;mpr     = mosfil ->tot_precip
;end if

;==============================================================================

dNam     = getvardims(msh)
tstp     = msh&$dNam(1)$
tstep    = max(tstp)+1
timestep = doubletoint(tstep)

sensht1@_FillValue = -9999
sensht1!0 = "timestep"
;sens_all           = (sensht1+sensht2+sensht3)/3
;sens_2             = (sensht2+sensht3)/2
;copy_VarMeta(sensht1,sens_all)
;copy_VarMeta(sensht1,sens_2)
latht1@_FillValue  = -9999
latht1!0 = "timestep"
;lh_all             = (latht1+latht2+latht3)/3
;lh_2               = (latht2+latht3)/2
;copy_VarMeta(latht1,lh_all)
;copy_VarMeta(latht1,lh_2)
;rad1@_FillValue    = -9999
;rad2@_FillValue    = -9999
;rad3@_FillValue    = -9999
;rad_all            = (rad1+rad2+rad3)/3
;;rad_2              = rad3 ;(rad2+rad3)/2
;copy_VarMeta(rad1,rad_all)
;;copy_VarMeta(rad1,rad_2)
temp1@_FillValue   = -9999
temp1!0 = "timestep"
;temp_all           = (temp1+temp2+temp3)/3 +273.16
;temp_2             = (temp2+temp3)/2 +273.16
;copy_VarMeta(temp1,temp_all)
;copy_VarMeta(temp1,temp_2)
;precip1@_FillValue = -9999
;precip2@_FillValue = -9999
;precip3@_FillValue = -9999
;pr_all             = (precip1+precip2+precip3)/3 *48
;;pr_2               = (precip2+precip3)/2 *48
;copy_VarMeta(precip1,pr_all)
;;copy_VarMeta(precip1,pr_2)

;qh_ta    = dim_avg_n_Wrap(qh,0)
;sens     = dim_avg_n_Wrap(sens_all,0)
;sens2    = dim_avg_n_Wrap(sens_2,0)
;qh_ta!0  = "month"
;qh_ta!1  = "timestep"
;sens!0   = "month"
;sens!1   = "timestep"
;sens2!0  = "month"
;sens2!1  = "timestep"
msh!0    = "month"
msh!1    = "timestep"
csh!0    = "month"
csh!1    = "timestep"

;qle_ta   = dim_avg_n_Wrap(qle,0)
;lats     = dim_avg_n_Wrap(lh_all,0)
;lats2    = dim_avg_n_Wrap(lh_2,0)
;qle_ta!0 = "month"
;qle_ta!1 = "timestep"
;lats!0   = "month"
;lats!1   = "timestep"
;lats2!0  = "month"
;lats2!1  = "timestep"
mlh!0    = "month"
mlh!1    = "timestep"
clh!0    = "month"
clh!1    = "timestep"

;tair_ta   = dim_avg_n_Wrap(tair,0)
;temp      = dim_avg_n_Wrap(temp_all,0)
;tem2      = dim_avg_n_Wrap(temp_2,0)
;tair_ta!0 = "month"
;tair_ta!1 = "timestep"
;temp!0    = "month"
;temp!1    = "timestep"
;tem2!0    = "month"
;tem2!1    = "timestep"
mta!0     = "month"
mta!1     = "timestep"
cta!0     = "month"
cta!1     = "timestep"

;rainf     = rainf*24*60*60
;mpr       = mpr *48
;cpr       = cpr *48
;rf_ta     = dim_avg_n_Wrap(rainf,0)
;prec      = dim_avg_n_Wrap(pr_all,0)
;;prec2     = dim_avg_n_Wrap(pr_2,0)
;rf_ta!0   = "month"
;rf_ta!1   = "timestep"
;prec!0    = "month"
;prec!1    = "timestep"
;;prec2!0   = "month"
;;prec2!1   = "timestep"
;mpr!0     = "month"
;mpr!1     = "timestep"
;cpr!0     = "month"
;cpr!1     = "timestep"

;;rad_ta    = dim_avg_n_Wrap(,0)
;rnet      = dim_avg_n_Wrap(rad_all,0)
;;rnet2     = dim_avg_n_Wrap(rad_2,0)
;;rad_ta!0  = "month"
;;rad_ta!1  = "timestep"
;rnet!0    = "month"
;rnet!1    = "timestep"
;;rnet2!0   = "month"
;;rnet2!1   = "timestep"
;mrad!0    = "month"
;mrad!1    = "timestep"
;crad!0    = "month"
;crad!1    = "timestep"

;----------------------------------------------------------------------------

y  = fspan(0,23.5,timestep)
y2 = ispan(0,23,1)

sh  = new(timestep,double)
lh  = new(timestep,double)
tmp = new(timestep,double)
do i=0,23
  sh(2*i+1)  = sensht1(i) 
  sh(2*i)  = sensht1(i)
  lh(2*i+1)  = latht1(i)
  lh(2*i)  = latht1(i)
  tmp(2*i+1) = temp1(i)
  tmp(2*i) = temp1(i)
end do
sh@_FillValue = -9999
lh@_FillValue = -9999
tmp@_FillValue = -9999
;print(tmp)

;----------------------------------------------------------------------------

; July
data      = new((/3,timestep/),double)
;data(0,:) = -9999 ;qh_ta(6,:)            ; igbp
data(0,:) = msh(6,:)              ; moses
data(1,:) = csh(6,:)              ; cable
data(2,:) = sh(:)             ; fluxnet

data1      = new((/3,timestep/),double)
;data1(0,:) = -9999 ;qle_ta(6,:)           ; igbp
data1(0,:) = mlh(6,:)              ; moses
data1(1,:) = clh(6,:)              ; cable
data1(2,:) = lh(:)             ; fluxnet

data7      = new((/3,timestep/),double)
;data7(0,:) = -9999 ;tair_ta(6,:)          ; igbp
data7(0,:) = mta(6,:)              ; moses
data7(1,:) = cta(6,:)              ; cable
data7(2,:) = tmp(:)             ; fluxnet

;data8      = new((/3,timestep/),double)
;;data8(0,:) = -9999 ;rad_ta(6,:)           ; igbp
;data8(0,:) = mrad(6,:)             ; moses
;data8(1,:) = crad(6,:)             ; cable
;data8(2,:) = rnet(6,:)             ; fluxnet

;data9      = new((/3,timestep/),double)
;;data9(0,:) = rf_ta(6,:)            ; igbp
;data9(0,:) = mpr(6,:)              ; moses
;data9(1,:) = cpr(6,:)              ; cable
;data9(2,:) = prec(6,:)             ; fluxnet

; Jan
;data2      = new((/3,timestep/),double)
;;data2(0,:) = qh_ta(0,:)            ; igbp
;data2(0,:) = msh(0,:)              ; moses
;data2(1,:) = csh(0,:)              ; cable
;data2(2,:) = sens(0,:)            ; fluxnet
;
;data3      = new((/3,timestep/),double)
;;data3(0,:) = qle_ta(0,:)           ; igbp
;data3(0,:) = mlh(0,:)              ; moses
;data3(1,:) = clh(0,:)              ; cable
;data3(2,:) = lats(0,:)            ; fluxnet
;
;data4      = new((/3,timestep/),double)
;;data4(0,:) = tair_ta(0,:)          ; igbp
;data4(0,:) = mta(0,:)              ; moses
;data4(1,:) = cta(0,:)              ; cable
;data4(2,:) = temp(0,:)             ; fluxnet
;
;data5      = new((/3,timestep/),double)
;;data5(0,:) = -9999 ;rad_ta(0,:)           ; igbp
;data5(0,:) = mrad(0,:)             ; moses
;data5(1,:) = crad(0,:)             ; cable
;data5(2,:) = rnet(0,:)            ; fluxnet
;
;data6      = new((/3,timestep/),double)
;;data6(0,:) = rf_ta(0,:)            ; igbp
;data6(0,:) = mpr(0,:)              ; moses
;data6(1,:) = cpr(0,:)              ; cable
;data6(2,:) = prec(0,:)            ; fluxnet

data0      = new((/3,timestep/),double)
data0      = 0 ;-9999

;Plot-------------------------------------------------------------------------

;ts_flux(LongName,data,data1,data7,data0)
;ts_rad(LongName,data5,data4,data8,data7)
;ts_prec(LongName,data6,data9)

pltName                    = "mmdc_hay_flux_"+numYr+"yrs"
wks                        = gsn_open_wks("ps","mmdc_hay_flux_"+numYr+"yrs")
;wks2                       = gsn_open_wks("ps","_ts_rad")
;wks3                       = gsn_open_wks("ps","_ts_prec")
res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
res@xyLineColors           = (/"blue","red","black"/)
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2/)
;res@xyLineThicknesses      = (/2,2,2,2/)
res@xyDashPattern          = 0.
res@tiYAxisString          = ""
res@tiXAxisString          = "Hours"
res@tiMainString           = "July Sensible Heat Flux (Wm-2)"
res@pmLegendDisplayMode    = "Always"            ; turn on legend
res@pmLegendSide           = "Top"               ; Change location of 
res@pmLegendParallelPosF   = .15                 ; move units right
res@pmLegendOrthogonalPosF = -0.3                ; move units down
res@pmLegendWidthF         = 0.09                ; Change width and
res@pmLegendHeightF        = 0.11                ; height of legend.
res@lgPerimOn              = False               ; turn off box around
res@lgLabelFontHeightF     = .018                ; label font height
res@xyExplicitLegendLabels = (/"moses","cable","fluxnet"/)         ; create explicit labels
;res@xyExplicitLegendLabels = (/"igbp","moses","cable","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

print("Plotting Hay MMDC with Fluxnet")
print("")

plot                       = new(4,graphic)
plot(0)                    = gsn_csm_xy(wks,y,data,res)
res@tiYAxisString          = ""
res@tiMainString           = "July Temperature (deg K)"
plot(1)                    = gsn_csm_xy(wks,y,data7,res) 
res@tiYAxisString          = ""
res@tiMainString           = "July Latent Heat Flux (Wm-2)"
plot(2)                    = gsn_csm_xy(wks,y,data1,res)
res@tiMainString           = ""
res@tiYAxisString          = "" 
plot(3)                    = gsn_csm_xy(wks,y,data0,res)

;radplt                     = new(4,graphic)
;res@tiYAxisString          = "Net Radiation (Wm-2)"
;res@tiMainString           = "January"
;radplt(0)                  = gsn_csm_xy(wks2,y,data5,res)
;res@tiYAxisString          = ""
;res@tiMainString           = "July"
;radplt(1)                  = gsn_csm_xy(wks2,y,data8,res)
;res@tiYAxisString          = "Temperature (deg K)"
;res@tiMainString           = ""
;radplt(2)                  = gsn_csm_xy(wks2,y,data4,res)
;res@tiYAxisString          = ""
;res@tiMainString           = ""
;radplt(3)                  = gsn_csm_xy(wks2,y,data7,res)
;
;prplot                     = new(2,graphic)
;res@tiYAxisString          = "Precipitation ()"
;res@tiMainString           = "January"
;prplot(0)                  = gsn_csm_xy(wks3,y,data6,res)
;res@tiYAxisString          = ""
;res@tiMainString           = "July"
;prplot(1)                  = gsn_csm_xy(wks3,y,data9,res)

; -----------------------------------------------------------------------------

pres                       = True
pres@gsnMaximize           = True
;pres@gsnOrientation        = "landscape"
;pres@PanelRowSpec          = True
gsn_panel(wks,plot,(/2,2/),pres)
;gsn_panel(wks2,radplt,(/2,2/),pres)
;gsn_panel(wks3,prplot,(/1,2/),pres)

if(jpeg) then
system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("convert -density 500 "+pltName+".eps "+pltName+".jpg")
;system("convert -density 1000 "+pltName+".ps "+pltName+".jpg")
;system("convert -density 1000 -geometry 1200x1600 "+pltName+".ps "+pltName+".jpg")
end if

end

