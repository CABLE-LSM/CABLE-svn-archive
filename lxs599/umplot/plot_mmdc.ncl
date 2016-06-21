;====================================================================
; Lauren Stevens 8 July 2011
;====================================================================
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;====================================================================

; CONTAINS: ts_flux, ts_rad, ts_prec, ts_prec4, ts_swlw

procedure ts_flux(FlxName:string,ShJan,LhJan,ShJul,LhJul,jpeg) ; :numeric

begin

;  pltName = Dir + "/" + varName + "_" + Yr + "yrs"
;  wks = gsn_open_wks(pltType,pltName)

numYr  = systemfunc("echo $YR")
tstep    = systemfunc("echo $TSTEP")
y = fspan(0,23.5,tstep)

dims = dimsizes(ShJan)

;-----------------------------------------------------------------------------

pltName                    = "mmdc_"+FlxName+"_flux_"+numYr+"yrs"
wks                        = gsn_open_wks("ps","mmdc_"+FlxName+"_flux_"+numYr+"yrs")

res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
if (dims(0).eq.5) then
res@xyLineColors           = (/"darkorange","forestgreen","blue","red","black"/)
else
res@xyLineColors           = (/"blue","red","black"/)
end if
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2,2/)
res@xyDashPattern          = 0.
res@tiYAxisString          = "Sensible Heat Flux (Wm~S~-2~N~)"
res@tiXAxisString          = "Hours"
res@gsnRightString         = "January"
res@gsnLeftString          = ""+FlxName+""
;res@tiMainString           = "January"
res@pmLegendDisplayMode    = "Always"            ; turn on legend
res@pmLegendSide           = "Top"               ; Change location of 
res@pmLegendParallelPosF   = .15                 ; move units right
res@pmLegendOrthogonalPosF = -0.3                ; move units down
res@pmLegendWidthF         = 0.09                ; Change width and
res@pmLegendHeightF        = 0.11                ; height of legend.
res@lgPerimOn              = False               ; turn off box around
res@lgLabelFontHeightF     = .016                ; label font height
if (dims(0).eq.5) then
res@xyExplicitLegendLabels = (/" a10"," a11_r0"," a11_r4"," a11_r5"," fluxnet"/)         ; create explicit labels
else
res@xyExplicitLegendLabels = (/"A1.0","A1.1","Obs"/)         ; create explicit labels
end if
;res@xyExplicitLegendLabels = (/"igbp","A1.0","A1.1","fluxnet"/)         ; create explicit labels

mxjan = (/max(ShJan),max(LhJan)/)
mnjan = (/min(ShJan),min(LhJan)/)
mxjul = (/max(ShJul),max(LhJul)/)
mnjul = (/min(ShJul),min(LhJul)/)
;if (FlxName .eq. "Boreas") then
; res@trYMinF = floor(min(mnjan))-10
; res@trYMaxF =  ceil(max(mxjan))+10
;end if
;if (FlxName .eq. "Boreas") then
; res@trYMinF = floor(min(mnjul))-10
; res@trYMaxF =  ceil(max(mxjul))+10
;end if

;------------------------------------------------------------------------------

plot                       = new(4,graphic)
if (FlxName .eq. "Boreas") then
 res@trYMinF = -20 ;floor(min(mnjan))-10
 res@trYMaxF =  60 ; ceil(max(mxjan))+10
end if
if (FlxName .eq. "Loobos") then
 res@trYMinF = -40 ;floor(min(mnjan))-10
 res@trYMaxF =  50 ; ceil(max(mxjan))+10
end if
if (FlxName .eq. "Mopane") then
 res@trYMinF = -50 ;floor(min(mnjan))-10
 res@trYMaxF = 450 ; ceil(max(mxjan))+10
end if
plot(0)                    = gsn_csm_xy(wks,y,ShJan,res) ; data2
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString        = "July"
;res@gsnLeftString          = ""
;res@tiMainString           = "July"
if (FlxName .eq. "Boreas") then
 res@trYMinF = -50 ;floor(min(mnjul))-10
 res@trYMaxF = 300 ; ceil(max(mxjul))+10
end if
if (FlxName .eq. "Loobos") then
 res@trYMinF = -40 ;floor(min(mnjul))-10
 res@trYMaxF = 280 ; ceil(max(mxjul))+10
end if
if (FlxName .eq. "Mopane") then
 res@trYMinF = -50 ;floor(min(mnjul))-10
 res@trYMaxF = 300 ; ceil(max(mxjul))+10
end if
plot(1)                    = gsn_csm_xy(wks,y,ShJul,res) ; data
res@tiYAxisOn              = True
res@tiYAxisString          = "Latent Heat Flux (Wm~S~-2~N~)"
res@gsnRightString        = "January"
;res@tiMainString           = ""
if (FlxName .eq. "Boreas") then
 res@trYMinF = -20 ;floor(min(mnjan))-10
 res@trYMaxF =  60 ; ceil(max(mxjan))+10
end if
if (FlxName .eq. "Loobos") then
 res@trYMinF = -40 ;floor(min(mnjan))-10
 res@trYMaxF =  50 ; ceil(max(mxjan))+10
end if
if (FlxName .eq. "Mopane") then
 res@trYMinF = -50 ;floor(min(mnjan))-10
 res@trYMaxF = 450 ; ceil(max(mxjan))+10
end if
plot(2)                    = gsn_csm_xy(wks,y,LhJan,res) ; data3
res@gsnRightString        = "July"
;res@tiMainString           = ""
;res@tiYAxisString          = "" 
res@tiYAxisOn              = False
if (FlxName .eq. "Boreas") then
 res@trYMinF = -50 ;floor(min(mnjul))-10
 res@trYMaxF = 300 ; ceil(max(mxjul))+10
end if
if (FlxName .eq. "Loobos") then
 res@trYMinF = -40 ;floor(min(mnjul))-10
 res@trYMaxF = 280 ; ceil(max(mxjul))+10
end if
if (FlxName .eq. "Mopane") then
 res@trYMinF = -50 ;floor(min(mnjul))-10
 res@trYMaxF = 300 ; ceil(max(mxjul))+10
end if
plot(3)                    = gsn_csm_xy(wks,y,LhJul,res) ; data1

; -----------------------------------------------------------------------------

pres                       = True
pres@gsnMaximize           = True
;pres@gsnOrientation        = "landscape"
;pres@PanelRowSpec          = True
pres@gsnPanelYWhiteSpacePercent = 5
gsn_panel(wks,plot,(/2,2/),pres)

 if (jpeg) then
 system("convert -density 500 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 ;system("convert -density 500 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 end if


end

;===============================================================================

procedure ts_rad(FlxName:string,RnJan,TaJan,RnJul,TaJul,jpeg) ; :numeric

begin

numYr  = systemfunc("echo $YR")
y = fspan(0,23.5,tstep)
dims = dimsizes(RnJan)

;-----------------------------------------------------------------------------

pltName                    = "mmdc_"+FlxName+"_rad_"+numYr+"yrs"
wks2                       = gsn_open_wks("ps","mmdc_"+FlxName+"_rad_"+numYr+"yrs")

res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
if (dims(0).eq.5) then
res@xyLineColors           = (/"darkorange","forestgreen","blue","red","black"/)
else
res@xyLineColors           = (/"blue","red","black"/)
end if
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2,2/)
res@xyDashPattern          = 0.
res@tiYAxisString          = "Net Radiation (Wm~S~-2~N~)"
res@tiXAxisString          = "Hours"
res@gsnRightString         = "January"
res@gsnLeftString          = ""+FlxName+""
;res@tiMainString           = "January"
res@pmLegendDisplayMode    = "Always"            ; turn on legend
res@pmLegendSide           = "Top"               ; Change location of
res@pmLegendParallelPosF   = .15                 ; move units right
res@pmLegendOrthogonalPosF = -0.3                ; move units down
res@pmLegendWidthF         = 0.09                ; Change width and
res@pmLegendHeightF        = 0.11                ; height of legend.
res@lgPerimOn              = False               ; turn off box around
res@lgLabelFontHeightF     = .016                ; label font height
if (dims(0).eq.5) then
res@xyExplicitLegendLabels = (/" a10"," a11_r0"," a11_r4"," a11_r5"," fluxnet"/)         ; create explicit labels
else
res@xyExplicitLegendLabels = (/"A1.0","A1.1","Obs"/)         ; create explicit labels
end if
;res@xyExplicitLegendLabels = (/"igbp","A1.0","A1.1","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

radplt                     = new(4,graphic)
radplt(0)                  = gsn_csm_xy(wks2,y,RnJan,res) ; data5
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
;res@gsnLeftString          = ""
;res@tiMainString           = "July"
radplt(1)                  = gsn_csm_xy(wks2,y,RnJul,res) ; data8
res@tiYAxisOn              = True
res@tiYAxisString          = "Temperature (deg K)"
res@gsnRightString         = "January"
;res@tiMainString           = ""
radplt(2)                  = gsn_csm_xy(wks2,y,TaJan,res) ; data4
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString         = "July"
;res@tiMainString           = ""
radplt(3)                  = gsn_csm_xy(wks2,y,TaJul,res) ; data7

; -----------------------------------------------------------------------------

pres                       = True
pres@gsnMaximize           = True
;pres@gsnOrientation        = "landscape"
;pres@PanelRowSpec          = True
pres@gsnPanelYWhiteSpacePercent = 5
gsn_panel(wks2,radplt,(/2,2/),pres)

 if (jpeg) then
 system("convert -density 500 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 ;system("convert -density 500 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 end if


end

;===============================================================================

procedure ts_prec4(FlxName:string,PrJan,ClJan,PrJul,ClJul,jpeg) ; :numeric

begin

numYr  = systemfunc("echo $YR")
y = fspan(0,23.5,tstep)

dims = dimsizes(PrJan)

;-----------------------------------------------------------------------------

pltName                    = "mmdc_"+FlxName+"_prec_"+numYr+"yrs"
wks3                       = gsn_open_wks("ps","mmdc_"+FlxName+"_prec_"+numYr+"yrs")

res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
if (dims(0).eq.5) then
res@xyLineColors           = (/"darkorange","forestgreen","blue","red","black"/)
else
res@xyLineColors           = (/"blue","red","black"/)
end if
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2,2/)
res@xyDashPattern          = 0.
res@tiYAxisString          = "Precipitation (mm/day)"
res@tiXAxisString          = "Hours"
res@gsnRightString         = "January"
res@gsnLeftString          = ""+FlxName+""
;res@tiMainString           = "January"
res@pmLegendDisplayMode    = "Always"            ; turn on legend
res@pmLegendSide           = "Top"               ; Change location of
res@pmLegendParallelPosF   = .15                 ; move units right
res@pmLegendOrthogonalPosF = -0.3                ; move units down
res@pmLegendWidthF         = 0.09                ; Change width and
res@pmLegendHeightF        = 0.11                ; height of legend.
res@lgPerimOn              = False               ; turn off box around
res@lgLabelFontHeightF     = .016                ; label font height
if (dims(0).eq.5) then
res@xyExplicitLegendLabels = (/" a10"," a11_r0"," a11_r4"," a11_r5"," fluxnet"/)         ; create explicit labels
else
res@xyExplicitLegendLabels = (/"A1.0","A1.1","Obs"/)         ; create explicit labels
end if
;res@xyExplicitLegendLabels = (/"igbp","A1.0","A1.1","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

prplot                     = new(4,graphic)
if (FlxName .eq. "Boreas") then
 res@trYMinF = 0. ;floor(min(mnjul))-10
 res@trYMaxF = 1.5; ceil(max(mxjul))+10
end if
prplot(0)                  = gsn_csm_xy(wks3,y,PrJan,res) ; data6
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
;res@gsnLeftString          = ""
;res@tiMainString           = "July"
if (FlxName .eq. "Boreas") then
 res@trYMinF =  0. ;floor(min(mnjul))-10
 res@trYMaxF = 10. ; ceil(max(mxjul))+10
end if
prplot(1)                  = gsn_csm_xy(wks3,y,PrJul,res) ; data9
;res@tiYAxisString          = ""
res@tiYAxisOn              = True
res@tiYAxisString          = "Total Cloud Fraction"
res@gsnRightString           = "January"
;res@gsnLeftString            = ""
if (FlxName .eq. "Boreas") then
 res@trYMinF = 0.56 ;floor(min(mnjul))-10
 res@trYMaxF = 0.80 ; ceil(max(mxjul))+10
end if
prplot(2)                  = gsn_csm_xy(wks3,y,ClJan,res) ; data6
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
;res@gsnLeftString            = ""
if (FlxName .eq. "Boreas") then
 res@trYMinF = 0.30 ;floor(min(mnjul))-10
 res@trYMaxF = 1.00 ; ceil(max(mxjul))+10
end if
prplot(3)                  = gsn_csm_xy(wks3,y,ClJul,res) ; data9

; -----------------------------------------------------------------------------

pres                       = True
;pres@gsnMaximize           = True
;pres@gsnOrientation        = "landscape"
;pres@PanelRowSpec          = True
pres@gsnPanelYWhiteSpacePercent = 5
gsn_panel(wks3,prplot,(/2,2/),pres)

 if (jpeg) then
 system("ps2epsi "+pltName+".ps "+pltName+".eps")
 system("convert -density 500 -rotate 270 "+pltName+".eps "+pltName+".jpg")
 ;system("convert -density 500 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 ;system("convert -density 500 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 end if

end

;===============================================================================

procedure ts_prec(FlxName:string,PrJan,PrJul,jpeg) ; :numeric

begin

numYr  = systemfunc("echo $YR")
y = fspan(0,23.5,tstep)

dims = dimsizes(PrJan)

;-----------------------------------------------------------------------------

pltName                    = "mmdc_"+FlxName+"_prec_"+numYr+"yrs"
wks3                       = gsn_open_wks("ps","mmdc_"+FlxName+"_prec_"+numYr+"yrs")

res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
if (dims(0).eq.5) then
res@xyLineColors           = (/"darkorange","forestgreen","blue","red","black"/)
else
res@xyLineColors           = (/"blue","red","black"/)
end if
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2,2/)
res@xyDashPattern          = 0.
res@tiYAxisString          = "Precipitation (mm/day)"
res@tiXAxisString          = "Hours"
res@gsnRightString         = "January"
res@gsnLeftString          = ""+FlxName+""
;res@tiMainString           = "January"
res@pmLegendDisplayMode    = "Always"            ; turn on legend
res@pmLegendSide           = "Top"               ; Change location of
res@pmLegendParallelPosF   = .15                 ; move units right
res@pmLegendOrthogonalPosF = -0.3                ; move units down
res@pmLegendWidthF         = 0.09                ; Change width and
res@pmLegendHeightF        = 0.11                ; height of legend.
res@lgPerimOn              = False               ; turn off box around
res@lgLabelFontHeightF     = .016                ; label font height
if (dims(0).eq.5) then
res@xyExplicitLegendLabels = (/" a10"," a11_r0"," a11_r4"," a11_r5"," fluxnet"/)         ; create explicit labels
else
res@xyExplicitLegendLabels = (/"A1.0","A1.1","Obs"/)         ; create explicit labels
end if
;res@xyExplicitLegendLabels = (/"igbp","A1.0","A1.1","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

prplot                     = new(2,graphic)
prplot(0)                  = gsn_csm_xy(wks3,y,PrJan,res) ; data6
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
;res@gsnLeftString          = ""
;res@tiMainString           = "July"
prplot(1)                  = gsn_csm_xy(wks3,y,PrJul,res) ; data9

; -----------------------------------------------------------------------------

pres                       = True
pres@gsnMaximize           = True
;pres@gsnOrientation        = "landscape"
;pres@PanelRowSpec          = True
pres@gsnPanelYWhiteSpacePercent = 5
gsn_panel(wks3,prplot,(/1,2/),pres)

 if (jpeg) then
 system("ps2epsi "+pltName+".ps "+pltName+".eps")
 system("convert -density 500 -rotate 270 "+pltName+".eps "+pltName+".jpg")
 ;system("convert -density 1000 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 ;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 end if

end

;===============================================================================

procedure ts_swlw(FlxName:string,SwJan,LwJan,SwJul,LwJul,jpeg) ; :numeric

begin

numYr  = systemfunc("echo $YR")
y = fspan(0,23.5,tstep)

dims = dimsizes(SwJan)

;-----------------------------------------------------------------------------

pltName                    = "mmdc_"+FlxName+"_swlw_"+numYr+"yrs"
wks3                       = gsn_open_wks("ps","mmdc_"+FlxName+"_swlw_"+numYr+"yrs")

res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
if (dims(0).eq.5) then
res@xyLineColors           = (/"darkorange","forestgreen","blue","red","black"/)
else
res@xyLineColors           = (/"blue","red","black"/)
end if
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2,2/)
res@xyDashPattern          = 0.
res@tiYAxisString          = "Incoming SW Rad (Wm~S~-2~N~)"
res@tiXAxisString          = "Hours"
res@gsnRightString         = "January"
res@gsnLeftString          = ""+FlxName+""
;res@tiMainString           = "January"
res@pmLegendDisplayMode    = "Always"            ; turn on legend
res@pmLegendSide           = "Top"               ; Change location of
res@pmLegendParallelPosF   = .15                 ; move units right
res@pmLegendOrthogonalPosF = -0.3                ; move units down
res@pmLegendWidthF         = 0.09                ; Change width and
res@pmLegendHeightF        = 0.11                ; height of legend.
res@lgPerimOn              = False               ; turn off box around
res@lgLabelFontHeightF     = .016                ; label font height
if (dims(0).eq.5) then
res@xyExplicitLegendLabels = (/" a10"," a11_r0"," a11_r4"," a11_r5"," fluxnet"/)         ; create explicit labels
else
res@xyExplicitLegendLabels = (/"A1.0","A1.1","Obs"/)         ; create explicit labels
end if
;res@xyExplicitLegendLabels = (/"igbp","A1.0","A1.1","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

prplot                     = new(4,graphic)
prplot(0)                  = gsn_csm_xy(wks3,y,SwJan,res) ; data6
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
;res@gsnLeftString          = ""
;res@tiMainString           = "July"
prplot(1)                  = gsn_csm_xy(wks3,y,SwJul,res) ; data9
;res@tiYAxisString          = ""
res@tiYAxisOn              = True
res@tiYAxisString          = "Incoming LW Rad (Wm~S~-2~N~)"
res@gsnRightString           = "January"
;res@gsnLeftString          = ""
prplot(2)                  = gsn_csm_xy(wks3,y,LwJan,res) ; data6
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
;res@gsnLeftString          = ""
prplot(3)                  = gsn_csm_xy(wks3,y,LwJul,res) ; data9

; -----------------------------------------------------------------------------

pres                       = True
;pres@gsnMaximize           = True
;pres@gsnOrientation        = "landscape"
;pres@PanelRowSpec          = True
pres@gsnPanelYWhiteSpacePercent = 5
gsn_panel(wks3,prplot,(/2,2/),pres)

 if (jpeg) then
 system("ps2epsi "+pltName+".ps "+pltName+".eps")
 system("convert -density 500 -rotate 270 "+pltName+".eps "+pltName+".jpg")
 ;system("convert -density 500 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 end if

end

