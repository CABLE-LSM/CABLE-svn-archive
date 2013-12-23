;====================================================================
; Lauren Stevens 8 July 2011
;====================================================================
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;====================================================================

procedure ts_flux(FlxName:string,ShJan,LhJan,ShJul,LhJul,jpeg) ; :numeric

begin

;  pltName = Dir + "/" + varName + "_" + Yr + "yrs"
;  wks = gsn_open_wks(pltType,pltName)

numYr  = systemfunc("echo $YR")
y = fspan(0,23.5,48)

;-----------------------------------------------------------------------------

pltName                    = "mmdc_"+FlxName+"_flux_"+numYr+"yrs"
wks                        = gsn_open_wks("ps","mmdc_"+FlxName+"_flux_"+numYr+"yrs")

res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
res@xyLineColors           = (/"blue","red","black"/)
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2/)
;res@xyLineThicknesses      = (/2,2,2,2/)
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
res@lgLabelFontHeightF     = .018                ; label font height
res@xyExplicitLegendLabels = (/"A1.0","A1.1","fluxnet"/)         ; create explicit labels
;res@xyExplicitLegendLabels = (/"igbp","A1.0","A1.1","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

plot                       = new(4,graphic)
plot(0)                    = gsn_csm_xy(wks,y,ShJan,res) ; data2
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString        = "July"
res@gsnLeftString          = ""
;res@tiMainString           = "July"
plot(1)                    = gsn_csm_xy(wks,y,ShJul,res) ; data
res@tiYAxisOn              = True
res@tiYAxisString          = "Latent Heat Flux (Wm~S~-2~N~)"
res@gsnRightString        = "January"
;res@tiMainString           = ""
plot(2)                    = gsn_csm_xy(wks,y,LhJan,res) ; data3
res@gsnRightString        = "July"
;res@tiMainString           = ""
;res@tiYAxisString          = "" 
res@tiYAxisOn              = False
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
y = fspan(0,23.5,48)

;-----------------------------------------------------------------------------

pltName                    = "mmdc_"+FlxName+"_rad_"+numYr+"yrs"
wks2                       = gsn_open_wks("ps","mmdc_"+FlxName+"_rad_"+numYr+"yrs")

res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
res@xyLineColors           = (/"blue","red","black"/)
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2/)
;res@xyLineThicknesses      = (/2,2,2,2/)
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
res@lgLabelFontHeightF     = .018                ; label font height
res@xyExplicitLegendLabels = (/"A1.0","A1.1","fluxnet"/)         ; create explicit labels
;res@xyExplicitLegendLabels = (/"igbp","A1.0","A1.1","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

radplt                     = new(4,graphic)
radplt(0)                  = gsn_csm_xy(wks2,y,RnJan,res) ; data5
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
res@gsnLeftString          = ""
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
y = fspan(0,23.5,48)

;-----------------------------------------------------------------------------

pltName                    = "mmdc_"+FlxName+"_prec_"+numYr+"yrs"
wks3                       = gsn_open_wks("ps","mmdc_"+FlxName+"_prec_"+numYr+"yrs")

res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
res@xyLineColors           = (/"blue","red","black"/)
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2/)
;res@xyLineThicknesses      = (/2,2,2,2/)
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
res@lgLabelFontHeightF     = .018                ; label font height
res@xyExplicitLegendLabels = (/"A1.0","A1.1","fluxnet"/)         ; create explicit labels
;res@xyExplicitLegendLabels = (/"igbp","A1.0","A1.1","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

prplot                     = new(4,graphic)
prplot(0)                  = gsn_csm_xy(wks3,y,PrJan,res) ; data6
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
res@gsnLeftString          = ""
;res@tiMainString           = "July"
prplot(1)                  = gsn_csm_xy(wks3,y,PrJul,res) ; data9
;res@tiYAxisString          = ""
res@tiYAxisOn              = True
res@tiYAxisString          = "Total Cloud Fraction"
res@gsnRightString           = "January"
res@gsnLeftString            = ""
prplot(2)                  = gsn_csm_xy(wks3,y,ClJan,res) ; data6
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
res@gsnLeftString            = ""
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
y = fspan(0,23.5,48)

;-----------------------------------------------------------------------------

pltName                    = "mmdc_"+FlxName+"_prec_"+numYr+"yrs"
wks3                       = gsn_open_wks("ps","mmdc_"+FlxName+"_prec_"+numYr+"yrs")

res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
res@xyLineColors           = (/"blue","red","black"/)
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2/)
;res@xyLineThicknesses      = (/2,2,2,2/)
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
res@lgLabelFontHeightF     = .018                ; label font height
res@xyExplicitLegendLabels = (/"A1.0","A1.1","fluxnet"/)         ; create explicit labels
;res@xyExplicitLegendLabels = (/"igbp","A1.0","A1.1","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

prplot                     = new(2,graphic)
prplot(0)                  = gsn_csm_xy(wks3,y,PrJan,res) ; data6
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
res@gsnLeftString          = ""
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
y = fspan(0,23.5,48)

;-----------------------------------------------------------------------------

pltName                    = "mmdc_"+FlxName+"_swlw_"+numYr+"yrs"
wks3                       = gsn_open_wks("ps","mmdc_"+FlxName+"_swlw_"+numYr+"yrs")

res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
res@xyLineColors           = (/"blue","red","black"/)
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2/)
;res@xyLineThicknesses      = (/2,2,2,2/)
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
res@lgLabelFontHeightF     = .018                ; label font height
res@xyExplicitLegendLabels = (/"A1.0","A1.1","fluxnet"/)         ; create explicit labels
;res@xyExplicitLegendLabels = (/"igbp","A1.0","A1.1","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

prplot                     = new(4,graphic)
prplot(0)                  = gsn_csm_xy(wks3,y,SwJan,res) ; data6
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
res@gsnLeftString          = ""
;res@tiMainString           = "July"
prplot(1)                  = gsn_csm_xy(wks3,y,SwJul,res) ; data9
;res@tiYAxisString          = ""
res@tiYAxisOn              = True
res@tiYAxisString          = "Incoming LW Rad (Wm~S~-2~N~)"
res@gsnRightString           = "January"
res@gsnLeftString          = ""
prplot(2)                  = gsn_csm_xy(wks3,y,LwJan,res) ; data6
;res@tiYAxisString          = ""
res@tiYAxisOn              = False
res@gsnRightString           = "July"
res@gsnLeftString          = ""
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

