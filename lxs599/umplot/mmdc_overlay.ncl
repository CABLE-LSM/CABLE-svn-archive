;==============================================================================
; Lauren Stevens 8 July 2011
;==============================================================================
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;==============================================================================

procedure tso_flux(FlxName:string,tstep,ShJan,LhJan,ShJul,LhJul) ; :numeric

begin

;  pltName = Dir + "/" + varName + "_" + Yr + "yrs"
;  wks = gsn_open_wks(pltType,pltName)

numYr  = systemfunc("echo $YR")
hmax   = 24 - (24/tstep)
y      = fspan(0,hmax,tstep)
;y2     = fspan(0,23+2/3,72)

;-----------------------------------------------------------------------------

wks                        = gsn_open_wks("ps","mmdc_"+FlxName+"_flux_"+numYr+"yrs")
colors                     = (/"green","blue","red","black"/)
labels                     = (/"igbp","moses","cable","fluxnet"/)
res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2,2/)
res@xyDashPattern          = 0.
res@tiYAxisString          = "Sensible Heat Flux (Wm-2)"
res@tiXAxisString          = "Hours"
res@tiMainString           = "January"
;res@pmLegendDisplayMode    = "Always"            ; turn on legend
;res@pmLegendSide           = "Top"               ; Change location of 
;res@pmLegendParallelPosF   = .15                 ; move units right
;res@pmLegendOrthogonalPosF = -0.3                ; move units down
;res@pmLegendWidthF         = 0.09                ; Change width and
;res@pmLegendHeightF        = 0.11                ; height of legend.
;res@lgPerimOn              = False               ; turn off box around
;res@lgLabelFontHeightF     = .018                ; label font height
;res@xyExplicitLegendLabels = (/"igbp","moses","cable","fluxnet"/)         ; create explicit labels
;res@trYMinF                = min((/sens,msh,csh,qh_ta/)) or data(0),data(1)etc
;res@trYMaxF                = max((/sens,msh,csh,qh_ta/))

;------------------------------------------------------------------------------

npanels          = 4
nplots_per_panel = 4
plot                       = new((/npanels,nplots_per_panel/),graphic)

;do i = 0,npanels-1,1
res@trYMinF                = min(ShJan)-10
res@trYMaxF                = max(ShJan)+10
res@xyLineColors           = colors(0)
plot(0,0)                  = gsn_csm_xy(wks,y,ShJan(0,:),res) ; data2
res@xyLineColors           = colors(1)
plot(0,1)                  = gsn_csm_xy(wks,y,ShJan(1,:),res) ; data2
res@xyLineColors           = colors(2)
plot(0,2)                  = gsn_csm_xy(wks,y,ShJan(2,:),res) ; data2
res@xyLineColors           = colors(3)
plot(0,3)                  = gsn_csm_xy(wks,y,ShJan(3,:),res) ; data2

res@tiYAxisString          = ""
res@tiMainString           = "July"
res@trYMinF                = min(ShJul)-10
res@trYMaxF                = max(ShJul)+10
res@xyLineColors           = colors(0)
plot(1,0)                  = gsn_csm_xy(wks,y,ShJul(0,:),res) ; data
res@xyLineColors           = colors(1)
plot(1,1)                  = gsn_csm_xy(wks,y,ShJul(1,:),res) ; data
res@xyLineColors           = colors(2)
plot(1,2)                  = gsn_csm_xy(wks,y,ShJul(2,:),res) ; data
res@xyLineColors           = colors(3)
plot(1,3)                  = gsn_csm_xy(wks,y,ShJul(3,:),res) ; data

res@tiYAxisString          = "Latent Heat Flux (Wm-2)"
res@tiMainString           = ""
res@trYMinF                = min(LhJan)-10
res@trYMaxF                = max(LhJan)+10
res@xyLineColors           = colors(0)
plot(2,0)                  = gsn_csm_xy(wks,y,LhJan(0,:),res) ; data3
res@xyLineColors           = colors(1)
plot(2,1)                  = gsn_csm_xy(wks,y,LhJan(1,:),res) ; data3
res@xyLineColors           = colors(2)
plot(2,2)                  = gsn_csm_xy(wks,y,LhJan(2,:),res) ; data3
res@xyLineColors           = colors(3)
plot(2,3)                  = gsn_csm_xy(wks,y,LhJan(3,:),res) ; data3

res@tiMainString           = ""
res@tiYAxisString          = "" 
res@trYMinF                = min(LhJul)-10
res@trYMaxF                = max(LhJul)+10
res@xyLineColors           = colors(0)
plot(3,0)                  = gsn_csm_xy(wks,y,LhJul(0,:),res) ; data1
res@xyLineColors           = colors(1)
plot(3,1)                  = gsn_csm_xy(wks,y,LhJul(1,:),res) ; data1
res@xyLineColors           = colors(2)
plot(3,2)                  = gsn_csm_xy(wks,y,LhJul(2,:),res) ; data1
res@xyLineColors           = colors(3)
plot(3,3)                  = gsn_csm_xy(wks,y,LhJul(3,:),res) ; data1

do i = 0,npanels-1,1
 do j=1,3
  overlay(plot(i,0),plot(i,j))

  lgres                    = True
  lgres@lgLineColors       = colors
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.13           ; width of legend (NDC)
  lgres@vpHeightF          = 0.10           ; height of legend (NDC)
  lgres@lgPerimOn          = False
  ;lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0

  legend = gsn_create_legend (wks, 4, labels,lgres)

  amres                    = True
  amres@amJust             = "TopLeft"    ; Use bottom right corner of box
                                            ; for determining its location.
  amres@amParallelPosF     = -0.46              ; Move legend to right
  amres@amOrthogonalPosF   = -0.47              ; Move legend down.

  annoid = gsn_add_annotation(plot(i,0),legend,amres)

 end do
end do

; -----------------------------------------------------------------------------

pres                       = True
pres@gsnMaximize           = True
;pres@gsnOrientation        = "landscape"
;pres@PanelRowSpec          = True
gsn_panel(wks,plot(:,0),(/2,2/),pres)

end

;===============================================================================

procedure tso_rad(FlxName:string,RnJan,TaJan,RnJul,TaJul) ; :numeric

begin

numYr  = systemfunc("echo $YR")
y      = fspan(0,23.5,48)
;y2     = fspan(0,23.43,72)


;-----------------------------------------------------------------------------

wks2                       = gsn_open_wks("ps","mmdc_"+FlxName+"_rad_"+numYr+"yrs")

colors                     = (/"green","blue","red","black"/)
labels                     = (/"igbp","moses","cable","fluxnet"/)
res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2,2/)
res@xyDashPattern          = 0.
res@tiYAxisString          = "Net Radiation (Wm-2)"
res@tiXAxisString          = "Hours"
res@tiMainString           = "January"
;res@pmLegendDisplayMode    = "Always"            ; turn on legend
;res@pmLegendSide           = "Top"               ; Change location of
;res@pmLegendParallelPosF   = .15                 ; move units right
;res@pmLegendOrthogonalPosF = -0.3                ; move units down
;res@pmLegendWidthF         = 0.09                ; Change width and
;res@pmLegendHeightF        = 0.11                ; height of legend.
;res@lgPerimOn              = False               ; turn off box around
;res@lgLabelFontHeightF     = .018                ; label font height
;res@xyExplicitLegendLabels = (/"igbp","moses","cable","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

npanels          = 4
nplots_per_panel = 4

radplt                     = new((/npanels,nplots_per_panel/),graphic)

;do i =0,npanels-1,1
res@trYMinF                = min(RnJan)-10
res@trYMaxF                = max(RnJan)+10
res@xyLineColors           = colors(0)
radplt(0,0)                = gsn_csm_xy(wks2,y,RnJan(0,:),res) ; data5
res@xyLineColors           = colors(1)
radplt(0,1)                = gsn_csm_xy(wks2,y,RnJan(1,:),res) ; data5
res@xyLineColors           = colors(2)
radplt(0,2)                = gsn_csm_xy(wks2,y,RnJan(2,:),res) ; data5
res@xyLineColors           = colors(3)
radplt(0,3)                = gsn_csm_xy(wks2,y,RnJan(3,:),res) ; data5

res@tiYAxisString          = ""
res@tiMainString           = "July"
res@trYMinF                = min(RnJul)-10
res@trYMaxF                = max(RnJul)+10
res@xyLineColors           = colors(0)
radplt(1,0)                = gsn_csm_xy(wks2,y,RnJul(0,:),res) ; data8
res@xyLineColors           = colors(1)
radplt(1,1)                = gsn_csm_xy(wks2,y,RnJul(1,:),res) ; data8
res@xyLineColors           = colors(2)
radplt(1,2)                = gsn_csm_xy(wks2,y,RnJul(2,:),res) ; data8
res@xyLineColors           = colors(3)
radplt(1,3)                = gsn_csm_xy(wks2,y,RnJul(3,:),res) ; data8

res@tiYAxisString          = "Temperature (deg K)"
res@tiMainString           = ""
res@trYMinF                = min(TaJan)-2
res@trYMaxF                = max(TaJan)+2
res@xyLineColors           = colors(0)
radplt(2,0)                = gsn_csm_xy(wks2,y,TaJan(0,:),res) ; data4
res@xyLineColors           = colors(1)
radplt(2,1)                = gsn_csm_xy(wks2,y,TaJan(1,:),res) ; data4
res@xyLineColors           = colors(2)
radplt(2,2)                = gsn_csm_xy(wks2,y,TaJan(2,:),res) ; data4
res@xyLineColors           = colors(3)
radplt(2,3)                = gsn_csm_xy(wks2,y,TaJan(3,:),res) ; data4

res@tiYAxisString          = ""
res@tiMainString           = ""
res@trYMinF                = min(TaJul)-2
res@trYMaxF                = max(TaJul)+2
res@xyLineColors           = colors(0)
radplt(3,0)                = gsn_csm_xy(wks2,y,TaJul(0,:),res) ; data7
res@xyLineColors           = colors(1)
radplt(3,1)                = gsn_csm_xy(wks2,y,TaJul(1,:),res) ; data7
res@xyLineColors           = colors(2)
radplt(3,2)                = gsn_csm_xy(wks2,y,TaJul(2,:),res) ; data7
res@xyLineColors           = colors(3)
radplt(3,3)                = gsn_csm_xy(wks2,y,TaJul(3,:),res) ; data7

do i =0,npanels-1,1
  do j=1,3
   overlay(radplt(i,0),radplt(i,j))

  lgres                    = True
  lgres@lgLineColors       = colors
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.13           ; width of legend (NDC)
  lgres@vpHeightF          = 0.10           ; height of legend (NDC)
  lgres@lgPerimOn          = False
  ;lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0

  legend = gsn_create_legend (wks2, 4, labels,lgres)

  amres                    = True
  amres@amJust             = "TopLeft"    ; Use bottom right corner of box
                                            ; for determining its location.
  amres@amParallelPosF     = -0.46              ; Move legend to right
  amres@amOrthogonalPosF   = -0.47              ; Move legend down.

  annoid = gsn_add_annotation(radplt(i,0),legend,amres)

  end do
end do

; -----------------------------------------------------------------------------

pres                       = True
pres@gsnMaximize           = True
;pres@gsnOrientation        = "landscape"
;pres@PanelRowSpec          = True
gsn_panel(wks2,radplt(:,0),(/2,2/),pres)

end

;===============================================================================

procedure tso_prec(FlxName:string,PrJan,PrJul) ; :numeric

begin

numYr  = systemfunc("echo $YR")
y      = fspan(0,23.5,48)
;y2     = fspan(0,23.43,72)

;-----------------------------------------------------------------------------

wks3                       = gsn_open_wks("ps","mmdc_"+FlxName+"_prec_"+numYr+"yrs")

colors                     = (/"green","blue","red","black"/)
labels                     = (/"igbp","moses","cable","fluxnet"/)
res                        = True
res@gsnDraw                = False
res@gsnFrame               = False
;res@xyLineColors           = (/"green","blue","red","black"/)
res@xyLineThicknesses      = (/2,2,2,2/)
res@xyDashPattern          = 0.
res@tiYAxisString          = "Precipitation (mm/day)"
res@tiXAxisString          = "Hours"
res@tiMainString           = "January"
;res@pmLegendDisplayMode    = "Always"            ; turn on legend
;res@pmLegendSide           = "Top"               ; Change location of
;res@pmLegendParallelPosF   = .15                 ; move units right
;res@pmLegendOrthogonalPosF = -0.3                ; move units down
;res@pmLegendWidthF         = 0.09                ; Change width and
;res@pmLegendHeightF        = 0.11                ; height of legend.
;res@lgPerimOn              = False               ; turn off box around
;res@lgLabelFontHeightF     = .018                ; label font height
;res@xyExplicitLegendLabels = (/"igbp","moses","cable","fluxnet"/)         ; create explicit labels

;------------------------------------------------------------------------------

npanels          = 2
nplots_per_panel = 4

prplot                     = new((/npanels,nplots_per_panel/),graphic)

res@trYMinF                = min(PrJan)-2
res@trYMaxF                = max(PrJan)+2
res@xyLineColors           = colors(0)
prplot(0,0)                = gsn_csm_xy(wks3,y,PrJan(0,:),res) ; data6
res@xyLineColors           = colors(1)
prplot(0,1)                = gsn_csm_xy(wks3,y,PrJan(1,:),res) ; data6
res@xyLineColors           = colors(2)
prplot(0,2)                = gsn_csm_xy(wks3,y,PrJan(2,:),res) ; data6
res@xyLineColors           = colors(3)
prplot(0,3)                = gsn_csm_xy(wks3,y,PrJan(3,:),res) ; data6

res@tiYAxisString          = ""
res@tiMainString           = "July"
res@trYMinF                = min(PrJul)-2
res@trYMaxF                = max(PrJul)+2
res@xyLineColors           = colors(0)
prplot(1,0)                = gsn_csm_xy(wks3,y,PrJul(0,:),res) ; data9
res@xyLineColors           = colors(1)
prplot(1,1)                = gsn_csm_xy(wks3,y,PrJul(1,:),res) ; data9
res@xyLineColors           = colors(2)
prplot(1,2)                = gsn_csm_xy(wks3,y,PrJul(2,:),res) ; data9
res@xyLineColors           = colors(3)
prplot(1,3)                = gsn_csm_xy(wks3,y,PrJul(3,:),res) ; data9

do i = 0,npanels-1,1
  do j=1,3
   overlay(prplot(i,0),prplot(i,j))

  lgres                    = True
  lgres@lgLineColors       = colors
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.13           ; width of legend (NDC)
  lgres@vpHeightF          = 0.10           ; height of legend (NDC)
  lgres@lgPerimOn          = False
  ;lgres@lgPerimThicknessF  = 2.0            ; thicken the box perimeter
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0

  legend = gsn_create_legend (wks3, 4, labels,lgres)

  amres                    = True
  amres@amJust             = "TopLeft"    ; Use bottom right corner of box
                                            ; for determining its location.
  amres@amParallelPosF     = -0.46              ; Move legend to right
  amres@amOrthogonalPosF   = -0.47              ; Move legend down.

  annoid = gsn_add_annotation(prplot(i,0),legend,amres)

  end do
end do

; -----------------------------------------------------------------------------

pres                       = True
pres@gsnMaximize           = True
;pres@gsnOrientation        = "landscape"
;pres@PanelRowSpec          = True
gsn_panel(wks3,prplot(:,0),(/1,2/),pres)

end

