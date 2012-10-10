; -----------------------------------------------------------------------------
; Lauren Stevens 31 Aug 2011
; Bar Chart plot
; NCL - XY Plots
; -----------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Set-Up ----------------------------------------------------------------------

procedure BarPlot(Title:string,data:numeric,pltName:string,pltType:string,Panel:logical)

begin

;print(data)

if(Panel)then
nDim   = dimsizes(data)
nBars  = nDim(1)/7
else
nBars  = dimsizes(data)/7
end if

if(nBars.eq.2) then
x      = (/1,2,3.5,4.5,6,7,8.5,9.5,11,12,13.5,14.5,16,17/)
;pos    = (/.75,3.25,5.75,8.25,10.75,13.25,15.75/)
models = (/"MOSES-ERAi","CABLE-ERAi"/)
labels = (/"Global","hNLat","mNlat","Trop","mSLat","hSLat","Land"/)
colors = (/"cornflowerblue","darkorchid4"/)
else
x      = (/1,2,3,4.5,5.5,6.5,8,9,10,11.5,12.5,13.5,15,16,17,18.5,19.5,20.5,22,23,24/)
models = (/"MOSES","CABLE","ERAi"/)
labels = (/"Global","hNLat","mNlat","Trop","mSLat","hSLat","b1500m"/)
colors = (/"cornflowerblue","darkorchid4","seagreen"/)
end if
nlab   = dimsizes(labels)

wks    = gsn_open_wks(pltType,pltName)

res                          = True
if(Panel)then
res@gsnDraw                  = False
end if
res@gsnFrame                 = False
res@gsnXYBarChart            = True
;res@gsnXYBarChartWidth       = 0.75
res@gsnXYBarChartColors2     = colors
;res@tmXBOn                   = False
res@trYMinF                  = floor(min(data))
res@trYMaxF                  = ceil(max(data))
res@trXMinF                  = 0
if(nBars.eq.2)then
res@trXMaxF                  = 18
else
res@trXMaxF                  = 25
end if
;res@tiMainString             = "Screen Temperature deg C"
;res@gsnRightString           = "Temperature at 1.5m (deg C)"
res@gsnYRefLine              = 0
;res@gsnAboveYRefLineBarColors = colors 
;res@gsnMaximize              = True
res@tmXBMode                 = "Explicit"
if(nBars.eq.2)then
res@tmXBValues               = (/1.5,4,6.5,9,11.5,14,16.5/)
else
res@tmXBValues               = (/2,5.5,9,12.5,16,19.5,23/)
end if
res@tmXBLabels               = labels

res@tmXBLabelFontHeightF        = 0.015
;res@gsnPanelBottom              = 0.13
res@tmXMajorGrid                = True
res@tmXMajorGridLineColor       = 1
res@tmXMajorGridLineDashPattern = 1

;res@tmYLMajorLengthF            = 0
;res@tmYLMajorOutwardLengthF     = 0.  
;res@tmXBMajorLengthF            = 0.01
;res@tmXBMajorOutwardLengthF     = 0.01
;res@tmXTOn                      = False 
;res@tmYMajorGrid                = True
;res@tmYMajorGridLineColor       = 1
;res@tmYMajorGridLineDashPattern = 1
;res@xyCurveDrawOrder            = "PostDraw"

if(Panel)then
plot              = new(4,graphic)
;res@tiMainString  = ""
res@tmXTOn        = False
res@gsnLeftString = "DJF"
plot(0) = gsn_csm_xy(wks,x,data(0,:),res)
res@gsnLeftString = "MAM"
plot(1) = gsn_csm_xy(wks,x,data(1,:),res)
res@gsnLeftString = "JJA"
plot(2) = gsn_csm_xy(wks,x,data(2,:),res)
res@gsnLeftString = "SON"
plot(3) = gsn_csm_xy(wks,x,data(3,:),res)

pres                       = True
pres@gsnMaximize           = True
pres@txString              = Title
;pres@gsnPanelBottom        = 0.18   ; PUT IN LATER WHEN HAVE LEGEND
;pres@txFontHeightF        = 0.03
;pres@pmLegendDisplayMode  = "Always"
;pres@gsnPanelFigureStrings = (/"DJF","MAM","JJA","SON"/)
;pres@gsnPanelFigureStringsBackgroundFillColor = "Transparent"
;pres@ = 0.15
;pres@gsnOrientation       = "Landscape"
;pres@gsnPanelRowSpec      = True
gsn_panel(wks,plot,(/2,2/),pres)

else
res@tiMainString                = Title
res@gsnMaximize                 = True
res@gsnPanelBottom              = 0.13
plot = gsn_csm_xy(wks,x,data,res)

lbres                    = True          ; labelbar only resources
lbres@vpWidthF           = 0.2 ;* vpw     ; labelbar width
lbres@vpHeightF          = 0.1 ;* vph     ; labelbar height
lbres@lbBoxMajorExtentF  = 0.15          ; puts space between color boxes
;lbres@lbFillColors       = res@gsnXYBarChartColors ; labelbar colors
lbres@lbMonoFillPattern  = True          ; Solid fill pattern
lbres@lbLabelFontHeightF = 0.02          ; font height. default is small
lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
lbres@lbPerimOn          = False
 
 do i=0,nBars-1,1
  lbres@lbFillColors = colors(i)  
  if (nBars.eq.2) then
  xpos = (/0.25,0.55/)
  else
  xpos = (/0.18,0.45,0.70/)
  end if
  ;ypos = (/0.22,0.26/)
  gsn_labelbar_ndc(wks,1,models(i),xpos(i),0.08,lbres)
  ;gsn_labelbar_ndc(wks,1,models(i),0.20,ypos(i),lbres)
end do
;draw(plot)
frame(wks)

end if

end

