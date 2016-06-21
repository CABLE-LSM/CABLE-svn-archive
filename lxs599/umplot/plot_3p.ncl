; -----------------------------------------------------------------------------
; Lauren Stevens 16 Aug 2011
; 3 panel bias plot
; NCL - contour plots
; -----------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;hdir = systemfunc("echo $hdir")
;usr = systemfunc("echo $USERID")
;loadscript(hdir+"/"+usr+"/umplot/bias_3panel.ncl")

procedure Panel3Plot(pltName:string,top:numeric,middle:numeric,bottom:numeric,jpeg:logical)

begin

a14 = False ;True
hdir = systemfunc("echo $hdir")
usr = systemfunc("echo $USERID")
numYr = systemfunc("echo $YR")
;jpg = systemfunc("echo $CHFMT")
;if(jpg.eq."y")then
;jpeg = True
;else
;jpeg = False
;end if

PltNam  = pltName+"_"+numYr+"yrs"
wks     = gsn_open_wks("ps",pltName+"_"+numYr+"yrs")
if (pltName.eq."bias_ann_prec_3p") then
 gsn_define_colormap(wks,"BlueDarkRed18")
 gsn_reverse_colormap(wks)
else
 cmap = RGBtoCmap(hdir+"/"+usr+"/umplot/colormap/BlueOrangeRed.rgb")
 gsn_define_colormap(wks,cmap)
 ;gsn_define_colormap(wks,"BlueDarkRed18");pltColor)
end if

res               = True
res@gsnDraw       = False
res@gsnFrame      = False

res@cnFillOn             = True              ; Turn on contour fill.
res@cnFillMode           = "RasterFill"      ; Pixels
res@cnLineLabelsOn       = False             ; Turn off line labels.
res@cnInfoLabelOn        = False             ; Turn off info label.
 ; if(lmask)
 ; res@cnFillDrawOrder      = "PreDraw"         ; Predraw cn so mask will appear
 ;                                              ; (if land transparent ow land will also be masked)
 ; end if
res@cnLinesOn            = False             ; Turn off Lines ?

res@gsnSpreadColors      = True              ;
res@lbLabelBarOn         = False             ;
res@lbLabelAutoStride    = True              ;

res@gsnLeftString        = ""                ;
res@gsnRightString       = ""                ;
res@mpCenterLonF         = -180. 
res@cnLevelSelectionMode = "ExplicitLevels"  ; Define your own
if (pltName.eq."bias_ann_prec_3p") then
 res@cnLevels = (/-8,-4,-2,-1,-0.5,0,0.5,1,2,4,8/)
else
 res@cnLevels = ispan(-9,9,1) ;cnLevels
 ;res@cnLevels = (/-8,-4,-2,-1,0,1,2,4,8/) ;cnLevels
end if

latS = -90
latN =  90
lonL =   0
lonR = 360

res@mpMinLatF  = latS
res@mpMaxLatF  = latN
res@mpMinLonF  = lonL
res@mpMaxLonF  = lonR

plot              = new(3,graphic)

;res@tiYAxisString = "Temperature (deg K)"
;res@tiXAxisString = ""
plot(0)           = gsn_csm_contour_map(wks,top,res)
plot(1)           = gsn_csm_contour_map(wks,middle,res)
plot(2)           = gsn_csm_contour_map(wks,bottom,res)

; Panel Plot -----------------------------------------------------------------

pres                      = True
pres@gsnMaximize          = True
pres@txString             = ""
if (a14) then
pres@gsnPanelFigureStrings = (/"A1.4","A1.3","A1.4-A1.3"/)
else
;pres@gsnPanelFigureStrings = (/"cable","moses","cable-moses"/)
pres@gsnPanelFigureStrings = (/"new (c)","old (m)","new-old"/)
end if
pres@gsnPanelLabelBar      = True
pres@gsnPanelBottom       = 0.18
;pres@txFontHeightF       = 0.03
;pres@pmLegendDisplayMode = "Always"
gsn_panel(wks,plot,(/3,1/),pres)

if(jpeg)then
 system("ps2epsi "+PltNam+".ps "+PltNam+".eps")
 system("convert -density 500 -rotate 270 "+PltNam+".eps "+PltNam+".jpg")
 ;system("convert -density 1000 -rotate 270 "+PltNam+".ps "+PltNam+".jpg")
 ;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+PltNam+".ps "+PltNam+".jpg")
end if

end

