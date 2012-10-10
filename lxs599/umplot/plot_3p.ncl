; -----------------------------------------------------------------------------
; Lauren Stevens 16 Aug 2011
; 3 panel bias plot
; NCL - contour plots
; -----------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "/home/cmar/ste69f/umplot/bias_3panel.ncl"

procedure Panel3Plot(pltName:string,top:numeric,middle:numeric,bottom:numeric,jpeg:logical)

begin

numYr = systemfunc("echo $YR")
;jpg = systemfunc("echo $JPEG")
;if(jpg.eq."y")then
;jpeg = True
;else
;jpeg = False
;end if

PltNam            = pltName+"_"+numYr+"yrs"
wks               = gsn_open_wks("ps",pltName+"_"+numYr+"yrs")
cmap = RGBtoCmap("/home/cmar/ste69f/umplot/colormap/BlueOrangeRed.rgb")
gsn_define_colormap(wks,cmap)
;gsn_define_colormap(wks,"BlueDarkRed18");pltColor)

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
    res@cnLevels = ispan(-9,9,1) ;cnLevels   
    ;res@cnLevels = (/-8,-4,-2,-1,0,1,2,4,8/) ;cnLevels   

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

panel_res                      = True
panel_res@gsnMaximize          = True
panel_res@txString             = ""
panel_res@gsnPanelFigureStrings = (/"cable","moses","cable-moses"/)
panel_res@gsnPanelLabelBar      = True
panel_res@gsnPanelBottom       = 0.18
;panel_res@txFontHeightF       = 0.03
;panel_res@pmLegendDisplayMode = "Always"
gsn_panel(wks,plot,(/3,1/),panel_res)

if(jpeg)then
system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+PltNam+".ps "+PltNam+".jpg")
end if

end

