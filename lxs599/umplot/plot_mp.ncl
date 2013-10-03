load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; Procedure Code ==================================================
; created by Lauren Stevens 11-21 April 2011
; during the NCL cource at BMTC.
; =================================================================

procedure Polar_Plot(Dir:string,varName:string,Yr:string,pltType:string,pltColor:string,\
                     field:numeric,region:logical,smthFlag:logical,cnLevels:numeric, \
                     Hem:string,jpeg:logical)

begin

  pltName = Dir + "/" + varName + "_" + Yr + "yrs"

  wks = gsn_open_wks(pltType,pltName)

; Color Map ===========================================================

  gsn_define_colormap(wks,pltColor)

; Resources ========================================================

  res = True
  
  res@gsnFrame             = False             ; Don't advance the frame
  res@gsnDraw              = False             ; only for Panel plots
  ;                                            ; switch off for individual plots

  res@gsnPolar = "NH"                          ; Default
  ;res@gsnPolar = "SH"

  res@cnFillOn             = True              ; Turn on contour fill.
  res@cnFillMode           = "RasterFill"      ; Pixels
  res@cnLineLabelsOn       = False             ; Turn off line labels.
  res@cnInfoLabelOn        = False             ; Turn off info label.
  ;res@cnFillDrawOrder      = "PreDraw"         ; Predraw cn so mask will appear
                                               ; (if land transparent ow land will also be masked)
  res@cnLinesOn            = False             ; Turn off Lines ?
   
  res@gsnSpreadColors      = True              ;
  res@lbLabelBarOn         = False             ;
  res@lbLabelAutoStride    = True              ;
  
  res@gsnLeftString        = ""                ;
  res@gsnRightString       = ""                ;

  ;if(sysFlag) then
  ; ;print("label bar scale")
  ; ;symMinMaxPlt(field,12,False,res)
  ;end if

;  if(cnLevels.eq."0") then
;    res@cnLevelSelectionMode = "AutomaticLevels"
;  else
    res@cnLevelSelectionMode = "ExplicitLevels"  ; Define your own
    res@cnLevels = cnLevels 
;  end if

  if (varName.eq."bias_precip".or.varName.eq."bias_evap") then
   gsn_reverse_colormap(wks)
  end if

  ;if (varName.eq."canht")
  ; res@cnLevels = (/0:28:2/)        ; ERROR
   ;mnlvl = 0                        ; Minimum contour level.
   ;mxlvl = 28                       ; Maximum contour level.
   ;spcng = 2                        ; Contour level spacing.
   ;ncn   = (mxlvl-mnlvl)/spcng + 1  ; Number of contour levels.
  ;end if
  
  ;res@cnLevelSelectionMode = "ManualLevels" ; Define your own
  ;res@cnMinLevelValF       = mnlvl          ; contour levels.
  ;res@cnMaxLevelValF       = mxlvl
  ;res@cnLevelSpacingF      = spcng
  
;region -----------------------------------------------------------------

  ;if (.not.region) then
  ; Area = "Global"
  ; latS = -90
  ; latN =  90 
  ; lonL =   0
  ; lonR = 360
  ;else
  ; Area = "Australia"
  ; latS = -45
  ; latN = -10 
  ; lonL = 110 
  ; lonR = 160
  ;end if 
  ; what about map projections - Arctic and Antarctic

  ;res@mpMinLatF  = latS
  ;res@mpMaxLatF  = latN
  ;res@mpMinLonF  = lonL
  ;res@mpMaxLonF  = lonR
  
  ;if(region) then
  ;res@gsnAddCyclic = False
  ;end if

  res@gsnPolar = Hem

  res@mpCenterLonF         = -180.             ; center on dateline

;mask -------------------------------------------------------------------
  
  ;;res@mpAreaMaskingOn      = True
  ;;res@mpMaskAreaSpecifiers = "Water"
  ;;res@mpGridMaskMode       = "MaskOcean"
  ;res@mpFillOn             = True
  ;res@mpOceanFillColor     = "white"
  ;res@mpLandFillColor      = "transparent"
  
;titles -----------------------------------------------------------------  
  
  res@tiMainOn      = False
  ;res@tiXAxisString = "Longitude"
  ;res@tiYAxisString = "Latitude"
  
; Plot ------------------------------------------------------------------
       
  plot = new(4,graphic)
    
  res@tmXBLabelsOn = False
  ;plot(0) = gsn_csm_contour_map_polar(wks,field(0,{latS:latN},{lonL:lonR}),res)
  plot(0) = gsn_csm_contour_map_polar(wks,field(0,:,:),res)
  res@tmYLLabelsOn = False
  ;plot(1) = gsn_csm_contour_map_polar(wks,field(1,{latS:latN},{lonL:lonR}),res)
  plot(1) = gsn_csm_contour_map_polar(wks,field(1,:,:),res)
  res@tmXBLabelsOn = True
  res@tmYLLabelsOn = True
  ;plot(2) = gsn_csm_contour_map_polar(wks,field(2,{latS:latN},{lonL:lonR}),res)
  plot(2) = gsn_csm_contour_map_polar(wks,field(2,:,:),res)
  res@tmYLLabelsOn = False
  ;plot(3) = gsn_csm_contour_map_polar(wks,field(3,{latS:latN},{lonL:lonR}),res)
  plot(3) = gsn_csm_contour_map_polar(wks,field(3,:,:),res)  

; Paneled Plot ====================================================

  pres                       = True
  pres@gsnMaximize           = True
  pres@gsnOrientation        = "landscape"
  pres@gsnPanelLabelBar      = True
  pres@gsnPanelFigureStrings = (/"DJF","MAM","JJA","SON"/)
  pres@gsnPanelRowSpec       = True
  pres@pmLabelBarWidthF      = 0.8         ; makes label bar bigger
  pres@txString              = (field@long_name + " (" + field@units + ")")
  
  gsn_panel(wks,plot,(/2,2/),pres)

 if(jpeg) then
 system("convert -density 1000 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 ;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 end if

end

