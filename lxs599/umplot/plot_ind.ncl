; Procedure Code ==================================================
; created by Lauren Stevens 12 Aug 2011
; =================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; =================================================================

procedure IndPlot(Dir:string,varName:string,Yr:string,pltType:string,pltColor:string,\
                        field:numeric,region:logical,sysFlag:logical,cnLevels:numeric,  \
                        lmask:logical,jpeg:logical)

begin

; Wks ===========================================================

  pltName = Dir + "/" + varName + "_" + Yr + "yrs"

  wks = gsn_open_wks(pltType,pltName)

; Color Map ===========================================================

  gsn_define_colormap(wks,pltColor)

; Resources ========================================================

  res = True
  
  ;res@gsnFrame             = False             ; Don't advance the frame
  ;res@gsnDraw              = False             ; only for Panel plots
  ;                                             ; switch off for individual plots

  res@cnFillOn             = True              ; Turn on contour fill.
  res@cnFillMode           = "RasterFill"      ; Pixels
if (varName.eq."austbias_tscrn".or.varName.eq."austbias_precip".or. \
    varName.eq."aust_tscrn".or.varName.eq."aust_precip" \
.or.varName.eq."apr_djf".or.varName.eq."apr_mam".or.varName.eq."apr_jja" \
.or.varName.eq."apr_son".or.varName.eq."bom_pr_djf".or.varName.eq."bom_pr_mam" \
.or.varName.eq."bom_pr_jja".or.varName.eq."bom_pr_son") then
  res@cnRasterSmoothingOn  = True
end if
  res@cnLineLabelsOn       = False             ; Turn off line labels.
  res@cnInfoLabelOn        = False             ; Turn off info label.
  if(lmask)
  res@cnFillDrawOrder      = "PreDraw"         ; Predraw cn so mask will appear
                                               ; (if land transparent ow land will also be masked)
  end if
  res@cnLinesOn            = False             ; Turn off Lines ?
   
  res@gsnSpreadColors      = True              ;
  ;res@lbLabelBarOn         = False             ;
  res@lbLabelAutoStride    = True              ;
  
  res@gsnLeftString        = ""                ;
  res@gsnRightString       = ""                ;

  if(sysFlag) then
   ;print("label bar scale")
   ;symMinMaxPlt(field,12,False,res)
  end if

;  if(cnLevels.eq."0") then
;    res@cnLevelSelectionMode = "AutomaticLevels"
;  else
    res@cnLevelSelectionMode = "ExplicitLevels"  ; Define your own
    res@cnLevels = cnLevels 
;  end if

  if (varName.eq."bias_precip".or.varName.eq."bias_evap" \
  .or.varName.eq."austbias_precip".or.varName.eq."austbias_precip_bom") then
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

  if (.not.region) then
   Area = "Global"
   latS = -90
   latN =  90 
   lonL =   0
   lonR = 360
  else
   Area = "Australia"
   latS = -45
   latN = -10 
   lonL = 110 
   lonR = 160
  end if 
  ; what about map projections - Arctic and Antarctic

  res@mpMinLatF  = latS
  res@mpMaxLatF  = latN
  res@mpMinLonF  = lonL
  res@mpMaxLonF  = lonR
  
  if(region) then
  res@gsnAddCyclic = False
  end if

  res@mpCenterLonF         = -180.             ; center on dateline

;mask -------------------------------------------------------------------
 
  if(lmask) 
  ;res@mpAreaMaskingOn      = True
  ;res@mpMaskAreaSpecifiers = "Water"
  ;res@mpGridMaskMode       = "MaskOcean"
  res@mpFillOn             = True
  res@mpOceanFillColor     = "white"
  res@mpLandFillColor      = "transparent"
  end if  

;titles -----------------------------------------------------------------  
  
  res@tiMainOn      = True
  res@tiMainString   = (field@long_name + " (" + field@units + ")")
  ;res@tiXAxisString = "Longitude"
  ;res@tiYAxisString = "Latitude"
 
  res@gsnMaximize           = True
  ;res@gsnOrientation        = "landscape"
  ;res@gsnPanelLabelBar      = True
  ;res@gsnPanelRowSpec       = True
  ;res@pmLabelBarWidthF      = 0.8         ; makes label bar bigger
  res@txString              = (field@long_name + " (" + field@units + ")")
 
; Plot ------------------------------------------------------------------
       
  plot = gsn_csm_contour_map(wks,field({latS:latN},{lonL:lonR}),res)

 if(jpeg) then  
 system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 end if

end
