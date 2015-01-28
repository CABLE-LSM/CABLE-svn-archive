load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; Procedure Code ==================================================
; created by Lauren Stevens 11-21 April 2011
; during the NCL cource at BMTC.
; Modified 21 May 2012
; =================================================================

procedure MmnPanel(Dir:string,varName:string,Yr:string,pltType:string,pltColor:string,\
                        field:numeric,region:logical,smthFlag:logical,cnLevels:numeric,\
                        lmask:logical,jpeg:logical)

begin

 if (region) then
  regno = 4
 end if

  pltName = Dir + "/" + varName + "_" + Yr + "yrs"

  wks = gsn_open_wks(pltType,pltName)

; Color Map ===========================================================

cmap = (/(/1.00, 1.00, 1.00/), (/0.00, 0.00, 0.00/), \
          (/.76,.64,.83/), \
          (/.84,.89,.93/), \
          (/.71,.79,1.0/), \
          (/.56,.70,1.0/), \
          (/.50,.59,1.0/), \
          (/.67,.81,.39/), \
          (/.91,.96,.62/), \
          (/1.0,.98,.08/), \
          (/1.0,.82,.13/), \
          (/1.0,.64,.04/), \
          (/1.0,.30,0.0/), \
          (/.73,0.0,.21/) /)

  ;gsn_define_colormap(wks,cmap)
  gsn_define_colormap(wks,pltColor)
  ;gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

; Resources ========================================================

  res = True
  
  res@gsnFrame             = False             ; Don't advance the frame
  res@gsnDraw              = False             ; only for Panel plots
                                               ; switch off for individual plots

  res@cnFillOn             = True              ; Turn on contour fill.
  res@cnFillMode           = "RasterFill"      ; Pixels
  res@cnLineLabelsOn       = False             ; Turn off line labels.
  res@cnInfoLabelOn        = False             ; Turn off info label.
if (lmask) then
  res@cnFillDrawOrder      = "PreDraw"         ; Predraw cn so mask will appear
                                               ; (if land transparent ow land will also be masked)
end if
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

  if (region) then
   if (regno.eq.1) then
    Area = "Global"
    latS = -90
    latN =  90 
    lonL =   0
    lonR = 360
   end if 
   if (regno.eq.2) then
    Area = "Australia"
    latS = -45
    latN = -10 
    lonL = 110 
    lonR = 160
   end if 
   if (regno.eq.3) then
    Area = "NorthAmerica"
    latS = 14
    latN = 85
    lonL = 190
    lonR = 308
   end if
   if (regno.eq.4) then
    Area = "45degN"
    latS = 0
    latN = 90
    lonL = 0
    lonR = 360
   end if
  else
   Area = "Global"
   latS = -90
   latN =  90
   lonL =   0
   lonR = 360
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
  
if (lmask) then
  ;res@mpAreaMaskingOn      = True
  ;res@mpMaskAreaSpecifiers = "Water"
  ;res@mpGridMaskMode       = "MaskOcean"
  res@mpFillOn             = True
  res@mpOceanFillColor     = "white"
  res@mpLandFillColor      = "transparent"
end if

;titles -----------------------------------------------------------------  
  
  res@tiMainOn      = False
  ;res@tiXAxisString = "Longitude"
  ;res@tiYAxisString = "Latitude"

; Plot ------------------------------------------------------------------
       
  plot  = new(12,graphic)
    
do i =0,11,1
  plot(i)  = gsn_csm_contour_map(wks,field(i,{latS:latN},{lonL:lonR}),res)
end do

; Paneled Plot ====================================================

  pres                       = True
  pres@gsnMaximize           = True
  pres@gsnOrientation        = "landscape"
  pres@gsnPanelLabelBar      = True
  pres@gsnPanelFigureStrings = (/"JAN","FEB","MAR","APR","MAY","JUN",\
                                      "JUL","AUG","SEP","OCT","NOV","DEC"/)
  pres@gsnPanelRowSpec       = True
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@pmLabelBarWidthF      = 0.8         ; makes label bar bigger
  pres@txString              = (field@long_name + " (" + field@units + ")")
  
  gsn_panel(wks,plot,(/4,4,4/),pres)
  ;gsn_panel(wks,plot,(/3,3,3,3/),pres)

 if(jpeg) then
  system("ps2epsi "+pltName+".ps "+pltName+".eps")
  ;system("convert -density 500 -rotate 270 "+pltName+".eps "+pltName+".jpg")
  ;system("convert -density 1000 -rotate 270 "+pltName+".ps "+pltName+".jpg")
  ;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
 end if

end

