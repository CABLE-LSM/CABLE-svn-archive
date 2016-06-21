; Procedure Code ==================================================
; created by Lauren Stevens 9-11 May 2012
; during the NCL cource at BMTC.
; =================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
hdir = systemfunc("echo $hdir")
usr = systemfunc("echo $USERID")
loadscript(hdir+"/"+usr+"/umplot/regrid.ncl")
; =================================================================

begin

numYr = systemfunc("echo $YR")
cDir  = systemfunc("echo $CABLE")
mDir  = systemfunc("echo $MOSES")
cname = systemfunc("echo $CNM")
resol = systemfunc("echo $RES")
ObsDir = hdir+"/"+usr+"/umplot/obs/"

moses  = addfile(mDir+"/seasonal_means_"+numYr+"yrs.nc","r")
cable  = addfile(cDir+"/seasonal_means_"+numYr+"yrs.nc","r")
obstas = addfile(ObsDir+"ERAi_precip.nc","r")

if ( cname .eq. "A" ) then
 rtype = "amip"
 a14 = False
end if
if ( cname .eq. "H" ) then
 rtype = "hist"
 a14 = True
end if
if ( cname .eq. "PI" ) then
 rtype = "pi"
 a14 = True
end if

if (isfilevar(moses,"precip")) then
 tas_m = moses  ->precip
else
if (isfilevar(cable,"tot_precip")) then
 tas_m = moses  ->tot_precip
else
 tas_m = moses  ->pr
end if
end if
if (isfilevar(cable,"precip")) then
 tas_c = cable  ->precip
else
if (isfilevar(cable,"tot_precip")) then
 tas_c = cable  ->tot_precip
else
 tas_c = cable  ->pr
end if
end if
;tas_l = cab11  ->pr
dims = dimsizes(dimsizes(tas_m))
if (dims.eq.4) then
tas_l = tas_c(:,0,:,:) - tas_m(:,0,:,:)
copy_VarMeta(tas_m(:,0,:,:),tas_l)
else
tas_l = tas_c(:,0,:,:) - tas_m(:,:,:)
copy_VarMeta(tas_m(:,:,:),tas_l)
end if
tas_o = obstas ->a

tas_m              = tas_m*24*60*60
tas_m@long_name    = "Precipitation"
tas_m@units        = "mm/day"
tas_l              = tas_l*24*60*60
tas_l@long_name    = "Precipitation"
tas_l@units        = "mm/day"
tas_c              = tas_c*24*60*60
tas_c@long_name    = "Precipitation"
tas_c@units        = "mm/day"
tas_o@long_name    = "ERA Interim Precipitation"
tas_o@units        = "mm/day"

era_t              = regrid_bilin(tas_o,tas_m)
era_t@long_name    = "ERAi Precip regridded to UM N"+resol+""
era_t@units        = "mm/day"

if (dims.eq.4) then
bias_tasm         = tas_m(:,0,:,:)-era_t(:,:,:)
bias_tasl         = tas_l(:,:,:);-era_t(:,:,:)
bias_tasc         = tas_c(:,0,:,:)-era_t(:,:,:)
else
bias_tasm         = tas_m(:,:,:)-era_t(:,:,:)
bias_tasl         = tas_l(:,:,:);-era_t(:,:,:)
bias_tasc         = tas_c(:,0,:,:)-era_t(:,:,:)
end if
copy_VarMeta(era_t(:,:,:),bias_tasm)
;copy_VarMeta(tas_c(:,:,:),bias_tasl)
;copy_VarMeta(era_t(:,:,:),bias_tasl)
copy_VarMeta(era_t(:,:,:),bias_tasc)
bias_tasm@long_name  = "Precipitation Bias: UM-ERAi"
bias_tasm@units      = "mm/day"
bias_tasl@long_name  = "Precipitation Bias: UM-ERAi"
bias_tasl@units      = "mm/day"
bias_tasc@long_name  = "Precipitation Bias: UM-ERAi"
bias_tasc@units      = "mm/day"

;WKS===================================================================

region   = False
;cnLevels = (/-8,-7,-6,-5,-4,-3,-2,-1,-0.5,0,0.5,1,2,3,4,5,6,7,8/)
cnLevels = (/-8,-4,-2,-1,-0.5,0,0.5,1,2,4,8/)
;cnLevels = ispan(-9,9,1)
lmask    = False
jpeg     = False

pltName = "bias_precip_"+rtype+"6p_"+numYr+"yrs" ;Dir + "/" + varName + "_" + Yr + "yrs"
;pltName = "bias_precip_hist6p_27yrs" ;Dir + "/" + varName + "_" + Yr + "yrs"
pltType = "ps"

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

cmap_bT = RGBtoCmap(hdir+"/"+usr+"/umplot/colormap/BlueOrangeRed.rgb")

; if (varName.eq."bias_tscrn") then
  gsn_define_colormap(wks,"BlueDarkRed18")
  ;gsn_define_colormap(wks,cmap_bT)
; else
;  ;gsn_define_colormap(wks,cmap)
;  gsn_define_colormap(wks,pltColor)
; end if

; Resources ========================================================

  res = True
  
  res@gsnFrame             = False             ; Don't advance the frame
  res@gsnDraw              = False             ; only for Panel plots
                                               ; switch off for individual plots

  res@cnFillOn             = True              ; Turn on contour fill.
  res@cnFillMode           = "RasterFill"      ; Pixels
;if (varName.eq."austbias_tscrn".or.varName.eq."austbias_precip".or. \
;    varName.eq."aust_tscrn".or.varName.eq."aust_precip") then
;  res@cnRasterSmoothingOn  = True
;end if
  res@cnLineLabelsOn       = False             ; Turn off line labels.
  res@cnInfoLabelOn        = False             ; Turn off info label.
  if(lmask)
  res@cnFillDrawOrder      = "PreDraw"         ; Predraw cn so mask will appear
                                ; (if land transparent ow land will also be masked)
  end if
  res@cnLinesOn            = False             ; Turn off Lines ?
   
  res@gsnSpreadColors      = True              ;
  res@lbLabelBarOn         = False             ;
  res@lbLabelAutoStride    = True              ;
  
  ;res@gsnLeftString        = ""                ;
  ;res@gsnRightString       = ""                ;

;  if(sysFlag) then
;   ;print("label bar scale")
;   ;symMinMaxPlt(field,12,False,res)
;  end if

;  if(cnLevels.eq."0") then
;    res@cnLevelSelectionMode = "AutomaticLevels"
;  else
    res@cnLevelSelectionMode = "ExplicitLevels"  ; Define your own
    res@cnLevels = cnLevels 
;  end if

;  if (varName.eq."bias_precip".or.varName.eq."bias_evap" \
;  .or.varName.eq."austbias_precip".or.varName.eq."austbias_precip_bom") then
   gsn_reverse_colormap(wks)
;  end if
;  if ( varName.eq."seas_smc".or.varName.eq."seas_smc1"  \
;   .or.varName.eq."seas_smc2".or.varName.eq."seas_smc3" \
;   .or.varName.eq."seas_smc4".or.varName.eq."seas_smc5" \
;   .or.varName.eq."seas_smc6" ) then
;   gsn_reverse_colormap(wks)
;  end if

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
  
  res@tiMainOn      = False
  ;res@tiXAxisString = "Longitude"
  ;res@tiYAxisString = "Latitude"
  res@gsnRightStringFontHeightF = .02
  res@gsnLeftString        = ""
  res@gsnRightString       = ""
  ;res@gsnRightString       = "ACCESS 1.0"
  
; Plot ------------------------------------------------------------------
       
  print("")
  print("Plotting AMOJ 6 Panel Precip Biases")
  print("")
  plot = new(6,graphic)

  res@gsnLeftStringFontHeightF = .021
  res@gsnRightStringFontHeightF = .021

  res@tmXBLabelsOn = False
  res@gsnRightString = "DJF"
if (a14) then
  res@gsnLeftString = "ACCESS1.3 "+rtype
else
  res@gsnLeftString = "ACCESS1.0 "+rtype
end if
  plot(0) = gsn_csm_contour_map(wks,bias_tasm(0,{latS:latN},{lonL:lonR}),res)
  res@tmYLLabelsOn = False
  res@gsnRightString = "JJA"
  plot(1) = gsn_csm_contour_map(wks,bias_tasm(2,{latS:latN},{lonL:lonR}),res)
  ;res@tmXBLabelsOn = True
  res@tmYLLabelsOn = True
  res@gsnRightString = "DJF"
if (a14) then
  res@gsnLeftString = "ACCESS1.4 "+rtype
else
  res@gsnLeftString = "ACCESS1.1 "+rtype
end if
  plot(2) = gsn_csm_contour_map(wks,bias_tasc(0,{latS:latN},{lonL:lonR}),res)
  res@tmYLLabelsOn = False
  ;res@tmXBLabelsOn = True
  ;res@tmYLLabelsOn = True
  res@gsnRightString = "JJA"
  plot(3) = gsn_csm_contour_map(wks,bias_tasc(2,{latS:latN},{lonL:lonR}),res)
  res@tmYLLabelsOn = True
  res@tmXBLabelsOn = True
  res@gsnRightString = "DJF"
if (a14) then
  res@gsnLeftString = "ACCESS1.4 - ACCESS1.3 "+rtype
else
  res@gsnLeftString = "ACCESS1.1 - ACCESS1.0 "+rtype
end if
  plot(4) = gsn_csm_contour_map(wks,bias_tasl(0,{latS:latN},{lonL:lonR}),res)
  res@tmYLLabelsOn = False
  res@gsnRightString = "JJA"
  plot(5) = gsn_csm_contour_map(wks,bias_tasl(2,{latS:latN},{lonL:lonR}),res)
  
; Paneled Plot ====================================================

  pres                       = True
  ;pres@gsnMaximize           = True
  ;pres@gsnOrientation        = "landscape"
  pres@gsnPanelLabelBar      = True
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelFigureStrings = (/"a","b","c","d","e","f"/)
  ;pres@gsnPanelFigureStrings = (/"DJF","DJF","DJF","JJA","JJA","JJA"/)
  ;pres@gsnPanelFigureStrings = (/"DJF","MAM","JJA","SON"/)
  ;pres@gsnPanelRowSpec       = True
  pres@pmLabelBarWidthF      = 0.8         ; makes label bar bigger
  pres@txString              = "Seasonal Mean "+(tas_m@long_name + " Bias (" + tas_m@units + ")")
  
  gsn_panel(wks,plot,(/3,2/),pres)

  if(jpeg) then
  system("convert -density 1000 -rotate 270 "+pltName+".ps "+pltName+".jpg")
;  system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
  end if

end
