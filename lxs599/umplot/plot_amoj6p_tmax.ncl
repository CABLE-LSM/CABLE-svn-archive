; Procedure Code ==================================================
; created by Lauren Stevens 9 July 2013
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

moses  = addfile(mDir+"/Tseasonal_means_"+numYr+"yrs.nc","r")
cable  = addfile(cDir+"/Tseasonal_means_"+numYr+"yrs.nc","r")
obstas = addfile(ObsDir+"ERAi_tmax_tmin.nc","r")

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

if (isfilevar(moses,"tmax")) then
tas_m = moses  ->tmax
else
tas_m = moses  ->temp ;tasmax
end if
if (isfilevar(cable,"tmax")) then
tas_c = cable  ->tmax
else
tas_c = cable  ->temp ;tasmax
end if
tas_l = tas_c - tas_m
;tas_l = tas_c(:,0,:,:) - tas_m
tas_o = obstas ->mx2t ;a
copy_VarMeta(tas_m,tas_l)

tas_m@long_name    = "Max Screen Temperature"; at 1.5m"
tas_m@units        = "~S~o~N~K"
tas_l@long_name    = "Max Screen Temperature"; at 1.5m"
tas_l@units        = "~S~o~N~K"
tas_c@long_name    = "Max Screen Temperature"; at 1.5m"
tas_c@units        = "~S~o~N~K"
;tas_o              = tas_o+273.15
tas_o@long_name    = "ERA Interim Max Temperature at 2m"
tas_o@units        = "~S~o~N~K"

era_t              = regrid_bilin(tas_o,tas_m)
era_t@long_name    = "ERAi Max Tscrn regridded to UM N"+resol+""
era_t@units        = "~S~o~N~K"

era_T           = new((/4,145,192/),double)
era_T(0,:,:)    = (era_t(0,:,:)+era_t(1,:,:)+era_t(11,:,:))/3
era_T(1,:,:)    = (era_t(2,:,:)+era_t(3,:,:)+era_t(4,:,:))/3
era_T(2,:,:)    = (era_t(5,:,:)+era_t(6,:,:)+era_t(7,:,:))/3
era_T(3,:,:)    = (era_t(8,:,:)+era_t(9,:,:)+era_t(10,:,:))/3
;copy_VarMeta(era_t,era_T)
era_T@long_name = "ERA Interim Maximum Temperature at 2m"
era_T@units     = "deg K"
era_T!1         = "latitude"
era_T!2         = "longitude"
era_T&latitude  = era_t&latitude
era_T&longitude = era_t&longitude

bias_tasm         = tas_m(:,0,:,:)-era_T(:,:,:)
bias_tasl         = tas_l(:,0,:,:);-era_T(:,:,:)
bias_tasc         = tas_c(:,0,:,:)-era_T(:,:,:)
bias_tas          = tas_c(:,0,:,:)-tas_m(:,0,:,:)
copy_VarMeta(era_T(:,:,:),bias_tasm)
copy_VarMeta(era_T(:,:,:),bias_tasl)
copy_VarMeta(era_T(:,:,:),bias_tasc)
copy_VarMeta(tas_c(:,0,:,:),bias_tas)
bias_tasm@long_name  = "Max Screen Temperature Bias: UM-ERAi"
bias_tasm@units      = "~S~o~N~K"
bias_tasl@long_name  = "Max Screen Temperature Bias: UM-ERAi"
bias_tasl@units      = "~S~o~N~K"
bias_tasc@long_name  = "Max Screen Temperature Bias: UM-ERAi"
bias_tasc@units      = "~S~o~N~K"
if (a14) then
bias_tas@long_name   = "Max Screen Temperature Bias: A1.4-A1.3"
else
bias_tas@long_name   = "Max Screen Temperature Bias: A1.1-A1.0"
end if
bias_tas@units       = "~S~o~N~K"

;WKS===================================================================

region   = False
;cnLevels = fspan(-3,3,13)
cnLevels = ispan(-9,9,1)
lmask    = False
jpeg     = False

pltName = "bias_tmax_"+rtype+"6p_"+numYr+"yrs" ;Dir + "/" + varName + "_" + Yr + "yrs"
pltType = "ps"

wks = gsn_open_wks(pltType,pltName)

; Color Map ===========================================================

cmap_bT = RGBtoCmap(hdir+"/"+usr+"/umplot/colormap/BlueOrangeRed.rgb")

; if (varName.eq."bias_tscrn") then
  gsn_define_colormap(wks,cmap_bT)
; else
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
;   gsn_reverse_colormap(wks)
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
  print("Plotting AMOJ 6 Panel Tmax Biases")
  print("")
  plot = new(6,graphic)
    
;  res@tmXBLabelsOn = False
;  plot(0) = gsn_csm_contour_map(wks,bias_tasm(0,{latS:latN},{lonL:lonR}),res)
;  res@tmYLLabelsOn = False
;  ;res@gsnRightString       = "ACCESS 1.1"
;  plot(1) = gsn_csm_contour_map(wks,bias_tasl(0,{latS:latN},{lonL:lonR}),res)
;  ;res@tmXBLabelsOn = True
;  ;res@tmYLLabelsOn = True
;  res@gsnRightString       = ""
;  plot(2) = gsn_csm_contour_map(wks,bias_tasc(0,{latS:latN},{lonL:lonR}),res)
;  ;res@tmYLLabelsOn = False
;  res@tmXBLabelsOn = True
;  res@tmYLLabelsOn = True
;  res@gsnRightString       = ""
;  plot(3) = gsn_csm_contour_map(wks,bias_tasm(2,{latS:latN},{lonL:lonR}),res)
;  res@tmYLLabelsOn = False
;  plot(4) = gsn_csm_contour_map(wks,bias_tasl(2,{latS:latN},{lonL:lonR}),res)
;  res@tmYLLabelsOn = False
;  plot(5) = gsn_csm_contour_map(wks,bias_tasc(2,{latS:latN},{lonL:lonR}),res)

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
if (a14) then
  res@gsnLeftString = "ACCESS1.4 "+rtype
else
  res@gsnLeftString = "ACCESS1.1 "+rtype
end if
  res@gsnRightString = "DJF"
  plot(2) = gsn_csm_contour_map(wks,bias_tasc(0,{latS:latN},{lonL:lonR}),res)
  res@tmYLLabelsOn = False
  ;res@tmXBLabelsOn = True
  ;res@tmYLLabelsOn = True
  res@gsnRightString = "JJA"
  plot(3) = gsn_csm_contour_map(wks,bias_tasc(2,{latS:latN},{lonL:lonR}),res)
  res@tmYLLabelsOn = True
  res@tmXBLabelsOn = True
if (a14) then
  res@gsnLeftString = "ACCESS1.4 - ACCESS1.3 "+rtype
else
  res@gsnLeftString = "ACCESS1.1 - ACCESS1.0 "+rtype
end if
  res@gsnRightString = "DJF"
  plot(4) = gsn_csm_contour_map(wks,bias_tas(0,{latS:latN},{lonL:lonR}),res)
  res@tmYLLabelsOn = False
  res@tmXBLabelsOn = True
  res@gsnRightString = "JJA"
  plot(5) = gsn_csm_contour_map(wks,bias_tas(2,{latS:latN},{lonL:lonR}),res)
  
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
  pres@txString              = "Seasonal Mean "+(tas_m@long_name + " Bias (K)")
  ;pres@txString              = "Seasonal Mean "+(tas_m@long_name + " Bias (" + tas_m@units + ")")
  
  gsn_panel(wks,plot,(/3,2/),pres)

  if(jpeg) then
  system("convert -density 1000 -rotate 270 "+pltName+".ps "+pltName+".jpg")
;  system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
  end if

end
