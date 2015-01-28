; Procedure Code ==================================================
; created by Lauren Stevens 9-11 May 2012
; during the NCL cource at BMTC.
; =================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
usr = systemfunc("echo $USERID")
;loadscript("~"+usr+"/umplot/regrid.ncl")
; =================================================================

begin

numYr = systemfunc("echo $YR") 
cDir  = systemfunc("echo $CABLE")
mDir  = systemfunc("echo $MOSES")
ObsDir = "~"+usr+"/umplot/obs/"
cname = systemfunc("echo $CNM")

if ( cname .eq. "A" ) then
 rtype = "amip"
 a14 = False
end if
if ( cname.eq. "H" ) then
 rtype = "hist"
 a14 = True
end if
if ( cname.eq. "PI" ) then
 rtype = "pi"
 a14 = True
end if

moses  = addfile(mDir+"/seasonal_means_"+numYr+"yrs.nc","r")
cable  = addfile(cDir+"/seasonal_means_"+numYr+"yrs.nc","r")

if (isfilevar(moses,"field203")) then
 swd_m = moses  ->field203 ;rsus
else
 swd_m = moses  ->rsds
end if
if (isfilevar(cable,"field203")) then
 swd_c = cable  ->field203 ;rsus
else
 swd_c = cable  ->rsds
end if
if (isfilevar(moses,"solar")) then
sw_m  = moses  ->solar    ;rsds
else
sw_m  = moses  ->rss
end if
sw_c  = cable  ->solar    ;rsds
;printVarSummary(sw_m)
sw_c@_FillValue = -9999
sw_m    = where(sw_m.ne.0,sw_m,sw_m@_FillValue)
sw_c    = where(sw_c.ne.0,sw_c,sw_c@_FillValue)
swd_m   = where(swd_m.ne.0,swd_m,swd_m@_FillValue)
swd_c   = where(swd_c.ne.0,swd_c,swd_c@_FillValue)
sup_m = swd_m - sw_m
sup_c = swd_c - sw_c
copy_VarMeta(sw_m,sup_m)
copy_VarMeta(sw_c,sup_c)

;solar@long_name = "RSS: Net Down Surface SW Flux: SW TS"
;solar_u@long_name = "RSUS: Surface Upwelling Shortwave Radiation"
;sw@long_name = "RSDS: Surface Downwelling Shortwave Radiation"

tas_m = sup_m/swd_m
tas_c = sup_c/swd_c
dimc = dimsizes(dimsizes(tas_c))
dimm = dimsizes(dimsizes(tas_m))
if (dimc.eq.4.and.dimm.eq.4) then
 tas_l = tas_c(:,0,:,:) - tas_m(:,0,:,:)
else
 tas_l = tas_c(:,0,:,:) - tas_m
end if
copy_VarMeta(sw_m,tas_m)
copy_VarMeta(sw_c,tas_c)
copy_VarMeta(sw_c(:,0,:,:),tas_l)
;salb@long_name = "Surface Albedo as Calculated by method in UMui"

tas_m@long_name    = "Surface Albedo"
;tas_m@units        = ""
tas_c@long_name    = "Surface Albedo"
;tas_c@units        = ""
tas_l@long_name    = "Surface Albedo"

;WKS===================================================================

region   = False
;cnLevels = (/.02,.04,.06,.08,0.1,.12,.14,.16,.18,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
;cnLevels = fspan(0.1,.9,9)
lmask    = True ;False
jpeg     = False

pltName = "seas_salb_"+rtype+"6p_"+numYr+"yrs" ;Dir + "/" + varName + "_" + Yr + "yrs"
pltType = "ps"

wks = gsn_open_wks(pltType,pltName)

; Color Map ===========================================================

cmap_bT = RGBtoCmap("/home/cmar/ste69f/umplot/colormap/BlueOrangeRed.rgb")


gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;gsn_define_colormap(wks,cmap_bT)

; Resources ========================================================

  res = True
  
  res@gsnFrame             = False             ; Don't advance the frame
  res@gsnDraw              = False             ; only for Panel plots
                                               ; switch off for individual plots

  res@cnFillOn             = True              ; Turn on contour fill.
  res@cnFillMode           = "RasterFill"      ; Pixels
;  res@cnRasterSmoothingOn  = True
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

  res@cnLevelSelectionMode = "ExplicitLevels"  ; Define your own
  ;res@cnLevels = cnLevels 

;  gsn_reverse_colormap(wks)

;printMinMax(tas_l,True)
;printMinMax(tas_l(0,0,:,:),True)
;printMinMax(tas_l(2,0,:,:),True)

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
  res@mpFillOn             = True
  res@mpOceanFillColor     = "white"
  res@mpLandFillColor      = "transparent"
  end if  

;titles -----------------------------------------------------------------  
  
  res@tiMainOn      = False
  res@gsnRightStringFontHeightF = .021
  res@gsnLeftStringFontHeightF = .021

  res1 = res
  
; Plot 1 -----------------------------------------------------------------
       
  print("")
  print("Plotting AMOJ 6 Panel Albedo")
  print("")
  plot1 = new(4,graphic)
    
  cnLevels = (/.02,.04,.06,.08,0.1,.12,.14,.16,.18,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
  res@cnLevels = cnLevels
  res@lbLabelFontHeightF = .010

  res@tmXBLabelsOn = False
  res@gsnRightString      = "DJF"
if (a14) then
  res@gsnLeftString       = "ACCESS1.3 "+rtype
else
  res@gsnLeftString       = "ACCESS1.0 "+rtype
end if
  if (dimm.eq.4) then
  plot1(0) = gsn_csm_contour_map(wks,tas_m(0,0,{latS:latN},{lonL:lonR}),res)
  else
  plot1(0) = gsn_csm_contour_map(wks,tas_m(0,{latS:latN},{lonL:lonR}),res)
  end if
  res@tmYLLabelsOn = False
  res@tmXBLabelsOn = False
  res@gsnRightString      = "JJA"
if (a14) then
  res@gsnLeftString       = "ACCESS1.3 "+rtype
else
  res@gsnLeftString       = "ACCESS1.0 "+rtype
end if
  if (dimm.eq.4) then
  plot1(1) = gsn_csm_contour_map(wks,tas_m(2,0,{latS:latN},{lonL:lonR}),res)
  else
  plot1(1) = gsn_csm_contour_map(wks,tas_m(2,{latS:latN},{lonL:lonR}),res)
  end if
  res@tmYLLabelsOn = True
  res@tmXBLabelsOn = True
  res@gsnRightString      = "DJF"
if (a14) then
  res@gsnLeftString       = "ACCESS1.4 "+rtype
else
  res@gsnLeftString       = "ACCESS1.1 "+rtype
end if
  plot1(2) = gsn_csm_contour_map(wks,tas_c(0,0,{latS:latN},{lonL:lonR}),res)
  res@tmYLLabelsOn = False
  res@gsnRightString      = "JJA"
if (a14) then
  res@gsnLeftString       = "ACCESS1.4 "+rtype
else
  res@gsnLeftString       = "ACCESS1.1 "+rtype
end if
  plot1(3) = gsn_csm_contour_map(wks,tas_c(2,0,{latS:latN},{lonL:lonR}),res)

; Paneled Plot 1 ====================================================

  pres1                      = True
  pres1@gsnFrame             = False
  pres1@gsnMaximize           = True
  ;pres1@gsnOrientation        = "landscape"
  pres1@gsnPanelLabelBar      = True
  pres1@lbLabelFontHeightF    = .011
  pres1@gsnPanelYWhiteSpacePercent = 5
  pres1@gsnPanelFigureStrings = (/"a","b","c","d"/)
  pres1@gsnPanelRowSpec       = True
  pres1@pmLabelBarWidthF      = 0.7         ; makes label bar bigger
  ;pres1@txString              = (tas_m@long_name + " (" + tas_m@units + ")")
  pres1@txString              = "Seasonal Mean "+(tas_m@long_name );+ " (" + tas_m@units + ")")
  pres1@gsnPanelBottom        = 0.4885;532;563
  
  gsn_panel(wks,plot1,(/2,2/),pres1)

; Plot 2 -----------------------------------------------------------------

  gsn_define_colormap(wks,cmap_bT)
  ;gsn_reverse_colormap(wks)

  plot2 = new(2,graphic)
  ;res1@cnLevels = (/-.3,-.2,-.1,-.05,-.025,-0.01,0,0.01,.025,.05,.1,.2,.3/)
  res1@cnLevels = (/-.4,-.3,-.2,-.1,-.05,-.025,0,.025,.05,.1,.2,.3,.4/)
  ;res1@cnLevels = (/-.4,-.2,-.1,-.08,-.06,-.04,-.02,0,.02,.04,.06,.08,.1,.2,.4/)
  ;res1@lbLabelFontHeightF = .010
  res1@tmXBLabelsOn = True
  res1@tmYLLabelsOn = True
  res1@gsnRightString      = "DJF"
if (a14) then
  res1@gsnLeftString       = "ACCESS1.4 - ACCESS1.3 "+rtype
else
  res1@gsnLeftString       = "ACCESS1.1 - ACCESS1.0 "+rtype
end if
  plot2(0) = gsn_csm_contour_map(wks,tas_l(0,{latS:latN},{lonL:lonR}),res1)
  res1@tmYLLabelsOn = False
  res1@gsnRightString      = "JJA"
if (a14) then
  res1@gsnLeftString       = "ACCESS1.4 - ACCESS1.3 "+rtype
else
  res1@gsnLeftString       = "ACCESS1.1 - ACCESS1.0 "+rtype
end if
  plot2(1) = gsn_csm_contour_map(wks,tas_l(2,{latS:latN},{lonL:lonR}),res1)

; Paneled Plot 2 ====================================================

  pres2                      = True
  pres2@gsnFrame             = False
  pres2@gsnMaximize           = True
  ;pres2@gsnOrientation        = "landscape"
  pres2@gsnPanelLabelBar      = True
  pres2@lbLabelFontHeightF    = .011
  pres2@gsnPanelYWhiteSpacePercent = 5
  pres2@gsnPanelFigureStrings = (/"e","f"/)
  pres2@gsnPanelRowSpec       = True
  pres2@pmLabelBarWidthF      = 0.7         ; makes label bar bigger
  ;pres2@txString              = (tas_m@long_name + " (" + tas_m@units + ")")
  pres2@gsnPanelTop           = 0.4885;532;563
  pres2@gsnPanelBottom        = 0.01

  gsn_panel(wks,plot2,(/2,1/),pres2)

; -----------------------------------------------------------------------

  frame(wks)

; -----------------------------------------------------------------------

  if(jpeg) then
  system("convert -density 500 -rotate 270 "+pltName+".ps "+pltName+".jpg")
;  system("convert -density 500 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
  end if

end
