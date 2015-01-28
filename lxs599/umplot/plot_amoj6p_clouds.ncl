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
cname = systemfunc("echo $CNM")
ObsDir = "~"+usr+"/umplot/obs/"

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

moses  = addfile(mDir+"/seasonal_means_"+numYr+"yrs.nc","r")
cable  = addfile(cDir+"/seasonal_means_"+numYr+"yrs.nc","r")

;cable  = addfile(wDir+"access13/amip/clt_yseasmean.nc","r")
;cable  = addfile(wDir+"saaba_a1p1_27yrs/seasonal_means_27yrs.nc","r")
;obstas = addfile(ObsDir+"ERAi_clouds_yseasmean_invertlat.nc","r")

if (isfilevar(moses,"clt")) then
 tas_m = moses  ->clt
else
 tas_m = moses  ->field30_1
end if
if (max(tas_m).le.2.0) then
 tas_m = tas_m*100.
end if
if (isfilevar(cable,"clt")) then
 tas_c = cable  ->clt
else
 tas_c = cable  ->field30_1
 tas_c = tas_c*100.
end if
;tas_o = obstas ->tcc
dims = dimsizes(dimsizes(tas_m))
if (dims.eq.4) then
tas_l = tas_c-tas_m
else
tas_l = tas_c(:,0,:,:)-tas_m(:,:,:)
;tas_l = tas_c(:,0,:,:)-tas_m
end if
copy_VarMeta(tas_m,tas_l)

;printVarSummary(tas_l)  ;[t | 4] x [unspecified_1 | 1] x [latitude | 145] x [longitude | 192]
;print("nb DJF: "+tas_l(0,0,117,139)+"")
;print("bo DJF: "+tas_l(0,0,104,145)+"")
;print("hy DJF: "+tas_l(0,0,121,13)+"")
;print("th DJF: "+tas_l(0,0,113,7)+"")
;print("lw DJF: "+tas_l(0,0,100,140)+"")
;print("wb DJF: "+tas_l(0,0,101,147)+"")
;print("lo DJF: "+tas_l(0,0,114,4)+"")
;print("hv DJF: "+tas_l(0,0,106,153)+"")
;print("")
;print("nb JJA: "+tas_l(2,0,117,139)+"")
;print("bo JJA: "+tas_l(2,0,104,145)+"")
;print("hy JJA: "+tas_l(2,0,121,13)+"")
;print("th JJA: "+tas_l(2,0,113,7)+"")
;print("lw JJA: "+tas_l(2,0,100,140)+"")
;print("wb JJA: "+tas_l(2,0,101,147)+"")
;print("lo JJA: "+tas_l(2,0,114,4)+"")
;print("hv JJA: "+tas_l(2,0,106,153)+"")

tas_m@long_name    = "Total Cloud"
tas_m@units        = "%"
tas_c@long_name    = "Total Cloud"
tas_c@units        = "%"
;tas_o              = tas_o*100
;tas_o@long_name    = "ERA Interim Total Cloud"
;tas_o@units        = "%"
tas_l@long_name    = "Total Cloud"
tas_l@units        = "%"

;era_t              = regrid_bilin(tas_o,tas_m)
;era_t@long_name    = "ERAi Total Cloud regridded to UM N96"
;era_t@units        = "%"

;bias_tasm         = tas_m(:,:,:)-era_t(:,:,:)
;bias_tasc         = tas_c(:,0,:,:)-era_t(:,:,:)
;copy_VarMeta(era_t(:,:,:),bias_tasm)
;copy_VarMeta(era_t(:,:,:),bias_tasc)
;bias_tasm@long_name  = "Total Cloud Bias: UM-ERAi"
;bias_tasm@units      = "%"
;bias_tasc@long_name  = "Total Cloud Bias: UM-ERAi"
;bias_tasc@units      = "%"

;WKS===================================================================

region   = False
;;cnLevels = (/-8,-7,-6,-5,-4,-3,-2,-1,-0.5,0,0.5,1,2,3,4,5,6,7,8/)
;;cnLevels = (/-8,-4,-2,-1,-0.5,0,0.5,1,2,4,8/)
;cnLevels = ispan(-75,75,10)
lmask    = True
jpeg     = False

pltName = "seas_cloud_"+rtype+"6p_"+numYr+"yrs"
pltType = "ps"

wks = gsn_open_wks(pltType,pltName)

; Color Map ===========================================================

cmap_bT = RGBtoCmap("/home/cmar/ste69f/umplot/colormap/BlueOrangeRed.rgb")
;gsn_define_colormap(wks,"spread_15lev")
;gsn_define_colormap(wks,"rh_19lev")
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
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

  gsn_reverse_colormap(wks)

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
  res@gsnLeftStringFontHeightF = .021
  res@gsnRightStringFontHeightF = .021

  res1 = res
  
; Plot 1 -----------------------------------------------------------------

  print("")
  print("Plotting AMOJ 6 Panel Cloud Fraction")
  print("")
       
  plot1 = new(4,graphic)
    
  res@cnLevels = (/1,2,3,4,5,6,7,8,9,10,15,25,35,45,55,65,75,85,95/)
  ;res@cnLevels = (/2,4,6,8,10,12,14,16,18,20,25,35,45,65,75,85,95/)
  res@tmXBLabelsOn = False
if (a14) then
  res@gsnLeftString       = "ACCESS1.3 "+rtype
else
  res@gsnLeftString       = "ACCESS1.0 "+rtype
end if
  res@gsnRightString      = "DJF"
  if (dims.eq.4) then
  plot1(0) = gsn_csm_contour_map(wks,tas_m(0,0,{latS:latN},{lonL:lonR}),res)
  else
  plot1(0) = gsn_csm_contour_map(wks,tas_m(0,{latS:latN},{lonL:lonR}),res)
  end if
  res@tmXBLabelsOn = False
  res@tmYLLabelsOn = False
  res@gsnRightString      = "JJA"
if (a14) then
  res@gsnLeftString       = "ACCESS1.3 "+rtype
else
  res@gsnLeftString       = "ACCESS1.0 "+rtype
end if
  if (dims.eq.4) then
  plot1(1) = gsn_csm_contour_map(wks,tas_m(2,0,{latS:latN},{lonL:lonR}),res)
  else
  plot1(1) = gsn_csm_contour_map(wks,tas_m(2,{latS:latN},{lonL:lonR}),res)
  end if
  res@tmYLLabelsOn = True
  res@tmXBLabelsOn = True
if (a14) then
  res@gsnLeftString       = "ACCESS1.4 "+rtype
else
  res@gsnLeftString       = "ACCESS1.1 "+rtype
end if
  res@gsnRightString      = "DJF"
  ;plot1(2) = gsn_csm_contour_map(wks,tas_c(0,{latS:latN},{lonL:lonR}),res)
  plot1(2) = gsn_csm_contour_map(wks,tas_c(0,0,{latS:latN},{lonL:lonR}),res)
  res@tmYLLabelsOn = False
  res@gsnRightString      = "JJA"
if (a14) then
  res@gsnLeftString       = "ACCESS1.4 "+rtype
else
  res@gsnLeftString       = "ACCESS1.1 "+rtype
end if
  ;plot1(3) = gsn_csm_contour_map(wks,tas_c(2,{latS:latN},{lonL:lonR}),res)
  plot1(3) = gsn_csm_contour_map(wks,tas_c(2,0,{latS:latN},{lonL:lonR}),res)

; Paneled Plot 1 ====================================================

  pres1                      = True
  pres1@gsnFrame             = False
  pres1@gsnMaximize           = True
  ;pres1@gsnOrientation        = "landscape"
  pres1@gsnPanelLabelBar      = True
  pres1@gsnPanelYWhiteSpacePercent = 5
  pres1@gsnPanelFigureStrings = (/"a","b","c","d"/)
  pres1@gsnPanelRowSpec       = True
  pres1@pmLabelBarWidthF      = 0.7         ; makes label bar bigger
  ;pres1@txString              = (tas_m@long_name + " (" + tas_m@units + ")")
  pres1@txString              = "Seasonal Mean "+(tas_m@long_name + " (" + tas_m@units + ")")
  pres1@gsnPanelBottom        = 0.4885;532;563
  
  gsn_panel(wks,plot1,(/2,2/),pres1)

; Plot 2 -----------------------------------------------------------------

  gsn_define_colormap(wks,cmap_bT)
  gsn_reverse_colormap(wks)

  plot2 = new(2,graphic)
  res1@cnLevels = (/-35,-25,-15,-10,-8,-6,-4,-2,0,2,4,6,8,10,15,25,35/)
  res1@tmXBLabelsOn = True
  res1@tmYLLabelsOn = True
  res1@gsnRightString      = "DJF"
if (a14) then
  res1@gsnLeftString       = "ACCESS1.4 - ACCESS1.3 "+rtype
else
  res1@gsnLeftString       = "ACCESS1.1 - ACCESS1.0 "+rtype
end if
  if (dims.eq.4) then
  plot2(0) = gsn_csm_contour_map(wks,tas_l(0,0,{latS:latN},{lonL:lonR}),res1)
  else
  plot2(0) = gsn_csm_contour_map(wks,tas_l(0,{latS:latN},{lonL:lonR}),res1)
  end if
  res1@tmYLLabelsOn = False
  res1@gsnRightString      = "JJA"
if (a14) then
  res1@gsnLeftString       = "ACCESS1.4 - ACCESS1.3 "+rtype
else
  res1@gsnLeftString       = "ACCESS1.1 - ACCESS1.0 "+rtype
end if
  if (dims.eq.4) then
  plot2(1) = gsn_csm_contour_map(wks,tas_l(2,0,{latS:latN},{lonL:lonR}),res1)
  else
  plot2(1) = gsn_csm_contour_map(wks,tas_l(2,{latS:latN},{lonL:lonR}),res1)
  end if

; Paneled Plot 2 ====================================================

  pres2                      = True
  pres2@gsnFrame             = False
  pres2@gsnMaximize           = True
  ;pres2@gsnOrientation        = "landscape"
  pres2@gsnPanelLabelBar      = True
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
  end if

end
