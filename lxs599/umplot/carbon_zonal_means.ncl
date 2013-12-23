load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;procedure zonal_means(fieldname:string,varname:string)

begin

;var_name  = (/"gpp","npp","presp","lresp","sresp","co2atm","co2ppm"/)
;var_field = (/"field1388","field1389","field1390","field1499",     \
;              "field1523","field1563","field1564"/)
var_name  = (/"gpp","npp","presp","sresp"/)
;var_field = (/"gpp","npp","ra","rh"/)
var_field = (/"field1388","field1389","field1390","field1523"/)
; note lresp has lat,lon,13,time 
; co2ppm has 38 levels
; co2atmos has missing values ?

usr = systemfunc("echo $USERID")
numYr    = systemfunc("echo $YR")
filDir   = systemfunc("echo $PWD")
;filDir   = systemfunc("echo $DIRW")
resol    = systemfunc("echo $RES")
vMask    = systemfunc("echo $MASK")
jpg      = systemfunc("echo $CHFMT")
if(jpg.eq."y") then
jpeg = True
else
jpeg = False
end if

Dir  = "~"+usr+"/umplot/"

print("Calculating Zonal Means for Carbon Fields")
print("")

do m=0,dimsizes(var_name)-1,1

; load soil respiration (rsoil,sresp,field1523,rh)
   filin  = addfile(filDir+"/mm."+var_name(m)+"_"+numYr+"yrs.nc","r")
   if (m.eq.1.and.isfilevar(filin,"field1389_1")) then
    var_field(1) = "field1389_1"
   end if
   soil_d = filin->$var_field(m)$

if (var_name(m).eq."gpp") then
   filobs = addfile(Dir+"obs/GPP_1982-2008.nc","r")
   ogpp   = filobs->GPP
end if

; convert double to float
   soil = dble2flt(soil_d)

; load PFT fractions for each tile
;   if (resol.eq.48) then
;   Sfrac = addfile(Dir+"data/fractions.nc","r")
;   else
;   if (vMask.eq.2) then
;   Sfrac = addfile(Dir+"data/fractions.nc","r")
;   else
;   Sfrac = addfile(Dir+"data/fractions.nc","r")
;   end if
;   end if
;   Sfrac = addfile(filDir+"/seasonal_means_"+numYr+"yrs.nc","r")
;   frac  = Sfrac->field1893(0,:,:,:) ; field1391

; load land fractions for each grid cell
   if (resol.eq.48) then
   ; Find n48 landfrac
   fland   = addfile(Dir+"data/landfrac_n48.nc","r")
   else ; n96
   if (vMask.eq.2) then
   fland   = addfile(Dir+"data/landfrac_v2.nc","r")
   else ; v1 mask
   fland   = addfile(Dir+"data/landfrac_ACCESS_N96.nc","r")
   end if
   end if
   lfrac = fland->lsm
   land_frac = where(lfrac.eq.0,lfrac@_FillValue,lfrac)
   land_frac@_FillValue = lfrac@_FillValue

; load area for each grid cell
;   ;garea     = addfile(Dir+"data/grid_area.nc","r") ; Tilo's calc grid_area
;   ;grid_area = garea->grid_area
;   garea     = addfile(Dir+"data/areacella_fx_ACCESS1-0_historical_r0i0p0.nc","r")
;   grid_area = garea->areacella

   
; calculate soil respiration
   pft = 13 ;12
   dims = dimsizes(soil)
   nlat = dims(2) ; 145
   nlon = dims(3) ; 192
   ntim = dims(0) ; 12*numYr
   nyr  = ntim/12 ; $YR or numYr

if (var_name(m).eq."gpp") then
   odim = dimsizes(ogpp)
   olat = odim(1)
   olon = odim(2)
   otim = odim(0)
   oyr  = otim/12
end if

; conversion factor (kgC/s -> gC/yr)
;   gcpy = 60. * 60. * 24. * 30. * 12. * 1000.
   gcpy = 60. * 60. * 24. * 365. * 1000.
   pyr  = 365.25
 
; zonal monthly mean (3 different way to calculate)
   mean_soil_1 = new((/ntim,nlat/),float)
   mean_soil_2 = new((/ntim,nlat/),float)
   mean_soil_3 = new((/ntim,nlat/),float)
   mean_soil_4 = new((/ntim,nlat/),float)
   mean_soil_1@_FillValue = -999 ; soil@_FillValue
   mean_soil_2@_FillValue = -999 ; soil@_FillValue
   mean_soil_3@_FillValue = -999 ; soil@_FillValue
   mean_soil_4@_FillValue = -999 ; soil@_FillValue

; Les can use  zfield = zonalAve(field)
;   mean1 = new((/ntim,nlat,nlon/),float)
;   mean1@_FillValue = -999
;
;   do l = 0,ntim-1
;      do j = 0,nlat-1
;         do i = 0,nlon-1,1
;           if (land_frac(0,0,j,i).ne.0) then
;             mean1(l,j,i) = soil(l,0,j,i)/land_frac(0,0,j,i)
;           else
;             mean1(l,j,i) = -999
;           end if
;         end do
;      end do
;   end do

   do l = 0,ntim-1
      do j = 0,nlat-1

        ;if (ismissing(sum(soil(l,j,:))) .or. ismissing(sum(land_frac(0,0,j,:)))) then
        if (ismissing(sum(soil(l,0,j,:)))) then
           mean_soil_1(l,j) = -999 ; soil@_FillValue
           mean_soil_2(l,j) = -999 ; soil@_FillValue
        else 
           ;mean_soil_1(l,j) = avg( mean1(l,j,:) * gcpy)
           mean_soil_1(l,j) = avg( soil(l,0,j,:) / land_frac(0,0,j,:)  * gcpy) 
           mean_soil_2(l,j) = avg( soil(l,0,j,:) * gcpy) 
        end if

        ;if (sum(land_frac(0,0,j,:)).ne.0) then
          mean_soil_3(l,j) = sum(soil(l,0,j,:)*gcpy) / sum(land_frac(0,0,j,:))
          mean_soil_4(l,j) = sum(soil(l,0,j,:)*land_frac(0,0,j,:)*gcpy) / sum(land_frac(0,0,j,:))
        ;else
        ;  mean_soil_3(l,j) = -999
        ;  mean_soil_4(l,j) = -999
        ;end if

      end do
   end do

;  zonal annual mean 
   mean_soil_y_1 = new((/nyr,nlat/),float)
   mean_soil_y_2 = new((/nyr,nlat/),float)
   mean_soil_y_3 = new((/nyr,nlat/),float)
   mean_soil_y_4 = new((/nyr,nlat/),float)
   mean_soil_y_1@_FillValue = -999
   mean_soil_y_2@_FillValue = -999
   mean_soil_y_3@_FillValue = -999
   mean_soil_y_4@_FillValue = -999
   mean_soil_y_1!0 = "time"
   mean_soil_y_1!1 = "latitude"
   mean_soil_y_1&time     = ispan(1,stringtoint(numYr),1)
   mean_soil_y_1&latitude = soil&latitude
   copy_VarMeta(mean_soil_y_1,mean_soil_y_2)
   copy_VarMeta(mean_soil_y_1,mean_soil_y_3)
   copy_VarMeta(mean_soil_y_1,mean_soil_y_4)

   do j = 0,nyr-1
      do k = 0,nlat-1
         mean_soil_y_1(j,k) = sum( mean_soil_1(j*12:(j+1)*12-1,k) ) / 12
         mean_soil_y_2(j,k) = sum( mean_soil_2(j*12:(j+1)*12-1,k) ) / 12
         mean_soil_y_3(j,k) = sum( mean_soil_3(j*12:(j+1)*12-1,k) ) / 12
         mean_soil_y_4(j,k) = sum( mean_soil_4(j*12:(j+1)*12-1,k) ) / 12
      end do 
   end do

; Obs =========================================

if (var_name(m).eq."gpp") then
obs_mean = dim_avg_n_Wrap(dim_avg_Wrap(ogpp),0)
obs_mean = obs_mean*pyr
end if

; plot ========================================

; plots (global setting)
  res          = True
  res@gsnDraw  = False                          ; Don't draw plots
  res@gsnFrame = False                          ; Don't advance frames

  res@xyLineColor      = "blue"
  res@xyLineThicknessF = 0.5
  res@xyDashPattern    = "SolidLine"
  res@tiXAxisString    = "Latitude"
;  res@tiYAxisString    = soil@title ;"Soil Respiration [gC/m2/yr] Zonal Mean"

;  generate x vector (lat)
  x = soil&latitude
if (var_name(m).eq."gpp") then
  ox = ogpp&latitude
end if

  wks = gsn_open_wks("ps","zonal_"+var_name(m)+"_"+numYr+"yrs")
  ;wks = gsn_open_wks("ps","zonal_gpp_"+numYr+"yrs")
  ;wks = gsn_open_wks("ps","zonal_rsoil_"+numYr+"yrs")

  plots = new(4,graphic)

  res@xyCurveDrawOrder = "PostDraw"
  plots(0) = gsn_csm_xy(wks,x,mean_soil_y_1({0:nyr},:),res)
  plots(1) = gsn_csm_xy(wks,x,mean_soil_y_2({0:nyr},:),res)
  plots(2) = gsn_csm_xy(wks,x,mean_soil_y_3({0:nyr},:),res)
  plots(3) = gsn_csm_xy(wks,x,mean_soil_y_4({0:nyr},:),res)

if (var_name(m).eq."gpp") then
  res@xyLineColor      = "red"
  res@xyLineThicknessF = 4.1
  res@xyCurveDrawOrder = "Draw"
  obsplot  = gsn_csm_xy(wks,ox,obs_mean,res)

  overlay(plots(3),obsplot)
end if

  pres             = True         ; Set panel resources.
  pres@gsnMaximize = True         ; Maximize plots in panel.
  pres@gsnPanelFigureStrings = (/"avg w /lf","avg","sum","sum w *lf"/)
  pres@txString    = soil@title+" [gC/m2/yr] Zonal Mean" 
  ;pres@txString    = soil@title+" ["+soil@units+"] Zonal Mean" 
                   ;"Soil Respiration [gC/m2/yr] Zonal Mean"

  gsn_panel(wks,plots,(/2,2/),pres)

end do

end
