load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


; load soil carbon (microbial pool) per tile
   in = addfile("soil_carbon.nc","r")
; temp(t,pseudo,lat,lon)
   soil_mic = in->temp
   soil_slw = in->temp_1
   soil_pas = in->temp_2

; load PFT fractions for each tile
   in = addfile("fractions.nc","r")
   frac = in->field1391

; load land fractions for each grid cell
   in = addfile("land_fraction.nc","r")
   land_frac = in->lsm

; load area for each grid cell
   in = addfile("grid_area.nc","r")
   grid_area = in->grid_area
   
; calculate soil content for each PFT
   tiles = 17
   pft = 12
   t = 0 ; 3 months, monthly values [0 .. 2]
   lat = 145;
   lon = 192;

; conversion factor (gC?)
   c = 1.
 
; zonal mean for each PFT
   mean_soil_mic = new((/tiles,lat/),float)
   mean_soil_glb_mic = new((/tiles/),float)

   mean_soil_slw = new((/tiles,lat/),float)
   mean_soil_glb_slw = new((/tiles/),float)

   mean_soil_pas = new((/tiles,lat/),float)
   mean_soil_glb_pas = new((/tiles/),float)

   mean_soil = new((/tiles,lat/),float)


   do i = 0,tiles-1
      do k = 0,lat-1
         mean_soil_mic(i,k) =  sum( soil_mic(t,i,k,:) * frac(0,i,k,:) * land_frac(0,0,k,:) * c ) / sum( land_frac(0,0,k,:) )
         mean_soil_slw(i,k) =  sum( soil_slw(t,i,k,:) * frac(0,i,k,:) * land_frac(0,0,k,:) * c ) / sum( land_frac(0,0,k,:) )
         mean_soil_pas(i,k) =  sum( soil_pas(t,i,k,:) * frac(0,i,k,:) * land_frac(0,0,k,:) * c ) / sum( land_frac(0,0,k,:) )
;         mean_soil(i,k) = mean_soil_mic(i,k) + mean_soil_slw(i,k) + mean_soil_pas(i,k)
         mean_soil(i,k) = sum( (/mean_soil_mic(i,k),mean_soil_slw(i,k),mean_soil_pas(i,k)/) )
      end do
; global value for each PFT [PgC = GtC] and pool
      mean_soil_glb_mic(i) = sum( soil_mic(t,i,:,:) * frac(0,i,:,:) * land_frac(0,0,:,:) * grid_area(:,:) * c ) / 1e15
      mean_soil_glb_slw(i) = sum( soil_slw(t,i,:,:) * frac(0,i,:,:) * land_frac(0,0,:,:) * grid_area(:,:) * c ) / 1e15
      mean_soil_glb_pas(i) = sum( soil_pas(t,i,:,:) * frac(0,i,:,:) * land_frac(0,0,:,:) * grid_area(:,:) * c ) / 1e15
   end do

 
; zonal mean total (all PFTs per pool)
   mean_soil_tot_mic = new((/lat/),float)
   mean_soil_tot_slw = new((/lat/),float)
   mean_soil_tot_pas = new((/lat/),float)
   mean_soil_tot = new((/lat/),float)

   do k = 0,lat-1
      mean_soil_tot_mic(k) = sum( mean_soil_mic(:,k) )
      mean_soil_tot_slw(k) = sum( mean_soil_slw(:,k) )
      mean_soil_tot_pas(k) = sum( mean_soil_pas(:,k) )
;      mean_soil_tot(k) = sum(mean_soil_mic(:,k)) + sum(mean_soil_slw(:,k)) + sum(mean_soil_pas(:,k))
      mean_soil_tot(k) = sum( (/mean_soil_tot_mic(k),mean_soil_tot_slw(k),mean_soil_tot_pas(k)/) )
   end do

; global value for each pool and total
   print( sum(mean_soil_glb_mic) ) ; total carbon in microbial pool
   print( sum(mean_soil_glb_slw) ) ; total carbon in slow pool
   print( sum(mean_soil_glb_pas) ) ; total carbon in passive pool
   print( sum(mean_soil_glb_mic) + sum(mean_soil_glb_slw) + sum(mean_soil_glb_pas) ) ; total soil carbon


; absolute latitude
   mean_soil_abs = new((/tiles,73/),float)
   mean_soil_tot_abs = new((/73/),float)
 
   mean_soil_abs(:,0) = mean_soil(:,72)
   mean_soil_tot_abs(0) = mean_soil_tot(72)
   do i = 1,72
;      mean_soil_abs(:,i) = mean_soil(:,72-i) + mean_soil(:,72+i)
      do j = 0,tiles-1
         mean_soil_abs(j,i) = sum( (/mean_soil(j,72-i),mean_soil(j,72+i)/) )
      end do
;      mean_soil_tot_abs(i) = mean_soil_tot(72-i) + mean_soil_tot(72+i)
     mean_soil_tot_abs(i) = sum( (/mean_soil_tot(72-i),mean_soil_tot(72+i)/) )

   end do


; plots (global setting)
  res = True
  res@gsnDraw  = False                          ; Don't draw plots
  res@gsnFrame = False                          ; Don't advance frames

;  generate x vector (lat)
   x = new( (/lat/),float )
   dummy = -90.
   do i = 0,lat-1
      x(i) = dummy
      dummy = dummy + 1.25
   end do

; absolute x vector (lat)   
   x1 = x(72:144)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plot soil carbon over absolute latitude 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   wks = gsn_open_wks("ps","soil_carbon_zonal_mean_abs")
   plots = new(12,graphic)
   plots(0) = gsn_xy(wks,x1,mean_soil_tot_abs,res)

   do i = 1,11
      plots(i) = gsn_xy(wks,x1,mean_soil_abs(i-1,:),res)
   end do

   pres                  = True         ; Set panel resources.
   pres@gsnMaximize      = True         ; Maximize plots in panel.
   pres@txString         = "Soil carbon [gC/m2] zonal mean"
   gsn_panel(wks,plots,(/4,3/),pres)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plot soil carbon over latitude
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   wks = gsn_open_wks("ps","soil_carbon_zonal_mean")
   plots = new(12,graphic)
   plots(0) = gsn_xy(wks,x,mean_soil_tot,res)

   do i = 1,11
      plots(i) = gsn_xy(wks,x,mean_soil(i-1,:),res)
   end do

   pres                  = True         ; Set panel resources.
   pres@gsnMaximize      = True         ; Maximize plots in panel.
   pres@txString         = "Soil carbon [gC/m2] zonal mean"
   gsn_panel(wks,plots,(/4,3/),pres)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plot microbial soil carbon over latitude
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   wks = gsn_open_wks("ps","soil_mic_carbon_zonal_mean")
   plots = new(12,graphic)
   plots(0) = gsn_xy(wks,x,mean_soil_tot_mic,res)

   do i = 1,11
      plots(i) = gsn_xy(wks,x,mean_soil_mic(i-1,:),res)
   end do

   pres                  = True         ; Set panel resources.
   pres@gsnMaximize      = True         ; Maximize plots in panel.
   pres@txString         = "Microbial soil carbon [gC/m2] zonal mean"
   gsn_panel(wks,plots,(/4,3/),pres)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plot slow soil carbon over latitude
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   wks = gsn_open_wks("ps","soil_slw_carbon_zonal_mean")
   plots = new(12,graphic)
   plots(0) = gsn_xy(wks,x,mean_soil_tot_slw,res)

   do i = 1,11
      plots(i) = gsn_xy(wks,x,mean_soil_slw(i-1,:),res)
   end do

   pres                  = True         ; Set panel resources.
   pres@gsnMaximize      = True         ; Maximize plots in panel.
   pres@txString         = "Slow soil carbon [gC/m2] zonal mean"
   gsn_panel(wks,plots,(/4,3/),pres)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plot passive soil carbon over latitude
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   wks = gsn_open_wks("ps","soil_pas_carbon_zonal_mean")
   plots = new(12,graphic)
   plots(0) = gsn_xy(wks,x,mean_soil_tot_pas,res)

   do i = 1,11
      plots(i) = gsn_xy(wks,x,mean_soil_pas(i-1,:),res)
   end do

   pres                  = True         ; Set panel resources.
   pres@gsnMaximize      = True         ; Maximize plots in panel.
   pres@txString         = "Passive soil carbon [gC/m2] zonal mean"
   gsn_panel(wks,plots,(/4,3/),pres)


end
