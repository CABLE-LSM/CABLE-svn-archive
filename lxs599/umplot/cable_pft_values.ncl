load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; load soil respiration per tile
   in = addfile("rsoil_tile_18yrs.nc","r")
   soil_d = in->temp_2
   soil = dble2flt(soil_d)

; load PFT fractions for each tile
   in = addfile("fractions.nc","r")
   frac = in->field1391

; load land fractions for each grid cell
   in = addfile("land_fraction.nc","r")
   land_frac = in->lsm

; load area for each grid cell
   in = addfile("grid_area.nc","r")
   grid_area = in->grid_area
   
; calculate soil respiration per month and PFT
   tiles = 17
   pft = 13 ;12
   t = 12*18
   y = 18

; conversion factor (kgC/s -> PgC/yr)
   pcpy = 60. * 60. * 24. * 365. / 1e12
   ;pcpy = 60. * 60. * 24. * 30. * 12. / 1e12
 
   mean_soil = new((/t,tiles/),float)    
 

   do j = 0,t-1
      do i = 0,tiles-1
         mean_soil(j,i) = sum( soil(j,i,:,:) * frac(0,i,:,:) * land_frac(0,0,:,:) * grid_area * pcpy)
      end do
   end do

; global values per month
   
   mean_soil_glob_m = new(t,float)
   do j = 0,t-1
      mean_soil_glob_m(j) = sum (mean_soil(j,:))
   end do
  

; global mean per year
  mean_soil_ym = new((/y,tiles,12/),float) 
  mean_soil_glob_ym = new((/y,12/),float) 
  mean_soil_glob_y = new(y,float) 

   do j = 0,y-1
      do i = 0,tiles-1
         mean_soil_ym(j,i,:) = mean_soil(j*12:(j+1)*12-1,i)
      end do
      mean_soil_glob_y(j) = sum( mean_soil(j*12:(j+1)*12-1,:) ) / 12
      mean_soil_glob_ym(j,:) = mean_soil_glob_m(j*12:(j+1)*12-1)
   end do 

   

  print(mean_soil_glob_y)

; plots (global setting)
  res = True
  res@gsnDraw  = False                          ; Don't draw plots
  res@gsnFrame = False                          ; Don't advance frames


;;;;;;;;;;;;;;;;;;
; plot global mean
;;;;;;;;;;;;;;;;;;


   x1 = ispan(1,12,1)  

   wks = gsn_open_wks("ps","plot_flux_pft")

   plots = new(12,graphic)

   plots(0) = gsn_xy(wks,x1,mean_soil_glob_ym(1:17,:),res)

   do i = 1,11
      plots(i) = gsn_xy(wks,x1,mean_soil_ym(1:17,i-1,:),res)
   end do

   pres                  = True         ; Set panel resources.
   pres@gsnMaximize      = True         ; Maximize plots in panel.
   pres@txString         = "Soil respiration [PgC/yr] global seasonal cycle"

   gsn_panel(wks,plots,(/4,3/),pres)

  

end
