load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


; load npp and rsoil monthly time series from CABLE

   in = addfile("/home/cmar/zie022/access/analysis/rsoil_18yrs.nc","r")
   soil_d = in->field1523

   in = addfile("/home/cmar/zie022/access/analysis/npp_18yrs.nc","r")
   npp_d = in->field1389

; load area for each grid cell

   in = addfile("/home/cmar/zie022/access/analysis/grid_area.nc","r")
   grid_area = in->grid_area

; load land faction 

   in = addfile("/home/cmar/zie022/access/analysis/land_fraction.nc","r")
   lfrac = in->lsm

; need to convert from double to float
   npp = dble2flt(npp_d)
   soil = dble2flt(soil_d)

; conversion factor (kgC/s -> PgC/yr)
;   c = 60. * 60. * 24. * 30. * 12. / 1e12
   c = 60. * 60. * 24. * 365. / 1e12

; 18 years of data, monthly values
   t = 12*18 
   y = 18

   npp_m  = new(t,float)
   nep_m  = new(t,float)
   rso_m  = new(t,float)
  
; monthly global mean
   do j = 0,t-1
      rso_m(j)  = sum(soil(j,0,:,:) * grid_area * c) ; PgC per yr
      npp_m(j)  = sum(npp(j,0,:,:) * grid_area * c) ; PgC per yr
      nep_m(j)  = sum( (npp(j,0,:,:)-soil(j,0,:,:)) * grid_area * c) ; PgC per yr
   end do


; annual global mean 
   rso_y = new((/y/),float)
   npp_y = new((/y/),float)
   nep_y = new((/y/),float)

   rso_ym = new((/y,12/),float)
   npp_ym = new((/y,12/),float)
   nep_ym = new((/y,12/),float)


   do j = 0,y-1
      rso_y(j) = sum( rso_m(j*12:(j+1)*12-1) ) / 12
      npp_y(j) = sum( npp_m(j*12:(j+1)*12-1) ) / 12
      nep_y(j) = sum( nep_m(j*12:(j+1)*12-1) ) / 12

; for plot only
      rso_ym(j,:) = rso_m(j*12:(j+1)*12-1)
      npp_ym(j,:) = npp_m(j*12:(j+1)*12-1)
      nep_ym(j,:) = nep_m(j*12:(j+1)*12-1)
   end do

   print(npp_y)
   print(rso_y)
   print(nep_y)
   print(avg(npp_y(0:17)))
   print(avg(rso_y(0:17)))
   print(avg(nep_y(0:17)))

; plots (global setting)
  res = True
  res@gsnDraw  = False                          ; Don't draw plots
  res@gsnFrame = False                          ; Don't advance frames


  res@tiXAxisString   = "Months"


;  generate x vector (time)
   x = ispan(1,12,1)
 
  wks = gsn_open_wks("ps","plot_flux_global")

  plots = new(4,graphic)

  plots(0) = gsn_xy(wks,x,npp_ym,res)
  plots(1) = gsn_xy(wks,x,rso_ym,res)
  plots(2) = gsn_xy(wks,x,nep_ym,res)
 

   pres                  = True         ; Set panel resources.
   pres@gsnMaximize      = True         ; Maximize plots in panel.
   pres@txString         = "Flux [PgC/yr] global mean"

   gsn_panel(wks,plots,(/2,2/),pres)


end
