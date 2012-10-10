load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

numYr = systemfunc("echo $YR")
resol = systemfunc("echo $RES")
vMask = systemfunc("echo $MASK")
t = 12*stringtoint(numYr)
y = stringtoint(numYr)

if (isfilepresent("mm.sresp_"+numYr+"yrs.nc")) then

; load npp and rsoil monthly time series from CABLE

   slin = addfile("mm.sresp_"+numYr+"yrs.nc","r")
   soil_d = slin->field1523

   nin = addfile("mm.npp_"+numYr+"yrs.nc","r")
   npp_d = nin->field1389

   gin = addfile("mm.gpp_"+numYr+"yrs.nc","r")
   gpp_d = gin->field1388

   pin = addfile("mm.presp_"+numYr+"yrs.nc","r")
   presp_d = pin->field1390

;if (isfilepresent("mm.lresp_"+numYr+"yrs.nc")) then
   lin = addfile("mm.lresp_"+numYr+"yrs.nc","r")
   lresp_d = lin->field1499
;end if

; load area for each grid cell

   Dir = "/home/cmar/ste69f/umplot/"
   ain = addfile(Dir+"data/areacella_fx_ACCESS1-0_historical_r0i0p0.nc","r")
   grid_area = ain->areacella
   ;ain = addfile("grid_area.nc","r")
   ;grid_area = ain->grid_area

; load PFT fractions for each tile
   Sfrac = addfile("seasonal_means_"+numYr+"yrs.nc","r")
  if (resol.eq.96) then
   frac  = Sfrac->field1893(0,0:12,:,:)
  else
   frac  = Sfrac->field1893(0,0:4,:,:)
  end if

; load land faction 

;   lfin = addfile(Dir+"data/landfrac_ACCESS_N96.nc","r")
;   ;lfin = addfile("land_fraction.nc","r")
;   lfrac = lfin->lsm

   if (resol.eq.48) then
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
;printVarSummary(land_frac)

; need to convert from double to float
   npp   = dble2flt(npp_d)
   soil  = dble2flt(soil_d)
   gpp   = dble2flt(gpp_d)
   presp = dble2flt(presp_d)
;if (isvar("lresp_d")) then
   lrsp  = dble2flt(lresp_d)
;end if
   prsp  = dble2flt(presp_d)
   gpplf = dble2flt(gpp_d)
   do ts=0,t-1,1
;if (isvar("lrsp")) then
   lrsp(ts,:,:,:)  = lrsp(ts,:,:,:)*frac(:,:,:)
;end if
   prsp(ts,:,:,:)  = presp(ts,:,:,:)/land_frac(0,:,:,:)
   gpplf(ts,:,:,:) = gpp(ts,:,:,:)/land_frac(0,:,:,:)
   end do
;if (isvar("lrsp")) then
   lresp = dim_sum_n_Wrap(lrsp,1)
   do ts=0,t-1,1
   lresp(ts,:,:) = lresp(ts,:,:)*land_frac(0,0,:,:)
   end do
;end if

; conversion factor (kgC/s -> PgC/yr)
;   c = 60. * 60. * 24. * 30. * 12. / 1e12
   c = 60. * 60. * 24. * 365. / 1e12

; 18 years of data, monthly values
;   t = 12*stringtoint(numYr) 
;   y = stringtoint(numYr)

   npp_m  = new(t,float)
   nep_m  = new(t,float)
   rso_m  = new(t,float)
   gpp_m  = new(t,float)
   rpl_m  = new(t,float)
   rlf_m  = new(t,float)
  
print("Creating Global Means for Carbon Fields")

   nep_1 = npp-soil
   nep_2 = (soil(:,0,:,:)+presp(:,0,:,:)+lresp(:,:,:))-gpp(:,0,:,:)
   ;nep_2 = (soil(:,0,:,:)+prsp(:,0,:,:)+lresp)-gpplf(:,0,:,:)
   ;nep_2 = (soil(:,0,:,:)+presp(:,0,:,:)+lresp)-gpp(:,0,:,:)

; monthly global mean
   do j = 0,t-1
      rso_m(j)  = sum(soil(j,0,:,:) * grid_area * c)                ; PgC per yr
      npp_m(j)  = sum(npp(j,0,:,:) * grid_area * c)                 ; PgC per yr
      nep_m(j)  = sum(nep_1(j,0,:,:) * grid_area * c)               ; PgC per yr
      ;nep_m(j)  = sum(nep_2(j,:,:) * grid_area * c)                 ; PgC per yr
      ;nep_m(j)  = sum((npp(j,0,:,:)-soil(j,0,:,:)) * grid_area * c) ; PgC per yr
      gpp_m(j)  = sum(gpp(j,0,:,:) * grid_area * c)                 ; PgC per yr
      rpl_m(j)  = sum(presp(j,0,:,:) * grid_area * c)               ; PgC per yr
      rlf_m(j)  = sum(lresp(j,:,:) * grid_area * c)                 ; PgC per yr
   end do


; annual global mean 
   rso_y = new((/y/),float)
   npp_y = new((/y/),float)
   nep_y = new((/y/),float)
   gpp_y = new((/y/),float)
   rpl_y = new((/y/),float)
   rlf_y = new((/y/),float)

   rso_ym = new((/y,12/),float)
   npp_ym = new((/y,12/),float)
   nep_ym = new((/y,12/),float)
   gpp_ym = new((/y,12/),float)
   rpl_ym = new((/y,12/),float)
   rlf_ym = new((/y,12/),float)

   do j = 0,y-1
      rso_y(j) = sum( rso_m(j*12:(j+1)*12-1) ) / 12
      npp_y(j) = sum( npp_m(j*12:(j+1)*12-1) ) / 12
      nep_y(j) = sum( nep_m(j*12:(j+1)*12-1) ) / 12
      gpp_y(j) = sum( gpp_m(j*12:(j+1)*12-1) ) / 12
      rpl_y(j) = sum( rpl_m(j*12:(j+1)*12-1) ) / 12
      rlf_y(j) = sum( rlf_m(j*12:(j+1)*12-1) ) / 12

; for plot only
      rso_ym(j,:) = rso_m(j*12:(j+1)*12-1)
      npp_ym(j,:) = npp_m(j*12:(j+1)*12-1)
      nep_ym(j,:) = nep_m(j*12:(j+1)*12-1)
      gpp_ym(j,:) = gpp_m(j*12:(j+1)*12-1)
      rpl_ym(j,:) = rpl_m(j*12:(j+1)*12-1)
      rlf_ym(j,:) = rlf_m(j*12:(j+1)*12-1)
   end do

   ;print(npp_y)
   ;print(rso_y)
   ;print(nep_y)
   ;print(npp_ym(5,:))
   ;print(rso_ym(5,:))
   ;print(nep_ym(5,:))
   print(avg(npp_y))
   print(avg(gpp_y))
   print(avg(rlf_y))
   print(avg(rpl_y))
   print(avg(rso_y))
   print(avg(nep_y))

;Table

vals  = (/avg(npp_y),avg(gpp_y),avg(rlf_y),avg(rpl_y),avg(rso_y),\
   (avg(npp_y)-avg(rso_y)),((avg(rlf_y)+avg(rpl_y)+avg(rso_y))-avg(gpp_y))/)

;ncr1  = (/1,2/)
;x1    = (/.205,.700/)
;y1    = (/.9,.955/)
;text1 = (/"Field","Mean"/)

ncr2  = (/7,1/)
x2    = (/.205,.450/)
y2    = (/.115,.9/)
text2 = (/"NPP","GPP","LEAF*","PLANT","SOIL","NEP = npp-soil","NEP = resp-gpp"/)

ncr3        = (/7,1/)
x3          = (/.450,.700/)
y3          = (/.115,.9/)
text_values = vals

res1               = True
res1@txFontHeightF = 0.016
;res1@gsFillColor   = (/"deepskyblue2","lightgoldenrod1","slateblue2","lightseagreen","tomato","lightgoldenrod1","tomato","lightgoldenrod1","tomato"/)
res2               = True
res2@txFontHeightF = 0.016
;res2@gsFillColor   = (/"forestgreen","deepskyblue2","darkorchid2"/)
res3               = True
;res3@gsnMaximize   = True
res3@txFontHeightF = 0.016
;res3@gsFillColor   = (//)

pltName = "carbon_table_"+numYr+"yrs"
wksT = gsn_open_wks("ps","carbon_table_"+numYr+"yrs")
print("Creating Carbon Balance Table")
;gsn_table(wksT,ncr1,x1,y1,text1,res1)
gsn_table(wksT,ncr2,x2,y2,text2,res2)
gsn_table(wksT,ncr3,x3,y3,text_values,res3)

system("ps2epsi "+pltName+".ps "+pltName+".eps")

; plots (global setting)
  res = True
  res@gsnDraw  = False                          ; Don't draw plots
  res@gsnFrame = False                          ; Don't advance frames


  res@tiXAxisString   = "Months"


;  generate x vector (time)
   x = ispan(1,12,1)
 
  wks = gsn_open_wks("ps","mon_flux_global_"+numYr+"yrs")

  plots = new(6,graphic)

  plots(0) = gsn_xy(wks,x,npp_ym,res)
  plots(1) = gsn_xy(wks,x,rso_ym,res)
  plots(2) = gsn_xy(wks,x,nep_ym,res)
  plots(3) = gsn_xy(wks,x,gpp_ym,res)
  plots(4) = gsn_xy(wks,x,rpl_ym,res)
  plots(5) = gsn_xy(wks,x,rlf_ym,res)
 

   pres                  = True         ; Set panel resources.
   pres@gsnMaximize      = True         ; Maximize plots in panel.
   pres@gsnPanelFigureStrings = (/"NPP","RSO","NEP","GPP","RPL","RLF"/)
   pres@txString         = "Flux [PgC/yr] global mean"

   gsn_panel(wks,plots,(/3,3/),pres)

end if

end

