load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

hdir = systemfunc("echo $hdir")
usr = systemfunc("echo $USERID")
numYr = systemfunc("echo $YR")
resol = systemfunc("echo $RES")
vMask = systemfunc("echo $MASK")
vers  = stringtoint(systemfunc("echo $VN"))
t = 12*stringtoint(numYr)
y = stringtoint(numYr)

if (isfilepresent("mm.sresp_"+numYr+"yrs.nc")) then

; load npp and rsoil monthly time series from CABLE

   slin = addfile("mm.sresp_"+numYr+"yrs.nc","r")
   soil_d = slin->field1523

   nin = addfile("mm.npp_"+numYr+"yrs.nc","r")
   if (isfilevar(nin,"field1389")) then
    npp_d = nin->field1389
   else
    npp_d = nin->field1389_1
   end if

   gin = addfile("mm.gpp_"+numYr+"yrs.nc","r")
   gpp_d = gin->field1388

   pin = addfile("mm.presp_"+numYr+"yrs.nc","r")
   presp_d = pin->field1390

;if (isfilepresent("mm.lresp_"+numYr+"yrs.nc")) then
   lin = addfile("mm.lresp_"+numYr+"yrs.nc","r")
   lresp_d = lin->field1499
;end if

delete(slin)
delete(nin)
delete(gin)
delete(pin)
delete(lin)

; load area for each grid cell

   Dir = hdir+"/"+usr+"/umplot/"
   if (vers.ge.85) then
   ain = addfile(Dir+"data/areas_n"+resol+"e_v85.nc","r")
   grid_area = doubletofloat(ain->areacella)
   else
   ain = addfile(Dir+"data/areacella_N"+resol+".nc","r")
   grid_area = ain->areacella
   end if
   ;ain = addfile("grid_area.nc","r")
   ;grid_area = ain->grid_area
delete(ain)

; load PFT fractions for each tile
   Sfrac = addfile("seasonal_means_"+numYr+"yrs.nc","r")
  if (resol.eq.320) then
   ;frac  = Sfrac->field1893(0,:,:,:)
   frac1 = Sfrac->field1893 ;(0,0:12,:,:)
   dim   = dimsizes(frac1)
   frac  = frac1(0,0:dim(1)-5,:,:)
   delete(frac1)
  end if
  if (resol.eq.96) then
   ;frac  = Sfrac->field1893(0,:,:,:)
   frac1 = Sfrac->field1893 ;(0,0:12,:,:)
   dim   = dimsizes(frac1)
   frac  = frac1(0,0:dim(1)-5,:,:)
   delete(frac1)
  end if
  if (resol.eq.48) then
   frac  = Sfrac->field1893(0,0:4,:,:)
  end if
delete(Sfrac)

; load land faction 

;   lfin = addfile(Dir+"data/landfrac_ACCESS_N"+resol+".nc","r")
;   ;lfin = addfile("land_fraction.nc","r")
;   lfrac = lfin->lsm

   if (resol.eq.48) then
    fland   = addfile(Dir+"data/landfrac_n"+resol+".nc","r")
   end if
   if (resol.eq.96) then
   if (vers.ge.85) then
     fland   = addfile(Dir+"data/landfrac_N"+resol+"_v85.nc","r")
   else
    if (vMask.eq.2) then
     fland   = addfile(Dir+"data/landfrac_N"+resol+".nc","r")
    else ; v1 mask
     fland   = addfile(Dir+"data/landfrac_ACCESS_N"+resol+".nc","r")
    end if
    end if
   else 
     fland   = addfile(Dir+"data/landfrac_N"+resol+".nc","r")
   end if
   lfrac = fland->lsm
   land_frac = where(lfrac.eq.0,lfrac@_FillValue,lfrac)
   land_frac@_FillValue = lfrac@_FillValue
;printVarSummary(land_frac)
delete(fland)
delete(resol)
delete(vMask)
delete(lfrac)

; need to convert from double to float
   npp   = dble2flt(npp_d)
   soil  = dble2flt(soil_d)
   gpp   = dble2flt(gpp_d)
   presp = dble2flt(presp_d)
   ;prsp  = dble2flt(presp_d)
   ;gpplf = dble2flt(gpp_d)
delete(npp_d)
delete(soil_d)
delete(gpp_d)
delete(presp_d)
   ;if (isvar("lresp_d")) then
   ;lrsp  = dble2flt(lresp_d)
   lrsp  = new(dimsizes(lresp_d),float)
   ;end if
;delete(lresp_d)
   do ts=0,t-1,1
   ;if (isvar("lrsp")) then
    lrsp(ts,0:12,:,:)  = dble2flt(lresp_d(ts,0:12,:,:))*frac(:,:,:)
    ;lrsp(ts,:,:,:)  = lrsp(ts,:,:,:)*frac(:,:,:)
   ;end if
    ;prsp(ts,:,:,:)  = presp(ts,:,:,:)/land_frac(0,:,:,:)
    ;gpplf(ts,:,:,:) = gpp(ts,:,:,:)/land_frac(0,:,:,:)
   end do
delete(lresp_d)

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
  
print("")
print("Calculating Global Means for Carbon Fields")

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
;   print(avg(npp_y))
;   print(avg(gpp_y))
;   print(avg(rlf_y))
;   print(avg(rpl_y))
;   print(avg(rso_y))
;   print(avg(nep_y))

;Table

val   = new(7,string)
vals  = new(7,float)
vals  = (/avg(npp_y),avg(gpp_y),avg(rlf_y),avg(rpl_y),avg(rso_y),\
   (avg(npp_y)-avg(rso_y)),((avg(rlf_y)+avg(rpl_y)+avg(rso_y))-avg(gpp_y))/)
val = "" + decimalPlaces(vals,3,True)

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
text_values = val;s

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
if (y.le.1) then
print("")
end if

;gsn_table(wksT,ncr1,x1,y1,text1,res1)
gsn_table(wksT,ncr2,x2,y2,text2,res2)
gsn_table(wksT,ncr3,x3,y3,text_values,res3)

delete(wksT)
system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("rm "+pltName+".ps")

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

; Interannual variability
if (y.gt.1) then
  
  print("Plotting Inter-Annual Var of Carbon")
  print("")
  res4 = res
  res4@tiXAxisString   = "Years"
  xr = ispan(1,y,1)

  wks1 = gsn_open_wks("ps","interann_carbon_"+numYr+"yrs")

  plots1 = new(6,graphic)
  plots1(0) = gsn_xy(wks1,xr,npp_y,res4)
  plots1(1) = gsn_xy(wks1,xr,rso_y,res4)
  plots1(2) = gsn_xy(wks1,xr,nep_y,res4)
  plots1(3) = gsn_xy(wks1,xr,gpp_y,res4)
  plots1(4) = gsn_xy(wks1,xr,rpl_y,res4)
  plots1(5) = gsn_xy(wks1,xr,rlf_y,res4)

  pres1                  = True         ; Set panel resources.
  pres1@gsnMaximize      = True         ; Maximize plots in panel.
  pres1@gsnPanelFigureStrings = (/"NPP","RSO","NEP","GPP","RPL","RLF"/)
  pres1@txString         = "Flux [PgC/yr] global mean"

  gsn_panel(wks1,plots1,(/3,3/),pres1)
end if

else

print("")
print("No Carbon Table Produced, as mm.var_yrs.nc Files Don't Exist")
print("")

end if

end
