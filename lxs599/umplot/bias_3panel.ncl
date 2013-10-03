; -----------------------------------------------------------------------------
; Lauren Stevens 16 Aug 2011
; 3 panel bias plot
; NCL - contour plots
; -----------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/home/cmar/ste69f/umplot/plot_3p.ncl"
load "/home/cmar/ste69f/umplot/regrid.ncl"

; Set-Up ----------------------------------------------------------------------

;procedure 3panel(pltDir:string,:string,Yr:string,pltType:string,Mstp:string)

begin

cabDir = systemfunc("echo $CABLE")
mosDir = systemfunc("echo $MOSES")
;cabsea = systemfunc("echo $cabsea")
;mossea = systemfunc("echo $mossea")
;cabtem = systemfunc("echo $cabtem")
;mostem = systemfunc("echo $mostem")
numYr  = systemfunc("echo $YR")
jpg    = systemfunc("echo $JPEG")
if(jpg.eq."y") then
jpeg = True
else
jpeg = False
end if
tname =  systemfunc("echo $tname")

ObsDir = "/home/cmar/ste69f/umplot/obs/"
LsmDir = "/home/cmar/ste69f/umplot/data/"

tscrn_obs  = addfile(ObsDir+"ERAi_tp2m.nc","r")
tmxx_obs   = addfile(ObsDir+"ERAi_tmax_tmin.nc","r")

moses_file = addfile(mosDir+"/seasonal_means_"+numYr+"yrs.nc","r")
cable_file = addfile(cabDir+"/seasonal_means_"+numYr+"yrs.nc","r")
if (isfilepresent(mosDir+"/Tseasonal_means_"+numYr+"yrs.nc")) then
tmoses     = addfile(mosDir+"/Tseasonal_means_"+numYr+"yrs.nc","r")
end if
tcable     = addfile(cabDir+"/Tseasonal_means_"+numYr+"yrs.nc","r")
;moses_file = addfile(mossea,"r")
;cable_file = addfile(cabsea,"r")
;tmoses     = addfile(mostem,"r")
;tcable     = addfile(cabtem,"r")

; Variables ------------------------------------------------------------------

if (isfilevar(moses_file,"tas")) then
 tscrnM1         = moses_file ->tas
else
if (isfilevar(moses_file,"tscrn")) then
 tscrnM1         = moses_file ->tscrn
else
 tscrnM1         = moses_file ->$tname$
end if
end if
if (isfilevar(cable_file,"tscrn")) then
tscrnC1         = cable_file ->tscrn
else
tscrnC1         = cable_file ->$tname$
end if
era_t2m1        = tscrn_obs  ->a
if (isfilevar(moses_file,"tasmax")) then
 tmaxM1          = moses_file->tasmax
 tminM1          = moses_file->tasmin
else
 tmaxM1          = tmoses     ->temp
 tminM1          = tmoses     ->temp_1
end if
tmaxC1          = tcable     ->temp
tminC1          = tcable     ->temp_1
era_tmax1       = tmxx_obs   ->mx2t
era_tmin1       = tmxx_obs   ->mn2t

era_t2m1            = era_t2m1+273.15
era_t2m1@long_name  = "ERA Interim Temperature at 2m"
era_t2m1@units      = "deg K"

if (isfilevar(moses_file,"tasmax")) then
tmaxM1!2 = "longitude_1"
tmaxC1!2 = "longitude_1"
tminM1!2 = "longitude_1"
tminC1!2 = "longitude_1"
else
tmaxM1!3 = "longitude_1"
tmaxC1!3 = "longitude_1"
tminM1!3 = "longitude_1"
tminC1!3 = "longitude_1"
end if

eraT_max           = new((/4,121,240/),double)
eraT_max(0,:,:)    = (era_tmax1(0,:,:)+era_tmax1(1,:,:)+era_tmax1(11,:,:))/ 3
eraT_max(1,:,:)    = (era_tmax1(2,:,:)+era_tmax1(3,:,:)+era_tmax1(4,:,:))/ 3
eraT_max(2,:,:)    = (era_tmax1(5,:,:)+era_tmax1(6,:,:)+era_tmax1(7,:,:))/ 3
eraT_max(3,:,:)    = (era_tmax1(8,:,:)+era_tmax1(9,:,:)+era_tmax1(10,:,:))/ 3
copy_VarMeta(era_t2m1,eraT_max)
eraT_max@long_name = "ERA Interim Maximum Temperature at 2m"
eraT_max!2         = "longitude"

eraT_min = new((/4,121,240/),double)
eraT_min(0,:,:) = (era_tmin1(0,:,:)+era_tmin1(1,:,:)+era_tmin1(11,:,:))/3
eraT_min(1,:,:) = (era_tmin1(2,:,:)+era_tmin1(3,:,:)+era_tmin1(4,:,:))/3
eraT_min(2,:,:) = (era_tmin1(5,:,:)+era_tmin1(6,:,:)+era_tmin1(7,:,:))/3
eraT_min(3,:,:) = (era_tmin1(8,:,:)+era_tmin1(9,:,:)+era_tmin1(10,:,:))/3
copy_VarMeta(era_t2m1,eraT_min)
eraT_min!2 = "longitude"

era_t_rg              = regrid_bilin(era_t2m1,tscrnM1)
era_t_rg@long_name    = "ERAi Tscrn regridded to UM N96"
era_t_rg@units        = "deg K"

era_x_rg              = regrid_bilin(eraT_max,tmaxM1)
era_x_rg@long_name    = "ERAi Tmax regridded to UM N96"
era_x_rg@units        = "deg K"

era_n_rg              = regrid_bilin(eraT_min,tminM1)
era_n_rg@long_name    = "ERAi Tmin regridded to UM N96"
era_n_rg@units        = "deg K"

;tscrnM          = dim_avg_n_Wrap(tscrnM1,0)
;tscrnC          = dim_avg_n_Wrap(tscrnC1,0)
;era_t2m         = dim_avg_n_Wrap(era_t_rg,0)
;tmaxM           = dim_avg_n_Wrap(tmaxM1,0)
;tmaxC           = dim_avg_n_Wrap(tmaxC1,0)
;tminM           = dim_avg_n_Wrap(tminM1,0)
;tminC           = dim_avg_n_Wrap(tminC1,0)
;era_tmax        = dim_avg_n_Wrap(era_x_rg,0)
;era_tmin        = dim_avg_n_Wrap(era_n_rg,0)

; Biases --------------------------------------------------------------------

;Mtbias = tscrnM(0,:,:)-era_t2m
;Ctbias = tscrnC(0,:,:)-era_t2m
;Dtbias = tscrnC-tscrnM
;copy_VarMeta(tscrnM,Mtbias)
;copy_VarMeta(tscrnC,Ctbias)
;copy_VarMeta(tscrnC,Dtbias)

;Mxbias = tmaxM(0,:,:)-era_tmax
;Cxbias = tmaxC(0,:,:)-era_tmax
;Dxbias = tmaxC-tmaxM
;copy_VarMeta(tmaxM,Mxbias)
;copy_VarMeta(tmaxC,Cxbias)
;copy_VarMeta(tmaxC,Dxbias)

;Mnbias = tminM(0,:,:)-era_tmin
;Cnbias = tminC(0,:,:)-era_tmin
;Dnbias = tminC-tminM
;copy_VarMeta(tminM,Mnbias)
;copy_VarMeta(tminC,Cnbias)
;copy_VarMeta(tminC,Dnbias)

  dNames = getvardims(tscrnM1)
  nDims  = dimsizes(tscrnM1)
  nSize  = dimsizes(nDims)
  ;lon    = var_gd&$dNames(nSize-1)$
  ;lat    = var_gd&$dNames(nSize-2)$
if (isfilevar(moses_file,"tas")) then
Mtbias4 = tscrnM1(:,:,:)-era_t_rg
else
Mtbias4 = tscrnM1(:,0,:,:)-era_t_rg
end if
Ctbias4 = tscrnC1(:,0,:,:)-era_t_rg
if (isfilevar(moses_file,"tas")) then
Dtbias4 = tscrnC1(:,0,:,:)-tscrnM1
else
Dtbias4 = tscrnC1(:,0,:,:)-tscrnM1(:,0,:,:)
end if
Mtbias4!0 = "seasons"
Mtbias4!1 = "latitude"
Mtbias4!2 = "longitude_1"
if (isdim(tscrnM1,"time")) then
Mtbias4&seasons   = tscrnM1&time
else
Mtbias4&seasons   = tscrnM1&t;ime
end if
if (isdim(tscrnM1,"lat")) then
Mtbias4&latitude  = tscrnM1&lat
else
Mtbias4&latitude  = tscrnM1&latitude
end if
Mtbias4&longitude_1 = tscrnM1&$dNames(nSize-1)$ ;longitude_1
copy_VarMeta(Mtbias4,Ctbias4)
copy_VarMeta(tscrnC1(:,0,:,:),Dtbias4)

if (isfilevar(moses_file,"tas")) then
Mxbias4 = tmaxM1(:,:,:)-era_x_rg
else
Mxbias4 = tmaxM1(:,0,:,:)-era_x_rg
end if
Cxbias4 = tmaxC1(:,0,:,:)-era_x_rg
if (isfilevar(moses_file,"tas")) then
Dxbias4 = tmaxC1(:,0,:,:)-tmaxM1
else
Dxbias4 = tmaxC1(:,0,:,:)-tmaxM1(:,0,:,:)
end if
copy_VarMeta(Mtbias4,Mxbias4)
copy_VarMeta(Mtbias4,Cxbias4)
copy_VarMeta(tmaxC1(:,0,:,:),Dxbias4)

if (isfilevar(moses_file,"tas")) then
Mnbias4 = tminM1(:,:,:)-era_n_rg
else
Mnbias4 = tminM1(:,0,:,:)-era_n_rg
end if
Cnbias4 = tminC1(:,0,:,:)-era_n_rg
if (isfilevar(moses_file,"tas")) then
Dnbias4 = tminC1(:,0,:,:)-tminM1
else
Dnbias4 = tminC1(:,0,:,:)-tminM1(:,0,:,:)
end if
copy_VarMeta(Mtbias4,Mnbias4)
copy_VarMeta(Mtbias4,Cnbias4)
copy_VarMeta(tminC1(:,0,:,:),Dnbias4)

Mtbias = dim_avg_n_Wrap(Mtbias4,0)
Ctbias = dim_avg_n_Wrap(Ctbias4,0)
Dtbias = dim_avg_n_Wrap(Dtbias4,0)
Mxbias = dim_avg_n_Wrap(Mxbias4,0)
Cxbias = dim_avg_n_Wrap(Cxbias4,0)
Mnbias = dim_avg_n_Wrap(Mnbias4,0)
Cnbias = dim_avg_n_Wrap(Cnbias4,0)
Dxbias = dim_avg_n_Wrap(Dxbias4,0)
Dnbias = dim_avg_n_Wrap(Dnbias4,0)

; Plot ----------------------------------------------------------------------

print("Plotting 3 Panel Bias Plots")
print("")

;ann------------------------------

Panel3Plot("bias_ann_tscrn_3p",Ctbias,Mtbias,Dtbias(:,:),jpeg)
;Panel3Plot("bias_ann_tmax_3p",Cxbias,Mxbias,Dxbias(:,:),jpeg)
Panel3Plot("bias_ann_tmin_3p",Cnbias,Mnbias,Dnbias(:,:),jpeg)

;djf------------------------------

Panel3Plot("bias_djf_tscrn_3p",Ctbias4(0,:,:),Mtbias4(0,:,:),Dtbias4(0,:,:),jpeg)
;Panel3Plot("bias_djf_tmax_3p",Cxbias4(0,:,:),Mxbias4(0,:,:),Dxbias4(0,:,:),jpeg)
Panel3Plot("bias_djf_tmin_3p",Cnbias4(0,:,:),Mnbias4(0,:,:),Dnbias4(0,:,:),jpeg)

;jja------------------------------

Panel3Plot("bias_jja_tscrn_3p",Ctbias4(2,:,:),Mtbias4(2,:,:),Dtbias4(2,:,:),jpeg)
;Panel3Plot("bias_jja_tmax_3p",Cxbias4(2,:,:),Mxbias4(2,:,:),Dxbias4(2,:,:),jpeg)
Panel3Plot("bias_jja_tmin_3p",Cnbias4(2,:,:),Mnbias4(2,:,:),Dnbias4(2,:,:),jpeg)

end

