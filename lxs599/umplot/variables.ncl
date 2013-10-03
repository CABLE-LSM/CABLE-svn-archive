; Main Code ========================================================
; created by Lauren Stevens: 18-21 April 2011
; during the NCL course at BMTC. 
; ==================================================================

load "/home/cmar/ste69f/umplot/plot.ncl"
load "/home/cmar/ste69f/umplot/plot_mmn.ncl"
load "/home/cmar/ste69f/umplot/plot_mp.ncl"
load "/home/cmar/ste69f/umplot/plot_ann.ncl"
load "/home/cmar/ste69f/umplot/plot_ind.ncl"
load "/home/cmar/ste69f/umplot/regrid.ncl"

procedure run_ncl(fileDir:string,Runid:string,numYr:string,Tiles:string,\
              sLevs:string,Res:string,Split:logical,Taylor:logical,    \
              compDir:string,lcable:logical,jpeg:logical)

begin

; Reading Data =====================================================

; Set-up -----------------------------------------------------------  

  PLOT    = True
  
  if (PLOT) then
   pltDir = "./"     ; "./plots"    ; fileDir+"/plots"
   type   = "ps" 
  end if

; Switches etc -----------------------------------------------------

  xaank_t = "n" ;systemfunc("echo $T15")
  if (xaank_t .eq."y") then
   t15 = True
  else
   t15 = False
  end if

  if (Tiles.eq."17") then
   ntiles = 17
   npft   = 13
  else
   ntiles = 9
   npft   = 5
  end if

  if (sLevs.eq."6") then
   nsoil  = 6
  else
   nsoil  = 4
  end if

  if (Res.eq."n96")
   nres   = 96
   nx     = 192
   ny     = 145
  else
   nres   = 48
   nx     = 96
   ny     = 73
  end if

; Directories ------------------------------------------------------

  ObsDir = "/home/cmar/ste69f/umplot/obs/"
  ;ModDir = fileDir
  ModDir = systemfunc("echo $PWD")
  ;ModDir = systemfunc("echo $DIRW")
  tname =  systemfunc("echo $tname")

; Obs & Model netCDF ----------------------------------------------- 

  print("Dmget Files")
  system("dmget " + ModDir + "/seasonal_means_"  +numYr+"yrs.nc" \
            + " " + ModDir + "/monthly_means_"   +numYr+"yrs.nc" \
            + " " + ModDir + "/Tseasonal_means_" +numYr+"yrs.nc" \
            + " " + ModDir + "/Tmonthly_means_"  +numYr+"yrs.nc" \
  ;         + " " + ModDir + "/Rseasonal_means_" +numYr+"yrs.nc" \
  ;         + " " + ModDir + "/Rmonthly_means_"  +numYr+"yrs.nc" \
            + " " + ObsDir + "ERAi_precip.nc" \
            + " " + ObsDir + "ERAi_tp2m.nc" \
            + " " + ObsDir + "ERAi_evap_cdo.nc" \
            + " " + ObsDir + "ERAi_tmax_tmin.nc" \
            + " " + ObsDir + "bom_precip_monthly.nc" )
  ;         + " " + ObsDir + "Evap_Sea-Ann_Means_89-09_unpk.nc" ) ;\
  ;         + " " + ObsDir + "ERAi_Means_tmax_tmin_89-09_unpk.nc" )
  ;print("Dmget Files Done")

  print("Reading in NetCDFs")
  ;sea_file   = addfile(ModDir+"/season_"      +numYr+"yrs.nc","r")
  ;tmp_file   = addfile(ModDir+"/Temp_monthly_"+numYr+"yrs.nc","r")
  ;mon_file   = addfile(ModDir+"/mnmean_"      +numYr+"yrs.nc","r")

  sea_file   = addfile(ModDir+"/seasonal_means_" +numYr+"yrs.nc","r")
  mon_file   = addfile(ModDir+"/monthly_means_"  +numYr+"yrs.nc","r")
if(isfilepresent(ModDir+"/Tseasonal_means_"+numYr+"yrs.nc")) then
  tsea_file  = addfile(ModDir+"/Tseasonal_means_"+numYr+"yrs.nc","r")
  tmon_file  = addfile(ModDir+"/Tmonthly_means_" +numYr+"yrs.nc","r")
end if
if(isfilepresent(ModDir+"/Rseasonal_means_"+numYr+"yrs.nc")) then
  rsea_file  = addfile(ModDir+"/Rseasonal_means_"+numYr+"yrs.nc","r")
  rmon_file  = addfile(ModDir+"/Rmonthly_means_" +numYr+"yrs.nc","r")
end if
  obs_prec   = addfile(ObsDir+"ERAi_precip.nc","r")
  obs_tp2m   = addfile(ObsDir+"ERAi_tp2m.nc","r")
  obs_evp    = addfile(ObsDir+"ERAi_evap_cdo.nc","r")
  obs_tmxx   = addfile(ObsDir+"ERAi_tmax_tmin.nc","r")
  bomopr     = addfile(ObsDir+"bom_precip_monthly.nc","r")
  ;obs_evp   = addfile(ObsDir+"Evap_Sea-Ann_Means_89-09_unpk.nc","r")
  ;obs_tmxx  = addfile(ObsDir+"ERAi_Means_tmax_tmin_89-09_unpk.nc","r")
  ;print("Reading in NetCDFs Done")
  ;printVarSummary()

; Variables ========================================================

;sea_file ----------------------------------------------------------
if (isfilevar(sea_file,"tscrn")) then
  tscrn    = sea_file    ->tscrn       ;Screen Temperature degK (T at 1.5m)
else
  tscrn    = sea_file    ->$tname$     ;Screen Temperature degK (T at 1.5m)
end if
  precip   = sea_file    ->precip      ;field77     ; Total Precip Accum kg/m2
if(isfilevar(sea_file,"temp")) then
  tsurf    = sea_file    ->temp        ; Surface Temperature (deg K)
end if
  sh       = sea_file    ->sh          ; Sensible Heat
  lh       = sea_file    ->lh          ; Latent Heat
  solar    = sea_file    ->solar       ; Net short wave radiation flux
  longwave = sea_file    ->longwave    ; Net long wave radiation flux
  lw       = sea_file    ->ilr         ;
  sw       = sea_file    ->field203    ; Total downward surface solar fl
  csky     = sea_file    ->field208    ; Clear-sky flux (type II) solar
if(isfilevar(sea_file,"field1526")) then
  evap1    = sea_file    ->field1526   ; Evaporation from soil or canopy
  evap2    = sea_file    ->field1527   ; Evaporation from soil or canopy
end if
if(isfilevar(sea_file,"field1894")) then
  lai      = sea_file    ->field1894   ; Leaf Area Index (LAI)
  canht    = sea_file    ->field1895   ; Canopy height
  sfrac    = sea_file    ->field1893   ; Surface Fraction
end if
if(isfilevar(sea_file,"field322")) then
  lalb     = sea_file    ->field322    ; Mean Land Albedo
end if
if(isfilevar(sea_file,"snowdepth")) then
  sndep    = sea_file    ->snowdepth   ; Snow Depth
end if
if(isfilevar(sea_file,"temp_1")) then
  lsftp    = sea_file    ->temp_1      ; Land Surace Temperature (deg K)
end if
if(isfilevar(sea_file,"field202_1")) then
  nsradf   = sea_file    ->field202_1  ; Net Sruface Rad Flux (no lh or sh)
end if
if(isfilevar(sea_file,"field141")) then 
  snmelt   = sea_file    ->field141    ;
end if
if(isfilevar(sea_file,"field1532")) then
  srnf     = sea_file    ->field1532   ; Surf Runoff                     ; **** rsea_file
  drnf     = sea_file    ->field1533   ; Sub-Surf Runoff (Drainage)      ; **** rsea_file
end if
if(isfilepresent(ModDir+"/Rseasonal_means_"+numYr+"yrs.nc")) then
if(isfilevar(rsea_file,"field1532")) then
  srnf     = rsea_file    ->field1532   ; Surf Runoff                     ; **** sea_file
  drnf     = rsea_file    ->field1533   ; Sub-Surf Runoff (Drainage)      ; **** sea_file
end if
end if
if(isfilevar(sea_file,"field1386")) then
  frozf    = sea_file    ->field1386   ; Frozen SM Frac Layer            ; **** rsea_file
end if
if(isfilepresent(ModDir+"/Rseasonal_means_"+numYr+"yrs.nc")) then
if(isfilevar(rsea_file,"field1386")) then
  frozf    = rsea_file    ->field1386   ; Frozen SM Frac Layer            ; **** rsea_file
end if
end if
if(isfilevar(sea_file,"sm_1").and.isfilevar(sea_file,"sm")) then
  smcl     = sea_file    ->sm_1        ; Soil Moisture in Layer
  ;smc      = sea_file    ->sm          ; Soil Moisture
end if
if(.not.isfilevar(sea_file,"sm_1").and.isfilevar(sea_file,"sm")) then
  smcl     = sea_file    ->sm          ; Soil Moisture in Layer
end if
if(isfilevar(sea_file,"field1388")) then
  gpp      = sea_file    ->field1388   ; Gross Primary Production
  npp      = sea_file    ->field1389   ; Net Primary Production
  frp      = sea_file    ->field1390   ; Plant Respiration
  frs      = sea_file    ->field1523   ; Soil Respiration
end if
if (lcable.and.t15) then
  latT     = tscrn&latitude 
  lonT     = tscrn&longitude_1
  soilT    = new((/4,nsoil,ntiles,dimsizes(latT),dimsizes(lonT)/),float)
  soilT    = -9999
  soilT@_FillValue = -9999
  ;soilT@long_name  = "Soil Temperatures in Layers from CABLE"
  ;soilT@units      = "deg K"
  ;soilT!0          = "t"
  ;soilT!1          = "soillev"
  ;soilT!2          = "pseudo"
  ;soilT!3          = "latitude"
  ;soilT!4          = "longitude_1"
  ;soilT&t          = tscrn&t
  ;soilT&soillev    = ispan(1,nsoil,1)
  ;soilT&pseudo     = ispan(1,ntiles,1)
  ;soilT&longitude_1= tscrn&longitude_1
  ;soilT&latitude   = tscrn&latitude
  soilT(:,0,:,:,:) = sea_file ->temp_2 ;soilT1
  soilT(:,1,:,:,:) = sea_file ->temp_3 ;soilT2
  soilT(:,2,:,:,:) = sea_file ->temp_4 ;soilT3
  soilT(:,3,:,:,:) = sea_file ->temp_5 ;soilT4
  soilT!1          = "soillev"
  soilT&soillev    = ispan(1,nsoil,1)
 if (nsoil.eq.6) then
  soilT(:,4,:,:,:) = sea_file ->temp_6 ;soilT5
  soilT(:,5,:,:,:) = sea_file ->temp_7 ;soilT6
 end if
end if
;if (Runid.eq."xaank".or.Runid.eq."xaanv".and.lcable) then
if (isfilevar(sea_file,"soiltemp")) then
  deepST   = sea_file    ->soiltemp    ; 
end if
;tsea_file----------------------------------------------------------
if(isfilepresent(ModDir+"/Tseasonal_means_"+numYr+"yrs.nc")) then
  tmax     = tsea_file   ->temp        ; Maximum Screen Temp degK
  tmin     = tsea_file   ->temp_1      ; Minimum Screen temp degK
end if
;mon_file-----------------------------------------------------------
  prmn     = mon_file    ->precip      ;field77
if (isfilevar(mon_file,"tscrn")) then
  t2mn     = mon_file    ->tscrn       ;temp_3
else
  t2mn     = mon_file    ->$tname$     ;temp_3
end if
if(isfilevar(sea_file,"field1894")) then
  laimn    = mon_file    ->field1894   ;
  chtmn    = mon_file    ->field1895   ;
end if
if(isfilevar(sea_file,"temp_1")) then
  lstmn    = mon_file    ->temp_1      ;
end if
if(isfilevar(sea_file,"field322")) then
  albmn    = mon_file    ->field322    ;
end if
if(isfilevar(sea_file,"snowdepth")) then
  sndmn    = mon_file    ->snowdepth   ;
end if
if(isfilevar(mon_file,"field1532")) then
  srmn     = mon_file    ->field1532   ; **** rmon_file
  drmn     = mon_file    ->field1533   ; **** rmon_file
end if
if(isfilepresent(ModDir+"/Rseasonal_means_"+numYr+"yrs.nc")) then
if(isfilevar(rmon_file,"field1532")) then
  srmn     = rmon_file    ->field1532   ; **** mon_file
  drmn     = rmon_file    ->field1533   ; **** mon_file
end if
end if
;tmp_file-----------------------------------------------------------
if(isfilepresent(ModDir+"/Tseasonal_means_"+numYr+"yrs.nc")) then
  txmn     = tmon_file    ->temp        ;
  tnmn     = tmon_file    ->temp_1      ;
end if
;obs----------------------------------------------------------------
  era_prec = obs_prec    ->a            ; Era Int Prec  Obs
  era_tp2m = obs_tp2m    ->a            ; Era Int Tscrn Obs
  era_evap = obs_evp     ->e            ; Era Int Evap  Obs
  era_tmax = obs_tmxx    ->mx2t         ; Era Int Tmax  Obs
  era_tmin = obs_tmxx    ->mn2t         ; Era Int Tmin  Obs
  bom_pr   = bomopr      ->precip       ; Bureau Precip Obs
;MetaData-----------------------------------------------------------

  precip              = precip*24*60*60  ;*48
  precip@long_name    = "Precipitation"
  precip@units        = "mm/day"
  prmn                = prmn*24*60*60    ;*48
  prmn@long_name      = "Precipitation"
  prmn@units          = "mm/day"

  tscrn@long_name     = "Temperature at 1.5m"
  tscrn@units         = "deg K"
if (isvar("tmax")) then
  tmax@long_name      = "Maximum Temperature at 1.5m"
  tmax@units          = "deg K"
  tmax!3              = "longitude_1"
  tmin@long_name      = "Minimum Temperature at 1.5m"
  tmin@units          = "deg K"
  tmin!3              = "longitude_1"

  txmn!3              = "longitude_1"
  tnmn!3              = "longitude_1"

  tamp                = tmax - tmin
  copy_VarMeta(tmax,tamp)
  tamp@long_name      = "Amplitude of Screen Temperature at 1.5m"
  tamp@units          = "deg K"
  tamn                = txmn - tnmn
  copy_VarMeta(txmn,tamn)
  tamn@long_name      = "Amplitude of Screen Temperature at 1.5m"
  tamn@units          = "deg K"
end if

if(isvar("evap1")) then
  evap                = ( evap1 + evap2 ) *365 *24 *60 *60
  copy_VarMeta(evap1,evap)
  evap@units          = "mm/yr"
  evap@long_name      = "Evaporation - soil+canopy"
end if

if(isvar("lai")) then  
  canht@long_name = "Canopy Height"
  canht@units     = "m"
  do i= 0,npft-1,1
   lai(:,i,:,:) = where(sfrac(:,i,:,:).gt.0,lai(:,i,:,:),0)
  end do
  lai@long_name   = "Leaf Area Index"
  lai@units       = "LAI"
end if
if(isvar("laimn")) then
  chtmn@long_name = "Canopy Height"
  chtmn@units     = "m"
  do i= 0,npft-1,1
  do j= 0,11,1
   laimn(j,i,:,:) = where(sfrac(0,i,:,:).gt.0,laimn(j,i,:,:),0)
  end do
  end do
  laimn@long_name = "Leaf Area Index"
  laimn@units     = "LAI"
end if

  csky  = where(csky.ne.0,csky,-9999) ;csky@_FillValue)
  sw    = where(sw.le.0,sw@_FillValue,sw)     ;sw@_FillValue)
  ;sw    = where(sw.ne.0,sw,-9999)     ;sw@_FillValue)
  ;lw    = where(lw.ne.0,lw,-9999)
  sw@long_name = "Total Downward Surface SW Flux"
  lw@long_name = "Downward LW Rad Flux: Surface"
  csky@long_name = "Clear-Sky (II) Down Surface SW Flux"
  calb  = (1-(sw/csky))
  copy_VarMeta(csky,calb)
  calb@long_name = "Cloud Albedo as Calculated by Jaeger"
  surfalb  = (-1)*((solar-sw)/sw)
  copy_VarMeta(sw,surfalb)
  surfalb@long_name = "Surface Albedo as Calculated by method in UMui"
  swnet = csky*(1-(-1)*(solar-sw)/sw)*(1-(1-(sw/csky)))    ; SWnet = SWclear(1-Acloud)(1-Asurface) = solar
  ;swnet = csky*(1-calb)*(1-surfalb)
  copy_VarMeta(sw,swnet)
  swnet@long_name = "Net SW as Calculated by ..."
  ;csky = where(ismissing(csky),0,csky)
  ;sw    = where(ismissing(sw),0,sw)

if(isvar("snmelt")) then
  seb           = solar+longwave-sh-lh-snmelt
else
  seb           = solar+longwave-sh-lh
end if
  copy_VarMeta(solar,seb)
  seb@long_name = "Surface Energy Balance"
  seb@units     = "W/m2"

if(isvar("srnf")) then
  srnf                = srnf *24 *60 *60
  srnf@units          = "mm/day"
  drnf                = drnf *24 *60 *60
  drnf@units          = "mm/day"
  srmn                = srmn *24 *60 *60
  srmn@units          = "mm/day"
  drmn                = drmn *24 *60 *60
  drmn@units          = "mm/day"

  tot_rnf            = srnf + drnf
  copy_VarMeta(srnf,tot_rnf)
  tot_rnf@long_name  = "Total Runoff: Surface+Drainage"
  tot_rnf@units      = "mm/day"

if(isvar("evap")) then
  wbal                = precip - evap/365 - tot_rnf
  copy_VarMeta(precip,wbal)
  wbal@long_name      = "Water Balance"
  wbal@units          = "mm/day"
end if
end if

if(isvar("smcl")) then
  if (nsoil.eq.6) then
   slayers            = (/0.022,0.058,0.154,0.409,1.085,2.872/) ; = 4.6 m
  else
   if (nsoil.eq.4 .and. lcable) then
   ;slayers           = (/0.022,0.158,0.900,3.120/)             ; = 4.2 m
    slayers           = (/0.022,0.158,0.900,1.920/)             ; = 3 m
   else  
    ; moses 4 soil layers
    slayers           = (/0.100,0.250,0.650,2.000/)             ; = 3 m
   end if
  end if

if (Runid.eq."saaba") then
nsoil = 4
end if
  do i=0,nsoil-1,1
   smcl(:,i,:,:)       = smcl(:,i,:,:) / (slayers(i) *1000)
  end do
  smcl@units           = "wb"
  if (nsoil.eq.6) then
   smcl@long_name = "Soil moisture in each (of 6) levels"
  end if

  smc                  = dim_sum_n_Wrap(smcl,1)
  smc@units            = "wbtot"
  smc@long_name = "Total Soil Moisture"
; frozf ???
if (isvar("frozf")) then
  frozf@long_name = "Frozen soil moisture fraction"
end if

end if

  if (lcable.and.t15) then
   soilTemp = new((/4,nsoil,dimsizes(latT),dimsizes(lonT)/),float)
   do i = 0,nsoil-1,1
    do j = 0,3,1
     soilTemp(j,i,:,:) = dim_sum_n_Wrap(sfrac(j,:,:,:) * soilT(j,i,:,:,:),0) 
    end do
   end do
   soilTemp@long_name = "Soil Temperature"
   soilTemp@units     = "dek K"
   soilTemp!0 = "t"
   soilTemp!1 = "soillev"
   soilTemp!2 = "latitude"
   soilTemp!3 = "longitude_1"
   soilTemp&t           = soilT&t
   soilTemp&soillev     = soilT&soillev
   soilTemp&latitude    = soilT&latitude
   soilTemp&longitude_1 = soilT&longitude_1
  end if

; ERA Interim ------------------------------------------------------

  era_prec@units      = "mm/day"
  era_prec@long_name  = "ERA Interim Precipitation"
  era_tp2m            = era_tp2m+273.15
  era_tp2m@long_name  = "ERA Interim Temperature at 2m"
  era_tp2m@units      = "deg K"

  era_evp           = new((/4,121,240/),double)
  era_evp(0,:,:)    = era_evap(0,:,:)+era_evap(1,:,:)+era_evap(11,:,:)
  era_evp(1,:,:)    = era_evap(2,:,:)+era_evap(3,:,:)+era_evap(4,:,:)
  era_evp(2,:,:)    = era_evap(5,:,:)+era_evap(6,:,:)+era_evap(7,:,:)
  era_evp(3,:,:)    = era_evap(8,:,:)+era_evap(9,:,:)+era_evap(10,:,:)
  era_evp           = era_evp * (-1) * 2000 * 365 /3 ; 200
  copy_VarMeta(era_prec,era_evp)
  era_evp@long_name = "ERA Interim Evaporation"
  era_evp@units     = "mm/yr"
  ;era_evp!2         = "longitude"

  ;era_tmax           = era_tmax+273.16
  ;era_tmax@units     = "deg K"
  eraT_max           = new((/4,121,240/),double)
  eraT_max(0,:,:)    = (era_tmax(0,:,:)+era_tmax(1,:,:)+era_tmax(11,:,:))/3
  eraT_max(1,:,:)    = (era_tmax(2,:,:)+era_tmax(3,:,:)+era_tmax(4,:,:))/3
  eraT_max(2,:,:)    = (era_tmax(5,:,:)+era_tmax(6,:,:)+era_tmax(7,:,:))/3
  eraT_max(3,:,:)    = (era_tmax(8,:,:)+era_tmax(9,:,:)+era_tmax(10,:,:))/3
  ;copy_VarMeta(era_tp2m,eraT_max)
  eraT_max@long_name = "ERA Interim Maximum Temperature at 2m"
  eraT_max@units     = "deg K"
  eraT_max!1         = "latitude"
  eraT_max!2         = "longitude"
  eraT_max&latitude  = era_tmax&latitude
  eraT_max&longitude = era_tmax&longitude

  ;era_tmin = era_tmin+273.16
  ;era_tmin@units = "deg K"
  eraT_min = new((/4,121,240/),double)
  eraT_min(0,:,:) = (era_tmin(0,:,:)+era_tmin(1,:,:)+era_tmin(11,:,:))/3
  eraT_min(1,:,:) = (era_tmin(2,:,:)+era_tmin(3,:,:)+era_tmin(4,:,:))/3
  eraT_min(2,:,:) = (era_tmin(5,:,:)+era_tmin(6,:,:)+era_tmin(7,:,:))/3
  eraT_min(3,:,:) = (era_tmin(8,:,:)+era_tmin(9,:,:)+era_tmin(10,:,:))/3
  copy_VarMeta(era_tp2m,eraT_min)
  eraT_min!2 = "longitude"
  eraT_min@long_name = "ERA Interim Minimum Temperature at 2m"

  era_tamp            = eraT_max-eraT_min
  copy_VarMeta(eraT_max,era_tamp)
  era_tamp@long_name  = "ERAi Amplitude of Screen Temperature"
  era_tamp@units      = "deg K"
 
  ;sndep = sndep */  ; units: water eq kgm-2 (mm) -> ??? 

; Bureau -----------------------------------------------------------

  bompr          = dim_avg_n_Wrap(bom_pr,0) 
  bmpr           = new((/4,139,178/),float)
  bmpr(0,:,:)    = (bompr(0,:,:)+bompr(1,:,:)+bompr(11,:,:))/3
  bmpr(1,:,:)    = (bompr(2,:,:)+bompr(3,:,:)+bompr(4,:,:))/3
  bmpr(2,:,:)    = (bompr(5,:,:)+bompr(6,:,:)+bompr(7,:,:))/3
  bmpr(3,:,:)    = (bompr(8,:,:)+bompr(9,:,:)+bompr(10,:,:))/3
  bmpr!1         = "latitude"
  bmpr!2         = "longitude"
  bmpr&latitude  = bompr&lat
  bmpr&longitude = bompr&lon
  bmpr@long_name = "Bureau Australian Precipitation"
  bmpr@units     = "mm/day"

; regrid obs -------------------------------------------------------

  print("Regridding Observations")
  era_t_rg              = regrid_bilin(era_tp2m,tscrn)
  era_t_rg@long_name    = "ERAi Tscrn regridded to UM N96"
  era_t_rg@units        = "deg K"

  era_p_rg              = regrid_bilin(era_prec,precip)
  era_p_rg@long_name    = "ERAi Precip regridded to UM N96"
  era_p_rg@units        = "mm/day"

if(isvar("tmax")) then
  era_x_rg              = regrid_bilin(eraT_max,tmax)
  era_x_rg@long_name    = "ERAi Tmax regridded to UM N96"
  era_x_rg@units        = "deg K"

  era_n_rg              = regrid_bilin(eraT_min,tmin)
  era_n_rg@long_name    = "ERAi Tmin regridded to UM N96"
  era_n_rg@units        = "deg K"

  era_d_rg              = regrid_bilin(era_tamp,tamp)
  era_d_rg@long_name    = "ERAi Tdiff regridded to UM N96"
  era_d_rg@units        = "deg K"
end if

if(isvar("evap")) then
  era_e_rg              = regrid_bilin(era_evp,evap)
  era_e_rg@long_name    = "ERAi Evap regridded to UM N96"
  era_e_rg@units        = "mm/day"
end if

  bom_p_rg              = regrid_bilin(bmpr,precip)
  bom_p_rg@long_name    = "Bureau of Meteorlogy Precip"
  bom_p_rg@units        = "mm/day"
  ;print("Regridding Observations Done")

; biases -----------------------------------------------------------

  print("Calculating Biases")
  bias_tscrn            = tscrn(:,0,:,:)-era_t_rg(:,:,:)
  copy_VarMeta(era_t_rg(:,:,:),bias_tscrn)
  bias_tscrn@long_name  = "Screen Temperature Bias: UM-ERAi"
  bias_tscrn@units      = "deg K"
  bias_tscrn!0          = "t"
  bias_tscrn&t          = tscrn&t

  bias_precip           = precip(:,0,:,:)-era_p_rg(:,:,:)
  copy_VarMeta(era_p_rg(:,:,:),bias_precip)
  bias_precip@long_name = "Precipitation Bias: UM-ERAi"
  bias_precip@units     = "mm/day"
  bias_precip!0         = "t"
  bias_precip&t         = precip&t

if(isvar("tmax")) then
  bias_tmax             = tmax(:,0,:,:)-era_x_rg
  copy_VarMeta(tmax(:,0,:,:),bias_tmax)
  bias_tmax@long_name   = "Maximum Temperature Bias: UM-ERAi"

  bias_tmin             = tmin(:,0,:,:)-era_n_rg
  copy_VarMeta(tmin(:,0,:,:),bias_tmin)
  bias_tmin@long_name   = "Minimum Temperature Bias: UM-ERAi"

  bias_tamp             = tamp(:,0,:,:)-era_d_rg
  copy_VarMeta(tamp(:,0,:,:),bias_tamp)
  bias_tamp@long_name   = "Amplitude of Temperature Bias: UM-ERAi"
end if

if(isvar("evap")) then
  bias_evap             = evap(:,0,:,:)-era_e_rg
  copy_VarMeta(evap(:,0,:,:),bias_evap)
  bias_evap@long_name   = "Evaporation Bias: UM-ERAi"
  ;print("Calculating Biases Done")
end if

  bias_bmpr             = precip(:,0,:,:) - bom_p_rg
  copy_VarMeta(bom_p_rg(:,:,:),bias_bmpr)
  bias_bmpr@long_name   = "Precipitation Bias: UM-BoM"
  bias_bmpr@units       = "mm/day"

; Plot Seasons======================================================

;print("Start Plotting") 
print("Plotting Seasons")

; Levels for contours
lev    =  "0"   ; For AutomaticLevels ?
levtp  = (/240,250,260,270,275,280,285,290,295,300,305/)
levpr  = (/0.5,1,2,4,8,16/)
levbpr = (/0.5,1,2,3,4,5,6,7,8/)
levrnf = (/0.1,0.5,1,2,4,8,16/)
cmapr  = (/"white","black","cornsilk3","dodgerblue","chartreuse3","gold","darkorange","red","firebrick"/)
cmaprn = (/"white","black","cornsilk3","blueviolet","dodgerblue","chartreuse3","gold","darkorange","red","firebrick"/) ; mediumpurple
levdf  = (/2,4,6,8,10,12,14,16,18,20,22,24,26/)
levam  = ispan(2,16,1)
levlai = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8/)
levalb = (/0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14/)
;levsnw = (/5,10,15,20,25,30,35,40,45,50,55,60,80,100,120,140,160,180/)
levsnw = (/5,10,20,40,60,80,100,150,200,400,600,800,1000,1500/)
cmapsn = (/"white","black","cornsilk3","lightslateblue","lightblue","dodgerblue","blue", \
           "olivedrab1","chartreuse3","forestgreen","khaki1","gold","darkgoldenrod2", \
           "darkorange","tomato1","red","firebrick"/)
levevp = (/-250,0,250,500,750,1000,1250,1500,1750,2000/)
levrad = (/-40,-20,0,20,40,60,80,100,120,140,160,180,200,220,240/)
levsal = (/.02,.04,.06,.08,0.1,.12,.14,.16,.18,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
levsfr = (/0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1/)
levsh  = (/-60,-40,-20,0,20,40,60,80,100,120,140/)
levlh  = (/-20,0,20,40,60,80,100,120,140,160,180/)
levwb  = (/-4,-2,-1,-0.5,-0.1,-0.05,-0.01,0,0.01,0.05,0.1,0.5,1,2,4/)
levht  = (/0,0.55,0.57,0.6,15.5,17,20,35/)
levsm  = (/0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5/)
cmapsm = "percent_11lev"
;levsc  = (/-200,-150,-100,-50,0,50,100,150,200,250,300/)
levsc  = (/0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/)
levsw  = (/0,50,100,150,200,250,300,350/)
levseb = (/-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6/)
levlw1 = (/50,100,150,200,250,300,350,400/)
levlw2 = (/-220,-200,-180,-160,-140,-120,-100,-80,-60,-40,-20,0,20/)
levsw1 = (/0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300/)
;levsw2 = (//)
levmlt = (/0.5,1,2,4,6,8,12,16,24/)

;BlAqGrYeOrReVi200
;perc2_9lev
Panel_Plot(pltDir,"seas_tscrn" ,numYr,type,"BlAqGrYeOrReVi200",tscrn(:,0,:,:) ,False,True,levtp,False,jpeg)
AnnPlot(pltDir,"ann_tscrn"     ,numYr,type,"BlAqGrYeOrReVi200",tscrn(:,0,:,:) ,False,True,levtp,False,jpeg)
Panel_Plot(pltDir,"seas_precip",numYr,type,cmapr,precip(:,0,:,:),False,True,levpr,False,jpeg)
AnnPlot(pltDir,"ann_precip"    ,numYr,type,cmapr,precip(:,0,:,:),False,True,levpr,False,jpeg)
if(isvar("tmax")) then
Panel_Plot(pltDir,"seas_tmax"  ,numYr,type,"BlAqGrYeOrReVi200",tmax(:,0,:,:)  ,False,True,levtp,True,jpeg) 
Panel_Plot(pltDir,"seas_tmin"  ,numYr,type,"BlAqGrYeOrReVi200",tmin(:,0,:,:)  ,False,True,levtp,True,jpeg)
Panel_Plot(pltDir,"seas_tamp"  ,numYr,type,"BlAqGrYeOrReVi200",tamp(:,0,:,:)  ,False,True,levam,True,jpeg)
end if

;cmapr2  = (/"white","black","cornsilk3","skyblue","dodgerblue","darkolivegreen1","chartreuse3","gold","darkorange","red","firebrick","darkorchid4"/)
cmapr2  = (/"white","black","cornsilk3","skyblue","dodgerblue","darkolivegreen2","chartreuse3","lightgoldenrod1","gold","darkorange","red","firebrick"/)
;cmapr2  = (/"white","black","cornsilk3","skyblue","dodgerblue","darkolivegreen2","chartreuse3","lightgoldenrod1","gold","orange","darkorange2","red","firebrick"/)
Panel_Plot(pltDir,"erai_tscrn" ,"20",type,"BlAqGrYeOrReVi200",era_tp2m ,False,True,levtp ,False,jpeg)
Panel_Plot(pltDir,"erai_precip","20",type, cmapr2            ,era_prec ,False,True,levbpr ,False,jpeg)
;Panel_Plot(pltDir,"erai_tmax"  ,"20",type,"BlAqGrYeOrReVi200",eraT_max ,False,True,levtp ,True,jpeg)
;Panel_Plot(pltDir,"erai_tmin"  ,"20",type,"BlAqGrYeOrReVi200",eraT_min ,False,True,levtp ,True,jpeg)
;Panel_Plot(pltDir,"erai_tamp"  ,"20",type,"BlAqGrYeOrReVi200",era_tamp ,False,True,levam ,True,jpeg)
;Panel_Plot(pltDir,"erai_evap"  ,"20",type, cmapsn            ,era_evp  ,False,True,levevp,True,jpeg)
AnnPlot(pltDir,"erai_atscrn","20",type,"BlAqGrYeOrReVi200",era_tp2m,False,True,levtp,False,jpeg)
AnnPlot(pltDir,"erai_aprec" ,"20",type,cmapr2             ,era_prec,False,True,levbpr,False,jpeg)

if(isvar("lai")) then
do n = 0,npft-1
b = n+1
lai@long_name   = "Leaf Area Index for PFT "+b
Panel_Plot(pltDir,"seas_lai"+b,numYr,type,"BlAqGrYeOrReVi200",lai(:,n,:,:),False,True,levlai,True,False)
end do
;if (.not.lcable) then
do n = 0,npft-1
c = n + 1
canht@long_name = "Canopy Height for PFT "+c
Panel_Plot(pltDir,"seas_canht"+c,numYr,type,"BlAqGrYeOrReVi200",canht(:,n,:,:),False,True,levdf,True,False)
end do
end if

if(isvar("lalb")) then
Panel_Plot(pltDir,"seas_lalb"  ,numYr,type,"BlAqGrYeOrReVi200",lalb(:,0,:,:)  ,False,True,levalb,True,jpeg)
end if
if(isvar("lsftp")) then
Panel_Plot(pltDir,"seas_lsftp" ,numYr,type,"BlAqGrYeOrReVi200",lsftp(:,0,:,:) ,False,True,levtp,True,jpeg)
end if
if(isvar("sfrac")) then
do n = 0,ntiles-1
d = n + 1
Panel_Plot(pltDir,"seas_sfrac"+d,numYr,type,"BlAqGrYeOrReVi200",sfrac(:,n,:,:),False,True,levsfr,True,False)
end do
end if

if(isvar("sndep")) then
Panel_Plot(pltDir,"seas_sndep" ,numYr,type,cmapsn,sndep(:,0,:,:) ,False,True,levsnw,True,jpeg)
end if
if(isvar("nsradf")) then
Panel_Plot(pltDir,"seas_nsradf",numYr,type,"BlAqGrYeOrReVi200",nsradf(:,0,:,:),False,True,levrad,True,jpeg)
end if
if(isvar("evap")) then
Panel_Plot(pltDir,"seas_evap"  ,numYr,type,cmapsn,evap(:,0,:,:)  ,False,True,levevp,True,jpeg)
end if
if(isvar("tsurf")) then
Panel_Plot(pltDir,"seas_tsurf" ,numYr,type,"BlAqGrYeOrReVi200",tsurf(:,0,:,:) ,False,True,levtp,True,jpeg)
end if
Panel_Plot(pltDir,"seas_solar",numYr,type,"BlAqGrYeOrReVi200",solar(:,0,:,:),False,True,levsw1,False,jpeg)
Panel_Plot(pltDir,"seas_swdown",numYr,type,"BlAqGrYeOrReVi200",sw(:,0,:,:),False,True,levsw1,False,jpeg)
Panel_Plot(pltDir,"seas_longwave",numYr,type,"BlAqGrYeOrReVi200",longwave(:,0,:,:),False,True,levlw2,False,jpeg)
Panel_Plot(pltDir,"seas_lwdown",numYr,type,"BlAqGrYeOrReVi200",lw(:,0,:,:),False,True,levlw1,False,jpeg)

Panel_Plot(pltDir,"seas_sh"      ,numYr,type,"BlAqGrYeOrReVi200",sh(:,0,:,:) ,False,True,levsh,True,jpeg)
Panel_Plot(pltDir,"seas_lh"      ,numYr,type,"BlAqGrYeOrReVi200",lh(:,0,:,:) ,False,True,levlh,True,jpeg)
if(isvar("snmelt")) then
Panel_Plot(pltDir,"seas_snmlt" ,numYr,type,"BlAqGrYeOrReVi200",snmelt(:,0,:,:) ,False,True,levmlt,True,jpeg)
end if
if(isvar("srnf")) then
Panel_Plot(pltDir,"seas_tot_rnf",numYr,type,cmaprn,tot_rnf(:,0,:,:),False,True,levrnf,True,jpeg)
if(isvar("wbal"))then
Panel_Plot(pltDir,"seas_wbal"    ,numYr,type,"BlAqGrYeOrReVi200",wbal(:,0,:,:) ,False,True,levwb,True,jpeg)
AnnPlot(pltDir,"ann_wbal"        ,numYr,type,"BlueDarkRed18",wbal(:,0,:,:) ,False,True,levwb,False,jpeg)
end if
Panel_Plot(pltDir,"seas_srnf"    ,numYr,type,cmaprn,srnf(:,0,:,:) ,False,True,levrnf,True,jpeg)
end if

if(isvar("frozf")) then
do n = 0,nsoil-1
e = n + 1
frozf@long_name = "Frozen soil moisture fraction: Level "+e
Panel_Plot(pltDir,"seas_frozf"+e,numYr,type,"BlAqGrYeOrReVi200",frozf(:,n,:,:),False,True,levsfr,True,False)
end do
end if

if(isvar("smcl")) then
do n = 0,nsoil-1
f = n + 1
smcl@long_name = "Soil moisture in each (of "+nsoil+") levels: Level "+f
Panel_Plot(pltDir,"seas_smc"+f  ,numYr,type,cmapsm,smcl(:,n,:,:) ,False,True,levsm,True,False)
end do
Panel_Plot(pltDir,"seas_smc" ,numYr,type,cmapsm,smc(:,:,:) ,False,True,levsc,True,False)
end if 

if(isvar("gpp")) then
levg = (/0,0.00000002,0.00000004,0.00000006,0.00000008,0.00000010,.00000012/)
levn = (/-.00000002,-.00000001,0,.00000001,.00000002,.00000003,.00000004,.00000005,.00000006,.00000007/)
levp = (/0,.000000004,.000000008,.000000012,.000000016,.000000020,.000000024,.000000028,.000000032/)
levs = (/0,.000000004,.000000008,.000000012,.000000016,.000000020,.000000024,.000000028,.000000032,.000000036,.00000004/)
Panel_Plot(pltDir,"seas_gpp" ,numYr,type,"BlAqGrYeOrReVi200",gpp(:,0,:,:) ,False,True,levg,True,False)
Panel_Plot(pltDir,"seas_npp" ,numYr,type,"BlAqGrYeOrReVi200",npp(:,0,:,:) ,False,True,levn,True,False)
Panel_Plot(pltDir,"seas_frp" ,numYr,type,"BlAqGrYeOrReVi200",frp(:,0,:,:) ,False,True,levp,True,False)
Panel_Plot(pltDir,"seas_frs" ,numYr,type,"BlAqGrYeOrReVi200",frs(:,0,:,:) ,False,True,levs,True,False)
end if

Panel_Plot(pltDir,"seas_cloudalb",numYr,type,"BlAqGrYeOrReVi200",calb(:,0,:,:),False,True,levsfr,False,jpeg)
Panel_Plot(pltDir,"seas_surfalb" ,numYr,type,"BlAqGrYeOrReVi200",surfalb(:,0,:,:),False,True,levsal,False,jpeg)
Panel_Plot(pltDir,"seas_swnet"   ,numYr,type,"BlAqGrYeOrReVi200",swnet(:,0,:,:),False,True,levsw,False,jpeg)
Panel_Plot(pltDir,"seas_seb" ,numYr,type,"BlueDarkRed18",seb(:,0,:,:)  ,False,True,levseb,True,jpeg)
; annual SEB ----------
AnnPlot(pltDir,"ann_seb"      ,numYr,type,"BlueDarkRed18",seb(:,0,:,:)  ,False,True,levseb,True,jpeg)

; soil temp -----------
;if (Runid.eq."xaank".or.Runid.eq."xaanv".and.lcable) then
if (isvar("deepST")) then
do n = 0,nsoil-1,1
f = n + 1
deepST@long_name = "Soil temperature: Level "+f
Panel_Plot(pltDir,"seas_tsl"+f,numYr,type,"BlAqGrYeOrReVi200",deepST(:,n,:,:),False,True,levtp,False,jpeg)
end do
end if
if (lcable.and.t15) then
do n = 0,nsoil-1,1
f = n + 1
Panel_Plot(pltDir,"seas_soilTcab"+f,numYr,type,"BlAqGrYeOrReVi200",soilTemp(:,n,:,:),False,True,levtp,False,jpeg)
end do
end if

; Plot Biases ======================================================

print("Plotting Biases") 
; color map -> BlwhRe

lev_bT   = ispan(-9,9,1)
;lev_bT   = (/-8,-4,-2,-1,0,1,2,4,8/)
lev_bP   = (/-8,-4,-2,-1,-0.5,0,0.5,1,2,4,8/)
lev_bE   = (/-1500,-1000,-500,-400,-300,-200,-100,-50,0,50,100,200,300,400,500,1000,1500/)

;cmap_bT = RGBtoCmap("/home/cmar/ste69f/umplot/colormap/BlueOrangeRed.rgb")

Panel_Plot(pltDir,"bias_tscrn" ,numYr,type,"BlueDarkRed18",bias_tscrn ,False,True,lev_bT,False,jpeg)
;Panel_Plot(pltDir,"bias_tscrn" ,numYr,type,"BlueDarkRed18",bias_tscrn ,False,True,lev_bT,False,jpeg)
AnnPlot(pltDir,"annbias_tscrn" ,numYr,type,"BlueDarkRed18",bias_tscrn ,False,True,lev_bT,False,jpeg)
Panel_Plot(pltDir,"bias_precip",numYr,type,"BlueDarkRed18",bias_precip,False,True,lev_bP,False,jpeg)
AnnPlot(pltDir,"annbias_precip",numYr,type,"BlueDarkRed18",bias_precip,False,True,lev_bP,False,jpeg)
if(isvar("tmax")) then
Panel_Plot(pltDir,"bias_tmax"  ,numYr,type,"BlueDarkRed18",bias_tmax  ,False,True,lev_bT,True,jpeg)
Panel_Plot(pltDir,"bias_tmin"  ,numYr,type,"BlueDarkRed18",bias_tmin  ,False,True,lev_bT,True,jpeg)
Panel_Plot(pltDir,"bias_tamp"  ,numYr,type,"BlueDarkRed18",bias_tamp  ,False,True,lev_bT,True,jpeg)
end if
if(isvar("evap")) then
Panel_Plot(pltDir,"bias_evap"  ,numYr,type,"BlueDarkRed18",bias_evap  ,False,True,lev_bE,True,jpeg)
end if

; Australian Region ===============================================

print("Plotting Australian Region")

lev_aT   = (/240,250,260,270,275,280,285,290,295,300,305,310,315/)
lev_aSnw = (/0,5,10,20,40,60,80,100,150,200,400,800,1000,1500/) 

Panel_Plot(pltDir,"aust_tscrn" ,numYr,type,"BlAqGrYeOrReVi200",tscrn(:,0,:,:) ,True,True,levtp,True,jpeg)
Panel_Plot(pltDir,"aust_precip",numYr,type,cmapr2,precip(:,0,:,:),True,True,levbpr,True,jpeg)

IndPlot(pltDir,"apr_djf",numYr,type,cmapr2,precip(0,0,:,:),True,True,levbpr,True,jpeg)
IndPlot(pltDir,"apr_mam",numYr,type,cmapr2,precip(1,0,:,:),True,True,levbpr,True,jpeg)
IndPlot(pltDir,"apr_jja",numYr,type,cmapr2,precip(2,0,:,:),True,True,levbpr,True,jpeg)
IndPlot(pltDir,"apr_son",numYr,type,cmapr2,precip(3,0,:,:),True,True,levbpr,True,jpeg)

if(isvar("tmax")) then
Panel_Plot(pltDir,"aust_tmax"  ,numYr,type,"BlAqGrYeOrReVi200",tmax(:,0,:,:)  ,True,True,levtp,True,jpeg)
Panel_Plot(pltDir,"aust_tmin"  ,numYr,type,"BlAqGrYeOrReVi200",tmin(:,0,:,:)  ,True,True,levtp,True,jpeg)
Panel_Plot(pltDir,"aust_tamp"  ,numYr,type,"BlAqGrYeOrReVi200",tamp(:,0,:,:)  ,True,True,levam,True,jpeg)
end if

if(isvar("lai")) then
;precip2_17lev
do n =0,npft-1
b = n+1
lai@long_name   = "Leaf Area Index for PFT "+b
Panel_Plot(pltDir,"aust_lai"+b,numYr,type,"BlAqGrYeOrReVi200",lai(:,n,:,:),True,True,levlai,True,False)
end do

levht1 = (/2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34/)
levht2 = (/0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0/)
levht3 = (/0.5,1.0,1.5,2.0,2.5,3.0,3.5/)

if (npft.eq.5) then
canht@long_name = "Canopy Height for PFT 1"
Panel_Plot(pltDir,"aust_canht1" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,0,:,:) ,True,True,levdf,True,False)
canht@long_name = "Canopy Height for PFT 2"
Panel_Plot(pltDir,"aust_canht2" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,1,:,:) ,True,True,levdf,True,False)
canht@long_name = "Canopy Height for PFT 3"
Panel_Plot(pltDir,"aust_canht3" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,2,:,:) ,True,True,levht2,True,False)
canht@long_name = "Canopy Height for PFT 4"
Panel_Plot(pltDir,"aust_canht4" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,3,:,:) ,True,True,levht2,True,False)
canht@long_name = "Canopy Height for PFT 5"
Panel_Plot(pltDir,"aust_canht5" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,4,:,:) ,True,True,levht3,True,False)
end if
if (npft.eq.13) then
canht@long_name = "Canopy Height for PFT 1"
Panel_Plot(pltDir,"aust_canht1" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,0,:,:) ,True,True,levht1,True,False)
canht@long_name = "Canopy Height for PFT 2"
Panel_Plot(pltDir,"aust_canht2" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,1,:,:) ,True,True,levht1,True,False)
canht@long_name = "Canopy Height for PFT 3"
Panel_Plot(pltDir,"aust_canht3" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,2,:,:) ,True,True,levht1,True,False)
canht@long_name = "Canopy Height for PFT 4"
Panel_Plot(pltDir,"aust_canht4" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,3,:,:) ,True,True,levht1,True,False)
canht@long_name = "Canopy Height for PFT 5"
Panel_Plot(pltDir,"aust_canht5" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,4,:,:) ,True,True,levht2,True,False)
canht@long_name = "Canopy Height for PFT 6"
Panel_Plot(pltDir,"aust_canht6" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,5,:,:) ,True,True,levht2,True,False)
canht@long_name = "Canopy Height for PFT 7"
Panel_Plot(pltDir,"aust_canht7" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,6,:,:) ,True,True,levht2,True,False)
canht@long_name = "Canopy Height for PFT 8"
Panel_Plot(pltDir,"aust_canht8" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,7,:,:) ,True,True,levht2,True,False)
canht@long_name = "Canopy Height for PFT 9"
Panel_Plot(pltDir,"aust_canht9" ,numYr,type,"BlAqGrYeOrReVi200",canht(:,8,:,:) ,True,True,levht2,True,False)
canht@long_name = "Canopy Height for PFT 10"
Panel_Plot(pltDir,"aust_canht10",numYr,type,"BlAqGrYeOrReVi200",canht(:,9,:,:) ,True,True,levht2,True,False)
canht@long_name = "Canopy Height for PFT 11"
Panel_Plot(pltDir,"aust_canht11",numYr,type,"BlAqGrYeOrReVi200",canht(:,10,:,:),True,True,levht2,True,False)
canht@long_name = "Canopy Height for PFT 12"
Panel_Plot(pltDir,"aust_canht12",numYr,type,"BlAqGrYeOrReVi200",canht(:,11,:,:),True,True,levht2,True,False)
canht@long_name = "Canopy Height for PFT 13"
Panel_Plot(pltDir,"aust_canht13",numYr,type,"BlAqGrYeOrReVi200",canht(:,12,:,:),True,True,levht2,True,False)
end if
end if ; if lai
 
if (isvar("lalb")) then
Panel_Plot(pltDir,"aust_lalb"  ,numYr,type,"BlAqGrYeOrReVi200",lalb(:,0,:,:)  ,True,True,levalb,True,jpeg)
end if
if(isvar("lsftp")) then
Panel_Plot(pltDir,"aust_lsftp" ,numYr,type,"BlAqGrYeOrReVi200",lsftp(:,0,:,:) ,True,True,levtp,True,jpeg)
end if

if(isvar("sfrac")) then
do n = 0,ntiles-1
d = n + 1
Panel_Plot(pltDir,"aust_sfrac"+d,numYr,type,"BlAqGrYeOrReVi200",sfrac(:,n,:,:),True,True,levsfr,True,False)
end do
end if

if(isvar("sndep")) then
Panel_Plot(pltDir,"aust_sndep" ,numYr,type,"BlAqGrYeOrReVi200",sndep(:,0,:,:) ,True,True,lev_aSnw,True,jpeg)
end if
if(isvar("nsradf")) then
Panel_Plot(pltDir,"aust_nsradf",numYr,type,"BlAqGrYeOrReVi200",nsradf(:,0,:,:),True,True,levrad,True,jpeg)
end if

if(isvar("evap")) then
Panel_Plot(pltDir,"aust_evap"  ,numYr,type,cmapsn,evap(:,0,:,:)  ,True,True,levevp,True,jpeg)
end if
if(isvar("tsurf")) then
Panel_Plot(pltDir,"aust_tsurf" ,numYr,type,"BlAqGrYeOrReVi200",tsurf(:,0,:,:) ,True,True,levtp,True,jpeg)
end if
Panel_Plot(pltDir,"aust_sh"    ,numYr,type,"BlAqGrYeOrReVi200",sh(:,0,:,:)  ,True,True,levsh,True,jpeg)
Panel_Plot(pltDir,"aust_lh"    ,numYr,type,"BlAqGrYeOrReVi200",lh(:,0,:,:)  ,True,True,levlh,True,jpeg)
;Panel_Plot(pltDir,"aust_seb"   ,numYr,type,"BlAqGrYeOrReVi200",seb(:,0,:,:) ,True,True,levseb,True,jpeg)
; annual SEB
AnnPlot(pltDir,"austann_seb"   ,numYr,type,"BlueDarkRed18",seb(:,0,:,:) ,True,True,levseb,True,jpeg)
;AnnPlot(pltDir,"austann_wbal"  ,numYr,type,"BlueDarkRed18",wbal(:,0,:,:),True,True,levwb,True,jpeg)
;Panel_Plot(pltDir,"aust_precip_bom"   ,numYr,type,cmapr2,bmpr(:,:,:) ,True,True,levbpr,True,jpeg)

IndPlot(pltDir,"bom_pr_djf","101",type,cmapr2,bmpr(0,:,:),True,True,levbpr,True,jpeg)
IndPlot(pltDir,"bom_pr_mam","101",type,cmapr2,bmpr(1,:,:),True,True,levbpr,True,jpeg)
IndPlot(pltDir,"bom_pr_jja","101",type,cmapr2,bmpr(2,:,:),True,True,levbpr,True,jpeg)
IndPlot(pltDir,"bom_pr_son","101",type,cmapr2,bmpr(3,:,:),True,True,levbpr,True,jpeg)

; Australian Biases ===============================================

print("Plotting Australian Biases")
; color map -> BlwhRe

;lev_abT   = (/-8,-4,-2,-1,0,1,2,4,8/)
lev_bP   = (/-8,-4,-2,-1,-0.5,0,0.5,1,2,4,8/)
lev_bP2  = (/-8,-7,-6,-5,-4,-3,-2,-1,-0.5,0,0.5,1,2,3,4,5,6,7,8/)
lev_bE   = (/-1500,-1000,-500,-400,-300,-200,-100,-50,0,50,100,200,300,400,500,1000,1500/)

Panel_Plot(pltDir,"austbias_tscrn" ,numYr,type,"BlueDarkRed18",bias_tscrn ,True,True,lev_bT,True,jpeg)
Panel_Plot(pltDir,"austbias_precip",numYr,type,"BlueDarkRed18",bias_precip,True,True,lev_bP2,True,jpeg)
if(isvar("tmax")) then
Panel_Plot(pltDir,"austbias_tmax"  ,numYr,type,"BlueDarkRed18",bias_tmax  ,True,True,lev_bT,True,jpeg)
Panel_Plot(pltDir,"austbias_tmin"  ,numYr,type,"BlueDarkRed18",bias_tmin  ,True,True,lev_bT,True,jpeg)
Panel_Plot(pltDir,"austbias_tamp"  ,numYr,type,"BlueDarkRed18",bias_tamp  ,True,True,lev_bT,True,jpeg)
end if
if(isvar("evap")) then
Panel_Plot(pltDir,"austbias_evap"  ,numYr,type,"BlueDarkRed18",bias_evap  ,True,True,lev_bE,True,jpeg)
end if
Panel_Plot(pltDir,"austbias_precip_bom" ,numYr,type,"BlueDarkRed18",bias_bmpr ,True,True,lev_bP2,True,jpeg)

; Jan and Jul Monthly Means =======================================

print("Plotting Jan/Jul Means")

MmnPanel(pltDir,"mmn_precip" ,numYr,type,cmapr,prmn(:,0,:,:) ,False,True,levpr,False)
MmnPanel(pltDir,"mmn_tscrn"  ,numYr,type,"BlAqGrYeOrReVi200",t2mn(:,0,:,:) ,False,True,levtp,False)
if(isvar("tmax")) then
MmnPanel(pltDir,"mmn_tmax"   ,numYr,type,"BlAqGrYeOrReVi200",txmn(:,0,:,:) ,False,True,levtp,False)
MmnPanel(pltDir,"mmn_tmin"   ,numYr,type,"BlAqGrYeOrReVi200",tnmn(:,0,:,:) ,False,True,levtp,False)
MmnPanel(pltDir,"mmn_tamp"   ,numYr,type,"BlAqGrYeOrReVi200",tamn(:,0,:,:) ,False,True,levam,False)
end if
if(isvar("lalb")) then
MmnPanel(pltDir,"mmn_lalb"   ,numYr,type,"BlAqGrYeOrReVi200",albmn(:,0,:,:),False,True,levalb,False)
end if

if(isvar("laimn")) then

do n =0,npft-1
d = n+1
lai@long_name   = "Leaf Area Index for PFT "+d
MmnPanel(pltDir,"mmn_lai"+d  ,numYr,type,"BlAqGrYeOrReVi200",laimn(:,n,:,:),False,True,levlai,False)
end do

do n =0,npft-1
g = n+1
chtmn@long_name = "Canopy Height for PFT "+g
MmnPanel(pltDir,"mmn_canht"+g,numYr,type,"BlAqGrYeOrReVi200",chtmn(:,n,:,:),False,True,levdf,False)
end do
if(isvar("lstmn")) then
MmnPanel(pltDir,"mmn_lsftp"  ,numYr,type,"BlAqGrYeOrReVi200",lstmn(:,0,:,:),False,True,levtp,False)
end if
end if

if(isvar("sndmn")) then
MmnPanel(pltDir,"mmn_sndep"  ,numYr,type,cmapsn,sndmn(:,0,:,:),False,True,levsnw,False)
end if

if(isvar("srmn")) then
MmnPanel(pltDir,"mmn_srnf"   ,numYr,type,cmaprn,srmn(:,0,:,:) ,False,True,levrnf,False)
MmnPanel(pltDir,"mmn_drnf"   ,numYr,type,cmaprn,drmn(:,0,:,:) ,False,True,levrnf,False)
end if

; soil temp -----------
;;if (Runid.eq."xaank".and.lcable) then
;if (isvar("deepST")) then
;do i = 0,nsoil-1,1
;MmnPanel(pltDir,"mmn_tsl",numYr,type,"BlAqGrYeOrReVi200",deepST(:,n,:,:),False,True,levtp,False);,jpeg)
;end do
;end if
;if (lcable.and.t15) then
;do i = 0,nsoil-1,1
;MmnPanel(pltDir,"mmn_soilTcab",numYr,type,"BlAqGrYeOrReVi200",soilTemp(:,n,:,:),False,True,levtp,False);,jpeg)
;end do
;end if

; Map Projections =================================================

print("Plotting Map Projections")
print("")

nHem = "NH"
sHem = "SH"
Polar_Plot(pltDir,"mp_tscrn_"+nHem,numYr,type,"perc2_9lev",tscrn(:,0,:,:),False,True,levtp,nHem,False)
Polar_Plot(pltDir,"mp_tscrn_"+sHem,numYr,type,"perc2_9lev",tscrn(:,0,:,:),False,True,levtp,sHem,False)
Polar_Plot(pltDir,"mpbias_tscrn_"+nHem,numYr,type,"BlueDarkRed18",bias_tscrn,False,True,lev_bT,nHem,False)
Polar_Plot(pltDir,"mpbias_tscrn_"+sHem,numYr,type,"BlueDarkRed18",bias_tscrn,False,True,lev_bT,sHem,False)

if(isvar("tmax")) then
Polar_Plot(pltDir,"mp_tmax_"+nHem,numYr,type,"perc2_9lev",tmax(:,0,:,:),False,True,levtp,nHem,False)
Polar_Plot(pltDir,"mp_tmax_"+sHem,numYr,type,"perc2_9lev",tmax(:,0,:,:),False,True,levtp,sHem,False)
Polar_Plot(pltDir,"mpbias_tmax_"+nHem,numYr,type,"BlueDarkRed18",bias_tmax,False,True,lev_bT,nHem,False)
Polar_Plot(pltDir,"mpbias_tmax_"+sHem,numYr,type,"BlueDarkRed18",bias_tmax,False,True,lev_bT,sHem,False)

Polar_Plot(pltDir,"mp_tmin_"+nHem,numYr,type,"perc2_9lev",tmin(:,0,:,:),False,True,levtp,nHem,False)
Polar_Plot(pltDir,"mp_tmin_"+sHem,numYr,type,"perc2_9lev",tmin(:,0,:,:),False,True,levtp,sHem,False)
Polar_Plot(pltDir,"mpbias_tmin_"+nHem,numYr,type,"BlueDarkRed18",bias_tmin,False,True,lev_bT,nHem,False)
Polar_Plot(pltDir,"mpbias_tmin_"+sHem,numYr,type,"BlueDarkRed18",bias_tmin,False,True,lev_bT,sHem,False)
end if

Polar_Plot(pltDir,"mp_precip_"+nHem,numYr,type,cmapr,precip(:,0,:,:),False,True,levpr,nHem,False)
Polar_Plot(pltDir,"mp_precip_"+sHem,numYr,type,cmapr,precip(:,0,:,:),False,True,levpr,sHem,False)
Polar_Plot(pltDir,"mpbias_precip_"+nHem,numYr,type,"BlueDarkRed18",bias_precip,False,True,lev_bP,nHem,False)
Polar_Plot(pltDir,"mpbias_precip_"+sHem,numYr,type,"BlueDarkRed18",bias_precip,False,True,lev_bP,sHem,False)
; annual ?

; Global Wbal Table ===============================================

; MMDC ============================================================

; Jaeger Plots ====================================================

; =================================================================

;if (Taylor) then

; Taylor Diagram ==================================================

; MMDC and Fluxnet ================================================

; Bar Charts ======================================================

; Bias 3panel =====================================================

; InterAnnual Variability =========================================

;IntVar_Mon(pltDir,"Int",numYr,type,)
;IntVar_Yr(pltDir,"Int",numYr,type,)
;Aust_IntVar_Mon()
;Aust_IntVar_Yr()

;==================================================================

; end if ; Taylor 
 
;print("Finished Plotting")

if(.not.Split)then
system("mv *.ps ./plots")
if(jpeg)then
system("mv *.jpg ./plots_jpg")
end if
;if(eps)then
;system("mv *.eps ./plots_eps")
;end if
end if

end
