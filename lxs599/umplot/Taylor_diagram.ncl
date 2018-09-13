; Lauren Stevens 22 Aug 2011 (Ferret - based on a script ~IPCC/data/other/taylor* - Mark Collier)
; Taylor Diagram using NCL 
;================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
hdir = systemfunc("echo $hdir")
usr = systemfunc("echo $USERID")
loadscript(hdir+"/"+usr+"/umplot/taylor_diagram.ncl")
;================================================================

begin

;usr = systemfunc("echo $USERID")
fileDir =systemfunc("echo $DIR")
;Runid=systemfunc("echo $RUNID")
numYr=systemfunc("echo $YR")
;CABLE=systemfunc("echo $MID")
;MOSES=systemfunc("echo $TOP")
mossea=systemfunc("echo $mossea")
cabsea=systemfunc("echo $cabsea")
mostem=systemfunc("echo $mostem")
cabtem=systemfunc("echo $cabtem")
tname =  systemfunc("echo $tname")
txname =  systemfunc("echo $txname")
tiname =  systemfunc("echo $tiname")
tstep =  systemfunc("echo $TSTEP")

;FILES-----------------------------------------------------------------------------------------------------

cable_file = addfile(cabsea,"r")
moses_file = addfile(mossea,"r")
if (isfilepresent(""+cabtem+"")) then
cable_tmp  = addfile(cabtem,"r")
end if
if (isfilepresent(""+mostem+"")) then
moses_tmp  = addfile(mostem,"r")
end if

era_tp = addfile(hdir+"/"+usr+"/umplot/obs/ERAi_tp2m.nc","r")         
era_pr = addfile(hdir+"/"+usr+"/umplot/obs/ERAi_precip.nc","r")
era_tx = addfile(hdir+"/"+usr+"/umplot/obs/ERAi_tmax_tmin.nc","r")
era_sw = addfile(hdir+"/"+usr+"/umplot/obs/ERA_INT_rss_ann.nc","r")    ; sw era int obs
era_lw = addfile(hdir+"/"+usr+"/umplot/obs/ERA_INT_rls_ann.nc","r")    ; lw era int obs
era_lh = addfile(hdir+"/"+usr+"/umplot/obs/ERA_INT_hfls_ann.nc","r") ; lh era int obs
era_sh = addfile(hdir+"/"+usr+"/umplot/obs/ERA_INT_hfss_ann.nc","r") ; sh era int obs

;----------------------------------------------------------------------------------------------------------

if (isfilevar(cable_file,"tscrn")) then
ctemp = cable_file ->tscrn
else
if (isfilevar(cable_file,"tas")) then
ctemp = cable_file ->tas
else
ctemp = cable_file ->$tname$
end if
end if
if (isvar("cable_tmp")) then
if (isfilevar(cable_tmp,"tmax")) then
 ctmax = cable_tmp  ->tmax
 ctmin = cable_tmp  ->tmin
else
 ctmax = cable_tmp  ->$txname$ ;temp
 ctmin = cable_tmp  ->$tiname$ ;temp_1
end if
end if
cprec = cable_file ->precip
clh   = cable_file ->lh
csh   = cable_file ->sh
csw   = cable_file ->solar
clw   = cable_file ->longwave

if (dimsizes(dimsizes(ctmax)).eq.4) then
 ctmax!2 = "latitude"
 ctmax!3 = "longitude"
 ctmin!2 = "latitude"
 ctmin!3 = "longitude"
else
 ctmax!1 = "latitude"
 ctmax!2 = "longitude"
 ctmin!1 = "latitude"
 ctmin!2 = "longitude"
end if

if (isfilevar(moses_file,"tscrn")) then
mtemp = moses_file ->tscrn
else
if (isfilevar(moses_file,"tas")) then
mtemp = moses_file ->tas
else
mtemp = moses_file ->$tname$
end if
end if
if (isvar("moses_tmp")) then
if (isfilevar(moses_tmp,"tmax")) then
 mtmax = moses_tmp  ->tmax
 mtmin = moses_tmp  ->tmin
else
 mtmax = moses_tmp  ->$txname$ ;temp
 mtmin = moses_tmp  ->$tiname$ ;temp_1
end if
end if
if (isfilevar(moses_file,"precip")) then
mprec = moses_file ->precip
else
mprec = moses_file ->pr
end if
if (isfilevar(moses_file,"lh")) then
mlh   = moses_file ->lh
else
mlh   = moses_file ->hfls
end if
if (isfilevar(moses_file,"sh")) then
msh   = moses_file ->sh
else
msh   = moses_file ->hfss
end if
if (isfilevar(moses_file,"solar")) then
msw   = moses_file ->solar
else
msw   = moses_file ->rss
end if
if (isfilevar(moses_file,"longwave")) then
mlw   = moses_file ->longwave
else
mlw   = moses_file ->rls
end if

tp2m  = era_tp ->a
prec  = era_pr ->a
tmax  = era_tx ->mx2t
tmin  = era_tx ->mn2t
sh    = era_sh ->hfss_ann
lh    = era_lh ->hfls_ann
lw    = era_lw ->rls_ann
sw    = era_sw ->rss_ann

tp2m        = tp2m + 273.15
tp2m@units  = "deg K" 
cprec       = cprec*24*60*60
cprec@units = "mm/day"
mprec       = mprec*24*60*60
mprec@units = "mm/day"

pi       = 4.0*atan(1.0)         ; 3.141592654
if (isdim(mtemp,"latitude")) then
lat      = mtemp&latitude
else
lat      = mtemp&lat
end if
wt       = cos(lat*pi/180)
wgt      = conform(ctemp(0,0,:,:),wt,0)
;wgt      = new((/73,96/),double)
;wgt      = 1.0
;do i = 0,95,1
;wgt(:,i) = cos(lat*pi/180)
;end do
;print(wgt)

ratio            = new((/2,8/),double)
cc               = new((/2,8/),double)
ratio            = -9999
cc               = -9999
ratio@_FillValue = -9999
cc@_FillValue    = -9999

;INITIALISATION--------------------------------------------------------------------------------------------
;STATS-----------------------------------------------------------------------------------------------------

function taylor_stats(wgt,obs,mod)
begin
sumw    = sum(wgt)
sumwo   = sum(wgt*obs)
sumwm   = sum(wgt*mod)
wmean_o = sumwo/sumw
wmean_m = sumwm/sumw
;wvar_o  = sum(wgt*(obs-wmean_o)^2)/sumw                                    ; NCL example
;wvar_m  = sum(wgt*(mod-wmean_m)^2)/sumw                                    ; NCL example
wvar_o  = sum((obs-wmean_o)^2)/sumw                                         ; M.Collier ~IPCC/data/other
wvar_m  = sum((mod-wmean_m)^2)/sumw                                         ; M.Collier ~IPCC/data/other
;wcc     = (sum(wgt*obs*mod) - sumwo*sumwm/sumw )/    \                     ; NCL example
;          ((sum(wgt*obs^2) - sumwo^2/sumw) *         \                     ; NCL example
;          (sum(wgt*mod^2) - sumwm^2/sumw) )^0.5                            ; NCL example
wcc     = (sum((obs-wmean_o)^2)+sum((mod-wmean_m)^2)-(sum(((obs-wmean_o)- \ ; M.Collier ~IPCC/data/other
          (mod-wmean_m))^2)))/((4*sum((obs-wmean_o)^2)* \                   ; M.Collier ~IPCC/data/other
          sum((mod-wmean_m)^2))^0.5)                                        ; M.Collier ~IPCC/data/other
wratio  = (wvar_m/wvar_o)^0.5
stats   = (/wratio,wcc/)
;print(sumw)
;print(sumwo)
;print(sumwm)
;print(wmean_o)
;print(wmean_m)
;print(wvar_o)
;print(wvar_m)
;print(wcc)
;print(wratio)
return(stats)
end

;CABLE--------------------------------------------------------------------------------------

; SCREEN TEMP -------------------------
  dNames = getvardims(ctemp)
  nDims  = dimsizes(ctemp)
  nSize  = dimsizes(nDims)
  lon    = ctemp&$dNames(nSize-1)$
tobs   = linint2(tp2m&longitude,tp2m&latitude,tp2m,True,lon,ctemp&latitude,0)
tvar_o = dim_avg_n_Wrap(tobs,0)
tvar_m = dim_avg_n_Wrap(ctemp(:,0,:,:),0)
ctscrn = taylor_stats(wgt,tvar_o,tvar_m)
ratio(0,0) = ctscrn(0)
cc(0,0)    = ctscrn(1)

;cvar_c = variance(tvar_o)
;cvar_t = variance(tvar_m)
;cstd_c = stddev()
;cstd_t = stddev()
;ccorr  = pattern_cor(tvar_o,tvar_m,wgt,0) 
;crmse  = wgt_arearmse(tvar_o,tvar_m,wgt,1.0,1)
;oratio = (cvar_t/cvar_c)^0.5
;print(cvar_t)
;print(cvar_c)
;print(oratio)

; PRECIP-------------------------------
  pNames = getvardims(cprec)
  pDims  = dimsizes(cprec)
  pSize  = dimsizes(nDims)
  lonp   = cprec&$dNames(nSize-1)$
;var_o = a[d=2,gx=field77[d=4],gy=field77[d=4],l=1:4@sum]/4
;var_m = (precip[d=4,l=1:4@sum]/4)*24*60*60 ;(field77[d=4,l=1:4@sum]/4)*tstep
pobs   = linint2(prec&longitude,prec&latitude,prec,True,lonp,cprec&latitude,0)
pvar_o = dim_avg_n_Wrap(pobs,0)
pvar_m = dim_avg_n_Wrap(cprec(:,0,:,:),0)
cprecip = taylor_stats(wgt,pvar_o,pvar_m)
ratio(0,1) = cprecip(0)
cc(0,1)    = cprecip(1)

; TMAX--------------------------------

if (isvar("cable_tmp")) then

;var_o = mx2t[d=3,gx=temp[d=6],gy=temp[d=6],l=1:4@sum]/4
;var_m = temp[d=6,l=1:4@sum]/4
xobs   = linint2(tmax&longitude,tmax&latitude,tmax,True,ctmax&longitude,ctmax&latitude,0)
xvar_o = dim_avg_n_Wrap(xobs,0)
xvar_m = dim_avg_n_Wrap(ctmax(:,0,:,:),0)
c_tmax  = taylor_stats(wgt,xvar_o,xvar_m)
ratio(0,2) = c_tmax(0)
cc(0,2)    = c_tmax(1)

; TMIN--------------------------------
;var_o = mn2t[d=3,gx=temp_1[d=6],gy=temp_1[d=6],l=1:4@sum]/4
;var_m = temp_1[d=6,l=1:4@sum]/4
nobs   = linint2(tmin&longitude,tmin&latitude,tmin,True,ctmin&longitude,ctmin&latitude,0)
nvar_o = dim_avg_n_Wrap(nobs,0)
nvar_m = dim_avg_n_Wrap(ctmin(:,0,:,:),0)
c_tmin  = taylor_stats(wgt,nvar_o,nvar_m)
ratio(0,3) = c_tmin(0)
cc(0,3)    = c_tmin(1)

else

ratio(0,2:3) = -9999
cc(0,2:3)    = -9999

end if

; SW RAD------------------------------
  sNames = getvardims(csw)
  sDims  = dimsizes(csw)
  sSize  = dimsizes(nDims)
  lons   = csw&$dNames(nSize-1)$
;var_o = RSS_ANN[d=8,gx=SOLAR[d=4],gy=SOLAR[d=4],l=1];:4@ave]
;var_m = SOLAR[d=4,l=1:4@sum]/4
sobs   = linint2(sw&lon,sw&lat,sw,True,lons,csw&latitude,0)
svar_o = dim_avg_n_Wrap(sobs,0)
svar_m = dim_avg_n_Wrap(csw(:,0,:,:),0)
c_sw    = taylor_stats(wgt,svar_o,svar_m)
ratio(0,4) = c_sw(0)
cc(0,4)    = c_sw(1)

; LW RAD------------------------------
  lNames = getvardims(clw)
  lDims  = dimsizes(clw)
  lSize  = dimsizes(nDims)
  lonl   = clw&$dNames(nSize-1)$
;var_o = RLS_ANN[d=9,gx=LONGWAVE[d=4],gy=LONGWAVE[d=4],l=1];:4@ave]
;var_m = LONGWAVE[d=4,l=1:4@sum]/4
lobs   = linint2(lw&lon,lw&lat,lw,True,lonl,clw&latitude,0)
lvar_o = dim_avg_n_Wrap(lobs,0)
lvar_m = dim_avg_n_Wrap(clw(:,0,:,:),0)
c_lw    = taylor_stats(wgt,lvar_o,lvar_m)
ratio(0,5) = c_lw(0)
cc(0,5)    = c_lw(1)

; SENS HEAT---------------------------
  hNames = getvardims(csh)
  hDims  = dimsizes(csh)
  hSize  = dimsizes(nDims)
  lonh   = csh&$dNames(nSize-1)$
;var_o = HFSS_ANN[d=11,gx=SH[d=4],gy=SH[d=4],l=1];:4@ave]
;var_m = SH[d=4,l=1:4@sum]/4
shobs   = linint2(sh&lon,sh&lat,sh,True,lonh,csh&latitude,0)
shvar_o = dim_avg_n_Wrap(shobs,0)
shvar_m = dim_avg_n_Wrap(csh(:,0,:,:),0)
c_sh     = taylor_stats(wgt,shvar_o,shvar_m)
ratio(0,6) = c_sh(0)
cc(0,6)    = c_sh(1)

; LATENT HEAT-------------------------
  tNames = getvardims(clh)
  tDims  = dimsizes(clh)
  tSize  = dimsizes(nDims)
  lont   = clh&$dNames(nSize-1)$
;var_o = HFLS_ANN[d=10,gx=LH[d=4],gy=LH[d=4],l=1];:4@ave]
;var_m = LH[d=4,l=1:4@sum]/4
lhobs   = linint2(lh&lon,lh&lat,lh,True,lont,clh&latitude,0)
lhvar_o = dim_avg_n_Wrap(lhobs,0)
lhvar_m = dim_avg_n_Wrap(clh(:,0,:,:),0)
c_lh     = taylor_stats(wgt,lhvar_o,lhvar_m)
ratio(0,7) = c_lh(0)
cc(0,7)    = c_lh(1)

;MOSES---------------------------------------------------------------------------------------

; SCREEN TEMP ------------------------
if (isdim(mtemp,"latitude")) then
mtobs   = linint2(tp2m&longitude,tp2m&latitude,tp2m,True,lon,mtemp&latitude,0)
else
mtobs   = linint2(tp2m&longitude,tp2m&latitude,tp2m,True,lon,mtemp&lat,0)
end if
mtvar_o = dim_avg_n_Wrap(mtobs,0)
if (dimsizes(dimsizes(mtemp)).eq.4) then
mtvar_m = dim_avg_n_Wrap(mtemp(:,0,:,:),0) 
else
mtvar_m = dim_avg_n_Wrap(mtemp(:,:,:),0) 
end if
;mcorr  = pattern_cor(mtvar_o,mtvar_m,wgt,0)
;mrmse  = wgt_arearmse(mtvar_o,mtvar_m,wgt,1.0,1)
mtscrn = taylor_stats(wgt,mtvar_o,mtvar_m)
ratio(1,0) = mtscrn(0)
cc(1,0)    = mtscrn(1)

; PRECIP------------------------------
;var_o = a[d=2,gx=field77[d=5],gy=field77[d=5],l=1:4@sum]/4
;var_m = (precip[d=5,l=1:4@sum]/4)*24*60*60 ;(field77[d=5,l=1:4@sum]/4)*tstep
if (isdim(mprec,"latitude")) then
mpobs   = linint2(prec&longitude,prec&latitude,prec,True,lonp,mprec&latitude,0)
else
mpobs   = linint2(prec&longitude,prec&latitude,prec,True,lonp,mprec&lat,0)
end if
mpvar_o = dim_avg_n_Wrap(mpobs,0)
if (dimsizes(dimsizes(mprec)).eq.4) then
mpvar_m = dim_avg_n_Wrap(mprec(:,0,:,:),0)
else
mpvar_m = dim_avg_n_Wrap(mprec(:,:,:),0)
end if
mprecip = taylor_stats(wgt,mpvar_o,mpvar_m)
ratio(1,1) = mprecip(0)
cc(1,1)    = mprecip(1)

; TMAX-------------------------------

if (isvar("moses_tmp")) then

;var_o = mx2t[d=3,gx=temp[d=7],gy=temp[d=7],l=1:4@sum]/4
;var_m = temp[d=7,l=1:4@sum]/4
if (isdim(mtmax,"latitude")) then
mxobs   = linint2(tmax&longitude,tmax&latitude,tmax,True,mtmax&longitude,mtmax&latitude,0)
else
mxobs   = linint2(tmax&longitude,tmax&latitude,tmax,True,mtmax&longitude,mtmax&lat,0)
end if
mxvar_o = dim_avg_n_Wrap(mxobs,0)
mxvar_m = dim_avg_n_Wrap(mtmax(:,0,:,:),0)
m_tmax = taylor_stats(wgt,mxvar_o,mxvar_m)
ratio(1,2) = m_tmax(0)
cc(1,2)    = m_tmax(1)

; TMIN-------------------------------
;var_o = mn2t[d=3,gx=temp_1[d=7],gy=temp_1[d=7],l=1:4@sum]/4
;var_m = temp_1[d=7,l=1:4@sum]/4
if (isdim(mtmin,"latitude")) then
mnobs   = linint2(tmin&longitude,tmin&latitude,tmin,True,mtmin&longitude,mtmin&latitude,0)
else
mnobs   = linint2(tmin&longitude,tmin&latitude,tmin,True,mtmin&longitude,mtmin&lat,0)
end if
mnvar_o = dim_avg_n_Wrap(mnobs,0)
mnvar_m = dim_avg_n_Wrap(mtmin(:,0,:,:),0)
m_tmin = taylor_stats(wgt,mnvar_o,mnvar_m)
ratio(1,3) = m_tmin(0)
cc(1,3)    = m_tmin(1)

else

ratio(1,2:3) = -9999
cc(1,2:3)    = -9999

end if

; SW RAD-----------------------------
;var_o = RSS_ANN[d=8,gx=SOLAR[d=5],gy=SOLAR[d=5],l=1];:4@ave]
;var_m = SOLAR[d=5,l=1:4@sum]/4
if (isdim(msw,"latitude")) then
msobs   = linint2(sw&lon,sw&lat,sw,True,lons,msw&latitude,0)
else
msobs   = linint2(sw&lon,sw&lat,sw,True,lons,msw&lat,0)
end if
msvar_o = dim_avg_n_Wrap(msobs,0)
if (dimsizes(dimsizes(msw)).eq.4) then
msvar_m = dim_avg_n_Wrap(msw(:,0,:,:),0)
else
msvar_m = dim_avg_n_Wrap(msw(:,:,:),0)
end if
m_sw = taylor_stats(wgt,msvar_o,msvar_m)
ratio(1,4) = m_sw(0)
cc(1,4)    = m_sw(1)

; LW RAD----------------------------
;var_o = RLS_ANN[d=9,gx=LONGWAVE[d=5],gy=LONGWAVE[d=5],l=1];:4@ave]
;var_m = LONGWAVE[d=5,l=1:4@sum]/4
if (isdim(mlw,"latitude")) then
mlobs   = linint2(lw&lon,lw&lat,lw,True,lonl,mlw&latitude,0)
else
mlobs   = linint2(lw&lon,lw&lat,lw,True,lonl,mlw&lat,0)
end if
mlvar_o = dim_avg_n_Wrap(mlobs,0)
if (dimsizes(dimsizes(mlw)).eq.4) then
mlvar_m = dim_avg_n_Wrap(mlw(:,0,:,:),0)
else
mlvar_m = dim_avg_n_Wrap(mlw(:,:,:),0)
end if
m_lw = taylor_stats(wgt,mlvar_o,mlvar_m)
ratio(1,5) = m_lw(0)
cc(1,5)    = m_lw(1)

; SENS HEAT-------------------------
;var_o = HFSS_ANN[d=11,gx=SH[d=5],gy=SH[d=5],l=1];:4@ave]
;var_m = SH[d=5,l=1:4@sum]/4
if (isdim(msh,"latitude")) then
mshobs   = linint2(sh&lon,sh&lat,sh,True,lonh,msh&latitude,0)
else
mshobs   = linint2(sh&lon,sh&lat,sh,True,lonh,msh&lat,0)
end if
mshvar_o = dim_avg_n_Wrap(mshobs,0)
if (dimsizes(dimsizes(msh)).eq.4) then
mshvar_m = dim_avg_n_Wrap(msh(:,0,:,:),0)
else
mshvar_m = dim_avg_n_Wrap(msh(:,:,:),0)
end if
m_sh = taylor_stats(wgt,mshvar_o,mshvar_m)
ratio(1,6) = m_sh(0)
cc(1,6)    = m_sh(1)

; LATENT HEAT-----------------------
;var_o = HFLS_ANN[d=10,gx=LH[d=5],gy=LH[d=5],l=1];:4@ave]
;var_m = LH[d=5,l=1:4@sum]/4
if (isdim(mlh,"latitude")) then
mlhobs   = linint2(lh&lon,lh&lat,lh,True,lont,mlh&latitude,0)
else
mlhobs   = linint2(lh&lon,lh&lat,lh,True,lont,mlh&lat,0)
end if
mlhvar_o = dim_avg_n_Wrap(mlhobs,0)
if (dimsizes(dimsizes(mlh)).eq.4) then
mlhvar_m = dim_avg_n_Wrap(mlh(:,0,:,:),0)
else
mlhvar_m = dim_avg_n_Wrap(mlh(:,:,:),0)
end if
m_lh = taylor_stats(wgt,mlhvar_o,mlhvar_m)
ratio(1,7) = m_lh(0)
cc(1,7)    = m_lh(1)

;print((/ratio,cc/))

ratio1        = new((/4,8/),double)
cc1           = new((/4,8/),double)
ratio1        = -9999
cc1           = -9999
ratio1(0:1,:) = ratio
cc1(0:1,:)    = cc
ratio1(2,:)   = (/1.0531,1.1340,1.0097,1.0168,1.0492,1.36,0.9568,1.1232/)     ; CABLE
cc1(2,:)      = (/0.9965,0.8804,0.9761,0.9679,0.9845,0.9113,0.8783,0.9799/)
ratio1(3,:)   = (/1.0724,1.133,1.029,1.0371,1.046,1.2022,1.0199,1.133/)       ; MOSES
cc1(3,:)      = (/.9966,0.8731,0.9761,.9697,0.9852,.9051,0.8932,.9761/)

;PLOT-------------------------------------------------------------------------------------------------------

; "Normalized Standard Deviation (mod/obs)"

case     = (/"Cable","Moses"/)
nCase    = dimsizes(case)
varNames = (/"Tscrn","Precip","Tmax","Tmin","SW","LW","SH","LH"/)
nVars    = dimsizes(varNames)

wks      = gsn_open_wks("ps","Taylor_Diagram_"+numYr+"yrs") ; Taylor_diagram_`YR`yrs.ps

res                     = True
res@Markers             = (/14,9/)
res@Colors              = (/"Red","Blue"/)
res@varLabels           = varNames
res@varLabelsYloc       = 0.55
res@caseLabels          = case
res@stnRad              = (/0.5,1.5/)
res@ccRays              = (/0.8,0.9,0.95,0.99/)
res@centerDiffRMS       = True
res@ccRays_color        = "LightGray"
res@centerDiffRMS_color = "LightGray"
res@gsnMaximize         = True

print("")
print("Plotting Taylor Diagram")
print("")

;plot                    = taylor_diagram(wks,ratio1,cc1,res)
plot                    = taylor_diagram(wks,ratio,cc,res)

end

