; Lauren Stevens 8 Nov 2012 ==================================
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;usr   = systemfunc("echo $USERID")
; loadscript("~"+usr+"/umplot/plot_ind.ncl")
; ============================================================

begin

 usr   = systemfunc("echo $USERID")
 rid   = systemfunc("echo $RUNID")
 dir   = systemfunc("echo $DIR")
 numyr = systemfunc("echo $YR")
 resol = systemfunc("echo $RES")
 vMask = systemfunc("echo $MASK")
 Pdai  = systemfunc("echo $Pdaily")
 Ptem  = systemfunc("echo $Ptemps")

 nfil = stringtoint(numyr)*4
 ;print(nfil)

 DatDir = "~"+usr+"/umplot/data/"

 ; read in files----------

 fila = systemfunc("ls "+dir+"/"+rid+"a."+Pdai+"?????.nc | head -"+nfil)
 fa   = addfiles(fila,"r")
 ListSetType(fa,"cat")    ; "join" adds extra dim (5,12,nlat,nlon)

 filb = systemfunc("ls "+dir+"/"+rid+"a."+Ptem+"?????.nc | head -"+nfil)
 fb   = addfiles(filb,"r")
 ListSetType(fb,"cat")    ; while "cat" extends time dim (60,nlat,nlon)

 if (resol.eq.48) then
  lsm_file = addfile(DatDir+"LSM_UMn"+resol+".nc","r")
 else
  if (vMask.eq.2) then
   lsm_file = addfile(DatDir+"mask_accessV2_n"+resol+".nc","r")
  else
   lsm_file = addfile(DatDir+"mask_accessV1_n"+resol+".nc","r")
  end if
 end if

 ; read in vars----------
 tsc  = fa[:]->temp_2
 ;if (.not.isfilevar(fa[0],"temp_2")) then
 ;tas  = fa[:]->temp_1
 ;else
 ;tas  = fa[:]->temp_2
 tas  = fa[:]->temp_34
 ;end if
 smc  = fa[:]->sm
 if (isfilevar(fa[0],"soiltemp")) then
 tsl  = fa[:]->soiltemp
 else
 tsl = new(dimsizes(smc),typeof(smc))
 tsl = -9999
 copy_VarMeta(smc,tsl)
 end if
 ;if (isfilevar(fa[0],"lh")) then
 lh   = fa[:]->lh
 sh   = fa[:]->sh
 ;else
 ;lh   = fb[:]->lh
 ;sh   = fb[:]->sh
 ;end if
 ;if (isfilevar(fa[0],"precip")) then
 pr   = fa[:]->precip
 ;else
 ;pr   = new(dimsizes(lh),typeof(lh))
 ;pr   = 0.
 ;copy_VarMeta(lh,pr)
 ;end if
 if (isfilevar(fa[0],"field1526")) then
 evp1 = fa[:]->field1526
 evp2 = fa[:]->field1527
 else
 ;if (isfilevar(fa[0],"field1526")) then
 evp1 = fa[:]->field1526
 evp2 = fa[:]->field1527
 ;else
 ;evp1 = new(dimsizes(lh),typeof(lh))
 ;evp2 = new(dimsizes(lh),typeof(lh))
 ;evp1 = 0.
 ;evp2 = 0.
 ;copy_VarMeta(lh,evp1)
 ;copy_VarMeta(lh,evp2)
 ;end if
 end if
 ;printVarSummary(pr)
 
 lsm1 = lsm_file->lsm
 ;printVarSummary(lsm)

 ; working space----------
 evp       = (evp1+evp2)*24*60*60
 ;evp       = (evp1+evp2)*365*24*60*60
 copy_VarMeta(evp1,evp)
 evp@units = "mm/yr"

 pr        = pr*24*60*60
 pr@units  = "mm/day"

 dims = dimsizes(smc)
 ntim = dims(0)
 nlev = dims(1)
 nlat = dims(2)
 nlon = dims(3)
 ;print((/ntim,nlev,nlat,nlon/))

; lsm  = new((/ntim,nlev,nlat,nlon/),float)
 lsm  = new((/ntim,1,nlat,nlon/),float)
; do j = 0,nlev-1,1
 do i = 0,ntim-1,1
;  lsm(i,j,:,:) = lsm1(0,0,:,:)
  lsm(i,0,:,:) = lsm1(0,0,:,:)
 end do
; end do
 ;copy_VarMeta()

 ; required vars & constants----------
 rad          = 4.0*atan(1.0) / 180. 
 lat          = pr&latitude
 lat!0        = "latitude"
 lat&latitude = pr&latitude
 ;print(lat)

 tim          = pr&t
 tim!0        = "time"
 tim&time     = pr&t
 ;print(pr&t)

 ; mask out ocean----------
; msmc = new((/ntim,nlev,nlat,nlon/),float)
; mtsl = new((/ntim,nlev,nlat,nlon/),float)
 mpr  = where(lsm(:,0,:,:).gt.0,pr(:,0,:,:) ,pr@_FillValue)
 mlh  = where(lsm(:,0,:,:).gt.0,lh(:,0,:,:) ,lh@_FillValue)
 msh  = where(lsm(:,0,:,:).gt.0,sh(:,0,:,:) ,sh@_FillValue)
 mtas = where(lsm(:,0,:,:).gt.0,tas(:,0,:,:),tas@_FillValue)
 mtsc = where(lsm(:,0,:,:).gt.0,tsc(:,0,:,:),tsc@_FillValue)
 mevp = where(lsm(:,0,:,:).gt.0,evp(:,0,:,:),evp@_FillValue)
 msmc = where(lsm(:,0,:,:).gt.0,smc(:,0,:,:),smc@_FillValue)
 mtsl = where(lsm(:,0,:,:).gt.0,tsl(:,0,:,:),tsl@_FillValue)
; do j = 0,nlev-1,1
;  print(j)
;  msmc(:,j,:,:) = where(lsm(:,j,:,:).gt.0,smc(:,j,:,:),smc@_FillValue)
;  mtsl(:,j,:,:) = where(lsm(:,j,:,:).gt.0,tsl(:,j,:,:),tsl@_FillValue)
; end do
; copy_VarMeta(smc,msmc)
; copy_VarMeta(tsl,mtsl)
; printVarSummary(msmc)
; printVarSummary(lh)
; printVarSummary(mlh)
; copy_VarMeta(lh(:,0,:,:),mlh)

 ; weighted area-average----------
 wgt_pr  = wgt_areaave(mpr ,cos(lat*rad),1.0,0)
 wgt_tas = wgt_areaave(mtas,cos(lat*rad),1.0,0)
 wgt_tsc = wgt_areaave(mtsc,cos(lat*rad),1.0,0)
 wgt_sm1 = wgt_areaave(msmc,cos(lat*rad),1.0,0)
 wgt_ts1 = wgt_areaave(mtsl,cos(lat*rad),1.0,0)
 wgt_lh  = wgt_areaave(mlh ,cos(lat*rad),1.0,0)
 wgt_sh  = wgt_areaave(msh ,cos(lat*rad),1.0,0)
 wgt_evp = wgt_areaave(mevp,cos(lat*rad),1.0,0)
 ;printVarSummary(wgt_ts1)
 ;printVarSummary(wgt_sm1)
 ;printVarSummary(wgt_tas)
 ;printVarSummary(smc)
 ;printMinMax(wgt_lh,True)
 ;printMinMax(mlh,True)

 wgt_pr!0       = "t"
; wgt_pr!1       = "surface"
 wgt_pr&t       = pr&t
; wgt_pr&surface = pr&surface
 copy_VarMeta(wgt_pr,wgt_tas)
 copy_VarMeta(wgt_pr,wgt_tsc)
 copy_VarMeta(wgt_pr,wgt_lh)
 copy_VarMeta(wgt_pr,wgt_sh)
 copy_VarMeta(wgt_pr,wgt_evp)
 wgt_sm1!0       = "t"
; wgt_sm1!1       = "surface"
 wgt_sm1&t       = smc&t 
; wgt_sm1&surface = smc&level6
 copy_VarMeta(wgt_sm1,wgt_ts1)

; wgt_smc = wgt_sm1(surface|:,t|:)
; wgt_tsl = wgt_ts1(surface|:,t|:)
 wgt_smc = wgt_sm1
 wgt_tsl = wgt_ts1

 ; detrended-----------
 dt_wgt_pr  = dtrend_quadratic(wgt_pr ,0) ;True,True,0)
 dt_wgt_evp = dtrend_quadratic(wgt_evp,0) ;True,True,0)
 dt_wgt_tsl = dtrend_quadratic(wgt_tsl,0) ;True,True,0)
 dt_wgt_smc = dtrend_quadratic(wgt_smc,0) ;True,True,0)
 dt_wgt_lh  = dtrend_quadratic(wgt_lh ,0) ;True,True,0)
 dt_wgt_sh  = dtrend_quadratic(wgt_sh ,0) ;True,True,0)
 dt_wgt_tas = dtrend_quadratic(wgt_tas,0) ;True,True,0)
 dt_wgt_tsc = dtrend_quadratic(wgt_tsc,0) ;True,True,0)

 ;dt_wgt_pr  = rmAnnCycle1D(wgt_pr)
 ;dt_wgt_evp = rmAnnCycle1D(wgt_evp)
 ;dt_wgt_tsl = rmAnnCycle1D(wgt_tsl)
 ;dt_wgt_smc = rmAnnCycle1D(wgt_smc)
 ;dt_wgt_lh  = rmAnnCycle1D(wgt_lh)
 ;dt_wgt_sh  = rmAnnCycle1D(wgt_sh)
 ;dt_wgt_tas = rmAnnCycle1D(wgt_tas)
 ;dt_wgt_tsc = rmAnnCycle1D(wgt_tsc)
 
 ; plot----------
; levlh  = (/-20,0,20,40,60,80,100,120,140,160,180/)
; IndPlot("./","lath","5","ps","BlAqGrYeOrReVi200",mlh(0,:,:),False,False,levlh,False,False)

 print("Plotting Globally Averaged Daily Mean Timeseries")

 wks1 = gsn_open_wks("ps","daily_means_ts1_"+numyr+"yrs")
 wks2 = gsn_open_wks("ps","daily_means_ts2_"+numyr+"yrs")
 wks3 = gsn_open_wks("ps","daily_means_dtts1_"+numyr+"yrs")
 wks4 = gsn_open_wks("ps","daily_means_dtts2_"+numyr+"yrs")
 res  = True
 res@gsnFrame = False
 res@gsnDraw  = False
 res@txFontHeightF = 0.016
 res@tiXAxisString = "days"

 plot_1 = new(4,graphic)
 plot_2 = new(4,graphic)
 plot_3 = new(4,graphic)
 plot_4 = new(4,graphic)
 
; plot_1(0) = gsn_csm_xy(wks1,tim,wgt_sh(:,0) ,res) 
; plot_1(1) = gsn_csm_xy(wks1,tim,wgt_tas(:,0),res) 
; plot_1(2) = gsn_csm_xy(wks1,tim,wgt_lh(:,0) ,res) 
; plot_1(3) = gsn_csm_xy(wks1,tim,wgt_tsc(:,0),res)
; plot_2(0) = gsn_csm_xy(wks2,tim,wgt_smc(0,:),res) ; top layer
; ;plot_2(0) = gsn_csm_xy(wks2,tim,wgt_smc(:,:),res) ; all levels
; plot_2(1) = gsn_csm_xy(wks2,tim,wgt_pr(:,0) ,res)
; plot_2(2) = gsn_csm_xy(wks2,tim,wgt_tsl(0,:),res) ; top layer
; ;plot_2(2) = gsn_csm_xy(wks2,tim,wgt_tsl(:,:),res) ; all levels
; plot_2(3) = gsn_csm_xy(wks2,tim,wgt_evp(:,0),res) 

 res@tiYAxisString = "Wm-2"
 plot_1(0) = gsn_csm_xy(wks1,tim,wgt_sh(:) ,res)
 res@tiYAxisString = "K"
 plot_1(1) = gsn_csm_xy(wks1,tim,wgt_tas(:),res)
 res@tiYAxisString = "Wm-2"
 plot_1(2) = gsn_csm_xy(wks1,tim,wgt_lh(:) ,res)
 res@tiYAxisString = "K"
 plot_1(3) = gsn_csm_xy(wks1,tim,wgt_tsc(:),res)
 res@tiYAxisString = "kgm-2"
 plot_2(0) = gsn_csm_xy(wks2,tim,wgt_smc(:),res) ; top layer
 ;plot_2(0) = gsn_csm_xy(wks2,tim,wgt_smc(:,:),res) ; all levels
 res@tiYAxisString = "mm/day"
 plot_2(1) = gsn_csm_xy(wks2,tim,wgt_pr(:) ,res)
 res@tiYAxisString = "K"
 plot_2(2) = gsn_csm_xy(wks2,tim,wgt_tsl(:),res) ; top layer
 ;plot_2(2) = gsn_csm_xy(wks2,tim,wgt_tsl(:,:),res) ; all levels
 res@tiYAxisString = "mm/day"
 plot_2(3) = gsn_csm_xy(wks2,tim,wgt_evp(:),res)

 res@tiYAxisString = "Wm-2"
 plot_3(0) = gsn_csm_xy(wks3,tim,dt_wgt_sh(:) ,res)
 res@tiYAxisString = "K"
 plot_3(1) = gsn_csm_xy(wks3,tim,dt_wgt_tas(:),res)
 res@tiYAxisString = "Wm-2"
 plot_3(2) = gsn_csm_xy(wks3,tim,dt_wgt_lh(:) ,res)
 res@tiYAxisString = "K"
 plot_3(3) = gsn_csm_xy(wks3,tim,dt_wgt_tsc(:),res)
 res@tiYAxisString = "kgm-2"
 plot_4(0) = gsn_csm_xy(wks4,tim,dt_wgt_smc(:),res) ; top layer
 ;plot_4(0) = gsn_csm_xy(wks4,tim,dt_wgt_smc(:,:),res) ; all levels
 res@tiYAxisString = "mm/day"
 plot_4(1) = gsn_csm_xy(wks4,tim,dt_wgt_pr(:) ,res)
 res@tiYAxisString = "K"
 plot_4(2) = gsn_csm_xy(wks4,tim,dt_wgt_tsl(:),res) ; top layer
 ;plot_4(2) = gsn_csm_xy(wks4,tim,dt_wgt_tsl(:,:),res) ; all levels
 res@tiYAxisString = "mm/day"
 plot_4(3) = gsn_csm_xy(wks4,tim,dt_wgt_evp(:),res)
 
 ; panel plot----------
 pres                       = True
; pres@gsnMaximize           = True
; pres@gsnOrientation        = "landscape"
; pres@gsnPanelRowSpec       = True
 ;pres@gsnPanelLabelBar      = True
 ;pres@pmLabelBarWidthF      = 0.8         ; makes label bar bigger
 ;pres@txString              = (field@long_name + " (" + field@units + ")")

 pres2 = pres
 pres@gsnPanelFigureStrings  = (/"sh","tas","lh","tsc"/)
 pres2@gsnPanelFigureStrings = (/"smc","pr","tsl","evp"/)
 ;pres@gsnPanelFigureStrings  = (/"sh-W/m2","tas-K","lh-W/m2","tsc-K"/)
 ;pres2@gsnPanelFigureStrings = (/"smc-kg/m2","pr-mm/dy","tsl-K","evp-mm/dy"/)

 gsn_panel(wks1,plot_1,(/2,2/),pres)
 gsn_panel(wks2,plot_2,(/2,2/),pres2)
 gsn_panel(wks3,plot_3,(/2,2/),pres)
 gsn_panel(wks4,plot_4,(/2,2/),pres2)

end

