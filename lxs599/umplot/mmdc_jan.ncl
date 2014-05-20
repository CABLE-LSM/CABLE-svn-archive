;==============================================================================
; Lauren Stevens 19 July 2011
;==============================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;usr = systemfunc("echo $USERID")
;loadscript("~"+usr+"/umplot/plot_mmdc.ncl")
;==============================================================================

;procedure mmdc_flux(LongName:string,FlxName:string,mosDir:string,cabDir:string)

begin

usr      = systemfunc("echo $USERID")
numYr    = systemfunc("echo $YR")
modl     = systemfunc("echo $MODEL")
filDir   = systemfunc("echo $PWD")
;filDir   = systemfunc("echo $DIRW")
jpg      = systemfunc("echo $CHFMT")
if(jpg.eq."y") then
 jpeg = True
else
 jpeg = False
end if

tsDir    = "~"+usr+"/umplot/obs/"
home     = "~"+usr+"/"
flxobs   = "/home/cmar/ste69f/obsdata/fluxnet/PALS/flxnt_mon/"

;hayfil = addfile( tsDir + "hayobsjuly.nc"         ,"r")
manfil = addfile( tsDir + "mmdc_ma96_roll.nc"     ,"r")
tbfil  = addfile( tsDir + "mmdc_tumbarumba_fv.nc" ,"r")
vi97   = addfile( tsDir + "mmdc_vi97_roll.nc"     ,"r")
hy_new = addfile(flxobs + "mmdc_ts_Hyytiala.nc"   ,"r")
lo_new = addfile(flxobs + "mmdc_ts_Loobos.nc"     ,"r")
tb_new = addfile(flxobs + "mmdc_ts_Tumba.nc"      ,"r")
nb_new = addfile(flxobs + "mmdc_ts_Boreas.nc"     ,"r")
bv_new = addfile(flxobs + "mmdc_ts_Bondv.nc"      ,"r")
th_new = addfile(flxobs + "mmdc_ts_Tharandt.nc"   ,"r")
hv_new = addfile(flxobs + "mmdc_ts_Harvard.nc"    ,"r")

;==============================================================================

fils   = systemfunc("ls "+tsDir+"mmdc_??.nc")
; 1:BV 2:HV 3:HY 4:LO 5:LW 6:NB 7:TH 8:VI 9:WB
f      = addfiles(fils,"r")
ListSetType(f,"join")
sens   = f[:]->H
lat    = f[:]->LH
radnet = f[:]->Rn
;temp2m = f[:]->Ta
;preci  = f[:]->PPT

; sites,mons,tstp; 9,1,1,12,48 ?
;sens2   = sens(:,0,0,:,:)   ;dim_avg_n_Wrap(sens,1)
;lat2    = lat(:,0,0,:,:)    ;dim_avg_n_Wrap(lat,1)
;radnet2 = radnet(:,0,0,:,:) ;dim_avg_n_Wrap(radnet,1)
sens2   = dim_avg_n_Wrap(sens,1)
lat2    = dim_avg_n_Wrap(lat,1)
radnet2 = dim_avg_n_Wrap(radnet,1)

; 0 :
; 18: afN afS bv dw dh hp hv hay hy in litw loob ma nb th tb vi wb
;a31: afN afr afS alk ant1 ant2 ant ant3 bv dw dh grnl hp hv hay hy ind lw lo ma nb rus1 rus2 sd1 sd2 sd3 sd4 th tb vi wb
;f31: man,hay,dwaters,bats,india,hyy,afN,afS,tumb,thar,bond,ding,lwas,walk,loo,harv,bore,viel
;    ,ant1,ant2,ant3,rus1,alk,afr,grl,ant,sau1,sau2,sau3,rus2,sau4
;
fils2  = systemfunc("ls "+filDir+"/mmdc_*_roll_"+numYr+"yrs.nc")
f2     = addfiles(fils2,"r")
ListSetType(f2,"join")
if (isfilevar(f2[0],"hfss")) then
 msh  = f2[:]->hfss      ;sh [ncl_join | 18] x [axis_0 | 12] x [axis_1 | 48]
 mlh  = f2[:]->hfls      ;lh
 mrad = f2[:]->field3333 ;field202
else
 msh  = f2[:]->sh
 mlh  = f2[:]->lh
 mrad = f2[:]->field202
end if

if (isfilevar(f2[0],"solar")) then
 mrad_1= f2[:]->solar
 mrad_2= f2[:]->longwave

 mrad_3 = mrad_1+mrad_2
 copy_VarMeta(mrad_1,mrad_3)

 mrad_3&axis_1 = ispan(0,47,6)

 dName2 = getvardims(mlh)
 ;dName2 = getvardims(mlh1)
 ; 12,8,31
 mrad_4 = ftcurvp(mrad_3&axis_1,mrad_3,48,mlh&$dName2(1)$)
 ;mrad_4 = linint1_n_Wrap(mrad_3&axis_1,mrad_3,True,mlh&$dName2(1)$,0,1)
 ;;mrad_4 = linint1_n_Wrap(mrad_3&axis_1,mrad_3,True,mlh1&$dName2(1)$,0,1)
 mrad1 = new((/12,48,31/),double)
 copy_VarMeta(msh,mrad1)
 ;copy_VarMeta(msh1,mrad1)

 mrad1 = mrad_4

if (dimsizes(fils2).eq.31) then
; 31: afN afr afS alk ant1 ant2 ant ant3 bv dw dh grnl hp hv hay hy ind lw lo ma nb rus1 rus2 sd1 sd2 sd3 sd4 th tb vi wb
 mrad1(:,42:47,8)   = mrad_4(:,0:5 ,8)   ;bv - black
 mrad1(:,0:41 ,8)   = mrad_4(:,6:47,8)   ;bv - black
 mrad1(:,44:47,13)  = mrad_4(:,0:3 ,13)  ;hv - blue
 mrad1(:,0:43 ,13)  = mrad_4(:,4:47,13)  ;hv - blue
 mrad1(:,0:6,15)    = mrad_4(:,41:47,15) ;hy - red
 mrad1(:,7:47 ,15)  = mrad_4(:,0:40 ,15) ;hy - red
 mrad1(:,0:4 ,18)   = mrad_4(:,43:47,18) ;lo - green
 mrad1(:,5:47,18)   = mrad_4(:,0:42 ,18) ;lo - green
 mrad1(:,41:47 ,17) = mrad_4(:,0:6,17)   ;lw - pale pink
 mrad1(:,0:40,17)   = mrad_4(:,7:47,17)  ;lw - pale pink
 mrad1(:,40:47 ,20) = mrad_4(:,0:7,20)   ;nb - light blue
 mrad1(:,0:39,20)   = mrad_4(:,8:47 ,20) ;nb - light blue
 mrad1(:,0:4 ,27)   = mrad_4(:,43:47,27) ;th - hot pink
 mrad1(:,5:47,27)   = mrad_4(:,0:42 ,27) ;th - hot pink
 ;vi 29
 mrad1(:,43:47,30)  = mrad_4(:,0:4 ,30)  ;wb - orange
 mrad1(:,0:42 ,30)  = mrad_4(:,5:47 ,30) ;wb - orange
 ;ma 19
 ;tb 28
else
; 18: afN afS bv dw dh hp hv hay hy in litw loob ma nb th tb vi wb
 mrad1(:,42:47,2)   = mrad_4(:,0:5 ,2)   ;bv - black
 mrad1(:,0:41 ,2)   = mrad_4(:,6:47,2)   ;bv - black
 mrad1(:,44:47,6)   = mrad_4(:,0:3 ,6)   ;hv - blue
 mrad1(:,0:43 ,6)   = mrad_4(:,4:47,6)   ;hv - blue
 mrad1(:,0:6,8)     = mrad_4(:,41:47,8)  ;hy - red
 mrad1(:,7:47 ,8)   = mrad_4(:,0:40 ,8)  ;hy - red
 mrad1(:,0:4 ,11)   = mrad_4(:,43:47,11) ;lo - green
 mrad1(:,5:47,11)   = mrad_4(:,0:42 ,11) ;lo - green
 mrad1(:,41:47 ,10) = mrad_4(:,0:6,10)   ;lw - pale pink
 mrad1(:,0:40,10)   = mrad_4(:,7:47,10)  ;lw - pale pink
 mrad1(:,40:47 ,13) = mrad_4(:,0:7,13)   ;nb - light blue
 mrad1(:,0:39,13)   = mrad_4(:,8:47 ,13) ;nb - light blue
 mrad1(:,0:4 ,14)   = mrad_4(:,43:47,14) ;th - hot pink
 mrad1(:,5:47,14)   = mrad_4(:,0:42 ,14) ;th - hot pink
 ;vi 16
 mrad1(:,43:47,17)  = mrad_4(:,0:4 ,17)  ;wb - orange
 mrad1(:,0:42 ,17)  = mrad_4(:,5:47 ,17) ;wb - orange
 ;ma 12
 ;tb 15
end if

; lest 7/1/14 order due to order in stash/file setup not order of file names.
; mrad1(:,0:6,5)     = mrad_4(:,41:47,5)  ;hy - red
; mrad1(:,7:47 ,5)   = mrad_4(:,0:40 ,5)  ;hy - red
; mrad1(:,0:4 ,9)    = mrad_4(:,43:47,9)  ;th - hot pink
; mrad1(:,5:47,9)    = mrad_4(:,0:42 ,9)  ;th - hot pink
; mrad1(:,42:47,10)  = mrad_4(:,0:5 ,10)  ;bv - black
; mrad1(:,0:41 ,10)  = mrad_4(:,6:47,10)  ;bv - black
; mrad1(:,41:47 ,12) = mrad_4(:,0:6,12)   ;lw - pale pink 
; mrad1(:,0:40,12)   = mrad_4(:,7:47,12)  ;lw - pale pink 
; mrad1(:,43:47,13)  = mrad_4(:,0:4 ,13)  ;wb - orange    
; mrad1(:,0:42 ,13)  = mrad_4(:,5:47 ,13) ;wb - orange    
; mrad1(:,0:4 ,14)   = mrad_4(:,43:47,14) ;lo - green     
; mrad1(:,5:47,14)   = mrad_4(:,0:42 ,14) ;lo - green     
; mrad1(:,44:47,15)  = mrad_4(:,0:3 ,15)  ;hv - blue
; mrad1(:,0:43 ,15)  = mrad_4(:,4:47,15)  ;hv - blue
; mrad1(:,40:47 ,16) = mrad_4(:,0:7,16)   ;nb - light blue
; mrad1(:,0:39,16)   = mrad_4(:,8:47 ,16) ;nb - light blue
end if

if (isvar("mrad1") .and. sum(mrad).eq.0.) then ; if (modl.eq."m") then
 mrad = mrad1
end if

; obs ---
sens_ma   = manfil ->H
lat_ma    = manfil ->LH
rad_ma    = manfil ->Rn
;temp_ma   = manfil ->Ta
;prec_ma   = manfil ->PPT
sens_tb   = tbfil  ->Qh     ;H
lat_tb    = tbfil  ->Qle    ;LH
rad_tb    = tbfil  ->Rnet   ;Rn
;sens_tb = tb_new ->Qh    ;H
;lat_tb  = tb_new ->Qle   ;LH
;rad_tb  = tb_new ->Rnet  ;Rn
;;temp_tb   = tbfil  ->Tair   ;Ta
;;prec_tb   = tbfil  ->Rainf  ;PPT

rad_vi    = vi97   ->Rn
sens_vi   = vi97   ->H
;sens_hay  = hayfil ->h
;lat_hay   = hayfil ->lh

sens_hy = hy_new->Qh
lat_hy  = hy_new->Qle
rad_hy  = hy_new->Rnet
sens_hy2= sens_hy(:,:,0)
lat_hy2 = lat_hy(:,:,0)
rad_hy2 = rad_hy(:,:,0)

sens_lo = lo_new->Qh
lat_lo  = lo_new->Qle
rad_lo  = lo_new->Rnet
sens_lo2= sens_lo(:,:,0)
lat_lo2 = lat_lo(:,:,0)
rad_lo2 = rad_lo(:,:,0)

sens_nb = nb_new->Qh
lat_nb  = nb_new->Qle
rad_nb  = nb_new->Rnet
sens_nb2= sens_nb(:,:,0)
lat_nb2 = lat_nb(:,:,0)
rad_nb2 = rad_nb(:,:,0)

;sens_th = th_new->Qh
;lat_th  = th_new->Qle
;rad_th  = th_new->Rnet ; missing
;sens_th2= sens_th(:,:,0)
;lat_th2 = lat_th(:,:,0)
;rad_th2 = rad_th(:,:,0)
;
;sens_hv = hv_new->Qh
;lat_hv  = hv_new->Qle
;rad_hv  = hv_new->Rnet
;sens_hv2= sens_hv(:,:,0)
;lat_hv2 = lat_hv(:,:,0)
;rad_hv2 = rad_hv(:,:,0)
;
;sens_bv = bv_new->Qh
;lat_bv  = bv_new->Qle
;rad_bv  = bv_new->Rnet
;sens_bv2= sens_bv(:,:,0)
;lat_bv2 = lat_bv(:,:,0)
;rad_bv2 = rad_bv(:,:,0)

;printVarSummary()
sens_tb2 = dim_avg_n_Wrap(sens_tb,0)
lat_tb2  = dim_avg_n_Wrap(lat_tb,0)
rad_tb2  = dim_avg_n_Wrap(rad_tb,0)
;sens_tb2 = sens_tb(:,:,0) ;dim_avg_n_Wrap(sens_tb,0)
;lat_tb2  = lat_tb(:,:,0)  ;dim_avg_n_Wrap(lat_tb,0)
;rad_tb2  = rad_tb(:,:,0)  ;dim_avg_n_Wrap(rad_tb,0)

;printVarSummary()
sens_ma2 = dim_avg_n_Wrap(sens_ma,0)
lat_ma2  = dim_avg_n_Wrap(lat_ma,0)
rad_ma2  = dim_avg_n_Wrap(rad_ma,0)
;sens_ma2 = sens_ma(0,:,:) ;dim_avg_n_Wrap(sens_ma,0)
;lat_ma2  = lat_ma(0,:,:)  ;dim_avg_n_Wrap(lat_ma,0)
;rad_ma2  = rad_ma(0,:,:)  ;dim_avg_n_Wrap(rad_ma,0)

sh_tb = new((/12,48/),double)
lh_tb = new((/12,48/),double)
rn_tb = new((/12,48/),double)
;sh_hay = new(48,double)
;lh_hay = new(48,double)
do i = 0,23
sh_tb(:,2*i+1)=sens_tb2(:,i)
sh_tb(:,2*i)  =sens_tb2(:,i)
lh_tb(:,2*i+1)=lat_tb2(:,i)
lh_tb(:,2*i)  =lat_tb2(:,i)
rn_tb(:,2*i+1)=rad_tb2(:,i)
rn_tb(:,2*i)  =rad_tb2(:,i)
;sh_hay(2*i+1)  =sens_hay(i)
;sh_hay(2*i)    =sens_hay(i)
;lh_hay(2*i+1)  =lat_hay(i)
;lh_hay(2*i)    =lat_hay(i)
end do

; Lest 15/1/14 - interpolate
;printVarSummary(sens_tb2)
;sens_tb2& = ispan(1,48,2)
;lat_tb2&  = ispan(1,48,2)
;rad_tb2&  = ispan(1,48,2)
;sh_tb = ftcurvp(sens_tb2&,sens_tb2,48,msh&timestep)
;lh_tb = ftcurvp( lat_tb2&, lat_tb2,48,mlh&timestep)
;rn_tb = ftcurvp( rad_tb2&, rad_tb2,48,mrnt&timestep)
;sens_hay& = ispan(1,48,2) ; 0,47,2
;lat_hay&  = ispan(1,48,2) ; 0,47,2
;sh_hay = ftcurvp(sens_hay&,sens_hay,48,msh&timestep)
;lh_hay = ftcurvp( lat_hay&, lat_hay,48,mlh&timestep)
;sh_hay!0  = "timestep"
;lh_hay!0  = "timestep"
;sh_hay&timestep = msh&timestep
;lh_hay&timestep = mlh&timestep

;obs ---

;==============================================================================

sensht            = new((/11,48/),double)
latht             = new((/11,48/),double)
rad               = new((/11,48/),double)
;temp              = new((/11,48/),double)
;precip            = new((/11,48/),double)
sensht            = -9999
latht             = -9999
rad               = -9999
;temp              = -9999
;precip            = -9999
sensht@_FillValue = -9999
latht@_FillValue  = -9999
rad@_FillValue    = -9999
;temp@_FillValue   = -9999
;precip@_FillValue = -9999

do k=0,8,1
sensht(k,:)       = sens2(k,0,0,0:47)   ;H
latht(k,:)        = lat2(k,0,0,0:47)    ;LH
rad(k,:)          = radnet2(k,0,0,0:47) ;Rn
;sensht(k,:)       = sens2(k,0,0:47)   ;H
;latht(k,:)        = lat2(k,0,0:47)    ;LH
;rad(k,:)          = radnet2(k,0,0:47) ;Rn
;;temp(k,:)         = temp2m(k,0:47)      ;Ta
;;precip(k,:)       = preci(k,0:47)       ;PPT
end do
sensht!0 = "axis_1"
latht!0  = "axis_1"
rad!0    = "axis_1"
;sh_hay!0  = "axis_2"
;lh_hay!0  = "axis_2"
;sensht(9,:)     = sh_hay(0:47)
;latht(9,:)      = lh_hay(0:47)
sensht(9,:)     = sens_ma2(0,0:47) ;H
latht(9,:)      = lat_ma2(0,0:47)  ;LH
rad(9,:)        = rad_ma2(0,0:47)  ;Rn
;temp(9,:)       = temp_ma(0:47,0)  ;Ta
;precip(9,:)     = prec_ma(0:47,0)  ;PPT
sensht(10,:)    = sh_tb(0,0:47)    ;H
latht(10,:)     = lh_tb(0,0:47)    ;LH
rad(10,:)       = rn_tb(0,0:47)    ;Rn
;temp(10,:)      = temp_tb(0:47,0)  ;Ta
;precip(10,:)    = prec_tb(0:47,0)  ;PPT

rad(2,:)    = rad_hy2(0,:);rad(2,:)*3 +2*9999
;rad(2,:)    = rad(2,:)*3 +2*9999
sensht(2,:) = sens_hy2(0,:);(sensht(2,:)*3 + 9999)/2
;sensht(2,:) = (sensht(2,:)*3 + 9999)/2
latht(2,:)  = lat_hy2(0,:);(latht(2,:)*3 + 9999)/2
;latht(2,:)  = (latht(2,:)*3 + 9999)/2

rad(3,:)    = rad_lo2(0,:);(rad(3,:)*3 + 9999)/2
;rad(3,:)    = (rad(3,:)*3 + 9999)/2
sensht(3,:) = sens_lo2(0,:);(sensht(3,:)*3 + 9999)/2
;sensht(3,:) = (sensht(3,:)*3 + 9999)/2
latht(3,:)  = lat_lo2(0,:);(latht(3,:)*3 + 9999)/2
;latht(3,:)  = (latht(3,:)*3 + 9999)/2

;sensht(6,:) = (sensht(6,:)*4 - sens_th())/3 
;latht(6,:)  = (latht(6,:)*4 - lh_th())/3 

rad(7,:)    = rad(7,:)*3 +2*9999
sensht(7,:) = (sensht(7,:)*3 + 9999 - sens_vi(0,0,:))
latht(7,:)  = (latht(7,:)*3 + 9999)/2

rad(5,:)    = (rad(5,:)*5+9999)/4
sensht(5,:) = (sensht(5,:)*5+9999)/4
latht(5,:)  = (latht(5,:)*5 + 9999)/4
;rad(5,:)    = rad_nb2(0,:);(rad(3,:)*3 + 9999)/2
;sensht(5,:) = sens_nb2(0,:);(sensht(3,:)*3 + 9999)/2
;latht(5,:)  = lat_nb2(0,:);(latht(3,:)*3 + 9999)/2

;; lest 7/1/14 - AMOJ
;sensht(9,:)     = -9999 ;sens_ma2(0,0:47) ;H
;latht(9,:)      = -9999 ;lat_ma2(0,0:47)  ;LH
;rad(9,:)        = -9999 ;rad_ma2(0,0:47)  ;Rn
;;temp(9,:)       = -9999 ;temp_ma(0:47,0)  ;Ta
;;precip(9,:)     = -9999 ;prec_ma(0:47,0)  ;PPT
;sensht(10,:)    = -9999 ;sh_tb(0,0:47)    ;H
;latht(10,:)     = -9999 ;lh_tb(0,0:47)    ;LH
;rad(10,:)       = -9999 ;rn_tb(0,0:47)    ;Rn
;;temp(10,:)      = -9999 ;temp_tb(0:47,0)  ;Ta
;;precip(10,:)    = -9999 ;prec_tb(0:47,0)  ;PPT
;rad(7,:)        = -9999 ;rad(7,:)*3 +2*9999
;sensht(7,:)     = -9999 ;(sensht(7,:)*3 + 9999 - sens_vi(0,0,:))
;latht(7,:)      = -9999 ;(latht(7,:)*3 + 9999)/2
;; lest - AMOJ

; 31: afN afr afS alk ant1 ant2 ant ant3 bv dw dh grnl hp hv hay hy ind lw lo ma nb rus1 rus2 sd1 sd2 sd3 sd4 th tb vi wb
mrad_ro = new((/11,48/),double)
;if (dimsizes(f2[0]).eq.18) then
mrad_ro(0,:)  = mrad(2,0,:)  ; bv ; black
mrad_ro(1,:)  = mrad(6,0,:)  ; hv ; blue
mrad_ro(2,:)  = mrad(8,0,:)  ; hy ; red
mrad_ro(3,:)  = mrad(11,0,:) ; lo ; green
mrad_ro(4,:)  = mrad(10,0,:) ; lw ; darkviolet
mrad_ro(5,:)  = mrad(13,0,:) ; nb ; cyan
mrad_ro(6,:)  = mrad(14,0,:) ; th ; magenta
mrad_ro(7,:)  = mrad(16,0,:) ; vi ; firebrick
mrad_ro(8,:)  = mrad(17,0,:) ; wb ; orange
mrad_ro(9,:)  = mrad(12,0,:) ; ma ; SlateBlue
mrad_ro(10,:) = mrad(15,0,:) ; tb ; navy       and hay = forestgreen
;else
;mrad_ro(0,:)  = mrad(8,0,:)  ; bv ; black
;mrad_ro(1,:)  = mrad(13,0,:)  ; hv ; blue
;mrad_ro(2,:)  = mrad(15,0,:)  ; hy ; red
;mrad_ro(3,:)  = mrad(18,0,:) ; lo ; green
;mrad_ro(4,:)  = mrad(17,0,:) ; lw ; darkviolet
;mrad_ro(5,:)  = mrad(20,0,:) ; nb ; cyan
;mrad_ro(6,:)  = mrad(27,0,:) ; th ; magenta
;mrad_ro(7,:)  = mrad(29,0,:) ; vi ; firebrick
;mrad_ro(8,:)  = mrad(30,0,:) ; wb ; orange
;mrad_ro(9,:)  = mrad(19,0,:) ; ma ; SlateBlue
;mrad_ro(10,:) = mrad(28,0,:) ; tb ; navy       and hay = forestgreen
;end if
; lest 7/1/14 - order from file
; mrad_ro(0,:)  = mrad(10,0,:) ; bv ; black
; mrad_ro(1,:)  = mrad(15,0,:) ; hv ; blue
; ;mrad_ro(1,:)  = mrad_h(0,:);mrad(15,0,:) ; hv ; blue
; mrad_ro(2,:)  = mrad(5,0,:)  ; hy ; red
; mrad_ro(3,:)  = mrad(14,0,:) ; lo ; green
; ;mrad_ro(3,:)  = mrad_l(0,:) ;mrad(14,0,:) ; lo ; green
; mrad_ro(4,:)  = mrad(12,0,:) ; lw ; darkviolet
; mrad_ro(5,:)  = mrad(16,0,:) ; nb ; cyan
; mrad_ro(6,:)  = mrad(9,0,:)  ; th ; magenta
; mrad_ro(7,:)  = mrad(17,0,:) ; vi ; firebrick
; mrad_ro(8,:)  = mrad(13,0,:) ; wb ; orange
; mrad_ro(9,:)  = mrad(0,0,:)  ; ma ; SlateBlue
; mrad_ro(10,:) = mrad(8,0,:)  ; tb ; navy       and hay = forestgreen

msh_ro = new((/11,48/),double)
;if (dimsizes(f2[0]).eq.18) then
msh_ro(0,:)  = msh(2,0,:)  ; bv
msh_ro(1,:)  = msh(6,0,:)  ; hv
msh_ro(2,:)  = msh(8,0,:)  ; hy
msh_ro(3,:)  = msh(11,0,:) ; lo
msh_ro(4,:)  = msh(10,0,:) ; lw
msh_ro(5,:)  = msh(13,0,:) ; nb
msh_ro(6,:)  = msh(14,0,:) ; th
msh_ro(7,:)  = msh(16,0,:) ; vi
msh_ro(8,:)  = msh(17,0,:) ; wb
msh_ro(9,:)  = msh(12,0,:) ; ma
msh_ro(10,:) = msh(15,0,:) ; tb
;else
;msh_ro(0,:)  = msh(8,0,:)  ; bv ; black
;msh_ro(1,:)  = msh(13,0,:)  ; hv ; blue
;msh_ro(2,:)  = msh(15,0,:)  ; hy ; red
;msh_ro(3,:)  = msh(18,0,:) ; lo ; green
;msh_ro(4,:)  = msh(17,0,:) ; lw ; darkviolet
;msh_ro(5,:)  = msh(20,0,:) ; nb ; cyan
;msh_ro(6,:)  = msh(27,0,:) ; th ; magenta
;msh_ro(7,:)  = msh(29,0,:) ; vi ; firebrick
;msh_ro(8,:)  = msh(30,0,:) ; wb ; orange
;msh_ro(9,:)  = msh(19,0,:) ; ma ; SlateBlue
;msh_ro(10,:) = msh(28,0,:) ; tb ; navy       and hay = forestgreen
;end if

mlh_ro = new((/11,48/),double)
;if (dimsizes(f2[0]).eq.18) then
mlh_ro(0,:)  = mlh(2,0,:)  ; bv
mlh_ro(1,:)  = mlh(6,0,:)  ; hv
mlh_ro(2,:)  = mlh(8,0,:)  ; hy
mlh_ro(3,:)  = mlh(11,0,:) ; lo
mlh_ro(4,:)  = mlh(10,0,:) ; lw
mlh_ro(5,:)  = mlh(13,0,:) ; nb
mlh_ro(6,:)  = mlh(14,0,:) ; th
mlh_ro(7,:)  = mlh(16,0,:) ; vi
mlh_ro(8,:)  = mlh(17,0,:) ; wb
mlh_ro(9,:)  = mlh(12,0,:) ; ma
mlh_ro(10,:) = mlh(15,0,:) ; tb
;else
;mlh_ro(0,:)  = mlh(8,0,:)  ; bv ; black
;mlh_ro(1,:)  = mlh(13,0,:)  ; hv ; blue
;mlh_ro(2,:)  = mlh(15,0,:)  ; hy ; red
;mlh_ro(3,:)  = mlh(18,0,:) ; lo ; green
;mlh_ro(4,:)  = mlh(17,0,:) ; lw ; darkviolet
;mlh_ro(5,:)  = mlh(20,0,:) ; nb ; cyan
;mlh_ro(6,:)  = mlh(27,0,:) ; th ; magenta
;mlh_ro(7,:)  = mlh(29,0,:) ; vi ; firebrick
;mlh_ro(8,:)  = mlh(30,0,:) ; wb ; orange
;mlh_ro(9,:)  = mlh(19,0,:) ; ma ; SlateBlue
;mlh_ro(10,:) = mlh(28,0,:) ; tb ; navy       and hay = forestgreen
;end if

;==============================================================================

sh_ANN            = new((/11,48/),double)
sh_ANN            = msh_ro
sh_ANN            = sh_ANN - sensht
lh_ANN            = new((/11,48/),double)
lh_ANN            = mlh_ro
lh_ANN            = lh_ANN - latht
rad_ANN           = new((/11,48/),double)
rad_ANN           = mrad_ro
rad_ANN           = rad_ANN - rad 

;Plot-------------------------------------------------------------------------

print("")
print("Plotting Jan MMDC Biases for all Sites")

y = fspan(0,23.5,48)
;y = ispan(1,12,1)

pltName               = "mmdc_jan_"+numYr+"yrs"
wks                   = gsn_open_wks("ps",pltName)
res                   = True
res@gsnFrame          = False
res@gsnDraw           = False
res@xyLineColors      = (/"black","deepskyblue1","red","green","mediumorchid1","cyan","magenta","firebrick","orange","SlateBlue","navy","forestgreen"/)
res@xyDashPatterns    = (/0,0,0,0,0,0,0,0,0,0,1,1/)
res@xyLineThicknesses = (/2,2,2,2,2,2,2,2,2,2,2,2/)
res@gsnYRefLine       = 0.

;-----------------------------------------------------------------------------

;res@gsnMaximize    = True
plot               = new(3,graphic)
res@tiYAxisFontHeightF = 0.02
res@tiYAxisOffsetXF    = 0.01
res@tiYAxisString      = "W/m2"
res@tiXAxisString      = "hours"
if (modl.eq."m") then
 res@gsnLeftString = "January - A1.0"
else
 res@gsnLeftString = "January - A1.1"
end if
res@gsnLeftStringFontHeightF    =  0.02
res@gsnLeftStringOrthogonalPosF = -0.10
res@gsnLeftStringParallelPosF   =  0.05

mxf    = (/0.,0.,0./)
;mxf!0  = "axis_0"
;mxf&axis_0 = (/1,2,3/)
mxf(0) = dble2flt(max(rad_ANN))
mxf(1) = dble2flt(max(sh_ANN))
mxf(2) = dble2flt(max(lh_ANN))
mnf    = (/0.,0.,0./)
;mnf!0  = "axis_0"
;mnf&axis_0 = (/1,2,3/)
mnf(0) = dble2flt(min(rad_ANN))
mnf(1) = dble2flt(min(sh_ANN))
mnf(2) = dble2flt(min(lh_ANN))
res@trYMinF         = floor(min(mnf))-5 ;-200
res@trYMaxF         = ceil(max(mxf)) +5 ;200
;res@trYMinF        = -200
;res@trYMaxF        = 200
plot(0)            = gsn_csm_xy(wks,y,rad_ANN,res)
;res@tiYAxisString = ""
plot(1)            = gsn_csm_xy(wks,y,sh_ANN,res)
 res@pmLegendDisplayMode    = "Always"            ; turn on legend
 res@pmLegendSide           = "Bottom"            ; Change location of
 res@pmLegendParallelPosF   = .15                 ; move units right
 res@pmLegendOrthogonalPosF = -0.55               ; move units down
 res@pmLegendWidthF         = 0.16                ; Change width and
 res@pmLegendHeightF        = 0.22                ; height of legend.
 res@lgPerimOn              = False               ; turn off box around
 res@lgLabelFontHeightF     = .015                ; label font height
 res@xyExplicitLegendLabels = (/"bond","harv","hyyt","loob","lwas","nbor","thar","vie","wbra","man","tumb"/)
plot(2)            = gsn_csm_xy(wks,y,lh_ANN,res)

; bv hv hy lo lw nb th vi wb ma tb
; Legend ---------------------------------------------------------------------

;lgres                    = True
;lgres@lgLineColors       = (/"black","blue","red","green","purple","cyan","magenta","firebrick","orange","SlateBlue","navy","SkyBlue"/)    ; colors for legend lines
;;lgres@lgDashIndex       = 0.
;lgres@vpWidthF           = 0.2                         ; width of legend (NDC)
;lgres@vpHeightF          = 0.5                         ; height of legend (NDC)
;lgres@lgPerimOn          = False                       ; turn off perimeter
;lgres@lgDashIndexes      = (/0,0,0,0,0,0,0,0,0,0,1,1/)  ; dash indexes
;lgres@lgLineLabelStrings = (/"",""/)                   ; no labels in lines
;;lst = (/"bondville","harvardforest","hyytiala","loobos","littlewashita","nsaboreas","tharandt","vielsalm","walkerbranch","hay","manaus","tumbarumba"/)
;lst  = (/"bv","hv","hy","lo","lw","nb","th","vi","wb","ma","tb"/)
;xpos = 0.9 ;(/0,0,0,0.25,0.25,0.25,0.5,0.5,0.5,0.75,0.75/)
;ypos = 0.9 ;(/0,.1,.2,0,.1,.2,0,.1,.2,0,.1/)
;
;lgd                      = gsn_create_legend(wks,11,lst,lgres)
;;gsn_legend_ndc(wks,11,lst,xpos,ypos,lgres)
;addres                   = True
;addres@amParallelPosF    = 0.05
;addres@amOrthogonalPosF  = 1.05
;add                      = gsn_add_annotation(plot(1),lgd,addres)

;-----------------------------------------------------------------------------


pres                       = True
pres@gsnMaximize           = True
pres@gsnPanelFigureStrings = (/"Rnet","SH","LH"/)
pres@gsnPanelFigureStringsBackgroundFillColor = "Transparent"

gsn_panel(wks,plot,(/1,3/),pres)

;frame(wks)

colors = (/"black","deepskyblue1","red","green","mediumorchid1","cyan","magenta","firebrick","orange","SlateBlue","navy","forestgreen"/)
;labels = (/"bondville","harvardforest","hyytiala","loobos","littlewashita","nsaboreas","tharandt","vielsalm","walkerbranch","manaus","tumbarumba"/)
labels = (/"Tumbarumba","Manaus","WalkerBranch","Vielsalm","Tharandt","NSABoreas","LittleWashita","Loobos","Hyytiala","HarvardForest","Bondville"/)

lgres                    = True
lgres@vpWidthF           = 0.4                         ; width of legend (NDC)
lgres@vpHeightF          = 0.5                         ; height of legend (NDC)
;;lgres@lgDashIndex       = 0.
;lgres@lgPerimOn          = False                       ; turn off perimeter
lgres@lgDashIndexes     = (/0,0,0,0,0,0,0,0,0,0,1,1/)  ; dash indexes
lgres@lgLineLabelStrings = (/"","","","","","","","","","",""/)                   ; no labels in lines
;lgres@lgLineLabelStrings =  labels
lgres@lgItemOrder        = (/10,9,8,7,6,5,4,3,2,1,0/)
lgres@lgOrientation      = "Horizontal"
lgres@lgLabelAngleF      = 270
lgres@gsnMaximize         = True
;lgres@gsnPaperOrientation = "landscape"
lgres@lgBoxMinorExtentF      = 0.2
;lgres@wkOrientation = "landscape"

xpos = 0.28  ;(/0.05,0.05,0.05,0.25,0.25,0.25,0.5,0.5,0.5,0.75,0.75/)
ypos = 0.95 ;(/0.05,.1,.2,0.05,.1,.2,0.05,.1,.2,0.05,.1/)

;lgd                      = gsn_create_legend(wks,11,lst,lgres)

;do i=0,10,1
lgres@lgLineColors = colors;(i)
;lgres@lgLineLabelFontColor = colors
;lgres@lgMonoLineLabelFontColor = True
;gsn_legend_ndc(wks,1,labels(i),xpos(i),ypos(i),lgres)
gsn_legend_ndc(wks,11,labels,xpos,ypos,lgres)
;end do

;draw(legend)
frame(wks)

print(" ")

if (jpeg) then
system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("convert -density 500 "+pltName+".eps "+pltName+".jpg")
;system("convert -density 1000 "+pltName+".ps "+pltName+".jpg")
;system("convert -density 1000 -geometry 1200x1600 "+pltName+".ps "+pltName+".jpg")
end if

end

