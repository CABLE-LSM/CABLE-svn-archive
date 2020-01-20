;==============================================================================
; Lauren Stevens 8 July 2011
;==============================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
hdir = systemfunc("echo $hdir")
usr = systemfunc("echo $USERID")
loadscript(hdir+"/"+usr+"/umplot/plot_mmdc.ncl")
;==============================================================================

;procedure harv_ts()

begin

FlxName  = "harvard"
LongName = "harvard_forest"
;usr = systemfunc("echo $USERID")
numYr    = systemfunc("echo $YR")
tstep    = systemfunc("echo $TSTEP")
mossf    = systemfunc("echo $moshar")
cabsf    = systemfunc("echo $cabhar")
jpg = systemfunc("echo $CHFMT")
if(jpg.eq."y") then
jpeg = True
else
jpeg = False
end if

tsDir    = hdir+"/"+usr+"/umplot/obs/"
home     = hdir+"/"+usr+"/"

;system("dmget " + tsDir + "mmdc_*_roll.nc")
;flux1  = addfile(tsDir + "mmdc_hv92n_roll.nc" ,"r")
;flux2  = addfile(tsDir + "mmdc_hv93n_roll.nc" ,"r")
;flux3  = addfile(tsDir + "mmdc_hv94n_roll.nc" ,"r")
;flux4  = addfile(tsDir + "mmdc_hv95n_roll.nc" ,"r")
;flux5  = addfile(tsDir + "mmdc_hv96n_roll.nc" ,"r")
;flux6  = addfile(tsDir + "mmdc_hv97n_roll.nc" ,"r")
;flux7  = addfile(tsDir + "mmdc_hv98n_roll.nc" ,"r")
;flux8  = addfile(tsDir + "mmdc_hv99n_roll.nc" ,"r")
flux   = addfile(tsDir  + "fluxnet/mmdc_ts_Harvard.nc","r")
;flux   = addfile(home  + "obsdata/fluxnet/PALS/flxnt_mon/mmdc_ts_Harvard.nc","r")
mosfil = addfile(mossf,"r")
cabfil = addfile(cabsf,"r")
;igbp   = addfile(home  + "MMDC_"+ FlxName +"_roll.nc","r")

;==============================================================================

sensht1 = flux ->Qh
latht1  = flux ->Qle
rad1    = flux ->Rnet
temp1   = flux ->Tair
precip1 = flux ->Rainf

;sensht1 = flux1 ->H
;latht1  = flux1 ->LH
;rad1    = flux1 ->Rn
;temp1   = flux1 ->Ta
;precip1 = flux1 ->PPT
;sensht2 = flux2 ->H
;latht2  = flux2 ->LH
;rad2    = flux2 ->Rn
;temp2   = flux2 ->Ta
;precip2 = flux2 ->PPT
;sensht3 = flux3 ->H
;latht3  = flux3 ->LH
;rad3    = flux3 ->Rn
;temp3   = flux3 ->Ta
;precip3 = flux3 ->PPT
;sensht4 = flux4 ->H
;latht4  = flux4 ->LH
;rad4    = flux4 ->Rn
;temp4   = flux4 ->Ta
;precip4 = flux4 ->PPT
;sensht5 = flux5 ->H
;latht5  = flux5 ->LH
;rad5    = flux5 ->Rn
;temp5   = flux5 ->Ta
;precip5 = flux5 ->PPT
;sensht6 = flux6 ->H
;latht6  = flux6 ->LH
;rad6    = flux6 ->Rn
;temp6   = flux6 ->Ta
;precip6 = flux6 ->PPT
;sensht7 = flux7 ->H
;latht7  = flux7 ->LH
;rad7    = flux7 ->Rn
;temp7   = flux7 ->Ta
;precip7 = flux7 ->PPT
;sensht8 = flux8 ->H
;latht8  = flux8 ->LH
;rad8    = flux8 ->Rn
;temp8   = flux8 ->Ta
;precip8 = flux8 ->PPT
;---------------------
;qh      = igbp  ->Qh
;qle     = igbp  ->Qle
;tair    = igbp  ->Tair
;rainf   = igbp  ->Rainf
;----------------------
;csh     = cabfil ->hfss
;msh     = mosfil ->hfss
;clh     = cabfil ->hfls
;mlh     = mosfil ->hfls
;crad    = cabfil ->field3333
;mrad    = mosfil ->field3333
;cta     = cabfil ->tas
;mta     = mosfil ->tas
;cpr     = cabfil ->field5226
;mpr     = mosfil ->field5226

if (isfilevar(cabfil,"hfss")) then
 csh     = cabfil ->hfss
else
 csh     = cabfil ->sh
end if
if (isfilevar(mosfil,"hfss")) then
 msh     = mosfil ->hfss
else
 msh     = mosfil ->sh
end if
if (isfilevar(cabfil,"hfls")) then
 clh     = cabfil ->hfls
else
 clh     = cabfil ->lh
end if
if (isfilevar(mosfil,"hfls")) then
 mlh     = mosfil ->hfls
else
 mlh     = mosfil ->lh
end if
if (isfilevar(cabfil,"field3333")) then
 crad    = cabfil ->field3333
else
 crad    = cabfil ->field202
 cswr    = cabfil ->solar
 clwr    = cabfil ->longwave

 dNames  = getvardims(crad)
 dName2  = getvardims(cswr)
 dims    = dimsizes(crad)

 crn1    = cswr+clwr
 copy_VarMeta(cswr,crn1)
 crn1&$dName2(1)$ = ispan(0,(tstep-1),6)
 crn2    = ftcurvp(crn1&$dName2(1)$,crn1,tstep,crad&$dNames(1)$)

 crnt    = new((/dims(0),dims(1)/),typeof(cswr))
 crnt(:,44:47)  = crn2(:,0:3)
 crnt(:,0:43)   = crn2(:,4:47)
 crnt!0         = "month"
 crnt!1         = "timestep"
end if
if (isfilevar(mosfil,"field3333")) then
 mrad    = mosfil ->field3333
 ;radsum = sum(mrad)
 ;print((/radsum/))
else
 mrad    = mosfil ->field202
 mswr    = mosfil ->solar
 mlwr    = mosfil ->longwave

 dNames  = getvardims(mrad)
 dName2  = getvardims(mswr)
 dims    = dimsizes(mrad)

 mrn1    = mswr+mlwr
 copy_VarMeta(mswr,mrn1)
 mrn1&$dName2(1)$ = ispan(0,(tstep-1),6)
 mrn2    = ftcurvp(mrn1&$dName2(1)$,mrn1,tstep,mrad&$dNames(1)$)

 ; Roll 4ts for Hyytiala
 mrnt    = new((/dims(0),dims(1)/),typeof(mswr))
 mrnt(:,44:47)  = mrn2(:,0:3)
 mrnt(:,0:43)   = mrn2(:,4:47)
 mrnt!0         = "month"
 mrnt!1         = "timestep"
end if
if (isfilevar(cabfil,"tas")) then
 cta     = cabfil ->tas
else
 cta     = cabfil ->temp
end if
if (isfilevar(mosfil,"tas")) then
 mta     = mosfil ->tas
else
 mta     = mosfil ->temp
end if
if (isfilevar(cabfil,"field5226")) then
 cpr     = cabfil ->field5226
else
 cpr     = cabfil ->tot_precip
end if
if (isfilevar(mosfil,"field5226")) then
 mpr     = mosfil ->field5226
else
 mpr     = mosfil ->tot_precip
end if
if (isfilevar(cabfil,"field30")) then
 ctc     = cabfil ->field30
;else
 ;ctc     = cabfil ->
end if
if (isfilevar(mosfil,"field30")) then
 mtc     = mosfil ->field30
;else
 ;mtc     = mosfil ->
end if

;==============================================================================

dNam     = getvardims(msh)
tstp     = msh&$dNam(1)$
tstep    = max(tstp)+1
timestep = doubletoint(tstep)

sensht1@_FillValue = -9999
;sensht2@_FillValue = -9999
;sensht3@_FillValue = -9999
;sensht4@_FillValue = -9999
;sensht5@_FillValue = -9999
;sensht6@_FillValue = -9999
;sensht7@_FillValue = -9999
;sensht8@_FillValue = -9999
sens_all           = sensht1(:,:,0) ;(sensht1+sensht2+sensht3+sensht4+sensht5+sensht6+sensht7+sensht8)/8
;sens_2             = (sensht2+sensht3)/2
copy_VarMeta(sensht1,sens_all)
;copy_VarMeta(sensht1,sens_2)
latht1@_FillValue  = -9999
;latht2@_FillValue  = -9999
;latht3@_FillValue  = -9999
;latht4@_FillValue  = -9999
;latht5@_FillValue  = -9999
;latht6@_FillValue  = -9999
;latht7@_FillValue  = -9999
;latht8@_FillValue  = -9999
lh_all             = latht1(:,:,0) ;(latht1+latht2+latht3+latht4+latht5+latht6+latht7+latht8)/8
;lh_2               = (latht2+latht3)/2
copy_VarMeta(latht1,lh_all)
;copy_VarMeta(latht1,lh_2)
;print(sensht1(0,0,0:10))
rad1@_FillValue    = -9999
;rad2@_FillValue    = -9999
;rad3@_FillValue    = -9999
;rad4@_FillValue    = -9999
;rad5@_FillValue    = -9999
;rad6@_FillValue    = -9999
;rad7@_FillValue    = -9999
;rad8@_FillValue    = -9999
rad_all            = rad1(:,:,0) ;(rad1+rad2+rad3+rad4+rad5+rad6+rad7+rad8)/8
;rad_2              = rad3 ;(rad2+rad3)/2
copy_VarMeta(rad1,rad_all)
;copy_VarMeta(rad1,rad_2)
temp1@_FillValue   = -9999
;temp2@_FillValue   = -9999
;temp3@_FillValue   = -9999
;temp4@_FillValue   = -9999
;temp5@_FillValue   = -9999
;temp6@_FillValue   = -9999
;temp7@_FillValue   = -9999
;temp8@_FillValue   = -9999
temp_all           = temp1(:,:,0) ;(temp1+temp2+temp3+temp4+temp5+temp6+temp7+temp8)/8 +273.16
;temp_2             = (temp2+temp3)/2 +273.16
copy_VarMeta(temp1,temp_all)
;copy_VarMeta(temp1,temp_2)
precip1@_FillValue = -9999
;precip2@_FillValue = -9999
;precip3@_FillValue = -9999
;precip4@_FillValue = -9999
;precip5@_FillValue = -9999
;precip6@_FillValue = -9999
;precip7@_FillValue = -9999
;precip8@_FillValue = -9999
pr_all             = precip1(:,:,0)*24*3600. ;(precip1+precip2+precip3+precip4+precip5+precip6+precip7+precip8)/8 *tstep
pr_all             = -9999
;pr_2               = (precip2+precip3)/2 *tstep
copy_VarMeta(precip1,pr_all)
;copy_VarMeta(precip1,pr_2)

sens     = sens_all;dim_avg_n_Wrap(sens_all,0)
;sens2    = dim_avg_n_Wrap(sens_2,0)
sens!0   = "month"
sens!1   = "timestep"
;sens2!0  = "month"
;sens2!1  = "timestep"
;qh_ta    = dim_avg_n_Wrap(qh,0)
;qh_ta!0  = "month"
;qh_ta!1  = "timestep"
msh!0    = "month"
msh!1    = "timestep"
csh!0    = "month"
csh!1    = "timestep"

lats     = lh_all;dim_avg_n_Wrap(lh_all,0)
;lats2    = dim_avg_n_Wrap(lh_2,0)
lats!0   = "month"
lats!1   = "timestep"
;lats2!0  = "month"
;lats2!1  = "timestep"
;qle_ta   = dim_avg_n_Wrap(qle,0)
;qle_ta!0 = "month"
;qle_ta!1 = "timestep"
mlh!0    = "month"
mlh!1    = "timestep"
clh!0    = "month"
clh!1    = "timestep"

temp      = temp_all;dim_avg_n_Wrap(temp_all,0)
;tem2      = dim_avg_n_Wrap(temp_2,0)
temp!0    = "month"
temp!1    = "timestep"
;tem2!0    = "month"
;tem2!1    = "timestep"
;tair_ta   = dim_avg_n_Wrap(tair,0)
;tair_ta!0 = "month"
;tair_ta!1 = "timestep"
mta!0     = "month"
mta!1     = "timestep"
cta!0     = "month"
cta!1     = "timestep"

prec      = pr_all;dim_avg_n_Wrap(pr_all,0)
;prec2     = dim_avg_n_Wrap(pr_2,0)
prec!0    = "month"
prec!1    = "timestep"
;prec2!0   = "month"
;prec2!1   = "timestep"
;rf_ta     = dim_avg_n_Wrap(rainf,0)
;rf_ta!0   = "month"
;rf_ta!1   = "timestep"
;rainf     = rainf*24*60*60
mpr       = mpr *timestep
cpr       = cpr *timestep
mpr!0     = "month"
mpr!1     = "timestep"
cpr!0     = "month"
cpr!1     = "timestep"

rnet      = rad_all;dim_avg_n_Wrap(rad_all,0)
;rnet2     = dim_avg_n_Wrap(rad_2,0)
rnet!0    = "month"
rnet!1    = "timestep"
;rnet2!0   = "month"
;rnet2!1   = "timestep"
;rad_ta    = dim_avg_n_Wrap(,0)
;rad_ta!0  = "month"
;rad_ta!1  = "timestep"
mrad!0    = "month"
mrad!1    = "timestep"
crad!0    = "month"
crad!1    = "timestep"

;----------------------------------------------------------------------------

y = fspan(0,23.5,timestep)

; July
data      = new((/3,timestep/),double)
;data(0,:) = -9999 ;qh_ta(6,:)            ; igbp
data(0,:) = msh(6,:)              ; moses
data(1,:) = csh(6,:)              ; cable
data(2,:) = sens(6,:)             ; fluxnet

data1      = new((/3,timestep/),double)
;data1(0,:) = -9999 ;qle_ta(6,:)           ; igbp
data1(0,:) = mlh(6,:)              ; moses
data1(1,:) = clh(6,:)              ; cable
data1(2,:) = lats(6,:)             ; fluxnet

data7      = new((/3,timestep/),double)
;data7(0,:) = -9999 ;tair_ta(6,:)          ; igbp
data7(0,:) = mta(6,:)              ; moses
data7(1,:) = cta(6,:)              ; cable
data7(2,:) = temp(6,:)             ; fluxnet

data8      = new((/3,timestep/),double)
;data8(0,:) = -9999 ;rad_ta(6,:)           ; igbp
if (sum(mrad(6,:)).ne.0) then
data8(0,:) = mrad(6,:)             ; moses
else
data8(0,:) = mrnt(6,:)             ; moses
end if
data8(1,:) = crnt(6,:);crad(6,:)             ; cable
data8(2,:) = rnet(6,:)             ; fluxnet

data9      = new((/3,timestep/),double)
;data9(0,:) = -9999 ;rf_ta(6,:)            ; igbp
data9(0,:) = mpr(6,:)              ; moses
data9(1,:) = cpr(6,:)              ; cable
data9(2,:) = prec(6,:)             ; fluxnet

if (isvar("ctc").and.isvar("mtc")) then
data14      = new((/3,timestep/),double)
;data14      = -9999                 ; igbp
;data14@_FillValue = -9999           ;
data14(0,:) = mtc(6,:)              ; moses
data14(1,:) = ctc(6,:)              ; cable
;data14(2,:) = -9999                 ; fluxnet
end if

; Jan
data2      = new((/3,timestep/),double)
;data2(0,:) = -9999 ;qh_ta(0,:)            ; igbp
data2(0,:) = msh(0,:)              ; moses
data2(1,:) = csh(0,:)              ; cable
data2(2,:) = sens(0,:)            ; fluxnet

data3      = new((/3,timestep/),double)
;data3(0,:) = -9999 ;qle_ta(0,:)           ; igbp
data3(0,:) = mlh(0,:)              ; moses
data3(1,:) = clh(0,:)              ; cable
data3(2,:) = lats(0,:)            ; fluxnet

data4      = new((/3,timestep/),double)
;data4(0,:) = -9999 ;tair_ta(0,:)          ; igbp
data4(0,:) = mta(0,:)              ; moses
data4(1,:) = cta(0,:)              ; cable
data4(2,:) = temp(0,:)             ; fluxnet

data5      = new((/3,timestep/),double)
;data5(0,:) = -9999 ;rad_ta(0,:)           ; igbp
if (sum(mrad(0,:)).ne.0) then
data5(0,:) = mrad(0,:)             ; moses
else
data5(0,:) = mrnt(0,:)             ; moses
end if
data5(1,:) = crnt(0,:);crad(0,:)             ; cable
data5(2,:) = rnet(0,:)            ; fluxnet

data6      = new((/3,timestep/),double)
;data6(0,:) = -9999 ;rf_ta(0,:)            ; igbp
data6(0,:) = mpr(0,:)              ; moses
data6(1,:) = cpr(0,:)              ; cable
data6(2,:) = prec(0,:)            ; fluxnet

if (isvar("ctc").and.isvar("mtc")) then
data12      = new((/3,timestep/),double)
;data12      = -9999                 ; igbp
;data12@_FillValue = -9999           ;
data12(0,:) = mtc(0,:)              ; moses
data12(1,:) = ctc(0,:)              ; cable
;data12(2,:) = -9999                 ; fluxnet
end if

;Plot-------------------------------------------------------------------------

print("Plotting Harvard Forest MMDC with Fluxnet")
print("")

ts_flux(FlxName,data2,data3,data,data1,jpeg)
ts_rad(FlxName,data5,data4,data8,data7,jpeg)
if (isvar("mtc") .and. isvar("ctc")) then
 ts_prec4(FlxName,data6,data12,data9,data14,jpeg)
else
 ts_prec(FlxName,data6,data9,jpeg)
end if

end
