;==============================================================================
; Lauren Stevens 8 July 2011
;==============================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/home/cmar/ste69f/umplot/plot_mmdc.ncl"
;==============================================================================

;procedure mmdc_flux(LongName:string,FlxName:string,mosDir:string,cabDir:string)

begin

FlxName  = "tharandt"
LongName = "tharandt"
numYr = systemfunc("echo $YR")
mossf = systemfunc("echo $mostha") ; mossf = systemfunc("echo $"+mosDir+"")
cabsf = systemfunc("echo $cabtha") ; cabsf = systemfunc("echo $"+cabDir+"")
jpg = systemfunc("echo $CHFMT")
if(jpg.eq."y") then
jpeg = True
else
jpeg = False
end if

tsDir = "/home/cmar/ste69f/umplot/obs/"
home  = "/home/cmar/ste69f/"

system("dmget " + tsDir + "mmdc_*_roll.nc" )
;flux  = addfiles()
flux1  = addfile(tsDir + "mmdc_th96e_roll.nc" ,"r")
flux2  = addfile(tsDir + "mmdc_th97e_roll.nc" ,"r")
flux3  = addfile(tsDir + "mmdc_th98e_roll.nc" ,"r")
flux4  = addfile(tsDir + "mmdc_th99e_roll.nc" ,"r")
mosfil = addfile(mossf,"r")
cabfil = addfile(cabsf,"r")
;igbp   = addfile(tsDir + "MMDC_"+LongName+"_roll.nc","r")

;==============================================================================

sensht1 = flux1 ->H
latht1  = flux1 ->LH
rad1    = flux1 ->Rn
temp1   = flux1 ->Ta
precip1 = flux1 ->PPT
sensht2 = flux2 ->H
latht2  = flux2 ->LH
rad2    = flux2 ->Rn
temp2   = flux2 ->Ta
precip2 = flux2 ->PPT
sensht3 = flux3 ->H
latht3  = flux3 ->LH
rad3    = flux3 ->Rn
temp3   = flux3 ->Ta
precip3 = flux3 ->PPT
sensht4 = flux4 ->H
latht4  = flux4 ->LH
rad4    = flux4 ->Rn
temp4   = flux4 ->Ta
precip4 = flux4 ->PPT
;---------------------
;qh      = igbp  ->Qh
;qle     = igbp  ->Qle
;tair    = igbp  ->Tair
;rainf   = igbp  ->Rainf
;netrad  = igbp  ->Rnet
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
end if
if (isfilevar(mosfil,"field3333")) then
mrad    = mosfil ->field3333
else
 mrad    = mosfil ->field202
 mswr    = mosfil ->solar
 mlwr    = mosfil ->longwave

 dNames  = getvardims(mrad)
 dName2  = getvardims(mswr)
 dims    = dimsizes(mrad)

 mrn1    = mswr+mlwr
 copy_VarMeta(mswr,mrn1)
 mrn1&$dName2(1)$ = ispan(0,47,6)
 mrn2    = linint1_n_Wrap(mrn1&$dName2(1)$,mrn1,True,mrad&$dNames(1)$,0,1)

 ; Roll 4ts for Hyytiala
 mrnt    = new((/dims(0),dims(1)/),typeof(mswr))
 mrnt(:,0:4)  = mrn2(:,43:47)
 mrnt(:,5:47) = mrn2(:,0:42)
mrnt!0    = "month"
mrnt!1    = "timestep"
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

;==============================================================================

dNam     = getvardims(msh)
tstp     = msh&$dNam(1)$
tstep    = max(tstp)+1
timestep = doubletoint(tstep)

sensht1@_FillValue = -9999
sensht2@_FillValue = -9999
sensht3@_FillValue = -9999
sensht4@_FillValue = -9999
sens_all           = (sensht1+sensht2+sensht3+sensht4)/4
sens_2             = (sensht2+sensht3+sensht4)/3
copy_VarMeta(sensht1,sens_all)
copy_VarMeta(sensht1,sens_2)
latht1@_FillValue  = -9999
latht2@_FillValue  = -9999
latht3@_FillValue  = -9999
latht4@_FillValue  = -9999
lh_all             = (latht1+latht2+latht3+latht4)/4
lh_2               = (latht2+latht3+latht4)/3
copy_VarMeta(latht1,lh_all)
copy_VarMeta(latht1,lh_2)
rad1@_FillValue    = -9999
rad2@_FillValue    = -9999
rad3@_FillValue    = -9999
rad4@_FillValue    = -9999
rad_all            = (rad1+rad2+rad3+rad4)/4
rad_2              = (rad2+rad3+rad4)/3
copy_VarMeta(rad1,rad_all)
copy_VarMeta(rad1,rad_2)
temp1@_FillValue   = -9999
temp2@_FillValue   = -9999
temp3@_FillValue   = -9999
temp4@_FillValue   = -9999
temp_all           = (temp1+temp2+temp3+temp4)/4 +273.16
temp_2             = (temp2+temp3+temp4)/3 +273.16
copy_VarMeta(temp1,temp_all)
copy_VarMeta(temp1,temp_2)
precip1@_FillValue = -9999
precip2@_FillValue = -9999
precip3@_FillValue = -9999
precip4@_FillValue = -9999
pr_all             = (precip1+precip2+precip3+precip4)/4 *48
pr_2               = (precip3+precip4)/2 *48
copy_VarMeta(precip1,pr_all)
copy_VarMeta(precip1,pr_2)

;qh_ta    = dim_avg_n_Wrap(qh,0)
sens     = dim_avg_n_Wrap(sens_all,0)
sens2    = dim_avg_n_Wrap(sens_2,0)
;qh_ta!0  = "month"
;qh_ta!1  = "timestep"
sens!0   = "month"
sens!1   = "timestep"
sens2!0  = "month"
sens2!1  = "timestep"
msh!0    = "month"
msh!1    = "timestep"
csh!0    = "month"
csh!1    = "timestep"

;qle_ta   = dim_avg_n_Wrap(qle,0)
lats     = dim_avg_n_Wrap(lh_all,0)
lats2    = dim_avg_n_Wrap(lh_2,0)
;qle_ta!0 = "month"
;qle_ta!1 = "timestep"
lats!0   = "month"
lats!1   = "timestep"
lats2!0  = "month"
lats2!1  = "timestep"
mlh!0    = "month"
mlh!1    = "timestep"
clh!0    = "month"
clh!1    = "timestep"

;tair_ta   = dim_avg_n_Wrap(tair,0)
temp      = dim_avg_n_Wrap(temp_all,0)
tem2      = dim_avg_n_Wrap(temp_2,0)
;tair_ta!0 = "month"
;tair_ta!1 = "timestep"
temp!0    = "month"
temp!1    = "timestep"
tem2!0    = "month"
tem2!1    = "timestep"
mta!0     = "month"
mta!1     = "timestep"
cta!0     = "month"
cta!1     = "timestep"

;rainf     = rainf*24*60*60
mpr       = mpr *timestep
cpr       = cpr *timestep
;rf_ta     = dim_avg_n_Wrap(rainf,0)
;rf_ta2    = dim_avg_n_Wrap(rainf(1:3,:,:),0)
prec      = dim_avg_n_Wrap(pr_all,0)
prec2     = dim_avg_n_Wrap(pr_2,0)
;rf_ta!0   = "month"
;rf_ta!1   = "timestep"
;rf_ta2!0  = "month"
;rf_ta2!1  = "timestep"
prec!0    = "month"
prec!1    = "timestep"
prec2!0   = "month"
prec2!1   = "timestep"
mpr!0     = "month"
mpr!1     = "timestep"
cpr!0     = "month"
cpr!1     = "timestep"

;rad_ta    = dim_avg_n_Wrap(netrad,0)
rnet      = dim_avg_n_Wrap(rad_all,0)
rnet2     = dim_avg_n_Wrap(rad_2,0)
;rad_ta!0  = "month"
;rad_ta!1  = "timestep"
rnet!0    = "month"
rnet!1    = "timestep"
rnet2!0   = "month"
rnet2!1   = "timestep"
mrad!0    = "month"
mrad!1    = "timestep"
crad!0    = "month"
crad!1    = "timestep"

;----------------------------------------------------------------------------

y = fspan(0,23.5,timestep)

; July
data      = new((/3,timestep/),double)
;data(0,:) = qh_ta(6,:)            ; igbp
data(0,:) = msh(6,:)              ; moses
data(1,:) = csh(6,:)              ; cable
data(2,:) = sens(6,:)             ; fluxnet

data1      = new((/3,timestep/),double)
;data1(0,:) = qle_ta(6,:)           ; igbp
data1(0,:) = mlh(6,:)              ; moses
data1(1,:) = clh(6,:)              ; cable
data1(2,:) = lats(6,:)             ; fluxnet

data7      = new((/3,timestep/),double)
;data7(0,:) = tair_ta(6,:)          ; igbp
data7(0,:) = mta(6,:)              ; moses
data7(1,:) = cta(6,:)              ; cable
data7(2,:) = temp(6,:)             ; fluxnet

data8      = new((/3,timestep/),double)
;data8(0,:) = rad_ta(6,:)           ; igbp
if (sum(mrad(6,:)).ne.0) then
data8(0,:) = mrad(6,:)             ; moses
else
data8(0,:) = mrnt(6,:)             ; moses
end if
data8(1,:) = crad(6,:)             ; cable
data8(2,:) = rnet(6,:)             ; fluxnet

data9      = new((/3,timestep/),double)
;data9(0,:) = rf_ta(6,:)            ; igbp
data9(0,:) = mpr(6,:)              ; moses
data9(1,:) = cpr(6,:)              ; cable
data9(2,:) = prec(6,:)             ; fluxnet

; Jan
data2      = new((/3,timestep/),double)
;data2(0,:) = qh_ta(0,:)            ; igbp
data2(0,:) = msh(0,:)              ; moses
data2(1,:) = csh(0,:)              ; cable
data2(2,:) = sens2(0,:)            ; fluxnet

data3      = new((/3,timestep/),double)
;data3(0,:) = qle_ta(0,:)           ; igbp
data3(0,:) = mlh(0,:)              ; moses
data3(1,:) = clh(0,:)              ; cable
data3(2,:) = lats2(0,:)            ; fluxnet

data4      = new((/3,timestep/),double)
;data4(0,:) = tair_ta(0,:)          ; igbp
data4(0,:) = mta(0,:)              ; moses
data4(1,:) = cta(0,:)              ; cable
data4(2,:) = temp(0,:)             ; fluxnet

data5      = new((/3,timestep/),double)
;data5(0,:) = rad_ta(0,:)           ; igbp
if (sum(mrad(0,:)).ne.0) then
data5(0,:) = mrad(0,:)             ; moses
else
data5(0,:) = mrnt(0,:)             ; moses
end if
data5(1,:) = crad(0,:)             ; cable
data5(2,:) = rnet(0,:)            ; fluxnet

data6      = new((/3,timestep/),double)
;data6(0,:) = rf_ta2(0,:)            ; igbp
data6(0,:) = mpr(0,:)              ; moses
data6(1,:) = cpr(0,:)              ; cable
data6(2,:) = prec2(0,:)            ; fluxnet

;Plot-------------------------------------------------------------------------

print("Plotting Tharandt MMDC with Fluxnet")
print("")

ts_flux(LongName,data2,data3,data,data1,jpeg)
ts_rad(LongName,data5,data4,data8,data7,jpeg)
ts_prec(LongName,data6,data9,jpeg)

end

