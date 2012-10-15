;==============================================================================
; Lauren Stevens 19 July 2011
;==============================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "/home/cmar/ste69f/umplot/plot_mmdc.ncl"
;==============================================================================

;procedure mmdc_flux(LongName:string,FlxName:string,mosDir:string,cabDir:string)

begin

numYr    = systemfunc("echo $YR")
filDir   = systemfunc("echo $PWD")
jpg      = systemfunc("echo $JPEG")
if(jpg.eq."y") then
jpeg = True
else
jpeg = False
end if

tsDir    = "/home/cmar/ste69f/umplot/obs/"
home     = "/home/cmar/ste69f/"

;system("dmget " + tsDir + "mmdc*_roll.nc" )
;hayfil = addfile(tsDir + "hayobsjuly.nc"                ,"r")
manfil = addfile(tsDir + "mmdc_ma96_roll.nc"         ,"r")
tbfil  = addfile(tsDir + "MMDC_tumbarumba_fv.nc" ,"r")
vi97   = addfile(tsDir + "mmdc_vi97_roll.nc"       ,"r")

;igbp     = addfile(home  + "MMDC_"+ LongName +"_roll.nc","r")

;==============================================================================

fils   = systemfunc("ls "+tsDir+"mmdc_??.nc")
f      = addfiles(fils,"r")
ListSetType(f,"join")
sens   = f[:]->H
lat    = f[:]->LH
radnet = f[:]->Rn
;temp2m = f[:]->Ta
;preci  = f[:]->PPT

sens2   = dim_avg_n_Wrap(sens,1)
lat2    = dim_avg_n_Wrap(lat,1)
radnet2 = dim_avg_n_Wrap(radnet,1)

fils2  = systemfunc("ls "+filDir+"/mmdc_*_roll_"+numYr+"yrs.nc")
f2     = addfiles(fils2,"r")
ListSetType(f2,"join")
msh  = f2[:]->hfss
mlh  = f2[:]->hfls
mrad = f2[:]->field3333

sens_ma = manfil ->H
lat_ma  = manfil ->LH
rad_ma  = manfil ->Rn
;temp_ma = manfil ->Ta
;prec_ma = manfil ->PPT
sens_tb = tbfil  ->Qh    ;H
lat_tb  = tbfil  ->Qle   ;LH
rad_tb  = tbfil  ->Rnet  ;Rn
;temp_tb = tbfil  ->Tair  ;Ta
;prec_tb = tbfil  ->Rainf ;PPT
;sens_hy  = hayfil ->h
;lat_hy   = hayfil ->lh
rad_vi   = vi97  ->Rn
sens_vi  = vi97  ->H

sens_tb2 = dim_avg_n_Wrap(sens_tb,0)
lat_tb2  = dim_avg_n_Wrap(lat_tb,0)
rad_tb2  = dim_avg_n_Wrap(rad_tb,0)

sens_ma2 = dim_avg_n_Wrap(sens_ma,0)
lat_ma2  = dim_avg_n_Wrap(lat_ma,0)
rad_ma2  = dim_avg_n_Wrap(rad_ma,0)

sh_tb = new((/12,48/),double)
lh_tb = new((/12,48/),double)
rn_tb = new((/12,48/),double)
;sh_hy = new(48,double)
;lh_hy = new(48,double)
do i = 0,23
sh_tb(:,2*i+1)=sens_tb2(:,i)
sh_tb(:,2*i)  =sens_tb2(:,i)
lh_tb(:,2*i+1)=lat_tb2(:,i)
lh_tb(:,2*i)  =lat_tb2(:,i)
rn_tb(:,2*i+1)=rad_tb2(:,i)
rn_tb(:,2*i)  =rad_tb2(:,i)
;sh_hy(2*i+1)  =sens_hy(i)
;sh_hy(2*i)    =sens_hy(i)
;lh_hy(2*i+1)  =lat_hy(i)
;lh_hy(2*i)    =lat_hy(i)
end do

;==============================================================================

sensht            = new((/11,48/),double)
latht             = new((/11,48/),double)
rad               = new((/11,48/),double)
;temp              = new((/11,12/),double)
;precip            = new((/11,12/),double)
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
;temp(k,:)         = temp2m(k,0:47)      ;Ta
;precip(k,:)       = preci(k,0:47)       ;PPT
end do
sensht!0 = "axis_1"
latht!0  = "axis_1"
rad!0    = "axis_1"
;sh_hy!0  = "axis_2"
;lh_hy!0  = "axis_2"
;sensht(9,:)     = sh_hy(0:47)
;latht(9,:)      = lh_hy(0:47)
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

rad(2,:)    = rad(2,:)*3 +2*9999
rad(4,:)    = (rad(4,:)*3 + 9999)/2
rad(5,:)    = (rad(5,:)*5+9999)/4
rad(7,:)    = rad(7,:)*3 +2*9999
sensht(2,:) = (sensht(2,:)*3 + 9999)/2
sensht(4,:) = (sensht(4,:)*3 + 9999)/2
sensht(5,:) = (sensht(5,:)*5+9999)/4
;sensht(6,:) = (sensht(6,:)*4 - sens_th())/3 
sensht(7,:) = (sensht(7,:)*3 + 9999 - sens_vi(0,0,:))
latht(2,:)  = (latht(2,:)*3 + 9999)/2
latht(4,:)  = (latht(4,:)*3 + 9999)/2
latht(5,:)  = (latht(5,:)*5 + 9999)/4
;latht(6,:)  = (latht(6,:)*4 - lh_th())/3 
latht(7,:)  = (latht(7,:)*3 + 9999)/2

mrad_ro = new((/11,48/),double)
mrad_ro(0,:)  = mrad(2,0,:)  ; bv ; black
mrad_ro(1,:)  = mrad(6,0,:)  ; hv ; blue
mrad_ro(2,:)  = mrad(8,0,:)  ; hy ; red
mrad_ro(3,:)  = mrad(10,0,:) ; lw ; green
mrad_ro(4,:)  = mrad(11,0,:) ; lo ; darkviolet
mrad_ro(5,:)  = mrad(13,0,:) ; nb ; cyan
mrad_ro(6,:)  = mrad(14,0,:) ; th ; magenta
mrad_ro(7,:)  = mrad(16,0,:) ; vi ; firebrick
mrad_ro(8,:)  = mrad(17,0,:) ; wb ; orange
mrad_ro(9,:)  = mrad(12,0,:) ; ma ; SlateBlue
mrad_ro(10,:) = mrad(15,0,:) ; tb ; navy       and hay = forestgreen

msh_ro = new((/11,48/),double)
msh_ro(0,:)  = msh(2,0,:)  ; bv
msh_ro(1,:)  = msh(6,0,:)  ; hv
msh_ro(2,:)  = msh(8,0,:)  ; hy
msh_ro(3,:)  = msh(10,0,:) ; lw
msh_ro(4,:)  = msh(11,0,:) ; lo
msh_ro(5,:)  = msh(13,0,:) ; nb
msh_ro(6,:)  = msh(14,0,:) ; th
msh_ro(7,:)  = msh(16,0,:) ; vi
msh_ro(8,:)  = msh(17,0,:) ; wb
msh_ro(9,:)  = msh(12,0,:) ; ma
msh_ro(10,:) = msh(15,0,:) ; tb

mlh_ro = new((/11,48/),double)
mlh_ro(0,:)  = mlh(2,0,:)  ; bv
mlh_ro(1,:)  = mlh(6,0,:)  ; hv
mlh_ro(2,:)  = mlh(8,0,:)  ; hy
mlh_ro(3,:)  = mlh(10,0,:) ; lw
mlh_ro(4,:)  = mlh(11,0,:) ; lo
mlh_ro(5,:)  = mlh(13,0,:) ; nb
mlh_ro(6,:)  = mlh(14,0,:) ; th
mlh_ro(7,:)  = mlh(16,0,:) ; vi
mlh_ro(8,:)  = mlh(17,0,:) ; wb
mlh_ro(9,:)  = mlh(12,0,:) ; ma
mlh_ro(10,:) = mlh(15,0,:) ; tb

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

print("Plotting Jan MMDC Biases for all Sites")

y = fspan(0,23.5,48)
;y = ispan(1,12,1)

;ts_flux(LongName,data2,data3,data,data1)
;ts_rad(LongName,data5,data4,data8,data7)
;ts_prec(LongName,data6,data9)

pltName               = "mmdc_jan_"+numYr+"yrs"
wks                   = gsn_open_wks("ps",pltName)
res                   = True
res@gsnFrame          = False
res@gsnDraw           = False
res@xyLineColors      = (/"black","blue","red","green","darkviolet","cyan","magenta","firebrick","orange","SlateBlue","navy","forestgreen"/)
res@xyDashPatterns    = (/0,0,0,0,0,0,0,0,0,0,1,1/)
res@xyLineThicknesses = (/2,2,2,2,2,2,2,2,2,2,2,2/)
res@gsnYRefLine       = 0.

;-----------------------------------------------------------------------------

;res@gsnMaximize    = True
plot               = new(3,graphic)
res@tiYAxisString  = "W/m2"
res@trYMinF        = -200
res@trYMaxF        = 200
plot(0)            = gsn_csm_xy(wks,y,rad_ANN,res)
;res@tiYAxisString = ""
plot(1)            = gsn_csm_xy(wks,y,sh_ANN,res)
plot(2)            = gsn_csm_xy(wks,y,lh_ANN,res)

; bv hv hy lw lo nb th vi wb ma tb
; Legend ---------------------------------------------------------------------

;lgres                    = True
;lgres@lgLineColors       = (/"black","blue","red","green","purple","cyan","magenta","firebrick","orange","SlateBlue","navy","SkyBlue"/)    ; colors for legend lines
;;lgres@lgDashIndex       = 0.
;lgres@vpWidthF           = 0.2                         ; width of legend (NDC)
;lgres@vpHeightF          = 0.5                         ; height of legend (NDC)
;lgres@lgPerimOn          = False                       ; turn off perimeter
;lgres@lgDashIndexes      = (/0,0,0,0,0,0,0,0,0,0,1,1/)  ; dash indexes
;lgres@lgLineLabelStrings = (/"",""/)                   ; no labels in lines
;;lst = (/"bondville","harvardforest","hyytiala","littlewashita","loobos","nsaboreas","tharandt","vielsalm","walkerbranch","hay","manaus","tumbarumba"/)
;lst  = (/"bv","hv","hy","lw","lo","nb","th","vi","wb","ma","tb"/)
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

colors = (/"black","blue","red","green","darkviolet","cyan","magenta","firebrick","orange","SlateBlue","navy","forestgreen"/)
;labels = (/"bondville","harvardforest","hyytiala","littlewashita","loobos","nsaboreas","tharandt","vielsalm","walkerbranch","manaus","tumbarumba"/)
labels = (/"Tumbarumba","Manaus","WalkerBranch","Vielsalm","Tharandt","NSABoreas","Loobos","LittleWashita","Hyytiala","HarvardForest","Bondville"/)

lgres                    = True
lgres@vpWidthF           = 0.4                         ; width of legend (NDC)
lgres@vpHeightF          = 0.5                         ; height of legend (NDC)
;;lgres@lgLineColors       = (/"black","blue","red","green","purple","cyan","magenta","firebrick","orange","SlateBlue","navy","SkyBlue"/)    ; colors for legend lines
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

if (jpeg) then
system("convert -density 1000 -geometry 1200x1600 "+pltName+".ps "+pltName+".jpg")
end if

end

