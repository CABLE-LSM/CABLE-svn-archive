;==============================================================================
; Lauren Stevens 12 Aug 2011
;==============================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;==============================================================================

;procedure mmdc(LongName:string,FlxName:string,mosDir:string,cabDir:string)

begin

hdir      = systemfunc("echo $hdir")
usr      = systemfunc("echo $USERID")
numYr    = systemfunc("echo $YR")
tstep    = stringtoint(systemfunc("echo $TSTEP"))
filDir   = systemfunc("echo $PWD")
;filDir   = systemfunc("echo $DIRW")
jpg     = systemfunc("echo $CHFMT")
if(jpg.eq."y") then
jpeg = True
else
jpeg = False
end if

home     = hdir+"/"+usr+"/"

;system("dmget " + filDir + "/mmdc_*_roll_"+numYr+"yrs.nc" )

;==============================================================================

fils2 = systemfunc("ls "+filDir+"/mmdc_*_roll_"+numYr+"yrs.nc")
;printVarSummary(fils2)
;print(dimsizes(fils2))
f2    = addfiles(fils2,"r")
ListSetType(f2,"join")
if (isfilevar(f2[0],"hfss")) then
msh   = f2[:]->hfss
else
msh   = f2[:]->sh
end if
if (isfilevar(f2[0],"hfls")) then
mlh   = f2[:]->hfls
else
mlh   = f2[:]->lh
end if
if (isfilevar(f2[0],"field3333")) then
mrad  = f2[:]->field3333
else
mrad  = f2[:]->field202
end if
if (isfilevar(f2[0],"field5226")) then
mpr   = f2[:]->field5226
else
mpr   = f2[:]->tot_precip
end if
if (isfilevar(f2[0],"tas")) then
mtmp  = f2[:]->tas
else
mtmp  = f2[:]->temp
end if

mpr       = mpr*tstep
mpr@units = "mm/day"

;==============================================================================

if (dimsizes(fils2).eq.18) then
 ; total = 18 
 stations = (/"afN","afS","bv","dw","dh","hp","hv","hay","hy","ind", \
              "lw","lo","ma","nb","th","tb","vi","wb"/)
 stNames  = (/"Africa: North of the Equator","Africa: South of the Equator","Bondville", \
              "Daly Waters","Dinghushan","Hapex (Bats)","Harvard Forest", "Hay",         \
              "Hyytiala","A Gridpt in India","Little Washita","Loobos","Manaus",         \
              "NSA Boreas","Tharandt","Tumbarumba","Vielsalm","Walker Branch"           /)
else
if (dimsizes(fils2).eq.31) then
 ; total = 31
 stations =(/"afN","afr","afS","alk","ant1","ant2","ant","ant3","bv","dw","dh",\
             "grnl","hp","hv","hay","hy","ind","lw","lo","ma","nb","rus1","rus2",\
             "sd1","sd2","sd3","sd4","th","tb","vi","wb"/)
 stNames  = (/"Africa: North of the Equator","Africa","Africa: South of the Equator",\
              "Alaska","Ant1","Ant2","Antarctica","Ant3","Bondville",\
              "Daly Waters","Dinghushan","Greenland","Hapex (Bats)","Harvard Forest","Hay",\
              "Hyytiala","A Gridpt in India","Little Washita","Loobos","Manaus",\
              "NSA Boreas","Russia1","Rus2","Saudi1","Sd2","Sd3","Sd4",\
              "Tharandt","Tumbarumba","Vielsalm","Walker Branch"/)
else
 print("Error in mmdc.ncl: Number of Sites and Files Do Not Match")
end if
end if

month    = (/"jan","jul"/)
mon      = (/0,6/)

;Plot==============================================================================

print(" ")
y = fspan(0,23.5,tstep)

do j=0,1,1
do i=0,dimsizes(stations)-1,1

pltName = "mmdc_"+stations(i)+"_"+month(j)+"_"+numYr+"yrs"
wks                   = gsn_open_wks("ps","mmdc_"+stations(i)+"_"+month(j)+"_"+numYr+"yrs")
res                   = True
res@gsnFrame          = False
res@gsnDraw           = False
res@xyLineColors      = "black"
res@xyDashPatterns    = 0
res@xyLineThicknesses = 2

;-----------------------------------------------------------------------------

print("Plotting "+month(j)+" MMDC for "+stNames(i)+" ("+stations(i)+")")

plot               = new(5,graphic)
res@tiYAxisString  = "W/m2"
plot(0)            = gsn_csm_xy(wks,y,mrad(i,mon(j),:),res)
plot(1)            = gsn_csm_xy(wks,y,msh(i,mon(j),:),res)
plot(2)            = gsn_csm_xy(wks,y,mlh(i,mon(j),:),res)
res@tiYAxisString  = "mm/day"
plot(3)            = gsn_csm_xy(wks,y,mpr(i,mon(j),:),res)
res@tiYAxisString  = "deg K"
plot(4)            = gsn_csm_xy(wks,y,mtmp(i,mon(j),:),res)

pres                       = True
pres@gsnMaximize           = True
pres@amJust                = "TopLeft"
pres@gsnPanelFigureStrings = (/"Rnet","SH","LH","Pr","TAS"/)

gsn_panel(wks,plot,(/2,3/),pres)

if(jpeg) then
system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("convert -density 500 -rotate 270 "+pltName+".eps "+pltName+".jpg")
;system("convert -density 1000 -rotate 270 "+pltName+".ps "+pltName+".jpg")
;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")
end if

end do 
end do

print(" ")

end

