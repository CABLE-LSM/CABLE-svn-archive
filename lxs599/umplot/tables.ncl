;==================================================================
; Lauren Stevens 17 Jan 2013
;==================================================================

;ste69f/umplot> grep -il gsn_table *
;global_means_carbon.ncl
;global_means_wbal.ncl
;wbal.ncl

;ste69f/umplot> grep -il gsn_table */*
;dcc_all_scripts/amoj_ebal_b.ncl
;dcc_all_scripts/amoj_ebal.ncl
;dcc_all_scripts/amoj_ecomponents_b.ncl
;dcc_all_scripts/amoj_tmp_b.ncl
;dcc_all_scripts/amoj_wbal_b.ncl
;dcc_all_scripts/amoj_wbal.ncl
;dcc_all_scripts/era_wbal_gmeans.ncl
;dcc_all_scripts/global_means_carbon.ncl
;dcc_all_scripts/global_means_wbal.ncl
;dcc_all_scripts/tables.ncl
;dcc_all_scripts/wbal.ncl

;scripts_dcc/amoj_ebal_b.ncl
;scripts_dcc/amoj_ebal.ncl
;scripts_dcc/amoj_wbal_b.ncl
;scripts_dcc/amoj_wbal.ncl
;scripts_dcc/wbal.ncl

;==================================================================

;CONTAINS
;procedure ebal_tab(land_vals,ocean_vals,global_vals,numYr,region)
;procedure wbal_tab(land_vals,ocean_vals,global_vals,numYr,region)
;procedure temp_tab(land_vals,ocean_vals,global_vals,numYr,region)
;procedure ecomps_compr(land_vals,ocean_vals,global_vals,numYr,region)
;procedure temp_compr(land_vals,ocean_vals,global_vals,numYr,region)
;procedure ebal_compr(land_vals,ocean_vals,global_vals,numYr,region)
;procedure wbal_compr(land_vals,ocean_vals,global_vals,numYr,region)

;==================================================================

procedure ebal_tab(land_vals,ocean_vals,global_vals,numYr,region)

begin

; swd swu sh lh lwu lwd rn
; TEST FIRST

land  = "" + decimalPlaces(land_vals,2,True)
ocean = "" + decimalPlaces(ocean_vals,2,True)
globe = "" + decimalPlaces(global_vals,2,True)

; Header-------------------------

ncr1 = (/1,8/)
x1   = (/0.080,.995/)
y1   = (/0.800,0.855/)
text1= (/"SW down","SW up","SH","LH",\
         "LW up","LW down","Rnet","eBal"/)

;ncr1_1= (/1,1/)
;x1_1  = (/0.080,.2325/)
;y1_1  = (/0.800,0.855/)
;text11= (/"SW down Wm~S~-2~N~"/)

;ncr1_2= (/1,1/)
;x1_2  = (/0.2325,.3850/)
;y1_2  = (/0.800,0.855/)
;text12= (/"SW up Wm~S~-2~N~"/)

;ncr1_3= (/1,1/)
;x1_3  = (/0.3850,0.5375/)
;y1_3  = (/0.800,0.855/)
;text13= (/"SH Wm~S~-2~N~"/)

;ncr1_4= (/1,1/)
;x1_4  = (/0.5375,0.69/)
;y1_4  = (/0.800,0.855/)
;text14= (/"LH Wm~S~-2~N~"/)

;ncr1_5= (/1,1/)
;x1_5  = (/0.69,0.8425/)
;y1_5  = (/0.800,0.855/)
;text15= (/"LW up Wm~S~-2~N~"/)

;ncr1_6= (/1,1/)
;x1_6  = (/0.8425,0.995/)
;y1_6  = (/0.800,0.855/)
;text16= (/"LW down Wm~S~-2~N~"/)

;ncr1_7= (/1,1/)
;x1_7  = (/0.8425,0.995/)
;y1_7  = (/0.800,0.855/)
;text17= (/"Rnet Wm~S~-2~N~"/)

;ncr1_8= (/1,1/)
;x1_8  = (/0.8425,0.995/)
;y1_8  = (/0.800,0.855/)
;text18= (/"Bal Wm~S~-2~N~"/)

; First Column-------------------
ncr2  = (/3,1/)
x2    = (/0.005,0.080/)
y2    = (/0.415,0.800/)
text2 = (/"Land","Ocean","Global"/)

ncr8  = (/1,1/)
x8    = (/0.005,0.080/)
y8    = (/0.800,0.855/)
text8 = (/region/)

; Main Table Body----------------
ncr3        = (/3,8/)
x3          = (/0.080,0.995/)
y3          = (/0.415,0.800/)
text_values = (/land,ocean,globe/)

pltName = "ebal_table" + "_" + numYr + "yrs"
wks     = gsn_open_wks("ps",pltName)

res1               = True
res1@txFontHeightF = 0.012
;res1@gsFillColor   = (/"deepskyblue2","lightgoldenrod1","slateblue2","lightseagreen","tomato","lightgoldenrod1","tomato","lightgoldenrod1","tomato","grey","grey","grey"/)
res2               = True
res2@txFontHeightF = 0.012
;res2@gsFillColor   = (/"forestgreen","deepskyblue2","darkorchid2"/)
res3               = True
;res3@gsnMaximize   = True
res3@txFontHeightF = 0.012
;res3@gsFillColor   = (/ (/"white","white","white","white","white","gray","gray","gray","gray","white","white","white"/), \
;                        (/"white","gray","gray","gray","gray","gray","gray","gray","gray","white","white","white"/),     \
;                        (/"white","gray","white","white","gray","white","gray","gray","gray","white","white","white"/)   /)

print("Creating Energy Components Table")
gsn_table(wks,ncr1,x1,y1,text1,res1)
gsn_table(wks,ncr8,x8,y8,text8,res1)
;gsn_table(wks,ncr1_1,x1_1,y1_1,text11,res1)
;gsn_table(wks,ncr1_2,x1_2,y1_2,text12,res1)
;gsn_table(wks,ncr1_3,x1_3,y1_3,text13,res1)
;gsn_table(wks,ncr1_4,x1_4,y1_4,text14,res1)
;gsn_table(wks,ncr1_5,x1_5,y1_5,text15,res1)
;gsn_table(wks,ncr1_6,x1_6,y1_6,text16,res1)
gsn_table(wks,ncr2,x2,y2,text2,res2)
gsn_table(wks,ncr3,x3,y3,text_values,res3)

system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("rm "+pltName+".ps")
;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")


end

;==================================================================

procedure wbal_tab(land_vals,ocean_vals,global_vals,numYr,region)

begin

land  = "" + decimalPlaces(land_vals,3,True)
ocean = "" + decimalPlaces(ocean_vals,3,True)
globe = "" + decimalPlaces(global_vals,3,True)

; Header-------------------------
ncr1  = (/1,7/) ;9;12
x1    = (/0.085,0.995/) ;.104
y1    = (/0.800,0.855/)
text1 = (/"Precip","Evap","SurfRun","Drainage","wBal", \
          "Evp(LH)","wBal2"/);, \
          ;"E(LHsnw)","wBal3"/);,\
          ;"SH","LH","Rnet"/)

; First Column-------------------
ncr2  = (/3,1/)
x2    = (/0.005,0.085/)
y2    = (/0.415,0.800/)
text2 = (/"Land","Ocean","Global"/)

ncr8  = (/1,1/)
x8    = (/0.005,0.085/)
y8    = (/0.800,0.855/)
text8 = (/region/)

; Main Table Body----------------
ncr3        = (/3,7/) ;9;12
x3          = (/0.085,0.995/)
y3          = (/0.415,0.800/)
text_values = (/land,ocean,globe/)

pltName = "wBal_table" + "_" + numYr + "yrs"
wks     = gsn_open_wks("ps",pltName)

res1               = True
res1@txFontHeightF = 0.012
;res1@gsFillColor   = (/"deepskyblue2","lightgoldenrod1","slateblue2","lightseagreen","tomato","lightgoldenrod1","tomato","lightgoldenrod1","tomato"/);,"lightskyblue","lightskyblue","lightskyblue"/)
res2               = True
res2@txFontHeightF = 0.012
;res2@gsFillColor   = (/"forestgreen","deepskyblue2","darkorchid2"/)
res3               = True
;res3@gsnMaximize   = True
res3@txFontHeightF = 0.012
res3@gsFillColor   = (/ (/"white","gray","white","white","gray","white","white"/), \
                        (/"white","gray","gray","gray","gray","white","white"/),     \
                        (/"white","gray","gray","gray","gray","white","gray"/)   /)
                        ;(/"white","gray","gray","gray","gray","white","white","white","white"/)   /)
                        ;(/"white","gray","white","white","gray","white","white","white","white"/)   /)
;res3@gsFillColor   = (/ (/"white","white","white","white","white","white","gray","gray","gray"/), \
;                        (/"white","gray","gray","gray","gray","white","gray","gray","gray"/),     \
;                        (/"white","gray","gray","gray","gray","white","gray","gray","gray"/)   /)
;                        ;(/"white","gray","gray","gray","gray","white","white","white","white"/)   /)
;                        ;(/"white","gray","white","white","gray","white","white","white","white"/)   /)
;res3@gsFillColor   = (/ (/"white","white","white","white","white","gray","gray","gray","gray","white","white","white"/), \
;                        (/"white","gray","gray","gray","gray","gray","gray","gray","gray","white","white","white"/),     \
;                        (/"white","gray","gray","gray","gray","white","gray","gray","gray","white","white","white"/)   /)
;                        ;(/"white","gray","gray","gray","gray","white","white","white","white"/)   /)
;                        ;(/"white","gray","white","white","gray","white","white","white","white"/)   /)

print("Creating Water Balance Table")
gsn_table(wks,ncr1,x1,y1,text1,res1)
gsn_table(wks,ncr8,x8,y8,text8,res1)
gsn_table(wks,ncr2,x2,y2,text2,res2)
gsn_table(wks,ncr3,x3,y3,text_values,res3)

system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("rm "+pltName+".ps")
;system("convert -density 1000 -rotate 270 "+pltName+".ps "+pltName+".jpg")
;;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")

end

;==================================================================

procedure temp_tab(land_vals,ocean_vals,global_vals,numYr,region)

begin

; tscrn tmax tmin (tamp)
; TEST FIRST

land  = "" + decimalPlaces(land_vals,3,True)
ocean = "" + decimalPlaces(ocean_vals,3,True)
globe = "" + decimalPlaces(global_vals,3,True)

; Header-------------------------

ncr1 = (/1,3/)
x1   = (/0.085,.995/)
y1   = (/0.800,0.855/)
text1= (/"Screen Temp","Tmax","Tmin"/)

;ncr1_1= (/1,1/)
;x1_1  = (/0.085,.475/)
;y1_1  = (/0.800,0.855/)
;text11= (/"Screen Temp ~S~o~N~C"/)

;ncr1_2= (/1,1/)
;x1_2  = (/0.475,.735/)
;y1_2  = (/0.800,0.855/)
;text12= (/"Max. Screen Temp ~S~o~N~C"/)

;ncr1_3= (/1,1/)
;x1_3  = (/0.735,0.995/)
;y1_3  = (/0.800,0.855/)
;text13= (/"Min. Screen Temp ~S~o~N~C"/)

; First Column-------------------
ncr2  = (/3,1/)
x2    = (/0.005,0.085/)
y2    = (/0.415,0.800/)
text2 = (/"Land","Ocean","Global"/)

ncr8  = (/1,1/)
x8    = (/0.005,0.085/)
y8    = (/0.800,0.855/)
text8 = (/region/)

; Main Table Body----------------
ncr3        = (/3,3/)
x3          = (/0.085,0.995/)
y3          = (/0.415,0.800/)
text_values = (/land,ocean,globe/)

pltName = "temp_table" + "_" + numYr + "yrs"
wks     = gsn_open_wks("ps",pltName)

res1               = True
res1@txFontHeightF = 0.012
;res1@gsFillColor   = (/"deepskyblue2","lightgoldenrod1","slateblue2","lightseagreen","tomato","lightgoldenrod1","tomato","lightgoldenrod1","tomato","grey","grey","grey"/)
res2               = True
res2@txFontHeightF = 0.012
;res2@gsFillColor   = (/"forestgreen","deepskyblue2","darkorchid2"/)
res3               = True
;res3@gsnMaximize   = True
res3@txFontHeightF = 0.012
;res3@gsFillColor   = (/ (/"white","white","white","white","white","gray","gray","gray","gray","white","white","white"/), \
;                        (/"white","gray","gray","gray","gray","gray","gray","gray","gray","white","white","white"/),     \
;                        (/"white","gray","white","white","gray","white","gray","gray","gray","white","white","white"/)   /)

print("Creating Temp. Table")
gsn_table(wks,ncr1,x1,y1,text1,res1)
gsn_table(wks,ncr8,x8,y8,text8,res1)
;gsn_table(wks,ncr1_1,x1_1,y1_1,text11,res1)
;gsn_table(wks,ncr1_2,x1_2,y1_2,text12,res1)
;gsn_table(wks,ncr1_3,x1_3,y1_3,text13,res1)
;gsn_table(wks,ncr1_4,x1_4,y1_4,text14,res1)
;gsn_table(wks,ncr1_5,x1_5,y1_5,text15,res1)
gsn_table(wks,ncr2,x2,y2,text2,res2)
gsn_table(wks,ncr3,x3,y3,text_values,res3)

system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("rm "+pltName+".ps")
;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")

end

;==================================================================

procedure ecomps_compr(land_vals,ocean_vals,global_vals,numYr,region,modv)

begin

land  = "" + decimalPlaces(land_vals,1,True)
ocean = "" + decimalPlaces(ocean_vals,1,True)
globe = "" + decimalPlaces(global_vals,1,True)

if (modv.eq."1") then
source = (/"A1.0","A1.1","Est","A1.0","A1.1","Est",\
           "A1.0","A1.1","Est","A1.0","A1.1","Est", \
           "A1.0","A1.1","Est","A1.0","A1.1","Est"/)
else
source = (/"A1.0","A1.3","Est","A1.0","A1.3","Est",\
           "A1.0","A1.3","Est","A1.0","A1.3","Est", \
           "A1.0","A1.3","Est","A1.0","A1.3","Est"/)
end if

; Header-------------------------

ncr1_1= (/1,1/)
x1_1  = (/0.080,.2325/)
y1_1  = (/0.800,0.855/)
text11= (/"SW down Wm~S~-2~N~"/)

ncr1_2= (/1,1/)
x1_2  = (/0.2325,.3850/)
y1_2  = (/0.800,0.855/)
text12= (/"SW up Wm~S~-2~N~"/)

ncr1_3= (/1,1/)
x1_3  = (/0.3850,0.5375/)
y1_3  = (/0.800,0.855/)
text13= (/"SH Wm~S~-2~N~"/)

ncr1_4= (/1,1/)
x1_4  = (/0.5375,0.69/)
y1_4  = (/0.800,0.855/)
text14= (/"LH Wm~S~-2~N~"/)

ncr1_5= (/1,1/)
x1_5  = (/0.69,0.8425/)
y1_5  = (/0.800,0.855/)
text15= (/"LW up Wm~S~-2~N~"/)

ncr1_6= (/1,1/)
x1_6  = (/0.8425,0.995/)
y1_6  = (/0.800,0.855/)
text16= (/"LW down Wm~S~-2~N~"/)

; First Column-------------------
ncr2  = (/4,1/)
x2    = (/0.005,0.080/)
y2    = (/0.415,0.800/)
text2 = (/"Source","Land","Ocean","Global"/)

ncr8  = (/1,1/)
x8    = (/0.005,0.080/)
y8    = (/0.800,0.855/)
text8 = (/region/)

; Main Table Body----------------
ncr3        = (/4,18/)
x3          = (/0.080,0.995/)
y3          = (/0.415,0.800/)
text_values = (/source,land,ocean,globe/)

pltName = "ecomps_ctable" + "_" + numYr + "yrs"
wks     = gsn_open_wks("ps",pltName)

res1               = True
res1@txFontHeightF = 0.012
;res1@gsFillColor   = (/"deepskyblue2","lightgoldenrod1","slateblue2","lightseagreen","tomato","lightgoldenrod1","tomato","lightgoldenrod1","tomato","grey","grey","grey"/)
res2               = True
res2@txFontHeightF = 0.012
;res2@gsFillColor   = (/"forestgreen","deepskyblue2","darkorchid2"/)
res3               = True
;res3@gsnMaximize   = True
res3@txFontHeightF = 0.012
;res3@gsFillColor   = (/ (/"white","white","white","white","white","gray","gray","gray","gray","white","white","white"/), \
;                        (/"white","gray","gray","gray","gray","gray","gray","gray","gray","white","white","white"/),     \
;                        (/"white","gray","white","white","gray","white","gray","gray","gray","white","white","white"/)   /)

print("Creating Energy Components Table")
;gsn_table(wks,ncr1,x1,y1,text1,res1)
gsn_table(wks,ncr1_1,x1_1,y1_1,text11,res1)
gsn_table(wks,ncr1_2,x1_2,y1_2,text12,res1)
gsn_table(wks,ncr1_3,x1_3,y1_3,text13,res1)
gsn_table(wks,ncr1_4,x1_4,y1_4,text14,res1)
gsn_table(wks,ncr1_5,x1_5,y1_5,text15,res1)
gsn_table(wks,ncr1_6,x1_6,y1_6,text16,res1)
gsn_table(wks,ncr8,x8,y8,text8,res2)
gsn_table(wks,ncr2,x2,y2,text2,res2)
gsn_table(wks,ncr3,x3,y3,text_values,res3)

system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("rm "+pltName+".ps")
;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")

end

;==================================================================

procedure temp_compr(land_vals,ocean_vals,global_vals,numYr,region,modv)

begin

land  = "" + decimalPlaces(land_vals,2,True)
ocean = "" + decimalPlaces(ocean_vals,2,True)
globe = "" + decimalPlaces(global_vals,2,True)

if (modv.eq."1") then
source = (/"A1.0","A1.1","NOAA","A1.0","A1.1",\
           "A1.0","A1.1"/)
else
source = (/"A1.0","A1.3","NOAA","A1.0","A1.3",\
           "A1.0","A1.3"/)
end if

; Header-------------------------

ncr1_1= (/1,1/)
x1_1  = (/0.085,.475/)
y1_1  = (/0.800,0.855/)
text11= (/"Screen Temp ~S~o~N~C"/)

ncr1_2= (/1,1/)
x1_2  = (/0.475,.735/)
y1_2  = (/0.800,0.855/)
text12= (/"Max. Screen Temp ~S~o~N~C"/)

ncr1_3= (/1,1/)
x1_3  = (/0.735,0.995/)
y1_3  = (/0.800,0.855/)
text13= (/"Min. Screen Temp ~S~o~N~C"/)

; First Column-------------------
ncr2  = (/4,1/)
x2    = (/0.005,0.085/)
y2    = (/0.415,0.800/)
text2 = (/"Source","Land","Ocean","Global"/)

ncr8  = (/1,1/)
x8    = (/0.005,0.085/)
y8    = (/0.800,0.855/)
text8 = (/region/)

; Main Table Body----------------
ncr3        = (/4,7/)
x3          = (/0.085,0.995/)
y3          = (/0.415,0.800/)
text_values = (/source,land,ocean,globe/)

pltName = "temp_ctable" + "_" + numYr + "yrs"
wks     = gsn_open_wks("ps",pltName)

res1               = True
res1@txFontHeightF = 0.012
;res1@gsFillColor   = (/"deepskyblue2","lightgoldenrod1","slateblue2","lightseagreen","tomato","lightgoldenrod1","tomato","lightgoldenrod1","tomato","grey","grey","grey"/)
res2               = True
res2@txFontHeightF = 0.012
;res2@gsFillColor   = (/"forestgreen","deepskyblue2","darkorchid2"/)
res3               = True
;res3@gsnMaximize   = True
res3@txFontHeightF = 0.012
;res3@gsFillColor   = (/ (/"white","white","white","white","white","gray","gray","gray","gray","white","white","white"/), \
;                        (/"white","gray","gray","gray","gray","gray","gray","gray","gray","white","white","white"/),     \
;                        (/"white","gray","white","white","gray","white","gray","gray","gray","white","white","white"/)   /)

print("Creating Temp. Table")
;gsn_table(wks,ncr1,x1,y1,text1,res1)
gsn_table(wks,ncr1_1,x1_1,y1_1,text11,res1)
gsn_table(wks,ncr1_2,x1_2,y1_2,text12,res1)
gsn_table(wks,ncr1_3,x1_3,y1_3,text13,res1)
;gsn_table(wks,ncr1_4,x1_4,y1_4,text14,res1)
;gsn_table(wks,ncr1_5,x1_5,y1_5,text15,res1)
gsn_table(wks,ncr8,x8,y8,text8,res2)
gsn_table(wks,ncr2,x2,y2,text2,res2)
gsn_table(wks,ncr3,x3,y3,text_values,res3)

system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("rm "+pltName+".ps")
;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")

end

;==================================================================

procedure ebal_compr(land_vals,ocean_vals,global_vals,numYr,region,modv)

begin

land  = "" + decimalPlaces(land_vals,2,True)
ocean = "" + decimalPlaces(ocean_vals,2,True)
globe = "" + decimalPlaces(global_vals,2,True)

if (modv.eq."1") then
source = (/"A1.0","A1.1","NOAA","A1.0","A1.1","Bud",\
           "Henn","A1.0","A1.1","Bud","Henn", \
           "A1.0","A1.1","Bud","Henn"/)
else
source = (/"A1.0","A1.3","NOAA","A1.0","A1.3","Bud",\
           "Henn","A1.0","A1.3","Bud","Henn", \
           "A1.0","A1.3","Bud","Henn"/)
end if

; Header-------------------------

ncr1_1= (/1,1/)
x1_1  = (/0.085,.267/)
y1_1  = (/0.800,0.855/)
text11= (/"Screen Temp ~S~o~N~C"/)

ncr1_2= (/1,1/)
x1_2  = (/0.267,.5097/)
y1_2  = (/0.800,0.855/)
text12= (/"Sensible Heat Wm~S~-2~N~"/)

ncr1_3= (/1,1/)
x1_3  = (/0.5097,0.7523/)
y1_3  = (/0.800,0.855/)
text13= (/"Latent Heat Wm~S~-2~N~"/)

ncr1_4= (/1,1/)
x1_4  = (/0.7523,0.995/)
y1_4  = (/0.800,0.855/)
text14= (/"Net Radiation Wm~S~-2~N~"/)

; First Column-------------------
ncr2  = (/4,1/)
x2    = (/0.005,0.085/)
y2    = (/0.415,0.800/)
text2 = (/"Source","Land","Ocean","Global"/)

ncr8  = (/1,1/)
x8    = (/0.005,0.085/)
y8    = (/0.800,0.855/)
text8 = (/region/)

; Main Table Body----------------
ncr3        = (/4,15/)
x3          = (/0.085,0.995/)
y3          = (/0.415,0.800/)
text_values = (/source,land,ocean,globe/)

pltName = "ebal_ctable" + "_" + numYr + "yrs"
wks     = gsn_open_wks("ps",pltName)

res1               = True
res1@txFontHeightF = 0.012
;res1@gsFillColor   = (/"deepskyblue2","lightgoldenrod1","slateblue2","lightseagreen","tomato","lightgoldenrod1","tomato","lightgoldenrod1","tomato","grey","grey","grey"/)
res2               = True
res2@txFontHeightF = 0.012
;res2@gsFillColor   = (/"forestgreen","deepskyblue2","darkorchid2"/)
res3               = True
;res3@gsnMaximize   = True
res3@txFontHeightF = 0.012
;res3@gsFillColor   = (/ (/"white","white","white","white","white","gray","gray","gray","gray","white","white","white"/), \
;                        (/"white","gray","gray","gray","gray","gray","gray","gray","gray","white","white","white"/),     \
;                        (/"white","gray","white","white","gray","white","gray","gray","gray","white","white","white"/)   /)

print("Creating Energy Balance Table")
;gsn_table(wks,ncr1,x1,y1,text1,res1)
gsn_table(wks,ncr1_1,x1_1,y1_1,text11,res1)
gsn_table(wks,ncr1_2,x1_2,y1_2,text12,res1)
gsn_table(wks,ncr1_3,x1_3,y1_3,text13,res1)
gsn_table(wks,ncr1_4,x1_4,y1_4,text14,res1)
;gsn_table(wks,ncr1_5,x1_5,y1_5,text15,res1)
gsn_table(wks,ncr8,x8,y8,text8,res2)
gsn_table(wks,ncr2,x2,y2,text2,res2)
gsn_table(wks,ncr3,x3,y3,text_values,res3)

system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("rm "+pltName+".ps")
;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")

end

;==================================================================

procedure wbal_compr(land_vals,ocean_vals,global_vals,numYr,region,modv)

begin

land  = "" + decimalPlaces(land_vals,2,True)
ocean = "" + decimalPlaces(ocean_vals,2,True)
globe = "" + decimalPlaces(global_vals,2,True)
;print(land)
;print(ocean)
;print(globe)

if (modv.eq."1") then
source = (/"A1.0","A1.1","B+R","L+W","A1.0","A1.1", \
           "B+R","A1.0","A1.1","A1.0","A1.1","B+R"/)
else
source = (/"A1.0","A1.3","B+R","L+W","A1.0","A1.3", \
           "B+R","A1.0","A1.3","A1.0","A1.3","B+R"/)
end if

; Header-------------------------

;ncr1  = (/1,1/)
;x1    = (/0.005,.085/)
;y1    = (/0.800,0.855/)
;text1 = (/"mm/d"/)

ncr1_1= (/1,1/)
x1_1  = (/0.085,.3883/)
y1_1  = (/0.800,0.855/)
text11= (/"Precipitation"/)

ncr1_2= (/1,1/)
x1_2  = (/0.3883,.6158/)
y1_2  = (/0.800,0.855/)
text12= (/"Evaporation"/)

ncr1_3= (/1,1/)
x1_3  = (/0.6158,0.7675/)
y1_3  = (/0.800,0.855/)
text13= (/"Surface Runoff"/)

ncr1_4= (/1,1/)
x1_4  = (/0.7675,0.9192/)
y1_4  = (/0.800,0.855/)
text14= (/"Drainage"/)

ncr1_5= (/1,1/)
x1_5  = (/0.9192,0.995/)
y1_5  = (/0.800,0.855/)
text15= (/"TotRun"/)

; First Column-------------------
ncr2  = (/4,1/)
x2    = (/0.005,0.085/)
y2    = (/0.415,0.800/)
text2 = (/"Source","Land","Ocean","Global"/)

ncr8  = (/1,1/)
x8    = (/0.005,0.085/)
y8    = (/0.800,0.855/)
text8 = (/region/)

; Main Table Body----------------
ncr3        = (/4,12/)
x3          = (/0.085,0.995/)
y3          = (/0.415,0.800/)
text_values = (/source,land,ocean,globe/)

pltName = "wbal_ctable" + "_" + numYr + "yrs"
wks     = gsn_open_wks("ps",pltName)

res1               = True
res1@txFontHeightF = 0.012
;res1@gsFillColor   = (/"deepskyblue2","lightgoldenrod1","slateblue2","lightseagreen","tomato","lightgoldenrod1","tomato","lightgoldenrod1","tomato","grey","grey","grey"/)
res2               = True
res2@txFontHeightF = 0.012
;res2@gsFillColor   = (/"forestgreen","deepskyblue2","darkorchid2"/)
res3               = True
;res3@gsnMaximize   = True
res3@txFontHeightF = 0.012
;res3@gsFillColor   = (/ (/"white","white","white","white","white","gray","gray","gray","gray","white","white","white"/), \
;                        (/"white","gray","gray","gray","gray","gray","gray","gray","gray","white","white","white"/),     \
;                        (/"white","gray","white","white","gray","white","gray","gray","gray","white","white","white"/)   /)

print("Creating Water Balance Table")
;gsn_table(wks,ncr1,x1,y1,text1,res1)
gsn_table(wks,ncr1_1,x1_1,y1_1,text11,res1)
gsn_table(wks,ncr1_2,x1_2,y1_2,text12,res1)
gsn_table(wks,ncr1_3,x1_3,y1_3,text13,res1)
gsn_table(wks,ncr1_4,x1_4,y1_4,text14,res1)
gsn_table(wks,ncr1_5,x1_5,y1_5,text15,res1)
gsn_table(wks,ncr8,x8,y8,text8,res2)
gsn_table(wks,ncr2,x2,y2,text2,res2)
gsn_table(wks,ncr3,x3,y3,text_values,res3)

system("ps2epsi "+pltName+".ps "+pltName+".eps")
system("rm "+pltName+".ps")
;system("convert -density 1000 -geometry 1200x1600 -rotate 270 "+pltName+".ps "+pltName+".jpg")

end
