;*************************************************
; Based on leg_7.ncl
;
; Concepts illustrated:
;   - Drawing a scatter plot on a map
;   - Changing the marker color and size in a map plot
;   - Plotting station locations using markers
;   - Manually creating a legend using markers and text
;   - Adding text to a plot
;   - Generating dummy data using "random_uniform"
;
;*************************************************
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
;-------Options--------

  colors = (/1,2,8,10,12,15,17,20,22,25/);,26,30,35,38,43,48,52,56,61,66,74,79,84,89,94,98/)  
                                        ; marker colors, dimsizes must 
                                        ; be equal to dimsizes()+1
  labels = new(10,string)  ; Labels for legend.
  ;labels = new(26,string)  ; Labels for legend.
;---------------------------
  npts = 10 ;26 ;30                                  ; Number of points.
  ;lat  = (/-2.06,-34.5,-16.26,43.68,-80,-80,-80,23.12,61.84,5,-5,58.3,58.59,-14.2,-35.7,50.9,40.11,\
  ;         23.2,34.9,35.96,52.17,42.54,55.88,50.3,65,-75/) 
  ;lon  = (/-60.21,144.93,133.38,-0.1,25,135,-105,83.2,24.3,20,20,-134.42,49.66,23.28,148.2,13.56,\
  ;         -88.37,112.5,-97.97,-84.28,5.74,-72.17,-98.48,5.99,-45,-90/) 
  lat  = (/61.3 ,27.5,27.5,22.5,-75,65 ,-13.8,58.8,58.8,21.3/)
  lon  = (/110.6,46.9,43.1,43.1,-90,-45, 22.5,-135,48.8,48.8/)

  R   = random_uniform(-1.2,35.,npts)   ; This is dummy data for determining
                                        ; how to color the markers.

;------------------------------ 
; Labels
;
  ;labels(0)  = "Man" ; "Manaus"
  ;labels(1)  = "Hay"
  ;labels(2)  = "DaW" ; "Daly Waters"
  ;labels(3)  = "Cau"
  ;labels(4)  = "Hap" ; "Bats (Hapex)"
  ;labels(5)  = "Aus" ; "Antarctica - Aus"
  ;labels(6)  = "Afr" ; "Antarctica - Afr"
  ;labels(7)  = "SAm" ; "Antarctica - SAm"
  ;labels(8)  = "Ind" ; "India"
  ;labels(9)  = "Hyy" ; "Hyytiala"
  ;labels(10) = "AfN" ; "Africa N Eq"
  ;labels(11) = "AfS" ; "Africa S Eq"
  ;labels(12) = "Alk" ; "Alaska"
  ;labels(13) = "Rus" ; "Russia"
  ;labels(14) = "Afr" ; "Africa"
  ;labels(15) = "Tum" ; "Tumbarumba"
  ;labels(16) = "Tha" ; "Tharandt"
  ;labels(17) = "Bon" ; "Bondville"
  ;labels(18) = "Din" ; "Dinghushan"
  ;labels(19) = "LiW" ; "Litte Washita"
  ;labels(20) = "WBr" ; "Walker Branch"
  ;labels(21) = "EK1"
  ;labels(22) = "EK2"
  ;labels(23) = "EK3"
  ;labels(24) = "Lob" ; "Loobos"
  ;labels(25) = "Har" ; "Harvard Forest"
  ;labels(26) = "Bor" ; "Boreas NSA"
  ;labels(27) = "Vie" ; "Vielsalm"
  ;labels(28) = "GrL" ; "Greenland"
  ;labels(29) = "Ant" ; "Antarctica"

  ;labels(0)  = "Man" ; "Manaus"
  ;labels(1)  = "Hay"
  ;labels(2)  = "DaW" ; "Daly Waters"
  ;labels(3)  = "Hap" ; "Bats (Hapex)"
  ;labels(4)  = "Aus" ; "Antarctica - Aus"
  ;labels(5)  = "Afr" ; "Antarctica - Afr"
  ;labels(6)  = "SAm" ; "Antarctica - SAm"   ;;;;
  ;labels(7)  = "Ind" ; "India"
  ;labels(8)  = "Hyy" ; "Hyytiala"
  ;labels(9)  = "AfN" ; "Africa N Eq"
  ;labels(10) = "AfS" ; "Africa S Eq"
  ;labels(11) = "Alk" ; "Alaska"
  ;labels(12) = "Rus" ; "Russia"
  ;labels(13) = "Afr" ; "Africa"             ;;;;
  ;labels(14) = "Tum" ; "Tumbarumba"
  ;labels(15) = "Tha" ; "Tharandt"
  ;labels(16) = "Bon" ; "Bondville"
  ;labels(17) = "Din" ; "Dinghushan"
  ;labels(18) = "LiW" ; "Litte Washita"
  ;labels(19) = "WBr" ; "Walker Branch"
  ;labels(20) = "Lob" ; "Loobos"             ;;;;
  ;labels(21) = "Har" ; "Harvard Forest"
  ;labels(22) = "Bor" ; "Boreas NSA"
  ;labels(23) = "Vie" ; "Vielsalm"
  ;labels(24) = "GrL" ; "Greenland"
  ;labels(25) = "Ant" ; "Antarctica"

  labels(0) = "RusW" ; "Russia"
  labels(1) = "SdA1" ; "Saudi A"
  labels(2) = "SdA2" ; "Saudi A"
  labels(3) = "SdA3" ; "Saudi A"
  labels(4) = "Ant"  ; "Antarctica"
  labels(5) = "GrL"  ; "Greenland"
  labels(6) = "Afr"  ; "Africa"
  labels(7) = "Alk"  ; "Alaska"
  labels(8) = "RusE" ; "Russia"
  labels(9) = "SdA4" ; "Saudi A"

;===========================================================================
; Begin plotting section.
;
  wks = gsn_open_wks("ps","map_2")                 ; Open a workstation and
  gsn_define_colormap(wks,"WhViBlGrYeOrRe")      ; define a different colormap.

  nc1 = NhlNewColor(wks,.8,.8,.8)        ; Add light gray to colormap,
                                         ; for continents.
;
; Set up some map resources.
;
  mpres                = True
  mpres@gsnMaximize    = True             ; Maximize plot in frame.
  mpres@gsnPanelBottom = 0.18
  mpres@gsnFrame       = False            ; Don't advance the frame
;
  mpres@mpFillColors = (/-1,-1,nc1,-1/)    ;assign light gray to land masses

;  mpres@tiMainString = "Data colored and sized according to range of values"
  map = gsn_csm_map(wks,mpres)
;
  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.

  txres               = True
  txres@txFontHeightF = 0.015
;
; Loop through each grouping of markers, and draw them one set at
; a time, assigning the proper color and size with gsn_marker.
;
; At the same time, draw a legend showing the meaning of the markers.
;

  xleg = (/0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.32,0.32,0.32,0.32,0.32,0.32,0.32,0.57,0.57,0.57,0.57,0.57,0.57,0.57,0.82,0.82,0.82,0.82,0.82,0.82,0.82/)                              ; Location of
  xtxt = (/0.16,0.16,0.16,0.16,0.16,0.16,0.16,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.66,0.66,0.66,0.66,0.66,0.66,0.66,0.91,0.91,0.91,0.91,0.91,0.91,0.91/)                              ; legend markers
  ;yleg = (/0.22,0.17,0.22,0.17,0.22,0.17,0.22,0.17/)                      ; and text
  yleg = (/0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.22,0.20,0.18,0.16,0.14,0.12,0.10/)
  ;ytxt = (/0.22,0.17,0.22,0.17,0.22,0.17,0.22,0.17/)                      ; strings.
  ytxt = (/0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.22,0.20,0.18,0.16,0.14,0.12,0.10/)

  do i = 0,npts-1 ;num_distinct_markers-1
      gsres@gsMarkerColor      = colors(i)
      gsres@gsMarkerThicknessF = 1.7 ;*(i+1)
      gsn_polymarker(wks,map,lon(i),lat(i),gsres)
;
; Add marker and text for the legend.
;
      gsn_polymarker_ndc(wks,          xleg(i),yleg(i),gsres)
      gsn_text_ndc      (wks,labels(i),xtxt(i),ytxt(i),txres)
  end do

  frame(wks)    ; Advance the frame.
end     
