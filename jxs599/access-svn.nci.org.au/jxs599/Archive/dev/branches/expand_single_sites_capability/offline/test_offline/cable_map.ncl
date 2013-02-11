;*************************************************
;
;
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


;**************************************************
; Read data 
;**************************************************

    varfile = "ebal_tot_cncheck.nc"
    f  = addfile( varfile, "r" )

    var = f->ebal_tot_cncheck
    dt = f ->time 

    varfile = "ebal_tot.nc"
    
    j  = addfile( varfile, "r" )

    var2 = j->ebal_tot
 
;************************************************
; to plot multiple lines, you must put them into
; a mulidimensional array
;************************************************
  data      = new((/2,dimsizes(dt)/),float) 
  data(0,:) = var 
  data(1,:) = var2 

;************************************************
; plotting 
;************************************************
 wks   = gsn_open_wks ("pdf",varfile)                  ; open workstation

 res                   = True                      ; plot mods desired
 res@tiMainString      = "Mulitple XY plot"        ; add title

; note a similiar resource is xyLineThicknessF, which will effect all
; the lines in the array.
 res@xyLineThicknesses = (/1.0,2.0/)               ; make 2nd lines thicker

 res@xyLineColors      = (/"blue","red"/)          ; change line color

 plot  = gsn_csm_xy (wks,dt,data,res) ; create plot

end
