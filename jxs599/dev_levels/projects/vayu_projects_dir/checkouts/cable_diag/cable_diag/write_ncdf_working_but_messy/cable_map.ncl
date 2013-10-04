;*************************************************
;Plot CABLE monthly output,example
;Q.Zhang 01/08/2011
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;**************************************************
; Read data
;**************************************************
    f  = addfile("nc.nc","r")

    var = f->var
    lon = f->lon 
    lat = f->lat
    var@_FillValue = -2000
    
    
    ;frac = f->patchfrac


    avg_var = dim_avg_n(var,0)
    ;yr     = dim_sum_n(avg_mn*frac,0)*12*365*86400*1.0E-6 ;convert umol/m2/s to gC/m2/yr
    ;yr     = dim_sum_n(avg_mn,0)*12*365*86400*1.0E-6 ;convert umol/m2/s to gC/m2/yr


    yr = dim_sum_n(avg_var,0)

    yr!0 = "lat" 
    yr!1 = "lon" 
    yr&lat = lat
    yr&lon = lon
    yr&lat@units         = "degrees_north"
    yr&lon@units         = "degrees_east"
    yr@units         = "xxxx"
    yr@_FillValue = -2000
    ;yr = mask(yr, yr.gt.-100,True)

;**************************************************
; plot
;**************************************************
    wks  = gsn_open_wks("x11","var")
    gsn_define_colormap(wks,"BlGrYeOrReVi200")  ; choose colormap

    res                      = True
    res@gsnSpreadColors      = True
    res@gsnSpreadColorEnd    = 180

    res@cnFillOn             = True
    res@cnFillMode           = "RasterFill"
    res@cnLinesOn            = False
  
;    res@cnLevelSelectionMode= "ManualLevels"   ; manual set levels
;    res@cnMaxLevelValF       = 3500
;    res@cnMinLevelValF       = 500
;    res@cnLevelSpacingF      = 500
  
    res@lbLabelStride        = 2
    res@lbRasterFillOn       = True

;    res@lbOrientation        = "Vertical"
;    res@lbLabelBarOn        = False            ; turn off individual label bars
     
     res@gsnAddCyclic = False 
    
    res@tiMainString        = "var"
    plot = gsn_csm_contour_map_ce(wks,yr,res)

end

