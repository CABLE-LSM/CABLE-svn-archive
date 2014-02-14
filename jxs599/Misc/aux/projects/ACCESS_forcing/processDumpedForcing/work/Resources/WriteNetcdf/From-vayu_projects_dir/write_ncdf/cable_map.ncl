;*************************************************
;
;
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


;**************************************************
; Read data -  USERS
;**************************************************

    varfile = "tgg.nc"
    f  = addfile( varfile, "r" )

    var = f->tgg 
    lon = f->lon 
    lat = f->lat
    tile= f->tile
    time= f->time
    var@_FillValue = -2000000000
    
    ;yr = var 
    ;yr!0 = "time"
    ;yr!1 = "tile"
    ;yr!2 = "lat" 
    ;yr!3 = "lon" 
    ;yr&time= time
    ;yr&tile= tile
    ;yr&lat = lat
    ;yr&lon = lon
    
    ;avg_var = dim_avg_n(var,0)
    ;yr = avg_var 
    ;yr!0 = "tile"
    ;yr!1 = "lat" 
    ;yr!2 = "lon" 
    ;yr&tile= tile
    ;yr&lat = lat
    ;yr&lon = lon
   
    

    avg_var = dim_avg_n(var,0)
    yr = dim_sum_n(avg_var,0)
    yr!0 = "lat" 
    yr!1 = "lon" 
    yr&lat = lat
    yr&lon = lon
   
   
   
    
    
    yr&lat@units         = "degrees_north"
    yr&lon@units         = "degrees_east"
    yr@units         = "xxxx"
    yr@_FillValue = -2000000000

; Read user data into a one-dimensional int array of length 14:
    ;uservar = asciiread("uservar.txt",-l,"string")
    ;npts = dimsizes(data)   ; should be 14
    ;print(data)     ; Print the values

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
  
    res@cnMaximize       = True

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
    ;plot = gsn_csm_contour_map_ce(wks,yr(2,:,:),res)
    plot = gsn_csm_contour_map_ce(wks,yr(:,:),res)
    ;plot = gsn_csm_contour_map_ce(wks,yr(1,16,:,:),res)

end

