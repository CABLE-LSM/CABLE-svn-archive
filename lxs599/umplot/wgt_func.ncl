; -------------------------------------------------------
; Lauren Stevens 25 Feb 2013
; -------------------------------------------------------
; Get the area weights for global, land & ocean averaging
; -------------------------------------------------------
; load --------------------------------------------------
;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; -------------------------------------------------------

;undef("get_wgt")
function get_wgt(var)

begin

dNames = getvardims(var)
dims   = dimsizes(var)
nsiz   = dimsizes(dims)
nlat   = dims(nsiz-2)
mlon   = dims(nsiz-1)

;print(dNames)
lat    = var&$dNames(nsiz-2)$
lon    = var&$dNames(nsiz-1)$

rad    = 4.0*atan(1.0) / 180.
re     = 6.37122e06
con    = re * rad
clat   = cos(lat*rad)

dlon   = (lon(2)-lon(1))
dlat   = (lat(2)-lat(1))
dx     = con * dlon * clat
dy     = con * dlat
dydx   = dy * dx

wgt    = new((/nlat,mlon/),double)
wgt    = conform(wgt,dydx,0)

;--land frac o-ocean and l-land--
;wgto   = wgt*(1-lfrac(0,0,:,:))
;wgtl   = wgt*lfrac(0,0,:,:)

;--land mask o-ocean and l-land--
;wgtlm  = wgt*(lsm(0,0,:,:))
;wgtom  = wgt*(1-lsm(0,0,:,:))

return(wgt)

end

