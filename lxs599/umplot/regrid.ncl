; Function Regrid_bilinear ======================================
; Lauren Stevens 29 April 2011
; Regirds a variable using linint2 (bilinear interpolation)
; ===============================================================

;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

function regrid_bilin(var:numeric,var_gd:numeric)

begin

  ;newlat = fspan(-90.,90.,145)
  ;print(newlat)
  ;newlon = fspan(0.,358.,193)
  ;print(newlon)
  ;newlat@units = "degrees_north"
  ;newlon@units = "degrees_east"

  dNames = getvardims(var_gd)
  nDims  = dimsizes(var)
  nSize  = dimsizes(nDims)
  lon    = var_gd&$dNames(3)$
  ;newvar = linint2(var&longitude,var&latitude,var,True,newlon,newlat,0)
  newvar = linint2(var&longitude,var&latitude,var,True,lon,var_gd&latitude,0)

  ;printVarSummary(newvar)

  ;copy_VarMeta(var_gd(:,0,:,:),var)
if (nSize.ge.3) then
  newvar!1         = "latitude"
  newvar!2         = "longitude"
else
  newvar!0         = "latitude"
  newvar!1         = "longitude"
end if
  newvar&latitude  = var_gd&latitude    ;newlat
  newvar&longitude = var_gd&$dNames(3)$ ;newlon

  return(newvar)

end

