; Function Regrid_bilinear ======================================
; Lauren Stevens 29 April 2011
; Regirds a variable using linint2 (bilinear interpolation)
; ===============================================================

;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;undef("regrid_bilin")
function regrid_bilin(var:numeric,var_gd:numeric)
;local dNames,nDims,nSize,lon,lat,dNam,ndim,nsiz,lonv,latv,newvar

begin

  ;newlat = fspan(-90.,90.,145)
  ;print(newlat)
  ;newlon = fspan(0.,358.,193)
  ;print(newlon)
  ;newlat@units = "degrees_north"
  ;newlon@units = "degrees_east"

  dNames = getvardims(var_gd)
  nDims  = dimsizes(var_gd)
  nSize  = dimsizes(nDims)
  lon    = var_gd&$dNames(nSize-1)$
  lat    = var_gd&$dNames(nSize-2)$

  dNam   = getvardims(var)
  ndim   = dimsizes(var)
  nsiz   = dimsizes(ndim)
  lonv   = var&$dNam(nsiz-1)$
  latv   = var&$dNam(nsiz-2)$

  ;newvar = linint2(var&longitude,var&latitude,var,True,newlon,newlat,0)
  ;newvar = linint2(var&longitude,var&latitude,var,True,lon,var_gd&latitude,0)
  newvar = linint2(lonv,latv,var,True,lon,lat,0)

  ;printVarSummary(var)
  ;printVarSummary(var_gd)
  ;printVarSummary(newvar)
  ;copy_VarMeta(var_gd,var)

if (nsiz.eq.4) then
  newvar!2         = ""+dNam(2)+""
  newvar!3         = ""+dNam(3)+""
else
 if (nsiz.eq.3) then
  newvar!1         = ""+dNam(1)+""
  newvar!2         = ""+dNam(2)+""
 else
  newvar!0         = ""+dNam(2)+""
  newvar!1         = ""+dNam(3)+""
 end if
end if
  newvar&$dNam(nsiz-2)$ = var_gd&$dNames(nSize-2)$ ;newlat
  newvar&$dNam(nsiz-1)$ = var_gd&$dNames(nSize-1)$ ;newlon

  return(newvar)

end

