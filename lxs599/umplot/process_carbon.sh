#!/bin/csh

set npp0=False
set npp1=False
set co2atm=False

# === Extract fields ===
# =======================================================================
set varn=`cdo showname Mmonthly_means_${YR}yrs.nc`
foreach nvr ( $varn )
 if ($nvr == field1389) then
  set npp0=True
  set npp1=False
 endif
 if ($nvr == field1389_1) then
  set npp0=False
  set npp1=True
 endif
 if ($nvr == field1563) then
  set co2atm=True
 endif
end
if ( $npp0 == True ) then
 ncks -v field1389 Mmonthly_means_${YR}yrs.nc mm.npp_${YR}yrs.nc       # NPP
endif
if ( $npp1 == True ) then
 ncks -v field1389_1 Mmonthly_means_${YR}yrs.nc mm.npp_${YR}yrs.nc     # NPP
endif
#ncks -v $nname Mmonthly_means_${YR}yrs.nc mm.npp_${YR}yrs.nc           # NPP

ncks -v field1388 Mmonthly_means_${YR}yrs.nc mm.gpp_${YR}yrs.nc        # GPP
ncks -v field1390 Mmonthly_means_${YR}yrs.nc mm.presp_${YR}yrs.nc      # Presp
ncks -v field1523 Mmonthly_means_${YR}yrs.nc mm.sresp_${YR}yrs.nc      # Sresp
ncks -v field1519 Mmonthly_means_${YR}yrs.nc mm.gppt_${YR}yrs.nc       # GPP on Tiles
ncks -v field1521 Mmonthly_means_${YR}yrs.nc mm.nppt_${YR}yrs.nc       # NPP on Tiles
ncks -v field1522 Mmonthly_means_${YR}yrs.nc mm.prespt_${YR}yrs.nc     # Presp on Tiles
ncks -v field1499 Mmonthly_means_${YR}yrs.nc mm.lresp_${YR}yrs.nc      # frday (leaf turnover)
if ( $co2atm == True ) then
 ncks -v field1563 Mmonthly_means_${YR}yrs.nc mm.co2atmos_${YR}yrs.nc  # CO2 totalflux to Atmos
endif
ncks -v field1564 Mmonthly_means_${YR}yrs.nc mm.co2_${YR}yrs.nc        # CO2 3D Tracer Mass Mixing Ratio

# =======================================================================

# Convert CO2 units to ppm
cdo mulc,659090.91 mm.co2_${YR}yrs.nc mm.co2ppm_${YR}yrs.nc

## Convert kgC/m2/s to gC/m2/yr (assume 365 days) = 1000*365*24*60*60
#cdo mulc,31536000000 mm.gpp.nc mm.gpp_gcpy.nc
#cdo mulc,31536000000 mm.npp.nc mm.npp_gcpy.nc
#cdo mulc,31536000000 mm.presp.nc mm.presp_gcpy.nc
#cdo mulc,31536000000 mm.sresp.nc mm.sresp_gcpy.nc

## Multiply(mul)/Divide(div) by land fraction
#cdo div /home/cmar/ste69f/umplot/landfrac_ACCESS_N96.nc mm.gpp_gpy.nc mm.gpp_lfr.nc
#cdo div /home/cmar/ste69f/umplot/landfrac_ACCESS_N96.nc mm.npp_gpy.nc mm.npp_lfr.nc
#cdo div /home/cmar/ste69f/umplot/landfrac_ACCESS_N96.nc mm.presp_gpy.nc mm.presp_lfr.nc
#cdo div /home/cmar/ste69f/umplot/landfrac_ACCESS_N96.nc mm.sresp_gpy.nc mm.sresp_lfr.nc

# =======================================================================

exit

