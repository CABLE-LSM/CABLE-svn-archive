#!/bin/csh

# Extract fields
ncks -v field1388 Mmonthly_means_${YR}yrs.nc mm.gpp_${YR}yrs.nc       # GPP
ncks -v field1389 Mmonthly_means_${YR}yrs.nc mm.npp_${YR}yrs.nc       # NPP
#ncks -v field1389_1 Mmonthly_means_${YR}yrs.nc mm.npp_${YR}yrs.nc       # NPP
#ncks -v $nname Mmonthly_means_${YR}yrs.nc mm.npp_${YR}yrs.nc       # NPP
ncks -v field1390 Mmonthly_means_${YR}yrs.nc mm.presp_${YR}yrs.nc     # Presp
ncks -v field1523 Mmonthly_means_${YR}yrs.nc mm.sresp_${YR}yrs.nc     # Sresp

ncks -v field1519 Mmonthly_means_${YR}yrs.nc mm.gppt_${YR}yrs.nc      # GPP on Tiles
ncks -v field1521 Mmonthly_means_${YR}yrs.nc mm.nppt_${YR}yrs.nc      # NPP on Tiles
ncks -v field1522 Mmonthly_means_${YR}yrs.nc mm.prespt_${YR}yrs.nc    # Presp on Tiles

#set varnames=`cdo showname Mmonthly_means_${YR}yrs.nc`
#setenv T15 n
#foreach name ( $varnames )
##if ($name == temp_15 && $T15 == n) then
#if ($tname == temp_15) then
#setenv T15 y
#endif
#if (temp_2) then 
# temporary removal - Les 20apr12 - due to xaanx having tscrn as temp_2
#ncks -v temp_2    Mmonthly_means_${YR}yrs.nc mm.srespt_${YR}yrs.nc    # Sresp on Tiles
#endif

ncks -v field1499 Mmonthly_means_${YR}yrs.nc mm.lresp_${YR}yrs.nc # frday (leaf turnover)
ncks -v field1563 Mmonthly_means_${YR}yrs.nc mm.co2atmos_${YR}yrs.nc  # CO2 totalflux to Atmos
ncks -v field1564 Mmonthly_means_${YR}yrs.nc mm.co2_${YR}yrs.nc       # CO2 3D Tracer Mass Mixing Ratio

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
# NHem j=90:145
# Trop j=57:89
# SHem j= 1:56

exit

