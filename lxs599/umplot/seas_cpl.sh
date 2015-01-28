#!/bin/csh

set a = a
@ noy = ($YR * 4)  # number of files for no. of years 
@ nom = ($YR * 12) # number of files for no. of years and months

# Monthly ------------------------------------------------------------

set pmlist=`ls $DIR/$RUNID.$Pmonth-??????????.nc | head -$nom` 
dmget $pmlist

set pmjan=`ls $DIR/$RUNID.$Pmonth-????001???.nc | head -$YR`
set pmfeb=`ls $DIR/$RUNID.$Pmonth-????002???.nc | head -$YR`
set pmmar=`ls $DIR/$RUNID.$Pmonth-????003???.nc | head -$YR`
set pmapr=`ls $DIR/$RUNID.$Pmonth-????004???.nc | head -$YR`
set pmmay=`ls $DIR/$RUNID.$Pmonth-????005???.nc | head -$YR`
set pmjun=`ls $DIR/$RUNID.$Pmonth-????006???.nc | head -$YR`
set pmjul=`ls $DIR/$RUNID.$Pmonth-????007???.nc | head -$YR`
set pmaug=`ls $DIR/$RUNID.$Pmonth-????008???.nc | head -$YR`
set pmsep=`ls $DIR/$RUNID.$Pmonth-????009???.nc | head -$YR`
set pmoct=`ls $DIR/$RUNID.$Pmonth-????010???.nc | head -$YR`
set pmnov=`ls $DIR/$RUNID.$Pmonth-????011???.nc | head -$YR`
set pmdec=`ls $DIR/$RUNID.$Pmonth-????012???.nc | head -$YR`

ncra -O $pmjan avejan.nc
ncra -O $pmfeb avefeb.nc
ncra -O $pmmar avemar.nc
ncra -O $pmapr aveapr.nc
ncra -O $pmmay avemay.nc
ncra -O $pmjun avejun.nc
ncra -O $pmjul avejul.nc
ncra -O $pmaug aveaug.nc
ncra -O $pmsep avesep.nc
ncra -O $pmoct aveoct.nc
ncra -O $pmnov avenov.nc
ncra -O $pmdec avedec.nc

ncrcat -O avejan.nc avefeb.nc avemar.nc aveapr.nc avemay.nc avejun.nc avejul.nc aveaug.nc avesep.nc aveoct.nc avenov.nc avedec.nc monthly_means_${YR}yrs.nc

# Seasonal -----------------------------------------------------------

ncra -O $pmdec $pmjan $pmfeb avedjf.nc
ncra -O $pmmar $pmapr $pmmay avemam.nc
ncra -O $pmjun $pmjul $pmaug avejja.nc
ncra -O $pmsep $pmoct $pmnov aveson.nc
ncrcat -O avedjf.nc avemam.nc avejja.nc aveson.nc seasonal_means_${YR}yrs.nc

#cdo mergetime $pmlist Mmonthly_means_${YR}yrs.nc

#if ( tas ) then
# $PLOT/dcc_all_scripts/remap_dcc_varnames.sh
#endif

# Timeseries ---------------------------------------------------------

dmget $DIR/$RUNID.$Ptimes-??????????.nc $DIR/$RUNID.$Ptemps-??????????.nc

# single points ----------
if ($REINIT == 1) then
 set pelist=`ls $DIR/$RUNID.$Ptimes-??????????.nc | head -$nom`
 set pflist=`ls $DIR/$RUNID.$Ptimes-??????????_*swlw*.nc | head -$nom`
 set pclist=`ls $DIR/$RUNID.pc-??????????_*swlw*.nc | head -$nom`
else
 set pelist=`ls $DIR/$RUNID.$Ptimes-??????????.nc | head -$noy`
 set pflist=`ls $DIR/$RUNID.$Ptimes-??????????_*swlw*.nc | head -$noy`
 set pclist=`ls $DIR/$RUNID.pc-??????????_*swlw*.nc | head -$noy`
endif
if ( ${#pelist} > 0 ) then
 cdo mergetime $pelist Timeseries_${YR}yrs.nc
 if ( $Ptimes == pf && ${#pflist} > 0 ) then
  cdo mergetime $DIR/$RUNID.$Ptimes-??????????\_noswlw.nc Timeseries_${YR}yrs_noswlw.nc
  cdo mergetime $DIR/$RUNID.$Ptimes-??????????\_swlw.nc Timeseries_${YR}yrs_swlw.nc
 endif
endif
if ( ${#pclist} > 0 ) then
  cdo mergetime $DIR/$RUNID.pc-??????????\_noswlw.nc PALS_ts_${YR}yrs_noswlw.nc
  cdo mergetime $DIR/$RUNID.pc-??????????\_swlw.nc PALS_ts_${YR}yrs_swlw.nc
endif
# single points ----------

# tmax and tmin ----------
if ($REINIT == 1) then
 set pblist=`ls $DIR/$RUNID.$Ptemps-??????????.nc | head -$nom`
else
 set pblist=`ls $DIR/$RUNID.$Ptemps-??????????.nc | head -$noy`
endif
if ( ${#pblist} > 0 ) then
 cdo mergetime $pblist Tempseries_${YR}yrs.nc
 cdo yseasmean Tempseries_${YR}yrs.nc Tseasonal_means_${YR}yrs.nc
 cdo ymonmean Tempseries_${YR}yrs.nc Tmonthly_means_${YR}yrs.nc
endif
# tmax and tmin ----------

rm ave*

exit

