#!/bin/csh

set a = a
@ noy = ($YR * 4)  # number of files for no. of years 
@ nom = ($YR * 12) # number of files for no. of years and months

# Monthly ------------------------------------------------------------

set pmlist=`/bin/ls $DIR/$RUNID.$Pmonth-??????????.nc | head -$nom` 
dmget $pmlist

set pmjan=`/bin/ls $DIR/$RUNID.$Pmonth-????001???.nc | head -$YR`
ncra -O $pmjan avejan.nc
set pmfeb=`/bin/ls $DIR/$RUNID.$Pmonth-????002???.nc | head -$YR`
ncra -O $pmfeb avefeb.nc
set pmmar=`/bin/ls $DIR/$RUNID.$Pmonth-????003???.nc | head -$YR`
ncra -O $pmmar avemar.nc
set pmapr=`/bin/ls $DIR/$RUNID.$Pmonth-????004???.nc | head -$YR`
ncra -O $pmapr aveapr.nc
set pmmay=`/bin/ls $DIR/$RUNID.$Pmonth-????005???.nc | head -$YR`
ncra -O $pmmay avemay.nc
set pmjun=`/bin/ls $DIR/$RUNID.$Pmonth-????006???.nc | head -$YR`
ncra -O $pmjun avejun.nc
set pmjul=`/bin/ls $DIR/$RUNID.$Pmonth-????007???.nc | head -$YR`
ncra -O $pmjul avejul.nc
set pmaug=`/bin/ls $DIR/$RUNID.$Pmonth-????008???.nc | head -$YR`
ncra -O $pmaug aveaug.nc
set pmsep=`/bin/ls $DIR/$RUNID.$Pmonth-????009???.nc | head -$YR`
ncra -O $pmsep avesep.nc
set pmoct=`/bin/ls $DIR/$RUNID.$Pmonth-????010???.nc | head -$YR`
ncra -O $pmoct aveoct.nc
set pmnov=`/bin/ls $DIR/$RUNID.$Pmonth-????011???.nc | head -$YR`
ncra -O $pmnov avenov.nc
set pmdec=`/bin/ls $DIR/$RUNID.$Pmonth-????012???.nc | head -$YR`
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
 set pelist=`/bin/ls $DIR/$RUNID.$Ptimes-??????????.nc | head -$nom`
else
 set pelist=`/bin/ls $DIR/$RUNID.$Ptimes-??????????.nc | head -$noy`
endif
if ( ${#pelist} > 0 ) then
 cdo mergetime $pelist Timeseries_${YR}yrs.nc

 if ( $Ptimes == pf ) then
  set month=( 001 002 003 004 005 006 007 008 009 010 011 012 )
  foreach mon ( $month )
   cdo mergetime $DIR/$RUNID.$Ptimes-????$mon???\_noswlw.nc $RUNID.$Ptimes\_$mon.nc
   cdo mergetime $DIR/$RUNID.$Ptimes-????$mon???\_swlw.nc $RUNID.$Ptimes\_$mon\_2.nc
  end
 endif
endif
# single points ----------

# tmax and tmin ----------
if ($REINIT == 1) then
 set pblist=`/bin/ls $DIR/$RUNID.$Ptemps-??????????.nc | head -$nom`
else
 set pblist=`/bin/ls $DIR/$RUNID.$Ptemps-??????????.nc | head -$noy`
endif
if ( ${#pblist} > 0 ) then
 cdo mergetime $pblist Tempseries_${YR}yrs.nc
 cdo yseasmean Tempseries_${YR}yrs.nc Tseasonal_means_${YR}yrs.nc
 cdo ymonmean Tempseries_${YR}yrs.nc Tmonthly_means_${YR}yrs.nc
endif
# tmax and tmin ----------

rm ave*

exit

