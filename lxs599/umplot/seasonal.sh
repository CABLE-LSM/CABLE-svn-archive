#!/bin/csh

set a = a
@ noy = ($YR * 4)  # number of files for no. of years 
@ nom = ($YR * 12) # number of files for no. of years and months
#@ djf = ($YR - 1)  # minus one year

# Monthly ------------------------------------------------------------

set pmlist=`ls $DIR/$RUNID$a.$Pmonth?????.nc | head -$nom` 
dmget $pmlist

set pmjan=`ls $DIR/$RUNID$a.$Pmonth??001.nc | head -$YR`
ncra -O $pmjan avejan.nc
set pmfeb=`ls $DIR/$RUNID$a.$Pmonth??002.nc | head -$YR`
ncra -O $pmfeb avefeb.nc
set pmmar=`ls $DIR/$RUNID$a.$Pmonth??003.nc | head -$YR`
ncra -O $pmmar avemar.nc
set pmapr=`ls $DIR/$RUNID$a.$Pmonth??004.nc | head -$YR`
ncra -O $pmapr aveapr.nc
set pmmay=`ls $DIR/$RUNID$a.$Pmonth??005.nc | head -$YR`
ncra -O $pmmay avemay.nc
set pmjun=`ls $DIR/$RUNID$a.$Pmonth??006.nc | head -$YR`
ncra -O $pmjun avejun.nc
set pmjul=`ls $DIR/$RUNID$a.$Pmonth??007.nc | head -$YR`
ncra -O $pmjul avejul.nc
set pmaug=`ls $DIR/$RUNID$a.$Pmonth??008.nc | head -$YR`
ncra -O $pmaug aveaug.nc
set pmsep=`ls $DIR/$RUNID$a.$Pmonth??009.nc | head -$YR`
ncra -O $pmsep avesep.nc
set pmoct=`ls $DIR/$RUNID$a.$Pmonth??010.nc | head -$YR`
ncra -O $pmoct aveoct.nc
set pmnov=`ls $DIR/$RUNID$a.$Pmonth??011.nc | head -$YR`
ncra -O $pmnov avenov.nc
set pmdec=`ls $DIR/$RUNID$a.$Pmonth??012.nc | head -$YR`
ncra -O $pmdec avedec.nc

#ncrcat -O avejan.nc avejul.nc mnmean_${YR}yrs.nc
ncrcat -O avejan.nc avefeb.nc avemar.nc aveapr.nc avemay.nc avejun.nc avejul.nc aveaug.nc avesep.nc aveoct.nc avenov.nc avedec.nc monthly_${YR}yrs.nc
#ncrcat -O avejan.nc avefeb.nc avemar.nc aveapr.nc avemay.nc avejun.nc avejul.nc aveaug.nc avesep.nc aveoct.nc avenov.nc avedec.nc monthly_means_${YR}yrs.nc

# Seasonal -----------------------------------------------------------

ncra -O $pmdec $pmjan $pmfeb avedjf.nc
ncra -O $pmmar $pmapr $pmmay avemam.nc
ncra -O $pmjun $pmjul $pmaug avejja.nc
ncra -O $pmsep $pmoct $pmnov aveson.nc
#ncrcat -O aveson.nc avedjf.nc avemam.nc avejja.nc season_${YR}yrs.nc
ncrcat -O avedjf.nc avemam.nc avejja.nc aveson.nc seasonal_${YR}yrs.nc
#ncrcat -O avedjf.nc avemam.nc avejja.nc aveson.nc seasonal_means_${YR}yrs.nc

# Timeseries ---------------------------------------------------------

dmget $DIR/$RUNID$a.$Ptimes?????.nc $DIR/$RUNID$a.$Ptemps?????.nc
#module load hdf5/1.8.1 netcdf/4.0 cdo/1.4.0

# single points ----------
#~ste69f/umplot/pe2nc.sh
if ($REINIT == 1) then
 set pelist=`ls $DIR/$RUNID$a.$Ptimes?????.nc | head -$nom`
else
 set pelist=`ls $DIR/$RUNID$a.$Ptimes?????.nc | head -$noy`
endif
cdo mergetime $pelist Timeseries_${YR}yrs.nc
#cdo copy $pelist Timeseries_${YR}yrs.nc
# single points ----------

# tmax and tmin ----------
if ($REINIT == 1) then
 set pblist=`ls $DIR/$RUNID$a.$Ptemps?????.nc | head -$nom`
else
 set pblist=`ls $DIR/$RUNID$a.$Ptemps?????.nc | head -$noy`
endif
cdo mergetime $pblist Tempseries_${YR}yrs.nc
#cdo copy $pblist Tempseries_${YR}yrs.nc
# tmax and tmin ----------

 cdo yseasavg Tempseries_${YR}yrs.nc Tseasonal_means_${YR}yrs.nc
 cdo ymonavg Tempseries_${YR}yrs.nc Tmonthly_means_${YR}yrs.nc

 mv seasonal_${YR}yrs.nc seasonal_means_${YR}yrs.nc
 mv monthly_${YR}yrs.nc monthly_means_${YR}yrs.nc

#rm *.sub *.sub.nc Trange_${YR}yrs.nc *.avg.nc
rm ave*

exit

