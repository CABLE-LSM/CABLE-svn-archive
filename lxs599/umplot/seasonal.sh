#!/bin/csh

set a = a
@ noy = ($YR * 4)  # number of files for no. of years 
@ nom = ($YR * 12) # number of files for no. of years and months
#@ djf = ($YR - 1)  # minus one year

# Monthly ------------------------------------------------------------

if ($RUNID != xagpd) then

set pmlist=`ls $DIR/$RUNID$a.pm*.nc | head -$nom` 
dmget $pmlist
set pmjan=`ls $DIR/$RUNID$a.pm*jan.nc | head -$YR`
ncra -O $pmjan avejan.nc
set pmfeb=`ls $DIR/$RUNID$a.pm*feb.nc | head -$YR`
ncra -O $pmfeb avefeb.nc
set pmmar=`ls $DIR/$RUNID$a.pm*mar.nc | head -$YR`
ncra -O $pmmar avemar.nc
set pmapr=`ls $DIR/$RUNID$a.pm*apr.nc | head -$YR`
ncra -O $pmapr aveapr.nc
set pmmay=`ls $DIR/$RUNID$a.pm*may.nc | head -$YR`
ncra -O $pmmay avemay.nc
set pmjun=`ls $DIR/$RUNID$a.pm*jun.nc | head -$YR`
ncra -O $pmjun avejun.nc
set pmjul=`ls $DIR/$RUNID$a.pm*jul.nc | head -$YR`
ncra -O $pmjul avejul.nc
set pmaug=`ls $DIR/$RUNID$a.pm*aug.nc | head -$YR`
ncra -O $pmaug aveaug.nc
set pmsep=`ls $DIR/$RUNID$a.pm*sep.nc | head -$YR`
ncra -O $pmsep avesep.nc
set pmoct=`ls $DIR/$RUNID$a.pm*oct.nc | head -$YR`
ncra -O $pmoct aveoct.nc
set pmnov=`ls $DIR/$RUNID$a.pm*nov.nc | head -$YR`
ncra -O $pmnov avenov.nc
set pmdec=`ls $DIR/$RUNID$a.pm*dec.nc | head -$YR`
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

endif

# Timeseries ---------------------------------------------------------

dmget $DIR/$RUNID$a.pe*0.nc $DIR/$RUNID$a.pb*0*
#module load hdf5/1.8.1 netcdf/4.0 cdo/1.4.0

#~ste69f/umplot/pe2nc.sh
if ($REINIT == 1) then
set pelist=`ls $DIR/$RUNID$a.pe*0.nc | head -$nom`
else
set pelist=`ls $DIR/$RUNID$a.pe*0.nc | head -$noy`
endif
cdo copy $pelist Timeseries_${YR}yrs.nc

# tmax and tmin
if ($REINIT == 1) then
set pblist=`ls $DIR/$RUNID$a.pb*0.nc | head -$nom`
else
set pblist=`ls $DIR/$RUNID$a.pb*0.nc | head -$noy`
endif
cdo copy $pblist Tempseries_${YR}yrs.nc

#dmget $DIR/$RUNID$a.pb*
#set pbum=`ls $DIR/$RUNID$a.pb*0 | head -$noy`

#foreach file ($pbum)
#if (-e $file.sub.nc) then
#python ~dix043/src/python/um/um_fields_subset.py -i $file -o $file.sub -v 3236
#~ste69f/umutils/conv2nc.tcl -i $file.sub -o $file.sub.nc
#cdo monavg $file.sub.nc $file.avg.nc
#endif
#end

#if (! -e $DIR/$RUNID$a.pb*0.avg.nc) then
cdo yseasavg Tempseries_${YR}yrs.nc Tseasonal_means_${YR}yrs.nc
cdo ymonavg Tempseries_${YR}yrs.nc Tmonthly_means_${YR}yrs.nc
#else
#set pbsub=`ls $DIR/$RUNID$a.pb*0.avg.nc | head -$noy`
#cdo copy $pbsub Trange_${YR}yrs.nc
#cdo yseasavg Trange_${YR}yrs.nc Tseasonal_means_${YR}yrs.nc
#cdo ymonavg Trange_${YR}yrs.nc Tmonthly_means_${YR}yrs.nc
#endif

# Runoff and SMC and FrozFrac
if ($RUNID == xaiya) then
set pasub=`ls $DIR/$RUNID$a.pa*0.avg.nc | head -$noy`
cdo copy $pasub TotRunoff_${YR}yrs.nc
cdo yseasavg TotRunoff_${YR}yrs.nc Rseasonal_means_${YR}yrs.nc
cdo ymonavg TotRunoff_${YR}yrs.nc Rmonthly_means_${YR}yrs.nc
ncrename -v sm,sm_1 seasonal_${YR}yrs.nc
ncrename -v sm,sm_1 monthly_${YR}yrs.nc
cdo merge seasonal_${YR}yrs.nc Rseasonal_means_${YR}yrs.nc seasonal_means_${YR}yrs.nc 
cdo merge monthly_${YR}yrs.nc Rmonthly_means_${YR}yrs.nc monthly_means_${YR}yrs.nc
else
mv seasonal_${YR}yrs.nc seasonal_means_${YR}yrs.nc
mv monthly_${YR}yrs.nc monthly_means_${YR}yrs.nc
endif

#rm *.sub *.sub.nc Trange_${YR}yrs.nc *.avg.nc
rm ave*

exit

