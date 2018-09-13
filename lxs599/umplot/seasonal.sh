#!/bin/csh

set a = a
@ noy = ($YR * 4)  # number of files for no. of years 
@ nom = ($YR * 12) # number of files for no. of years and months
@ nod = ($YR * 12 * 31) # number of files for no. of years, months & days

# Monthly ------------------------------------------------------------

set pmlist=`ls $DIR/$RUNID$a.$Pmonth?????.nc | head -$nom` 
#dmget $pmlist

set pmjan=`ls $DIR/$RUNID$a.$Pmonth??001.nc | head -$YR`
set pmfeb=`ls $DIR/$RUNID$a.$Pmonth??002.nc | head -$YR`
set pmmar=`ls $DIR/$RUNID$a.$Pmonth??003.nc | head -$YR`
set pmapr=`ls $DIR/$RUNID$a.$Pmonth??004.nc | head -$YR`
set pmmay=`ls $DIR/$RUNID$a.$Pmonth??005.nc | head -$YR`
set pmjun=`ls $DIR/$RUNID$a.$Pmonth??006.nc | head -$YR`
set pmjul=`ls $DIR/$RUNID$a.$Pmonth??007.nc | head -$YR`
set pmaug=`ls $DIR/$RUNID$a.$Pmonth??008.nc | head -$YR`
set pmsep=`ls $DIR/$RUNID$a.$Pmonth??009.nc | head -$YR`
set pmoct=`ls $DIR/$RUNID$a.$Pmonth??010.nc | head -$YR`
set pmnov=`ls $DIR/$RUNID$a.$Pmonth??011.nc | head -$YR`
set pmdec=`ls $DIR/$RUNID$a.$Pmonth??012.nc | head -$YR`

#if (${#pmjan} == 0) then
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
#else
#exit()
#endif

#ncrcat -O avejan.nc avefeb.nc avemar.nc aveapr.nc avemay.nc avejun.nc avejul.nc aveaug.nc avesep.nc aveoct.nc avenov.nc avedec.nc monthly_${YR}yrs.nc
ncrcat -O avejan.nc avefeb.nc avemar.nc aveapr.nc avemay.nc avejun.nc avejul.nc aveaug.nc avesep.nc aveoct.nc avenov.nc avedec.nc monthly_means_${YR}yrs.nc

# Seasonal -----------------------------------------------------------

ncra -O $pmdec $pmjan $pmfeb avedjf.nc
ncra -O $pmmar $pmapr $pmmay avemam.nc
ncra -O $pmjun $pmjul $pmaug avejja.nc
ncra -O $pmsep $pmoct $pmnov aveson.nc
#ncrcat -O aveson.nc avedjf.nc avemam.nc avejja.nc season_${YR}yrs.nc
#ncrcat -O avedjf.nc avemam.nc avejja.nc aveson.nc seasonal_${YR}yrs.nc
ncrcat -O avedjf.nc avemam.nc avejja.nc aveson.nc seasonal_means_${YR}yrs.nc

# Timeseries ---------------------------------------------------------

#dmget $DIR/$RUNID$a.$Ptimes?????.nc $DIR/$RUNID$a.$Ptemps?????.nc

# single points ----------
if ($REINIT == 1) then
 if ( $TSTEP == 288 ) then
  set pelist=`ls $DIR/$RUNID$a.$Ptimes?????.nc | head -$nod`
  set pflist=`ls $DIR/$RUNID$a.$Ptimes?????_*swlw*.nc | head -$nod`
 else
  set pelist=`ls $DIR/$RUNID$a.$Ptimes?????.nc | head -$nom`
  set pflist=`ls $DIR/$RUNID$a.$Ptimes?????_*swlw*.nc | head -$nom`
 endif
 set pclist=`ls $DIR/$RUNID$a.pc?????_*swlw*.nc | head -$nom`
else
 set pelist=`ls $DIR/$RUNID$a.$Ptimes?????.nc | head -$noy`
 set pflist=`ls $DIR/$RUNID$a.$Ptimes?????_*swlw*.nc | head -$noy`
 set pclist=`ls $DIR/$RUNID$a.pc?????_*swlw*.nc | head -$noy`
endif
if ( ${#pelist} > 0 ) then
 cdo mergetime $pelist Timeseries_${YR}yrs.nc
 #if ( $Ptimes == pf && ${#pflist} > 0 ) then
  cdo mergetime $DIR/$RUNID$a.$Ptimes?????_noswlw.nc Timeseries_${YR}yrs_noswlw.nc
  cdo mergetime $DIR/$RUNID$a.$Ptimes?????_swlw.nc Timeseries_${YR}yrs_swlw.nc
 #endif
set hname=`echo $HOST | head -c4`
if ($hname == chip || $hname == raij)  then
 ncrename -d t,time Timeseries_${YR}yrs.nc 
 ncrename -v t,time Timeseries_${YR}yrs.nc
 ncrename -d t,time Timeseries_${YR}yrs_swlw.nc 
 ncrename -v t,time Timeseries_${YR}yrs_swlw.nc
 ncrename -d t,time Timeseries_${YR}yrs_noswlw.nc 
 ncrename -v t,time Timeseries_${YR}yrs_noswlw.nc
endif
endif
if ( ${#pclist} > 0 ) then
  cdo mergetime $DIR/$RUNID$a.pc?????_noswlw.nc ts_pc_${YR}yrs_noswlw.nc
  cdo mergetime $DIR/$RUNID$a.pc?????_swlw.nc ts_pc_${YR}yrs_swlw.nc
endif
# single points ----------

# tmax and tmin ----------
if ($REINIT == 1) then
 set pblist=`ls $DIR/$RUNID$a.$Ptemps?????.nc | head -$nom`
else
 set pblist=`ls $DIR/$RUNID$a.$Ptemps?????.nc | head -$noy`
endif
if ( ${#pblist} > 0 ) then
 cdo mergetime $pblist Tempseries_${YR}yrs.nc
 cdo yseasmean Tempseries_${YR}yrs.nc Tseasonal_means_${YR}yrs.nc
 cdo ymonmean Tempseries_${YR}yrs.nc Tmonthly_means_${YR}yrs.nc
endif
# tmax and tmin ----------

#mv seasonal_${YR}yrs.nc seasonal_means_${YR}yrs.nc
#mv monthly_${YR}yrs.nc monthly_means_${YR}yrs.nc

rm ave*

exit

