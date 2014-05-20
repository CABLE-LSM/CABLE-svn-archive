#!/bin/csh

set a = a
@ noy = ($YR * 4)  # number of files for no. of years 
@ nom = ($YR * 12) # number of files for no. of years and months

# Monthly ------------------------------------------------------------

set pmlist=`ls $DIR/$RUNID$a.$Pmonth?????.nc | head -$nom` 
dmget $pmlist

set pmjan=`ls $DIR/$RUNID$a.$Pmonth??001.nc | head -$YR`
#if (${#pmjan} == 0) then
ncra -O $pmjan avejan.nc
#else
#exit()
#endif
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

dmget $DIR/$RUNID$a.$Ptimes?????.nc $DIR/$RUNID$a.$Ptemps?????.nc

# single points ----------
if ($REINIT == 1) then
 set pelist=`ls $DIR/$RUNID$a.$Ptimes?????.nc | head -$nom`
else
 set pelist=`ls $DIR/$RUNID$a.$Ptimes?????.nc | head -$noy`
endif
if ( ${#pelist} > 0 ) then
 cdo mergetime $pelist Timeseries_${YR}yrs.nc
 if ( $Ptimes == pf ) then
  cdo mergetime $DIR/$RUNID$a.$Ptimes?????_noswlw.nc Timeseries_${YR}yrs_noswlw.nc
  cdo mergetime $DIR/$RUNID$a.$Ptimes?????_swlw.nc Timeseries_${YR}yrs_swlw.nc
  set month=( 001 002 003 004 005 006 007 008 009 010 011 012 )
  foreach mon ( $month )
   cdo mergetime $DIR/$RUNID$a.$Ptimes??$mon\_noswlw.nc $RUNID$a.$Ptimes\_$mon\_${YR}yrs.nc
   cdo mergetime $DIR/$RUNID$a.$Ptimes??$mon\_swlw.nc $RUNID$a.$Ptimes\_$mon\_2_${YR}yrs.nc
  end
 endif
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

