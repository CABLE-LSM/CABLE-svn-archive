#!/bin/csh

set a=a
set date=`date`
set DIRW=$PWD
#module load hdf5/1.8.1 netcdf/4.0 cdo/1.4.0
#dmget $DIR/$RUNID$a.ps*.nc 
dmget $DIR/$RUNID$a.pa*.nc
dmget $DIR/$RUNID$a.pb*.nc 
dmget $DIR/$RUNID$a.pe*.nc 
dmget $DIR/$RUNID$a.pm*.nc

#if ($YR == 20) then
# setenv block '1 2 3 4'
#else if ($YR == 15) then
# setenv block '1 2 3'
#else if ($YR == 10) then
# setenv block '1 2'
#endif

#foreach blk ($block)
#setenv BLOCK $blk

@ noy = ($FullYrs * 4)
@ nom = ($FullYrs * 12)
@ Numyr= $FullYrs + 1
set count_all=`ls $DIRW/h[0123456789]*.nc $DIRW/i[0123456789]*.nc $DIRW/j[0123456789]*.nc $DIRW/k[0123456789]*.nc | wc -l`
@ countm1 = $count_all - 1
@ countm2 = $count_all - 2

if (! -d $DIRW/block${BLOCK}_5yrs) then
 mkdir $DIRW/block${BLOCK}_5yrs
 cd $DIRW/block${BLOCK}_5yrs

 set pelist=`ls $DIR/$RUNID$a.pe*0.nc | head -${noy}`
 set pblist=`ls $DIR/$RUNID$a.pb*0.nc | head -${noy}`
 if ($RUNID == xaiya) then
 set palist=`ls $DIR/$RUNID$a.pa*0.avg.nc | head -${noy}`
 endif
 if ($count_all == $Numyr) then
 echo "(1)     Warning: Number of Years = YR+1. Please Check."
 set monlist=`ls $DIRW/h[0123456789]*.nc $DIRW/i[0123456789]*.nc $DIRW/j[0123456789]*.nc $DIRW/k[0123456789]*.nc | head -${countm1} | tail -${countm2}` ### head -$NumYrs ???
 endif

if ($BLOCK == 1) then
 set pmjan=`ls $DIR/$RUNID$a.pm*jan.nc | head -5`
 set pmfeb=`ls $DIR/$RUNID$a.pm*feb.nc | head -5`
 set pmmar=`ls $DIR/$RUNID$a.pm*mar.nc | head -5`
 set pmapr=`ls $DIR/$RUNID$a.pm*apr.nc | head -5`
 set pmmay=`ls $DIR/$RUNID$a.pm*may.nc | head -5`
 set pmjun=`ls $DIR/$RUNID$a.pm*jun.nc | head -5`
 set pmjul=`ls $DIR/$RUNID$a.pm*jul.nc | head -5`
 set pmaug=`ls $DIR/$RUNID$a.pm*aug.nc | head -5`
 set pmsep=`ls $DIR/$RUNID$a.pm*sep.nc | head -5`
 set pmoct=`ls $DIR/$RUNID$a.pm*oct.nc | head -5`
 set pmnov=`ls $DIR/$RUNID$a.pm*nov.nc | head -5`
 set pmdec=`ls $DIR/$RUNID$a.pm*dec.nc | head -5`

 set pels=`ls $pelist | head -20`
if ($RUNID == uabab) then
 set pbls=`ls $pblist | head -60`
else
 set pbls=`ls $pblist | head -20`
endif
 if ($RUNID == xaiya) then
 set pals=`ls $palist | head -20`
 endif

 if ($count_all == $Numyr) then
 set cdo_blk=`ls $monlist | head -5`
 else
 set cdo_blk=`ls $DIRW/h[0123456789]*.nc $DIRW/i[0123456789]*.nc $DIRW/j[0123456789]*.nc $DIRW/k[0123456789]*.nc | head -5`
 endif

#if (-e $DIR/$RUNID$a.pb*0.avg.nc) then
#set pbsub=`ls $DIR/$RUNID$a.pb*0.avg.nc | head -20`
#endif

else if ($BLOCK == 2) then
 set pmjan=`ls $DIR/$RUNID$a.pm*jan.nc | head -10 | tail -5` # tail +6 | head -5
 set pmfeb=`ls $DIR/$RUNID$a.pm*feb.nc | head -10 | tail -5`
 set pmmar=`ls $DIR/$RUNID$a.pm*mar.nc | head -10 | tail -5`
 set pmapr=`ls $DIR/$RUNID$a.pm*apr.nc | head -10 | tail -5`
 set pmmay=`ls $DIR/$RUNID$a.pm*may.nc | head -10 | tail -5`
 set pmjun=`ls $DIR/$RUNID$a.pm*jun.nc | head -10 | tail -5`
 set pmjul=`ls $DIR/$RUNID$a.pm*jul.nc | head -10 | tail -5`
 set pmaug=`ls $DIR/$RUNID$a.pm*aug.nc | head -10 | tail -5`
 set pmsep=`ls $DIR/$RUNID$a.pm*sep.nc | head -10 | tail -5`
 set pmoct=`ls $DIR/$RUNID$a.pm*oct.nc | head -10 | tail -5`
 set pmnov=`ls $DIR/$RUNID$a.pm*nov.nc | head -10 | tail -5`
 set pmdec=`ls $DIR/$RUNID$a.pm*dec.nc | head -10 | tail -5`

 set pels=`ls $pelist | head -40 | tail -20` #tail +21 | head -20
if ($RUNID == uabab) then
 set pbls=`ls $pblist | head -120 | tail -60`
else
 set pbls=`ls $pblist | head -40 | tail -20`
endif
 if ($RUNID == xaiya) then
 set pals=`ls $palist | head -40 | tail -20`
 endif

 if ($count_all == $Numyr) then
 set cdo_blk=`ls $monlist | head -10 | tail -5`
 else
 set cdo_blk=`ls $DIRW/h[0123456789]*.nc $DIRW/i[0123456789]*.nc $DIRW/j[0123456789]*.nc $DIRW/k[0123456789]*.nc | head -10 | tail -5`
 endif

#if (-e $DIR/$RUNID$a.pb*0.avg.nc) then
#set pbsub=`ls $DIR/$RUNID$a.pb*0.avg.nc | tail +21 | head -20`
#endif

else if ($BLOCK == 3) then
 set pmjan=`ls $DIR/$RUNID$a.pm*jan.nc | head -15 | tail -5` # tail +11 | head -5
 set pmfeb=`ls $DIR/$RUNID$a.pm*feb.nc | head -15 | tail -5`
 set pmmar=`ls $DIR/$RUNID$a.pm*mar.nc | head -15 | tail -5`
 set pmapr=`ls $DIR/$RUNID$a.pm*apr.nc | head -15 | tail -5`
 set pmmay=`ls $DIR/$RUNID$a.pm*may.nc | head -15 | tail -5`
 set pmjun=`ls $DIR/$RUNID$a.pm*jun.nc | head -15 | tail -5`
 set pmjul=`ls $DIR/$RUNID$a.pm*jul.nc | head -15 | tail -5`
 set pmaug=`ls $DIR/$RUNID$a.pm*aug.nc | head -15 | tail -5`
 set pmsep=`ls $DIR/$RUNID$a.pm*sep.nc | head -15 | tail -5`
 set pmoct=`ls $DIR/$RUNID$a.pm*oct.nc | head -15 | tail -5`
 set pmnov=`ls $DIR/$RUNID$a.pm*nov.nc | head -15 | tail -5`
 set pmdec=`ls $DIR/$RUNID$a.pm*dec.nc | head -15 | tail -5`

 set pels=`ls $pelist | head -60 | tail -20` # tail +41 | head -20
if ($RUNID == uabab) then
 set pbls=`ls $pblist | head -180 | tail -60`
else
 set pbls=`ls $pblist | head -60 | tail -20`
endif
 if ($RUNID == xaiya) then
 set pals=`ls $palist | head -60 | tail -20`
 endif

 if ($count_all == $Numyr) then
 set cdo_blk=`ls $monlist | head -15 | tail -5`
 else
 set cdo_blk=`ls $DIRW/h[0123456789]*.nc $DIRW/i[0123456789]*.nc $DIRW/j[0123456789]*.nc $DIRW/k[0123456789]*.nc | head -15 | tail -5`
 endif

#if (-e $DIR/$RUNID$a.pb*0.avg.nc) then
#set pbsub=`ls $DIR/$RUNID$a.pb*0.avg.nc | tail +41 | head -20`
#endif

else if ($BLOCK == 4) then
 set pmjan=`ls $DIR/$RUNID$a.pm*jan.nc | head -20 | tail -5`
 set pmfeb=`ls $DIR/$RUNID$a.pm*feb.nc | head -20 | tail -5`
 set pmmar=`ls $DIR/$RUNID$a.pm*mar.nc | head -20 | tail -5`
 set pmapr=`ls $DIR/$RUNID$a.pm*apr.nc | head -20 | tail -5`
 set pmmay=`ls $DIR/$RUNID$a.pm*may.nc | head -20 | tail -5`
 set pmjun=`ls $DIR/$RUNID$a.pm*jun.nc | head -20 | tail -5`
 set pmjul=`ls $DIR/$RUNID$a.pm*jul.nc | head -20 | tail -5`
 set pmaug=`ls $DIR/$RUNID$a.pm*aug.nc | head -20 | tail -5`
 set pmsep=`ls $DIR/$RUNID$a.pm*sep.nc | head -20 | tail -5`
 set pmoct=`ls $DIR/$RUNID$a.pm*oct.nc | head -20 | tail -5`
 set pmnov=`ls $DIR/$RUNID$a.pm*nov.nc | head -20 | tail -5`
 set pmdec=`ls $DIR/$RUNID$a.pm*dec.nc | head -20 | tail -5`

 set pels=`ls $pelist | tail -20`
if ($RUNID == uabab) then
 set pbls=`ls $pblist | tail -60`
else
 set pbls=`ls $pblist | tail -20`
endif
 if ($RUNID == xaiya) then
 set pals=`ls $palist | tail -20`
 endif

 if ($count_all == $Numyr) then
 set cdo_blk=`ls $monlist | tail -5`
 else
 set cdo_blk=`ls $DIRW/h[0123456789]*.nc $DIRW/i[0123456789]*.nc $DIRW/j[0123456789]*.nc $DIRW/k[0123456789]*.nc | head -20 | tail -5`
 endif

#if (-e $DIR/$RUNID$a.pb*0.avg.nc) then
#set pbsub=`ls $DIR/$RUNID$a.pb*0.avg.nc | tail -20`
#endif

endif

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
ncra -O $pmsep $pmoct $pmnov aveson.nc
ncra -O $pmdec $pmjan $pmfeb avedjf.nc
ncra -O $pmmar $pmapr $pmmay avemam.nc
ncra -O $pmjun $pmjul $pmaug avejja.nc
ncrcat -O avejan.nc avefeb.nc avemar.nc aveapr.nc avemay.nc avejun.nc avejul.nc aveaug.nc avesep.nc aveoct.nc avenov.nc avedec.nc monthly_5yrs.nc
#ncrcat -O avejan.nc avefeb.nc avemar.nc aveapr.nc avemay.nc avejun.nc avejul.nc aveaug.nc avesep.nc aveoct.nc avenov.nc avedec.nc monthly_means_5yrs.nc
ncrcat -O avedjf.nc avemam.nc avejja.nc aveson.nc seasonal_5yrs.nc
#ncrcat -O avedjf.nc avemam.nc avejja.nc aveson.nc seasonal_means_5yrs.nc

cdo copy $pels Timeseries_5yrs.nc
cdo copy $pbls Tempseries_5yrs.nc
if ($RUNID == xaiya) then
cdo copy $pals Runfseries_5yrs.nc
endif

cdo copy $cdo_blk Mmonthly_means_5yrs.nc
cdo yearavg Mmonthly_means_5yrs.nc yearly_means_5yrs.nc

#if (-e $DIR/$RUNID$a.pb*0.avg.nc) then
#cdo copy $pbsub Trange_5yrs.nc
#cdo yseasavg Trange_5yrs.nc Tseasonal_means_5yrs.nc
#cdo ymonavg Trange_5yrs.nc Tmonthly_means_5yrs.nc
#else
cdo yseasavg Tempseries_5yrs.nc Tseasonal_means_5yrs.nc
cdo ymonavg Tempseries_5yrs.nc Tmonthly_means_5yrs.nc
#endif

if ($RUNID == xaiya) then
cdo yseasavg Runfseries_5yrs.nc Rseasonal_means_5yrs.nc
cdo ymonavg Runfseries_5yrs.nc Rmonthly_means_5yrs.nc
ncrename -v sm,sm_1 seasonal_5yrs.nc
ncrename -v sm,sm_1 monthly_5yrs.nc
cdo merge seasonal_5yrs.nc Rseasonal_means_5yrs.nc seasonal_means_5yrs.nc
cdo merge monthly_5yrs.nc Rmonthly_means_5yrs.nc monthly_means_5yrs.nc
else
mv seasonal_5yrs.nc seasonal_means_5yrs.nc
mv monthly_5yrs.nc monthly_means_5yrs.nc
endif

# cd $DIRW

else
 echo "(1)     $date - Directory already Exists"
endif

rm ave*.nc

#end

exit

