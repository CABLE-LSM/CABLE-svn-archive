#!/bin/csh 

unalias ls
#alias ls ls

#----------------------------
set a=a
set date=`date`
set DIRW=$PWD
#----------------------------

#module load hdf5/1.8.1 netcdf/4.0 cdo/1.4.0
#dmget $DIR/$RUNID$a.ps?????.nc 
echo " "
dmget $DIR/$RUNID$a.$Pdaily?????.nc
dmget $DIR/$RUNID$a.$Ptemps?????.nc 
dmget $DIR/$RUNID$a.$Ptimes?????.nc 
dmget $DIR/$RUNID$a.$Pmonth?????.nc

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
set count_all=`ls $DIRW/h[0123456789].nc $DIRW/i[0123456789].nc $DIRW/j[0123456789].nc $DIRW/k[0123456789].nc | wc -l`
@ countm1 = $count_all - 1
@ countm2 = $count_all - 2

if (! -d $DIRW/block${BLOCK}_5yrs) then
 mkdir $DIRW/block${BLOCK}_5yrs
 cd $DIRW/block${BLOCK}_5yrs

if ($REINIT == 1) then
 set pelist=`ls $DIR/$RUNID$a.$Ptimes?????.nc | head -${nom}`
 set pblist=`ls $DIR/$RUNID$a.$Ptemps?????.nc | head -${nom}`
else
 set pelist=`ls $DIR/$RUNID$a.$Ptimes?????.nc | head -${noy}`
 set pblist=`ls $DIR/$RUNID$a.$Ptemps?????.nc | head -${noy}`
endif
 if ($count_all == $Numyr) then
  echo "(1)     Warning: Number of Years = YR+1. Please Check."
  set monlist=`ls $DIRW/h[0123456789].nc $DIRW/i[0123456789].nc $DIRW/j[0123456789].nc $DIRW/k[0123456789].nc | head -${countm1} | tail -${countm2}` ### head -$NumYrs ???
 endif

if ($BLOCK == 1) then
 set pmjan=`ls $DIR/$RUNID$a.$Pmonth??001.nc | head -5`
 set pmfeb=`ls $DIR/$RUNID$a.$Pmonth??002.nc | head -5`
 set pmmar=`ls $DIR/$RUNID$a.$Pmonth??003.nc | head -5`
 set pmapr=`ls $DIR/$RUNID$a.$Pmonth??004.nc | head -5`
 set pmmay=`ls $DIR/$RUNID$a.$Pmonth??005.nc | head -5`
 set pmjun=`ls $DIR/$RUNID$a.$Pmonth??006.nc | head -5`
 set pmjul=`ls $DIR/$RUNID$a.$Pmonth??007.nc | head -5`
 set pmaug=`ls $DIR/$RUNID$a.$Pmonth??008.nc | head -5`
 set pmsep=`ls $DIR/$RUNID$a.$Pmonth??009.nc | head -5`
 set pmoct=`ls $DIR/$RUNID$a.$Pmonth??010.nc | head -5`
 set pmnov=`ls $DIR/$RUNID$a.$Pmonth??011.nc | head -5`
 set pmdec=`ls $DIR/$RUNID$a.$Pmonth??012.nc | head -5`

if ($REINIT == 1) then
 set pels=`ls $pelist | head -60`
 set pbls=`ls $pblist | head -60`
else
 set pels=`ls $pelist | head -20`
 set pbls=`ls $pblist | head -20`
endif

 if ($count_all == $Numyr) then
  set cdo_blk=`ls $monlist | head -5`
 else
  set cdo_blk=`ls $DIRW/h[0123456789].nc $DIRW/i[0123456789].nc $DIRW/j[0123456789].nc $DIRW/k[0123456789].nc | head -5`
 endif

else if ($BLOCK == 2) then
 set pmjan=`ls $DIR/$RUNID$a.$Pmonth??001.nc | head -10 | tail -5` # tail +6 | head -5
 set pmfeb=`ls $DIR/$RUNID$a.$Pmonth??002.nc | head -10 | tail -5`
 set pmmar=`ls $DIR/$RUNID$a.$Pmonth??003.nc | head -10 | tail -5`
 set pmapr=`ls $DIR/$RUNID$a.$Pmonth??004.nc | head -10 | tail -5`
 set pmmay=`ls $DIR/$RUNID$a.$Pmonth??005.nc | head -10 | tail -5`
 set pmjun=`ls $DIR/$RUNID$a.$Pmonth??006.nc | head -10 | tail -5`
 set pmjul=`ls $DIR/$RUNID$a.$Pmonth??007.nc | head -10 | tail -5`
 set pmaug=`ls $DIR/$RUNID$a.$Pmonth??008.nc | head -10 | tail -5`
 set pmsep=`ls $DIR/$RUNID$a.$Pmonth??009.nc | head -10 | tail -5`
 set pmoct=`ls $DIR/$RUNID$a.$Pmonth??010.nc | head -10 | tail -5`
 set pmnov=`ls $DIR/$RUNID$a.$Pmonth??011.nc | head -10 | tail -5`
 set pmdec=`ls $DIR/$RUNID$a.$Pmonth??012.nc | head -10 | tail -5`

if ($REINIT == 1) then
 set pels=`ls $pelist | head -120 | tail -60` #tail +21 | head -20
 set pbls=`ls $pblist | head -120 | tail -60`
else
 set pels=`ls $pelist | head -40 | tail -20` #tail +21 | head -20
 set pbls=`ls $pblist | head -40 | tail -20`
endif

 if ($count_all == $Numyr) then
  set cdo_blk=`ls $monlist | head -10 | tail -5`
 else
  set cdo_blk=`ls $DIRW/h[0123456789].nc $DIRW/i[0123456789].nc $DIRW/j[0123456789].nc $DIRW/k[0123456789].nc | head -10 | tail -5`
 endif

else if ($BLOCK == 3) then
 set pmjan=`ls $DIR/$RUNID$a.$Pmonth??001.nc | head -15 | tail -5` # tail +11 | head -5
 set pmfeb=`ls $DIR/$RUNID$a.$Pmonth??002.nc | head -15 | tail -5`
 set pmmar=`ls $DIR/$RUNID$a.$Pmonth??003.nc | head -15 | tail -5`
 set pmapr=`ls $DIR/$RUNID$a.$Pmonth??004.nc | head -15 | tail -5`
 set pmmay=`ls $DIR/$RUNID$a.$Pmonth??005.nc | head -15 | tail -5`
 set pmjun=`ls $DIR/$RUNID$a.$Pmonth??006.nc | head -15 | tail -5`
 set pmjul=`ls $DIR/$RUNID$a.$Pmonth??007.nc | head -15 | tail -5`
 set pmaug=`ls $DIR/$RUNID$a.$Pmonth??008.nc | head -15 | tail -5`
 set pmsep=`ls $DIR/$RUNID$a.$Pmonth??009.nc | head -15 | tail -5`
 set pmoct=`ls $DIR/$RUNID$a.$Pmonth??010.nc | head -15 | tail -5`
 set pmnov=`ls $DIR/$RUNID$a.$Pmonth??011.nc | head -15 | tail -5`
 set pmdec=`ls $DIR/$RUNID$a.$Pmonth??012.nc | head -15 | tail -5`

if ($REINIT == 1) then
 set pels=`ls $pelist | head -180 | tail -60` # tail +41 | head -20
 set pbls=`ls $pblist | head -180 | tail -60`
else
 set pels=`ls $pelist | head -60 | tail -20` # tail +41 | head -20
 set pbls=`ls $pblist | head -60 | tail -20`
endif

 if ($count_all == $Numyr) then
  set cdo_blk=`ls $monlist | head -15 | tail -5`
 else
  set cdo_blk=`ls $DIRW/h[0123456789].nc $DIRW/i[0123456789].nc $DIRW/j[0123456789].nc $DIRW/k[0123456789].nc | head -15 | tail -5`
 endif

else if ($BLOCK == 4) then
 set pmjan=`ls $DIR/$RUNID$a.$Pmonth??001.nc | head -20 | tail -5`
 set pmfeb=`ls $DIR/$RUNID$a.$Pmonth??002.nc | head -20 | tail -5`
 set pmmar=`ls $DIR/$RUNID$a.$Pmonth??003.nc | head -20 | tail -5`
 set pmapr=`ls $DIR/$RUNID$a.$Pmonth??004.nc | head -20 | tail -5`
 set pmmay=`ls $DIR/$RUNID$a.$Pmonth??005.nc | head -20 | tail -5`
 set pmjun=`ls $DIR/$RUNID$a.$Pmonth??006.nc | head -20 | tail -5`
 set pmjul=`ls $DIR/$RUNID$a.$Pmonth??007.nc | head -20 | tail -5`
 set pmaug=`ls $DIR/$RUNID$a.$Pmonth??008.nc | head -20 | tail -5`
 set pmsep=`ls $DIR/$RUNID$a.$Pmonth??009.nc | head -20 | tail -5`
 set pmoct=`ls $DIR/$RUNID$a.$Pmonth??010.nc | head -20 | tail -5`
 set pmnov=`ls $DIR/$RUNID$a.$Pmonth??011.nc | head -20 | tail -5`
 set pmdec=`ls $DIR/$RUNID$a.$Pmonth??012.nc | head -20 | tail -5`

if ($REINIT == 1) then
 set pels=`ls $pelist | tail -60`
 set pbls=`ls $pblist | tail -60`
else
 set pels=`ls $pelist | tail -20`
 set pbls=`ls $pblist | tail -20`
endif

 if ($count_all == $Numyr) then
  set cdo_blk=`ls $monlist | tail -5`
 else
  set cdo_blk=`ls $DIRW/h[0123456789].nc $DIRW/i[0123456789].nc $DIRW/j[0123456789].nc $DIRW/k[0123456789].nc | head -20 | tail -5`
 endif

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

cdo mergetime $pels Timeseries_5yrs.nc
cdo mergetime $pbls Tempseries_5yrs.nc
#cdo copy $pels Timeseries_5yrs.nc
#cdo copy $pbls Tempseries_5yrs.nc

cdo copy $cdo_blk Mmonthly_means_5yrs.nc
cdo yearavg Mmonthly_means_5yrs.nc yearly_means_5yrs.nc

 cdo yseasavg Tempseries_5yrs.nc Tseasonal_means_5yrs.nc
 cdo ymonavg Tempseries_5yrs.nc Tmonthly_means_5yrs.nc

 mv seasonal_5yrs.nc seasonal_means_5yrs.nc
 mv monthly_5yrs.nc monthly_means_5yrs.nc

# cd $DIRW

else
 echo "(1)     $date - Block Directory Already Exists"
 echo "(1)     Please Delete Directory and Run Again"
 exit (1)
endif

rm ave*.nc

#end

exit

