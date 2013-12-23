#!/bin/csh 

unalias ls

#----------------------------
set date=`date`
set DIRW=$PWD
#----------------------------

echo " "
dmget $DIR/$RUNID.$Pdaily-??????????.nc
dmget $DIR/$RUNID.$Ptemps-??????????.nc 
dmget $DIR/$RUNID.$Ptimes-??????????.nc 
dmget $DIR/$RUNID.$Pmonth-??????????.nc

set bblk=  5*${BLOCK}
set mblk= 60*${BLOCK}
set yblk= 20*${BLOCK}
@ noy  = ($FullYrs * 4)
@ nom  = ($FullYrs * 12)
@ Numyr= ($FullYrs + 1)

#set count_all=`/bin/ls $DIRW/h[0123456789].nc $DIRW/i[0123456789].nc $DIRW/j[0123456789].nc $DIRW/k[0123456789].nc | wc -l`
#@ countm1 = $count_all - 1
#@ countm2 = $count_all - 2
#if ( $count_all > $FullYrs ) then # need to  limit cdo_merge
#
#  echo "(1)     Warning: Number of Years = YR+1. Please Check."
#  echo "(1)     Warning: Or Year Doesn't Start in January. Please Check."
#if ($Numyr > $FullYrs) then
#  echo "(1)     Warning: Number of Years Greater Than Full Number of Years"
#  echo "(1)     Please Check"
#  echo "(1)     Hint: Move Undesired Files to a Temporary Folder"
#endif

if (! -d $DIRW/block${BLOCK}_5yrs) then
 mkdir $DIRW/block${BLOCK}_5yrs
 cd $DIRW/block${BLOCK}_5yrs

if ($REINIT == 1) then
 set pelist=`/bin/ls $DIR/$RUNID.$Ptimes-??????????.nc | head -${nom}`
 set pblist=`/bin/ls $DIR/$RUNID.$Ptemps-??????????.nc | head -${nom}`
 set pmlist=`/bin/ls $DIR/$RUNID.$Pmonth-??????????.nc | head -${nom}`
else
 set pelist=`/bin/ls $DIR/$RUNID.$Ptimes-??????????.nc | head -${noy}`
 set pblist=`/bin/ls $DIR/$RUNID.$Ptemps-??????????.nc | head -${noy}`
 set pmlist=`/bin/ls $DIR/$RUNID.$Pmonth-??????????.nc | head -${noy}`
endif

 echo $bblk
 set pmjan=`/bin/ls $DIR/$RUNID.$Pmonth-????001???.nc | head -${bblk} | tail -5`
 set pmfeb=`/bin/ls $DIR/$RUNID.$Pmonth-????002???.nc | head -${bblk} | tail -5`
 set pmmar=`/bin/ls $DIR/$RUNID.$Pmonth-????003???.nc | head -${bblk} | tail -5`
 set pmapr=`/bin/ls $DIR/$RUNID.$Pmonth-????004???.nc | head -${bblk} | tail -5`
 set pmmay=`/bin/ls $DIR/$RUNID.$Pmonth-????005???.nc | head -${bblk} | tail -5`
 set pmjun=`/bin/ls $DIR/$RUNID.$Pmonth-????006???.nc | head -${bblk} | tail -5`
 set pmjul=`/bin/ls $DIR/$RUNID.$Pmonth-????007???.nc | head -${bblk} | tail -5`
 set pmaug=`/bin/ls $DIR/$RUNID.$Pmonth-????008???.nc | head -${bblk} | tail -5`
 set pmsep=`/bin/ls $DIR/$RUNID.$Pmonth-????009???.nc | head -${bblk} | tail -5`
 set pmoct=`/bin/ls $DIR/$RUNID.$Pmonth-????010???.nc | head -${bblk} | tail -5`
 set pmnov=`/bin/ls $DIR/$RUNID.$Pmonth-????011???.nc | head -${bblk} | tail -5`
 set pmdec=`/bin/ls $DIR/$RUNID.$Pmonth-????012???.nc | head -${bblk} | tail -5`

if ($REINIT == 1) then
 set pels=`/bin/ls $pelist | head -${mblk} | tail -60`
 set pbls=`/bin/ls $pblist | head -${mblk} | tail -60`
else
 set pels=`/bin/ls $pelist | head -${yblk} | tail -20`
 set pbls=`/bin/ls $pblist | head -${yblk} | tail -20`
endif
 set pmls=`/bin/ls $pmlist | head -${mblk} | tail -60`

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
 ncrcat -O avejan.nc avefeb.nc avemar.nc aveapr.nc avemay.nc avejun.nc avejul.nc aveaug.nc avesep.nc aveoct.nc avenov.nc avedec.nc monthly_means_5yrs.nc
 ncrcat -O avedjf.nc avemam.nc avejja.nc aveson.nc seasonal_means_5yrs.nc

 cdo mergetime $pels Timeseries_5yrs.nc
 cdo mergetime $pbls Tempseries_5yrs.nc
 cdo mergetime $pmls Mmonthly_means_5yrs.nc
 cdo yearmean Mmonthly_means_5yrs.nc yearly_means_5yrs.nc
 cdo yseasmean Tempseries_5yrs.nc Tseasonal_means_5yrs.nc
 cdo ymonmean Tempseries_5yrs.nc Tmonthly_means_5yrs.nc

else
 echo "(1)     $date - Block Directory Already Exists"
 echo "(1)     Please Delete Directory and Run Again"
 exit (1)
endif

rm ave*.nc

exit

