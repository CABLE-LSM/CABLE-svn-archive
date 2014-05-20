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

#set count_all=`ls $DIRW/h[0123456789].nc $DIRW/i[0123456789].nc $DIRW/j[0123456789].nc $DIRW/k[0123456789].nc | wc -l`
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
 set pelist=`ls $DIR/$RUNID.$Ptimes-??????????.nc | head -${nom}`
 #set peswlw=`ls $DIR/$RUNID.$Ptimes-??????????_swlw.nc | head -${nom}`
 #set penosw=`ls $DIR/$RUNID.$Ptimes-??????????_noswlw.nc | head -${nom}`
 set pblist=`ls $DIR/$RUNID.$Ptemps-??????????.nc | head -${nom}`
 set pmlist=`ls $DIR/$RUNID.$Pmonth-??????????.nc | head -${nom}`
else
 set pelist=`ls $DIR/$RUNID.$Ptimes-??????????.nc | head -${noy}`
 #set peswlw=`ls $DIR/$RUNID.$Ptimes-??????????_swlw.nc | head -${noy}`
 #set penosw=`ls $DIR/$RUNID.$Ptimes-??????????_noswlw.nc | head -${noy}`
 set pblist=`ls $DIR/$RUNID.$Ptemps-??????????.nc | head -${noy}`
 set pmlist=`ls $DIR/$RUNID.$Pmonth-??????????.nc | head -${noy}`
endif

 echo $bblk
 set pmjan=`ls $DIR/$RUNID.$Pmonth-????001???.nc | head -${bblk} | tail -5`
 set pmfeb=`ls $DIR/$RUNID.$Pmonth-????002???.nc | head -${bblk} | tail -5`
 set pmmar=`ls $DIR/$RUNID.$Pmonth-????003???.nc | head -${bblk} | tail -5`
 set pmapr=`ls $DIR/$RUNID.$Pmonth-????004???.nc | head -${bblk} | tail -5`
 set pmmay=`ls $DIR/$RUNID.$Pmonth-????005???.nc | head -${bblk} | tail -5`
 set pmjun=`ls $DIR/$RUNID.$Pmonth-????006???.nc | head -${bblk} | tail -5`
 set pmjul=`ls $DIR/$RUNID.$Pmonth-????007???.nc | head -${bblk} | tail -5`
 set pmaug=`ls $DIR/$RUNID.$Pmonth-????008???.nc | head -${bblk} | tail -5`
 set pmsep=`ls $DIR/$RUNID.$Pmonth-????009???.nc | head -${bblk} | tail -5`
 set pmoct=`ls $DIR/$RUNID.$Pmonth-????010???.nc | head -${bblk} | tail -5`
 set pmnov=`ls $DIR/$RUNID.$Pmonth-????011???.nc | head -${bblk} | tail -5`
 set pmdec=`ls $DIR/$RUNID.$Pmonth-????012???.nc | head -${bblk} | tail -5`

if ($REINIT == 1) then
 set pels=`ls $pelist | head -${mblk} | tail -60`
 #set pesw=`ls $peswlw | head -${mblk} | tail -60`
 #set peno=`ls $penosw | head -${mblk} | tail -60`
 set pbls=`ls $pblist | head -${mblk} | tail -60`
else
 set pels=`ls $pelist | head -${yblk} | tail -20`
 #set pesw=`ls $peswlw | head -${yblk} | tail -20`
 #set peno=`ls $penosw | head -${yblk} | tail -20`
 set pbls=`ls $pblist | head -${yblk} | tail -20`
endif
 set pmls=`ls $pmlist | head -${mblk} | tail -60`

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
 #cdo mergetime $pesw Timeseries_5yrs_swlw.nc
 #cdo mergetime $peno Timeseries_5yrs_noswlw.nc
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

