#!/bin/csh
#---------------------------------------------------------------------------
# UMPLOT v2.1.2
# Created by Lauren Stevens, email: Lauren.Stevens@csiro.au
# 2008-2014
#---------------------------------------------------------------------------

#*******Edit********
if ($HOSTNAME == ruby) then
 setenv USERID ste69f  # ruby/Burnet/Shine/Pearcey
else
 setenv USERID lxs599  # NCI- Raijin/Accessdev
endif
#*******************

echo ""
echo "=================================================================="
echo "                          UMPLOT v2.1.2                           "
echo "                   Developed by Lauren Stevens                    "
echo "                 Email: Lauren.Stevens@csiro.au                   "
echo "=================================================================="
echo ""

# SET UP - see Chapter 2 in Documentation.
# You must have certain modules loaded before running UMPLOT. E.g.
module load ncl/6.1.2 #/6.1.0-beta or later
#module load nco
module load cdo/1.6.1
#module load python
#module load cdat

if ( ! -e ~/.hluresfile ) then
 echo " "
 echo "(1)     You Do Not Have a ~/.hluresfile, Please Copy One to $HOME"
 echo "(1)     It Suggests That You May Be Unfamiliar with NCL"
 echo "(1)     Please Make Sure NCL is Loaded Properly, For Further Details"
 echo "(1)     Check the Website http://www.ncl.ucar.edu/"
 echo " "
 exit (1)
endif

# Namelist replaces Questions - see Chapter 3 in Documentation.
# Output Files Generally Located in $HOME/access.

if (! -d $HOME/umplot) then
 mkdir $HOME/umplot
 echo " "
 echo "(0)     Making Directory $HOME/umplot"
endif
if (! -d $HOME/umplot/nml) then
 mkdir $HOME/umplot/nml
 if ($HOSTNAME == ruby) then
  cp ~ste69f/umplot/nml/example_*.nml $HOME/umplot/nml
 else
  cp ~lxs599/umplot/nml/example_*.nml $HOME/umplot/nml
 endif
 echo " "
 echo "(0)     Making Directory $HOME/umplot/nml"
 echo "(0)     Copying Example NML to $HOME/umplot/nml"
 echo "(0)     Please Copy an Example NML to $HOME/umplot/nml/umplot.nml and Modify"
endif

#=================
setenv DIRW $PWD
#=================

if ( -e $HOME/umplot/nml/umplot.nml ) then
  source $HOME/umplot/nml/umplot.nml
else
 echo " "
 echo "(1)     Namelist File Does Not Exist"
 echo "(1)     Please Create Namelist File in $HOME/umplot/nml"
 echo "(1)     Example Namelist File: ruby:~ste69f/umplot/nml/umplot.nml"
 echo "(1)     and/or Email Lauren.Stevens@csiro.au"
 echo " "
 exit (1)
endif

set rfile=`ls $DIR/ave*a.p* | wc -l`
if ($rfile > 0) then
 echo " "
 echo "(1)     RUNID is similar to processed files ave*.nc"
 echo "(1)     WARNING: Therefore necessary files may be deleted"
 echo "(1)     Please Email Lauren.Stevens@csiro.au For Help"
 echo " "
 exit (1)
endif
echo " "

#if ($SELYR == y) then
# @ YR = 1 + $EYR - $SYR
#endif
#=================
setenv FullYrs $YR
#=================

ncl ~$USERID/umplot/mod_func.ncl
# WARNING: Can't put anything between this line and status
if ( $status == 38 ) then
 @ noblks = ${FullYrs} / 5
else
 @ noblks = 1
endif
echo ""
if ( $noblks > 1 ) then
 if ($SPLIT == y) then
  setenv YR 5
 endif
else
 setenv SPLIT n
endif

if ($MODEL == c && $RES >= 96) then
 setenv TILE 17
 setenv SOIL 6
else
 setenv TILE 9
 setenv SOIL 4
endif

if ($RES == 320) then
 setenv TSTEP 288
else
 setenv TSTEP 48
endif

# NCI Transfer =============================================
#set nsdate=`date`
#if (${?NCI_CP} == 1) then
if ($NCI_CP == y) then
 echo "(0)     Transferring Fields Files to ruby"
 echo ""
 ~$USERID/umplot/tools/rsync_raijin.sh
 echo ""
else
 echo "(0)     No Fields Files Transferred from Raijin"
 echo ""
endif
#else
# echo "(1)     NCI_CP not set. Please update umplot.nml"
# echo "(1)     Contact Lauren.Stevens@csiro.au"
#endif
#set nedate=`date`

# if ~/access doesn't exist ?
if ($HOSTNAME == ruby) then
 if ($DIRW == ~/access) then
  setenv DIRW ~/access/$RUNID
  cd $DIRW
 endif
#else
#/short/p66/$NCI_ID/
endif
# NCI Transfer =============================================

# Nc Convert And Check =====================================
set a=a
if (${?CPL} == 0) then
 setenv CPL n
else
 echo "(1)     NOTE: CPL is Set and Set to" $CPL
 echo ""
endif
if ($CPL == n) then
 set pnc=`ls $DIR/$RUNID$a.p??????.nc | wc -l`
else
 set pnc=`ls $DIR/$RUNID.p?-??????????.nc | wc -l`
endif
if (${?CNM} == 0) then
 setenv CNM A
endif

#Check - so you don't edit someone else's $DIR
if ( $CNV2NC == y ) then
 if ($HOSTNAME == ruby) then
  set nchar=17
 else
  set nchar=16 # for /home
  #set ncharS=17 # for /short
  #set usid=`echo "$PWD" | head -c3` #/sh /ho
 endif
 set udir1=`echo $PWD | head -c$nchar | tail -c6`
 cd $DIR
 set udir2=`echo $PWD | head -c$nchar | tail -c6`
 cd $DIRW
 if ($udir1 != $udir2) then
  echo "Directories:" $udir1 $udir2
  setenv CNV2NC n
  echo "(1)     You Were About to Edit Someone Else's Directory"
  echo "(1)     Check How DIR is Set in umplot.nml"
  #exit(1)
 endif
endif

if ($CNV2NC == y) then

  #set csdate=`date`
   #echo ""
   echo "(0)     Converting Fields Files to NetCDF"
   echo ""
   cd $DIR
   ~$USERID/umplot/conv_um2nc.sh
  #set cedate=`date`
 if ($CPL == n) then
   echo " "
   echo "(0)     Extracting Tmax and Tmin to" $Ptemps "File"
   ~$USERID/umplot/tools/extract_pb.sh
 endif
   cd $DIRW

  #echo " "
  echo "(0)     Converted Fields Files to NetCDF before Running UMPLOT"

else #CNV2NC == n

 if ( $pnc == 0 && ! -e $DIRW/seasonal_means_${YR}yrs.nc ) then
 #if ( || $SPLIT == y && ! -e $DIRW/block1_5yrs/seasonal_means_${YR}yrs.nc ) then
  echo " "
  echo "(1)     Fields Files Not Converted to NetCDF"
  echo "(1)     Cannot Run UMPLOT until NetCDFs are Present"
  echo "(1)     Either Convert Offline or Switch CNV2NC to y in Namelist"
  echo " "
  exit (1)
 endif #pnc

 echo " "
 echo "(0)     Converted Fields Files to NetCDF before Running UMPLOT"

endif #CNV2NC

#cd $DIR
#if ( $SELYR == y ) then
# need to do for both SELYR cases?
# ~$USERID/umplot/{workdir/selyr/}umformat2years.sh
# setenv CPL y
#endif
#cd $DIRW

if ($CPL == n) then
 set pnc2=`ls $DIR/$RUNID$a.p??????.nc | wc -l`
 set chk=`ls -s $DIR/$RUNID$a.p??????.nc`
else
 set pnc2=`ls $DIR/$RUNID.p?-??????????.nc | wc -l`
 set chk=`ls -s $DIR/$RUNID.p?-??????????.nc`
endif
#set cnt=`find -size 4 | wc -l`
set i = 0
 while ($i < $pnc2)
  @ k = 2 * $i + 1  # 2*$i-1 for i=1
  @ m = 2 * $i + 2  # 2*$i   for i=1
   if ( $chk[$k] == "4") then
    echo " "
    echo "(1)     Not all NetCDF Have Been Converted Properly"
    echo " "
    echo "(1)     E.g. See File" $chk[$m]
    echo " "
    echo "(1)     If You See Netcdf (i.e. $RUNID$a.p*.nc) with SIZE 32"
    echo "(1)     You'll Need to Delete the UM File(s) and the Netcdf(s)"
    echo "(1)     Or Set CNV2NC = n Once the Netcdf(s) Have been Removed"
    echo "(1)     WARNING: There Could be More than one File with Size 32"
    exit (1)
   endif
  @ i ++
 end
# Nc Convert And Check =====================================

# Lest 12.11.13 will move to umplot.nml
if (${?UMPLT} == 0) then
 setenv UMPLT y
else
 echo " "
 echo "(1)     NOTE: UMPLT is Set and Set to" $UMPLT
endif

if (${?REGION} == 0) then
 setenv REGION 2
else
 echo " "
 echo "(1)     NOTE: REGION is Set and Set to" $REGION
endif

if ($UMPLT == y) then

echo " "
echo "========================= Running UMPLOT ========================="
echo " "
echo "You have Set the Parameters of this UMPLOT Run as:"
echo "Jobid:" $RUNID "with files in Directory:" $DIR
echo "Running for" $FullYrs "years with" $SPLIT "5yr Splitting," $BMRK "Extra Plots, and" $CHFMT $FMT "files"
if ($BMRK == y) then
echo "You are Comparing results with" $BDIR
endif

echo " "
#set a=a
set sdate=`date`
echo "Start Time:" $sdate
echo " "
echo "=================================================================="
echo " "

if ($DIRW == $DIR) then
 echo "(0)     Working directory is the same as Data directory"
else
 echo "(1)     Working directory is NOT the same as Data directory"

 # in $DIRW
 if ($CPL == n) then
  set rid=`ls ?????a.p?????? | head -1 | head -c5`
 else
  set rid=`ls *.p?-?????????? | head -1 | head -c5`
 endif
 #set rid2=`ls $DIR/?????a.p?????? | head -1 | head -c5`
 setenv RID $rid
 #setenv RID2 $rid2

 if ( ${#rid} > 0 ) then
  #if ( $RID != $RUNID || $RID2 != $RUNID ) then
  if ( $RID != $RUNID ) then
   echo " "
   echo "(1)     WARNING: Different RUNIDs Present"
   echo "(1)     WARNING: You May Be Running A Different Job"
   echo "(1)     WARNING: Check that you have the Correct NML or DIR"
   exit(1)
  endif
 endif

endif

# Nc Check =================================================
if ($CPL == n) then
 set pmnc=`ls $DIR/$RUNID$a.$Pmonth?????.nc | wc -l`
else
 set pmnc=`ls $DIR/$RUNID.$Pmonth-??????????.nc | wc -l`
endif
 if ( $pmnc == 0 && ! -e $DIRW/seasonal_means_${YR}yrs.nc ) then
    echo " "
    echo "(1)     If you get a message like: 'please enter files interactively'"
    echo "(1)     Make Sure Monthly Netcdf Files are in Directory"
    echo "(1)     i.e. $RUNID$a.$Pmonth?????.nc"
 endif
# Nc Check =================================================

echo " "

# ==================================================================================

if ($SPLIT == y) then

if ( -d $DIRW/block1_5yrs ) then
 echo "(1)     $sdate - Block Directory Already Exists"
 echo "(1)     Please Move or Delete Directory and Run Again"
 exit (1)
endif

# Processing ===============================================
if ($CPL == n) then
~$USERID/umplot/cdo_merge.sh
else
~$USERID/umplot/cdo_merge_cpl.sh
#~$USERID/umplot/{workdir/selyr/}process_years.sh
endif

#foreach blk ( $block )
#setenv BLOCK $blk

@ j = 1
while ( $j <= $noblks )
setenv BLOCK ${j}

source ~$USERID/umplot/nml/envars.sh
echo " "
echo "(0)     Post-Processing NetCDF Files"

if ($CPL == n) then
 ~$USERID/umplot/split_seasonal.sh
else
 ~$USERID/umplot/split_seas_cpl.sh
endif

cd $DIRW/block${BLOCK}_5yrs

set flen=`cdo ntime $DIRW/block${BLOCK}_5yrs/Mmonthly_means_${YR}yrs.nc`
@ nom = ($YR * 12)
@ nof = ($FullYrs * 12)
if ($flen != $nom || $pmnc < $nof) then
echo " "
echo "(1)     Number of Months is Incorrect"
echo "(1)     Please Check Length of Mmonthly_means_${YR}yrs.nc"
echo "(1)     It May Be That Not All .p Files were converted to NetCDF"
echo "(1)     Or not all Fields Files were Transferred"
#echo "(1)     "
#echo "(1)     If This is the (Only) Issue, 'setenv RUNID' in $DIR"
#echo "(1)     and"
#echo "(1)     If REINIT=1, Run ~/umplot/umtonc_ss_reinit.sh on Command Line"
#echo "(1)     or"
#echo "(1)     If REINIT=3, Run ~/umplot/runallv.sh on Command Line"
#echo "(1)     Then Run umplot Again"
exit (1)
endif

ncrename -v $tname,tscrn seasonal_means_${YR}yrs.nc
ncrename -v $tname,tscrn monthly_means_${YR}yrs.nc
#ncrename -v $tname,tscrn Mmonthly_means_${YR}yrs.nc
ncrename -v $tname,tscrn yearly_means_${YR}yrs.nc

~$USERID/umplot/process_carbon.sh
# Processing ===============================================

echo " "
ncl ~$USERID/umplot/run.ncl
ncl ~$USERID/umplot/global_means_wbal.ncl
ncl ~$USERID/umplot/global_means_tables.ncl
if ($BMRK == y) then
 if ( -e ${CABLE}/Timeseries_${YR}yrs.nc && -e ${MOSES}/Timeseries_${YR}yrs.nc ) then
  ncl ~$USERID/umplot/global_means_ctables.ncl
 endif
endif
#if ( -e mm.*.nc ) then
 ncl ~$USERID/umplot/carbon_zonal_means.ncl
#endif
ncl ~$USERID/umplot/global_means_carbon.ncl
#ncl ~$USERID/umplot/daily_mean_timeseries.ncl

# Timeseries ===============================================
if( -e $DIRW/block${BLOCK}_5yrs/Timeseries_${YR}yrs.nc ) then
 if ($CAL == 360 ) then
  python ~$USERID/umplot/MMDC.py $YR
 else
  set years=`cdo showyear $DIRW/block${BLOCK}_5yrs/Timeseries_${YR}yrs.nc`
  echo ""
  python ~$USERID/umplot/MMDC_365cal.py $YR $years[1] $years[$YR] #[$#years]
 endif
 # non-Jan start: python $PLOT/MMDC_roll.py $YR
 mv MMDC_${YR}yrs.nc MeanMnthDailyCycles_${YR}yrs.nc
 python ~$USERID/umplot/mmdc_roll.py $YR

 if ($CAL == 360) then
  python ~$USERID/umplot/AnnualCycle.py $YR
 else
  set years=`cdo showyear $DIRW/block${BLOCK}_5yrs/Timeseries_${YR}yrs.nc`
  echo ""
  python ~$USERID/umplot/AnnualCycle_365cal.py $YR $years[1] $years[$YR] #[$#years]
 #cdo ymonmean Timeseries_${YR}yrs.nc AnnCycle_${YR}yrs.nc
 endif
 # non-Jan start: python $PLOT/AnnCycle_roll.py $YR
 mv AnnCycle_${YR}yrs.nc AnnualCycle_${YR}yrs.nc
 
 # New - getting Tested
 #ncl process_timeseries_1.ncl
 #ncl process_timeseries_1.ncl
 #python process_timeseries_2.py $YR $years[1] $years[$YR] False
 #python process_timeseries_2.py $YR $years[1] $years[$YR] True

 # Jaeger Figure6 ----------
 ncl ~$USERID/umplot/mmdc.ncl
 ncl ~$USERID/umplot/mmdc_jan.ncl
 ncl ~$USERID/umplot/mmdc_jul.ncl
 ncl ~$USERID/umplot/annualcycles.ncl
endif # if Timeseries
# Timeseries ===============================================

# Benchmarking =============================================
if ($BMRK == y) then
 if ( -e ${CABLE}/seasonal_means_${YR}yrs.nc && -e ${MOSES}/seasonal_means_${YR}yrs.nc) then
  ncl ~$USERID/umplot/Taylor_diagram.ncl
  ncl ~$USERID/umplot/Mon_IntAnnVar.ncl
  ncl ~$USERID/umplot/Mon_IntAnnVar_Aust.ncl
  if ($YR > 1) then
   ncl ~$USERID/umplot/Yr_IntAnnVar.ncl
   ncl ~$USERID/umplot/Yr_IntAnnVar_Aust.ncl
  endif
  ncl ~$USERID/umplot/bias_3panel.ncl
  ncl ~$USERID/umplot/zonal_pr.ncl
  ncl ~$USERID/umplot/zonal_clt.ncl
  ncl ~$USERID/umplot/plot_amoj6p_precip.ncl
  ncl ~$USERID/umplot/plot_amoj6p_tscrn.ncl
  if ( -e ${CABLE}/Tseasonal_means_${YR}yrs.nc && -e ${MOSES}/Tseasonal_means_${YR}yrs.nc) then
   ncl ~$USERID/umplot/plot_amoj6p_tmin.ncl
   ncl ~$USERID/umplot/plot_amoj6p_tmax.ncl
  endif
  ncl ~$USERID/umplot/plot_amoj6p_clouds.ncl
  ncl ~$USERID/umplot/plot_diff6p_surfalb.ncl
  ncl ~$USERID/umplot/plot_amoj6p_trunoff.ncl
  ncl ~$USERID/umplot/barchart.ncl
  #ncl chv=0 ~$USERID/umplot/barchart_amp.ncl
  #ncl chv=1 ~$USERID/umplot/barchart_amp.ncl
  #ncl ~$USERID/umplot/barchart_clt.ncl
 else
  echo "(1)     Cannot plot Taylor Diagram and Panel Plots - Necessary Files Don't Exist"
  echo "(1)     e.g. seasonal_means_${YR}yrs.nc"
  echo ""
 endif
 if ( -e ${CABLE}/Timeseries_${YR}yrs.nc && -e ${MOSES}/Timeseries_${YR}yrs.nc) then
  ncl ~$USERID/umplot/mmdc_hyytiala.ncl
  ncl ~$USERID/umplot/mmdc_bondville.ncl
  ncl ~$USERID/umplot/mmdc_hay.ncl
  ncl ~$USERID/umplot/mmdc_walker_branch.ncl
  ncl ~$USERID/umplot/mmdc_tharandt.ncl
  ncl ~$USERID/umplot/mmdc_tumbarumba.ncl
  ncl ~$USERID/umplot/mmdc_little_washita.ncl
  ncl ~$USERID/umplot/mmdc_vielsalm.ncl
  ncl ~$USERID/umplot/mmdc_nsaboreas.ncl
  ncl ~$USERID/umplot/mmdc_harvard.ncl
  ncl ~$USERID/umplot/mmdc_loobos.ncl
  ncl ~$USERID/umplot/mmdc_manaus.ncl
 else
  echo "(1)     Cannot plot Fluxnet with Timeseries - Necessary Files Don't Exist"
  echo "(1)     e.g. Timeseries_${YR}yrs.nc"
  echo ""
 endif
endif
# Benchmarking =============================================

#rm Timeseries_$YR\yrs.nc Tempseries_$YR\yrs.nc
cd $DIRW 

 @ j++
end

if ($CPL == n) then
rm h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc
else
rm y????.nc
endif

# ==================================================================================

else if ($SPLIT == n) then

#mkdir netcdf
#mkdir dumps

source ~$USERID/umplot/nml/envars.sh

# Processing ===============================================
if (! -e $DIRW/seasonal_means_${YR}yrs.nc ) then
echo "(0)     Post-Processing NetCDF Files"
echo " "

if ($CPL == n) then
 ~$USERID/umplot/seasonal.sh
 ~$USERID/umplot/cdo_merge.sh
else
 ~$USERID/umplot/seas_cpl.sh
 ~$USERID/umplot/cdo_merge_cpl.sh
endif

set flen=`cdo ntime $DIRW/Mmonthly_means_${YR}yrs.nc`
@ nom = ($YR * 12)
if ($flen != $nom) then
echo " "
echo "(1)     Number of Months is Incorrect"
echo "(1)     Please Check Length of Mmonthly_means_${YR}yrs.nc"
echo "(1)     It May Be That Not All .p Files were converted to NetCDF"
echo "(1)     Or not all Fields Files were Transferred"
#echo "(1)     "
#echo "(1)     If This is the (Only) Issue, 'setenv RUNID' in $DIR"
#echo "(1)     and"
#echo "(1)     If REINIT=1, Run ~/umplot/umtonc_ss_reinit.sh on Command Line"
#echo "(1)     or"
#echo "(1)     If REINIT=3, Run ~/umplot/runallv.sh on Command Line"
#echo "(1)     Then Run umplot Again"
exit (1)
endif

ncrename -v $tname,tscrn seasonal_means_${YR}yrs.nc
ncrename -v $tname,tscrn monthly_means_${YR}yrs.nc
#ncrename -v $tname,tscrn Mmonthly_means_${YR}yrs.nc
ncrename -v $tname,tscrn yearly_means_${YR}yrs.nc

~$USERID/umplot/process_carbon.sh

else
echo " "
echo "(0)     NetCDF Files Already Exist. Will Use These Files to Plot."
echo "(1)     Suggestion: Remove Created _${YR}yrs.nc Files if You Want to Start Again."

endif # if ( ! -e $DIRW/season_${YR}yrs.nc )
# Processing ===============================================

echo " "
# Interface bw Shell and NCL - runs major plotting routines
ncl ~$USERID/umplot/run.ncl
ncl ~$USERID/umplot/global_means_wbal.ncl
ncl ~$USERID/umplot/global_means_tables.ncl
if ($BMRK == y) then
 if ( -e ${CABLE}/Timeseries_${YR}yrs.nc && -e ${MOSES}/Timeseries_${YR}yrs.nc ) then
  ncl ~$USERID/umplot/global_means_ctables.ncl
 endif
endif
#if ( -e mm.*.nc ) then
 ncl ~$USERID/umplot/carbon_zonal_means.ncl
#endif
ncl ~$USERID/umplot/global_means_carbon.ncl
#ncl ~$USERID/umplot/daily_mean_timeseries.ncl

# Timeseries ===============================================
if ( -e $DIRW/Timeseries_${YR}yrs.nc ) then
 if ($CAL == 360) then
  python ~$USERID/umplot/MMDC.py $YR
 else
  set years=`cdo showyear $DIRW/Timeseries_${YR}yrs.nc`
  echo ""
  python ~$USERID/umplot/MMDC_365cal.py $YR $years[1] $years[$YR] #[$#years]
 endif
 # non-Jan start: python $PLOT/MMDC_roll.py $YR
 mv MMDC_${YR}yrs.nc MeanMnthDailyCycles_${YR}yrs.nc
 python ~$USERID/umplot/mmdc_roll.py $YR

 if ($CAL == 360) then
  python ~$USERID/umplot/AnnualCycle.py $YR
 else
  set years=`cdo showyear $DIRW/Timeseries_${YR}yrs.nc`
  echo ""
  python ~$USERID/umplot/AnnualCycle_365cal.py $YR $years[1] $years[$YR] #[$#years]
 #cdo ymonmean Timeseries_${YR}yrs.nc AnnCycle_${YR}yrs.nc
 endif
 #non-Jan start: python $PLOT/AnnCycle_roll.py $YR
 mv AnnCycle_${YR}yrs.nc AnnualCycle_${YR}yrs.nc

 # New - getting Tested
 #ncl process_timeseries_1.ncl
 #ncl process_timeseries_1.ncl
 #python process_timeseries_2.py $YR $years[1] $years[$YR] False
 #python process_timeseries_2.py $YR $years[1] $years[$YR] True

 # Jaeger Figure6 ----------
 ncl ~$USERID/umplot/mmdc.ncl
 ncl ~$USERID/umplot/mmdc_jan.ncl
 ncl ~$USERID/umplot/mmdc_jul.ncl
 ncl ~$USERID/umplot/annualcycles.ncl
endif # if Timeseries
# Timeseries ===============================================

# Benchmarking =============================================
if ($BMRK == y) then
 if ( -e ${CABLE}/seasonal_means_${YR}yrs.nc && -e ${MOSES}/seasonal_means_${YR}yrs.nc) then
  ncl ~$USERID/umplot/Taylor_diagram.ncl
  ncl ~$USERID/umplot/Mon_IntAnnVar.ncl
  ncl ~$USERID/umplot/Mon_IntAnnVar_Aust.ncl
  if ($YR > 1) then
   ncl ~$USERID/umplot/Yr_IntAnnVar.ncl
   ncl ~$USERID/umplot/Yr_IntAnnVar_Aust.ncl
  endif
  ncl ~$USERID/umplot/bias_3panel.ncl
  ncl ~$USERID/umplot/zonal_pr.ncl
  ncl ~$USERID/umplot/zonal_clt.ncl
  ncl ~$USERID/umplot/plot_amoj6p_precip.ncl
  ncl ~$USERID/umplot/plot_amoj6p_tscrn.ncl
  if ( -e ${CABLE}/Tseasonal_means_${YR}yrs.nc && -e ${MOSES}/Tseasonal_means_${YR}yrs.nc) then
   ncl ~$USERID/umplot/plot_amoj6p_tmin.ncl
   ncl ~$USERID/umplot/plot_amoj6p_tmax.ncl
  endif
  ncl ~$USERID/umplot/plot_amoj6p_clouds.ncl
  ncl ~$USERID/umplot/plot_diff6p_surfalb.ncl
  ncl ~$USERID/umplot/plot_amoj6p_trunoff.ncl
  ncl ~$USERID/umplot/barchart.ncl
  #ncl chv=0 ~$USERID/umplot/barchart_amp.ncl
  #ncl chv=1 ~$USERID/umplot/barchart_amp.ncl
  #ncl ~$USERID/umplot/barchart_clt.ncl
 else
  echo "(1)     Cannot plot Taylor Diagram and Panel Plots - Necessary Files Don't Exist"
  echo "(1)     e.g. seasonal_means_${YR}yrs.nc"
  echo ""
 endif
 if ( -e ${CABLE}/Timeseries_${YR}yrs.nc && -e ${MOSES}/Timeseries_${YR}yrs.nc) then
  ncl ~$USERID/umplot/mmdc_hyytiala.ncl
  ncl ~$USERID/umplot/mmdc_bondville.ncl
  ncl ~$USERID/umplot/mmdc_hay.ncl
  ncl ~$USERID/umplot/mmdc_walker_branch.ncl
  ncl ~$USERID/umplot/mmdc_tharandt.ncl
  ncl ~$USERID/umplot/mmdc_tumbarumba.ncl
  ncl ~$USERID/umplot/mmdc_little_washita.ncl
  ncl ~$USERID/umplot/mmdc_vielsalm.ncl
  ncl ~$USERID/umplot/mmdc_nsaboreas.ncl
  ncl ~$USERID/umplot/mmdc_harvard.ncl
  ncl ~$USERID/umplot/mmdc_loobos.ncl
  ncl ~$USERID/umplot/mmdc_manaus.ncl
 else
  echo "(1)     Cannot plot Fluxnet with Timeseries - Necessary Files Don't Exist"
  echo "(1)     e.g. Timeseries_${YR}yrs.nc"
  echo ""
 endif
endif
# Benchmarking =============================================

#rm Timeseries_$YR\yrs.nc Tempseries_$YR\yrs.nc

mv *.ps plots/
mv *.eps plots/
#mv *.ps plots_${YR}yrs/
#mv *.eps plots_${YR}yrs/
if ($CHFMT == y) then
 mv *.jpg plots_jpg/
 #mv *.$FMT plots_$FMT/
endif
#if ($EPS == y) then
 #mv *.eps plots_eps/
#endif

# Need to first create dirs
#mv *_${YR}yrs.nc netcdf/
#mv $RUNID$a.p*.nc netcdf/
#mv $RUNID$a.da* dumps/

endif # split

#---------------------------------------------------------------------------

echo " "
echo "========================= Finished UMPLOT ========================="
echo " "
echo "You have Set the Parameters of this UMPLOT Run as:"
echo "Jobid:" $RUNID "with files in Directory:" $DIR
echo "Running for" $FullYrs "years with" $SPLIT "5yr Splitting," $BMRK "Extra Plots, and" $CHFMT $FMT "files"
if ($BMRK == y) then
echo "You are Comparing results with" $BDIR
endif

echo " "
echo "==================================================================="
echo " "

echo "Start Time:" $sdate
set edate=`date`
echo "End   Time:" $edate
echo " "

#endif # UMPLT

#if ($UMPLT == y) then

setenv wLog n
if ($wLog == y) then
 echo " "                                                          > log_`date +%d.%m.%y`
 echo "You have Set the Parameters of this UMPLOT Run as:"        >> log_`date +%d.%m.%y`
 echo "Jobid:" $RUNID "with files in Directory:" $DIR             >> log_`date +%d.%m.%y`
 echo "Running for" $FullYrs "years with" $SPLIT "5yr Splitting," $BMRK "Extra Plots, and" $CHFMT $FMT "files" >> log_`date +%d.%m.%y`
 if ($BMRK == y) then
  echo "You are Comparing results with" $BDIR                     >> log_`date +%d.%m.%y`
 endif
 echo " "                                                         >> log_`date +%d.%m.%y`
 echo "Start Time:" $sdate                                        >> log_`date +%d.%m.%y`
 echo " "                                                         >> log_`date +%d.%m.%y`
 echo "  End Time:" $edate                                        >> log_`date +%d.%m.%y`
 echo " "                                                         >> log_`date +%d.%m.%y`
 mv log_`date +%d.%m.%y` log_`date +%d.%m.%y_%H.%M.%S`
endif # wlog

#else
#endif # UMPLT

endif # UMPLT

exit

