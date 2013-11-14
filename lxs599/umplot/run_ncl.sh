#!/bin/csh
# UMPLOT v2.0
# Created by Lauren Stevens, email: lauren.stevens@csiro.au
# 2008-2013

# SET UP - see Chapter 2 in Documentation.
# You must have certain modules loaded before running UMPLOT. E.g.
#module load ncl
#module load nco
#module load cdo
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
 cp ~ste69f/umplot/nml/example_*.nml $HOME/umplot/nml
 echo " "
 echo "(0)     Making Directory $HOME/umplot/nml"
endif

setenv DIRW $PWD

if ( -e $HOME/umplot/nml/umplot.nml ) then
  source $HOME/umplot/nml/umplot.nml
else
 echo " "
 echo "(1)     Namelist File Does Not Exist"
 echo "(1)     Please Create Namelist File in $HOME/umplot/nml"
 echo "(1)     Example Namelist File: ~ste69f/umplot/nml/umplot.nml"
 echo "(1)     and/or Email Lauren.Stevens@csiro.au"
 echo " "
 exit (1)
endif

setenv FullYrs $YR

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

if ($YR == 20 || $YR == 15 || $YR == 10) then
# lxs 4oct13: use mod func ?
 if ($SPLIT == y) then
  if ($YR == 20) then
   setenv block '1 2 3 4'
  else if ($YR == 15) then
   setenv block '1 2 3'
  else if ($YR == 10) then
   setenv block '1 2'
  endif
  setenv YR 5
 endif
else
 setenv SPLIT n
endif

if ($MODEL == c && $RES == 96) then
 setenv TILE 17
 setenv SOIL 6
else
 setenv TILE 9
 setenv SOIL 4
endif

#set date=`date`
if (${?NCI_CP} == 1) then
if ($NCI_CP == y) then
 ~ste69f/umplot/tools/rsync_raijin.sh
else
 echo "(0)     No Fields Files Transferred from Raijin"
 echo ""
endif
else
 echo "(1)     NCI_CP not set. Please update umplot.nml"
 echo "(1)     Contact Lauren.Stevens@csiro.au"
endif
#set date=`date`

# if ~/access doesn't exist ?
if ($DIRW == ~/access) then
 setenv DIRW ~/access/$RUNID
 cd $DIRW
endif

set a=a
set pnc=`ls $DIR/$RUNID$a.p??????.nc | wc -l`

if ($CNV2NC == y) then

  #set date=`date`
   cd $DIR
   ~ste69f/umplot/conv_um2nc.sh
  #set date=`date`
   echo " "
   echo "(0)     Extracting Tmax and Tmin to" $Ptemps "File"
   ~ste69f/umplot/tools/extract_pb.sh
   cd $DIRW

  echo " "
  echo "(0)     Converted Fields Files to NetCDF before Running UMPLOT"

else #CNV2NC == n

 if ( $pnc == 0 ) then
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

set pnc2=`ls $DIR/$RUNID$a.p??????.nc | wc -l`
set chk=`ls -s $DIR/$RUNID$a.p??????.nc`
set i = 0
 while ($i < $pnc2)
  @ k = 2 * $i + 1
  @ m = 2 * $i + 2
   if ( $chk[$k] == "4") then
    echo " "
    echo "(1)     Not all NetCDF Have Been Converted Properly"
    echo " "
    echo "(1)     See File" $chk[$m]
    echo " "
    echo "(1)     If You See Netcdf (i.e. $RUNID$a.p*.nc) with SIZE 32"
    echo "(1)     You'll Need to Delete the UM File(s) and the Netcdf(s)"
    echo "(1)     Or Set CNV2NC = n Once the Netcdf(s) Have been Removed"
    echo "(1)     WARNING: There Could be More than one File with Size 32"
    exit (1)
   endif
  @ i ++
 end

# Lest 12.11.13 will move to umplot.nml
if (${?UMPLT} == 0) then
setenv UMPLT y
else
echo "(1)     NOTE: UMPLT is Set and set to" $UMPLT
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
set date=`date`
echo "Start Time:" $date
echo " "
echo "====================================================================="
echo " "

if ($DIRW == $DIR) then
 echo "(0)     Working directory is the same as Data directory"
else
 echo "(1)     Working directory is NOT the same as Data directory"

 # in $DIRW
 set rid=`ls ?????a.p?????? | head -1 | head -c5`
 setenv RID $rid

 if ( ${#rid} > 0 ) then
  if ( $RID != $RUNID ) then
   echo " "
   echo "(1)     WARNING: Different RUNIDs Present"
   echo "(1)     WARNING: You May Be Running A Different Job"
   echo "(1)     WARNING: Check that you have the Correct NML or DIR"
   exit(1)
  endif
 endif

endif

set pmnc=`ls $DIR/$RUNID$a.$Pmonth?????.nc | wc -l`
 if ( $pmnc == 0 ) then
    echo " "
    echo "(1)     If you get a message like: 'please enter files interactively'"
    echo "(1)     Make Sure Monthly Netcdf Files are in Directory"
    echo "(1)     i.e. $RUNID$a.$Pmonth?????.nc"
 endif

echo " "

# ==================================================================================

if ($SPLIT == y) then

if ( -d $DIRW/block1_5yrs ) then
 echo "(1)     $date - Block Directory Already Exists"
 echo "(1)     Please Delete Directory and Run Again"
 exit (1)
endif

/home/cmar/ste69f/umplot/cdo_merge.sh

foreach blk ( $block )
setenv BLOCK $blk

source ~ste69f/umplot/nml/envars.sh
echo " "
echo "(0)     Post-Processing NetCDF Files"

/home/cmar/ste69f/umplot/split_seasonal.sh

cd $DIRW/block${BLOCK}_5yrs

set flen=`cdo ntime $DIRW/block${BLOCK}_5yrs/Mmonthly_means_${YR}yrs.nc`
@ nom = ($YR * 12)
if ($flen != $nom) then
echo " "
echo "(1)     Number of Months is Incorrect"
echo "(1)     Please Check Length of Mmonthly_means_${YR}yrs.nc"
echo "(1)     It May Be That Not All .p Files are NetCDF"
echo "(1)     "
echo "(1)     If This is the (Only) Issue, 'setenv RUNID' in $DIR"
echo "(1)     and"
echo "(1)     If REINIT=1, Run ~/umplot/umtonc_ss_reinit.sh on Command Line"
echo "(1)     or"
echo "(1)     If REINIT=3, Run ~/umplot/runallv.sh on Command Line"
echo "(1)     Then Run umplot Again"
exit (1)
endif

ncrename -v $tname,tscrn seasonal_means_${YR}yrs.nc
ncrename -v $tname,tscrn monthly_means_${YR}yrs.nc
#ncrename -v $tname,tscrn Mmonthly_means_${YR}yrs.nc
ncrename -v $tname,tscrn yearly_means_${YR}yrs.nc

/home/cmar/ste69f/umplot/process_carbon.sh

echo " "
ncl /home/cmar/ste69f/umplot/run.ncl
ncl /home/cmar/ste69f/umplot/global_means_wbal.ncl
ncl /home/cmar/ste69f/umplot/global_means_tables.ncl
if ($BMRK == y) then
ncl /home/cmar/ste69f/umplot/global_means_ctables.ncl
endif
ncl /home/cmar/ste69f/umplot/carbon_zonal_means.ncl
ncl /home/cmar/ste69f/umplot/global_means_carbon.ncl
#ncl /home/cmar/ste69f/umplot/daily_mean_timeseries.ncl

if( -e $DIRW/block${BLOCK}_5yrs/Timeseries_${YR}yrs.nc ) then

if ($CAL == 360 ) then
 python /home/cmar/ste69f/umplot/MMDC.py $YR
else
 set years=`cdo showyear $DIRW/block${BLOCK}_5yrs/Timeseries_${YR}yrs.nc`
 python /home/cmar/ste69f/umplot/MMDC_365cal.py $YR $years[1] $years[$YR] #[$#years]
endif
# non-Jan start: python $PLOT/MMDC_roll.py $YR
mv MMDC_${YR}yrs.nc MeanMnthDailyCycles_${YR}yrs.nc
python /home/cmar/ste69f/umplot/mmdc_roll.py $YR

ncl /home/cmar/ste69f/umplot/mmdc.ncl

# Figure 6 Jaeger ===================================================

if ($CAL == 360) then
 python /home/cmar/ste69f/umplot/AnnualCycle.py $YR
else
 set years=`cdo showyear $DIRW/block${BLOCK}_5yrs/Timeseries_${YR}yrs.nc`
 python /home/cmar/ste69f/umplot/AnnualCycle_365cal.py $YR $years[1] $years[$YR] #[$#years]
#cdo ymonavg Timeseries_${YR}yrs.nc AnnCycle_${YR}yrs.nc
endif
# non-Jan start: python $PLOT/AnnCycle_roll.py $YR
mv AnnCycle_${YR}yrs.nc AnnualCycle_${YR}yrs.nc

# figure 6 with 6/9 plots - find way to join
ncl /home/cmar/ste69f/umplot/mmdc_jul.ncl
ncl /home/cmar/ste69f/umplot/mmdc_jan.ncl
ncl /home/cmar/ste69f/umplot/annualcycles.ncl

endif

if ($BMRK == y) then
if ( -e ${CABLE}/seasonal_means_${YR}yrs.nc && -e ${MOSES}/seasonal_means_${YR}yrs.nc) then
ncl /home/cmar/ste69f/umplot/Taylor_diagram.ncl
if ( -e ${CABLE}/Timeseries_${YR}yrs.nc && -e ${MOSES}/Timeseries_${YR}yrs.nc) then
ncl /home/cmar/ste69f/umplot/mmdc_hyytiala.ncl
ncl /home/cmar/ste69f/umplot/mmdc_bondville.ncl
ncl /home/cmar/ste69f/umplot/mmdc_hay.ncl
ncl /home/cmar/ste69f/umplot/mmdc_walker_branch.ncl
ncl /home/cmar/ste69f/umplot/mmdc_tharandt.ncl
ncl /home/cmar/ste69f/umplot/mmdc_tumbarumba.ncl
ncl /home/cmar/ste69f/umplot/mmdc_little_washita.ncl
ncl /home/cmar/ste69f/umplot/mmdc_vielsalm.ncl
ncl /home/cmar/ste69f/umplot/mmdc_nsaboreas.ncl
ncl /home/cmar/ste69f/umplot/mmdc_harvard.ncl
ncl /home/cmar/ste69f/umplot/mmdc_loobos.ncl
ncl /home/cmar/ste69f/umplot/mmdc_manaus.ncl
endif
ncl /home/cmar/ste69f/umplot/Mon_IntAnnVar.ncl
ncl /home/cmar/ste69f/umplot/Mon_IntAnnVar_Aust.ncl
if ($YR > 1) then
ncl /home/cmar/ste69f/umplot/Yr_IntAnnVar.ncl
ncl /home/cmar/ste69f/umplot/Yr_IntAnnVar_Aust.ncl
endif
ncl /home/cmar/ste69f/umplot/bias_3panel.ncl
ncl /home/cmar/ste69f/umplot/barchart.ncl
#ncl chv=0 /home/cmar/ste69f/umplot/barchart_amp.ncl
#ncl chv=1 /home/cmar/ste69f/umplot/barchart_amp.ncl
#ncl /home/cmar/ste69f/umplot/barchart_clt.ncl
else
echo "(1)     Cannot plot Fluxnet with Timeseries - Necessary Files Don't Exist"
endif
endif

#rm Timeseries_$YR\yrs.nc Tempseries_$YR\yrs.nc
cd $DIRW 

end

rm h[0123456789].nc i[0123456789].nc j[0123456789].nc k[0123456789].nc

# ==================================================================================

else if ($SPLIT == n) then

#mkdir netcdf
#mkdir dumps

source ~ste69f/umplot/nml/envars.sh

if (! -e $DIRW/seasonal_means_${YR}yrs.nc ) then
echo "(0)     Post-Processing NetCDF Files"
echo " "

/home/cmar/ste69f/umplot/seasonal.sh
/home/cmar/ste69f/umplot/cdo_merge.sh

set flen=`cdo ntime $DIRW/Mmonthly_means_${YR}yrs.nc`
@ nom = ($YR * 12)
if ($flen != $nom) then
echo " "
echo "(1)     Number of Months is Incorrect"
echo "(1)     Please Check Length of Mmonthly_means_${YR}yrs.nc"
echo "(1)     It May Be That Not All .p Files are in NetCDF"
echo "(1)     "
echo "(1)     If This is the (Only) Issue, 'setenv RUNID' in $DIR"
echo "(1)     and"
echo "(1)     If REINIT=1, Run ~/umplot/umtonc_ss_reinit.sh on Command Line"
echo "(1)     or"
echo "(1)     If REINIT=3, Run ~/umplot/runallv.sh on Command Line"
echo "(1)     Then Run umplot Again"
exit (1)
endif

ncrename -v $tname,tscrn seasonal_means_${YR}yrs.nc
ncrename -v $tname,tscrn monthly_means_${YR}yrs.nc
#ncrename -v $tname,tscrn Mmonthly_means_${YR}yrs.nc
ncrename -v $tname,tscrn yearly_means_${YR}yrs.nc

/home/cmar/ste69f/umplot/process_carbon.sh

else
echo " "
echo "(0)     NetCDF Files Already Exist. Will Use These Files to Plot."
echo "(1)     Suggestion: Remove Created _${YR}yrs.nc Files if You Want to Start Again."

endif # if ( ! -e $DIRW/season_${YR}yrs.nc )

echo " "
# Interface bw Shell and NCL - runs major plotting routines
ncl /home/cmar/ste69f/umplot/run.ncl
ncl /home/cmar/ste69f/umplot/global_means_wbal.ncl
ncl /home/cmar/ste69f/umplot/global_means_tables.ncl
if ($BMRK == y) then
ncl /home/cmar/ste69f/umplot/global_means_ctables.ncl
endif
ncl /home/cmar/ste69f/umplot/carbon_zonal_means.ncl
ncl /home/cmar/ste69f/umplot/global_means_carbon.ncl
#ncl /home/cmar/ste69f/umplot/daily_mean_timeseries.ncl

if ( -e $DIRW/Timeseries_${YR}yrs.nc ) then

if ($CAL == 360) then
 python /home/cmar/ste69f/umplot/MMDC.py $YR
else
 set years=`cdo showyear $DIRW/Timeseries_${YR}yrs.nc`
 python /home/cmar/ste69f/umplot/MMDC_365cal.py $YR $years[1] $years[$YR] #[$#years]
endif
# non-Jan start: python $PLOT/MMDC_roll.py $YR
mv MMDC_${YR}yrs.nc MeanMnthDailyCycles_${YR}yrs.nc
python /home/cmar/ste69f/umplot/mmdc_roll.py $YR

ncl /home/cmar/ste69f/umplot/mmdc.ncl

# Figure 6 Jaeger ===================================================

if ($CAL == 360) then
 python /home/cmar/ste69f/umplot/AnnualCycle.py $YR
else
 set years=`cdo showyear $DIRW/Timeseries_${YR}yrs.nc`
 python /home/cmar/ste69f/umplot/AnnualCycle_365cal.py $YR $years[1] $years[$YR] #[$#years]
#cdo ymonavg Timeseries_${YR}yrs.nc AnnCycle_${YR}yrs.nc
endif
#non-Jan start: python $PLOT/AnnCycle_roll.py $YR
mv AnnCycle_${YR}yrs.nc AnnualCycle_${YR}yrs.nc

# figure 6 with 6/9 plots - find way to join
ncl /home/cmar/ste69f/umplot/mmdc_jul.ncl
ncl /home/cmar/ste69f/umplot/mmdc_jan.ncl
ncl /home/cmar/ste69f/umplot/annualcycles.ncl

endif # if Timeseries

if ($BMRK == y) then
if ( -e ${CABLE}/seasonal_means_${YR}yrs.nc && -e ${MOSES}/seasonal_means_${YR}yrs.nc) then
ncl /home/cmar/ste69f/umplot/Taylor_diagram.ncl
if ( -e ${CABLE}/Timeseries_${YR}yrs.nc && -e ${MOSES}/Timeseries_${YR}yrs.nc) then
ncl /home/cmar/ste69f/umplot/mmdc_hyytiala.ncl
ncl /home/cmar/ste69f/umplot/mmdc_bondville.ncl
ncl /home/cmar/ste69f/umplot/mmdc_hay.ncl
ncl /home/cmar/ste69f/umplot/mmdc_walker_branch.ncl
ncl /home/cmar/ste69f/umplot/mmdc_tharandt.ncl
ncl /home/cmar/ste69f/umplot/mmdc_tumbarumba.ncl
ncl /home/cmar/ste69f/umplot/mmdc_little_washita.ncl
ncl /home/cmar/ste69f/umplot/mmdc_vielsalm.ncl
ncl /home/cmar/ste69f/umplot/mmdc_nsaboreas.ncl
ncl /home/cmar/ste69f/umplot/mmdc_harvard.ncl
ncl /home/cmar/ste69f/umplot/mmdc_loobos.ncl
ncl /home/cmar/ste69f/umplot/mmdc_manaus.ncl
endif
ncl /home/cmar/ste69f/umplot/Mon_IntAnnVar.ncl
ncl /home/cmar/ste69f/umplot/Mon_IntAnnVar_Aust.ncl
if ($YR > 1) then
ncl /home/cmar/ste69f/umplot/Yr_IntAnnVar.ncl
ncl /home/cmar/ste69f/umplot/Yr_IntAnnVar_Aust.ncl
endif
ncl /home/cmar/ste69f/umplot/bias_3panel.ncl
ncl /home/cmar/ste69f/umplot/barchart.ncl
#ncl chv=0 /home/cmar/ste69f/umplot/barchart_amp.ncl
#ncl chv=1 /home/cmar/ste69f/umplot/barchart_amp.ncl
#ncl /home/cmar/ste69f/umplot/barchart_clt.ncl
else
echo "(1)     Cannot plot Fluxnet with Timeseries - Necessary Files Don't Exist"
endif
endif

#rm Timeseries_$YR\yrs.nc Tempseries_$YR\yrs.nc

mv *.ps plots/
mv *.eps plots/
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

endif

#---------------------------------------------------------------------------

echo " "
echo "========================= Finished UMPLOT ==========================="
echo " "
echo "You have Set the Parameters of this UMPLOT Run as:"
echo "Jobid:" $RUNID "with files in Directory:" $DIR
echo "Running for" $FullYrs "years with" $SPLIT "5yr Splitting," $BMRK "Extra Plots, and" $CHFMT $FMT "files"
if ($BMRK == y) then
echo "You are Comparing results with" $BDIR
endif

echo " "
echo "====================================================================="
echo " "

echo "Start Time:" $date
set date=`date`
echo "End   Time:" $date
echo " "

endif # UMPLT

exit

