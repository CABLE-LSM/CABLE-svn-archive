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

setenv DIRW $PWD

if (! -d $HOME/umplot) then
 mkdir $HOME/umplot
 echo " "
 echo "(0)     Making Directory $HOME/umplot"
endif
if (! -d $HOME/umplot/nml) then
 mkdir $HOME/umplot/nml
 echo " "
 echo "(0)     Making Directory $HOME/umplot/nml"
endif

if ( -e $HOME/umplot/nml/umplot.nml ) then
# if ( -e $HOME/umplot/nml/umplot.nml.$RUNID ) then
#  source $HOME/umplot/nml/umplot.nml.$RUNID
# else
  source $HOME/umplot/nml/umplot.nml
# endif
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

if ($YR == 20 || $YR == 15 || $YR == 10) then
 if ($SPLIT == y || $SPLIT == "") then
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

#echo " "
#if ($DIRW == $DIR) then
# echo "(0)     Working directory is the same as Data directory"
#else
# echo "(1)     Working directory is NOT the same as Data directory"
#endif

set a=a
set pnc=`ls $DIR/$RUNID$a.p*.nc | wc -l`
set pmnc=`ls $DIR/$RUNID$a.pm*.nc | wc -l`
set panc=`ls $DIR/$RUNID$a.pa*.nc | wc -l`
set chk=`ls -s $DIR/$RUNID$a.p*.nc`
if ( $pnc == 0 ) then
 if ($CNV2NC == y || $CNV2NC == "") then
  if ($REINIT == 1) then
   cd $DIR
   ~ste69f/umplot/umtonc_ss_reinit.sh
   cd $DIRW
  else
   ~ste69f/umplot/runallv.sh
    #if ( $pmnc == 0 && $panc != 0 ) then #|| $pmnc < ( $YR * 12 ) && $panc < ( $YR * 12 ) ) then
    if ( ${?Pmonth} == 1 ) then
     if ( ${Pmonth} == pa ) then
      ~ste69f/umplot/tools/pa2pm.sh
     endif
    else
     echo "(1)     Please Set Pmonth - New Environment Variable"
     echo "(1)     Contact Lauren - lauren.stevens@csiro.au"
     exit (1)
    endif
    if ( ${?Pdaily} == 1 ) then
     if ( ${Pdaily} == pe ) then
      ~ste69f/umplot/tools/extract_pb.sh
      #~ste69f/umplot/tools/pe2pb.sh
     endif
    endif
  endif
  echo " "
  echo "(0)     Converted Fields Files to NetCDF before Running UMPLOT"
 else
  echo " "
  echo "(1)     Fields Files Not Converted to NetCDF"
  echo "(1)     Cannot Run UMPLOT until NetCDFs are Present"
  echo "(1)     Either Convert Offline or Switch CNV2NC to y in Namelist"
  echo " "
#if (! -e seasonal_means_5yrs.nc) then
  exit (1)
#endif
 endif
else
# echo " "
# echo "(0)     You Either Manually Converted Files to NetCDF or"
# echo "(0)     Runallv Was Used to Create $RUNID$a.p*.nc Fields Files"
# echo "        OR"
# echo "(1)     NetCDF Haven't Been Converted Properly"
# echo "(1)     If you get a message like: 'please enter files interactively'"
# echo "(1)     Try Removing $RUNID$a.p*.nc and/or Set CNV2NC to y and Run Again"

 set i = 0
 while ($i < $pnc)
  @ k = 2 * $i + 1
   if ( $chk[$k] == "4") then
    echo " "
    echo "(1)     NetCDF Haven't Been Converted Properly"
    echo "(1)     If you get a message like: 'please enter files interactively'"
    echo "(1)     Or If You See Netcdf (i.e. $RUNID$a.p*.nc) with SIZE 32"
    echo "(1)     Try Removing $RUNID$a.p*.nc and Set CNV2NC to y and Run Again"
    exit (1)
   endif
  @ i ++
 end

# if ($BYPS == y) then
# if ($CNV2NC == y || $CNV2NC == "") then
#  if ($REINIT == 1) then
#   cd $DIR
#   ~ste69f/umplot/umtonc_ss_reinit.sh
#   cd $DIRW
#  else
#   ~ste69f/umplot/runallv.sh
#  endif
#  echo " "
#  echo "(0)     Converted Fields Files to NetCDF before Running UMPLOT"
# else
#  echo " "
#  echo "(2)     Fields Files Not Converted to NetCDF"
#  echo "(2)     Cannot Run UMPLOT until NetCDFs are Present"
#  echo "(2)     Either Convert Offline or Switch CNV2NC to y in Namelist"
#  echo " "
#  exit (2)
# endif
# endif

 #if ( $pmnc == 0 && $panc != 0 ) then # || $pmnc < ($YR * 12) && $panc < ($YR * 12) ) then
 if ( ${?Pmonth} == 1 ) then
  if ( ${Pmonth} == pa ) then
   ~ste69f/umplot/tools/pa2pm.sh
  endif
 else
  echo "(1)    Please Set Pmonth - New Environment Variable"
  echo "(1)    Contact Lauren - lauren.stevens@csiro.au"
  exit (1)
 endif

 if ( ${?Pdaily} == 1 ) then
  if ( ${Pdaily} == pe ) then
   ~ste69f/umplot/tools/extract_pb.sh
   #~ste69f/umplot/tools/pe2pb.sh
  endif
 endif

 echo " "
 echo "(0)     You Either Manually Converted UM Files to NetCDF"
 echo "(0)     OR"
 echo "(0)     Runallv Was Used to Create $RUNID$a.p*.nc Fields Files"

endif

echo " "
echo "========================= Running UMPLOT ========================="
echo " "
echo "You have Set the Parameters of this Run as:"
echo "Jobid:" $RUNID "with files in Directory:" $DIR
echo "Running for" $FullYrs "years with" $SPLIT "5yr Splitting," $TAY "Extra Plots, and" $JPEG "Jpeg files"
if ($TAY == y || $TAY == "") then
echo "You are Comparing results with" $LDIR
endif

#set a=a
set date=`date`
echo "Start Time:" $date
echo " "

if ($DIRW == $DIR) then
 echo "(0)     Working directory is the same as Data directory"
else
 echo "(1)     Working directory is NOT the same as Data directory"

 #ncl /home/cmar/ste69f/umplot/get_runid.ncl

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
echo " "

# ==================================================================================

if ($SPLIT == y || $SPLIT == "") then

if ( -d $DIRW/block1_5yrs ) then
 echo "(1)     $date - Block Directory Already Exists"
 echo "(1)     Please Delete Directory and Run Again"
 exit (1)
endif

#/home/cmar/ste69f/umplot/cdo_copy.sh
#/home/cmar/ste69f/umplot/cdo_copy2.sh
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

#cdo chname,temp_26,tscrn ifile ofile
#mv ofile ifile
ncrename -v $tname,tscrn seasonal_means_${YR}yrs.nc
ncrename -v $tname,tscrn monthly_means_${YR}yrs.nc
#ncrename -v $tname,tscrn Mmonthly_means_${YR}yrs.nc
ncrename -v $tname,tscrn yearly_means_${YR}yrs.nc

/home/cmar/ste69f/umplot/process_carbon.sh

echo " "
ncl /home/cmar/ste69f/umplot/run.ncl
ncl /home/cmar/ste69f/umplot/global_means_wbal.ncl
ncl /home/cmar/ste69f/umplot/global_means_tables.ncl
if ($TAY == y) then
ncl /home/cmar/ste69f/umplot/global_means_ctables.ncl
endif
ncl /home/cmar/ste69f/umplot/carbon_zonal_means.ncl
ncl /home/cmar/ste69f/umplot/global_means_carbon.ncl
#qsub -I -lnodes=1:ppn=4,vmem=4gb,walltime=5:00 -d $DIRW/block${BLOCK}_5yrs -x 'source $HOME/umplot/nml/umplot.nml; ncl $PLOT/global_means_carbon.ncl'
ncl /home/cmar/ste69f/umplot/daily_mean_timeseries.ncl

if( -e $DIRW/block${BLOCK}_5yrs/Timeseries_${YR}yrs.nc ) then

if ($CAL == 360 ) then
 python /home/cmar/ste69f/umplot/MMDC.py $YR
else
 set years=`cdo showyear $DIRW/Timeseries_${YR}yrs.nc`
 python /home/cmar/ste69f/umplot/MMDC_365cal.py $YR $years[1] $years[$YR] #[$#years]
endif
# for a start date on Oct 1978
if ($RUNID == xagpb) then
# if the model starts from a month other than Jan
python /home/cmar/ste69f/umplot/MMDC_roll.py $YR
else
mv MMDC_${YR}yrs.nc MeanMnthDailyCycles_${YR}yrs.nc
endif
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
if ($RUNID == xagpb) then
python /home/cmar/ste69f/umplot/AnnCycle_roll.py $YR
else
mv AnnCycle_${YR}yrs.nc AnnualCycle_${YR}yrs.nc
endif

# figure 6 with 6/9 plots - find way to join
ncl /home/cmar/ste69f/umplot/mmdc_jul.ncl
ncl /home/cmar/ste69f/umplot/mmdc_jan.ncl
ncl /home/cmar/ste69f/umplot/annualcycles.ncl

endif

if ($TAY == y || $TAY == "") then
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
#/home/cmar/ste69f/umplot/cdo_copy.sh
#/home/cmar/ste69f/umplot/cdo_copy2.sh
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

#cdo chname,temp_26,tscrn ifile ofile
#mv ofile ifile
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
if ($TAY == y) then
ncl /home/cmar/ste69f/umplot/global_means_ctables.ncl
endif
ncl /home/cmar/ste69f/umplot/carbon_zonal_means.ncl
ncl /home/cmar/ste69f/umplot/global_means_carbon.ncl
#qsub -I -lnodes=1:ppn=4,vmem=4gb,walltime=5:00 -d $DIRW -x 'source $HOME/umplot/nml/umplot.nml; ncl $PLOT/global_means_carbon.ncl'
ncl /home/cmar/ste69f/umplot/daily_mean_timeseries.ncl

if ( -e $DIRW/Timeseries_${YR}yrs.nc ) then

if ($CAL == 360) then
 python /home/cmar/ste69f/umplot/MMDC.py $YR
else
 set years=`cdo showyear $DIRW/Timeseries_${YR}yrs.nc`
 python /home/cmar/ste69f/umplot/MMDC_365cal.py $YR $years[1] $years[$YR] #[$#years]
endif
# for a start date on Oct 1978
if ($RUNID == xagpb) then
# if the model starts from a month other than Jan
python /home/cmar/ste69f/umplot/MMDC_roll.py $YR
else
mv MMDC_${YR}yrs.nc MeanMnthDailyCycles_${YR}yrs.nc
endif
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
if ($RUNID == xagpb) then
python /home/cmar/ste69f/umplot/AnnCycle_roll.py $YR
else
mv AnnCycle_${YR}yrs.nc AnnualCycle_${YR}yrs.nc
endif

# figure 6 with 6/9 plots - find way to join
ncl /home/cmar/ste69f/umplot/mmdc_jul.ncl
ncl /home/cmar/ste69f/umplot/mmdc_jan.ncl
ncl /home/cmar/ste69f/umplot/annualcycles.ncl

endif # if Timeseries

if ($TAY == y || $TAY == "") then
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
if ($JPEG == y || $JPEG == "") then
mv *.jpg plots_jpg/
endif
#if ($EPS == y || $EPS == "") then
#mv *.eps plots_eps/
#endif

# Need to first create dirs
#mv *_${YR}yrs.nc netcdf/
#mv $RUNID$a.p*.nc netcdf/
#mv $RUNID$a.da* dumps/

endif

#--------------------------------------------------------

echo "(0)     Finished Plotting"
echo " "
echo "You have Set the Parameters of this Run as:"
echo "Jobid:" $RUNID "with files in Directory:" $DIR
echo "Running for" $FullYrs "years with" $SPLIT "5yr Splitting," $TAY "Extra Plots, and" $JPEG "Jpeg files"
if ($TAY == y || $TAY == "") then
echo "You are Comparing results with" $LDIR
endif

echo "Start Time:" $date
set date=`date`
echo "End   Time:" $date
echo " "

exit

