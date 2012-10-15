#!/bin/csh
# UMPLOT
# Created by Lauren Stevens, email: lauren.stevens@csiro.au
# 2008-2012

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
 echo "(1)     It May Mean You Are Unfamiliar with NCL"
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
 echo "(0)     Making Directory $HOME/umplot"
endif
if (! -d $HOME/umplot/nml) then
 mkdir $HOME/umplot/nml
 echo "(0)     Making Directory $HOME/umplot/nml"
endif

if ( -e $HOME/umplot/nml/umplot.nml ) then
 source $HOME/umplot/nml/umplot.nml
else
 echo "(1)     Namelist File Does Not Exist"
 echo "(1)     Please Create Namelist File in $HOME/umplot/nml"
 echo "(1)     Example Namelist File: ~ste69f/umplot/nml/umplot.nml"
 echo "(1)     and/or Email Lauren.Stevens@csiro.au"
 echo " "
 exit (1)
endif

setenv FullYrs $YR

#set rfile=`ls $DIR/ave*a.p*.nc | wc -l`
#if ($rfile > 0) then
# echo "(1)     RUNID is similar to processed files ave*.nc"
# echo "(1)     WARNING: Therefore necessary files may be deleted"
# exit (1)
#endif

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

echo " "
echo "You have Set the Parameters of this Run as:"
echo "Jobid:" $RUNID "with files in Directory:" $DIR
echo "Running for" $FullYrs "years with" $SPLIT  "Splitting," $TAY "Extra Plots, and" $JPEG "Jpeg files"
if ($TAY == y || $TAY == "") then
echo "You are Comparing results with" $LDIR
endif

set a=a
set date=`date`
echo "Start Time:" $date
echo " "

if ($DIRW == $DIR) then
 echo "(0)     Working directory is the same as Data directory"
else
 echo "(1)     Working directory is NOT the same as Data directory"
endif
echo " "

#set pnc=`ls $DIR/$RUNID$a.p*.nc | wc -l`
#if ( $pnc == 0 ) then
# if ($CNV2NC == y || $CNV2NC == "") then
#  echo "(0)     Converting Fields Files to NetCDF before Running"
#  if ($REINIT == 1)
#   $PLOT/new*/umtonc_ss* 
#  else
#   source # runallv ?
#  endif
# else
#  echo "(1)     Fields Files Not Converted to NetCDF"
#  echo "(1)     Cannot Run UMPLOT until NetCDFs are Present"
#  echo "(1)     Either Convert Offline or Switch CNV2NC to y in Namelist"
#  echo " "
#  exit (1)
# endif
#else
# echo "(0)     You Either Manually Converted to NetCDF or"
# echo "(0)     Runallv Was Used to Create $RUNID$a.p*.nc Fields Files"
# echo " "
#endif

# ==================================================================================

if ($SPLIT == y || $SPLIT == "") then

#/home/cmar/ste69f/umplot/cdo_copy.sh
/home/cmar/ste69f/umplot/cdo_copy2.sh

foreach blk ( $block )
setenv BLOCK $blk

source ~ste69f/umplot/nml/envars.sh
echo " "
echo "(0)     Creating NetCDF Files"

/home/cmar/ste69f/umplot/split_seasonal.sh

cd $DIRW/block${BLOCK}_5yrs

ncrename -v $tname,tscrn seasonal_means_${YR}yrs.nc
ncrename -v $tname,tscrn monthly_means_${YR}yrs.nc
ncrename -v $tname,tscrn Mmonthly_means_${YR}yrs.nc
ncrename -v $tname,tscrn yearly_means_${YR}yrs.nc

/home/cmar/ste69f/umplot/process_carbon.sh

echo " "
ncl /home/cmar/ste69f/umplot/run.ncl
ncl /home/cmar/ste69f/umplot/global_means_wbal.ncl
ncl /home/cmar/ste69f/umplot/carbon_zonal_means.ncl
ncl /home/cmar/ste69f/umplot/global_means_carbon.ncl

if( -e $DIRW/block${BLOCK}_5yrs/Timeseries_${YR}yrs.nc ) then

python /home/cmar/ste69f/umplot/MMDC.py $YR
# for a start date on Oct 1978
if ($RUNID == xagpb) then
python /home/cmar/ste69f/umplot/MMDC_roll.py $YR
else
mv MMDC_${YR}yrs.nc MeanMnthDailyCycles_${YR}yrs.nc
endif
python /home/cmar/ste69f/umplot/mmdc_roll.py $YR

ncl /home/cmar/ste69f/umplot/mmdc.ncl

# Figure 6 Jaeger ===================================================

python /home/cmar/ste69f/umplot/AnnualCycle.py $YR
#cdo ymonavg Timeseries_${YR}yrs.nc AnnualCycle_${YR}yrs.nc
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
if ( -e $CABLE/seasonal_means_${YR}yrs.nc && -e ${MOSES}/seasonal_means_${YR}yrs.nc) then
ncl /home/u/cmar/ste69f/umplot/Taylor_diagram.ncl
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
ncl /home/cmar/ste69f/umplot/Mon_IntAnnVar.ncl
ncl /home/cmar/ste69f/umplot/Yr_IntAnnVar.ncl
ncl /home/cmar/ste69f/umplot/Mon_IntAnnVar_Aust.ncl
ncl /home/cmar/ste69f/umplot/Yr_IntAnnVar_Aust.ncl
ncl /home/cmar/ste69f/umplot/bias_3panel.ncl
ncl /home/cmar/ste69f/umplot/barchart.ncl
else
echo "(1)     Cannot plot Fluxnet with Timeseries - Files Don't Exist"
endif
endif

#rm Timeseries_$YR\yrs.nc Tempseries_$YR\yrs.nc
cd $DIRW 

end

rm h[0123456789]*.nc i[0123456789]*.nc j[0123456789]*.nc k[0123456789]*.nc

# ==================================================================================

else if ($SPLIT == n) then

#mkdir netcdf
#mkdir dumps

source ~ste69f/umplot/nml/envars.sh

if (! -e $DIRW/seasonal_means_${YR}yrs.nc ) then
echo "(0)     Creating NetCDF Files"
echo " "

/home/cmar/ste69f/umplot/seasonal.sh
#/home/cmar/ste69f/umplot/cdo_copy.sh
/home/cmar/ste69f/umplot/cdo_copy2.sh

ncrename -v $tname,tscrn seasonal_means_${YR}yrs.nc
ncrename -v $tname,tscrn monthly_means_${YR}yrs.nc
ncrename -v $tname,tscrn Mmonthly_means_${YR}yrs.nc
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
ncl /home/cmar/ste69f/umplot/carbon_zonal_means.ncl
ncl /home/cmar/ste69f/umplot/global_means_carbon.ncl

if ( -e $DIRW/Timeseries_${YR}yrs.nc ) then

python /home/cmar/ste69f/umplot/MMDC.py $YR
# for a start date on Oct 1978
if ($RUNID == xagpb) then
python /home/cmar/ste69f/umplot/MMDC_roll.py $YR
else
mv MMDC_${YR}yrs.nc MeanMnthDailyCycles_${YR}yrs.nc
endif
python /home/cmar/ste69f/umplot/mmdc_roll.py $YR

ncl /home/cmar/ste69f/umplot/mmdc.ncl

# Figure 6 Jaeger ===================================================

python /home/cmar/ste69f/umplot/AnnualCycle.py $YR
#cdo ymonavg Timeseries_${YR}yrs.nc AnnualCycle_${YR}yrs.nc
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
ncl /home/cmar/ste69f/umplot/Mon_IntAnnVar.ncl
ncl /home/cmar/ste69f/umplot/Yr_IntAnnVar.ncl
ncl /home/cmar/ste69f/umplot/Mon_IntAnnVar_Aust.ncl
ncl /home/cmar/ste69f/umplot/Yr_IntAnnVar_Aust.ncl
ncl /home/cmar/ste69f/umplot/bias_3panel.ncl
ncl /home/cmar/ste69f/umplot/barchart.ncl
else
echo "(1)     Cannot plot Fluxnet with Timeseries - Files Don't Exist"
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

echo "(0)     Finished Plotting"
echo " "
echo "You have Set the Parameters of this Run as:"
echo "Jobid:" $RUNID "with files in Directory:" $DIR
echo "Running for" $FullYrs "years with" $SPLIT  "Splitting," $TAY "Extra Plots, and" $JPEG "Jpeg files"
if ($TAY == y || $TAY == "") then
echo "You are Comparing results with" $LDIR
endif

echo "Start Time:" $date
set date=`date`
echo "End   Time:" $date
echo " "

exit

