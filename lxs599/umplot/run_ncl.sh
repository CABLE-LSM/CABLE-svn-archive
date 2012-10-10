#!/bin/csh
# Lauren Stevens, email: lauren.stevens@csiro.au
# 2008-2011

#echo " "
#echo "Please Make Sure NCL is Loaded, If You Are Unfamiliar With This"
#echo "Check the Website http://www.ncl.ucar.edu/"
#echo " "

#module load ncl
#if ( ! -e ~/.hluresfile ) then
# echo "You Do Not Have a ~/.hluresfile, Please Copy One to $HOME"
#endif
#module load nco
#module load cdo
#module load python

echo "Please Enter Directory e.g. ~kow014/access/xaank, ~ste69f/access/xagpa :"
setenv DIR $<
setenv DIRW $PWD
echo "Please Enter Run Id e.g. xaank, xaiip, xagpa :"
setenv RUNID $<
echo "Please Enter No. of Years e.g. 1, 3, 5, 10, 15, 20 :"
setenv YR $<
setenv FullYrs $YR

#if ($RUNID == ave*) then
#echo "RUNID is similar to process files ave*.nc"
#echo "Therefore files may be deleted"
#exit (1)
#endif

#Maybe have nml file in local directory ?
if ( -e ~ste69f/umplot/nml/$RUNID.sh ) then
source ~ste69f/umplot/nml/$RUNID.sh
else
echo "(1)     Namelist File Does Not Exist."
echo "(1)     Please Create New Namelist File and/or"
echo "(1)     Email Lauren.Stevens@csiro.au"
echo " "
exit (1)
endif

if ($RES == 96) then
# 17Feb2012 Temporary Question for CABLE n96 - v1 or v2 mask used
echo "Temporary Question: Which Mask ? e.g. 1, 2 :"
setenv MASK $<
endif

if ($YR == 20 || $YR == 15 || $YR == 10) then
 echo "Do you want to split the run into 5 year blocks ? y or n"
 setenv SPLIT $< 
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

echo "Do you want to Compare Runs: Taylor Plot, Barcharts and Fluxnet Plots ? y or n"
setenv TAY $<
if ($TAY == y || $TAY == "") then
 echo "Please Enter Another Model Directory for Comparison e.g. MOSES: ~ste69f/access/xaiyl_m48_20yrs .OR. ~ste69f/access/xagpb_n96_m_20yrs :"
 setenv LDIR $<
endif

if ($SPLIT == n) then
echo "Do you want plots for PowerPoint (jpeg) ? y or n"
setenv JPEG $<
else
setenv JPEG n
endif

#if ($SPLIT == n) then
#if ($JPEG == n) then
#echo "Do you want plots for PowerPoint (eps) ? y or n"
#setenv EPS $<
#else
#setenv EPS n
#endif
#endif

echo " "
if ($DIRW == $DIR) then
 echo "(0)     Working directory is the same as Data directory"
else
 echo "(1)     Working directory is NOT the same as Data directory"
endif

echo " "
echo "You have Set the Parameters of this Run as:"
echo "Jobid:" $RUNID "with files in Directory:" $DIR
echo "Running for" $FullYrs "years with" $SPLIT  "Splitting," $TAY "Extra Plots, and" $JPEG "Jpeg files"
if ($TAY == y || $TAY == "") then
echo "You are Comparing results with" $LDIR
endif

#set a=a
set date=`date`
echo "Start Time:" $date
echo " "

# ==================================================================================

if ($SPLIT == y || $SPLIT == "") then

#if ( ! -e $DIR/$RUNID$a.pm*.nc ) then
#echo "(1)     Please Convert Fields Files to NetCDF before Running"
#echo "(0)     Or Do You Want to Convert Fields Files to NetCDF Now ? y or n"
#setenv CONV2NC $<
#if ($CONV2NC == y || $CONV2NC == "")
#source 
#else
#exit (1)
#endif
#endif

#echo "(0)     Runallv Was Used to Create .nc Fields Files"
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

setenv T15 n
if ($RUNID == xaank) then
#set varnames=`cdo showname seasonal_means_${YR}yrs.nc`
##setenv T15 n
#foreach name ( $varnames )
#if ($name == temp_15 && $T15 == n) then
if ($tname == temp_15) then
setenv T15 y
endif
#end
endif

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

source ~ste69f/umplot/nml/envars.sh

#if ( ! -e $DIR/$RUNID$a.pm*.nc ) then
#echo "(1)     Please Convert Fields Files to NetCDF before Running"
#echo "(0)     Or Do You Want to Convert Fields Files to NetCDF Now ? y or n"
#setenv CONV2NC $<
#if ($CONV2NC == y || $CONV2NC == "")
#source
#else
#exit (1)
#endif
#endif

#echo "(0)     Runallv Was Used to Create .nc Fields Files"
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

setenv T15 n
if ($RUNID == xaank) then
#set varnames=`cdo showname seasonal_means_${YR}yrs.nc`
##setenv T15 n
#foreach name ( $varnames )
#if ($name == temp_15 && $T15 == n) then
if ($tname == temp_15 && $T15 == n) then
setenv T15 y
endif
#end
endif

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

