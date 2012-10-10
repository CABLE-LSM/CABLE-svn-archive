#!/bin/csh
# Lauren Stevens 24 Oct 2011
# ============================== DON'T NEED TO MODIFY ==============================
# ==================================================================================

if ($TAY == y) then

if ($SPLIT == y) then
 if ($MODEL == c) then # if you are running CABLE
#  setenv moses $LDIR
  setenv MOSES $LDIR/block${BLOCK}_5yrs
  setenv CABLE $DIRW/block${BLOCK}_5yrs
 else                  # if you are running MOSES
#  setenv cable $LDIR
  setenv MOSES $DIRW/block${BLOCK}_5yrs
  setenv CABLE $LDIR/block${BLOCK}_5yrs
 endif
else
 if ($MODEL == c) then
  setenv MOSES $LDIR
  setenv CABLE $DIRW  # cable for split == n
 else
  setenv MOSES $DIRW  # moses for split == n
  setenv CABLE $LDIR
 endif
endif

setenv cabsea ${CABLE}/seasonal_means_${YR}yrs.nc
setenv mossea ${MOSES}/seasonal_means_${YR}yrs.nc
setenv cabtem ${CABLE}/Tseasonal_means_${YR}yrs.nc
setenv mostem ${MOSES}/Tseasonal_means_${YR}yrs.nc
setenv moshyy ${MOSES}/mmdc_hyytiala_roll_${YR}yrs.nc
setenv cabhyy ${CABLE}/mmdc_hyytiala_roll_${YR}yrs.nc
setenv mostum ${MOSES}/mmdc_tumbarum_roll_${YR}yrs.nc
setenv cabtum ${CABLE}/mmdc_tumbarum_roll_${YR}yrs.nc
setenv moshay ${MOSES}/mmdc_hay_roll_${YR}yrs.nc
setenv cabhay ${CABLE}/mmdc_hay_roll_${YR}yrs.nc
setenv mostha ${MOSES}/mmdc_tharandt_roll_${YR}yrs.nc
setenv cabtha ${CABLE}/mmdc_tharandt_roll_${YR}yrs.nc
setenv mosman ${MOSES}/mmdc_manaus_roll_${YR}yrs.nc
setenv cabman ${CABLE}/mmdc_manaus_roll_${YR}yrs.nc
setenv moslit ${MOSES}/mmdc_litwashi_roll_${YR}yrs.nc
setenv cablit ${CABLE}/mmdc_litwashi_roll_${YR}yrs.nc
setenv moswal ${MOSES}/mmdc_walkerbr_roll_${YR}yrs.nc
setenv cabwal ${CABLE}/mmdc_walkerbr_roll_${YR}yrs.nc
setenv mosbon ${MOSES}/mmdc_bondvill_roll_${YR}yrs.nc
setenv cabbon ${CABLE}/mmdc_bondvill_roll_${YR}yrs.nc
setenv moshar ${MOSES}/mmdc_harvard_roll_${YR}yrs.nc
setenv cabhar ${CABLE}/mmdc_harvard_roll_${YR}yrs.nc
setenv moslob ${MOSES}/mmdc_loobos_roll_${YR}yrs.nc
setenv cablob ${CABLE}/mmdc_loobos_roll_${YR}yrs.nc
setenv mosvie ${MOSES}/mmdc_vielsalm_roll_${YR}yrs.nc
setenv cabvie ${CABLE}/mmdc_vielsalm_roll_${YR}yrs.nc
setenv mosnsa ${MOSES}/mmdc_nsabor_roll_${YR}yrs.nc
setenv cabnsa ${CABLE}/mmdc_nsabor_roll_${YR}yrs.nc

endif

exit

