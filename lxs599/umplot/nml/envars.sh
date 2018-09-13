#!/bin/csh
# Lauren Stevens 24 Oct 2011
# ============================== DON'T NEED TO MODIFY ==============================
# ==================================================================================

if ($BMRK == y) then

if ($SPLIT == y) then
 if ($MODEL == c) then # if you are running CABLE
#  setenv moses $BDIR
  setenv TOP $BDIR/block${BLOCK}_5yrs
  setenv MID $DIRW/block${BLOCK}_5yrs
 else                  # if you are running MOSES
#  setenv cable $BDIR
  setenv TOP $DIRW/block${BLOCK}_5yrs
  setenv MID $BDIR/block${BLOCK}_5yrs
 endif
else
 if ($MODEL == c) then
  setenv TOP $BDIR
  setenv MID $DIRW  # cable for split == n
 else
  setenv TOP $DIRW  # moses for split == n
  setenv MID $BDIR
 endif
endif

setenv cabsea ${MID}/seasonal_means_${YR}yrs.nc
setenv cabtem ${MID}/Tseasonal_means_${YR}yrs.nc
setenv cabhyy ${MID}/mmdc_hyytiala_roll_${YR}yrs.nc
setenv cabtum ${MID}/mmdc_tumbarum_roll_${YR}yrs.nc
setenv cabhay ${MID}/mmdc_hay_roll_${YR}yrs.nc
setenv cabtha ${MID}/mmdc_tharandt_roll_${YR}yrs.nc
setenv cabman ${MID}/mmdc_manaus_roll_${YR}yrs.nc
setenv cablit ${MID}/mmdc_litwashi_roll_${YR}yrs.nc
setenv cabwal ${MID}/mmdc_walkerbr_roll_${YR}yrs.nc
setenv cabbon ${MID}/mmdc_bondvill_roll_${YR}yrs.nc
setenv cabhar ${MID}/mmdc_harvard_roll_${YR}yrs.nc
setenv cablob ${MID}/mmdc_loobos_roll_${YR}yrs.nc
setenv cabvie ${MID}/mmdc_vielsalm_roll_${YR}yrs.nc
setenv cabnsa ${MID}/mmdc_nsabor_roll_${YR}yrs.nc

setenv mossea ${TOP}/seasonal_means_${YR}yrs.nc
setenv mostem ${TOP}/Tseasonal_means_${YR}yrs.nc
setenv moshyy ${TOP}/mmdc_hyytiala_roll_${YR}yrs.nc
setenv mostum ${TOP}/mmdc_tumbarum_roll_${YR}yrs.nc
setenv moshay ${TOP}/mmdc_hay_roll_${YR}yrs.nc
setenv mostha ${TOP}/mmdc_tharandt_roll_${YR}yrs.nc
setenv mosman ${TOP}/mmdc_manaus_roll_${YR}yrs.nc
setenv moslit ${TOP}/mmdc_litwashi_roll_${YR}yrs.nc
setenv moswal ${TOP}/mmdc_walkerbr_roll_${YR}yrs.nc
setenv mosbon ${TOP}/mmdc_bondvill_roll_${YR}yrs.nc
setenv moshar ${TOP}/mmdc_harvard_roll_${YR}yrs.nc
setenv moslob ${TOP}/mmdc_loobos_roll_${YR}yrs.nc
setenv mosvie ${TOP}/mmdc_vielsalm_roll_${YR}yrs.nc
setenv mosnsa ${TOP}/mmdc_nsabor_roll_${YR}yrs.nc

endif

exit

