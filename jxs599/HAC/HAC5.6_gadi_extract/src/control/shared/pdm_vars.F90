! Module containing the variables for PDM.

MODULE pdm_vars

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

REAL(KIND=real_jlslsm), ALLOCATABLE :: slope_gb(:)
  ! Terrain slope: to be used in spatial varying b_pdm calculation

END MODULE pdm_vars
