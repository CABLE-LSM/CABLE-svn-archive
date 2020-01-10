! File containing variables for ozone implementation

MODULE ozone_vars

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Ozone forcing
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  o3_gb(:)        ! Surface ozone concentration (ppb).


! Ozone diagnostics
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  flux_o3_pft(:,:)                                                            &
               ! Flux of O3 to stomata (nmol O3/m2/s).
 ,fo3_pft(:,:)  ! Ozone exposure factor.

END MODULE
