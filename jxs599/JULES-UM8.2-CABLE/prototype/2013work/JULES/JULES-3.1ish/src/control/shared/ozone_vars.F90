! File containing variables for ozone implementation

MODULE ozone_vars

#if !defined(UM_JULES)
! Ozone forcing - only required for standalone
  REAL, ALLOCATABLE ::         &
    o3(:)        ! Surface ozone concentration (ppb).
#endif

! Ozone diagnostics
  REAL, ALLOCATABLE ::         &
    flux_o3_ft(:,:)            &
                 ! Flux of O3 to stomata (nmol O3/m2/s).
   ,fo3_ft(:,:)  ! Ozone exposure factor.

END MODULE
