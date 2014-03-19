#if !defined(UM_JULES)
! Module containing the variables for Topmodel and PDM.

! DBC Arguably sthzw and zw should be stored in PROGNOSTICS since they
! are indeed prognostics. fsat needs to persist between timesteps (but
! can be initialised (recalculated) from soil moisture).


MODULE top_pdm

  IMPLICIT NONE

  REAL, DIMENSION(:), ALLOCATABLE :: FEXP
!                                      ! Decay factor in Sat. Conductivity in water table layer
  REAL, DIMENSION(:), ALLOCATABLE :: GAMTOT
!                                  ! Integrated complete Gamma function
!    DBC gamtot doesn't need to be in a module in this version, but left there for
!        now for compatability.
  REAL, DIMENSION(:), ALLOCATABLE :: TI_MEAN
!                                  ! Mean topographic index
  REAL, DIMENSION(:), ALLOCATABLE :: TI_SIG
!                                  ! Standard dev. of topographic index
  REAL, DIMENSION(:), ALLOCATABLE :: FSAT
!                                  ! Surface saturation fraction
  REAL, DIMENSION(:), ALLOCATABLE :: FWETL
!                                  ! Wetland fraction
  REAL, DIMENSION(:), ALLOCATABLE :: ZW
!                                  ! Water table depth (m)
  REAL, DIMENSION(:), ALLOCATABLE :: DRAIN
!                                  ! Drainage out of bottom (nshyd) soil layer (kg/m2/s)
  REAL, DIMENSION(:), ALLOCATABLE :: DUN_ROFF
!                                  ! Dunne part of sfc runoff (kg/m2/s)
  REAL, DIMENSION(:), ALLOCATABLE :: QBASE
!                                  ! Base flow (kg/m2/s)
  REAL, DIMENSION(:), ALLOCATABLE :: QBASE_ZW
!                                  ! Base flow from ZW layer (kg/m2/s)
  REAL, DIMENSION(:), ALLOCATABLE :: FCH4_WETL
!                                  ! Scaled wetland methane flux (10^-9 kg C/m2/s)

  REAL, ALLOCATABLE :: INLANDOUT_ATM(:)  ! TRIP inland basin outflow (for land points only)(kg/m2/s)
  REAL, ALLOCATABLE :: STHZW(:)     ! soil moist fraction in deep (water table) layer.
  REAL, ALLOCATABLE :: A_FSAT(:)    ! Fitting parameter for Fsat in LSH model
  REAL, ALLOCATABLE :: C_FSAT(:)    ! Fitting parameter for Fsat in LSH model
  REAL, ALLOCATABLE :: A_FWET(:)    ! Fitting parameter for Fwet in LSH model
  REAL, ALLOCATABLE :: C_FWET(:)    ! Fitting parameter for Fwet in LSH model

END MODULE top_pdm
#endif
