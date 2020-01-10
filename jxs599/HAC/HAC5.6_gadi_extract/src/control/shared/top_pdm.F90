! Module containing the variables for Topmodel and PDM.

! DBC Arguably sthzw and zw should be stored in PROGNOSTICS since they
! are indeed prognostics. fsat needs to persist between timesteps (but
! can be initialised (recalculated) from soil moisture).

MODULE top_pdm

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

REAL(KIND=real_jlslsm), ALLOCATABLE :: fexp_soilt(:,:)
  ! Decay factor in Sat. Conductivity in deep LSH/TOPMODEL layer
REAL(KIND=real_jlslsm), ALLOCATABLE :: gamtot_soilt(:,:)
  ! Integrated complete Gamma function
  ! DBC gamtot doesn't need to be in a module in this version, but left there
  !for now for compatability.
REAL(KIND=real_jlslsm), ALLOCATABLE :: ti_mean_soilt(:,:)
  ! Mean topographic index
REAL(KIND=real_jlslsm), ALLOCATABLE :: ti_sig_soilt(:,:)
  ! Standard dev. of topographic index
REAL(KIND=real_jlslsm), ALLOCATABLE :: fsat_soilt(:,:)
  ! Surface saturation fraction
REAL(KIND=real_jlslsm), ALLOCATABLE :: fwetl_soilt(:,:)
  ! Wetland fraction
REAL(KIND=real_jlslsm), ALLOCATABLE :: zw_soilt(:,:)
  ! Water table depth (m)
REAL(KIND=real_jlslsm), ALLOCATABLE :: drain_soilt(:,:)
  ! Drainage out of bottom (nshyd) soil layer (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: dun_roff_soilt(:,:)
  ! Dunne part of sfc runoff (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: qbase_soilt(:,:)
  ! Base flow (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: qbase_zw_soilt(:,:)
  ! Base flow from deep LSH/TOPMODEL layer (kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: fch4_wetl_soilt(:,:)
  ! Scaled wetland methane flux, as
  ! used in atmospheric chemistry.
  ! The substrate is set by parameter ch4_substrate.
  ! (Note different units: 10^-9 kg C/m2/s).
REAL(KIND=real_jlslsm), ALLOCATABLE :: fch4_wetl_cs_soilt(:,:)
  ! Scaled wetland methane flux using
  ! soil carbon as substrate (kg C/m2/s).
REAL(KIND=real_jlslsm), ALLOCATABLE :: fch4_wetl_npp_soilt(:,:)
  ! Scaled wetland methane flux using
  ! NPP as substrate (kg C/m2/s).
REAL(KIND=real_jlslsm), ALLOCATABLE :: fch4_wetl_resps_soilt(:,:)
  ! Scaled wetland methane flux using
  ! soil respiration as substrate (kg C/m2/s).
REAL(KIND=real_jlslsm), ALLOCATABLE :: fch4_wetl_acc_soilt(:,:)
  ! Accum scaled wetland methane flux (kg C/m2)

REAL(KIND=real_jlslsm), ALLOCATABLE :: inlandout_atm_gb(:)
  ! TRIP inland basin outflow (for land points only)(kg/m2/s)
REAL(KIND=real_jlslsm), ALLOCATABLE :: sthzw_soilt(:,:)
  ! soil moist fraction in deep LSH/TOPMODEL layer.
REAL(KIND=real_jlslsm), ALLOCATABLE :: a_fsat_soilt(:,:)
  ! Fitting parameter for Fsat in LSH model
REAL(KIND=real_jlslsm), ALLOCATABLE :: c_fsat_soilt(:,:)
  ! Fitting parameter for Fsat in LSH model
REAL(KIND=real_jlslsm), ALLOCATABLE :: a_fwet_soilt(:,:)
  ! Fitting parameter for Fwet in LSH model
REAL(KIND=real_jlslsm), ALLOCATABLE :: c_fwet_soilt(:,:)
  ! Fitting parameter for Fwet in LSH model

END MODULE top_pdm
