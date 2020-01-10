MODULE trifctl

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module containing the variables for TRIFFID and plant phenology
! (not parameter values)
!-----------------------------------------------------------------------------

INTEGER :: asteps_since_triffid
    ! Number of atmospheric timesteps since last call to TRIFFID

REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: g_leaf_acc_pft
  ! Accumulated leaf turnover rate (/360days).
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: npp_acc_pft
  ! Accumulated NPP_FT (kg m-2).
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: g_leaf_phen_acc_pft
  ! Accumulated leaf turnover rate including phenology (/360days).
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: resp_w_acc_pft
  ! Accumulated RESP_W_FT (kg m-2).
REAL(KIND=real_jlslsm), DIMENSION(:,:,:,:), ALLOCATABLE :: resp_s_acc_soilt
  ! Accumulated RESP_S (kg m-2).

REAL(KIND=real_jlslsm), DIMENSION(:),   ALLOCATABLE :: gpp_gb
  ! Gross primary productivity (kg C/m2/s)
REAL(KIND=real_jlslsm), DIMENSION(:),   ALLOCATABLE :: npp_gb
  ! Net primary productivity (kg C/m2/s)
REAL(KIND=real_jlslsm), DIMENSION(:),   ALLOCATABLE :: resp_p_gb
  ! Plant respiration (kg C/m2/s)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: g_leaf_pft
  ! Leaf turnover rate (/360days)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: g_leaf_phen_pft
  ! Mean leaf turnover rate over phenology period (/360days)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: gpp_pft
  ! Gross primary productivity on PFTs (kg C/m2/s)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: npp_pft
  ! Net primary productivity on PFTs (kg C/m2/s)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: resp_p_pft
  ! Plant respiration on PFTs (kg C/m2/s)
REAL(KIND=real_jlslsm), DIMENSION(:,:,:,:), ALLOCATABLE :: resp_s_soilt
  ! Soil respiration (kg C/m2/s)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: resp_w_pft
  ! Wood maintenance respiration (kg C/m2/s)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: lai_phen_pft
  ! LAI of PFTs after phenology.
  ! Required as separate variable for top-level argument list matching
  !with VEG_IC2A
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: c_veg_pft
  ! Total carbon content of the vegetation (kg C/m2)
REAL(KIND=real_jlslsm), DIMENSION(:),   ALLOCATABLE :: cv_gb
  ! Gridbox mean vegetation carbon (kg C/m2)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: g_leaf_day_pft
  ! Mean leaf turnover rate for input to PHENOL (/360days)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: g_leaf_dr_out_pft
  ! Mean leaf turnover rate for driving TRIFFID (/360days)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: lit_c_pft
  ! Carbon Litter (kg C/m2/360days)
REAL(KIND=real_jlslsm), DIMENSION(:),   ALLOCATABLE :: lit_c_mn_gb
  ! Gridbox mean carbon litter (kg C/m2/360days)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: npp_dr_out_pft
  ! Mean NPP for driving TRIFFID (kg C/m2/360days)
REAL(KIND=real_jlslsm), DIMENSION(:,:), ALLOCATABLE :: resp_w_dr_out_pft
  ! Mean wood respiration for driving TRIFFID (kg C/m2/360days)
REAL(KIND=real_jlslsm), DIMENSION(:,:,:), ALLOCATABLE :: resp_s_dr_out_gb
  ! Mean soil respiration for driving TRIFFID (kg C/m2/360days)
REAL(KIND=real_jlslsm), DIMENSION(:),   ALLOCATABLE :: frac_agr_gb
  ! Fraction of agriculture

END MODULE trifctl
