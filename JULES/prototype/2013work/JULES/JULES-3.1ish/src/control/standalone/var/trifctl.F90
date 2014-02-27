#if !defined(UM_JULES)
!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/control/standalone/var/trifctl.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

MODULE trifctl

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module containing the variables for TRIFFID and plant phenology
! (not parameter values)
!-----------------------------------------------------------------------------

  INTEGER ::                                                                  &
    ASTEPS_SINCE_TRIFFID,                                                     &
      ! Number of atmospheric timesteps since last call to TRIFFID
    PHENOL_PERIOD,                                                            &
      ! Phenology period (days)
    TRIFFID_PERIOD
      ! TRIFFID period (days)

  REAL, DIMENSION(:,:), ALLOCATABLE :: G_LEAF_ACC
!                                    ! Accumulated leaf turnover rate
  REAL, DIMENSION(:,:), ALLOCATABLE :: NPP_FT_ACC
!                                    ! Accumulated NPP_FT
  REAL, DIMENSION(:,:), ALLOCATABLE :: G_LEAF_PHEN_ACC
!                                    ! Accumulated leaf turnover rate including phenology
  REAL, DIMENSION(:,:), ALLOCATABLE :: RESP_W_FT_ACC
!                                    ! Accum RESP_W_FT
  REAL, DIMENSION(:,:), ALLOCATABLE :: RESP_S_ACC
!                                    ! Accumulated RESP_S
  REAL, DIMENSION(:),   ALLOCATABLE :: GPP
!                                    ! Gross primary productivity (kg C/m2/s)
  REAL, DIMENSION(:),   ALLOCATABLE :: NPP
!                                    ! Net primary productivity (kg C/m2/s)
  REAL, DIMENSION(:),   ALLOCATABLE :: RESP_P
!                                    ! Plant respiration (kg C/m2/s)
  REAL, DIMENSION(:,:), ALLOCATABLE :: G_LEAF
!                                    ! Leaf turnover rate (/360days)
  REAL, DIMENSION(:,:), ALLOCATABLE :: G_LEAF_PHEN
!                                    ! Mean leaf turnover rate over phenology period(/360days)
  REAL, DIMENSION(:,:), ALLOCATABLE :: GPP_FT
!                                    ! Gross primary productivity on PFTs (kg C/m2/s)
  REAL, DIMENSION(:,:), ALLOCATABLE :: NPP_FT
!                                    ! Net primary productivity on PFTs (kg C/m2/s)
  REAL, DIMENSION(:,:), ALLOCATABLE :: RESP_P_FT
!                                    ! Plant respiration on PFTs (kg C/m2/s)
  REAL, DIMENSION(:,:), ALLOCATABLE :: RESP_S
!                                    ! Soil respiration (kg C/m2/s)
  REAL, DIMENSION(:,:), ALLOCATABLE :: RESP_W_FT
!                                    ! Wood maintenance respiration (kg C/m2/s)
  REAL, DIMENSION(:,:), ALLOCATABLE :: LAI_PHEN
!                                    ! LAI of PFTs after phenology.
!                                    !   Required as separate variable for top-level argument list matching with VEG_IC2A
  REAL, DIMENSION(:,:), ALLOCATABLE :: C_VEG
!                                    ! Total carbon content of the vegetation (kg C/m2)
  REAL, DIMENSION(:),   ALLOCATABLE :: CV
!                                    ! Gridbox mean vegetation carbon (kg C/m2)
  REAL, DIMENSION(:,:), ALLOCATABLE :: G_LEAF_DAY
!                                    ! Mean leaf turnover rate for input to PHENOL (/360days)
  REAL, DIMENSION(:,:), ALLOCATABLE :: G_LEAF_DR_OUT
!                                    ! Mean leaf turnover rate for driving TRIFFID (/360days)
  REAL, DIMENSION(:,:), ALLOCATABLE :: LIT_C
!                                    ! Carbon Litter (kg C/m2/360days)
  REAL, DIMENSION(:),   ALLOCATABLE :: LIT_C_MN
!                                    ! Gridbox mean carbon litter (kg C/m2/360days)
  REAL, DIMENSION(:,:), ALLOCATABLE :: NPP_DR_OUT
!                                    ! Mean NPP for driving TRIFFID (kg C/m2/360days)
  REAL, DIMENSION(:,:), ALLOCATABLE :: RESP_W_DR_OUT
!                                    ! Mean wood respiration for driving TRIFFID (kg C/m2/360days)
  REAL, DIMENSION(:,:), ALLOCATABLE :: RESP_S_DR_OUT
!                                    ! Mean soil respiration for driving TRIFFID (kg C/m2/360days)
  REAL, DIMENSION(:),   ALLOCATABLE :: FRAC_AGR
!                                    ! Fraction of agriculture

END MODULE trifctl
#endif
