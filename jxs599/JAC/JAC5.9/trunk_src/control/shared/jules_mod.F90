! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  module of declarations and initialisation routines for JULES/UM

MODULE jules_mod

  ! Description:
  !   Module containing declarations of JULES-related variables
  !   and initialisation subroutines.
  !
  ! Method:
  !   Module contains declarations of variables needed for
  !   implementation of JULES routines in the UM, and also
  !   subroutines required to initialise these variables.
  !
  !   Since these variables will be in the UM whether JULES routines
  !   are implemented or not, the idea is to collect them here
  !   and then have different versions of this module for different
  !   physics versions. Specifically, the new variables will have
  !   different dimensions when JULES routines are not used.
  !
  ! Code Owner: Please refer to ModuleLeaders.txt
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !
  ! Declarations:

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!========================
! extra JULES variables
!========================

! snow variable
!--------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  snowdep_surft(:,:)
                  ! Depth of snow (=snowdepth) for all
                  ! surfaces except those using the
                  ! snow canopy, for which it is the
                  ! depth of snow in the canopy (m)
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  albobs_scaling_surft(:,:,:)
                  ! scaling factor applied to match observed
                  ! albedo, within limits, for each tile, for
                  ! SW radn or VIS and NIR
                  ! Allocated size of last dimension depends on l_spec_albedo
                  ! rad_nband is set in check_jules_radiation 
                  ! in jules_radiation_mod and stored in ancil_info

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_MOD'

CONTAINS

SUBROUTINE jules_alloc(land_pts,ntype,nsurft,rad_nband,                       &
                l_albedo_obs)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts,ntype,nsurft,rad_nband

LOGICAL, INTENT(IN) :: l_albedo_obs

!Local variables
INTEGER :: temp_size, temp_tiles, temp_layers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='JULES_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( snowdep_surft(land_pts,nsurft))
snowdep_surft(:,:)          = 0.0

IF ( l_albedo_obs ) THEN
  ALLOCATE( albobs_scaling_surft(land_pts,ntype,rad_nband))
  albobs_scaling_surft(:,:,:) = 0.0
ELSE
  ALLOCATE( albobs_scaling_surft(1,1,1))
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE jules_alloc
END MODULE jules_mod

