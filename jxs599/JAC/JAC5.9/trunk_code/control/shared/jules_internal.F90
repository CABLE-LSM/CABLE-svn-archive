! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE jules_internal

!-----------------------------------------------------------------------------
! Description:
!   Contains internal arrays that need to be passed between different
!   areas of JULES (eg. Surface exchange to snow).
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
      unload_backgrnd_pft(:,:)
                  ! Background unloading rate for canopy (s-1)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_INTERNAL'

CONTAINS

SUBROUTINE jules_internal_alloc(land_pts,npft,                                &
                                cansnowtile)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts,npft

LOGICAL, INTENT(IN) :: cansnowtile(npft)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='JULES_INTERNAL_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ANY(cansnowtile(1:npft)) .EQV. .TRUE.) THEN
  ALLOCATE( unload_backgrnd_pft(land_pts,npft))
ELSE
  ALLOCATE( unload_backgrnd_pft(1,1))
END IF
unload_backgrnd_pft(:,:) = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE jules_internal_alloc

END MODULE jules_internal
