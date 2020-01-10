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

END MODULE jules_internal
