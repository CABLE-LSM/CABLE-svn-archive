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

END MODULE jules_mod

