#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/control/standalone/spinup/spinup_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

MODULE spinup_mod

  USE model_interface_mod, ONLY : IDENTIFIER_LEN

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
  INTEGER, PARAMETER :: MAX_SPINUP_VARS = 2

!-----------------------------------------------------------------------------
! Module types
!-----------------------------------------------------------------------------
  TYPE spinup_field

    CHARACTER(len=IDENTIFIER_LEN) :: identifier
                   ! The identifier of the model variable associated with this
                   ! spinup field

    LOGICAL :: use_percent = .FALSE.  ! T - tolerance is specified in %
                                      ! F - tolerance is specified as an
                                      !     absolute value
    REAL :: tolerance = 0.0  ! The tolerance to use for this field

    REAL, POINTER :: data(:,:) => NULL()  ! The data from the last comparison

  END TYPE spinup_field

!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
  INTEGER :: nvars = 0
  TYPE(spinup_field), SAVE :: spinup_vars(MAX_SPINUP_VARS)

CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "spinup_init.inc"
#include "spinup_check.inc"

END MODULE spinup_mod
#endif
