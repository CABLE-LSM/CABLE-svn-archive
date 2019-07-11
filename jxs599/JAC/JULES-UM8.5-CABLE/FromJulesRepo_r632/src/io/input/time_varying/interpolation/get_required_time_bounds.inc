#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE get_required_time_bounds(interp_flags, lower_bound, upper_bound)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given a list of interpolation flags, returns the common lower bound and
!   upper bound required to use all the given interpolation schemes in the
!   same file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  CHARACTER(len=*), INTENT(IN) :: interp_flags(:)

! Return types
  INTEGER, INTENT(OUT) :: lower_bound, upper_bound


! Work variables
  LOGICAL :: required_times(-1:2)  ! Array indicating the required times
  INTEGER :: i  ! Loop variable

!-----------------------------------------------------------------------------

  required_times(:) = .FALSE.

  DO i = 1,SIZE(interp_flags)

! If all the times are already required, we can exit early
    IF ( ALL(required_times) ) EXIT

! Otherwise, set the required times based on the interpolation flag
    SELECT CASE ( interp_flags(i) )
      CASE ( INTERP_AVE_BACKWARD )
        required_times(0:2) = .TRUE.

      CASE ( INTERP_AVE_CENTRED )
        required_times(-1:2) = .TRUE.

      CASE ( INTERP_AVE_FORWARD )
        required_times(-1:1) = .TRUE.

      CASE ( INTERP_INSTANT )
        required_times(0:1) = .TRUE.

      CASE ( NO_INTERP_END )
        required_times(1) = .TRUE.

      CASE ( NO_INTERP_CENTRED )
        required_times(0:1) = .TRUE.

      CASE ( NO_INTERP_START )
        required_times(0) = .TRUE.

      CASE DEFAULT
        CALL log_fatal("get_required_times",                                  &
                       "Unrecognised interpolation flag - " // interp_flags(i))
    END SELECT

  END DO

! The required lower bound is the first element for which required_times=T
  DO i = -1,2
    IF ( required_times(i) ) THEN
      lower_bound = i
      EXIT
    END IF
  END DO

! The required upper bound is the last element for which required_times=T
! So step backwards through the required_times array to detect it
  DO i = 2,-1,-1
    IF ( required_times(i) ) THEN
      upper_bound = i
      EXIT
    END IF
  END DO

  RETURN

END SUBROUTINE get_required_time_bounds
#endif
