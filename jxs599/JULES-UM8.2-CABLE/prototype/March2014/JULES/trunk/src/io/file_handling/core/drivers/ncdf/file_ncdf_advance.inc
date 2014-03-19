#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_advance(file)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Advances the file by one record
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(INOUT) :: file  ! The file to advance


!-----------------------------------------------------------------------------

  IF ( file%record_dim < 0 )                                                  &
    CALL log_fatal("file_ncdf_advance",                                       &
                   "Cannot advance - no record dimension defined")

  file%current_record = file%current_record + 1

  RETURN

END SUBROUTINE file_ncdf_advance
#endif
