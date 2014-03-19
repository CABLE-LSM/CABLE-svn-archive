#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_enddef(file)

  USE io_constants, ONLY : MODE_WRITE

  USE netcdf, ONLY : nf90_enddef

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes the file out of definition mode - no more dimensions or variables
!   may be defined after this
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(INOUT) :: file

  INTEGER :: error      ! Error indicator


!-----------------------------------------------------------------------------

! We are only genuinely in define mode if the file is open in write mode
  IF ( file%mode == MODE_WRITE ) THEN
    error = nf90_enddef(file%id)
    IF ( error /= NF90_NOERR )                                                &
      CALL log_fatal_ncdf("file_ncdf_enddef",                                 &
                          "Error exiting define mode", error)
  END IF

  RETURN

END SUBROUTINE file_ncdf_enddef
#endif
