#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_close(file)

  USE driver_ascii_mod, ONLY : file_ascii_close
  USE driver_ncdf_mod, ONLY : file_ncdf_close

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Closes and frees any resources consumed by the given file object
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_handle), INTENT(INOUT) :: file ! The file to close


!-----------------------------------------------------------------------------

  SELECT CASE ( file%driver )
    CASE ( DRIVER_ASCII )
      CALL file_ascii_close(file%ascii)

    CASE ( DRIVER_NCDF )
      CALL file_ncdf_close(file%ncdf)

    CASE DEFAULT
      CALL log_fatal("file_close",                                            &
                     "Unrecognised driver - see docs for available drivers")
  END SELECT

  RETURN

END SUBROUTINE file_close
#endif
