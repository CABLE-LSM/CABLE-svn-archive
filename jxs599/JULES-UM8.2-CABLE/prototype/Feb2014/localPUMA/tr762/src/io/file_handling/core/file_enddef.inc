#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_enddef(file)

  USE driver_ascii_mod, ONLY : file_ascii_enddef
  USE driver_ncdf_mod, ONLY : file_ncdf_enddef

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
  TYPE(file_handle), INTENT(INOUT) :: file


!-----------------------------------------------------------------------------

  SELECT CASE ( file%driver )
    CASE ( DRIVER_ASCII )
      CALL file_ascii_enddef(file%ascii)

    CASE ( DRIVER_NCDF )
      CALL file_ncdf_enddef(file%ncdf)

    CASE DEFAULT
      CALL log_fatal("file_enddef",                                           &
                     "Unrecognised driver - see docs for available drivers")
  END SELECT

  RETURN

END SUBROUTINE file_enddef
#endif
