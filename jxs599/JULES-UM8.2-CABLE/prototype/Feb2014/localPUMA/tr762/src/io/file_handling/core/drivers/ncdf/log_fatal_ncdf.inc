#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE log_fatal_ncdf(originating_proc, message, ncdf_err)

  USE netcdf, ONLY : nf90_strerror

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Logs a NetCDF error
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  CHARACTER(len=*), INTENT(IN) :: originating_proc  ! The procedure where
                                                    ! the error occured
  CHARACTER(len=*), INTENT(IN) :: message  ! Message to accompany error
  INTEGER, INTENT(IN) :: ncdf_err  ! The NetCDF error code


!-----------------------------------------------------------------------------

  CALL log_fatal(originating_proc, TRIM(message) //                           &
                 " (NetCDF error - " // TRIM(nf90_strerror(ncdf_err)) // ")")

  RETURN

END SUBROUTINE log_fatal_ncdf
#endif
