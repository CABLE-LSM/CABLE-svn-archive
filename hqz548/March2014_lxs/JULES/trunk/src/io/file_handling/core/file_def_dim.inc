#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_def_dim(file, dim_name, dim_len) RESULT(dim_id)

  USE io_constants, ONLY : MAX_SDF_NAME_LEN

  USE driver_ascii_mod, ONLY : file_ascii_def_dim
  USE driver_ncdf_mod, ONLY : file_ncdf_def_dim

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines a dimension on the given file, returning the dimension id
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_handle), INTENT(INOUT) :: file
                                  ! The file to define the dimension on
  CHARACTER(len=*), INTENT(IN) :: dim_name
                                  ! The name of the dimension
  INTEGER, INTENT(IN) :: dim_len  ! The length of the dimension

! Return type
  INTEGER :: dim_id               ! The dimension id


!-----------------------------------------------------------------------------

  SELECT CASE ( file%driver )
    CASE ( DRIVER_ASCII )
      dim_id = file_ascii_def_dim(file%ascii, dim_name, dim_len)

    CASE ( DRIVER_NCDF )
      dim_id = file_ncdf_def_dim(file%ncdf, dim_name, dim_len)

    CASE DEFAULT
      CALL log_fatal("file_def_dim",                                          &
                     "Unrecognised driver - see docs for available drivers")
  END SELECT

  RETURN

END FUNCTION file_def_dim
#endif
