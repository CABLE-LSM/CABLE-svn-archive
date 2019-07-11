#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_inquire_var(file, var_name, var_id, ndims, dim_ids, is_record)

  USE io_constants, ONLY : MAX_DIM_FILE

  USE driver_ascii_mod, ONLY : file_ascii_inquire_var
  USE driver_ncdf_mod, ONLY : file_ncdf_inquire_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given a variable name, returns its id and information about dimensions
!   If the returned id < 0, the variable doesn't exist
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_handle), INTENT(IN) :: file  ! The file to inspect
  CHARACTER(len=*), INTENT(IN) :: var_name  ! The name of the variable to
                                            ! inquire about
  INTEGER, INTENT(OUT) :: var_id  ! The id of the variable in the file
  INTEGER, INTENT(OUT), OPTIONAL :: ndims  ! The number of dimensions that
                                           ! the variable has
  INTEGER, INTENT(OUT), OPTIONAL :: dim_ids(:)  ! The dimension ids
  LOGICAL, INTENT(OUT), OPTIONAL :: is_record  ! Indicates if the variable
                                               ! uses the record dimension


! Work variables
! Local versions of the optional arguments to pass to the underlying
! implementations
  INTEGER :: ndims_local
  INTEGER :: dim_ids_local(MAX_DIM_FILE)
  LOGICAL :: is_record_local


!-----------------------------------------------------------------------------

  SELECT CASE ( file%driver )
    CASE ( DRIVER_ASCII )
      CALL file_ascii_inquire_var(                                            &
        file%ascii, var_name, var_id,                                         &
        ndims_local, dim_ids_local, is_record_local                           &
      )

    CASE ( DRIVER_NCDF )
      CALL file_ncdf_inquire_var(                                             &
        file%ncdf, var_name, var_id,                                          &
        ndims_local, dim_ids_local, is_record_local                           &
      )

    CASE DEFAULT
      CALL log_fatal("file_inquire_var",                                      &
                     "Unrecognised driver - see docs for available drivers")
  END SELECT

! Copy the returned values into the optional arguments if present
  IF ( PRESENT(ndims) ) ndims = ndims_local
  IF ( PRESENT(dim_ids) ) dim_ids(1:ndims) = dim_ids_local(1:ndims)
  IF ( PRESENT(is_record) ) is_record = is_record_local

  RETURN

END SUBROUTINE file_inquire_var
#endif
