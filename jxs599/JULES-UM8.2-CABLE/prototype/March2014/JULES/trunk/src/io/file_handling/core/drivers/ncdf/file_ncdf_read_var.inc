#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ncdf_read_var_scalar(file, var_id, value, start)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE netcdf, ONLY : nf90_get_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a scalar value from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, INTENT(OUT) :: value  ! The value read from file
  INTEGER, INTENT(IN) :: start(MAX_DIM_VAR)
                              ! The point to start reading from
                              ! This should contain one value for each
                              ! NON-RECORD dimension of the variable, and
                              ! any unused slots should be one

! Work variables
  TYPE(var_ncdf) :: var  ! Structure containing information about the
                         ! variable

  INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                         ! has, including record dimension

  INTEGER, ALLOCATABLE :: local_start(:)

  INTEGER :: error      ! Error indicator


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
  var = file%vars(var_id)

  ndims_tot = var%ndims
  IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
  ALLOCATE(local_start(ndims_tot))
  local_start(1:var%ndims) = start(1:var%ndims)
  IF ( var%is_record ) local_start(ndims_tot) = file%current_record

! Try to read the variable
  error = nf90_get_var(file%id, var%id, value, local_start)
  IF ( error /= NF90_NOERR )                                                  &
    CALL log_fatal_ncdf("file_ncdf_read_var_scalar",                          &
                        "Error reading variable", error)

  DEALLOCATE(local_start)

  RETURN

END SUBROUTINE file_ncdf_read_var_scalar


SUBROUTINE file_ncdf_read_var_1d(file, var_id, values, start, count)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE netcdf, ONLY : nf90_get_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 1d array from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, DIMENSION(:), INTENT(OUT) :: values  ! The values read from file
  INTEGER, INTENT(IN) :: start(MAX_DIM_VAR)  ! The point to start reading from
                                   ! (one value for each NON-RECORD
                                   ! dimension of the variable)
  INTEGER, INTENT(IN) :: count(MAX_DIM_VAR)  ! The number of points to read
                                   ! in each dimension of the variable


! Work variables
  TYPE(var_ncdf) :: var  ! Structure containing information about the
                         ! variable

  INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                         ! has, including record dimension

  INTEGER, ALLOCATABLE :: local_start(:)
  INTEGER, ALLOCATABLE :: local_count(:)

  INTEGER :: error      ! Error indicator


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
  var = file%vars(var_id)

  ndims_tot = var%ndims
  IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
  ALLOCATE(local_start(ndims_tot))
  local_start(1:var%ndims) = start(1:var%ndims)
  IF ( var%is_record ) local_start(ndims_tot) = file%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given count array, plus a value so that
! we read from only the correct record, if required
  ALLOCATE(local_count(ndims_tot))
  local_count(1:var%ndims) = count(1:var%ndims)
  IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
  error = nf90_get_var(file%id, var%id, values, local_start, local_count)
  IF ( error /= NF90_NOERR )                                                  &
    CALL log_fatal_ncdf("file_ncdf_read_var_1d",                              &
                        "Error reading variable", error)

  DEALLOCATE(local_start)
  DEALLOCATE(local_count)

  RETURN

END SUBROUTINE file_ncdf_read_var_1d


SUBROUTINE file_ncdf_read_var_2d(file, var_id, values, start, count)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE netcdf, ONLY : nf90_get_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 2d array from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, DIMENSION(:,:), INTENT(OUT) :: values  ! The values read from file
  INTEGER, INTENT(IN) :: start(MAX_DIM_VAR)  ! The point to start reading from
                                   ! (one value for each NON-RECORD
                                   ! dimension of the variable)
  INTEGER, INTENT(IN) :: count(MAX_DIM_VAR)  ! The number of points to read
                                   ! in each dimension of the variable


! Work variables
  TYPE(var_ncdf) :: var  ! Structure containing information about the
                         ! variable

  INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                         ! has, including record dimension

  INTEGER, ALLOCATABLE :: local_start(:)
  INTEGER, ALLOCATABLE :: local_count(:)

  INTEGER :: error      ! Error indicator


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
  var = file%vars(var_id)

  ndims_tot = var%ndims
  IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
  ALLOCATE(local_start(ndims_tot))
  local_start(1:var%ndims) = start(1:var%ndims)
  IF ( var%is_record ) local_start(ndims_tot) = file%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given count array, plus a value so that
! we read from only the correct record, if required
  ALLOCATE(local_count(ndims_tot))
  local_count(1:var%ndims) = count(1:var%ndims)
  IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
  error = nf90_get_var(file%id, var%id, values, local_start, local_count)
  IF ( error /= NF90_NOERR )                                                  &
    CALL log_fatal_ncdf("file_ncdf_read_var_2d",                              &
                        "Error reading variable", error)

  DEALLOCATE(local_start)
  DEALLOCATE(local_count)

  RETURN

END SUBROUTINE file_ncdf_read_var_2d


SUBROUTINE file_ncdf_read_var_3d(file, var_id, values, start, count)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE netcdf, ONLY : nf90_get_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 3d array from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, DIMENSION(:,:,:), INTENT(OUT) :: values  ! The values read from file
  INTEGER, INTENT(IN) :: start(MAX_DIM_VAR)  ! The point to start reading from
                                   ! (one value for each NON-RECORD
                                   ! dimension of the variable)
  INTEGER, INTENT(IN) :: count(MAX_DIM_VAR)  ! The number of points to read
                                   ! in each dimension of the variable


! Work variables
  TYPE(var_ncdf) :: var  ! Structure containing information about the
                         ! variable

  INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                         ! has, including record dimension

  INTEGER, ALLOCATABLE :: local_start(:)
  INTEGER, ALLOCATABLE :: local_count(:)

  INTEGER :: error      ! Error indicator


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
  var = file%vars(var_id)

  ndims_tot = var%ndims
  IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
  ALLOCATE(local_start(ndims_tot))
  local_start(1:var%ndims) = start(1:var%ndims)
  IF ( var%is_record ) local_start(ndims_tot) = file%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given count array, plus a value so that
! we read from only the correct record, if required
  ALLOCATE(local_count(ndims_tot))
  local_count(1:var%ndims) = count(1:var%ndims)
  IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
  error = nf90_get_var(file%id, var%id, values, local_start, local_count)
  IF ( error /= NF90_NOERR )                                                  &
    CALL log_fatal_ncdf("file_ncdf_read_var_3d",                              &
                        "Error reading variable", error)

  DEALLOCATE(local_start)
  DEALLOCATE(local_count)

  RETURN

END SUBROUTINE file_ncdf_read_var_3d


SUBROUTINE file_ncdf_read_var_4d(file, var_id, values, start, count)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE netcdf, ONLY : nf90_get_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 4d array from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: values  ! The values read from file
  INTEGER, INTENT(IN) :: start(MAX_DIM_VAR)  ! The point to start reading from
                                   ! (one value for each NON-RECORD
                                   ! dimension of the variable)
  INTEGER, INTENT(IN) :: count(MAX_DIM_VAR)  ! The number of points to read
                                   ! in each dimension of the variable


! Work variables
  TYPE(var_ncdf) :: var  ! Structure containing information about the
                         ! variable

  INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                         ! has, including record dimension

  INTEGER, ALLOCATABLE :: local_start(:)
  INTEGER, ALLOCATABLE :: local_count(:)

  INTEGER :: error      ! Error indicator


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
  var = file%vars(var_id)

  ndims_tot = var%ndims
  IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
  ALLOCATE(local_start(ndims_tot))
  local_start(1:var%ndims) = start(1:var%ndims)
  IF ( var%is_record ) local_start(ndims_tot) = file%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given count array, plus a value so that
! we read from only the correct record, if required
  ALLOCATE(local_count(ndims_tot))
  local_count(1:var%ndims) = count(1:var%ndims)
  IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
  error = nf90_get_var(file%id, var%id, values, local_start, local_count)
  IF ( error /= NF90_NOERR )                                                  &
    CALL log_fatal_ncdf("file_ncdf_read_var_4d",                              &
                        "Error reading variable", error)

  DEALLOCATE(local_start)
  DEALLOCATE(local_count)

  RETURN

END SUBROUTINE file_ncdf_read_var_4d


SUBROUTINE file_ncdf_read_var_5d(file, var_id, values, start, count)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE netcdf, ONLY : nf90_get_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 5d array from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, DIMENSION(:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file
  INTEGER, INTENT(IN) :: start(MAX_DIM_VAR)  ! The point to start reading from
                                   ! (one value for each NON-RECORD
                                   ! dimension of the variable)
  INTEGER, INTENT(IN) :: count(MAX_DIM_VAR)  ! The number of points to read
                                   ! in each dimension of the variable


! Work variables
  TYPE(var_ncdf) :: var  ! Structure containing information about the
                         ! variable

  INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                         ! has, including record dimension

  INTEGER, ALLOCATABLE :: local_start(:)
  INTEGER, ALLOCATABLE :: local_count(:)

  INTEGER :: error      ! Error indicator


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
  var = file%vars(var_id)

  ndims_tot = var%ndims
  IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
  ALLOCATE(local_start(ndims_tot))
  local_start(1:var%ndims) = start(1:var%ndims)
  IF ( var%is_record ) local_start(ndims_tot) = file%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given count array, plus a value so that
! we read from only the correct record, if required
  ALLOCATE(local_count(ndims_tot))
  local_count(1:var%ndims) = count(1:var%ndims)
  IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
  error = nf90_get_var(file%id, var%id, values, local_start, local_count)
  IF ( error /= NF90_NOERR )                                                  &
    CALL log_fatal_ncdf("file_ncdf_read_var_5d",                              &
                        "Error reading variable", error)

  DEALLOCATE(local_start)
  DEALLOCATE(local_count)

  RETURN

END SUBROUTINE file_ncdf_read_var_5d


SUBROUTINE file_ncdf_read_var_6d(file, var_id, values, start, count)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE netcdf, ONLY : nf90_get_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 6d array from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file
  INTEGER, INTENT(IN) :: start(MAX_DIM_VAR)  ! The point to start reading from
                                   ! (one value for each NON-RECORD
                                   ! dimension of the variable)
  INTEGER, INTENT(IN) :: count(MAX_DIM_VAR)  ! The number of points to read
                                   ! in each dimension of the variable


! Work variables
  TYPE(var_ncdf) :: var  ! Structure containing information about the
                         ! variable

  INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                         ! has, including record dimension

  INTEGER, ALLOCATABLE :: local_start(:)
  INTEGER, ALLOCATABLE :: local_count(:)

  INTEGER :: error      ! Error indicator


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
  var = file%vars(var_id)

  ndims_tot = var%ndims
  IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
  ALLOCATE(local_start(ndims_tot))
  local_start(1:var%ndims) = start(1:var%ndims)
  IF ( var%is_record ) local_start(ndims_tot) = file%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given count array, plus a value so that
! we read from only the correct record, if required
  ALLOCATE(local_count(ndims_tot))
  local_count(1:var%ndims) = count(1:var%ndims)
  IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
  error = nf90_get_var(file%id, var%id, values, local_start, local_count)
  IF ( error /= NF90_NOERR )                                                  &
    CALL log_fatal_ncdf("file_ncdf_read_var_6d",                              &
                        "Error reading variable", error)

  DEALLOCATE(local_start)
  DEALLOCATE(local_count)

  RETURN

END SUBROUTINE file_ncdf_read_var_6d


SUBROUTINE file_ncdf_read_var_7d(file, var_id, values, start, count)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE netcdf, ONLY : nf90_get_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 7d array from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, DIMENSION(:,:,:,:,:,:,:), INTENT(OUT) :: values  ! The values read from file
  INTEGER, INTENT(IN) :: start(MAX_DIM_VAR)  ! The point to start reading from
                                   ! (one value for each NON-RECORD
                                   ! dimension of the variable)
  INTEGER, INTENT(IN) :: count(MAX_DIM_VAR)  ! The number of points to read
                                   ! in each dimension of the variable


! Work variables
  TYPE(var_ncdf) :: var  ! Structure containing information about the
                         ! variable

  INTEGER :: ndims_tot  ! The total number of dimensions that the variable
                         ! has, including record dimension

  INTEGER, ALLOCATABLE :: local_start(:)
  INTEGER, ALLOCATABLE :: local_count(:)

  INTEGER :: error      ! Error indicator


!-----------------------------------------------------------------------------

! Copy values into a local variable for convenience
  var = file%vars(var_id)

  ndims_tot = var%ndims
  IF ( var%is_record ) ndims_tot = ndims_tot + 1

! Set up the local_start array - this consists of one value for each
! dimension of the variable from the given start array, plus a value so that
! we read from the correct record, if required
  ALLOCATE(local_start(ndims_tot))
  local_start(1:var%ndims) = start(1:var%ndims)
  IF ( var%is_record ) local_start(ndims_tot) = file%current_record

! Set up the local_count array - this consists of one value for each
! dimension of the variable from the given count array, plus a value so that
! we read from only the correct record, if required
  ALLOCATE(local_count(ndims_tot))
  local_count(1:var%ndims) = count(1:var%ndims)
  IF ( var%is_record ) local_count(ndims_tot) = 1

! Try to read the variable
  error = nf90_get_var(file%id, var%id, values, local_start, local_count)
  IF ( error /= NF90_NOERR )                                                  &
    CALL log_fatal_ncdf("file_ncdf_read_var_7d",                              &
                        "Error reading variable", error)

  DEALLOCATE(local_start)
  DEALLOCATE(local_count)

  RETURN

END SUBROUTINE file_ncdf_read_var_7d
#endif
