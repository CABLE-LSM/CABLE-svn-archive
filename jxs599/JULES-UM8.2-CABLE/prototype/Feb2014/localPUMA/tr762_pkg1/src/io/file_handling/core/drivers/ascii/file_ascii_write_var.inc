#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_write_var_scalar(file, var_id, value, start)

  USE io_constants, ONLY : MAX_DIM_VAR

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a scalar value to a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ascii), INTENT(INOUT) :: file  ! The file to write to
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
  REAL, INTENT(IN) :: value  ! The value to write
  INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                   ! (one value for each dimension
                                   ! of the variable in the file)

! Work variables
  REAL :: local_values(1)  ! An array version of value
  INTEGER :: local_count(MAX_DIM_VAR)  ! A version of count, which will be
                                       ! all 1s

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result
  local_count(:) = 1
  local_values(1) = value

  CALL file_ascii_write_var_1d(file, var_id, local_values, start, local_count)

  RETURN

END SUBROUTINE file_ascii_write_var_scalar


SUBROUTINE file_ascii_write_var_1d(file, var_id, values, start, count_in)

  USE io_constants, ONLY : MODE_WRITE

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 1d array to a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ascii), INTENT(INOUT) :: file  ! The file to write to
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
  REAL, DIMENSION(:), INTENT(IN) :: values  ! The values to write
  INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                   ! (one value for each dimension
                                   ! of the variable in the file)
  INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                      ! in each dimension of the variable

! Work variables
  INTEGER :: ndims  ! The number of dimensions that the variable has in file

  INTEGER :: offset  ! The offset to use in the buffer


!-----------------------------------------------------------------------------

! We can't write data if the file is in define mode
  IF ( file%define_mode )                                                     &
    CALL log_fatal("file_ascii_write_var_1d",                                 &
                   "Cannot write data - file is still in define mode")

! We can't write data unless we are in write mode...
  IF ( file%mode /= MODE_WRITE )                                              &
    CALL log_fatal("file_ascii_write_var_1d",                                 &
                   "Can only write data if file is opened in write mode")

  ndims = file%var_ndims(var_id)

! ASCII files have the restriction that we can only write the whole variable
  IF ( ANY(start(1:ndims) /= 1) )                                             &
    CALL log_fatal("file_ascii_write_var_1d",                                 &
                   "start must be 1 for all dimensions - writing part of " // &
                   "a variable in an ASCII file is not supported")

  IF ( SIZE(values) /= PRODUCT(count_in(1:ndims)) .AND.                       &
       SIZE(values) /= file%var_sizes(var_id) )                               &
    CALL log_fatal("file_ascii_write_var_1d",                                 &
                   "values must have the same number of elements as the " //  &
                   "variable in file - writing part of a variable in an " //  &
                   "ASCII file is not supported")

! Now we know we can write to the buffer correctly
  offset = file%var_offsets(var_id)
  file%buffer(offset:(offset + file%var_sizes(var_id) - 1)) = values(:)

! Mark the buffer as dirty (i.e. has writes that have not been committed to file
  file%buffer_is_dirty = .TRUE.

  RETURN

END SUBROUTINE file_ascii_write_var_1d


SUBROUTINE file_ascii_write_var_2d(file, var_id, values, start, count_in)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 2d array to a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ascii), INTENT(INOUT) :: file  ! The file to write to
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
  REAL, DIMENSION(:,:), INTENT(IN) :: values  ! The values to write
  INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                   ! (one value for each dimension
                                   ! of the variable in the file)
  INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                      ! in each dimension of the variable


! Work variables
  REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

  ALLOCATE(local_values(PRODUCT(SHAPE(values))))

  local_values = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

  CALL file_ascii_write_var_1d(file, var_id, local_values, start, count_in)

  DEALLOCATE(local_values)

  RETURN

END SUBROUTINE file_ascii_write_var_2d


SUBROUTINE file_ascii_write_var_3d(file, var_id, values, start, count_in)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 3d array to a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ascii), INTENT(INOUT) :: file  ! The file to write to
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
  REAL, DIMENSION(:,:,:), INTENT(IN) :: values  ! The values to write
  INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                   ! (one value for each dimension
                                   ! of the variable in the file)
  INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                      ! in each dimension of the variable


! Work variables
  REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

  ALLOCATE(local_values(PRODUCT(SHAPE(values))))

  local_values = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

  CALL file_ascii_write_var_1d(file, var_id, local_values, start, count_in)

  DEALLOCATE(local_values)

  RETURN

END SUBROUTINE file_ascii_write_var_3d


SUBROUTINE file_ascii_write_var_4d(file, var_id, values, start, count_in)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 4d array to a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ascii), INTENT(INOUT) :: file  ! The file to write to
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
  REAL, DIMENSION(:,:,:,:), INTENT(IN) :: values  ! The values to write
  INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                   ! (one value for each dimension
                                   ! of the variable in the file)
  INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                      ! in each dimension of the variable


! Work variables
  REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

  ALLOCATE(local_values(PRODUCT(SHAPE(values))))

  local_values = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

  CALL file_ascii_write_var_1d(file, var_id, local_values, start, count_in)

  DEALLOCATE(local_values)

  RETURN

END SUBROUTINE file_ascii_write_var_4d


SUBROUTINE file_ascii_write_var_5d(file, var_id, values, start, count_in)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 5d array to a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ascii), INTENT(INOUT) :: file  ! The file to write to
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
  REAL, DIMENSION(:,:,:,:,:), INTENT(IN) :: values  ! The values to write
  INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                   ! (one value for each dimension
                                   ! of the variable in the file)
  INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                      ! in each dimension of the variable


! Work variables
  REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

  ALLOCATE(local_values(PRODUCT(SHAPE(values))))

  local_values = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

  CALL file_ascii_write_var_1d(file, var_id, local_values, start, count_in)

  DEALLOCATE(local_values)

  RETURN

END SUBROUTINE file_ascii_write_var_5d


SUBROUTINE file_ascii_write_var_6d(file, var_id, values, start, count_in)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 6d array to a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ascii), INTENT(INOUT) :: file  ! The file to write to
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
  REAL, DIMENSION(:,:,:,:,:,:), INTENT(IN) :: values  ! The values to write
  INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                   ! (one value for each dimension
                                   ! of the variable in the file)
  INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                      ! in each dimension of the variable


! Work variables
  REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

  ALLOCATE(local_values(PRODUCT(SHAPE(values))))

  local_values = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

  CALL file_ascii_write_var_1d(file, var_id, local_values, start, count_in)

  DEALLOCATE(local_values)

  RETURN

END SUBROUTINE file_ascii_write_var_6d


SUBROUTINE file_ascii_write_var_7d(file, var_id, values, start, count_in)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes a 7d array to a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ascii), INTENT(INOUT) :: file  ! The file to write to
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to write to
  REAL, DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: values  ! The values to write
  INTEGER, INTENT(IN) :: start(:)  ! The point to start writing at
                                   ! (one value for each dimension
                                   ! of the variable in the file)
  INTEGER, INTENT(IN) :: count_in(:)  ! The number of points to write
                                      ! in each dimension of the variable


! Work variables
  REAL, ALLOCATABLE :: local_values(:)  ! A 1d version of values

!-----------------------------------------------------------------------------

! To implement this routine, we just delegate to the 1d routine and reshape
! the result

  ALLOCATE(local_values(PRODUCT(SHAPE(values))))

  local_values = RESHAPE(values, (/ PRODUCT(SHAPE(values)) /))

  CALL file_ascii_write_var_1d(file, var_id, local_values, start, count_in)

  DEALLOCATE(local_values)

  RETURN

END SUBROUTINE file_ascii_write_var_7d
#endif
