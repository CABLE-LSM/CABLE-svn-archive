#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ascii_def_var(file, var_name, dims, is_record) RESULT(var_id)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines a variable in the given file, returning the variable id
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ascii), INTENT(INOUT) :: file
                                    ! The file to define the variable in
  CHARACTER(len=*), INTENT(IN) :: var_name
                                    ! The name of the variable
  INTEGER, INTENT(IN) :: dims(:)    ! The ids of the NON-RECORD dimensions of
                                    ! the variable
  LOGICAL, INTENT(IN) :: is_record  ! Indicates whether the variable uses the
                                    ! record dimension

! Return type
  INTEGER :: var_id                 ! The variable id - this is the index
                                    ! in the vars array of the file_ascii object
                                    ! of the variable


! Work variables
  INTEGER :: dim_sizes(SIZE(dims))  ! The size of each specified dimension

  INTEGER :: i  ! Loop variable


!-----------------------------------------------------------------------------


! If we are not in define mode, error out
  IF ( .NOT. file%define_mode )                                               &
    CALL log_fatal("file_ascii_def_var",                                      &
                   "Cannot define variable - file is not in define mode")

! If adding another variable will cause us to have too many variables,
! error out
  IF ( file%nvars >= MAX_VAR_FILE )                                           &
    CALL log_fatal("file_ascii_def_var",                                      &
                   "Too many variables in file - try increasing MAX_VAR_FILE")

! The record dimension needs to be defined before any variables can use it
  IF ( is_record .AND. .NOT. file%has_record_dim )                            &
    CALL log_fatal("file_ascii_def_var",                                      &
                   "Record dimension must be defined before defining a " //   &
                   "variable as using it")

! If the record dimension is defined, all variables must use it
  IF ( file%has_record_dim .AND. .NOT. is_record )                            &
    CALL log_fatal("file_ascii_def_var",                                      &
                   "Record dimension has been defined, so all variables " //  &
                   "must use it")

! Check that all the given dimensions are valid
  IF ( ANY(dims < 1) )                                                        &
    CALL log_fatal("file_ascii_def_var",                                      &
                   "Negative dimension id given - note that use of the " //   &
                   "record dimension should be specified using 'is_record'")

! Get the dimension sizes from the given ids
  DO i = 1,SIZE(dims)
    dim_sizes(i) = file%dim_sizes( dims(i) )
  END DO


!-----------------------------------------------------------------------------
! Save information about the variable for later use
!-----------------------------------------------------------------------------
  file%nvars = file%nvars + 1

! The return value is the index of the variable in the vars array on the
! file_ascii object
  var_id = file%nvars

  file%var_names(var_id) = var_name
  file%var_sizes(var_id) = PRODUCT(dim_sizes)

  file%var_ndims(var_id) = SIZE(dims)
  file%var_dims(var_id,1:SIZE(dims)) = dims(:)

  RETURN

END FUNCTION file_ascii_def_var
#endif
