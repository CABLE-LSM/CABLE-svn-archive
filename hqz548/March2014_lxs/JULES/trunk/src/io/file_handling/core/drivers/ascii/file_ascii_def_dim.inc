#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ascii_def_dim(file, dim_name, dim_len) RESULT(dim_id)

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
  TYPE(file_ascii), INTENT(INOUT) :: file
                                  ! The file to define the dimension on
  CHARACTER(len=*), INTENT(IN) :: dim_name
                                  ! The name of the dimension
  INTEGER, INTENT(IN) :: dim_len  ! The length of the dimension

! Return type
  INTEGER :: dim_id               ! The dimension id


!-----------------------------------------------------------------------------

  IF ( .NOT. file%define_mode )                                               &
    CALL log_fatal("file_ascii_def_dim",                                      &
                   "Cannot define dimension - file is not in define mode")

! If adding another dimension will cause us to have too many dimensions,
! error out
  IF ( file%ndims >= MAX_DIM_FILE )                                           &
    CALL log_fatal("file_ascii_def_dim",                                      &
                   "Too many dimensions in file - try increasing MAX_DIM_FILE")

! Check the dimension has a sensible size (i.e. > 0)
  IF ( dim_len < 1 )                                                          &
    CALL log_fatal("file_ascii_def_dim",                                      &
                   "Dimension size should be >= 1 - size of " //              &
                   TRIM(to_string(dim_len)) // " given for dimension " //     &
                   TRIM(dim_name))

!-----------------------------------------------------------------------------
! Store the dimension info to use later
!-----------------------------------------------------------------------------
  file%ndims = file%ndims + 1

! The returned dimension id is just the index in the dims array on the file object
  dim_id = file%ndims

  file%dim_names(dim_id) = dim_name
  file%dim_sizes(dim_id) = dim_len

  RETURN

END FUNCTION file_ascii_def_dim
#endif
