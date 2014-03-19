#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ts_def_time_dim(file, dim_name) RESULT(dim_id)

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
  TYPE(file_ts), INTENT(INOUT) :: file
                                  ! The file to define the dimension on
  CHARACTER(len=*), INTENT(IN) :: dim_name
                                  ! The name of the dimension

! Return type
  INTEGER :: dim_id               ! The dimension id


!-----------------------------------------------------------------------------

! If we are not in define mode, error out
  IF ( .NOT. file%define_mode )                                               &
    CALL log_fatal("file_ts_def_time_dim",                                    &
                   "Cannot define time dimension - file is not in define mode")

! A file can only have one time dimension
  IF ( file%has_time_dim )                                                    &
    CALL log_fatal("file_ts_def_time_dim",                                    &
                   "Time dimension has already been defined")

!-----------------------------------------------------------------------------
! Store the dimension attributes so that they can be used later to define
! the dimension on actual file(s)
!-----------------------------------------------------------------------------
  file%has_time_dim = .TRUE.

  file%time_dim%name   = dim_name
  file%time_dim%length = -1  ! The time dimension does not have a length

! Return -1 as the dimension id, since the return value should never be used
! and -1 can never be a valid dimension id for a non-time dimension
  dim_id = -1

  RETURN

END FUNCTION file_ts_def_time_dim
#endif
