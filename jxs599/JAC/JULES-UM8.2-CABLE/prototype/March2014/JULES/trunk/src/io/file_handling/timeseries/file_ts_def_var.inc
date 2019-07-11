#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ts_def_var(file, var_name, levels_dims, use_time) RESULT(var_id)

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
  TYPE(file_ts), INTENT(INOUT) :: file
                                  ! The file to define the variable in
  CHARACTER(len=*), INTENT(IN) :: var_name
                                    ! The name of the variable
  INTEGER, INTENT(IN), OPTIONAL :: levels_dims(:)
                                    ! The ids of the dimensions to use for the
                                    ! vertical levels of the variable, if
                                    ! required
                                    ! If not given or an array of 0 length is
                                    ! given, no vertical levels are used
  LOGICAL, INTENT(IN) :: use_time   ! Indicates whether the variable uses the
                                    ! time dimension

! Return type
  INTEGER :: var_id               ! The variable id


! Work variables


!-----------------------------------------------------------------------------

! If we are not in define mode, error out
  IF ( .NOT. file%define_mode )                                               &
    CALL log_fatal("file_ts_def_var",                                         &
                   "Cannot define variable - file is not in define mode")

! If adding another variable will cause us to have too many variables,
! error out
  IF ( file%nvars >= MAX_VAR_FILE )                                           &
    CALL log_fatal("file_ts_def_var",                                         &
                   "Too many variables in file - try increasing MAX_VAR_FILE")


!-----------------------------------------------------------------------------
! Store the variable attributes so that they can be used later to define
! variables on actual file(s)
!-----------------------------------------------------------------------------

  file%nvars = file%nvars + 1

! The returned variable id is just the index in the vars array on the file object
  var_id = file%nvars

  file%vars(var_id)%name       = var_name
  file%vars(var_id)%use_time   = use_time

! This exploits the fact that providing a levels_dims array of size 0 to
! file_gridded_def_var is the same as not providing it
  IF ( PRESENT(levels_dims) ) THEN
    ALLOCATE(file%vars(var_id)%levels_dims(SIZE(levels_dims)))
    file%vars(var_id)%levels_dims(:) = levels_dims(:)
  ELSE
    ALLOCATE(file%vars(var_id)%levels_dims(0))
  END IF

  RETURN

END FUNCTION file_ts_def_var
#endif
