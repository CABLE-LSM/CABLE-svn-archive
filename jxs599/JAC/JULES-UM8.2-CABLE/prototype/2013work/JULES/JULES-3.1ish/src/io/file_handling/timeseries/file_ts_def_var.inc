#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/timeseries/file_ts_def_var.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

FUNCTION file_ts_def_var(file, var_name, levels_dim, use_time) RESULT(var_id)

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
  INTEGER, INTENT(IN), OPTIONAL :: levels_dim
                                    ! The id of the dimension to use for the
                                    ! vertical levels of the variable, if
                                    ! required
                                    ! If a value of levels_dim < 1 is given,
                                    ! that is the same as not providing it
  LOGICAL, INTENT(IN) :: use_time   ! Indicates whether the variable uses the
                                    ! time dimension

! Return type
  INTEGER :: var_id               ! The variable id


! Work variables
  INTEGER :: levels_dim_local  ! Value of levels dim to use


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
! This exploits the fact that providing a value of levels_dim < 1 to
! file_gridded_def_var is the same as not providing it
  levels_dim_local = -1
  IF ( PRESENT(levels_dim) ) levels_dim_local = levels_dim

  file%nvars = file%nvars + 1

! The returned variable id is just the index in the vars array on the file object
  var_id = file%nvars

  file%vars(var_id)%name       = var_name
  file%vars(var_id)%has_levels = levels_dim_local > 0
  file%vars(var_id)%levels_dim = levels_dim_local
  file%vars(var_id)%use_time   = use_time

  RETURN

END FUNCTION file_ts_def_var
#endif
