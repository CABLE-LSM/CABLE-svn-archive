#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/gridded/file_gridded_def_var.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

FUNCTION file_gridded_def_var(file, var_name, levels_dim, is_record)          &
                                                                RESULT(var_id)

  USE file_mod, ONLY : file_def_var

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
  TYPE(file_gridded), INTENT(INOUT) :: file
                                  ! The file to define the variable in
  CHARACTER(len=*), INTENT(IN) :: var_name
                                    ! The name of the variable
  INTEGER, INTENT(IN), OPTIONAL :: levels_dim
                                    ! The id of the dimension to use for the
                                    ! vertical levels of the variable, if
                                    ! required
                                    ! If a value of levels_dim < 1 is given,
                                    ! that is the same as not providing it
  LOGICAL, INTENT(IN) :: is_record  ! Indicates whether the variable uses the
                                    ! record dimension

! Return type
  INTEGER :: var_id               ! The variable id


! Work variables
  INTEGER :: levels_dim_local  ! Value of levels dim to use
  INTEGER :: ndims       ! The number of dimensions that the variable has
  INTEGER :: dim_ids(3)  ! The ids in the file of the dimensions of the
                         ! variable
                         ! The variable can have at most 3 non-record
                         ! dimensions (1 or 2 for grid + 0 or 1 for levels)


!-----------------------------------------------------------------------------


! Check if we already have the maximum number of variables
  IF ( file%nvars >= MAX_VAR_FILE )                                           &
    CALL log_fatal("file_gridded_def_var",                                    &
                   "Too many variables in file - try increasing MAX_VAR_FILE")

!-----------------------------------------------------------------------------
! Construct the array of dimension ids for the variable
!-----------------------------------------------------------------------------
  levels_dim_local = -1
  IF ( PRESENT(levels_dim) ) levels_dim_local = levels_dim

  dim_ids(:) = -1
! Grid dimensions come first
  IF ( file%grid_is_1d ) THEN
    ndims = 1
    dim_ids(1) = file%grid_dim_id
  ELSE
    ndims = 2
    dim_ids(1:2) = (/ file%grid_x_dim_id, file%grid_y_dim_id /)
  END IF

  IF ( levels_dim_local > 0 ) THEN
    ndims = ndims + 1
! levels_dim is an index to file%dim_ids as returned by file_gridded_def_dim
    dim_ids(ndims) = file%dim_ids(levels_dim_local)
  END IF

!-----------------------------------------------------------------------------
! Define the variable in the underlying file and store details for future use
!-----------------------------------------------------------------------------
  file%nvars = file%nvars + 1

! The returned variable id is just the index in the var_ids array of the file
  var_id = file%nvars

  file%var_ids(var_id) = file_def_var(                                        &
    file%fh, var_name, dim_ids(1:ndims), is_record                            &
  )
  file%var_has_levs(var_id) = .FALSE.
  IF ( levels_dim_local > 0 ) THEN
    file%var_has_levs(var_id) = .TRUE.
    file%var_nlevs(var_id)    = file%dim_sizes(levels_dim_local)
  END IF

  RETURN

END FUNCTION file_gridded_def_var
#endif
