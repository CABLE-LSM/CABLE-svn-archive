#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_def_grid(file, grid)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines the grid used by variables in the file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ts), INTENT(INOUT) :: file
      ! The file to define the grid on

  TYPE(grid_info), INTENT(IN) :: grid  ! The grid to define


!-----------------------------------------------------------------------------


! If we are not in define mode, error out
  IF ( .NOT. file%define_mode )                                               &
    CALL log_fatal("file_ts_def_grid",                                        &
                   "Cannot define grid - file is not in define mode")

  file%grid = grid


  RETURN

END SUBROUTINE file_ts_def_grid
#endif
