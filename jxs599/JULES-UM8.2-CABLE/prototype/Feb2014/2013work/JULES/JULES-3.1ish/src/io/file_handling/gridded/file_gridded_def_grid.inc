#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/gridded/file_gridded_def_grid.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE file_gridded_def_grid(file, is_1d, grid_name, grid_npoints,        &
                                 grid_x_name, grid_nx, grid_y_name, grid_ny)

  USE file_mod, ONLY : file_def_dim

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
  TYPE(file_gridded), INTENT(INOUT) :: file
      ! The file to define the grid on

  LOGICAL, INTENT(IN) :: is_1d
      ! T - define a 1d grid (i.e. a vector)
      ! F - define a 2d grid

  CHARACTER(len=*), INTENT(IN), OPTIONAL :: grid_name
      ! ONLY USED IF is_1d=T
      ! The name of the single grid dimension
  INTEGER, INTENT(IN), OPTIONAL :: grid_npoints
      ! ONLY USED IF is_1d=T
      ! The size of the single grid dimension

  CHARACTER(len=*), INTENT(IN), OPTIONAL :: grid_x_name, grid_y_name
      ! ONLY USED IF is_1d=F
      ! The names of the x and y dimensions of the grid respectively
  INTEGER, INTENT(IN), OPTIONAL :: grid_nx, grid_ny
      ! ONLY USED IF is_1d=F
      ! The sizes of the x and y dimensions respectively


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that we have the right combination of arguments to do something
!-----------------------------------------------------------------------------
  IF ( is_1d ) THEN
! For a 1d grid, we need grid_name and grid_npoints
    IF ( .NOT. PRESENT(grid_name) .OR. .NOT. PRESENT(grid_npoints) )          &
      CALL log_fatal("file_gridded_def_grid",                                 &
                     "If 1d grid is specified, both grid_name and " //        &
                     "grid_npoints must be given")
  ELSE
! For a 2d grid, we need information about the x and y dimensions
    IF ( .NOT. PRESENT(grid_x_name) .OR. .NOT. PRESENT(grid_nx) .OR.          &
         .NOT. PRESENT(grid_y_name) .OR. .NOT. PRESENT(grid_ny) )             &
      CALL log_fatal("file_gridded_def_grid",                                 &
                     "If 2d grid is specified, grid_x_name, grid_nx, " //     &
                     "grid_y_name and grid_ny must all be given")
  END IF

!-----------------------------------------------------------------------------
! Define the dimensions for the grid on the underlying file and store the
! dimension sizes and ids for later use
!-----------------------------------------------------------------------------
  file%grid_is_1d = is_1d

  IF ( is_1d ) THEN
    file%grid_nx      = grid_npoints
    file%grid_ny      = 1
    file%grid_dim_id  = file_def_dim(file%fh, grid_name, grid_npoints)
  ELSE
    file%grid_nx       = grid_nx
    file%grid_ny       = grid_ny
    file%grid_x_dim_id = file_def_dim(file%fh, grid_x_name, grid_nx)
    file%grid_y_dim_id = file_def_dim(file%fh, grid_y_name, grid_ny)
  END IF

  RETURN

END SUBROUTINE file_gridded_def_grid
#endif
