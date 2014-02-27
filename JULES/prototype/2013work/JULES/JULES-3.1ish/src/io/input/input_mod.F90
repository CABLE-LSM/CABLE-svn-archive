#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/input/input_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

MODULE input_mod

  USE io_constants, ONLY : MAX_SDF_NAME_LEN

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   The module contains information about the grid and levels dimensions
!   to use for input, and routines for reading non-time-varying input
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Definition of the input grid
!-----------------------------------------------------------------------------
  LOGICAL :: grid_is_1d = .FALSE.
      ! T => input files contain data on a 1d grid (vector)
      ! F => input files contain data on a 2d grid

! USED ONLY IF is_1d=T
  CHARACTER(len=MAX_SDF_NAME_LEN) :: grid_dim_name = "land"
      ! The dimension name of the single grid dimension

! USED ONLY IF is_1d=F
  CHARACTER(len=MAX_SDF_NAME_LEN) :: x_dim_name = "x", y_dim_name = "y"

  INTEGER :: nx = 0, ny = 0
      ! The number of points in the x and y dimensions respectively
      ! If is_1d=T, nx is the number of points and ny=1
      ! Default values are set so that we get an error if no grid is
      ! specified


!-----------------------------------------------------------------------------
! Definition of the model grid relative to the input grid
!-----------------------------------------------------------------------------
  LOGICAL :: use_subgrid = .FALSE.
      ! T => the model grid is a subset of the input grid
      ! F => the model grid is the input grid

! Used if use_subgrid = T
  INTEGER :: subgrid_npoints = 0  ! The number of points in the subgrid
  INTEGER, ALLOCATABLE ::                                                     &
    points(:),                                                                &
        ! USED IF grid_is_1d=T - The indices of the points in the subgrid
    x_coords(:), y_coords(:)
        ! USED IF grid_is_1d=F - The indices of the selected points in the
        !                        x and y directions respectively


CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "fill_variables_from_file.inc"

END MODULE input_mod
#endif
