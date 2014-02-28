#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/gridded/file_gridded_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

MODULE file_gridded_mod

  USE io_constants, ONLY : MAX_DIM_FILE, MAX_VAR_FILE

  USE file_mod, ONLY : file_handle

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE


! Type definitions
  TYPE file_gridded

    TYPE(file_handle) :: fh  ! Handle to the underlying file

!-----------------------------------------------------------------------------
! Information about the grid to use
!-----------------------------------------------------------------------------
    LOGICAL :: grid_is_1d = .TRUE.  ! T - file is using a 1d grid
                                    ! F - file is using a 2d grid

    INTEGER :: grid_dim_id = -1  ! Dimension id of the single grid dimension
                                 ! if using a 1d grid
    INTEGER :: grid_x_dim_id = -1  ! Dimension ids of the x and y dimensions
    INTEGER :: grid_y_dim_id = -1  ! if using a 2d grid

    INTEGER :: grid_nx = 0  ! The sizes of the x and y dimensions if using a
    INTEGER :: grid_ny = 0  ! 2d grid

!-----------------------------------------------------------------------------
! Information about the 'levels' dimensions that have been defined
! (i.e. dimensions that define the number of vertical levels a variable has)
! Each variable can have at most one 'levels' dimension
!-----------------------------------------------------------------------------
    INTEGER :: ndims = 0  ! The number of 'levels' dimensions defined
    INTEGER :: dim_ids(MAX_DIM_FILE)  ! The ids of the 'levels' dimensions
                                      ! in the file
    INTEGER :: dim_sizes(MAX_DIM_FILE)  ! The sizes of the 'levels' dimensions

!-----------------------------------------------------------------------------
! Information about the number of levels that variable in the file have
!-----------------------------------------------------------------------------
    INTEGER :: nvars = 0  ! The number of variables in the file
    INTEGER :: var_ids(MAX_VAR_FILE)  ! The ids of the variables in the file
    LOGICAL :: var_has_levs(MAX_VAR_FILE)  ! T - variable has vertical levels
                                           ! F - variable has no vertical levels
    INTEGER :: var_nlevs(MAX_VAR_FILE)  ! The number of vertical levels the
                                        ! variable has

  END TYPE file_gridded


! Overloads for file_gridded_def_attr
  INTERFACE file_gridded_def_attr
    MODULE PROCEDURE file_gridded_def_attr_real, file_gridded_def_attr_int,   &
                     file_gridded_def_attr_char
  END INTERFACE file_gridded_def_attr


CONTAINS

! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "file_gridded_open.inc"
#include "file_gridded_def_grid.inc"
#include "file_gridded_def_dim.inc"
#include "file_gridded_def_record_dim.inc"
#include "file_gridded_def_var.inc"
#include "file_gridded_def_attr.inc"
#include "file_gridded_enddef.inc"
#include "file_gridded_seek.inc"
#include "file_gridded_advance.inc"
#include "file_gridded_read_var.inc"
#include "file_gridded_write_var.inc"
#include "file_gridded_close.inc"

END MODULE file_gridded_mod
#endif
