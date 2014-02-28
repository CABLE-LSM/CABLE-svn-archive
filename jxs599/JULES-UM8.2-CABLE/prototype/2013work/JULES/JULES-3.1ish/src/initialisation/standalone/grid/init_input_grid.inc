#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/grid/init_input_grid.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_input_grid()

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE input_mod, ONLY : grid_is_1d, grid_dim_name, x_dim_name, y_dim_name,    &
                        nx, ny

  USE model_interface_mod, ONLY : pft_dim_name, nvg_dim_name, type_dim_name,  &
                                  tile_dim_name, snow_dim_name, soil_dim_name,&
                                  scpool_dim_name

  USE time_varying_input_mod, ONLY : time_dim_name

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about the grid for input data and checks for
!   consistency
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the jules_input_grid namelist - this combines variables
! from input_mod with local variables
!-----------------------------------------------------------------------------
  INTEGER :: npoints
  NAMELIST /jules_input_grid/ grid_is_1d,                                     &
! Required options for input on a 1d grid
                              grid_dim_name, npoints,                         &
! Required options for input on a 2d grid
                              x_dim_name, y_dim_name, nx, ny,                 &
! Levels dimension names
                              pft_dim_name, nvg_dim_name, type_dim_name,      &
                              tile_dim_name, snow_dim_name, soil_dim_name,    &
                              scpool_dim_name,                                &
! Name of the time dimension for time-varying files
                              time_dim_name


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  npoints = 0

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_input_grid", "Reading JULES_INPUT_GRID namelist...")

! First, we read the namelist
  READ(NAMELIST_UNIT, nml=jules_input_grid, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_input_grid",                                         &
                   "Error reading namelist JULES_INPUT_GRID " //              &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Set values derived from namelist and verify for consistency
!-----------------------------------------------------------------------------
! Reset nx and ny if we are using a 1d grid
  IF ( grid_is_1d ) THEN
    nx = npoints
    ny = 1
  END IF

! Check that a grid of non-zero size has been specified
  IF ( nx * ny <= 0 )                                                         &
    CALL log_fatal("init_input_grid",                                         &
                   "Invalid dimensions specified for grid - nx =  " //        &
                   TRIM(to_string(nx)) // "; ny = " // TRIM(to_string(ny)))

  CALL log_info("init_input_grid",                                            &
                "Size of input grid - " // TRIM(to_string(nx)) // " x " //    &
                TRIM(to_string(ny)))

  RETURN

END SUBROUTINE init_input_grid
#endif
