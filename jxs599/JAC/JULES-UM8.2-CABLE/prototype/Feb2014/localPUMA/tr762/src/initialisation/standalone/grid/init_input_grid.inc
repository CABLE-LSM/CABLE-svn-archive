#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_input_grid()

  USE io_constants, ONLY : NAMELIST_UNIT, MAX_SDF_NAME_LEN

  USE grid_utils_mod, ONLY : grid_create

  USE string_utils_mod, ONLY : to_string

  USE input_mod, ONLY : input_grid => grid

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
  LOGICAL :: grid_is_1d  ! T - grid is a vector
                         ! F - grid is 2D (e.g. lat/lon, but not necessarily)

! Used if grid_is_1d = T - the name and size of the single grid dimension
  CHARACTER(len=MAX_SDF_NAME_LEN) :: grid_dim_name
  INTEGER :: npoints

! Used if grid_is_1d = F - the names and sizes of the x and y dimensions
  CHARACTER(len=MAX_SDF_NAME_LEN) :: x_dim_name, y_dim_name
  INTEGER :: nx, ny

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
  grid_is_1d = .FALSE.

  grid_dim_name = "land"
  npoints       = 0

  x_dim_name = "x"
  y_dim_name = "y"
  nx         = 0
  ny         = 0

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
! Build the input grid object from the namelist values
!-----------------------------------------------------------------------------
  input_grid = grid_create(                                                   &
    grid_is_1d, grid_dim_name, npoints, x_dim_name, nx, y_dim_name, ny        &
  )

! Check that a grid of non-zero size has been specified
  IF ( input_grid%nx * input_grid%ny <= 0 )                                   &
    CALL log_fatal("init_input_grid",                                         &
                   "Invalid dimensions specified for grid - " //              &
                   "nx =  " // TRIM(to_string(input_grid%nx)) // "; " //      &
                   "ny = " // TRIM(to_string(input_grid%ny)))

  CALL log_info("init_input_grid",                                            &
                "Size of input grid - " // TRIM(to_string(input_grid%nx)) //  &
                " x " // TRIM(to_string(input_grid%ny)))

  RETURN

END SUBROUTINE init_input_grid
#endif
