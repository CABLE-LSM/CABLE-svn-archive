#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_model_grid()

USE allocate_jules_arrays_mod, ONLY: allocate_jules_arrays
USE allocate_cable_progs_mod, ONLY: allocate_cable_progs

USE io_constants, ONLY: max_file_name_len, namelist_unit, points_file_unit

USE string_utils_mod, ONLY: to_string

USE grid_utils_mod, ONLY: subgrid_info, grid_create, subgrid_create,          &
                           subgrid_extract, subgrid_restrict

USE input_mod, ONLY: grid_in => grid,                                         &
! Variables to do with the extraction of subgrids if specified
! use_subgrid does not have an _in so that it can appear as use_subgrid in the
! namelist
                        use_subgrid, subgrid_in => subgrid

USE output_mod, ONLY: grid_out => grid,                                       &
! Variables to do with writing subgrids if multiple tasks are used
                         use_subgrid_out => use_subgrid,                      &
                         subgrid_out => subgrid

USE parallel_mod, ONLY: decompose_domain

USE switches, ONLY: l_co2_interactive

USE jules_soil_biogeochem_mod, ONLY: soil_model_1pool, soil_model_ecosse,     &
                                      soil_model_rothc, l_layeredC,           &
                                      soil_bgc_model

USE jules_surface_mod, ONLY: l_point_data, l_aggregate

USE jules_vegetation_mod, ONLY: l_triffid, l_crop

USE model_grid_mod, ONLY: model_grid, global_land_pts, global_land_mask,      &
                           latitude, longitude,                               &
                           latitude_of_land_pts, longitude_of_land_pts

USE coastal, ONLY: fland, flandg

USE atm_fields_bounds_mod, ONLY: atm_fields_bounds_init
USE theta_field_sizes, ONLY: t_i_length, t_j_length

USE ancil_info, ONLY: land_pts, land_mask, land_index, co2_dim_len,           &
                       co2_dim_row, dim_cs1, dim_cs2, land_pts_trif,          &
                       npft_trif, nsurft, dim_cslayer,                        &
                       row_length, rows, n_rows, nsoilt

USE jules_deposition_mod, ONLY: ndry_dep_species

USE jules_surface_types_mod, ONLY: npft, ncpft, nnvg, ntype

USE jules_snow_mod, ONLY: nsmax

USE jules_sea_seaice_mod, ONLY: nice, nice_use

USE jules_soil_mod, ONLY: sm_levels, ns_deep, l_tile_soil

USE model_interface_mod, ONLY: bl_level_dim_size, bedrock_dim_size,           &
                                pft_dim_size, cpft_dim_size, nvg_dim_size,    &
                                type_dim_size, tile_dim_size, snow_dim_size,  &
                                soil_dim_size, scpool_dim_size,               &
                                sclayer_dim_size, tracer_dim_size

USE nlsizes_namelist_mod, ONLY: bl_levels, jules_nlsizes

USE errormessagelength_mod, ONLY: errormessagelength

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about how the model grid is defined, sets it up and
!   checks for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
LOGICAL :: input_mask(grid_in%nx, grid_in%ny)
                                      ! Defined on the full input grid
                                      !   T - the point will be modelled
                                      !   F - the point will not be modelled
REAL(KIND=real_jlslsm) :: point_lat, point_lon
                              ! Latitude and longitude of points read from
                              ! points file

TYPE(subgrid_info) :: task_subgrid  ! The subgrid that the current MPI task
                                    ! will model, as a subgrid of the full
                                    ! model grid

REAL(KIND=real_jlslsm)                                                        &
    , ALLOCATABLE :: global_lat(:,:), global_lon(:,:), global_land_frac(:,:)
                             ! Latitude, longitude and land fraction on
                             ! the full model grid
                             ! These are not required outside of this
                             ! routine

INTEGER :: error, error_sum  ! Error indicators
CHARACTER(LEN=errormessagelength) :: iomessage

INTEGER :: i, j, l  ! Index variables

!-----------------------------------------------------------------------------
! Definition of the jules_model_grid namelist - this combines variables
! from input_mod with some local variables
!-----------------------------------------------------------------------------
LOGICAL :: land_only      ! T - only model land points
                          ! F - model all specified points


LOGICAL :: force_1d_grid  ! T - model grid is 1D
                          ! F - model grid is 1D unless input grid is 2D
                          !     and model grid is the whole input grid

LOGICAL :: latlon_region  ! T - subgrid is to be selected with latitude and
                          !     longitude bounds
                          ! F - subgrid is to be selected using a list
                          !     of latitudes and longitudes
REAL(KIND=real_jlslsm) :: lat_bounds(2)   ! USED IF latlon_region=T
REAL(KIND=real_jlslsm) :: lon_bounds(2)
                        ! Upper and lower bounds for latitude and longitude

INTEGER :: npoints = 0 ! The number of points to read from file
CHARACTER(LEN=max_file_name_len) :: points_file
                        ! The file to read latitudes and longitudes for
                        ! specified points

NAMELIST  / jules_model_grid/ land_only, use_subgrid, latlon_region,          &
                            lat_bounds, lon_bounds, npoints, points_file,     &
                            force_1d_grid


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
land_only = .TRUE.  ! Default is to model only land points from whatever
                    ! points are selected
force_1d_grid = .FALSE. ! Default is model grid is 1D unless input grid is 2D
                        ! and model grid is the whole input grid
points_file = ''    ! Empty file name.

!-----------------------------------------------------------------------------
! Read the init_model_grid namelist.
!-----------------------------------------------------------------------------
CALL log_info("init_model_grid", "Reading JULES_MODEL_GRID namelist...")

READ(namelist_unit, NML = jules_model_grid, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_model_grid",                                           &
                 "Error reading namelist JULES_MODEL_GRID " //                &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

!-----------------------------------------------------------------------------
! Read the jules_nlsizes namelist.
!-----------------------------------------------------------------------------
CALL log_info("init_model_grid", "Reading JULES_NLSIZES namelist...")

READ(namelist_unit, NML = jules_nlsizes, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal("init_model_grid",                                           &
                 "Error reading namelist JULES_NLSIZES " //                   &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

!-----------------------------------------------------------------------------
! Set up the model grid variables
!-----------------------------------------------------------------------------
CALL log_info("init_model_grid", "Setting up model grid variables...")

!*****************************************************************************
!*****************************************************************************
! Work out what points in the input grid will make up the model grid
!
! By the time we get here, we know that we have latitude, longitude and
! land fraction on the full input grid already read in
! Use these in combination with the values in the jules_model_grid namelist
! to work out what points to use
!*****************************************************************************
!*****************************************************************************
! To start with, assume we are using all the points
input_mask(:,:) = .TRUE.

IF ( use_subgrid ) THEN
  IF ( latlon_region ) THEN
    !-----------------------------------------------------------------------------
    ! The subgrid will be selected using latitude and longitude bounds
    !-----------------------------------------------------------------------------
    CALL log_info("init_model_grid",                                          &
                  "Subgrid will be selected using latitude and " //           &
                  "longitude bounds")
    CALL log_info("init_model_grid",                                          &
                  "Latitude range - " // TRIM(to_string(lat_bounds(1))) //    &
                  " to " // TRIM(to_string(lat_bounds(2))))
    CALL log_info("init_model_grid",                                          &
                  "Longitude range - " // TRIM(to_string(lon_bounds(1))) //   &
                  " to " // TRIM(to_string(lon_bounds(2))))

    input_mask = ( lat_bounds(1) <= latitude ) .AND.                          &
                 ( latitude <= lat_bounds(2) ) .AND.                          &
                 ( lon_bounds(1) <= longitude ) .AND.                         &
                 ( longitude <= lon_bounds(2) )
  ELSE
    !-----------------------------------------------------------------------------
    ! The subgrid will be selected using a list of latitudes and longitudes
    ! from the given file
    !-----------------------------------------------------------------------------
    CALL log_info("init_model_grid",                                          &
                  "Subgrid will be selected using a list of points " //       &
                  "from " // TRIM(points_file))

    !     Check that a file name was provided.
    IF ( LEN_TRIM(points_file) == 0 )                                         &
      CALL log_fatal("init_model_grid", "No points file provided" )

    ! No points have been specified yet, so set input_mask accordingly
    input_mask(:,:) = .FALSE.

    !-----------------------------------------------------------------------------
    ! Read the point latitudes and longitudes from file
    !-----------------------------------------------------------------------------
    OPEN(points_file_unit, FILE=points_file,                                  &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
    IF ( error /= 0 )                                                         &
    CALL log_fatal("init_model_grid",                                         &
                   "Error opening points file " //                            &
                   "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //       &
                   TRIM(iomessage) // ")")

    DO i = 1,npoints
      ! Try to read the next latitude/longitude pair
      READ(points_file_unit, *, IOSTAT = error, IOMSG = iomessage) point_lat, &
                                                               point_lon
      IF ( error /= 0 )                                                       &
        CALL log_fatal("init_model_grid",                                     &
                       "Error reading lat/lon pair from points file " //      &
                       "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //   &
                       TRIM(iomessage) // ")")

      ! Set input_mask to .TRUE. at that point
      WHERE ( ABS(latitude - point_lat) < EPSILON(1.0) .AND.                  &
              ABS(longitude - point_lon) < EPSILON(1.0) )
        input_mask = .TRUE.
      END WHERE
    END DO

    CLOSE(points_file_unit, IOSTAT = error, IOMSG = iomessage)
    IF ( error /= 0 )                                                         &
      CALL log_fatal("init_model_grid",                                       &
                     "Error closing points file " //                          &
                     "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //     &
                     TRIM(iomessage) // ")")
  END IF  ! latlon_region
END IF  ! use_subgrid

IF ( land_only ) THEN
  !-----------------------------------------------------------------------------
  ! If requested, select land points only from the points we have left
  !-----------------------------------------------------------------------------
  CALL log_info("init_model_grid",                                            &
                "From the points specified, only land points will be " //     &
                "modelled")

  input_mask = input_mask .AND. ( flandg > 0.0 )
END IF


!*****************************************************************************
!-----------------------------------------------------------------------------
! At this point, the .TRUE. points in input_mask will comprise the FULL
! model grid
! Issue an error if all points have been excluded
!-----------------------------------------------------------------------------
IF ( .NOT. ANY(input_mask) )                                                  &
  CALL log_fatal("init_model_grid",                                           &
                 "All points in input grid have been excluded from " //       &
                 "model grid")
!*****************************************************************************

!-----------------------------------------------------------------------------
! Build a subgrid from the input grid and mask
!-----------------------------------------------------------------------------
subgrid_in = subgrid_create(grid_in, input_mask, force_1d_grid)

! Construct the model grid to be the same size as the input subgrid
! It doesn't matter about dimension names as it will not be used for I/O
! It is always created as a 2D grid
model_grid = grid_create(                                                     &
  .FALSE., x_name = "", y_name = "", nx = subgrid_in%nx, ny = subgrid_in%ny   &
)

!-----------------------------------------------------------------------------
! Store the latitude, longitude and land fraction on the full model grid
!-----------------------------------------------------------------------------
ALLOCATE(global_lat(model_grid%nx,model_grid%ny))
global_lat(:,:) = subgrid_extract(subgrid_in, latitude)

ALLOCATE(global_lon(model_grid%nx,model_grid%ny))
global_lon(:,:) = subgrid_extract(subgrid_in, longitude)

ALLOCATE(global_land_frac(model_grid%nx,model_grid%ny))
global_land_frac(:,:) = subgrid_extract(subgrid_in, flandg)

! Work out how many land points are in the full model grid and set the land
! mask
ALLOCATE(global_land_mask(model_grid%nx,model_grid%ny))
global_land_mask(:,:) = ( global_land_frac(:,:) > EPSILON(1.0) )
global_land_pts = COUNT(global_land_mask)

!-----------------------------------------------------------------------------
! Decompose the model domain
!
! Returns the subgrid of the full model grid that the current task will model
!-----------------------------------------------------------------------------
task_subgrid = decompose_domain(model_grid)

!-----------------------------------------------------------------------------
! Set up the extraction of the input grid
!-----------------------------------------------------------------------------
! First, we combine the task subgrid with the subgrid generated from input_mask
! to get the input subgrid for this task
subgrid_in = subgrid_restrict(subgrid_in, task_subgrid)

! Do we actually have a subgrid, or does the subgrid cover the whole input grid?
use_subgrid = ( subgrid_in%nx /= grid_in%nx ) .OR.                            &
              ( subgrid_in%ny /= grid_in%ny )

IF ( use_subgrid ) THEN
  ! Regrid the latitude and longitude onto the model grid for this task
  DEALLOCATE(latitude)
  ALLOCATE(latitude(task_subgrid%nx,task_subgrid%ny))
  latitude(:,:) = subgrid_extract(task_subgrid, global_lat)

  DEALLOCATE(longitude)
  ALLOCATE(longitude(task_subgrid%nx,task_subgrid%ny))
  longitude(:,:) = subgrid_extract(task_subgrid, global_lon)

  DEALLOCATE(flandg)
  ALLOCATE(flandg(task_subgrid%nx,task_subgrid%ny))
  flandg(:,:) = subgrid_extract(task_subgrid, global_land_frac)
END IF

! Deallocate some arrays that are no longer required
DEALLOCATE(global_lat)
DEALLOCATE(global_lon)
DEALLOCATE(global_land_frac)

!-----------------------------------------------------------------------------
! Set up the output grid
!
! The full output grid is the same as the full model grid
! The subgrid that this task is responsible for outputting is the same as
! the subgrid it is responsible for modelling
!-----------------------------------------------------------------------------
grid_out    = model_grid
subgrid_out = task_subgrid

! Do we actually have a subgrid, or does the subgrid cover the whole output grid?
use_subgrid_out = ( subgrid_out%nx /= grid_out%nx ) .OR.                      &
                  ( subgrid_out%ny /= grid_out%ny )

!-----------------------------------------------------------------------------
! Set up the model grid for this task, i.e. the grid that the science code
! knows about. Only control code has to deal with the concept of a 'task
! model grid' and a 'full model grid'
!-----------------------------------------------------------------------------
t_i_length = task_subgrid%nx
t_j_length = task_subgrid%ny

! At the moment, JULES uses no halos, and the t, p, u and v grids are the same
CALL atm_fields_bounds_init(0, 0, 0, 0, t_i_length, t_j_length, t_j_length, 1,&
                            bl_levels_opt = bl_levels)
! Copy the values into row_length, rows and n_rows for legacy purposes
row_length = t_i_length
rows   = t_j_length
n_rows = t_j_length

CALL log_info("init_model_grid",                                              &
              "Size of model grid - " // TRIM(to_string(t_i_length)) //       &
              " x " // TRIM(to_string(t_j_length)))

! Calculate the number of land points
land_pts = COUNT(flandg > EPSILON(1.0))
CALL log_info("init_model_grid",                                              &
              "Selected grid contains " // TRIM(to_string(land_pts)) //       &
              " land points")

CALL log_info("init_model_grid",                                              &
              "bl_levels = " // TRIM(to_string(bl_levels)))

! Calculate number of surface types based on whether we are using the aggregate
! surface type
IF ( l_aggregate ) THEN
  nsurft = 1
ELSE
  nsurft = ntype
END IF

! Set up the number of soil tiles
IF ( l_tile_soil ) THEN
  nsoilt = nsurft
ELSE
  nsoilt = 1
  !This is also the default value, but reaffirming it here for clarity
END IF

! Set up the ice categories - there is only one
nice = 1
! nice_use MUST EQUAL NICE FOR STANDALONE JULES
nice_use = nice

! Set dimensions for TRIFFID arrays.
IF ( l_triffid .OR. l_crop ) THEN
  land_pts_trif = land_pts
  npft_trif     = npft
ELSE
  land_pts_trif = 1
  npft_trif     = 1
END IF

! Set dimensions for soil carbon model.
SELECT CASE ( soil_bgc_model )
CASE ( soil_model_1pool )
  dim_cs1     = 1  !  1 pool
  dim_cs2     = 1  !  1 as this model is not used with TRIFFID
  IF (l_layeredC) THEN
    dim_cslayer = sm_levels
  ELSE
    dim_cslayer = 1
  END IF
CASE ( soil_model_rothc )
  dim_cs1 = 4   !  for 4 soil pools
  dim_cs2 = land_pts
  IF (l_layeredC) THEN
    dim_cslayer = sm_levels
  ELSE
    dim_cslayer = 1
  END IF
CASE ( soil_model_ecosse )
  dim_cs1 = 4   !  4 soil pools: DPM, RPM, BIO, HUM
  dim_cs2 = land_pts
  ! dim_cslayer was set in jules_soil_biogeochem_mod.
END SELECT

IF ( l_co2_interactive ) THEN
  co2_dim_len = t_i_length
  co2_dim_row = t_j_length
ELSE
  co2_dim_len = 1
  co2_dim_row = 1
END IF

! Warn if using l_point_data for a non-single-point run
IF ( l_point_data .AND. t_i_length * t_j_length > 1 )                         &
  CALL log_warn("init_model_grid",                                            &
                "l_point_data is selected but there is more than one point")

!-----------------------------------------------------------------------------
! Set up the mapping of land points to the full model grid
!-----------------------------------------------------------------------------
ALLOCATE(land_mask(t_i_length,t_j_length), stat = error)
error_sum = error
ALLOCATE(land_index(land_pts), stat = error)
error_sum = error_sum + error
ALLOCATE(fland(land_pts), stat = error)
error_sum = error_sum + error
ALLOCATE(latitude_of_land_pts(land_pts), stat = error)
error_sum = error_sum + error
ALLOCATE(longitude_of_land_pts(land_pts), stat = error)
error_sum = error_sum + error
IF ( error_sum /= 0 )                                                         &
  CALL log_fatal("init_model_grid", "Error allocating land arrays")

land_mask(:,:) = ( flandg(:,:) > EPSILON(1.0) )

l = 0
DO j = 1,t_j_length
  DO i = 1,t_i_length
    IF ( land_mask(i,j) ) THEN
      l = l + 1
      land_index(l) = (j-1) * t_i_length + i
      fland(l) = flandg(i,j)
      latitude_of_land_pts(l)  = latitude(i,j)
      longitude_of_land_pts(l) = longitude(i,j)
    END IF
  END DO
END DO

!-----------------------------------------------------------------------------
! Allocate the rest of the JULES model arrays
!
! Note: CABLE model arrays are only allocated later on in init.F90 due to a
!       dependency on surft_pts, which has not yet been initialised
!-----------------------------------------------------------------------------
CALL allocate_jules_arrays()
CALL allocate_cable_progs() ! CABLE_LSM

!-----------------------------------------------------------------------------
! Now we know all the required info, set up the z-dimension sizes for IO
!-----------------------------------------------------------------------------
bl_level_dim_size = bl_levels
pft_dim_size    = npft
cpft_dim_size   = ncpft
nvg_dim_size    = nnvg
type_dim_size   = ntype
tile_dim_size   = nsurft
soil_dim_size   = sm_levels
scpool_dim_size = dim_cs1
sclayer_dim_size = dim_cslayer
snow_dim_size   = nsmax
tracer_dim_size = ndry_dep_species
bedrock_dim_size= ns_deep

! CABLE_LSM { including this from 8.5 but dont remember why - except that we do siilar in UM
land_pts_trif = land_pts
dim_cs2 = land_pts
!}
RETURN

END SUBROUTINE init_model_grid
#endif
