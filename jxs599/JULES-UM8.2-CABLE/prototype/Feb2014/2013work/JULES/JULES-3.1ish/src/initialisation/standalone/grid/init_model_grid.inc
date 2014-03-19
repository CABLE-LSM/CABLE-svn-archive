#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/grid/init_model_grid.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_model_grid()

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, NAMELIST_UNIT, POINTS_FILE_UNIT

  USE string_utils_mod, ONLY : to_string

  USE input_mod, ONLY : grid_is_1d, nx_in => nx, ny_in => ny,                 &
! Variables to do with the extraction of subgrids if specified
                        use_subgrid, subgrid_npoints, points, x_coords, y_coords

  USE switches, ONLY : l_co2_interactive, l_point_data, l_triffid

  USE latlon_mod, ONLY : latitude, longitude

  USE coastal, ONLY : fland, flandg

  USE atm_fields_bounds_mod, ONLY : atm_fields_bounds_init
  USE theta_field_sizes, ONLY : t_i_length, t_j_length

  USE ancil_info, ONLY : land_pts, land_mask, land_index, co2_dim_len,        &
                         co2_dim_row, dim_cs2, land_pts_trif, row_length,     &
                         rows, n_rows

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about how the model grid is defined, sets it up and
!   checks for consistency
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  LOGICAL :: points_mask(nx_in, ny_in)  ! Defined on the full input grid
                                        !   T - the point will be modelled
                                        !   F - the point will not be modelled
  REAL :: point_lat, point_lon
                                ! Latitude and longitude of points read from
                                ! points file

  REAL, ALLOCATABLE :: subgrid_data(:,:)  ! Used when regridding latitude,
                                          ! longitude and land fraction if
                                          ! required

  INTEGER :: error, error_sum  ! Error indicators

  INTEGER :: i, j, l, p  ! Index variables

!-----------------------------------------------------------------------------
! Definition of the jules_model_grid namelist - this combines variables
! from input_mod with some local variables
!-----------------------------------------------------------------------------
  LOGICAL :: land_only      ! T - only model land points
                            ! F - model all specified points

  LOGICAL :: latlon_region  ! T - subgrid is to be selected with latitude and
                            !     longitude bounds
                            ! F - subgrid is to be selected using a list
                            !     of latitudes and longitudes
  REAL :: lat_bounds(2)   ! USED IF latlon_region=T
  REAL :: lon_bounds(2)   ! Upper and lower bounds for latitude and longitude

  INTEGER :: npoints = 0 ! The number of points to read from file
  CHARACTER(len=MAX_FILE_NAME_LEN) :: points_file
                          ! The file to read latitudes and longitudes for
                          ! specified points

  NAMELIST /jules_model_grid/ land_only, use_subgrid, latlon_region,          &
                              lat_bounds, lon_bounds, npoints, points_file


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  land_only = .TRUE.  ! Default is to model only land points from whatever
                      ! points are selected

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_model_grid", "Reading JULES_MODEL_GRID namelist...")

! First, we read the namelist
  READ(NAMELIST_UNIT, nml=jules_model_grid, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_model_grid",                                         &
                   "Error reading namelist JULES_MODEL_GRID " //              &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Set up the model grid variables
!-----------------------------------------------------------------------------
  CALL log_info("init_model_grid", "Setting up model grid variables...")

!*****************************************************************************
! Work out what points in the input grid will make up the model grid
!
! By the time we get here, we know that we have latitude, longitude and
! land fraction on the full input grid already read in
! Use these in combination with the values in the jules_model_grid namelist
! to work out what points to use
!*****************************************************************************
! To start with, assume we are using all the points
  points_mask(:,:) = .TRUE.

  IF ( use_subgrid ) THEN
    IF ( latlon_region ) THEN
!-----------------------------------------------------------------------------
! The subgrid will be selected using latitude and longitude bounds
!-----------------------------------------------------------------------------
      CALL log_info("init_model_grid",                                        &
                    "Subgrid will be selected using latitude and " //         &
                    "longitude bounds")
      CALL log_info("init_model_grid",                                        &
                    "Latitude range - " // TRIM(to_string(lat_bounds(1))) //  &
                    " to " // TRIM(to_string(lat_bounds(2))))
      CALL log_info("init_model_grid",                                        &
                    "Longitude range - " // TRIM(to_string(lon_bounds(1))) // &
                    " to " // TRIM(to_string(lon_bounds(2))))

      points_mask = ( lat_bounds(1) <= latitude ) .AND.                       &
                    ( latitude <= lat_bounds(2) ) .AND.                       &
                    ( lon_bounds(1) <= longitude ) .AND.                      &
                    ( longitude <= lon_bounds(2) )
    ELSE
!-----------------------------------------------------------------------------
! The subgrid will be selected using a list of latitudes and longitudes
! from the given file
!-----------------------------------------------------------------------------
      CALL log_info("init_model_grid",                                        &
                    "Subgrid will be selected using a list of points " //     &
                    "from " // TRIM(points_file))

! No points have been specified yet, so set points_mask accordingly
      points_mask(:,:) = .FALSE.

!-----------------------------------------------------------------------------
! Read the point latitudes and longitudes from file
!-----------------------------------------------------------------------------
      OPEN(POINTS_FILE_UNIT, FILE=points_file,                                &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
      IF ( error /= 0 )                                                       &
      CALL log_fatal("init_model_grid",                                       &
                     "Error opening points file " //                          &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")

      DO i = 1,npoints
! Try to read the next latitude/longitude pair
        READ(POINTS_FILE_UNIT, *, IOSTAT=error) point_lat, point_lon
        IF ( error /= 0 )                                                     &
          CALL log_fatal("init_model_grid",                                   &
                         "Error reading lat/lon pair from points file " //    &
                         "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Set points_mask to .TRUE. at that point
        WHERE ( ABS(latitude - point_lat) < EPSILON(1.0) .AND.                &
                ABS(longitude - point_lon) < EPSILON(1.0) )
          points_mask = .TRUE.
        END WHERE
      END DO

      CLOSE(POINTS_FILE_UNIT, IOSTAT=error)
      IF ( error /= 0 )                                                       &
        CALL log_fatal("init_model_grid",                                     &
                       "Error closing points file " //                        &
                       "(IOSTAT=" // TRIM(to_string(error)) // ")")
    END IF  ! latlon_region
  END IF  ! use_subgrid

  IF ( land_only ) THEN
!-----------------------------------------------------------------------------
! If requested, select land points only from the points we have left
!-----------------------------------------------------------------------------
    CALL log_info("init_model_grid",                                          &
                  "From the points specified, only land points will be " //   &
                  "modelled")

    points_mask = points_mask .AND. ( flandg > 0.0 )
  END IF


!-----------------------------------------------------------------------------
! The model grid is the points that are left - configure the variables in
! input_mod so that it will be extracted correctly
!-----------------------------------------------------------------------------
  subgrid_npoints = COUNT(points_mask)

  IF ( subgrid_npoints <= 0 )                                                 &
    CALL log_fatal("init_model_grid",                                         &
                   "All points in input grid have been excluded from " //     &
                   "model grid")

! If the number of .TRUE. points in points_mask is < the number of points in
! the input grid, then we are using a subgrid
  use_subgrid = ( subgrid_npoints < (nx_in * ny_in) )

  IF ( use_subgrid ) THEN
! Set up points/x_coords/y_coords depending on whether the input grid is 1d
! or 2d
    ALLOCATE(x_coords(subgrid_npoints), STAT=error)
    error_sum = error
    ALLOCATE(y_coords(subgrid_npoints), STAT=error)
    error_sum = error_sum + error
    IF ( error_sum /= 0 )                                                     &
      CALL log_fatal("init_model_grid",                                       &
                     "Error allocating subgrid coordinate arrays")

    p = 0
    DO j = 1,ny_in
      DO i = 1,nx_in
        IF ( points_mask(i,j) ) THEN
          p = p + 1
          x_coords(p) = i
          y_coords(p) = j
        END IF
      END DO
    END DO

! Regrid the latitude, longitude and land fraction onto the subgrid
    ALLOCATE(subgrid_data(subgrid_npoints,1), STAT=error)
    error_sum = error

! Regrid the latitude
    DO p = 1,subgrid_npoints
      subgrid_data(p,1) = latitude(x_coords(p),y_coords(p))
    END DO
    DEALLOCATE(latitude, STAT=error)
    error_sum = error_sum + error
    ALLOCATE(latitude(subgrid_npoints,1), STAT=error)
    error_sum = error_sum + error
    latitude(:,:) = subgrid_data(:,:)

! Regrid the longitude
    DO p = 1,subgrid_npoints
      subgrid_data(p,1) = longitude(x_coords(p),y_coords(p))
    END DO
    DEALLOCATE(longitude, STAT=error)
    error_sum = error_sum + error
    ALLOCATE(longitude(subgrid_npoints,1), STAT=error)
    error_sum = error_sum + error
    longitude(:,:) = subgrid_data(:,:)

! Regrid the land fraction
    DO p = 1,subgrid_npoints
      subgrid_data(p,1) = flandg(x_coords(p),y_coords(p))
    END DO
    DEALLOCATE(flandg, STAT=error)
    error_sum = error_sum + error
    ALLOCATE(flandg(subgrid_npoints,1), STAT=error)
    error_sum = error_sum + error
    flandg(:,:) = subgrid_data(:,:)

    DEALLOCATE(subgrid_data, STAT=error)
    error_sum = error_sum + error

! Since a 1d grid is just a 2d grid with nx=npoints and ny=1, we can just
! swap points for x_coords if using a 1d grid
    IF ( grid_is_1d ) THEN
      ALLOCATE(points(subgrid_npoints), STAT=error)
      error_sum = error_sum + error
      points(:) = x_coords(:)
      DEALLOCATE(x_coords, STAT=error)
      error_sum = error_sum + error
      DEALLOCATE(y_coords, STAT=error)
      error_sum = error_sum + error
    END IF

    IF ( error_sum /= 0 )                                                     &
      CALL log_fatal("init_model_grid", "Error reallocating arrays to subgrid")
  END IF  ! use_subgrid


!-----------------------------------------------------------------------------
! Set up the grid sizes for the model grid and output grid
!-----------------------------------------------------------------------------
! Set up the structures containing the dimension sizes for the model variables
  IF ( use_subgrid ) THEN
    t_i_length = subgrid_npoints
    t_j_length = 1
  ELSE
    t_i_length = nx_in
    t_j_length = ny_in
  END IF
! At the moment, JULES uses no halos, and the t, p, u and v grids are the same
  CALL atm_fields_bounds_init(0, 0, 0, 0, t_i_length, t_j_length, t_j_length, 1)
! Copy the values into row_length, rows and n_rows for legacy purposes
  row_length = t_i_length
  rows   = t_j_length
  n_rows = t_j_length

  CALL log_info("init_model_grid",                                            &
                "Size of model grid - " // TRIM(to_string(t_i_length)) //     &
                " x " // TRIM(to_string(t_j_length)))

! Calculate the number of land points
  land_pts = COUNT(flandg > 0.0)
  CALL log_info("init_model_grid",                                            &
                "Selected grid contains " // TRIM(to_string(land_pts)) //     &
                " land points")

! Set up other model grid information
  IF ( l_triffid ) THEN
    land_pts_trif = land_pts
    dim_cs2 = land_pts
  ELSE
    land_pts_trif = 1
    dim_cs2 = 1
  ENDIF

  IF( l_co2_interactive ) THEN
    co2_dim_len = t_i_length
    co2_dim_row = t_j_length
  ELSE
    co2_dim_len = 1
    co2_dim_row = 1
  ENDIF

! Warn if using l_point_data for a non-single-point run
  IF ( l_point_data .AND. t_i_length * t_j_length > 1 )                       &
    CALL log_warn("init_model_grid",                                          &
                  "l_point_data is selected but there is more than one point")

!-----------------------------------------------------------------------------
! Set up the mapping of land points to the full model grid
!-----------------------------------------------------------------------------
  ALLOCATE(land_mask(t_i_length,t_j_length), STAT=error)
  error_sum = error
  ALLOCATE(land_index(land_pts), STAT=error)
  error_sum = error_sum + error
  ALLOCATE(fland(land_pts), STAT=error)
  error_sum = error_sum + error
  IF ( error_sum /= 0 )                                                       &
    CALL log_fatal("init_model_grid", "Error allocating land arrays")

  land_mask = ( flandg > 0.0 )

  l = 0
  DO j = 1,t_j_length
    DO i = 1,t_i_length
      IF ( land_mask(i,j) ) THEN
        l = l + 1
        land_index(l) = (j - 1) * t_i_length + i
        fland(l) = flandg(i,j)
      END IF
    END DO
  END DO

!-----------------------------------------------------------------------------
! Allocate the rest of the model arrays
!-----------------------------------------------------------------------------
  CALL allocate_jules_arrays()

  RETURN

END SUBROUTINE init_model_grid
#endif
