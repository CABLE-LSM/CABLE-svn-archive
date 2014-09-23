#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION internal_open_output_file(                                           &
  file_start, file_end, output_period, is_time_template, file_name, fields    &
) RESULT(file)

  USE io_constants, ONLY : MODE_WRITE, MAX_DIM_FILE, ATTR_GLOBAL

  USE dictionary_mod, ONLY : dict, dict_create, dict_free

  USE data_cube_mod, ONLY : cube_get_data, cube_free

  USE file_ts_mod, ONLY : file_ts_open, file_ts_def_grid, file_ts_def_time_dim,&
                          file_ts_def_attr, file_ts_enddef, file_ts_write_var

  USE model_interface_mod, ONLY : get_var_id, extract_var

  IMPLICIT NONE

  INCLUDE "mpif.h"

!-----------------------------------------------------------------------------
! Description:
!   Opens a file_ts object for a new output file using the supplied information
!   and returns it
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(datetime), INTENT(IN) :: file_start  ! The start time of the file
  TYPE(datetime), INTENT(IN) :: file_end    ! The end time of the file

  INTEGER, INTENT(IN) :: output_period  ! The output period of the file

  LOGICAL, INTENT(IN) :: is_time_template  ! Indicates if file_name includes
                                           ! time templating strings, or
                                           ! whether it is a single file name
  CHARACTER(len=MAX_FILE_NAME_LEN), INTENT(IN) :: file_name
                                           ! The name or time template to use
                                           ! for the file

  TYPE(output_field), INTENT(INOUT) :: fields(:)
                                           ! The fields we are defining in the
                                           ! file
                                           ! These are INOUT because the
                                           ! file_id field will be set to the
                                           ! id of the variable in the file

! Return type
  TYPE(file_ts) :: file


! Work variables
  INTEGER :: dummy  ! Dummy variable to receive dimension id of time dim

  TYPE(dict) :: defined_dims  ! Dictionary of dimensions defined in the file
                              ! so far

  TYPE(data_cube) :: lat_cube, lon_cube  ! Cubes for lat/lon data
  REAL :: point_lat, point_lon  ! Real values for lat/lon of single point
                                ! for ASCII files
  INTEGER :: lat_id, lon_id  ! Id of lat/lon variable in file for NetCDF files

  CHARACTER(len=15) :: cell_methods_val  ! Value of cell_methods attribute
                                         ! for the current variable

  INTEGER :: i  ! Loop variable


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Open a new file using the properties we were given
!-----------------------------------------------------------------------------
  file = file_ts_open(MODE_WRITE, file_start, file_end,                       &
                                  output_period, .FALSE.,                     &
! If is_time_template = T, file_name will be used as a time template
                                  is_time_template, file_name,                &
! If is_time_template = F, file_names and file_times will be used
! So we give a single file name starting at the start of the file (i.e. one
! file for all output)
                                  (/ file_name /), (/ file_start /),          &
! Pass MPI variables to file_ts_open so that parallel I/O will be used
                                  MPI_COMM_WORLD, MPI_INFO_NULL)


!-----------------------------------------------------------------------------
! Define the grid
!
! Before passing the grid_info object from output_mod.F90, overide the
! dimension names that will be used
!-----------------------------------------------------------------------------
  grid%dim_name = GRID_DIM_NAME
  grid%x_name   = X_DIM_NAME
  grid%y_name   = Y_DIM_NAME
  CALL file_ts_def_grid(file, grid)


!-----------------------------------------------------------------------------
! Define the time dimension
!-----------------------------------------------------------------------------
  dummy = file_ts_def_time_dim(file, TIME_DIM_NAME)


!-----------------------------------------------------------------------------
! Create a dictionary to store the levels dimensions defined so far
!    dim_name => dim_id
!-----------------------------------------------------------------------------
  defined_dims = dict_create(MAX_DIM_FILE, INT(1))


!-----------------------------------------------------------------------------
! Get the latitude and longitude data cubes
!-----------------------------------------------------------------------------
  lat_cube = extract_var(get_var_id('latitude'))
  lon_cube = extract_var(get_var_id('longitude'))


!-----------------------------------------------------------------------------
! Do the things that need to be done differently for ASCII and NetCDF
!
! In particular, if we are using ASCII files every variable must have a time
! dimension. However, the grid is restricted to 1 x 1 - so we use global
! attributes for latitude and longitude
!-----------------------------------------------------------------------------
  SELECT CASE ( OUTPUT_FORMAT )
    CASE ( FORMAT_ASCII )
! For ASCII, populate the lat/lon attributes
      CALL cube_get_data(lat_cube, point_lat)
      CALL cube_get_data(lon_cube, point_lon)

      CALL file_ts_def_attr(file, ATTR_GLOBAL, 'latitude', point_lat)
      CALL file_ts_def_attr(file, ATTR_GLOBAL, 'longitude', point_lon)

    CASE ( FORMAT_NCDF )
! For NetCDF, create non-time-varying lat/lon variables
      CALL internal_define_var(                                                 &
        file, defined_dims, get_var_id('latitude'), 'latitude', .FALSE., lat_id &
      )

      CALL internal_define_var(                                                   &
        file, defined_dims, get_var_id('longitude'), 'longitude', .FALSE., lon_id &
      )

! No default case, so that using a format other than those defined is a definite
! error
  END SELECT


!-----------------------------------------------------------------------------
! Set up the output variables in the file
!-----------------------------------------------------------------------------
  DO i = 1,SIZE(fields)

! Define the variable in the file (inc. dimensions and attributes)
    CALL internal_define_var(                                                 &
      file, defined_dims,                                                     &
      fields(i)%var_id, fields(i)%output_name, .TRUE., fields(i)%file_id      &
    )

! Add the CF convention "cell_methods" attribute to indicate whether the field
! is a snapshot, mean or accumulation
    SELECT CASE ( fields(i)%type )
      CASE ( OUTPUT_SNAPSHOT )
        cell_methods_val = "time : point"

      CASE ( OUTPUT_ACCUM )
        cell_methods_val = "time : sum"

      CASE ( OUTPUT_MEAN )
        cell_methods_val = "time : mean"

    END SELECT

    CALL file_ts_def_attr(                                                    &
      file, fields(i)%file_id, "cell_methods", cell_methods_val               &
    )

  END DO

! We have finished defining things on the file handle
  CALL file_ts_enddef(file)

! We no longer need the defined_dims dictionary
  CALL dict_free(defined_dims)


!-----------------------------------------------------------------------------
! Fill the latitude and longitude if we need to
!-----------------------------------------------------------------------------
  SELECT CASE ( OUTPUT_FORMAT )
    CASE ( FORMAT_ASCII )
! For ASCII files, we do nothing - we use a select statement so we get an
! error if someone has added a file type but not considered it's consequences
! here

    CASE ( FORMAT_NCDF )
      CALL file_ts_write_var(                                                 &
        file, lat_id, lat_cube,                                               &
! Subgrid information (for writing a slab of the output grid in parallel mode)
        use_subgrid, subgrid                                                  &
      )
      CALL file_ts_write_var(                                                 &
        file, lon_id, lon_cube,                                               &
! Subgrid information (for writing a slab of the output grid in parallel mode)
        use_subgrid, subgrid                                                  &
      )

! No default case, so that using a format other than those defined is a definite
! error
  END SELECT


!-----------------------------------------------------------------------------
! Deallocate the lat/lon data cubes
!-----------------------------------------------------------------------------
  CALL cube_free(lat_cube)
  CALL cube_free(lon_cube)


  RETURN

END FUNCTION internal_open_output_file
#endif
