#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/input/time_varying/register_input_file.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE register_input_file(data_start, data_end, data_period,             &
                               is_climatology, use_time_template, template,   &
                               file_names, file_times, identifiers,           &
                               sdf_names, interp_flags)

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, MODE_READ, MAX_SDF_NAME_LEN,    &
                           MAX_DIM_FILE

  USE model_time_mod, ONLY : run_min_time, run_max_time, timestep_len

  USE datetime_mod, ONLY : PERIOD_MONTH, PERIOD_YEAR, datetime_advance

  USE dictionary_mod, ONLY : dict, dict_create, dict_get, dict_set, dict_free

  USE file_ts_mod, ONLY : file_ts_open, file_ts_def_grid, file_ts_def_dim,    &
                          file_ts_def_time_dim, file_ts_def_var,              &
                          file_ts_enddef

  USE model_interface_mod, ONLY : get_var_id, get_var_levs_dim

  USE interpolation_mod, ONLY : get_required_time_bounds

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Registers a file or group of files as providing time-varying data for the
!   model variables specified by the given identifiers
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Argument types
!-----------------------------------------------------------------------------
  TYPE(datetime), INTENT(IN) :: data_start
                              ! The date and time of the first data
  TYPE(datetime), INTENT(IN) :: data_end
                              ! The date and time of the last data
  INTEGER, INTENT(IN) :: data_period
                              ! The period of the data
                              ! (in seconds or a 'special' period)
  LOGICAL, INTENT(IN) :: is_climatology
                              ! .TRUE. - the data is a climatology
                              ! .FALSE. - the data is not a climatology

  LOGICAL, INTENT(IN) :: use_time_template
                              ! .TRUE. - use time templating
                              ! .FALSE. - use lists of file names and times
                              !           of first data in each file

! With time templating
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: template
                              ! The time template to use

! With a file list
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: file_names(:)
                              ! List of file names
  TYPE(datetime), INTENT(IN), OPTIONAL :: file_times(:)
                              ! Time of first data for each file

  CHARACTER(len=*), INTENT(IN) :: identifiers(:)
                              ! The model identifiers of the variables
                              ! this file(s) provides data for
  CHARACTER(len=*), INTENT(IN) :: sdf_names(:)
                              ! The name of each variable in the file(s)
  CHARACTER(len=*), INTENT(IN) :: interp_flags(:)
                              ! The type of interpolation to use for each
                              ! variable


!-----------------------------------------------------------------------------
! Work variables
!-----------------------------------------------------------------------------
! Local variables passed to file_ts_open
! In the case where optional variables are not given, these are set to values
! that will cause errors if they are required, but will be ignored if they
! are not required
  CHARACTER(len=MAX_FILE_NAME_LEN) :: template_local
  CHARACTER(len=MAX_FILE_NAME_LEN), ALLOCATABLE :: file_names_local(:)
  TYPE(datetime), ALLOCATABLE :: file_times_local(:)

  INTEGER :: nvars  ! The number of variables in the file

  TYPE(dict) :: file_dims  ! Dictionary containing the dimensions to define
                           ! in the file
                           ! Maps dim_name => dim_size

  TYPE(dict) :: file_dim_ids  ! Dictionary containing the dimension ids in file
                              ! Maps dim_name => dim_id

  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_name  ! Temporary variables used to
  INTEGER :: dim_size                          ! populate dimension dictionaries
  INTEGER :: dim_id
  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_name_ignore  ! Used as a dummy argument

  INTEGER :: levels_dim_id  ! The id in the file of the vertical levels
                            ! dimension for the variable being processed
                            ! - used when defining variables


! Variables used to set up input_field objects
  INTEGER :: data_nx, data_ny  ! The number of points in the x and y dimension
                               ! of the model grid
  INTEGER :: var_ids(SIZE(identifiers))
                               ! The variable id from model_interface_mod
                               ! for each requested variable
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_levs_dims(SIZE(identifiers))
                               ! The name of the levels dimension for each
                               ! requested variable
  INTEGER :: var_nlevs(SIZE(identifiers))
                               ! The number of vertical levels for each
                               ! variable
  INTEGER :: var_file_ids(SIZE(identifiers))
                               ! The id of each variable in the file
  INTEGER :: times_lbound, times_ubound
                               ! The lower and upper bounds for the time
                               ! dimension, depending on what interpolation
                               ! is required


  INTEGER :: dummy  ! This is used to store the result from defining the time
                    ! dimension - we don't need to keep it

  INTEGER :: i  ! Loop counter


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  file_dims    = dict_create(MAX_DIM_FILE, INT(1))
  file_dim_ids = dict_create(MAX_DIM_FILE, INT(1))

!-----------------------------------------------------------------------------
! Check arguments
!-----------------------------------------------------------------------------
! Check that we have space to register another file
  IF ( nfiles >= NFILES_MAX )                                                 &
    CALL log_fatal("register_input_file",                                     &
                   "Too many files registered - try increasing NFILES_MAX")

! Check that we have an sdf_name and an interpolation flag for each variable
  nvars = SIZE(identifiers)
  IF ( nvars /= SIZE(sdf_names) )                                             &
    CALL log_fatal("register_input_file",                                     &
                   "identifiers and sdf_names must have the same number " //  &
                   "of elements")
  IF ( nvars /= SIZE(interp_flags) )                                          &
    CALL log_fatal("register_input_file",                                     &
                   "identifiers and interp_flags must have the same " //      &
                   "number of elements")

!-----------------------------------------------------------------------------
! Get the integer indentifiers from model_interface_mod for each requested
! string identifier
!-----------------------------------------------------------------------------
  DO i = 1,nvars
    var_ids(i) = get_var_id(identifiers(i))
  END DO

!-----------------------------------------------------------------------------
! Check that data will be provided for the whole run at an appropriate
! timestep
!-----------------------------------------------------------------------------
! If not using a climatology, data must be provided for the whole run
! If using a climatology, the data must cover a whole year, but this is checked
! by file_ts_open
  IF ( .NOT. is_climatology .AND.                                             &
       ( data_start > run_min_time .OR. data_end < run_max_time ) )           &
    CALL log_fatal("register_input_file",                                     &
                   "Each input file must provide data for the entire run")

! Unless data is given on a 'special' period, the data period must be a whole
! number of model timesteps
  IF ( data_period /= PERIOD_MONTH .AND. data_period /= PERIOD_YEAR .AND.     &
       MOD(data_period, timestep_len) /= 0 )                                  &
    CALL log_fatal("register_input_file",                                     &
                   "Data period must be a special period or a multiple " //   &
                   "of timestep length")

!-----------------------------------------------------------------------------
! Set up the optional arguments depending on what is available
!-----------------------------------------------------------------------------
  IF ( PRESENT(template) ) THEN
    template_local = template
  ELSE
! If template is not provided, we provide an empty string
! This will cause an error unless a file list is specified
    template_local = ""
  END IF

  IF ( PRESENT(file_names) ) THEN
! If file_names is present, copy it's values into the local counterpart
! that has fixed length strings
    ALLOCATE(file_names_local(SIZE(file_names)))
    file_names_local(:) = file_names(:)
  ELSE
! If file_names is not provided, we provide one file that doesn't exist
! This will cause an error unless time templating is specified
    ALLOCATE(file_names_local(1))
    file_names_local(1) = ""
  END IF

  IF ( PRESENT(file_times) ) THEN
! If file_times is present, copy it's values into the local counterpart
    ALLOCATE(file_times_local(SIZE(file_times)))
    file_times_local(:) = file_times(:)
  ELSE
! If file_times is not provided, we provide one time that is not equal to
! data_start
! This will cause an error unless time templating is specified
    ALLOCATE(file_times_local(1))
    file_times_local(1) = datetime_advance(data_start, 1)
  END IF

!-----------------------------------------------------------------------------
! Open the file handle
!-----------------------------------------------------------------------------
  nfiles = nfiles + 1

  files(nfiles)%fh = file_ts_open(MODE_READ, data_start, data_end,            &
                                  data_period, is_climatology,                &
                                  use_time_template, template_local,          &
                                  file_names_local, file_times_local)

! We have finished with the local pointers
  DEALLOCATE(file_names_local)
  DEALLOCATE(file_times_local)

!-----------------------------------------------------------------------------
! Define the grid
!-----------------------------------------------------------------------------
  CALL file_ts_def_grid(                                                      &
    files(nfiles)%fh, grid_is_1d, grid_dim_name, nx,                          &
                                  x_dim_name, nx, y_dim_name, ny              &
  )

!-----------------------------------------------------------------------------
! Define the required levels dimensions
!-----------------------------------------------------------------------------
! Make a pass over the requested variables to work out what levels dimensions
! we expect to have in the file
  DO i = 1,nvars
! We use the input name and ignore the output name for dimensions
    CALL get_var_levs_dim(var_ids(i), dim_name, dim_name_ignore, dim_size)

! Store the levels dim name and size for this variable for use later
    var_levs_dims(i) = dim_name
    var_nlevs(i) = dim_size

! For a variable with no levels dimension, the returned dim_name will be the
! empty string - we don't want to add this to the dictionary
    IF ( LEN_TRIM(dim_name) < 1 ) CYCLE

! Otherwise store it in the dictionary
! Duplicate keys will overwritten, so we don't have to worry about duplication
    CALL dict_set(file_dims, dim_name, dim_size)
  END DO

! Define the required dimensions
  DO i = 1,file_dims%length
! The dimension names are the dictionary keys, and the sizes are the values
    dim_name = file_dims%keys(i)
    CALL dict_get(file_dims, dim_name, dim_size)

! Define the dimension
    dim_id = file_ts_def_dim(files(nfiles)%fh, dim_name, dim_size)

! Store the id
    CALL dict_set(file_dim_ids, dim_name, dim_id)
  END DO

!-----------------------------------------------------------------------------
! Define the time dimension
!-----------------------------------------------------------------------------
  dummy = file_ts_def_time_dim(files(nfiles)%fh, time_dim_name)

!-----------------------------------------------------------------------------
! Define the variables and save their ids for later
!-----------------------------------------------------------------------------
  DO i = 1,nvars
! Get the dimension id to use for the levels dim for this variable
    IF ( LEN_TRIM(var_levs_dims(i)) < 1 ) THEN
! If an empty string was returned as the levels dim name, then the variable
! has no levels dimension
! We exploit the fact that giving -1 as the dimension id is equivalent to not
! giving a dimension id
      levels_dim_id = -1
    ELSE
! Otherwise retrieve the dimension id from the dictionary
      CALL dict_get(file_dim_ids, var_levs_dims(i), levels_dim_id)
    END IF

! Define the variable and store its id
    var_file_ids(i) = file_ts_def_var(                                        &
      files(nfiles)%fh, sdf_names(i), levels_dim_id, .TRUE.                   &
    )
  END DO

! We have finished defining things on the file handle
  CALL file_ts_enddef(files(nfiles)%fh)

  CALL dict_free(file_dims)
  CALL dict_free(file_dim_ids)


!-----------------------------------------------------------------------------
! Set up the input_file object and its corresponding input_field objects
!-----------------------------------------------------------------------------
! Calculate values required to dimension the data arrays of the input fields
! Get the size of the grid dimensions
  data_nx = nx
  data_ny = ny
  IF ( use_subgrid ) THEN
    data_nx = subgrid_npoints
    data_ny = 1
  END IF
! Get the upper and lower bounds to use for the time dimension based on what
! interpolation each variable is using
  CALL get_required_time_bounds(interp_flags, times_lbound, times_ubound)
  files(nfiles)%times_lbound = times_lbound
  files(nfiles)%times_ubound = times_ubound

! Set up input fields
! Allocate space on the input_file object for the fields
  files(nfiles)%nfields = nvars
  ALLOCATE(files(nfiles)%fields(nvars))

  DO i = 1,nvars
! Allocate space for the data
    ALLOCATE(files(nfiles)%fields(i)%data(data_nx, data_ny, var_nlevs(i),     &
                                          times_lbound:times_ubound))

! Store other information required later
    files(nfiles)%fields(i)%var_id      = var_ids(i)
    files(nfiles)%fields(i)%file_id     = var_file_ids(i)
    files(nfiles)%fields(i)%interp_flag = interp_flags(i)
  END DO

  RETURN

END SUBROUTINE register_input_file
#endif
