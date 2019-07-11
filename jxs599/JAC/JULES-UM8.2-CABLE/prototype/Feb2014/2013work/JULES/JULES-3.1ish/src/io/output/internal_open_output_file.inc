#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/output/internal_open_output_file.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE internal_open_output_file(profile)

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, FORMAT_ASCII, FORMAT_NCDF,      &
                           MODE_WRITE, ATTR_GLOBAL

  USE datetime_mod, ONLY : PERIOD_MONTH, PERIOD_YEAR, datetime_advance

  USE templating_mod, ONLY : TPL_YR_4DIGIT, TPL_MON_2DIGIT

  USE model_time_mod, ONLY : is_spinup, spinup_start, spinup_end, spinup_cycle

  USE model_interface_mod, ONLY : extract_var, get_var_id


  USE file_ts_mod, ONLY : file_ts_open, file_ts_def_attr, file_ts_enddef,     &
                          file_ts_write_var, file_ts_close

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a file_ts object representing a new output file with name(s) based
!   on the current model time and returns it
!   This routine does not do any checks to see if the profile will use the
!   file - these are done in the calling routine
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(output_profile), INTENT(INOUT) :: profile  ! The profile to open a new file for

! Work variables
! Local variables passed to file_ts_open
! We will use time templating, except when one file contains output for the
! whole period, in which case we specify no time templating with one file name
! and one file time
  TYPE(datetime) :: file_start
  TYPE(datetime) :: file_end

  LOGICAL :: use_template
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file_name

! Variables used to set up latitude and longitude
  REAL :: lat_data(nx,ny,1)  ! Latitudes on the output grid
  REAL :: lon_data(nx,ny,1)  ! Longitudes on the output grid
  INTEGER :: lat_var_id, lon_var_id  ! The variable ids in the defined file
                                     ! for latitude and longitude

! Variables used to define variables in file and set up output fields
  INTEGER :: var_ids(profile%nfields)
                     ! Stores the variable ids of the variables to create
                     ! in the file
                     ! See model_interface_mod
  INTEGER :: var_ids_tmp(profile%nfields+2)
                     ! Temporary version of var_ids for when lat and
                     ! lon ids are added
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_names(profile%nfields)
                     ! Stores the name to use in output files for each
                     ! variable
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_names_tmp(profile%nfields+2)
                     ! Temporary version of var_names for when lat and
                     ! lon names are added
  INTEGER :: var_file_ids(profile%nfields)
                    ! Stores the variable id in the file for each variable
  INTEGER :: var_file_ids_tmp(profile%nfields+2)
                    ! Temporarily stores the output from the call to
                    ! internal_define_vars until it can be split out into
                    ! lat_var_id and lon_var_id (if required) and var_file_ids
  INTEGER :: var_nlevs(profile%nfields)
                    ! Stores the number of vertical levels for each variable
  INTEGER :: var_nlevs_tmp(profile%nfields+2)
                    ! Analog of var_ids_tmp for var_nlevs


  INTEGER :: i  ! Loop counters


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Close the current file if we need to
!-----------------------------------------------------------------------------
  IF ( profile%has_open_file ) CALL file_ts_close(profile%fh)

!-----------------------------------------------------------------------------
! Work out when the data in the file will start and end
!-----------------------------------------------------------------------------
  IF ( is_spinup ) THEN
    file_start = spinup_start
    file_end   = spinup_end
  ELSE
    file_start = profile%output_start
    file_end   = profile%output_end
  END IF

!-----------------------------------------------------------------------------
! Work out what file name/template we want to use
!-----------------------------------------------------------------------------
! All file names start with the run id followed by the profile name
  file_name = TRIM(run_id) // "." // TRIM(profile%profile_name)

  IF ( is_spinup ) THEN
! If we are in a spinup cycle, add that to the file name
    file_name = TRIM(file_name) // "." // "spin" // TRIM(to_string(spinup_cycle))
  END IF

! Indicate whether we will be using time templating or not and add the required
! time template to the file name
  use_template = .FALSE.

  IF ( profile%file_period == PERIOD_YEAR .OR.                                &
       profile%file_period == PERIOD_MONTH ) THEN
    use_template = .TRUE.

    file_name = TRIM(file_name) // "." // TPL_YR_4DIGIT
    IF ( profile%file_period == PERIOD_MONTH )                                &
      file_name = TRIM(file_name) // TPL_MON_2DIGIT
  END IF

! Add an extension based on the output format
  SELECT CASE ( OUTPUT_FORMAT )
    CASE ( FORMAT_ASCII )
      file_name = TRIM(file_name) // ".asc"

    CASE ( FORMAT_NCDF )
      file_name = TRIM(file_name) // ".nc"

    CASE DEFAULT
      CALL log_fatal("internal_open_output_file",                             &
                     "Unrecognised file format - " // TRIM(OUTPUT_FORMAT))
  END SELECT

! Prepend the output directory
  file_name = TRIM(output_dir) // "/" // TRIM(file_name)

!-----------------------------------------------------------------------------
! Open a new file using the properties we have gathered
!-----------------------------------------------------------------------------
  profile%fh = file_ts_open(MODE_WRITE, file_start, file_end,                 &
                                        profile%output_period, .FALSE.,       &
                                        use_template, file_name,              &
! For file_names and file_times we give a single file starting at the output
! start time
! These will only be used in the case of no time templating
                                        (/ file_name /), (/ file_start /))

!-----------------------------------------------------------------------------
! Define the variables (and dimensions/attributes)
!-----------------------------------------------------------------------------
! Build the var_ids and var_names lists from the output profiles
  var_ids(:)   = profile%fields%var_id
  var_names(:) = profile%fields%output_name

  SELECT CASE ( OUTPUT_FORMAT )
    CASE ( FORMAT_ASCII )
! If we are using ASCII files, every variable must have a time dimension but
! the grid is restricted to 1 x 1 - so we use global attributes
! for latitude and longitude
      CALL extract_var(get_var_id('latitude'), lat_data)
      CALL extract_var(get_var_id('longitude'), lon_data)

      CALL file_ts_def_attr(profile%fh, ATTR_GLOBAL, 'latitude', lat_data(1,1,1))
      CALL file_ts_def_attr(profile%fh, ATTR_GLOBAL, 'longitude', lon_data(1,1,1))

! Define the rest of the variables normally, all as using the time dimension
      CALL internal_define_vars(                                                &
        profile%fh,                                                             &
        var_ids, var_names, (/ (.TRUE., i = 1,profile%nfields) /),              &
        var_file_ids, var_nlevs                                                 &
      )

    CASE ( FORMAT_NCDF )
! For a NetCDF file, a gridded format, latitude and longitude are defined as
! variables with no time dimension
! We define all the variables together by adding items for latitude and longitude
! to temporary arrays
      var_ids_tmp(1) = get_var_id('latitude')
      var_ids_tmp(2) = get_var_id('longitude')
      var_ids_tmp(3:) = var_ids(:)

      var_names_tmp(1) = 'latitude'
      var_names_tmp(2) = 'longitude'
      var_names_tmp(3:) = var_names(:)

      CALL internal_define_vars(                                                &
        profile%fh, var_ids_tmp, var_names_tmp,                                 &
! lat and lon don't use time, the rest do
        (/ (.FALSE., i = 1,2), (.TRUE., i = 1,profile%nfields) /),              &
        var_file_ids_tmp, var_nlevs_tmp                                         &
      )

! Extract the ids and levels for latitude, longitude and the rest of the
! variables
      lat_var_id = var_file_ids_tmp(1)
      lon_var_id = var_file_ids_tmp(2)
      var_file_ids(:) = var_file_ids_tmp(3:)
      var_nlevs(:) = var_nlevs_tmp(3:)

! No default case, so that using a format other than those defined is a definite
! error
  END SELECT

! We have finished defining things on the file handle
  CALL file_ts_enddef(profile%fh)

!-----------------------------------------------------------------------------
! Fill the latitude and longitude if we need to
!-----------------------------------------------------------------------------
  SELECT CASE ( OUTPUT_FORMAT )
    CASE ( FORMAT_ASCII )
! For ASCII files, we do nothing - we use a select statement so we get an
! error if someone has added a file type but not considered it's consequences
! here

    CASE ( FORMAT_NCDF )
      CALL extract_var(get_var_id('latitude'), lat_data)
      CALL extract_var(get_var_id('longitude'), lon_data)

! The final argument is given as F since we don't use subgrids for output
! (at the moment)
      CALL file_ts_write_var(profile%fh, lat_var_id, lat_data, .FALSE.)
      CALL file_ts_write_var(profile%fh, lon_var_id, lon_data, .FALSE.)

! No default case, so that using a format other than those defined is a definite
! error
  END SELECT

!-----------------------------------------------------------------------------
! Initialise the profile so that the sampling and output periods start at
! the beginning of the newly opened file
!-----------------------------------------------------------------------------
  profile%next_sample_time    = file_start
  profile%current_output_time = file_start
  profile%next_output_time    = datetime_advance(file_start, profile%output_period)
  profile%samples_in_period   = 0

!-----------------------------------------------------------------------------
! Allocate the data fields on the output_field objects
! We do this here so they don't get allocated unless they are used
! Also set the variable id that each field uses in the current file
!-----------------------------------------------------------------------------
  DO i = 1,profile%nfields
! Allocate space for the data
    IF ( .NOT. ASSOCIATED(profile%fields(i)%data) )                           &
      ALLOCATE(profile%fields(i)%data(nx, ny, var_nlevs(i)))

    profile%fields(i)%data(:,:,:) = 0

    profile%fields(i)%file_id = var_file_ids(i)
  END DO

! Indicate that the profile now has an open file
  profile%has_open_file = .TRUE.

  RETURN

END SUBROUTINE internal_open_output_file
#endif
