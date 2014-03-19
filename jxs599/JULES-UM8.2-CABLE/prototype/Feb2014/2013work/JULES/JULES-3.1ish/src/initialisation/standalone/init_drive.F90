#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/init_drive.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

SUBROUTINE init_drive()

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN, NAMELIST_UNIT

  USE model_interface_mod, ONLY : IDENTIFIER_LEN

  USE datetime_mod, ONLY : DATETIME_STR_LEN, datetime, datetime_from_string

  USE string_utils_mod, ONLY : to_string

  USE templating_mod, ONLY : tpl_detect_period, tpl_has_var_name,             &
                             tpl_substitute_var

  USE update_mod, ONLY : io_precip_type, io_rad_type, io_wind_speed,          &
                         use_diff_rad, diff_frac_const, t_for_con_rain,       &
                         t_for_snow

  USE time_varying_input_mod, ONLY : register_input_file

  USE switches, ONLY : l_imogen, l_point_data

  USE ancil_info, ONLY : z1_uv, z1_tq

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about driving data and initialises it
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER, PARAMETER :: MAX_DRIVE_VARS = 20

  INTEGER :: nvars_required      ! The number of variables that are
                                 ! required in the current configuration
  CHARACTER(len=IDENTIFIER_LEN) :: required_vars(MAX_DRIVE_VARS)
                                 ! The variable identifiers of the required
                                 ! variables

  INTEGER :: nvars_file  ! The number of files to be read from file

  TYPE(datetime) :: data_start_dt, data_end_dt
                              ! Start and end times for data as datetime objects

  LOGICAL :: use_time_tpl     ! Determines whether we will be using time
                              ! templating or whether we will pass lists
                              ! of files and times to register_input_file

  LOGICAL :: use_var_name_tpl  ! Determines whether we will be using variable
                               ! name templating or not

  CHARACTER(len=MAX_FILE_NAME_LEN), ALLOCATABLE :: file_names(:)
                              ! If use_time_tpl=T, the list of file names
  TYPE(datetime), ALLOCATABLE :: file_times(:)
                              ! If use_time_tpl=T, the list of file start times
  CHARACTER(len=DATETIME_STR_LEN) :: file_time_str
                              ! The read file time as a string until it can
                              ! be converted

  INTEGER :: i,j  ! Index variables

  INTEGER :: error, error_sum  ! Error indicators

!-----------------------------------------------------------------------------
! Definition of the jules_drive namelist
!-----------------------------------------------------------------------------
! Information about the data in the file
  CHARACTER(len=DATETIME_STR_LEN) :: data_start, data_end
                             ! Start and end times for driving data as strings
  INTEGER :: data_period     ! The period of the driving data

! Properties that determine how dependent variables are updated
  REAL :: z1_uv_in, z1_tq_in
      ! values of z1_uv and z1_tq that are used at all points on grid

! Information about the files to use
  LOGICAL :: read_list       ! T - the given file contains a list of file
                             !     names and times of first data
                             !     These files may contain variable name
                             !     templating but not time templating
                             ! F - driving data should be read directly from
                             !     the given file
                             !     The file may contain variable and time
                             !     templating
  INTEGER :: nfiles          ! The number of files/file times to read from the
                             ! list file
                             ! ONLY USED IF read_list=T
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
                             ! The file to use for whatever read_list
                             ! determines

! Information about the variables contained in the file
  INTEGER :: nvars
  CHARACTER(len=IDENTIFIER_LEN) :: var(MAX_DRIVE_VARS)
                        ! The variable identifiers of the variables
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_name(MAX_DRIVE_VARS)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: tpl_name(MAX_DRIVE_VARS)
                        ! The name to substitute in a template for each
                        ! variable
  CHARACTER(len=2) :: interp(MAX_DRIVE_VARS)
                        ! Flag indicating the type of interpolation to use
                        ! for each variable


  NAMELIST /jules_drive/ t_for_snow, t_for_con_rain, diff_frac_const,         &
                         z1_uv_in, z1_tq_in,                                  &
                         data_start, data_end, data_period,                   &
                         read_list, nfiles, file,                             &
                         nvars, var, var_name, tpl_name, interp


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  nvars_required = 0
  nvars_file     = 0
  read_list      = .FALSE.  ! Default is to assume that the given file name
                            ! either provides the data or is a template
  nfiles         = 0
  nvars          = 0
  z1_uv_in       = 10.0
  z1_tq_in       = 10.0

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_drive", "Reading JULES_DRIVE namelist...")

  OPEN(NAMELIST_UNIT, FILE='drive.nml',                                       &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_drive",                                              &
                   "Error opening namelist file drive.nml " //                &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  READ(NAMELIST_UNIT, nml=jules_drive, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_drive",                                              &
                   "Error reading namelist JULES_DRIVE " //                   &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_drive",                                              &
                   "Error closing namelist file drive.nml " //                &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Set the reference level heights for wind and air temperature/humidity to
! user-specified constant values.
!-----------------------------------------------------------------------------
  z1_uv(:,:) = z1_uv_in
  z1_tq(:,:) = z1_tq_in

!-----------------------------------------------------------------------------
! Check if IMOGEN will be used to update the main driving variables
!-----------------------------------------------------------------------------
  IF ( l_imogen ) THEN
    CALL log_info("init_drive", "Driving data will be provided by IMOGEN")

! Using IMOGEN to update driving variables is conceptually the same as having
! sw_down and lw_down provided (io_rad_type = 1), all 4 precip components
! provided (io_precip_type = 4) and both components of wind provided
! (io_wind_speed = F)
! We also use a constant diffuse fraction
    io_rad_type = 1
    io_precip_type = 4
    io_wind_speed = .FALSE.
    use_diff_rad = .FALSE.

    RETURN
  END IF

!-----------------------------------------------------------------------------
! Convert data_start and data_end to datetime objects
!-----------------------------------------------------------------------------
  data_start_dt = datetime_from_string(data_start)
  data_end_dt   = datetime_from_string(data_end)

!-----------------------------------------------------------------------------
! Work out what files we will be using to read driving data
!-----------------------------------------------------------------------------
  IF ( read_list ) THEN
    IF ( nfiles < 1 )                                                         &
      CALL log_fatal("init_drive",                                            &
                     "If reading a list of file names and file times, at " // &
                     "least one file must be given")

! If we are reading a list of files, then we will definitely not be using
! time templating
    use_time_tpl = .FALSE.

!-----------------------------------------------------------------------------
! Read the list of file names and times from the file
!-----------------------------------------------------------------------------
    CALL log_info("init_drive",                                               &
                  "Reading list of drive file names and start times...")

! First allocate space for the elements
    ALLOCATE(file_names(nfiles), STAT=error)
    error_sum = error
    ALLOCATE(file_times(nfiles), STAT=error)
    error_sum = error_sum + error
    IF ( error_sum /= 0 )                                                     &
      CALL log_fatal("init_drive",                                            &
                     "Error allocating arrays for file names and times")

! Open the file
    OPEN(NAMELIST_UNIT, FILE=file, STATUS='old', POSITION='rewind',           &
                                                 ACTION='read', IOSTAT=error)
    IF ( error /= 0 )                                                         &
      CALL log_fatal("init_drive",                                            &
                     "Error opening file " // TRIM(file) // " " //            &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Read the number of files that we have been told to read
    DO i = 1,nfiles
      READ(NAMELIST_UNIT, *, IOSTAT=error) file_names(i), file_time_str
      IF ( error /= 0 )                                                       &
        CALL log_fatal("init_drive",                                          &
                       "Error reading file name/time pair from file " //      &
                       "(IOSTAT=" // TRIM(to_string(error)) // ")")
! Convert the file time as a string into a datetime object
      file_times(i) = datetime_from_string(file_time_str)
    END DO

! Close the file
    CLOSE(NAMELIST_UNIT, IOSTAT=error)
    IF ( error /= 0 )                                                         &
      CALL log_fatal("init_drive",                                            &
                     "Error closing file " // TRIM(file) // " " //            &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")
  ELSE
! If we are not using a list of files, then we detect if we are using time
! templating or not
    use_time_tpl = ( tpl_detect_period(file) < 0 )

    IF ( use_time_tpl ) THEN
      CALL log_info("init_drive",                                             &
                    "Using time templating to get drive file names")
    ELSE
      CALL log_info("init_drive",                                             &
                    "Using single file for all driving data times")
    END IF

! If we are not using templating, then we give lists indicating one file that
! starts at the data start time to register_input_file
    ALLOCATE(file_names(1), STAT=error)
    error_sum = error
    ALLOCATE(file_times(1), STAT=error)
    error_sum = error_sum + error
    IF ( error_sum /= 0 )                                                     &
      CALL log_fatal("init_drive",                                            &
                     "Error allocating arrays for file names and times")

    file_names(1) = file
    file_times(1) = data_start_dt
  END IF

! Check if we are going to be using variable name templating or not
  IF ( use_time_tpl ) THEN
! If using time templating, just check the given file name for variable templating
    use_var_name_tpl = tpl_has_var_name(file)
  ELSE
! If using a list of files, check that they either all have variable templating
! or all don't

! Check the first file first
    use_var_name_tpl = tpl_has_var_name(file_names(1))

    DO i = 2,nfiles
! Check that the rest of the file names match
      IF ( use_var_name_tpl .NEQV. tpl_has_var_name(file_names(i)) )          &
        CALL log_fatal("init_drive",                                          &
                       "If providing a list of files, either they must " //   &
                       "all use variable name templating or all NOT use " //  &
                       "variable name templating")
    END DO
  END IF

  IF ( use_var_name_tpl )                                                     &
    CALL log_info("init_drive",                                               &
                  "Using variable name templating to get drive file names")


!-----------------------------------------------------------------------------
! Build the list of required driving variables
!-----------------------------------------------------------------------------
! First, variables that are always read in directly - they can't be derived
  nvars_required = 3
  required_vars(1:3) = (/ 'pstar', 'q    ', 't    ' /)

!-----------------------------------------------------------------------------
! Radiation
! Currently we always require downward fluxes, but these can be
! derived from net fluxes. This case is rather unusual in that the
! net fluxes are stored in the variables used for downward fluxes, until
! the downward fluxes are calculated in subroutine CONTROL.
! So although we require the downward fluxes, we may not save them
! at this level of code - waits until CONTROL.
! CONTROL uses the value of io_rad_type to decide what to do with the fluxes
! We can detect what value of io_rad_type to use by looking at the given
! variables
!-----------------------------------------------------------------------------
  IF ( ANY(var(1:nvars) == 'rad_net') ) THEN
    CALL log_info("init_drive",                                               &
                  "Downward SW and downward LW radiation will be derived " // &
                  "from net downward all-wavelength and downward SW radiation")

! If net downward all-wavelength radiation is given, we need to use io_rad_type=2
    io_rad_type = 2

! io_rad_type = 2 requires net downward radiation and sw down radiation
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'rad_net'

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'sw_down'

  ELSE IF ( ANY(var(1:nvars) == 'lw_net') ) THEN
    CALL log_info("init_drive",                                               &
                  "Downward SW and downward LW radiation will be derived " // &
                  "from net downward LW and net downward SW radiation")

! If net downward longwave radiation is given, we need to use io_rad_type = 3
    io_rad_type = 3

! io_rad_type = 3 requires net downward lw and net downward sw
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'lw_net'

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'sw_net'

  ELSE
    CALL log_info("init_drive",                                               &
                  "Downward LW and downward SW radiation are both provided directly")

! The default is to use io_rad_type=1 - downward lw and sw are given directly
    io_rad_type = 1

! io_rad_type = 1 requires downward longwave and downward shortwave radiation
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'lw_down'

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'sw_down'
  END IF

! Check if diffuse radiation will be supplied from file
! Otherwise we will use the given constant
  IF ( ANY(var(1:nvars) == 'diff_rad') ) THEN
    CALL log_info("init_drive", "Diffuse radiation will be read from file")

    use_diff_rad = .TRUE.

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'diff_rad'
  ELSE
    CALL log_info("init_drive", "Diffuse radiation will be set as a constant")
  END IF

!-----------------------------------------------------------------------------
! Precipitation
! Currently we always require four components (large scale/convective
! solid/liquid precip), but these can be derived from totals
!-----------------------------------------------------------------------------
  IF ( ANY(var(1:nvars) == 'precip') ) THEN
    CALL log_info("init_drive",                                              &
                  "Precipitation components will be derived from total " //  &
                  "precipitation")

! If total precip is given, we need io_precip_type=1
    io_precip_type = 1

! The only required variable is the total precip
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'precip'

  ELSE IF ( ANY(var(1:nvars) == 'tot_rain') ) THEN
    CALL log_info("init_drive",                                              &
                  "Precipitation components will be derived from total " //  &
                  "rainfall and total snowfall")

! If total rainfall is given, we need io_precip_type=2
    io_precip_type = 2

! We require total rain and total snow to be given
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'tot_rain'

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'tot_snow'

  ELSE IF ( ANY(var(1:nvars) == 'tot_snow') ) THEN
    CALL log_info("init_drive",                                              &
                  "Precipitation components will be derived from " //        &
                  "convective rainfall, large-scale rainfall and total " //  &
                  "snowfall")

! If we are given total snow without total rain (checked above), we need
! io_precip_type = 3
    io_precip_type = 3

! This requires rainfall components to be given separately but total snowfall
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'con_rain'

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'ls_rain'

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'tot_snow'

  ELSE
    CALL log_info("init_drive",                                              &
                  "All precipitation components are provided directly")

! The default is to assume all four components are input directly
! This is io_precip_type=4
    io_precip_type = 4

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'con_rain'

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'ls_rain'

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'con_snow'

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'ls_snow'
  END IF

!-------------------------------------------------------------------------------
! Wind
! Currently we always require two components, but these can be set from the
! total
!-------------------------------------------------------------------------------
  IF ( ANY(var(1:nvars) == 'wind') ) THEN
    CALL log_info("init_drive",                                               &
                  "Horizontal components of wind will be derived from wind speed")

! If total wind is provided, then use it
    io_wind_speed = .TRUE.

! The total wind speed is the only required variable
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'wind'

  ELSE
    CALL log_info("init_drive",                                               &
                  "Horizontal components of wind given directly")

! Otherwise we require both components to be given
    io_wind_speed = .FALSE.

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'u'

    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'v'
  END IF


!-----------------------------------------------------------------------------
! Check deduced values for consistency
!-----------------------------------------------------------------------------
  IF ( .NOT. l_point_data .AND.                                               &
       ( io_precip_type == 1 .OR. io_precip_type == 2) ) THEN
    IF ( t_for_con_rain <= t_for_snow )                                       &
      CALL log_fatal("init_drive",                                            &
                     "t_for_con_rain <= t_for_snow - precip must be " //      &
                     "liquid at t_for_con_rain")
  ENDIF


!-----------------------------------------------------------------------------
! Check that all the required variables have been given
!-----------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                        &
      CALL log_fatal("init_drive",                                            &
                     "Could not find required variable '" //                  &
                     TRIM(required_vars(i)) // "' in list")
  END DO

!-----------------------------------------------------------------------------
! Check which given variables we will be using
!-----------------------------------------------------------------------------
  DO i = 1,nvars
!-----------------------------------------------------------------------------
! If the variable is one of the required vars, then we will be using it
!-----------------------------------------------------------------------------
    IF ( ANY(required_vars(1:nvars_required) == var(i)) ) THEN
      CALL log_info("init_drive",                                             &
                    "'" // TRIM(var(i)) // "' will be read from file")

! If the variable will be filled from file, register it here
      nvars_file = nvars_file + 1
! Since nvars_file <= i (so we will not overwrite unprocessed values)
! and we do not need the values from these arrays for any non-file variables
! from now on, we can just compress them down onto variables that are in the file
      var(nvars_file)      = var(i)
      var_name(nvars_file) = var_name(i)
      tpl_name(nvars_file) = tpl_name(i)
      interp(nvars_file)   = interp(i)

    ELSE
! If the variable is not a required variable, warn about not using it
      CALL log_warn("init_drive",                                             &
                    "Provided variable '" // TRIM(var(i)) //                  &
                    "' is not required, so will be ignored")
    END IF
  END DO

!-----------------------------------------------------------------------------
! Register the input file(s)
!-----------------------------------------------------------------------------
  IF ( use_var_name_tpl ) THEN
! If we are using a variable name template, we must loop through the variables
! and register one file (set of files) per variable
    DO i = 1,nvars_file
      CALL register_input_file(data_start_dt, data_end_dt, data_period,       &
! Driving data should not be a climatology
                               .FALSE.,                                       &
                               use_time_tpl,                                  &
! Substitute the variable name into the file
                               tpl_substitute_var(file, tpl_name(i)),         &
! Build the list of file names using an array comprehension by substituting
! the variable name into each file name
                               (/ (tpl_substitute_var(                        &
                                     file_names(j), tpl_name(i)               &
                               ), j = 1,SIZE(file_names)) /),                 &
                               file_times,    &
! Use array constructors to give the single values related to the variable
                               (/ var(i) /), (/ var_name(i) /),               &
                               (/ interp(i) /))
    END DO

  ELSE
! We are not using variable name templating, so register the same file as
! providing all variables
    CALL register_input_file(data_start_dt, data_end_dt, data_period,         &
! Driving data should not be a climatology
                             .FALSE.,                                         &
                             use_time_tpl, file, file_names, file_times,      &
                             var(1:nvars_file), var_name(1:nvars_file),       &
                             interp(1:nvars_file))
  END IF

  DEALLOCATE(file_names)
  DEALLOCATE(file_times)

  RETURN

END SUBROUTINE init_drive
#endif
