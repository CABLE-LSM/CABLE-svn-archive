#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/init_prescribed_data.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

SUBROUTINE init_prescribed_data()

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN, MAX_VAR_FILE, &
                           NAMELIST_UNIT, FILE_LIST_UNIT

  USE model_interface_mod, ONLY : IDENTIFIER_LEN

  USE datetime_mod, ONLY : DATETIME_STR_LEN, datetime, datetime_from_string

  USE string_utils_mod, ONLY : to_string

  USE templating_mod, ONLY : tpl_detect_period, tpl_has_var_name,             &
                             tpl_substitute_var

  USE time_varying_input_mod, ONLY : register_input_file

  USE update_mod, ONLY : have_prescribed_veg

  USE switches, ONLY : l_phenol, l_triffid, l_o3_damage

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about what variables are prescribed on what timestep
!   from what files and initialises the files
!   We allow the user to specify any variables they like, so long as the code
!   exists in model_interface_mod to populate them from file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  TYPE(datetime) :: data_start_dt, data_end_dt
                         ! Datetime objects created from strings given in
                         ! output profiles

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

  LOGICAL :: have_ozone  ! Used to check for the presence of ozone
                         ! as a prescribed variable, since this
                         ! is required if ozone damage is switched
                         ! on

  LOGICAL :: have_lai    ! Used to check for the presence of veg
  LOGICAL :: have_canht  ! variables as prescribed variables, since
                         ! this is only allowed in certain model
                         ! configurations and requires special
                         ! treatment in update_mod

  INTEGER :: i,j,k  ! Index variables

  INTEGER :: error, error_sum  ! Error indicators

!-----------------------------------------------------------------------------
! Definition of the jules_prescribed namelist
!-----------------------------------------------------------------------------
  INTEGER :: n_datasets      ! The number of datasets containing prescribed
                             ! variables
  NAMELIST /jules_prescribed/ n_datasets

!-----------------------------------------------------------------------------
! Definition of the jules_prescribed_dataset namelist
!-----------------------------------------------------------------------------
! Information about the data in the file
  CHARACTER(len=DATETIME_STR_LEN) :: data_start, data_end
                             ! Start and end times for data as strings
  INTEGER :: data_period     ! The period of the driving data
  LOGICAL :: is_climatology  ! T - use the file to provide a climatology
                             !     for the specified variables
                             ! F - the file is not used as a climatology

! Information about the files to use
  LOGICAL :: read_list       ! T - the given file contains a list of file
                             !     names and times of first data
                             !     These files may contain variable name
                             !     templating but not time templating
                             ! F - data should be read directly from
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
  CHARACTER(len=IDENTIFIER_LEN) :: var(MAX_VAR_FILE)
                        ! The variable identifiers of the variables
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_name(MAX_VAR_FILE)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: tpl_name(MAX_VAR_FILE)
                        ! The name to substitute in a template for each
                        ! variable
  CHARACTER(len=2) :: interp(MAX_VAR_FILE)
                        ! Flag indicating the type of interpolation to use
                        ! for each variable


  NAMELIST /jules_prescribed_dataset/ data_start, data_end, data_period,      &
                                      is_climatology, read_list, nfiles, file,&
                                      nvars, var, var_name, tpl_name, interp


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  have_ozone = .FALSE.
  have_lai   = .FALSE.
  have_canht = .FALSE.
  n_datasets = 0

!-----------------------------------------------------------------------------
! Read namelists
!-----------------------------------------------------------------------------
  OPEN(NAMELIST_UNIT, FILE='prescribed_data.nml',                             &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_prescribed_data",                                    &
                   "Error opening namelist file prescribed_data.nml " //      &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Read the JULES_PRESCRIBED namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_prescribed_data", "Reading JULES_PRESCRIBED namelist...")
  READ(NAMELIST_UNIT, nml=jules_prescribed, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_prescribed_data",                                    &
                   "Error reading namelist JULES_PRESCRIBED " //              &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Process the information from the JULES_PRESCRIBED namelist
!-----------------------------------------------------------------------------
  IF ( n_datasets < 1 ) THEN
    CALL log_info("init_prescribed_data",                                     &
                  "No data (other than driving data) is prescribed for " //   &
                  "this run")
  END IF

!-----------------------------------------------------------------------------
! Read and process information about each input file in turn
!-----------------------------------------------------------------------------
  DO i = 1,n_datasets
!-----------------------------------------------------------------------------
! Set namelist values to their defaults before reading about the next file
!-----------------------------------------------------------------------------
    data_start     = ''
    data_end       = ''
    data_period    = 0
    is_climatology = .FALSE.
    read_list      = .FALSE.
    nfiles         = 0
    file           = ''
    nvars          = 0

!-----------------------------------------------------------------------------
! Read the namelist
!-----------------------------------------------------------------------------
    CALL log_info("init_prescribed_data",                                     &
                  "Reading JULES_PRESCRIBED_DATASET namelist...")
    READ(NAMELIST_UNIT, nml=jules_prescribed_dataset, IOSTAT=error)
    IF ( error /= 0 )                                                         &
      CALL log_fatal("init_prescribed_data",                                  &
                     "Error reading namelist JULES_PRESCRIBED_DATASET " //    &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! If the file is not providing any variables, we can skip it
!-----------------------------------------------------------------------------
    IF ( nvars < 1 ) THEN
      CALL log_error("init_prescribed_data",                                  &
                     "File " // TRIM(file) // " is not providing any " //     &
                     "variables - ignoring")
      CYCLE
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
      IF ( nfiles < 1 )                                                       &
        CALL log_fatal("init_prescribed_data",                                &
                       "If reading a list of file names and file times, " //  &
                       "at least one file must be given")

! If we are reading a list of files, then we will definitely not be using
! time templating
      use_time_tpl = .FALSE.

!-----------------------------------------------------------------------------
! Read the list of file names and times from the file
!-----------------------------------------------------------------------------
      CALL log_info("init_prescribed_data",                                   &
                    "Reading list of data file names and start times " //     &
                    "from " // TRIM(file) // "...")

! First allocate space for the elements
      ALLOCATE(file_names(nfiles), STAT=error)
      error_sum = error
      ALLOCATE(file_times(nfiles), STAT=error)
      error_sum = error_sum + error
      IF ( error_sum /= 0 )                                                   &
        CALL log_fatal("init_prescribed_data",                                &
                       "Error allocating arrays for file names and times")

! Open the file
      OPEN(FILE_LIST_UNIT, FILE=file, STATUS='old', POSITION='rewind',        &
                                                 ACTION='read', IOSTAT=error)
      IF ( error /= 0 )                                                       &
        CALL log_fatal("init_prescribed_data",                                &
                       "Error opening file " // TRIM(file) // " " //          &
                       "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Read the number of files that we have been told to read
      DO j = 1,nfiles
        READ(FILE_LIST_UNIT, *, IOSTAT=error) file_names(j), file_time_str
        IF ( error /= 0 )                                                     &
          CALL log_fatal("init_prescribed_data",                              &
                         "Error reading file name/time pair from file " //    &
                         "(IOSTAT=" // TRIM(to_string(error)) // ")")
! Convert the file time as a string into a datetime object
        file_times(j) = datetime_from_string(file_time_str)
      END DO

! Close the file
      CLOSE(FILE_LIST_UNIT, IOSTAT=error)
      IF ( error /= 0 )                                                       &
        CALL log_fatal("init_prescribed_data",                                &
                       "Error closing file " // TRIM(file) // " " //          &
                       "(IOSTAT=" // TRIM(to_string(error)) // ")")
    ELSE
!-----------------------------------------------------------------------------
! If we are not using a list of files, then we detect if we are using time
! templating or not
!-----------------------------------------------------------------------------
      use_time_tpl = ( tpl_detect_period(file) < 0 )

      IF ( use_time_tpl ) THEN
        CALL log_info("init_prescribed_data",                                 &
                      "Using time templating to get data file names")
      ELSE
        CALL log_info("init_prescribed_data",                                 &
                      "Using single file for all data times")
      END IF

!-----------------------------------------------------------------------------
! If we are not using templating, then we give lists indicating one file that
! starts at the data start time to register_input_file
!-----------------------------------------------------------------------------
      ALLOCATE(file_names(1), STAT=error)
      error_sum = error
      ALLOCATE(file_times(1), STAT=error)
      error_sum = error_sum + error
      IF ( error_sum /= 0 )                                                   &
        CALL log_fatal("init_prescribed_data",                                &
                       "Error allocating arrays for file names and times")

      file_names(1) = file
      file_times(1) = data_start_dt
    END IF

!-----------------------------------------------------------------------------
! Check if we are going to be using variable name templating or not
!-----------------------------------------------------------------------------
    IF ( use_time_tpl ) THEN
! If using time templating, just check the given file name for variable templating
      use_var_name_tpl = tpl_has_var_name(file)
    ELSE
! If using a list of files, check that they either all have variable templating
! or all don't

! Check the first file first
      use_var_name_tpl = tpl_has_var_name(file_names(1))

      DO j = 2,nfiles
! Check that the rest of the file names match
        IF ( use_var_name_tpl .NEQV. tpl_has_var_name(file_names(j)) )        &
          CALL log_fatal("init_prescribed_data",                              &
                         "If providing a list of files, either they must " // &
                         "all use variable name templating or all NOT use " //&
                         "variable name templating")
      END DO
    END IF

    IF ( use_var_name_tpl )                                                   &
      CALL log_info("init_prescribed_data",                                   &
                    "Using variable name templating to get data file names")

!-----------------------------------------------------------------------------
! Register the input file(s)
!-----------------------------------------------------------------------------
    IF ( use_var_name_tpl ) THEN
! If we are using a variable name template, we must loop through the variables
! and register one file (set of files) per variable
      DO j = 1,nvars
        CALL register_input_file(data_start_dt, data_end_dt, data_period,     &
                                 is_climatology, use_time_tpl,                &
! Substitute the variable name into the file
                                 tpl_substitute_var(file, tpl_name(j)),       &
! Build the list of file names using an array comprehension by substituting
! the variable name into each file name
                                 (/ (tpl_substitute_var(                      &
                                       file_names(k), tpl_name(j)             &
                                 ), k = 1,SIZE(file_names)) /),               &
                                 file_times,    &
! Use array constructors to give the single values related to the variable
                                 (/ var(j) /), (/ var_name(j) /),             &
                                 (/ interp(j) /))
      END DO

    ELSE
! We are not using variable name templating, so register the same file as
! providing all variables
      CALL register_input_file(data_start_dt, data_end_dt, data_period,       &
                               is_climatology, use_time_tpl,                  &
                               file, file_names, file_times,                  &
                               var(1:nvars), var_name(1:nvars),               &
                               interp(1:nvars))
    END IF

    DEALLOCATE(file_names)
    DEALLOCATE(file_times)

!-----------------------------------------------------------------------------
! Check if any of the specified variables require special treatment, either
! here or in update_derived_variables (once the prescribed variables
! themselves have been updated from file)
!-----------------------------------------------------------------------------
! We want to know if we have ozone, since we require it to be supplied if
! ozone damage is on
    have_ozone = have_ozone .OR. ANY(var(1:nvars) == 'ozone')
! We also want to know when we have veg variables
    have_lai   = have_lai   .OR. ANY(var(1:nvars) == 'lai')
    have_canht = have_canht .OR. ANY(var(1:nvars) == 'canht')

  END DO  ! datasets

! Can't have ozone damage without prescribing ozone
  IF ( l_o3_damage .AND. .NOT. have_ozone )                                   &
    CALL log_fatal("init_prescribed_data",                                    &
                   "When ozone damage is on, ozone must be prescribed")

! Can't prescribe lai if phenology is on, since it is prognostic
  IF ( l_phenol .AND. have_lai )                                              &
    CALL log_fatal("init_prescribed_data",                                    &
                   "When phenology is on, lai is prognostic and cannot " //   &
                   "be prescribed")

! Can't prescribe canopy height if TRIFFID is on, since it is prognostic
  IF ( l_triffid .AND. have_canht )                                           &
    CALL log_fatal("init_prescribed_data",                                    &
                   "When TRIFFID is on, canopy height is prognostic and " //  &
                   "cannot be prescribed")

! update_mod needs to know if we have a veg variable
  have_prescribed_veg = have_lai .OR. have_canht


  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_prescribed_data",                                    &
                   "Error closing namelist file prescribed_data.nml " //      &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  RETURN

END SUBROUTINE init_prescribed_data
#endif
