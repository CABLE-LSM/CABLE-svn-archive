#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ts_open(mode, data_start, data_end, data_period, is_climatology,&
                      use_time_template, template, files, file_times,         &
                      comm, info)                                             &
                              RESULT(file)

  USE io_constants, ONLY : MODE_READ, MODE_WRITE

  USE datetime_mod, ONLY : PERIOD_MONTH, PERIOD_YEAR, datetime_advance

  USE templating_mod, ONLY : tpl_detect_period

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a timeseries file and returns a file_ts object representing it
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  INTEGER, INTENT(IN) :: mode ! The mode to open the file
                              ! One of MODE_READ or MODE_WRITE

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
  CHARACTER(len=*), OPTIONAL, INTENT(IN) :: template
                              ! The time template to use

! With a file list
  CHARACTER(len=*), OPTIONAL, INTENT(IN) :: files(:)
                              ! List of file names
  TYPE(datetime), OPTIONAL, INTENT(IN) :: file_times(:)
                              ! Time of first data for each file

! MPI variables
  INTEGER, OPTIONAL, INTENT(IN) :: comm  ! MPI communicator to use for parallel
                                         ! I/O
                                         ! If not given, serial I/O is used
  INTEGER, OPTIONAL, INTENT(IN) :: info  ! MPI info object to use for parallel
                                         ! I/O
                                         ! If not given, serial I/O is used

! Return type
  TYPE(file_ts) :: file


! Work variables
  INTEGER :: nfiles ! The number of files in the files/file_times lists
  INTEGER :: tpl_period ! The detected templating period

  TYPE(datetime) :: prev_file_start  ! The start time of the previous file
                                     ! in the list
                                     ! Used when checking that files
                                     ! given in a list are in chronological
                                     ! order
  LOGICAL :: file_exists  ! Indicates if the file exists
                          ! Used when checking that files given in a list
                          ! exist

  INTEGER :: i ! Loop counter


!-----------------------------------------------------------------------------


! Initialise nfiles - this is not used when time-templating is enabled
  nfiles = -1

!*****************************************************************************
! Check that a valid combination of options has been supplied
!*****************************************************************************
!-----------------------------------------------------------------------------
! Check that the start and end times for the data make sense
!-----------------------------------------------------------------------------
! Obviously, the data must start before they end, so bail if not true
  IF ( data_start >= data_end )                                               &
    CALL log_fatal("file_ts_open",                                            &
                   "Data start must be strictly before data end.")

! Check that the data period is allowed
  SELECT CASE ( data_period )
    CASE ( 1:, PERIOD_MONTH, PERIOD_YEAR )
! These are the allowed data periods
      CALL log_info("file_ts_open",                                           &
                    "Opening time series with data_period=" //                &
                    TRIM(to_string(data_period)))

    CASE DEFAULT
      CALL log_fatal("file_ts_open",                                          &
                     "Data period must be > 0 or a 'special' period " //      &
                     "(supplied " // TRIM(to_string(data_period)) // ")")
  END SELECT

! Check that the start time is appropriate if using a special period
  IF ( data_period == PERIOD_MONTH .OR. data_period == PERIOD_YEAR ) THEN
    IF ( data_start%time /= 0 .OR. data_start%day /= 1 )                      &
      CALL log_fatal("file_ts_open",                                          &
                     "When using data_period=" //                             &
                     TRIM(to_string(data_period)) //                          &
                     ", data must start at 00:00:00 on 1st of month")

    IF ( data_period == PERIOD_YEAR .AND. data_start%month /= 1 )             &
      CALL log_fatal("file_ts_open",                                          &
                     "When using data_period=" //                             &
                     TRIM(to_string(data_period)) //                          &
                     ", data must start at 00:00:00 on 1st of January")
  END IF

!-----------------------------------------------------------------------------
! Check that we have appropriate arguments to be able to locate the
! correct files to use
!-----------------------------------------------------------------------------
  IF ( use_time_template ) THEN
! If using time templating, we need a template to use
    IF ( .NOT. PRESENT(template) )                                            &
      CALL log_fatal("file_ts_open",                                          &
                     "Time templating selected but no template given")

! Check that template is not an empty string
    IF ( LEN_TRIM(template) <= 0 )                                            &
      CALL log_fatal("file_ts_open",                                          &
                     "Time templating selected but template is empty string")

    tpl_period = tpl_detect_period(template)

! If the period is not one that has code written for it, abort
    SELECT CASE ( tpl_period )
      CASE ( PERIOD_MONTH, PERIOD_YEAR )
! Fine - just log some information
        CALL log_info("file_ts_open",                                         &
                      "Detected period=" // TRIM(to_string(tpl_period)) //    &
                      " for template " // TRIM(template))

      CASE DEFAULT
        CALL log_fatal("file_ts_open",                                        &
                       "Could not detect supported templating period " //     &
                       "for template " // TRIM(template))
    END SELECT
  ELSE
! If using lists of files and start times, check that they are present and
! consistent
    IF ( .NOT. PRESENT(files) )                                               &
      CALL log_fatal("file_ts_open",                                          &
                     "Time templating is not selected - a list of files " //  &
                     "must be given")

    IF ( .NOT. PRESENT(file_times) )                                          &
      CALL log_fatal("file_ts_open",                                          &
                     "Time templating is not selected - a list of file " //   &
                     "times must be given")

    nfiles = SIZE(files)

    IF ( nfiles <= 0 )                                                        &
      CALL log_fatal("file_ts_open",                                          &
                     "Time templating is not selected - list must contain" // &
                     " at least one file")

    IF ( SIZE(file_times) /= nfiles )                                         &
      CALL log_fatal("file_ts_open",                                          &
                     "'files' and 'file_times' must have the same " //        &
                     "number of entries")

! If using yearly data, there must be one file containing all the data
! We know from previous checks on file_times(1) and data_start that this file
! starts on the 1st Jan for some year
    IF ( data_period == PERIOD_YEAR .AND. nfiles > 1 )                        &
      CALL log_fatal("file_ts_open",                                          &
                     "Yearly data must be contained in a single file for " // &
                     "the entirety of the data")

! Do checks on the individual files
    DO i = 1,nfiles
      IF ( mode == MODE_READ ) THEN
! By this point, we know that we don't have variable name or time templating,
! so we should check if the given files exist
        INQUIRE(FILE=files(i), EXIST=file_exists)
        IF ( .NOT. file_exists )                                              &
          CALL log_fatal("file_ts_open",                                      &
                         "Given file '" // TRIM(files(i)) // "' does not exist")
      ELSE
! If we are in write mode, just check that a non-empty file name has been
! given
        IF ( LEN_TRIM(files(i)) <= 0 )                                        &
          CALL log_fatal("file_ts_open",                                      &
                         "List of files given, but one or more file names" // &
                         " are the empty string")
      END IF

! Check that the start times of the files are suitable
      IF ( i == 1 ) THEN
! Check that the first file starts at data start
        IF ( file_times(i) /= data_start )                                    &
          CALL log_fatal("file_ts_open",                                      &
                         "Start time for the first file must match data_start")

        prev_file_start = file_times(i)
      ELSE
! Check that all subsequent files start after the previous file
        IF ( file_times(i) <= prev_file_start )                               &
          CALL log_fatal("file_ts_open",                                      &
                         "Files must be given in chronological order")

        prev_file_start = file_times(i)
      END IF

! If using monthly data, files must start at midnight on the 1st of some
! month (but not necessarily only contain data for one month)
      IF ( data_period == PERIOD_MONTH ) THEN
        IF ( file_times(i)%time /= 0 .OR. file_times(i)%day /= 1 )            &
          CALL log_fatal("file_ts_open",                                      &
                         "When using monthly data with a list of files, " //  &
                         "all files must start at 00:00:00 on 1st of " //     &
                         "some month")
      END IF
    END DO
  END IF

!-----------------------------------------------------------------------------
! If we have a climatology, check that the given arguments and detected
! properties make sense
!-----------------------------------------------------------------------------
  IF ( is_climatology ) THEN
! It makes no sense to specify is_climatology if in write mode
    IF ( mode == MODE_WRITE )                                                 &
      CALL log_fatal("file_ts_open",                                          &
                     "Cannot open a file in write mode as a climatology.")

! The data must start at the beginning of the year
    IF ( data_start%time /= 0 .OR.                                            &
         data_start%day /= 1 .OR. data_start%month /= 1 )                     &
      CALL log_fatal("file_ts_open",                                          &
                     "When using data as a climatology, data must start " //  &
                     "at 00:00:00 on 1st of January")

! The data must apply for exactly one year
    IF ( datetime_advance(data_start, PERIOD_YEAR) /= data_end )              &
      CALL log_fatal("file_ts_open",                                          &
                     "When using data as a climatology, exactly one year " // &
                     "of data must be given.")

! Warn the user if they have supplied yearly data - this means the same data
! will be used every year, and is the same as having fixed data
    IF ( data_period == PERIOD_YEAR )                                         &
      CALL log_warn("file_ts_open",                                           &
                    "When using a climatology, a data period of a year " //   &
                    "is equivalent to having fixed data")
  END IF

!-----------------------------------------------------------------------------
! Check that the MPI variables are either specified together or not at all
!-----------------------------------------------------------------------------
  IF ( PRESENT(comm) .NEQV. PRESENT(info) )                                   &
    CALL log_fatal("file_ts_open",                                            &
                   "Only one of comm and info is present - either give a " // &
                   "value for both MPI variables for parallel access or " //  &
                   "omit both for serial access")


!-----------------------------------------------------------------------------
! Now we are happy that the arguments are consistent, we can set up the
! file_ts object
!-----------------------------------------------------------------------------
  file%mode = mode

! Store info about MPI-IO settings
  IF ( PRESENT(comm) ) THEN
    file%use_mpiio = .TRUE.
    file%comm = comm
    file%info = info
  END IF

  file%data_start  = data_start
  file%data_end    = data_end
  file%data_period = data_period
  file%is_climatology = is_climatology

  file%use_time_template = use_time_template

  IF ( use_time_template ) THEN
! With time templating
    file%template   = template
    file%tpl_period = tpl_period
  ELSE
! With a file list
    file%nfiles = nfiles

! Allocate space for the files and file_times arrays
    ALLOCATE(file%files(nfiles))
    ALLOCATE(file%file_times(nfiles))

    file%files(:)      = files(:)
    file%file_times(:) = file_times(:)
  END IF

  RETURN

END FUNCTION file_ts_open
#endif
