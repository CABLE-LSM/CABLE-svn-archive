#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/init_time.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

SUBROUTINE init_time()

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE datetime_mod, ONLY : DATETIME_STR_LEN, SECS_IN_DAY,                     &
                           datetime,                                          &
! Also import the comparison operators on the datetime type
                           OPERATOR(.eq.), OPERATOR(.ne.), OPERATOR(.lt.),    &
                           OPERATOR(.gt.), OPERATOR(.le.), OPERATOR(.ge.),    &
                           datetime_from_string, datetime_to_string,          &
                           datetime_advance, datetime_subtract

  USE model_interface_mod, ONLY : IDENTIFIER_LEN

  USE model_time_mod, ONLY :                                                  &
    run_min_time, run_max_time, current_time, timestep, timestep_len,         &
    is_spinup, max_spinup_cycles, spinup_cycle, terminate_on_spinup_fail,     &
! Alias some variables to different names
    main_run_start_dt => main_run_start, main_run_end_dt => main_run_end,     &
    spinup_start_dt => spinup_start, spinup_end_dt => spinup_end,             &
    start_of_year, end_of_year

  USE spinup_mod, ONLY : MAX_SPINUP_VARS, nvars, spinup_vars

  USE trifctl, ONLY : phenol_period, triffid_period, asteps_since_triffid

  USE switches, ONLY : l_triffid, l_trif_eq

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about the run times and spinup and checks for
!   consistency
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  TYPE(datetime) :: dt  ! Placeholder for a datetime in a calculation

  INTEGER :: error  ! Error indicator

  INTEGER :: i  ! Loop counter

!-----------------------------------------------------------------------------
! Definition of jules_time namelist - this brings together some variables
! defined in disparate modules, and also reads some local variables that
! are transformed in this routine
!-----------------------------------------------------------------------------
  CHARACTER(len=DATETIME_STR_LEN) :: main_run_start, main_run_end
                                     ! Character strings for times of main run
  NAMELIST /jules_time/ timestep_len, main_run_start, main_run_end,           &
                        phenol_period, triffid_period

!-----------------------------------------------------------------------------
! Definition of jules_spinup namelist - this reads some variables in spinup_mod
! and some local variables that are transformed in this routine
!-----------------------------------------------------------------------------
  CHARACTER(len=DATETIME_STR_LEN) :: spinup_start, spinup_end
                                ! Character strings for times of spinup cycles
  CHARACTER(len=IDENTIFIER_LEN) :: var(MAX_SPINUP_VARS)
                                ! Identifiers of the variables to use for
                                ! spinup
  LOGICAL :: use_percent(MAX_SPINUP_VARS)
                                ! T - use a percentage of the previous value
                                !     as a tolerance
                                ! F - use an absolute tolerance
  REAL :: tolerance(MAX_SPINUP_VARS)
                                ! The tolerance to use
  NAMELIST /jules_spinup/ max_spinup_cycles, spinup_start, spinup_end,        &
                          terminate_on_spinup_fail, nvars, var, use_percent,  &
                          tolerance


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  use_percent(:) = .FALSE.  ! Default is to use absolute tolerance for all
                            ! spinup variables

!-----------------------------------------------------------------------------
! Read namelists
!-----------------------------------------------------------------------------
  OPEN(NAMELIST_UNIT, FILE='timesteps.nml',                                   &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_time",                                               &
                   "Error opening namelist file timesteps.nml " //            &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! First read the time namelist
  CALL log_info("init_time", "Reading JULES_TIME namelist...")
  READ(NAMELIST_UNIT, nml=jules_time, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_time",                                               &
                   "Error reading namelist JULES_TIME " //                    &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Then read the spinup namelist
  CALL log_info("init_time", "Reading JULES_SPINUP namelist...")
  READ(NAMELIST_UNIT, nml=jules_spinup, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_time",                                               &
                   "Error reading namelist JULES_SPINUP " //                  &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_time",                                               &
                   "Error closing namelist file timesteps.nml " //            &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


!-----------------------------------------------------------------------------
! Set values derived from namelists and verify for consistency
!-----------------------------------------------------------------------------
! For simplicity, we enforce that timestep is a divisor of one day
  IF ( MOD(SECS_IN_DAY, timestep_len) /= 0 )                                  &
    CALL log_fatal("init_time",                                               &
                   "A day must contain a whole number of timesteps")
  CALL log_info("init_time",                                                  &
                "Timestep is " // TRIM(to_string(timestep_len)) // " seconds")

! Check that TRIFFID timestep (the coupling period) seems sensible
! In equilibrium mode, the coupling period should be sufficient to average
! out seasonal (and ideally interannual) variability, so a period >= 1yr
! and perhaps 5-10 yrs is recommended
! In transient mode, more frequent coupling is required - typically every
! 10 days or so
  IF ( l_trif_eq .AND. triffid_period < 360 ) THEN
    CALL log_error("init_time",                                               &
                   "triffid_period < 360 - in equilibrium mode a " //         &
                   "TRIFFID timestep of at least 1 year is advised")

!   Note: We would expect "normal usage" to be that the coupling period
!   is an integer number of years (or at least close to this if the
!   calendar includes leap years). However, there's generally nothing
!   intrinsically wrong with using partial years (as long as the period
!   is >> 1 year) - it's just a slightly odd choice. In any case, we are
!   not testing if the period is a number of years.

  ELSE IF ( l_triffid .AND. triffid_period > 30 ) THEN
    CALL log_error("init_time",                                               &
                   "triffid_period > 30 - in dynamic mode a TRIFFID " //      &
                   "timestep of <30 days is recommended; 10 days is " //      &
                   "often used")
  END IF

! Convert character strings for main run times into datetimes
  main_run_start_dt = datetime_from_string(main_run_start)
  main_run_end_dt = datetime_from_string(main_run_end)

! Check that the main run times are in chronological order
  IF ( main_run_end_dt <= main_run_start_dt )                                 &
    CALL log_fatal("init_time",                                               &
                   "Start time for main run must be before end time")

  CALL log_info("init_time",                                                  &
                "Main run start - " // datetime_to_string(main_run_start_dt))
  CALL log_info("init_time",                                                  &
                "Main run end - " // datetime_to_string(main_run_end_dt))

! Assume that there is no spinup for now
! I.e. the main run is the entire run, and the time for the first timestep
! is the start of the run
  run_min_time = main_run_start_dt
  run_max_time = main_run_end_dt
  current_time = main_run_start_dt

! Check for spinup
  is_spinup = max_spinup_cycles > 0

  IF ( is_spinup .AND. nvars < 1 ) THEN
! If no spinup variables were specified, that is the same as specifying no
! spinup
    CALL log_warn("init_time",                                                &
                  "Spinup requested (max_spinup_cycles > 0) but no " //       &
                  "variables specified - spinup will not be performed")

    is_spinup         = .FALSE.
    max_spinup_cycles = 0
  END IF

  IF ( is_spinup ) THEN
    CALL log_info("init_time",                                                &
                  "Spinup requested - maximum number of spinup cycles is " // &
                  TRIM(to_string(max_spinup_cycles)))

! Convert character strings for spinup times into datetimes
    spinup_start_dt = datetime_from_string(spinup_start)
    spinup_end_dt = datetime_from_string(spinup_end)

! Check that the spinup times are in chronological order
    IF ( spinup_end_dt <= spinup_start_dt )                                   &
      CALL log_fatal("init_time",                                             &
                     "Start time for spinup must be before end time")

    CALL log_info("init_time",                                                &
                  "Spinup start - " // datetime_to_string(spinup_start_dt))
    CALL log_info("init_time",                                                &
                  "Spinup end - " // datetime_to_string(spinup_end_dt))

! Adjust the run times as required
    IF ( spinup_start_dt < main_run_start_dt ) run_min_time = spinup_start_dt
    IF ( spinup_end_dt > main_run_end_dt ) run_max_time = spinup_end_dt
! The time for the first timestep is now the start of spinup
    current_time = spinup_start_dt

! We start in the first cycle of spinup
    spinup_cycle = 1

! Set up the spinup variables
    IF ( nvars > MAX_SPINUP_VARS )                                            &
      CALL log_fatal("init_time", "Too many spinup variables specified")

    DO i = 1,nvars
      spinup_vars(i)%identifier  = var(i)
      spinup_vars(i)%use_percent = use_percent(i)
      spinup_vars(i)%tolerance   = tolerance(i)
    END DO
  ELSE
    CALL log_info("init_time", "No spinup requested")
  END IF

! We are (obviously) starting in the first timestep
  timestep = 1
! Initialise the number of timesteps since triffid was called
  asteps_since_triffid = 1

! Set the variables that indicate if we are starting or ending a year
! If going back one timestep length causes the year to change, then we are in
! the first timestep of a year
  dt = datetime_subtract(current_time, timestep_len)
  start_of_year = dt%year /= current_time%year
! If advancing the current time causes the year to change, then it is the last
! timestep in a year
  dt = datetime_advance(current_time, timestep_len)
  end_of_year = dt%year /= current_time%year

  RETURN

END SUBROUTINE init_time
#endif
