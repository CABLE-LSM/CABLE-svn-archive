#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/output/register_output_profile.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE register_output_profile(profile_name,                              &
                                   output_spinup, output_main_run,            &
                                   output_start, output_end,                  &
                                   output_period, sample_period, file_period, &
                                   identifiers, var_names, output_types)

  USE datetime_mod, ONLY : PERIOD_MONTH, PERIOD_YEAR, datetime_to_string

  USE model_time_mod, ONLY : max_spinup_cycles, main_run_start, main_run_end, &
                             timestep_len

  USE model_interface_mod, ONLY : get_var_id

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Register an output profile to output given variables over the given
!   period
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  CHARACTER(len=*), INTENT(IN) :: profile_name
      ! The name of this profile - only used in filenames

  LOGICAL, INTENT(IN) :: output_spinup
      ! T - profile will provide output during spinup
      ! F - profile will not provide output during spinup

  LOGICAL, INTENT(IN) :: output_main_run
      ! T - profile will provide output for the specified portion of the
      !     main run only
      ! F - profile will not provide any output during the main run

  TYPE(datetime), INTENT(IN) :: output_start, output_end
      ! USED ONLY IF output_main_run=T
      ! Start and end dates for output if providing output for main run

  INTEGER, INTENT(IN) ::                                                      &
    output_period,                                                            &
      ! The output period for this profile - this should be a multiple of
      ! the model timestep or a 'special' period (monthly/yearly)
    sample_period,                                                            &
      ! The sampling period for this profile - this should be a multiple of
      ! the model timestep or a special period
    file_period
      ! This should be a special period for monthly or yearly files
      ! Any other value will result in a single file for the whole of the
      ! output period

  CHARACTER(len=*), INTENT(IN) :: identifiers(:)
      ! The model identifiers of the variables this profile provides output for

  CHARACTER(len=*), INTENT(IN) :: var_names(:)
      ! The name to use in output files for each variable
      ! Each value can be the empty string - if it is, the model identifier will
      ! be used

  CHARACTER(len=*), INTENT(IN) :: output_types(:)
      ! The type of output to use for each variable


! Work variables
  INTEGER :: nvars  ! The number of output fields in the profile

! Local versions of INTENT(IN) arguments that we can modify
  LOGICAL :: output_spinup_local
  TYPE(datetime) :: output_start_local, output_end_local
  CHARACTER(len=LEN(var_names)) :: var_names_local(SIZE(var_names))

  INTEGER :: i  ! Loop counters


!-----------------------------------------------------------------------------

! Check that we have space to register another profile
  IF ( nprofiles >= NPROFILES_MAX )                                           &
    CALL log_fatal("register_output_profile",                                 &
                   "Too many profiles registered - try increasing NPROFILES_MAX")

! Check that we have an output type for each variable
  nvars = SIZE(identifiers)

  IF ( nvars /= SIZE(var_names) )                                             &
    CALL log_fatal("register_output_profile",                                 &
                   "identifiers and var_names must have the same " //         &
                   "number of elements")

  IF ( nvars /= SIZE(output_types) )                                          &
    CALL log_fatal("register_output_profile",                                 &
                   "identifiers and output_types must have the same " //      &
                   "number of elements")

!-----------------------------------------------------------------------------
! Check that arguments make sense
!-----------------------------------------------------------------------------
! First copy INTENT(IN) arguments that we might want to modify into local
! versions
  output_spinup_local = output_spinup
  output_start_local  = output_start
  output_end_local    = output_end

! Check that we have been given a profile name
  IF ( LEN_TRIM(profile_name) == 0 )                                          &
    CALL log_fatal("register_output_profile", "No profile name given")

  IF ( output_spinup_local .AND. max_spinup_cycles <= 0 ) THEN
! If spinup output has been requested but there is no spinup, issue a warning
! and ignore it
    CALL log_error("register_output_profile",                                 &
                   "Model has no spinup - ignoring request for output " //    &
                   "during spinup from profile " // TRIM(profile_name))
    output_spinup_local = .FALSE.
  END IF

! Check if the above means that this profile has requested no output at all
  IF ( .NOT. output_spinup_local .AND. .NOT. output_main_run ) THEN
    CALL log_error("register_output_profile",                                 &
                   "Profile " // TRIM(profile_name) // " will provide no " // &
                   "output with current model setup - ignoring")
    RETURN
  END IF

  IF ( output_main_run ) THEN
! If output has been requested for times that are outside the main run, we
! issue a warning and truncate
    IF ( output_start_local < main_run_start ) THEN
      CALL log_error("register_output_profile",                               &
                     "Output has been requested for times before the " //     &
                     "start of the main run - output will start at the " //   &
                     "start of the main run")
      output_start_local = main_run_start
    END IF

    IF ( output_end_local > main_run_end ) THEN
      CALL log_error("register_output_profile",                               &
                     "Output has been requested for times after the " //      &
                     "end of the main run - output will end at the " //       &
                     "end of the main run")
      output_end_local = main_run_end
    END IF

! Check that the times make sense
    IF ( output_end_local <= output_start_local )                             &
      CALL log_fatal("register_output_profile",                               &
                     "Output cannot end before it has started")
  END IF

! Check that the given periods make sense
  IF ( file_period /= PERIOD_MONTH .AND. file_period /= PERIOD_YEAR )         &
! Warn that one file will be used for all output, since file_period is not
! a special period
    CALL log_info("register_output_profile",                                  &
                  "Since file_period is not a 'special' period, all " //      &
                  "output from each section will go into one file")

  IF ( output_period /= PERIOD_MONTH .AND. output_period /= PERIOD_YEAR .AND. &
       MOD(output_period, timestep_len) /= 0 )                                &
! If output period is not a special period, then it must be a multiple of
! the model timestep
    CALL log_fatal("register_output_profile",                                 &
                   "Output period must be a 'special' period or a " //        &
                   "multiple of model timestep")

  IF ( sample_period < 1 .OR. MOD(sample_period, timestep_len) /= 0 )         &
    CALL log_fatal("register_output_profile",                                 &
                   "Sample period must be a multiple of model timestep")

!-----------------------------------------------------------------------------
! Indicate what output we are providing
!-----------------------------------------------------------------------------
  IF ( output_spinup_local )                                                  &
    CALL log_info("register_output_profile",                                  &
                  "Profile with name " // TRIM(profile_name) // " " //        &
                  "registered to provide output during spinup")

  IF ( output_main_run )                                                      &
    CALL log_info("register_output_profile",                                  &
                  "Profile with name " // TRIM(profile_name) // " " //        &
                  "registered to provide output for main run from " //        &
                  datetime_to_string(output_start_local) // " to " //         &
                  datetime_to_string(output_end_local))

!-----------------------------------------------------------------------------
! Set up the variable names that will be used for output
! Use user-supplied output name if given, otherwise use the model identifier
!
! We also check for duplicates in the output variable names as we go
! Note that duplicate identifiers are allowed, as long different names have
! been specified for them to use in output files
!-----------------------------------------------------------------------------
  var_names_local(:) = var_names(:)
  DO i = 1,nvars
    IF ( LEN_TRIM(var_names(i)) <= 0 ) var_names_local(i) = identifiers(i)

! If the output name being processed matches any of those previously processed,
! issue a fatal error
! We issue a fatal error to force the user to resolve the problem and avoid any
! confusion over what has been output
    IF ( ANY(var_names_local(1:i-1) == var_names_local(i)) )                  &
      CALL log_fatal("register_output_profile",                               &
                     "Duplicate variable name for output files " //           &
                     "detected for profile " // TRIM(profile_name))
  END DO

!-----------------------------------------------------------------------------
! Set up the output_profile object and its corresponding output_field objects
!-----------------------------------------------------------------------------
! Store profile constant data
  nprofiles = nprofiles + 1
  profiles(nprofiles)%profile_name    = profile_name
  profiles(nprofiles)%output_spinup   = output_spinup_local
  profiles(nprofiles)%output_main_run = output_main_run
  profiles(nprofiles)%file_period     = file_period
  profiles(nprofiles)%output_period   = output_period
  profiles(nprofiles)%sample_period   = sample_period
  profiles(nprofiles)%output_start    = output_start_local
  profiles(nprofiles)%output_end      = output_end_local

! Allocate space on the output_profile object for the fields
  profiles(nprofiles)%nfields = nvars
  ALLOCATE(profiles(nprofiles)%fields(nvars))

  DO i = 1,nvars
! Store as much info about fields as we currently know for later
    profiles(nprofiles)%fields(i)%var_id      = get_var_id(identifiers(i))
    profiles(nprofiles)%fields(i)%output_name = var_names_local(i)
    profiles(nprofiles)%fields(i)%type        = output_types(i)
  END DO

  RETURN

END SUBROUTINE register_output_profile
#endif
