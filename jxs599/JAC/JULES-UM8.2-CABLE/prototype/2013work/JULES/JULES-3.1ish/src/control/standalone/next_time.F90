#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/control/standalone/next_time.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

SUBROUTINE next_time()

  USE string_utils_mod, ONLY : to_string

  USE datetime_mod, ONLY : datetime,                                          &
! Also import the comparison operators on the datetime type
                           OPERATOR(.eq.), OPERATOR(.ne.), OPERATOR(.lt.),    &
                           OPERATOR(.gt.), OPERATOR(.le.), OPERATOR(.ge.),    &
                           datetime_advance, datetime_subtract,               &
                           datetime_to_string

  USE time_varying_input_mod, ONLY : seek_all_to_current_datetime,            &
                                     advance_all

  USE spinup_mod, ONLY : spinup_check

  USE model_time_mod

  USE dump_mod, ONLY : write_dump

  USE trifctl, ONLY : asteps_since_triffid

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Advances the model to the next timestep. Returns logical indicating end
!   of run in the given variable
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  TYPE(datetime) :: next_tstep  ! Temporarily holds the value of the time at
                                ! timestep we are moving in to

  LOGICAL :: seek_required  !   T - the input files require a seek
                            !   F - the input files just need to be advanced
                            ! A seek of the input files is only required
                            ! when moving from spinup into the main run
                            ! Ramping between spinup cycles is handled by
                            ! time_varying_input_mod

  TYPE(datetime) :: dt  ! Placeholder for datetime in start/end of year
                        ! calculations


!-----------------------------------------------------------------------------


  next_tstep = datetime_advance(current_time, timestep_len)
  seek_required = .FALSE.  ! The default is to advance input files normally,
                           ! rather than seek to an out-of-sequence time

!-----------------------------------------------------------------------------
! The only conditions that can cause us to do something other than advance
! the time normally are:
!   - End of a spinup cycle
!   - End of main run
! So we check for these here
!-----------------------------------------------------------------------------
! Check if the current spinup cycle has ended
  IF ( is_spinup .AND. ( next_tstep >= spinup_end ) ) THEN

! Possible things to happen from this point are:
!   - Starting the next spinup cycle
!   - Terminating after a failed spinup
!   - Starting the main run
! All of these require a dump to be written for the timestep we are about to
! leave
    CALL write_dump()

    IF ( spinup_check() ) THEN
!-----------------------------------------------------------------------------
! All points are spun up, we can move on to the main run with no problems
!-----------------------------------------------------------------------------
      CALL log_info("next_time",                                              &
                    "All points spun up after " //                            &
                    TRIM(to_string(spinup_cycle)) // " cycles")
      CALL log_info("next_time", "Starting main run")

      is_spinup = .FALSE.
      next_tstep = main_run_start

! A seek of the input files is required when starting the main run
      seek_required = .TRUE.

    ELSE IF ( spinup_cycle >= max_spinup_cycles ) THEN
!-----------------------------------------------------------------------------
! We have reached the end of the last spinup cycle without all points being
! spun up
!-----------------------------------------------------------------------------
      CALL log_warn("next_time",                                              &
                    "Model is not fully spun up after maximum number " //     &
                    "of spinup cycles")

      IF ( terminate_on_spinup_fail ) THEN
        CALL log_fatal("next_time",                                           &
                       "Model failed to spin up - run terminated")
      ELSE
        CALL log_warn("next_time",                                            &
                      "Model failed to spin up - continuing with main run")

        is_spinup = .FALSE.
        next_tstep = main_run_start

! A seek of the input files is required when starting the main run
        seek_required = .TRUE.
      END IF

    ELSE
!-----------------------------------------------------------------------------
! Some points are still not spun up
!-----------------------------------------------------------------------------
      CALL log_info("next_time",                                              &
                    "Model is not fully spun up after " //                    &
                    TRIM(to_string(spinup_cycle)) // " cycles")
      CALL log_info("next_time",                                              &
                    "Starting another cycle of spinup")

      next_tstep = spinup_start
      spinup_cycle = spinup_cycle + 1

! A seek of the input files is not required for starting a new cycle of spinup
! advance_all handles ramping between spinup_end and spinup_start smoothly
    END IF

! Check for end of main run
  ELSE IF ( .NOT. is_spinup .AND. next_tstep >= main_run_end ) THEN
! We want to write a dump at the end of the last timestep
    CALL write_dump()
    end_of_run = .TRUE.
    RETURN

  ELSE IF ( end_of_year ) THEN
! We also want to write a dump if the timestep we are about to move out of
! was the last timestep in a year
! We do this as part of this IF statement to avoid writing dumps twice (i.e.
! if spinup cycles end on the last timestep of a year)
    CALL write_dump()
  END IF

!-----------------------------------------------------------------------------
! Now we have corrected the next time in light of spinup and stopped the run
! if required, we can update the necessary variables
!-----------------------------------------------------------------------------
! Set the current model time and advance the timestep
  current_time = next_tstep
  timestep = timestep + 1
  asteps_since_triffid = asteps_since_triffid + 1

! Seek the input data if required, otherwise just advance it normally
! A seek of the input files is only required when moving from spinup into
! the main run
! Ramping between spinup cycles is handled in advance_all
  IF ( seek_required ) THEN
    CALL seek_all_to_current_datetime()
  ELSE
    CALL advance_all()
  END IF

! Set the variables that indicate if we are starting or ending a year

! If going back one timestep length causes the year to change, then we are in
! the first timestep of a year
  dt = datetime_subtract(current_time, timestep_len)
  start_of_year = dt%year /= current_time%year
! If advancing the current time causes the year to change, then it is the last
! timestep in a year
  dt = datetime_advance(current_time, timestep_len)
  end_of_year = dt%year /= current_time%year


  CALL log_info("next_time",                                                  &
                "Timestep: " // TRIM(to_string(timestep)) // "; " //          &
                "Started at: " // datetime_to_string(current_time))
  IF ( is_spinup )                                                            &
    CALL log_info("next_time",                                                &
                  "Spinup cycle: " // TRIM(to_string(spinup_cycle)))

  RETURN

END SUBROUTINE next_time
#endif
