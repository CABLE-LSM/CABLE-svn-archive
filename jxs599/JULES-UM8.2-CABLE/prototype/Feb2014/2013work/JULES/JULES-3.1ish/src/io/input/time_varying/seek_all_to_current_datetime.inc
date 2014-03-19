#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/input/time_varying/seek_all_to_current_datetime.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE seek_all_to_current_datetime()

  USE datetime_mod, ONLY : PERIOD_MONTH, PERIOD_YEAR,                         &
                           datetime_clone, datetime_diff,                     &
                           datetime_advance, datetime_subtract

  USE model_time_mod, ONLY : current_time, timestep_len

  USE file_ts_mod, ONLY : file_ts_seek_to_datetime, file_ts_read_var,         &
                          file_ts_advance

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Sets up the input files so that they are able to provide data for the
!   current timestep next time they are asked for data
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  TYPE(datetime) :: data_dt  ! The datetime for data that applies for
                             ! current_dt
                             ! See step 1 in algorithm description below
  TYPE(datetime) :: seek_dt  ! The datetime to seek the underlying file to
                             ! See step 2 of algorithm below

  INTEGER :: diff_secs  ! Used in calculation of tsteps_in_data_period
                        ! The number of seconds between two times

  INTEGER :: i,j,k  ! Loop counters

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Overview of algorithm (for each file):
!
!   1. Find the data time for the current datetime (i.e. the closest
!      time going backwards from current_dt that data exists for)
!   2. Depending what times are needed for interpolation, seek the underlying
!      file to the appropriate time
!   3. Read the timesteps required for interpolation
!   4. Calculate the data times for all times from -1 to 2, even if they are
!      not required (and even if they are not actually part of the data)
!   5. Calculate the number of model timesteps in each data period and the
!      number of timesteps into the current interpolation period we are
!
!-----------------------------------------------------------------------------


  DO i = 1,nfiles
!-----------------------------------------------------------------------------
! Step 1 - find the data time for the requested datetime for this file
!-----------------------------------------------------------------------------
    SELECT CASE ( files(i)%fh%data_period )
      CASE ( 1: )
! In the case of a period of seconds >= 1, we need to know how many seconds
! past the closest data time current_dt is, and take that off current_dt
        data_dt = datetime_subtract(                                          &
          current_time,                                                       &
          MOD(                                                                &
            datetime_diff(files(i)%fh%data_start, current_time),              &
            files(i)%fh%data_period                                           &
          )                                                                   &
        )

      CASE ( PERIOD_MONTH )
! In the case of a monthly period, we know that the data applies from the
! start of the month
        data_dt = datetime_clone(current_time)
        data_dt%day  = 1
        data_dt%time = 0

      CASE ( PERIOD_YEAR )
! In the case of a yearly period, we know that the data will apply from the
! start of the year
        data_dt = datetime_clone(current_time)
        data_dt%month = 1
        data_dt%day   = 1
        data_dt%time  = 0
    END SELECT

! Update the time of the current data
    files(i)%data_times(0) = data_dt


!-----------------------------------------------------------------------------
! Step 2 - depending on the interpolation in use for this file, seek the
!          underlying timeseries file to the correct place
!-----------------------------------------------------------------------------
    seek_dt = datetime_clone(data_dt)
! Because data_period could be a special period, which we can't meaningfully
! multiply, we do this by looping
    DO j = 1,ABS(files(i)%times_lbound)
      IF ( files(i)%times_lbound < 0 ) THEN
        seek_dt = datetime_subtract(seek_dt, files(i)%fh%data_period)
      ELSE
        seek_dt = datetime_advance(seek_dt, files(i)%fh%data_period)
      END IF
    END DO

    CALL file_ts_seek_to_datetime(files(i)%fh, seek_dt)


!-----------------------------------------------------------------------------
! Step 3 - read the number of timesteps required for interpolation for
!          each field
!-----------------------------------------------------------------------------
    DO j = files(i)%times_lbound,files(i)%times_ubound
      DO k = 1,files(i)%nfields
! Read the slab of data for this time
        CALL file_ts_read_var(                                                &
          files(i)%fh,                                                        &
          files(i)%fields(k)%file_id,                                         &
          files(i)%fields(k)%data(:,:,:,j),                                   &
! Subgrid information from input_grid_mod
          use_subgrid, points, x_coords, y_coords                             &
        )
      END DO  ! fields

! If we need to read more data now, advance the file
      IF ( j /= files(i)%times_ubound ) CALL file_ts_advance(files(i)%fh)
    END DO  ! times


!-----------------------------------------------------------------------------
! Step 4 - Populate data_times for all times -1 to 2
!-----------------------------------------------------------------------------
! We know data_times(0), so populate the other times from that
    files(i)%data_times(-1) = datetime_subtract(                              &
      files(i)%data_times(0), files(i)%fh%data_period                         &
    )

    DO j = 1,2
      files(i)%data_times(j) = datetime_advance(                              &
        files(i)%data_times(j - 1), files(i)%fh%data_period                   &
      )
    END DO

!-----------------------------------------------------------------------------
! Step 5 - Populate tsteps_in_data_period and current_tstep
!-----------------------------------------------------------------------------
! Populate tsteps_in_data_period from the data times
    DO j = -1,1
      diff_secs = datetime_diff(                                              &
        files(i)%data_times(j), files(i)%data_times(j + 1)                    &
      )

      IF ( MOD(diff_secs, timestep_len) /= 0 )                                &
        CALL log_fatal("seek_all_to_current_datetime",                        &
                       "Data should be a whole number of model timesteps apart")

      files(i)%tsteps_in_data_period(j) = REAL(diff_secs) / REAL(timestep_len)
    END DO

! Work out how many model timesteps into the interpolation period we are
    diff_secs = datetime_diff(files(i)%data_times(0), current_time)

    IF ( MOD(diff_secs, timestep_len) /= 0 )                                  &
      CALL log_fatal("seek_all_to_current_datetime",                          &
                     "Data should be a whole number of model timesteps apart")

    files(i)%current_tstep = REAL(diff_secs) / REAL(timestep_len)

  END DO  ! files

  RETURN

END SUBROUTINE seek_all_to_current_datetime
#endif
