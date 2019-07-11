#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_advance(file)

  USE io_constants, ONLY : MODE_WRITE

  USE datetime_mod, ONLY : datetime_advance, datetime_to_string

  USE templating_mod, ONLY : tpl_substitute_datetime

  USE file_gridded_mod, ONLY : file_gridded_advance

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Advances the file by one timestep
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument type
  TYPE(file_ts), INTENT(INOUT) :: file


! Work variables
  TYPE(datetime) :: next_time  ! The datetime that we are trying to advance to

! These are used only if advancing requires opening a new file
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file_name  ! The name of the file
  INTEGER :: file_index  ! The index of the file in the files array of the
                         ! file_ts object
  TYPE(datetime) :: next_file_start  ! The start time of the next file after
                                     ! the one we need to open


!-----------------------------------------------------------------------------


! If file is still in define mode, error out
  IF ( file%define_mode )                                                     &
    CALL log_fatal("file_ts_advance",                                         &
                   "Cannot seek while file is in define mode")


  next_time = datetime_advance(                                               &
    file%current_datetime, file%data_period                                   &
  )

  IF ( next_time >= file%data_end ) THEN
    IF ( file%is_climatology ) THEN
! If we are using a climatology, then advancing past the end of the data
! just means we go back to the start
      CALL file_ts_seek_to_datetime(file, file%data_start)
      RETURN
    ELSE
! If we are not using a climatology, stop with an error as we cannot advance
! beyond the end of the data
      CALL log_fatal("file_ts_advance",                                       &
                     "Cannot advance to " // datetime_to_string(next_time) // &
                     " - data_end is " // datetime_to_string(file%data_end))
    END IF
  END IF

! Check if we need to open a new file
  IF ( next_time >= file%next_file_start ) THEN
!-----------------------------------------------------------------------------
! We need to get the name of the file to open and the start time of the next
! file so we can update the file_ts object
!
! We know that we will be using the first timestep in the file, so we don't
! need to calculate where in the file to start
!
! This is done differently depending on whether we are using time templating
! or not
!-----------------------------------------------------------------------------
! Initialise file_index - not used when time templating is in use
    file_index = 0

    IF ( file%use_time_template ) THEN
! Time templating is in use
! We already know that the next time is within the data range, so we get the
! file to open by just substituting values into the template
      file_name = tpl_substitute_datetime(file%template, next_time)
! We get the start time for the next file by skipping one templating period
! on from the previous value
      next_file_start = datetime_advance(                                     &
        file%next_file_start, file%tpl_period                                 &
      )
    ELSE
! No time templating - use the list of files instead

! We know that the file we are moving out of is not the last file, otherwise
! the test on file%data_end earlier would fail
      file_index = file%open_file_index + 1

! So we can easily get the name of the file now
      file_name = file%files(file_index)

! Get the start time for the next file
! If we are in the moving in to the last file (i.e. there is no next file),
! then set this to 1 second past the end of the data so it never gets hit
      IF ( file_index == file%nfiles ) THEN
        next_file_start = datetime_advance(file%data_end, 1)
      ELSE
        next_file_start = file%file_times(file_index + 1)
      END IF
    END IF

! Make sure we have the correct file open
    CALL file_ts_internal_open_file(file, file_name)

! Update open file related properties on the file_ts object
    file%open_file_index  = file_index
    file%next_file_start  = next_file_start

  ELSE
!-----------------------------------------------------------------------------
! No new file required - all we have to do is advance the current file
!-----------------------------------------------------------------------------
    CALL file_gridded_advance(file%open_file)
  END IF

! We always update the current datetime
  file%current_datetime = next_time


! Indicate that the timestamp for this file timestep has not been written yet
  file%timestamp_written = .FALSE.

  RETURN

END SUBROUTINE file_ts_advance
#endif
