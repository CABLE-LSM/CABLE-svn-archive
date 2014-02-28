#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE output_initial_data()

  USE datetime_mod, ONLY : datetime_subtract

  USE data_cube_mod, ONLY : cube_free

  USE model_time_mod, ONLY : is_spinup, spinup_start, main_run_start,         &
                             spinup_cycle, current_time

  USE file_ts_mod, ONLY : file_ts_write_var, file_ts_close

  USE model_interface_mod, ONLY : extract_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   For each profile, decide if we need to output a file containing initial
!   data for a section of the run
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  TYPE(file_ts) :: file  ! The initial data file

  CHARACTER(len=MAX_FILE_NAME_LEN) :: file_name
                                          ! The name for the initial data file

  TYPE(datetime) :: file_start  ! The start time for the file

  TYPE(data_cube) :: data_for_output  ! Cube to hold data for output

  INTEGER :: i,j  ! Loop counter


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! If we are not in the first timestep of a section, there is nothing to do
! for any profile
!-----------------------------------------------------------------------------
  IF ( is_spinup .AND. current_time /= spinup_start )                         &
    RETURN

  IF ( .NOT. is_spinup .AND. current_time /= main_run_start )                 &
    RETURN


!-----------------------------------------------------------------------------
! Check what we need to do for each profile
!-----------------------------------------------------------------------------
  DO i = 1,nprofiles

!-----------------------------------------------------------------------------
! Check if we need to output initial data for this profile
!-----------------------------------------------------------------------------
    IF ( .NOT. profiles(i)%output_initial )                                   &
! If initial state was not asked for, there is nothing to do
      CYCLE

    IF ( .NOT. ANY(                                                           &
            profiles(i)%fields(1:profiles(i)%nfields)%type == OUTPUT_SNAPSHOT &
       ))                                                                     &
! If there is no snapshot data in the profile, there is nothing to output here
      CYCLE

    IF ( is_spinup .AND. .NOT. profiles(i)%output_spinup )                    &
! If we are at the start of a spinup cycle but profile is not outputting during
! spinup, we don't need to do anything
      CYCLE

    IF ( .NOT. is_spinup ) THEN
! We are at the start of the main run
      IF ( .NOT. profiles(i)%output_main_run )                                &
! If the profile is not outputting for the main run, there is nothing to do
        CYCLE

      IF ( profiles(i)%output_start /= main_run_start )                       &
! If the profile does not start outputting until later in the run, there is nothing to do
        CYCLE
    END IF


!-----------------------------------------------------------------------------
! Open an initial data file for the profile
!-----------------------------------------------------------------------------
! All file names start with the run id followed by the profile name
    file_name = TRIM(run_id) // "." // TRIM(profiles(i)%profile_name)

    IF ( is_spinup ) THEN
! If we are in a spinup cycle, add that to the file name
      file_name = TRIM(file_name) // "." // "spin" // TRIM(to_string(spinup_cycle))
    END IF

! Indicate that it is an initial data file
    file_name = TRIM(file_name) // ".initial"

! Add an extension based on the output format
    SELECT CASE ( OUTPUT_FORMAT )
      CASE ( FORMAT_ASCII )
        file_name = TRIM(file_name) // ".asc"

      CASE ( FORMAT_NCDF )
        file_name = TRIM(file_name) // ".nc"

      CASE DEFAULT
        CALL log_fatal("output_initial_data",                                 &
                       "Unrecognised file format - " // TRIM(OUTPUT_FORMAT))
    END SELECT

! Prepend the output directory
    file_name = TRIM(output_dir) // "/" // TRIM(file_name)

! The timestamp in files is the end of the output period, so we create a file
! with exactly one output period ending at the current time
! We don't much care about time_bounds since the only real valid output is
! snapshot variables, for which the timestamp is the important thing
! So we create a file that ends at the current time, starts 1 second before
! and has an output period of 1
    file_start = datetime_subtract(current_time, 1)

    file = internal_open_output_file(                                         &
      file_start, current_time, 1, .FALSE., file_name,                        &
      profiles(i)%fields(1:profiles(i)%nfields)                               &
    )


!-----------------------------------------------------------------------------
! Populate any snapshot variables in the file
!-----------------------------------------------------------------------------
    DO j = 1,profiles(i)%nfields
! Ignore variables that are not snapshot variables
      IF ( profiles(i)%fields(j)%type /= OUTPUT_SNAPSHOT ) CYCLE

! Extract the data to output
      data_for_output = extract_var(profiles(i)%fields(j)%var_id)

! Output the data to file
      CALL file_ts_write_var(                                                 &
        file, profiles(i)%fields(j)%file_id, data_for_output,                 &
! Subgrid information (for writing a slab of the output grid in parallel mode)
        use_subgrid, subgrid                                                  &
      )

! Free the temporary cube
      CALL cube_free(data_for_output)
    END DO


!-----------------------------------------------------------------------------
! Close the file
!-----------------------------------------------------------------------------
    CALL file_ts_close(file)

  END DO

  RETURN

END SUBROUTINE output_initial_data
#endif
