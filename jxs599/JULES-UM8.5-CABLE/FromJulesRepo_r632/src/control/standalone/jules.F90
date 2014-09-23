#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


PROGRAM jules

  USE io_constants, ONLY : MAX_FILE_NAME_LEN

  USE time_varying_input_mod, ONLY :                                          &
    update_prescribed_variables => update_model_variables,                    &
    input_close_all => close_all

  USE update_mod, ONLY : update_derived_variables
  
  USE output_mod, ONLY : output_initial_data, sample_data, output_data,       &
                         output_close_all => close_all
  
  USE model_time_mod, ONLY : timestep, start_of_year, end_of_year, end_of_run

  USE switches, ONLY : l_imogen

  USE logging_mod, ONLY : log_init, log_shutdown

  IMPLICIT NONE
  
!-----------------------------------------------------------------------------
! Description:
!   This is the main program routine for JULES
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  CHARACTER(len=MAX_FILE_NAME_LEN) :: nml_dir  ! Directory containing namelists

  INTEGER :: error  ! Error indicator


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise the MPI environment
!-----------------------------------------------------------------------------
! We don't check the error since most (all?) MPI implementations will just
! fail if a call is unsuccessful
  CALL MPI_INIT(error)

!-----------------------------------------------------------------------------
! Try to read a single argument from the command line
!
! If present, that single argument will be the directory we try to read
! namelists from
! If not present, we use current working directory instead
!-----------------------------------------------------------------------------
  CALL GET_COMMAND_ARGUMENT(1, nml_dir)
! If no argument is given, GET_COMMAND_ARGUMENT returns a blank string
  IF ( LEN_TRIM(nml_dir) == 0 ) nml_dir = "."

!-----------------------------------------------------------------------------
! Now MPI_INIT has been called, we can initialise the logging
!-----------------------------------------------------------------------------
  CALL log_init(nml_dir)

!-----------------------------------------------------------------------------
! Initialise the model
!-----------------------------------------------------------------------------
  CALL init(nml_dir)

!-----------------------------------------------------------------------------
! Loop over timesteps.
! Note that the number of timesteps is of unknown length at the start of run,
! if the model is to determine when it has spun up.
!-----------------------------------------------------------------------------
  DO    !  timestep

!-----------------------------------------------------------------------------
! Update the IMOGEN climate variables if required
!-----------------------------------------------------------------------------
    IF ( l_imogen .AND. start_of_year ) CALL imogen_update_clim()

!-----------------------------------------------------------------------------
! The update of prescribed data is done in two phases
!  - Update variables provided by files
!  - Update variables that are derived from those updated in the first phase
!-----------------------------------------------------------------------------
    CALL update_prescribed_variables()
    CALL update_derived_variables()

!-----------------------------------------------------------------------------
! Check if this is a timestep that we need to output initial data for (i.e.
! start of spinup cycle or start of main run), and output that data if
! required
!-----------------------------------------------------------------------------
    CALL output_initial_data()

!-----------------------------------------------------------------------------
! Call the main model science routine
!-----------------------------------------------------------------------------
    CALL control(timestep)

!-----------------------------------------------------------------------------
! Update IMOGEN carbon if required
!-----------------------------------------------------------------------------
    IF ( l_imogen .AND. end_of_year ) CALL imogen_update_carb()

!-----------------------------------------------------------------------------
! Sample variables for output
!-----------------------------------------------------------------------------
    CALL sample_data()

!-----------------------------------------------------------------------------
! Output collected data if required
!-----------------------------------------------------------------------------
    CALL output_data()

!-----------------------------------------------------------------------------
! Move the model on to the next timestep
!-----------------------------------------------------------------------------
    CALL next_time()

    IF ( end_of_run ) EXIT

  ENDDO  !  timestep loop

!-----------------------------------------------------------------------------
! Clean up by closing all open files
!-----------------------------------------------------------------------------
  CALL input_close_all()
  CALL output_close_all()

!-----------------------------------------------------------------------------
! Shut down the logging
!-----------------------------------------------------------------------------
  CALL log_shutdown()

!-----------------------------------------------------------------------------
! Clean up the MPI environment
!-----------------------------------------------------------------------------
  CALL MPI_FINALIZE(error)


END PROGRAM jules
#endif
