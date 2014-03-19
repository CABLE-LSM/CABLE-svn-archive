#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/init.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

SUBROUTINE init()

  USE time_varying_input_mod, ONLY : seek_all_to_current_datetime

  USE model_time_mod, ONLY : is_spinup

  USE init_grid_mod, ONLY : init_grid

  USE init_ancillaries_mod, ONLY : init_ancillaries

  USE init_params_mod, ONLY : init_params

  USE initial_conditions_mod, ONLY : init_ic

  USE spinup_mod, ONLY : spinup_init

  USE dump_mod, ONLY : write_dump

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   The main initialisation routine - initialises the model by calling
!   specialised routines
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------


! Initialise the model switches
  CALL init_switches()

! Initialise the vertical levels for the model
  CALL init_model_levels()

! Intialise the times for the run
  CALL init_time()

! Initialise the input, model and output grids
  CALL init_grid()

! Initialise the model ancils
  CALL init_ancillaries()

! Initialise model parameters
  CALL init_params()

! Initialise urban schemes
  CALL init_urban()

! Initialise IMOGEN
  CALL init_imogen()

! Initialise meteorological forcing
  CALL init_drive()

! Initialise other prescribed data
  CALL init_prescribed_data()

! Initialise the model prognostics
  CALL init_ic()

! Initialise output
  CALL init_output()


!-----------------------------------------------------------------------------
! Other initialisation that does not depend on user input
!-----------------------------------------------------------------------------
! Temporary (i.e. for this version) initialisation of variables.
  CALL init_vars_tmp()

! Set index arrays and initialise other variables.
  CALL init_parms()

! Seek the input files to the start of the run
  CALL seek_all_to_current_datetime()

! Save initial state if spinning up. Arrays are allocated here.
  IF ( is_spinup ) CALL spinup_init()

! Write an initial dump
  CALL write_dump()


  CALL log_info("init", "Initialisation is complete")

  RETURN

END SUBROUTINE init
#endif
