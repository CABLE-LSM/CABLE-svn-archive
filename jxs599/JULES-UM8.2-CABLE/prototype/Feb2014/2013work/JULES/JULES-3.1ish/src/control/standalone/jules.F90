#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/control/standalone/jules.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

PROGRAM jules

  USE time_varying_input_mod, ONLY :                                          &
    update_prescribed_variables => update_model_variables,                    &
    input_close_all => close_all

  USE update_mod, ONLY : update_derived_variables
  
  USE output_mod, ONLY : sample_data, output_data,                            &
                         output_close_all => close_all
  
  USE model_time_mod, ONLY : timestep, start_of_year, end_of_year, end_of_run

  USE switches, ONLY : l_imogen

  IMPLICIT NONE
  

!-----------------------------------------------------------------------------
! Initialise the model
!-----------------------------------------------------------------------------
  CALL init()

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
! Call the main model science routine
!-----------------------------------------------------------------------------
    CALL control(timestep)

!-----------------------------------------------------------------------------
! Update IMOGEN carbon if required
!-----------------------------------------------------------------------------
    IF ( l_imogen .AND. end_of_year ) CALL imogen_update_carb()

!-----------------------------------------------------------------------------
! The gathering of data for output and the actual outputting of data are
! done in two different phases
!-----------------------------------------------------------------------------
    CALL sample_data()
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

  WRITE(*,"(/,a)")'End of run.'

END PROGRAM jules
#endif
