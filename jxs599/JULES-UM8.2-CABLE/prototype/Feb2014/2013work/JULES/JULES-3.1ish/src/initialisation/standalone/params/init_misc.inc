#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/params/init_misc.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_misc()

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE surf_param, ONLY : jules_surf_param

  USE csmin, ONLY : jules_csmin

  USE aero, ONLY : jules_aero

  USE seed, ONLY : jules_seed

  USE sigm, ONLY : jules_sigm

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the surface parameters that have not been initialised yet...
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: error  ! Error indicator


!-----------------------------------------------------------------------------

! Open the misc parameters namelist file
  OPEN(NAMELIST_UNIT, FILE='misc_params.nml',                                 &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_misc",                                               &
                   "Error opening namelist file misc_params.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


! There are five(?!) namelists to read from this file
  CALL log_info("init_misc", "Reading JULES_SURF_PARAM namelist...")
  READ(NAMELIST_UNIT, nml=jules_surf_param, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_misc",                                               &
                   "Error reading namelist JULES_SURF_PARAM " //              &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CALL log_info("init_misc", "Reading JULES_CSMIN namelist...")
  READ(NAMELIST_UNIT, nml=jules_csmin, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_misc",                                               &
                   "Error reading namelist JULES_CSMIN " //                   &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CALL log_info("init_misc", "Reading JULES_AERO namelist...")
  READ(NAMELIST_UNIT, nml=jules_aero, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_misc",                                               &
                   "Error reading namelist JULES_AERO " //                    &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CALL log_info("init_misc", "Reading JULES_SEED namelist...")
  READ(NAMELIST_UNIT, nml=jules_seed, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_misc",                                               &
                   "Error reading namelist JULES_SEED " //                    &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CALL log_info("init_misc", "Reading JULES_SIGM namelist...")
  READ(NAMELIST_UNIT, nml=jules_sigm, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_misc",                                               &
                   "Error reading namelist JULES_SIGM " //                    &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_misc",                                               &
                   "Error closing namelist file misc_params.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! No post processing required!!

  RETURN

END SUBROUTINE init_misc
#endif
