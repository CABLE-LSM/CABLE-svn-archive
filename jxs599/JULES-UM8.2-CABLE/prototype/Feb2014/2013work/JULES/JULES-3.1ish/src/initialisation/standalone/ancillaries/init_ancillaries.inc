#if !defined(UM_JULES)
#define CABLE_init T
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/ancillaries/init_ancillaries.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_ancillaries()

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE switches, ONLY : l_pdm, b_pdm, dz_pdm

  USE ancil_info, ONLY : frac

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about the model ancillaries
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
! The jules_pdm namelist
!-----------------------------------------------------------------------------
  NAMELIST /jules_pdm/ b_pdm, dz_pdm


!-----------------------------------------------------------------------------


! Open the ancillaries namelist file
  OPEN(NAMELIST_UNIT, FILE='ancillaries.nml',                                 &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_ancillaries",                                        &
                   "Error opening namelist file ancillaries.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Defer to specialist routines to process each namelist
  CALL init_frac()
  CALL init_soil()

! CABLE
#ifdef CABLE_init
  CALL init_cable()
#endif

! The pdm namelist is easy - just 2 properties - and we only need to read it
! if pdm is on
  IF ( l_pdm ) THEN
    READ(NAMELIST_UNIT, nml=jules_pdm, IOSTAT=error)
    IF ( error /= 0 )                                                         &
      CALL log_fatal("init_ancillaries",                                      &
                     "Error reading namelist JULES_PDM " //                   &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")
  END IF

  CALL init_top()
  CALL init_agric()

! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_ancillaries",                                        &
                   "Error closing namelist file ancillaries.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  RETURN

END SUBROUTINE init_ancillaries
#endif
