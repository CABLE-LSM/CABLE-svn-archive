#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_ancillaries(nml_dir)

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
! Arguments
  CHARACTER(len=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                           ! namelists

! Work variables
  INTEGER :: error  ! Error indicator


!-----------------------------------------------------------------------------
! The jules_pdm namelist
!-----------------------------------------------------------------------------
  NAMELIST /jules_pdm/ b_pdm, dz_pdm


!-----------------------------------------------------------------------------


! Open the ancillaries namelist file
  OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'ancillaries.nml'),       &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_ancillaries",                                        &
                   "Error opening namelist file ancillaries.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Defer to specialist routines to process each namelist
  CALL init_frac()
  CALL init_soil()

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
