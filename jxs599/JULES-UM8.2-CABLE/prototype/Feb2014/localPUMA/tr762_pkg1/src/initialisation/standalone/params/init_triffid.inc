#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_triffid(nml_dir)

  USE c_mdi, ONLY :                                                           &
!  imported scalar parameters
     rmdi

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE nstypes, ONLY : npft

  USE switches, ONLY : l_phenol, l_triffid

  USE trif

  USE trif_io

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the TRIFFID parameters
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


! Nothing to do if neither triffid or phenology are selected
  IF ( .NOT. l_triffid .AND. .NOT. l_phenol ) RETURN

!-----------------------------------------------------------------------------
! Initialise namelist values before reading them
!-----------------------------------------------------------------------------
  crop_io(:) = -9

  g_area_io(:)  = RMDI
  g_grow_io(:)  = RMDI
  g_root_io(:)  = RMDI
  g_wood_io(:)  = RMDI
  lai_max_io(:) = RMDI
  lai_min_io(:) = RMDI

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_triffid", "Reading JULES_TRIFFID namelist...")

! Open the pft parameters namelist file
  OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'triffid_params.nml'),    &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_triffid",                                            &
                   "Error opening namelist file triffid_params.nml " //       &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  READ(NAMELIST_UNIT, nml=jules_triffid, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_triffid",                                            &
                   "Error reading namelist JULES_TRIFFID " //                 &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_triffid",                                            &
                   "Error closing namelist file triffid_params.nml " //       &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


!-----------------------------------------------------------------------------
! Copy values from fixed length arrays used in namelist into allocated arrays
!-----------------------------------------------------------------------------
  crop(:)    = crop_io(1:npft)
  g_area(:)  = g_area_io(1:npft)
  g_grow(:)  = g_grow_io(1:npft)
  g_root(:)  = g_root_io(1:npft)
  g_wood(:)  = g_wood_io(1:npft)
  lai_max(:) = lai_max_io(1:npft)
  lai_min(:) = lai_min_io(1:npft)

!-----------------------------------------------------------------------------
! Check that all required variables were present in the namelist.
! The namelist variables were initialised to RMDI.
!-----------------------------------------------------------------------------
  error = 0
  IF ( ANY( crop(:) < 0 ) ) THEN  ! crop was initialised to < 0
    error = 1
    CALL log_error("init_triffid", "No value for crop")
  END IF
  IF ( ANY( ABS( g_area(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_triffid", "No value for g_area")
  END IF
  IF ( ANY( ABS( g_grow(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_triffid", "No value for g_grow")
  END IF
  IF ( ANY( ABS( g_root(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_triffid", "No value for g_root")
  END IF
  IF ( ANY( ABS( g_wood(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_triffid", "No value for g_wood")
  END IF
  IF ( ANY( ABS( lai_max(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_triffid", "No value for lai_max")
  END IF
  IF ( ANY( ABS( lai_min(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_triffid", "No value for lai_min")
  END IF

  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_triffid",                                            &
                   "Variable(s) missing from namelist - see earlier " //      &
                   "error message(s)")

  RETURN

END SUBROUTINE init_triffid
#endif
