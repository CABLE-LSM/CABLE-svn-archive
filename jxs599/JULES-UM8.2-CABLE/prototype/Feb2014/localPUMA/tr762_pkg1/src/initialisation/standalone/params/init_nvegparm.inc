#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_nvegparm(nml_dir)

  USE c_mdi, ONLY :                                                           &
!  imported scalar parameters
     rmdi

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE switches_urban, ONLY : l_urban2t

  USE nstypes, ONLY : nnvg, npft, urban_canyon, urban_roof

  USE c_z0h_z0m, ONLY : z0h_z0m,  z0h_z0m_classic

  USE nvegparm

  USE nvegparm_io

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the non-vegetation parameters
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


!-----------------------------------------------------------------------------
! Initialise namelist values before reading them
!-----------------------------------------------------------------------------
  albsnc_nvg_io(:)       = RMDI
  albsnf_nvgu_io(:)      = RMDI
  albsnf_nvg_io(:)       = RMDI
  albsnf_nvgl_io(:)      = RMDI
  catch_nvg_io(:)        = RMDI
  gs_nvg_io(:)           = RMDI
  infil_nvg_io(:)        = RMDI
  z0_nvg_io(:)           = RMDI
  ch_nvg_io(:)           = RMDI
  vf_nvg_io(:)           = RMDI
  emis_nvg_io(:)         = RMDI
  z0hm_nvg_io(:)         = RMDI
  z0hm_classic_nvg_io(:) = RMDI

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_nvegparm", "Reading JULES_NVEGPARM namelist...")

! Open the pft parameters namelist file
  OPEN(NAMELIST_UNIT, FILE=(TRIM(nml_dir) // '/' // 'nveg_params.nml'),       &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_nvegparm",                                           &
                   "Error opening namelist file nveg_params.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  READ(NAMELIST_UNIT, nml=jules_nvegparm, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_nvegparm",                                           &
                   "Error reading namelist JULES_NVEGPARM " //                &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_nvegparm",                                           &
                   "Error closing namelist file nveg_params.nml " //          &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


!-----------------------------------------------------------------------------
! Process the namelist values
!-----------------------------------------------------------------------------
! Copy values from fixed length arrays used in namelist into allocated arrays
  nvgname(:)    = nvgname_io(1:nnvg)
  albsnc_nvg(:) = albsnc_nvg_io(1:nnvg)
  albsnf_nvgu(:)= albsnf_nvgu_io(1:nnvg)
  albsnf_nvg(:) = albsnf_nvg_io(1:nnvg)
  albsnf_nvgl(:)= albsnf_nvgl_io(1:nnvg)
  catch_nvg(:)  = catch_nvg_io(1:nnvg)
  gs_nvg(:)     = gs_nvg_io(1:nnvg)
  infil_nvg(:)  = infil_nvg_io(1:nnvg)
  z0_nvg(:)     = z0_nvg_io(1:nnvg)
  ch_nvg(:)     = ch_nvg_io(1:nnvg)
  vf_nvg(:)     = vf_nvg_io(1:nnvg)
  emis_nvg(:)   = emis_nvg_io(1:nnvg)

  z0h_z0m(npft+1:) = z0hm_nvg_io(1:nnvg)
  z0h_z0m_classic(npft+1:) = z0hm_classic_nvg_io(1:nnvg)

!-----------------------------------------------------------------------------
! Check that all required variables were present in the namelist.
! The namelist variables were initialised to RMDI.
! Some configurations don't need all parameters (e.g. albsnc_nvg) but for now
! we insist on getting all parameters (and that there are not RMDI!).
!-----------------------------------------------------------------------------
  error = 0
  IF ( ANY( ABS( albsnc_nvg(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for albsnc_nvg")
  END IF
  IF ( ANY( ABS( albsnf_nvgu(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for albsnf_nvgu")
  END IF
  IF ( ANY( ABS( albsnf_nvg(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for albsnf_nvg")
  END IF
  IF ( ANY( ABS( albsnf_nvgl(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for albsnf_nvgl")
  END IF
  IF ( ANY( ABS( catch_nvg(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for catch_nvg")
  END IF
  IF ( ANY( ABS( gs_nvg(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for gs_nvg")
  END IF
  IF ( ANY( ABS( infil_nvg(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for infil_nvg")
  END IF
  IF ( ANY( ABS( z0_nvg(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for z0_nvg")
  END IF
  IF ( ANY( ABS( ch_nvg(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for ch_nvg")
  END IF
  IF ( ANY( ABS( vf_nvg(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for vf_nvg")
  END IF
  IF ( ANY( ABS( emis_nvg(:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for emis_nvg")
  END IF
  IF ( ANY( ABS( z0h_z0m(npft+1:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for z0hm_nvg")
  END IF
  IF ( ANY( ABS( z0h_z0m_classic(npft+1:) - RMDI ) < EPSILON(1.0) ) ) THEN
    error = 1
    CALL log_error("init_nvegparm", "No value for z0hm_classic_nvg")
  END IF

  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_nvegparm",                                           &
                   "Variable(s) missing from namelist - see earlier " //      &
                   "error message(s)")

  RETURN

END SUBROUTINE init_nvegparm
#endif
