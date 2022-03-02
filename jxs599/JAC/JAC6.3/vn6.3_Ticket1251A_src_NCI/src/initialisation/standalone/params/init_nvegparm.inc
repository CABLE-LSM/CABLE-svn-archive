#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_nvegparm_jules(nml_dir)

USE missing_data_mod, ONLY:                                                    &
!  imported scalar parameters
     rmdi

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE jules_surface_types_mod, ONLY: nnvg, npft

USE c_z0h_z0m, ONLY: z0h_z0m,  z0h_z0m_classic

USE nvegparm, ONLY: albsnc_nvg, albsnf_nvgu, albsnf_nvg, albsnf_nvgl,          &
                    catch_nvg, emis_nvg, gs_nvg, infil_nvg, z0_nvg, ch_nvg,    &
                    vf_nvg

USE nvegparm_io, ONLY: jules_nvegparm, init_nvegparm_allocated

USE jules_radiation_mod, ONLY: l_albedo_obs, l_snow_albedo

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the non-vegetation parameters
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists
! Work variables
INTEGER :: ERROR  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

CHARACTER(LEN=*), PARAMETER :: routinename='INIT_NVEGPARM_JULES'

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info(routinename, "Reading JULES_NVEGPARM namelist...")

! Open the nveg parameters namelist file
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'nveg_params.nml'),          &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = ERROR, &
               IOMSG = iomessage)
IF ( ERROR /= 0 )                                                              &
  CALL log_fatal(routinename,                                                  &
                 "Error opening namelist file nveg_params.nml " //             &
                 "(IOSTAT=" // TRIM(to_string(ERROR)) // " IOMSG=" //          &
                 TRIM(iomessage) // ")")

READ(namelist_unit, NML = jules_nvegparm, IOSTAT = ERROR, IOMSG = iomessage)
IF ( ERROR /= 0 )                                                              &
  CALL log_fatal(routinename,                                                  &
                 "Error reading namelist JULES_NVEGPARM " //                   &
                 "(IOSTAT=" // TRIM(to_string(ERROR)) // " IOMSG=" //          &
                 TRIM(iomessage) // ")")

! Close the namelist file
CLOSE(namelist_unit, IOSTAT = ERROR, IOMSG = iomessage)
IF ( ERROR /= 0 )                                                              &
  CALL log_fatal(routinename,                                                  &
                 "Error closing namelist file nveg_params.nml " //             &
                 "(IOSTAT=" // TRIM(to_string(ERROR)) // " IOMSG=" //          &
                 TRIM(iomessage) // ")")


!-----------------------------------------------------------------------------
! Process the namelist values
!-----------------------------------------------------------------------------
! Copy values from fixed length arrays used in namelist into allocated arrays
CALL init_nvegparm_allocated()

!-----------------------------------------------------------------------------
! Check that all required variables were present in the namelist.
! The namelist variables were initialised to rmdi.
! Some configurations don't need all parameters (e.g. albsnc_nvg) but for now
! we insist on getting all parameters (and that there are not rmdi!).
!-----------------------------------------------------------------------------
ERROR = 0
IF ( ANY( ABS( albsnf_nvg(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "No value for albsnf_nvg")
END IF
IF ( ANY( ABS( catch_nvg(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "No value for catch_nvg")
END IF
IF ( ANY( ABS( gs_nvg(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "No value for gs_nvg")
END IF
IF ( ANY( ABS( infil_nvg(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "No value for infil_nvg")
END IF
IF ( ANY( ABS( z0_nvg(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "No value for z0_nvg")
END IF
IF ( ANY( ABS( ch_nvg(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "No value for ch_nvg")
END IF
IF ( ANY( ABS( vf_nvg(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "No value for vf_nvg")
END IF
IF ( ANY( ABS( emis_nvg(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "No value for emis_nvg")
END IF
IF ( ANY( ABS( z0h_z0m(npft+1:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "No value for z0hm_nvg")
END IF
IF ( ANY( ABS( z0h_z0m_classic(npft+1:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "No value for z0hm_classic_nvg")
END IF

IF ( l_albedo_obs ) THEN
  IF ( ANY( ABS( albsnf_nvgl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL log_error(routinename, "No value for albsnf_nvgl")
  END IF
  IF ( ANY( ABS( albsnf_nvgu(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL log_error(routinename, "No value for albsnf_nvgu")
  END IF
END IF

IF ( .NOT. l_snow_albedo ) THEN
  IF ( ANY( ABS( albsnc_nvg(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL log_error(routinename, "No value for albsnc_nvg")
  END IF
END IF

IF ( ERROR /= 0 )                                                              &
  CALL log_fatal(routinename,                                                  &
                 "Variable(s) missing from namelist - see earlier " //         &
                 "error message(s)")

ERROR = 0
IF ( ANY( vf_nvg(:) > 1.0 ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "At least one value of vf_nvg > 1.0")
END IF
IF ( ANY( vf_nvg(:) < 0.0 ) ) THEN
  ERROR = 1
  CALL log_error(routinename, "At least one value of vf_nvg < 0.0")
END IF

IF ( ERROR /= 0 )                                                              &
   CALL log_fatal(routinename,                                                 &
                  "Variable(s) outside of allowed range - see earlier " //     &
                  "error message(s)")

RETURN

END SUBROUTINE init_nvegparm_jules

#endif
