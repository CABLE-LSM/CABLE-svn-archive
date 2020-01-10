#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE init_ancillaries(nml_dir)

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE ancil_info, ONLY: frac_surft

USE aero, ONLY: co2_mmr_mod => co2_mmr

USE dump_mod, ONLY: ancil_dump_read

USE errormessagelength_mod, ONLY: errormessagelength

USE jules_surface_mod, ONLY: l_urban2t

USE switches_urban, ONLY: l_moruses

USE lsm_switch_mod, ONLY: lsm_id, cable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in information about the model ancillaries
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

CHARACTER(LEN=errormessagelength) :: iomessage

! Work variables
INTEGER :: error  ! Error indicator

LOGICAL :: read_from_dump
REAL(KIND=real_jlslsm)    :: co2_mmr

CHARACTER(LEN=*), PARAMETER :: RoutineName='init_ancillaries'

NAMELIST  / jules_co2 / read_from_dump, co2_mmr

!-----------------------------------------------------------------------------
! Open the ancillaries namelist file
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'ancillaries.nml'),         &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,&
               IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error opening namelist file ancillaries.nml " //            &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

! Defer to specialist routines to process each namelist
CALL init_frac()
CALL init_soil_props()
CALL init_top()
CALL init_pdm()
CALL init_agric()
CALL init_crop_props()
CALL init_irrig()
CALL init_rivers_props()
!CABLE_LSM:
IF (lsm_id == cable) THEN
  CALL init_cable()
END IF
IF ( l_urban2t .OR. l_moruses ) CALL init_urban_props()

! Read the JULES_CO2 namelist
! This is so simple it doesn't really need its own subroutine
read_from_dump = .FALSE.

!Copy across the default value set up in the module.
co2_mmr = co2_mmr_mod

CALL log_info(RoutineName, "Reading JULES_CO2 namelist...")
READ(namelist_unit, NML = jules_co2, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error reading namelist JULES_CO2 " //                       &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

ancil_dump_read%co2 = read_from_dump

IF ( .NOT. ancil_dump_read%co2) THEN
  co2_mmr_mod = co2_mmr
ELSE
  ! We read from the dump file.
  CALL log_info(RoutineName,                                                  &
                "co2_mmr will be read from the dump file.  " //               &
                "Namelist value ignored")
END IF

! Close the namelist file
CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 ) THEN
  CALL log_fatal(RoutineName,                                                 &
                 "Error closing namelist file ancillaries.nml " //            &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")
END IF

RETURN

END SUBROUTINE init_ancillaries
#endif
