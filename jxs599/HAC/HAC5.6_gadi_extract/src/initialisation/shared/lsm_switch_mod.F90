! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE lsm_switch_mod

!-----------------------------------------------------------------------------
! Description:
!   CABLE_LSM:
!   Declare default and read namelist overide for LSM switch between
!   JULES & CABLE for coupled runs and standalone runs
!
! Current Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!
! NB: Checks have not been implemented at this time
!-----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER, PARAMETER :: unset = 0
INTEGER, PARAMETER :: jules = 1
INTEGER, PARAMETER :: cable = 2

INTEGER :: lsm_id = unset
CHARACTER(LEN=5), DIMENSION(2) :: lsm_name = (/ 'JULES', 'CABLE' /)

!-----------------------------------------------------------------------------
! Namelist definition
!-----------------------------------------------------------------------------
NAMELIST  / jules_lsm_switch/ lsm_id


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSM_SWITCHES_MOD'

CONTAINS

#if !defined(UM_JULES)
SUBROUTINE init_lsm(nml_dir)

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE errormessagelength_mod, ONLY: errormessagelength

USE logging_mod, ONLY: log_info, log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!  Read the LSM_SWITCHES namelist for standalone runs
!-----------------------------------------------------------------------------
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists


!-----------------------------------------------------------------------------
! Work variables
!-----------------------------------------------------------------------------
INTEGER :: error   ! Error indicators
CHARACTER(LEN=errormessagelength) :: iomessage
!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_lsm", "Reading JULES_LSM_SWITCH namelist...")

OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'jules_lsm_switch.nml'),    &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT  = error, &
               IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_lsm",                                                  &
                 "Error opening namelist file jules_lsm_switch.nml " //       &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

READ(namelist_unit, NML = jules_lsm_switch, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_lsm",                                                  &
                 "Error reading namelist JULES_LSM_SWITCH " //                &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_lsm",                                                  &
                 "Error closing namelist file jules_lsm_switch.nml " //       &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

!-----------------------------------------------------------------------------
! Check if land surface model id is valid
!-----------------------------------------------------------------------------

IF ( (lsm_id /= cable) .AND. (lsm_id /= jules) ) THEN
  CALL log_fatal("init_lsm", "Invalid LSM ID provided.")
ELSE
  CALL log_info("init_lsm", "Land surface model selected is " //              &
                lsm_name(lsm_id))
END IF

RETURN

END SUBROUTINE init_lsm
#endif

SUBROUTINE init_lsm_um()
! Sets the LSM ID for the UM. Currently only JULES.

IMPLICIT NONE

lsm_id = jules

END SUBROUTINE init_lsm_um

END MODULE lsm_switch_mod


