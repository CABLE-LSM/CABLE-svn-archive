#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
MODULE init_rivers_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE init_rivers(nml_dir)
!-----------------------------------------------------------------------------
! Description:
!   Reads in the river routing namelist items and checks them for consistency
!   Initialises river routing parameters
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
USE mpi, ONLY: mpi_comm_world

USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE jules_rivers_mod, ONLY: jules_rivers, check_jules_rivers,                 &
                            l_rivers, i_river_vn, nstep_rivers
USE logging_mod, ONLY: log_info, log_debug, log_warn, log_error, log_fatal
USE overbank_inundation_mod, ONLY: jules_overbank

IMPLICIT NONE


! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists

! Work variables
INTEGER :: error  ! Error indicator
INTEGER :: ntasks ! Parallel mode indicator

!-----------------------------------------------------------------------------
! Read river routing namelist
!----------------------------------------------------------------------------
CALL log_info("init_rivers", "Reading JULES_RIVERS namelist...")

! Open the river routing parameters namelist file
OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'jules_rivers.nml'),        &
               STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_rivers",                                               &
                 "Error opening namelist file jules_rivers.nml " //           &
                 "(IOSTAT=" // TRIM(to_string(error)) // ")")

READ(namelist_unit, NML = jules_rivers, IOSTAT = error)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_rivers",                                               &
                 "Error reading namelist JULES_RIVERS " //                    &
                 "(IOSTAT=" // TRIM(to_string(error)) // ")")

IF (l_rivers) THEN
  CALL log_info("init_rivers", "Reading JULES_OVERBANK namelist...")
  READ(namelist_unit, NML = jules_overbank, IOSTAT = error)
  IF ( error /= 0 ) THEN
    CALL log_fatal("init_overbank",                                           &
                   "Error reading namelist JULES_OVERBANK " //                &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")
  END IF
END IF

! Close the namelist file
CLOSE(namelist_unit, IOSTAT = error)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_rivers",                                               &
                 "Error closing namelist file jules_rivers.nml " //           &
                 "(IOSTAT=" // TRIM(to_string(error)) // ")")

CALL check_jules_rivers()

IF (l_rivers) THEN
  CALL log_info("init_rivers",                                                &
                "river routing" // TRIM(to_string(i_river_vn)) //             &
                "is selected")
  CALL log_info("init_rivers",                                                &
                "Integer river routing timestep = " //                        &
                 TRIM(to_string(nstep_rivers)))
ELSE
  CALL log_info("init_rivers", "No river routing selected")
END IF

END SUBROUTINE init_rivers

END MODULE init_rivers_mod
#endif
