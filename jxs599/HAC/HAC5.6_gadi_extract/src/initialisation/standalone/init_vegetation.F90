#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
MODULE init_vegetation_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE init_vegetation(nml_dir)
!-----------------------------------------------------------------------------
! Description:
!   Reads in the vegetation namelist items and checks them for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE jules_vegetation_mod

USE logging_mod, ONLY: log_info, log_error, log_fatal

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists
! Work variables
INTEGER :: error  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------
! First, we read the vegetation namelist
!-----------------------------------------------------------------------------
CALL log_info("init_vegetation", "Reading JULES_VEGETATION namelist...")

OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'jules_vegetation.nml'),    &
     STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = error,          &
     IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_vegetation",                                           &
                 "Error opening namelist file jules_vegetation.nml " //       &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

READ(namelist_unit, NML = jules_vegetation, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_vegetation",                                           &
                 "Error reading namelist JULES_VEGETATION " //                &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

CLOSE(namelist_unit, IOSTAT = error, IOMSG = iomessage)
IF ( error /= 0 )                                                             &
  CALL log_fatal("init_vegetation",                                           &
                 "Error closing namelist file vegetation.nml " //             &
                 "(IOSTAT=" // TRIM(to_string(error)) // " IOMSG=" //         &
                 TRIM(iomessage) // ")")

CALL check_jules_vegetation()

!-----------------------------------------------------------------------------
! Print some human friendly summary information about the selected options
!-----------------------------------------------------------------------------
IF ( l_phenol )                                                               &
  CALL log_info("init_vegetation", "Phenology is on")

IF ( l_triffid )                                                              &
  CALL log_info("init_vegetation", "TRIFFID is on")

IF ( l_veg_compete )                                                          &
  CALL log_info("init_vegetation", "Competing vegetation is on")
  
IF ( l_veg_compete .AND. ( .NOT. l_ht_compete ) )                             &
  CALL log_info("init_vegetation", "l_veg_compete=T, l_ht_compete=F " //      &
                "assumes that the " //                                        &
                "natural PFTs are BT, NT, C3, C4, SH (in that order).")

IF ( l_crop )                                                                 &
  CALL log_info("init_vegetation", "Crop model is on")

IF ( l_irrig_dmd )                                                            &
  CALL log_info("init_vegetation", "Irrigation demand model is on")

SELECT CASE ( irr_crop )
CASE ( 0 )
  CALL log_info("init_vegetation",                                            &
                "irr_crop = 0: continuous irrigation")
CASE ( 1 )
  CALL log_info("init_vegetation",                                            &
                "irr_crop = 1: following Doell & Siebert (2002)")
CASE ( 2 )
  CALL log_info("init_vegetation",                                            &
                "irr_crop = 2: irrigation triggered by DVI from crop model")
CASE DEFAULT
  CALL log_error("init_vegetation",                                           &
                 "irr_crop value not recognised")
END SELECT

IF ( l_irrig_limit )                                                          &
  CALL log_info("init_vegetation",                                            &
                "Irrigation is limited by water availability")

CALL log_info("init_vegetation",                                              &
              "Using can_rad_mod = " // TRIM(to_string(can_rad_mod)) //       &
              " and can_model = " // TRIM(to_string(can_model)))

SELECT CASE ( photo_model )
CASE ( photo_collatz )
  CALL log_info("init_vegetation",                                            &
                "C3 plants use the Collatz model of photosynthesis.")
CASE ( photo_farquhar )
  CALL log_info("init_vegetation",                                            &
                "C3 plants use the Farquhar model of photosynthesis.")
END SELECT
CALL log_info("init_vegetation",                                              &
              "C4 plants use the Collatz model of photosynthesis.")

SELECT CASE ( stomata_model )
CASE ( stomata_jacobs )
  CALL log_info("init_vegetation",                                            &
                "Using the original model of stomatal conductance " //        &
                "including the Jacobs closure.")
CASE ( stomata_medlyn )
  CALL log_info("init_vegetation",                                            &
                "Using the Medlyn et al. model of stomatal conductance.")
END SELECT

IF ( l_inferno ) THEN
  CALL log_info("init_vegetation",                                            &
                "Interactive fires and emissions (INFERNO) will be diagnosed")
  IF (ignition_method == 1 ) THEN
    CALL log_info("init_vegetation",                                          &
                  "Constant or ubiquitous ignitions (INFERNO)")
  ELSE IF (ignition_method == 2 ) THEN
    CALL log_info("init_vegetation",                                          &
                  "Constant human ignitions, varying lightning (INFERNO)")
  ELSE IF (ignition_method == 3 ) THEN
    CALL log_info("init_vegetation",                                          &
                  "Fully prescribed ignitions (INFERNO)")
  END IF
END IF


! Check that TRIFFID timestep (the coupling period) seems sensible
! In equilibrium mode, the coupling period should be sufficient to average
! out seasonal (and ideally interannual) variability, so a period >= 1yr
! and perhaps 5-10 yrs is recommended
! In transient mode, more frequent coupling is required - typically every
! 10 days or so
IF ( l_trif_eq .AND. triffid_period < 360 ) THEN
  CALL log_error("init_vegetation",                                           &
                 "triffid_period < 360 - in equilibrium mode a " //           &
                 "TRIFFID timestep of at least 1 year is advised")

  !   Note: We would expect "normal usage" to be that the coupling period
  !   is an integer number of years (or at least close to this if the
  !   calendar includes leap years). However, there's generally nothing
  !   intrinsically wrong with using partial years (as long as the period
  !   is >> 1 year) - it's just a slightly odd choice. In any case, we are
  !   not testing if the period is a number of years.

ELSE IF ( l_triffid .AND. triffid_period > 30 ) THEN
  CALL log_error("init_vegetation",                                           &
                 "triffid_period > 30 - in dynamic mode a TRIFFID " //        &
                 "timestep of <30 days is recommended; 10 days is " //        &
                 "often used")
END IF

RETURN
END SUBROUTINE init_vegetation

END MODULE init_vegetation_mod
#endif
