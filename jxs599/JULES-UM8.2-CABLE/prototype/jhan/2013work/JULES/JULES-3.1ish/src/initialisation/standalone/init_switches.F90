#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/init_switches.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

SUBROUTINE init_switches()

  USE io_constants, ONLY : NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE datetime_mod, ONLY : l_360_dt => l_360

  USE switches

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the model switches and checks them for consistency
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

!-----------------------------------------------------------------------------
! First, we read the switches namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_switches", "Reading JULES_SWITCHES namelist...")

  OPEN(NAMELIST_UNIT, FILE='switches.nml',                                    &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_switches",                                           &
                   "Error opening namelist file switches.nml " //             &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  READ(NAMELIST_UNIT, nml=jules_switches, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_switches",                                           &
                   "Error reading namelist JULES_SWITCHES " //                &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_switches",                                           &
                   "Error closing namelist file switches.nml " //             &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


!-----------------------------------------------------------------------------
! Check that the given combination of switches makes sense
!-----------------------------------------------------------------------------
  IF ( l_veg_compete .AND. ( .NOT. l_triffid ) )                              &
    CALL log_fatal("init_switches",                                           &
                   'Cannot use competing vegetation with l_triffid = false')

  IF ( l_aggregate .AND. (l_phenol .OR. l_triffid .OR. l_trif_eq) )           &
    CALL log_fatal("init_switches",                                           &
                   'Phenology or TRIFFID cannot be used with the ' //         &
                   'aggregated surface scheme (i.e. l_aggregate = true)')

  IF ( can_model<1 .OR. can_model>4 )                                         &
    CALL log_fatal("init_switches",                                           &
                   'can_model should be in range 1 to 4')

  IF ( can_model==4 .AND. l_aggregate )                                       &
    CALL log_fatal("init_switches",                                           &
                   'can_model=4 cannot be used with the aggregated ' //       &
                   'surface scheme (i.e. l_aggregate = true)')

  IF ( can_rad_mod < 1 .OR. can_rad_mod > 5 )                                 &
    CALL log_fatal("init_switches",                                           &
                   'can_rad_mod should be in range 1 to 5')

! If spectral albedo selected, use spectral snow albedo.
  l_snow_albedo = .FALSE.
  IF ( l_spec_albedo ) l_snow_albedo = .TRUE.

! Check that at most one of TOPMODEL and PDM have been selected.
  IF ( l_pdm .AND. l_top )                                                    &
    CALL log_fatal("init_switches",                                           &
                   'At most one of PDM and TOPMODEL can be selected')

  IF ( soilhc_method < 1 .OR. soilhc_method > 2 )                             &
    CALL log_fatal("init_switches", "soilhc_method must be 1 or 2")

! Set l_360 in datetime_mod to the same value as in switches
  l_360_dt = l_360

!------------------------------------------------------------------------------
! Set l_q10 to TRUE if not using TRIFFID.
! This is anyway done later in MICROBE, so do here too for ease.
!------------------------------------------------------------------------------
  IF ( .NOT. l_triffid ) l_q10 = .TRUE.

!-----------------------------------------------------------------------------
! Log some information about the main switches we are using
!-----------------------------------------------------------------------------
  IF ( l_aggregate )                                                          &
    CALL log_info("init_switches",                                            &
                  "Aggregate surface scheme has been selected")

  IF ( l_phenol )                                                             &
    CALL log_info("init_switches", "Phenology is on")

  IF ( l_triffid )                                                            &
    CALL log_info("init_switches", "TRIFFID is on")

  IF ( l_veg_compete )                                                        &
    CALL log_info("init_switches", "Competing vegetation is on")

  IF ( l_360 )                                                                &
    CALL log_info("init_switches", "Running with 360 day year")

  CALL log_info("init_switches",                                              &
                "Using can_rad_mod = " // TRIM(to_string(can_rad_mod)) //     &
                " and can_model = " // TRIM(to_string(can_model)))

  IF ( l_vg_soil )                                                            &
    CALL log_info("init_switches", "van Genuchten model will be used")

  IF ( l_soil_sat_down ) THEN
    CALL log_info("init_switches",                                            &
                  "l_soil_sat_down=T - excess water is pushed down")
  ELSE
    CALL log_info("init_switches",                                            &
                  "l_soil_sat_down=F - excess water is pushed up")
  ENDIF

  IF ( soilHc_method == 1 ) THEN
    CALL log_info("init_switches",                                            &
                  "soilHc_method=1:  following Cox et al (1999)")
  ELSE
    CALL log_info("init_switches",                                            &
                  "soilHc_method=2:  following simplified Johansen (1975)")
  ENDIF

  IF ( l_q10 ) THEN
    CALL log_info("init_switches",                                            &
                  "Q10 equation will be used for soil respiration")
  ELSE
     CALL log_info("init_switches",                                            &
                  "RothC equations will be used for soil respiration")
  ENDIF

  IF ( l_cable )                                                               &
    CALL log_info("init_switches", "CABLE is switched on as the LSM")

  RETURN

END SUBROUTINE init_switches
#endif
