#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE init_mod
CONTAINS
SUBROUTINE init(nml_dir)

USE init_jules_sf_diags_mod,      ONLY: allocate_sf_diags,                    &
                                        set_sf_diag_switches
USE lsm_switch_mod,               ONLY: init_lsm, lsm_id, cable
USE time_varying_input_mod,       ONLY: seek_all_to_current_datetime
USE model_time_mod,               ONLY: is_spinup
USE init_grid_mod,                ONLY: init_grid
USE init_ancillaries_mod,         ONLY: init_ancillaries
USE init_params_mod,              ONLY: init_params
USE initial_conditions_mod,       ONLY: init_ic
USE spinup_mod,                   ONLY: spinup_init
USE dump_mod,                     ONLY: write_dump
USE logging_mod,                  ONLY: log_info
USE init_output_mod,              ONLY: init_output
USE init_plant_n_uptake_mod,      ONLY: init_plant_n_uptake
USE init_vegetation_mod,          ONLY: init_vegetation
USE init_hydrology_mod,           ONLY: init_hydrology
USE init_model_environment_mod,   ONLY: init_model_environment
USE init_radiation_mod,           ONLY: init_radiation
USE init_snow_mod,                ONLY: init_snow
USE init_surface_types_mod,       ONLY: init_jules_surface_types,             &
                                        init_cable_surface_types
USE init_soil_biogeochem_mod,     ONLY: init_soil_biogeochem
USE init_soil_mod,                ONLY: init_soil
USE init_surface_mod,             ONLY: init_surface
USE init_rivers_mod,              ONLY: init_rivers
USE metstats_mod,                 ONLY: l_metstats
USE metstats_init_mod,            ONLY: metstats_init
USE init_parms_mod,               ONLY: init_parms
USE init_soil_ecosse_mod,         ONLY: init_soil_ecosse
USE init_check_compatibility_mod, ONLY: init_check_compatibility
USE init_cable_mod,               ONLY: init_cable_grid, init_cable_veg,      &
                                        init_cable_soil, init_cable_met,      &
                                        init_cable_rad, init_cable_canopy,    &
                                        init_cable_soilsnow
USE init_deposition_mod,          ONLY: init_deposition
USE allocate_cable_arrays_mod,    ONLY: allocate_cable_arrays
USE jules_science_fixes_mod,      ONLY: init_science_fixes

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   The main initialisation routine - initialises the model by calling
!   specialised routines
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


!-----------------------------------------------------------------------------

! Determine what environment JULES is being run in
CALL init_model_environment(nml_dir)

! Check which LSM is to be used
CALL init_lsm(nml_dir)

! Initialise the surface types
IF ( lsm_id /= cable ) CALL init_jules_surface_types(nml_dir)
IF ( lsm_id == cable ) CALL init_cable_surface_types(nml_dir)

! Initialise the surface options
CALL init_surface(nml_dir)

! Intialise the times for the run
CALL init_time(nml_dir)

! Initialise radiation options
CALL init_radiation(nml_dir)

! Initialise hydrology options
CALL init_hydrology(nml_dir)

! Initialise soil options
CALL init_soil(nml_dir)

! Initialise vegetation options
! Must be called after init_soil to check soil tiling flags
CALL init_vegetation(nml_dir)

! Initialise vegetation N uptake model options
CALL init_plant_n_uptake(nml_dir)

! Initialise soil biogeochemistry options.
CALL init_soil_biogeochem(nml_dir)

! Initialise ECOSSE soil biogeochemistry options.
CALL init_soil_ecosse(nml_dir)

! Initialise atmospheric deposition module.
CALL init_deposition(nml_dir)

! Initialise snow options
CALL init_snow(nml_dir)

! Initialise river routing parameters, ancils and grid
CALL init_rivers(nml_dir)

! Initialise two-tile urban schemes (not ancillary data) and must be called
! after init_radiation
CALL init_urban(nml_dir)

! Initialise science fixes
CALL init_science_fixes(nml_dir)

! Initialise the input, model and output grids **also allocates arrays**
CALL init_grid(nml_dir)

! Initialise the model ancils
CALL init_ancillaries(nml_dir)

! Initialise model parameters
CALL init_params(nml_dir)

! Initialise fire module
CALL init_fire(nml_dir)

! Initialise the metstats module (used by fire)
IF (l_metstats) CALL metstats_init()

! Initialise meteorological forcing
CALL init_drive(nml_dir)

! Initialise IMOGEN
CALL init_imogen(nml_dir)

! Initialise other prescribed data
CALL init_prescribed_data(nml_dir)

! Initialise the model prognostics
CALL init_ic(nml_dir)

! Initialise output
CALL init_output(nml_dir)

!-----------------------------------------------------------------------------
! Other initialisation that does not depend on further user input.
!-----------------------------------------------------------------------------

! Check that the enabled schemes are compatible
CALL init_check_compatibility()

! Set sf_diag switches that are required by the science configuration (not by
! the choice of diagnostics).
CALL set_sf_diag_switches

! Allocate for coupled model diagnostics (sf_diag).
CALL allocate_sf_diags

! Further initialisation of variables.
CALL init_vars_tmp()

! Set index arrays and initialise other variables.
CALL init_parms()

! Seek the input files to the start of the run
CALL seek_all_to_current_datetime()

! If the LSM is CABLE, initialise CABLE variables
IF ( lsm_id == cable ) THEN
  CALL allocate_cable_arrays()
  CALL init_cable_grid()
  CALL init_cable_veg()
!CABLE_LSM:hac::WE DONT NEED THESE 
!  CALL init_cable_soil()
!  CALL init_cable_met()
!  CALL init_cable_rad()
!  CALL init_cable_canopy()
!  CALL init_cable_soilsnow()
END IF

! Save initial state if spinning up. Arrays are allocated here.
IF ( is_spinup ) CALL spinup_init()

! Write an initial dump
CALL write_dump()

CALL log_info("init", "Initialisation is complete")

RETURN

END SUBROUTINE init
END MODULE init_mod
#endif
