#if !defined(UM_JULES)
!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************
 
MODULE deposition_species_mod

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Module containing species and parameters for the deposition model, and
!   code for checking their values.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in BIOGENIC FLUXES.
!
! Code Description:
!   Language: Fortran 90.
!-----------------------------------------------------------------------------
 
! Private scope by default.
PRIVATE

PUBLIC                                                                        &
! public variables
  ch4_mml, ch4_scaling, ch4dd_tundra, cuticle_o3, dd_ice_coeff,               &
  dep_species_name, diffusion_corr, diffusion_coeff,                          &
  h2dd_c, h2dd_m, h2dd_q,                                                     &
  r_tundra, r_wet_soil_o3, rsurf_std, species_name_len,                       &
! public subroutines
  check_jules_deposition_species

!-----------------------------------------------------------------------------
! Parameters.
!-----------------------------------------------------------------------------
INTEGER, PARAMETER ::                                                         &
  species_name_len = 10
    ! Length of species names.

REAL, PARAMETER ::                                                            &
  ch4_mml = 1.008e5
    ! Factor to convert methane flux (ug m-2 h-1) to dry deposition velocity
    ! (m s-1). ch4_mml = 3600 * 0.016 * 1.0E9 * 1.75E-6, where
    ! 0.016 = molar mass of methane (kg), 1.0E9 converts ug -> kg, and
    ! 1.75E-6 = assumed CH4 volumetric mixing ratio

!-----------------------------------------------------------------------------
! Variables that are set via namelists.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Parameters that are only used with a single species. These are added to
! the species-specific namelists but values are only used for a single
! species.
!-----------------------------------------------------------------------------
! Parameters for O3.
REAL ::                                                                       &
  cuticle_o3,                                                                 &
    ! Constant for calculation of cuticular resistance for ozone (s m-1).
  r_wet_soil_o3
    ! Wet soil surface resistance for ozone (s m-1).

! Parameters for CH4.
REAL ::                                                                       &
  ch4_scaling
    ! Scaling applied to CH4 soil uptake (dimensionless).
    ! Originally this was used to match the value from the IPCC TAR.

REAL ::                                                                       &
  ch4dd_tundra(4)
    ! Coefficients of cubic polynomial relating CH4 loss for tundra to
    ! temperature. N.B. Loss flux is in units of ug(CH4) m-2 s-1.

! Parameters for H2.
REAL, ALLOCATABLE ::                                                          &
  h2dd_c(:),                                                                  &
    ! Constant in quadratic function relating hydrogen deposition to soil
    ! moisture, for each surface type (s m-1).
  h2dd_m(:),                                                                  &
    ! Coefficient of first order term in quadratic function relating hydrogen
    ! deposition to soil moisture, for each surface type (s m-1).
  h2dd_q(:)
    ! Coefficient of second order term in quadratic function relating
    ! hydrogen deposition to soil moisture, for each surface type (s m-1).

!-----------------------------------------------------------------------------
! Parameters that depend on species.
!-----------------------------------------------------------------------------
REAL, ALLOCATABLE ::                                                          &
  dd_ice_coeff(:,:),                                                          &
    ! Coefficients in quadratic function relating dry deposition over ice to
    ! temperature. These are only used with selected species.
  diffusion_coeff(:),                                                         &
    ! Diffusion coefficient for each species (m2 s-1).
    ! A value < 0 is used to indicate that a species does not dry deposit -
    ! this was inherited from the UKCA code but will be reviewed in future
    ! when the parameterisation is implemented within JULES.
  diffusion_corr(:),                                                          &
    ! Diffusion correction for stomatal resistance for each species,
    ! accounting for the different diffusivities of water and other species
    ! (dimensionless). These are only used with selected species.
  rsurf_std(:,:),                                                             &
    ! Standard value of surface resistance for each surface type and tracer 
    ! species (s m-1).
  r_tundra(:)
    ! Surface resistance used in tundra region for each species (s m-1).
    ! These are only used with selected species.

CHARACTER(LEN=10), ALLOCATABLE ::                                             &
  dep_species_name(:)
    ! Names of atmospheric trace species.
  
CONTAINS

!-----------------------------------------------------------------------------

SUBROUTINE check_jules_deposition_species()

!-----------------------------------------------------------------------------
! Description:
!   Checks values from the JULES_DEPOSITION_SPECIES namelist.
!-----------------------------------------------------------------------------

USE ereport_mod, ONLY: ereport
USE jules_deposition_mod, ONLY: ndry_dep_species
USE missing_data_mod, ONLY: rmdi

IMPLICIT NONE

INTEGER ::                                                                    &
  errorstatus,                                                                &
  i  ! Loop counter.

CHARACTER(LEN=*), PARAMETER ::                                                &
  RoutineName = 'CHECK_JULES_DEPOSITION_SPECIES'   ! Name of this procedure.

!-----------------------------------------------------------------------------
!end of header

!-----------------------------------------------------------------------------
! Check that all required variables were provided by the namelist.
! Some of these variables are coded with a deafult value, but safer to check
! everything in case that ever changes.
!-----------------------------------------------------------------------------

! Check we have a species name.
IF ( ANY( dep_species_name(:) == 'unset' ) ) THEN
  errorstatus = 1   ! a fatal error
  CALL ereport( RoutineName, errorstatus,                                     &
                "Missing value for dep_species_name." )
END IF

! Check we don't have any duplicate species names.
DO i = 1, ndry_dep_species - 1
  IF ( ANY( dep_species_name(i+1:) == dep_species_name(i) ) ) THEN
    errorstatus = 1   ! a fatal error
    CALL ereport( RoutineName, errorstatus,                                   &
                  "Duplicate species names." )
  END IF
END DO

!-----------------------------------------------------------------------------
! Check variables for which we require a value for each species.
!-----------------------------------------------------------------------------
IF ( ANY( ABS( diffusion_coeff(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  errorstatus = 1   ! a fatal error
  CALL ereport( RoutineName, errorstatus,                                     &
                "Missing value for diffusion_coeff." )
END IF

IF ( ANY( ABS( rsurf_std(:,:) - rmdi ) < EPSILON(1.0) ) ) THEN
  errorstatus = 1   ! a fatal error
  CALL ereport( RoutineName, errorstatus,                                     &
                "Missing value for rsurf_std." )
END IF

!-----------------------------------------------------------------------------
! Check variables that are only required for one species.
!-----------------------------------------------------------------------------
DO i = 1, ndry_dep_species
  SELECT CASE ( dep_species_name(i) )


  CASE ( 'CH4' )
    IF ( ABS( ch4_scaling - rmdi ) < EPSILON(1.0) ) THEN
      errorstatus = 1   ! a fatal error
      CALL ereport( RoutineName, errorstatus,                                 &
                    "Missing value for ch4_scaling." )
    END IF
    IF ( ANY( ABS( ch4dd_tundra(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      errorstatus = 1   ! a fatal error
      CALL ereport( RoutineName, errorstatus,                                 &
                    "Missing value for ch4dd_tundra." )
    END IF

  CASE ( 'H2' )
    IF ( ANY( ABS( h2dd_c(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      errorstatus = 1   ! a fatal error
      CALL ereport( RoutineName, errorstatus,                                 &
                    "Missing value for h2dd_c." )
    END IF
    IF ( ANY( ABS( h2dd_m(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      errorstatus = 1   ! a fatal error
      CALL ereport( RoutineName, errorstatus,                                 &
                    "Missing value for h2dd_m." )
    END IF
    IF ( ANY( ABS( h2dd_q(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      errorstatus = 1   ! a fatal error
      CALL ereport( RoutineName, errorstatus,                                 &
                    "Missing value for h2dd_q." )
    END IF

  CASE ( 'O3' )
    IF ( ABS( cuticle_o3 - rmdi ) < EPSILON(1.0) ) THEN
      errorstatus = 1   ! a fatal error
      CALL ereport( RoutineName, errorstatus,                                 &
                    "Missing value for cuticle_o3." )
    END IF
    IF ( ABS( r_wet_soil_o3 - rmdi ) < EPSILON(1.0) ) THEN
      errorstatus = 1   ! a fatal error
      CALL ereport( RoutineName, errorstatus,                                 &
                    "Missing value for r_wet_soil_o3." )
    END IF

  END SELECT
END DO

!-----------------------------------------------------------------------------
! Check variables that are only required for certain species, but which are
! held in arrays that are sized with the number of species.
! Raise an error if a required value is not provided.
! Issue a warning if a value is not required but is provided.
! This is used to recreate the legacy UKCA code (which treats some species
! differently, requiring extra parameters) while providing some flexibility
! for possible future changes. Any changes to which parameters are used
! for each species should be reflected in updated checks here.
!----------------------------------------------------------------------------
DO i = 1, ndry_dep_species

  !--------------------------------------------------------------------------
  ! Check that dd_ice_coeff is provided when required.
  !--------------------------------------------------------------------------
  SELECT CASE ( dep_species_name(i) )

  CASE ( 'HNO3', 'HONO2', 'ISON' )
    ! A value is required for these species.
    IF ( ANY( ABS( dd_ice_coeff(i,:) - rmdi ) < EPSILON(1.0) ) ) THEN
      errorstatus = 1   ! a fatal error
      CALL ereport( RoutineName, errorstatus,                                 &
                    "Missing value for dd_ice_coeff for species=" //          &
                    TRIM(dep_species_name(i)) )
    END IF

  CASE DEFAULT
    ! A value is not required for these species.
    IF ( ANY( ABS( dd_ice_coeff(i,:) - rmdi ) > EPSILON(1.0) ) ) THEN
      errorstatus = -1   ! a warning
      CALL ereport( RoutineName, errorstatus,                                 &
                    "dd_ice_coeff is provided for a species for which " //    &
                    " it is not required. Species=" //                        &
                    TRIM(dep_species_name(i)) )
    END IF

  END SELECT

  !--------------------------------------------------------------------------
  ! Check that diffusion_corr is provided when required.
  !--------------------------------------------------------------------------
  SELECT CASE ( dep_species_name(i) )

  CASE ( 'NO2', 'O3', 'PAN', 'PPAN', 'MPAN', 'ONITU' )
    ! A value is required for these species.
    IF ( ABS( diffusion_corr(i) - rmdi ) < EPSILON(1.0) ) THEN
      errorstatus = 1   ! a fatal error
      CALL ereport( RoutineName, errorstatus,                                 &
                    "Missing value for diffusion_corr for species=" //        &
                    TRIM(dep_species_name(i)) )
    END IF

  CASE DEFAULT
    ! A value is not required for these species.
    IF ( ABS( diffusion_corr(i) - rmdi ) > EPSILON(1.0) ) THEN
      errorstatus = -1   ! a warning
      CALL ereport( RoutineName, errorstatus,                                 &
                    "diffusion_corr is provided for a species for which " //  &
                    " it is not required. Species=" //                        &
                    TRIM(dep_species_name(i)) )
    END IF

  END SELECT

  !--------------------------------------------------------------------------
  ! Check that r_tundra is provided when required.
  !--------------------------------------------------------------------------
  SELECT CASE ( dep_species_name(i) )

  CASE ( 'CO', 'NO2', 'O3', 'PAN', 'PPAN', 'MPAN', 'ONITU' )
    ! A value is required for these species.
    IF ( ABS( r_tundra(i) - rmdi ) < EPSILON(1.0) ) THEN
      errorstatus = 1   ! a fatal error
      CALL ereport( RoutineName, errorstatus,                                 &
                    "Missing value for r_tundra for species=" //              &
                    TRIM(dep_species_name(i)) )
    END IF

  CASE DEFAULT
    ! A value is not required for these species.
    IF ( ABS( r_tundra(i) - rmdi ) > EPSILON(1.0) ) THEN
      errorstatus = -1   ! a warning
      CALL ereport( RoutineName, errorstatus,                                 &
                    "r_tundra is provided for a species for which " //        &
                    " it is not required. Species=" //                        &
                    TRIM(dep_species_name(i)) )
    END IF

  END SELECT

END DO  !  species

END SUBROUTINE check_jules_deposition_species

END MODULE deposition_species_mod
#endif
