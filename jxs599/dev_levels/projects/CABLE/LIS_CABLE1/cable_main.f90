! cable_main.f90
!
! Source file containing main routine execution order for CABLE
!
! Development by Ying-Ping Wang, Eva Kowalczyk, Ray Leuning
! Gab Abramowitz, Martin Dix, Harvey Davies, Mike Raupach, 
!
! bugs to bernard.pak@csiro.au
!
! This file contains modules:
!   main_cable_module 
! The subroutines included are:
!   cable_main
!
! Most user-defined types (e.g. met%tk) are defined in cable_types module
! in cable_variables.f90

MODULE main_cable_module
  USE canopy_vh_module
  USE canopy_module
  USE carbon_module
  USE soil_snow_module, ONLY: soil_snow
  USE cable_types
  USE cable_dimensions, ONLY: i_d,r_1
  USE physical_constants, ONLY: emsoil, sboltz
  USE roughness_module
  USE radiation_module
  USE albedo_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC cable_main 
CONTAINS
  SUBROUTINE cable_main(ktau, kstart, kend, dels, air, bgc, &
       canopy, met, bal, rad, rough, soil, ssoil, sum_flux, &
       veg, nvegt, nsoilt, model_structure,prin)
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number
    INTEGER(i_d), INTENT(IN)            :: kstart ! starting value of ktau
    INTEGER(i_d), INTENT(IN)            :: kend ! total # timesteps in run
    REAL(r_1), INTENT(IN)               :: dels ! time setp size (s)
    TYPE (air_type), INTENT(INOUT)      :: air
    TYPE (bgc_pool_type), INTENT(INOUT) :: bgc
    TYPE (canopy_type), INTENT(INOUT)   :: canopy
    TYPE (met_type), INTENT(INOUT)      :: met
    TYPE (balances_type), INTENT(INOUT) :: bal
    TYPE (radiation_type), INTENT(INOUT)      :: rad
    TYPE (roughness_type), INTENT(INOUT)      :: rough
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE (sum_flux_type), INTENT(INOUT)       :: sum_flux
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE (model_structure_type), INTENT(IN)  :: model_structure
    INTEGER(i_d), INTENT(IN)            :: nvegt  ! Number of vegetation types
    INTEGER(i_d), INTENT(IN)            :: nsoilt ! Number of soil types
    
    LOGICAL :: prin

    ! Fix in-canopy turbulence scheme globally:
    veg%meth = 1

    CALL ruff_resist(veg, rough, ssoil, canopy)
    CALL init_radiation(met, rad, veg, canopy)
    CALL surface_albedo(ssoil, veg, met, rad, canopy)

    IF(ANY(model_structure%canopy=='canopy_vh')) &
       CALL define_canopy_vh(ktau,bal,rad,rough,air,met,dels,ssoil,soil,veg, &
            & bgc,canopy,model_structure)
    IF(ANY((model_structure%canopy=='default').OR.(model_structure%canopy=='hawkesbury'))) &
!         CALL define_canopy(ktau,bal,rad,rough,air,met,dels,ssoil,soil,veg, &
!            & bgc,canopy,model_structure)
         CALL define_canopy(ktau,bal,rad,rough,air,met,dels,ssoil,soil,veg, &
            & bgc,canopy,model_structure, prin)


    ! Calculate soil and snow variables:
    IF(ANY(model_structure%soil=='soilsnow')) &
       CALL soil_snow(dels, ktau, soil, ssoil, veg, canopy, met, rad,air,prin)
    IF(ANY(model_structure%soil=='sli')) &
       CALL sli_main(ktau,dels,veg,soil,ssoil,met,canopy,air)

    ! Call carbon routines:
    CALL sumcflux(ktau, kstart, kend, dels, bgc, canopy,  &
         & soil, ssoil, sum_flux, veg, met, nvegt, nsoilt)

  END SUBROUTINE cable_main

END MODULE main_cable_module
