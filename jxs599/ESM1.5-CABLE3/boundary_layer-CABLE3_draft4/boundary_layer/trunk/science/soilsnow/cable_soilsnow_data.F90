MODULE cable_soil_snow_data_mod

  USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
       veg_parameter_type, canopy_type, met_type,        &
       balances_type, r_2, ms, mp
!distribute these per sbr
USE cable_phys_constants_mod, ONLY : CTFRZ => TFRZ
USE cable_phys_constants_mod, ONLY : CHL => HL
USE cable_phys_constants_mod, ONLY : CHLF => HLF
USE cable_phys_constants_mod, ONLY : Cdensity_liq => density_liq
USE cable_phys_constants_mod, ONLY : Ccgsnow => cgsnow
USE cable_phys_constants_mod, ONLY : Ccswat => cswat
USE cable_phys_constants_mod, ONLY : Ccsice => csice
USE cable_phys_constants_mod, ONLY : Ccs_rho_wat => cs_rho_wat
USE cable_phys_constants_mod, ONLY : Ccs_rho_ice => cs_rho_ice

  USE cable_common_module, ONLY: cable_user,snow_ccnsw,snmin,&
       max_ssdn,max_sconds,frozen_limit,&
       max_glacier_snowd

   IMPLICIT NONE

    REAL, PARAMETER ::                                                          &
      cgsnow = 2090.0,     & ! specific heat capacity for snow
      csice = 2.100e3,     & ! specific heat capacity for ice
      cswat = 4.218e3,     & ! specific heat capacity for water
      rhowat = 1000.0 !,     & ! density of water
      !snmin = 1.,          & ! for 3-layer;
      !max_ssdn = 750.0,    & !
      !max_sconds = 2.51,   & !
      !frozen_limit = 0.85    ! EAK Feb2011 (could be 0.95)
   
!   REAL :: cp    ! specific heat capacity for air
!   
!   !jhan:make parameter
!   REAL :: max_glacier_snowd
!       max_glacier_snowd

END MODULE cable_soil_snow_data_mod

