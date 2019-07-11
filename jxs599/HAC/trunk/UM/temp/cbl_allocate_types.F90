module cbl_allocate_types_mod
   
USE cable_def_types_mod, ONLY :  air_type,  bgc_pool_type,  met_type,          &
                        balances_type, radiation_type, roughness_type,         &
                        soil_parameter_type, soil_snow_type, sum_flux_type,    &
                        veg_parameter_type, canopy_type

IMPLICIT NONE

TYPE(air_type), SAVE             :: air
TYPE(bgc_pool_type), SAVE        :: bgc
TYPE(met_type), SAVE             :: met
TYPE(balances_type), SAVE        :: bal
TYPE(radiation_type),  SAVE       :: rad
TYPE(roughness_type),  SAVE       :: rough
TYPE(soil_parameter_type), SAVE  :: soil       ! soil parameters
TYPE(soil_snow_type), SAVE       :: ssnow
TYPE(sum_flux_type), SAVE        :: sum_flux
TYPE(veg_parameter_type), SAVE   :: veg        ! vegetation parameters
TYPE(canopy_type), SAVE          :: canopy

Public :: alloc_cbl_types

contains

SUBROUTINE alloc_cbl_types( mp )

USE cable_def_types_mod, ONLY : alloc_cbm_var

implicit none

integer :: mp
logical, save :: first_call =.true.

if ( first_call  ) then
  CALL alloc_cbm_var(air, mp)
  CALL alloc_cbm_var(canopy, mp)
  CALL alloc_cbm_var(met, mp)
  CALL alloc_cbm_var(bal, mp)
  CALL alloc_cbm_var(rad, mp)
  CALL alloc_cbm_var(rough, mp)
  CALL alloc_cbm_var(soil, mp)
  CALL alloc_cbm_var(ssnow, mp)
  CALL alloc_cbm_var(sum_flux, mp)
  CALL alloc_cbm_var(veg, mp)
  CALL alloc_cbm_var(bgc, mp)
endif
first_call =.false.
END SUBROUTINE alloc_cbl_types

End module cbl_allocate_types_mod
