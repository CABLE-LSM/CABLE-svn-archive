
program cable_main

   use cable_driver_module
   use cable_def_types_mod
   use cable_Duplicate_types_mod, ONLY : TdupVars
   use cable_Buffer_types_mod, only : TBufVars
   USE casavariable,        ONLY: casafile, casa_biome, casa_pool, casa_flux,  &
                                  casa_met, casa_balance
   USE phenvariable,        ONLY: phen_variable
   
   ! CABLE variables
   TYPE (met_type), save       :: met     ! met input variables
   TYPE (air_type), save       :: air     ! air property variables
   TYPE (canopy_type), save    :: canopy  ! vegetation variables
   TYPE (radiation_type), save :: rad     ! radiation variables
   TYPE (roughness_type), save :: rough   ! roughness varibles
   TYPE (balances_type), save  :: bal     ! energy and water balance variables
   TYPE (soil_snow_type), save :: ssnow   ! soil and snow variables
   
   TYPE (sum_flux_type), save  :: sum_flux ! cumulative flux variables
   TYPE (bgc_pool_type), save  :: bgc  ! carbon pool variables
 
   ! CABLE parameters
   TYPE (soil_parameter_type), save :: soil ! soil parameters	
   TYPE (veg_parameter_type), save  :: veg  ! vegetation parameters	 
 
    ! CASA-CNP variables 
   TYPE (casa_biome), save     :: casabiome
   TYPE (casa_pool), save      :: casapool
   TYPE (casa_flux), save      :: casaflux
   TYPE (casa_met), save       :: casamet
   TYPE (casa_balance), save   :: casabal
   TYPE (phen_variable), save  :: phen 
   
   type (TdupVars), save :: dup
   type (TbufVars), save :: buf

   call cable_offline_driver( met, air, canopy, rad, rough, &
                          ssnow, soil, veg, bal,        &
                          sum_flux, bgc, &
                          casabiome, casapool, casaflux, &
                          casamet, casabal, phen, dup, buf ) 

end program cable_main
