
MODULE cbm_module
   
   USE canopy_module
   USE albedo_module
  
   IMPLICIT NONE
  
   PRIVATE
   PUBLIC cbm 

CONTAINS

   SUBROUTINE cbm( dels, air, bgc, canopy, met,                                &
                   bal, rad, rough, soil,                                      &
                   ssoil, sum_flux, veg )
    
    USE define_dimensions
    USE cable_common_module
    USE carbon_module
    USE soil_snow_module
    USE define_types
    USE physical_constants
    USE roughness_module
    USE radiation_module
    USE air_module
    USE cable_diag_module, only : cable_stat
    USE casadimension,     only : icycle ! used in casa_cnp

   ! CABLE model variables
   TYPE (air_type),       INTENT(INOUT) :: air
   TYPE (bgc_pool_type),  INTENT(INOUT) :: bgc
   TYPE (canopy_type),    INTENT(INOUT) :: canopy
   TYPE (met_type),       INTENT(INOUT) :: met
   TYPE (balances_type),  INTENT(INOUT) :: bal
   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (soil_snow_type), INTENT(INOUT) :: ssoil
   TYPE (sum_flux_type),  INTENT(INOUT) :: sum_flux
    
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil 
   TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg  

   REAL, INTENT(IN)               :: dels ! time setp size (s)
    
   INTEGER :: k,kk,j  

   IF( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         CALL cable_stat('cbm')

   !jhan:do this in iniotialization   
   ! Fix in-canopy turbulence scheme globally:
   !veg%meth = 1

   IF( cable_runtime%um ) THEN
      
      cable_runtime%um_radiation = .FALSE.
      
      IF( cable_runtime%um_explicit ) THEN
         CALL ruff_resist(veg, rough, ssoil, soil, met, canopy)
         met%tk = met%tk + grav/capp*(rough%zref_tq + 0.9*rough%z0m)
         !jhan:consult group b4 rm
         !met%tc = met%tk - tfrz
      ENDIF
      
      CALL define_air (met, air)
   
   ELSE
      call ruff_resist(veg, rough, ssoil, soil, met, canopy)
   ENDIF


   CALL init_radiation(met,rad,veg, canopy) ! need to be called at every dt

   IF( cable_runtime%um ) THEN
      
      IF( cable_runtime%um_explicit ) THEN
         CALL surface_albedo(ssoil, veg, met, rad, soil, canopy)
      ENDIF
   
   ELSE
      CALL surface_albedo(ssoil, veg, met, rad, soil, canopy)
   ENDIf
    
   ! Calculate canopy variables:
   CALL define_canopy(bal,rad,rough,air,met,dels,ssoil,soil,veg, canopy)

   ssoil%otss_0 = ssoil%otss
   ssoil%otss = ssoil%tss

   IF( cable_runtime%um ) THEN
      
      ssoil%owetfac = ssoil%wetfac
      IF( cable_runtime%um_hydrology ) THEN
         CALL soil_snow(dels, soil, ssoil, canopy, met, bal,veg)
      ENDIF

   ELSE
      call soil_snow(dels, soil, ssoil, canopy, met, bal,veg)
   ENDIF

   ssoil%deltss = ssoil%tss-ssoil%otss

   canopy%fhs = canopy%fhs + ( ssoil%tss-ssoil%otss ) * ssoil%dfh_dtg
   
   canopy%fhs_cor = canopy%fhs_cor + ( ssoil%tss-ssoil%otss ) * ssoil%dfh_dtg
   
    canopy%fh = canopy%fhv + canopy%fhs

   canopy%fes = canopy%fes + ( ssoil%tss-ssoil%otss ) *                        &
                ( ssoil%cls * ssoil%dfe_ddq * ssoil%ddq_dtg )
   
   canopy%fes_cor = canopy%fes_cor + ( ssoil%tss-ssoil%otss ) *                &
                    ( ssoil%cls * ssoil%dfe_ddq * ssoil%ddq_dtg )

   ! need to adjust fe after soilsnow
   canopy%fev  = canopy%fevc + canopy%fevw
  
   ! Calculate total latent heat flux:
   canopy%fe = canopy%fev + canopy%fes

   ! Calculate net radiation absorbed by soil + veg
   canopy%rnet = canopy%fns + canopy%fnv

   ! Calculate radiative/skin temperature:
   rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 +                             &
              rad%transd * ssoil%tss**4 )**0.25

   ! rml 17/1/11 move all plant resp and soil resp calculations here            
   ! from canopy. in UM only call on implicit step.
   ! put old and new soil resp calculations into soilcarb subroutine
   ! make new plantcarb subroutine
   IF (.not.cable_runtime%um_explicit .AND. icycle == 0) THEN

      !calculate canopy%frp
      CALL plantcarb(veg,bgc,met,canopy)
     
      !calculate canopy%frs
      CALL soilcarb(soil, ssoil, veg, bgc, met, canopy)

      CALL carbon_pl(dels, soil, ssoil, veg, canopy, bgc)

      canopy%fnpp = -1.0* canopy%fpn - canopy%frp
      canopy%fnee = canopy%fpn + canopy%frs + canopy%frp

   ENDIF

  
END SUBROUTINE cbm

END MODULE cbm_module


