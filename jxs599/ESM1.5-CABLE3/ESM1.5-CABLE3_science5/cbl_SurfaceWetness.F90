MODULE cbl_SurfaceWetness_module

IMPLICIT NONE

PUBLIC :: Surf_wetness_fact
PRIVATE

CONTAINS

SUBROUTINE Surf_wetness_fact( cansat, canopy, ssnow,veg, met, soil, dels )

    USE cable_common_module
    USE cable_def_types_mod
! physical constants
USE cable_phys_constants_mod, ONLY : CTFRZ   => TFRZ
    !H!USE cable_gw_hydro_module, ONLY : calc_srf_wet_fraction

    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
    TYPE (soil_snow_type), INTENT(inout):: ssnow
    TYPE (soil_parameter_type), INTENT(inout)   :: soil
    TYPE (canopy_type), INTENT(INOUT)   :: canopy
    TYPE (met_type), INTENT(INOUT)   :: met

    REAL, INTENT(IN) :: dels ! integration time setp (s)

    REAL,INTENT(IN), DIMENSION(:) :: cansat ! max canopy intercept. (mm)

    !local variables
    REAL, DIMENSION(mp)  :: lower_limit, upper_limit,ftemp

    INTEGER :: j, i



   ! Rainfall variable is limited so canopy interception is limited,
   ! used to stabilise latent fluxes.
   ! to avoid excessive direct canopy evaporation (EK nov2007, snow scheme)
   upper_limit = 4.0 * MIN(dels,1800.0) / (60.0 * 1440.0 ) 
   ftemp =MIN(met%precip-met%precip_sn, upper_limit )
   ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
   lower_limit = cansat - canopy%cansto
   upper_limit = max(lower_limit, 0.0) 
   canopy%wcint = MERGE( MIN( upper_limit, ftemp ), 0.0,                       &
                  ftemp > 0.0  .AND. met%tk > Ctfrz)  !EAK, 09/10

   ! Define canopy throughfall (100% of precip if temp < 0C, see above):
   canopy%through = met%precip_sn + MIN( met%precip - met%precip_sn ,          &
                    MAX( 0.0, met%precip - met%precip_sn - canopy%wcint) ) 

   ! Add canopy interception to canopy storage term:
   canopy%cansto = canopy%cansto + canopy%wcint

   ! Calculate fraction of canopy which is wet:
   canopy%fwet   = MAX( 0.0, MIN( 0.9, 0.8 * canopy%cansto /                   &
                   MAX( cansat, 0.01 ) ) )

   ssnow%wetfac = MAX( 1.e-6, MIN( 1.0,                                        &
                  ( REAL (ssnow%wb(:,1) ) - soil%swilt/ 2.0 )                  &
                  / ( soil%sfc - soil%swilt/2.0 ) ) )
  
   DO j=1,mp
   
      IF( ssnow%wbice(j,1) > 0. )                                              &
         ssnow%wetfac(j) = ssnow%wetfac(j) * MAX( 0.5, 1. - MIN( 0.2,          &
                           ( ssnow%wbice(j,1) / ssnow%wb(j,1) )**2 ) )

      IF( ssnow%snowd(j) > 0.1) ssnow%wetfac(j) = 0.9
      
      IF ( veg%iveg(j) == 16 .and. met%tk(j) >= Ctfrz + 5. )                  &
         ssnow%wetfac(j) = 1.0 ! lakes: hard-wired number to be removed
      
      IF( veg%iveg(j) == 16 .and. met%tk(j) < Ctfrz + 5. )                    &
         ssnow%wetfac(j) = 0.7 ! lakes: hard-wired number to be removed

   ENDDO 
      

   ! owetfac introduced to reduce sharp changes in dry regions,
   ! especially in offline runs in which there may be discrepancies b/n
   ! timing of precip and temperature change (EAK apr2009)
   ssnow%wetfac = 0.5*(ssnow%wetfac + ssnow%owetfac)




  END SUBROUTINE Surf_wetness_fact


END MODULE cbl_SurfaceWetness_module
