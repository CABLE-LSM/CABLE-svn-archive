MODULE cable_wetleaf_module

  IMPLICIT NONE

  PUBLIC wetleaf
  PRIVATE

CONTAINS
SUBROUTINE wetLeaf( dels, rad, rough, air, met, veg, canopy, cansat, tlfy,     &
                    gbhu, gbhf, ghwet, &
                    mp, CLAI_thresh, CCAPP, CRmair, & 
                    canopy_fevw,canopy_fevw_pot, canopy_fhvw )

USE cable_def_types_mod, ONLY : radiation_type, roughness_type, air_type,     &  
met_type, canopy_type, veg_parameter_type, r_2
   
   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (air_type),       INTENT(INOUT) :: air
   TYPE (met_type),       INTENT(INOUT) :: met
   TYPE (canopy_type),    INTENT(INOUT) :: canopy
   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg

INTEGER, INTENT(IN) :: mp
REAL, INTENT(IN) :: CLAI_thresh, CCAPP, CRmair
REAL, INTENT(INOUT) :: canopy_fevw(mp)
REAL, INTENT(INOUT) :: canopy_fevw_pot(mp)
REAL, INTENT(INOUT) :: canopy_fhvw(mp)

   REAL,INTENT(IN), DIMENSION(:) ::                                            &
      tlfy,          & ! leaf temp (K) - assCUMINg the temperature of 
                       ! wet leaf is equal that of dry leaf ="tlfy"
      cansat           ! max canopy intercept. (mm)

   REAL(r_2), INTENT(IN), DIMENSION(:,:) ::                                    &
      gbhu,          & ! forcedConvectionBndryLayerCond
      gbhf             ! freeConvectionBndryLayerCond

   REAL(r_2), INTENT(OUT), DIMENSION(:) ::                                     &
      ghwet            ! cond for heat for a wet canopy

   REAL, INTENT(IN)     :: dels ! integration time step (s)

   ! local variables  
   REAL, DIMENSION(mp) ::                                                      &
     ccfevw,        & ! limitation term for
     gwwet,         & ! cond for water for a wet canopy
     ghrwet           ! wet canopy cond: heat & thermal rad
   
   !i sums, terms of convenience/readability
   REAL, DIMENSION(mp) ::                                                      &
     sum_gbh, sum_rad_rniso, sum_rad_gradis, xx1

   INTEGER :: j
   
   ! END header

   ghwet = 1.0e-3
   gwwet = 1.0e-3
   ghrwet= 1.0e-3

   canopy_fevw = 0.0
   canopy_fhvw = 0.0

   sum_gbh = SUM((gbhu+gbhf),2)
   sum_rad_rniso = SUM(rad%rniso,2)
   sum_rad_gradis = SUM(rad%gradis,2)

   DO j=1,mp

      IF(canopy%vlaiw(j) > CLAI_THRESH) THEN

         ! VEG SENSIBLE & LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
         ! calculate total thermal resistance, rthv in s/m
         ghwet(j) = 2.0   * sum_gbh(j) 
         gwwet(j) = 1.075 * sum_gbh(j) 
         ghrwet(j) = sum_rad_gradis(j) + ghwet(j)
         
         ! Calculate fraction of canopy which is wet:
         canopy%fwet(j) = MAX( 0.0, MIN( 1.0,                                  &
                          0.8 * canopy%cansto(j) / MAX( cansat(j), 0.01 ) ) )
         
         ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
         ! to avoid excessive evaporation:
         ccfevw(j) = MIN(canopy%cansto(j) * air%rlam(j) / dels, &
                       2.0 / (1440.0 / (dels/60.0)) * air%rlam(j) )
   
         canopy_fevw(j) = MIN( canopy%fwet(j) * ( air%dsatdk(j) *              &
                         ( sum_rad_rniso(j)- CCAPP*Crmair*( met%tvair(j)     &
                         - met%tk(j) ) * sum_rad_gradis(j) )                   &
                         + CCAPP * Crmair * met%dva(j) * ghrwet(j) )         &
                         / ( air%dsatdk(j)+air%psyc(j)*ghrwet(j) / gwwet(j) )  &
                         , ccfevw(j) )

         canopy%fevw_pot(j) = ( air%dsatdk(j)* (sum_rad_rniso(j) -             &
                              CCAPP * Crmair * ( met%tvair(j) - met%tk(j) )  &
                              *sum_rad_gradis(j) )                             &
                              + CCAPP * Crmair * met%dva(j) * ghrwet(j))     &
                              / (air%dsatdk(j)+air%psyc(j)*ghrwet(j)/gwwet(j) )
          
         ! calculate sens heat from wet canopy:
         canopy_fhvw(j) = canopy%fwet(j) * ( sum_rad_rniso(j) -CCAPP * Crmair&
                          * ( tlfy(j) - met%tk(j) ) * sum_rad_gradis(j) )      &
                           - canopy_fevw(j)

      ENDIF
       
   ENDDO 
           
END SUBROUTINE wetLeaf

END MODULE cable_wetleaf_module
