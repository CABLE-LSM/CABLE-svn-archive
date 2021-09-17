MODULE cable_cbm_module
   
   USE cable_canopy_module
   USE cable_albedo_module
  
   IMPLICIT NONE
  
   PRIVATE
   PUBLIC cbm 

CONTAINS

   SUBROUTINE cbm( dels, air, bgc, canopy, met,                                &
                   bal, rad, rough, soil,                                      &
                   ssnow, sum_flux, veg,                                       &
                   xk, c1, rhoch )
    
   USE cable_common_module
   USE cable_carbon_module
   USE cable_soil_snow_module
   USE cable_def_types_mod
   USE cable_roughness_module
   USE cable_radiation_module
   USE cable_air_module
!CBL3 
!CBL3 USE cbl_albedo_mod, ONLY: albedo
USE cbl_masks_mod, ONLY: fveg_mask,  fsunlit_mask,  fsunlit_veg_mask
USE cbl_masks_mod, ONLY: veg_mask,  sunlit_mask,  sunlit_veg_mask
!jhan:pass these !data
USE cable_other_constants_mod, ONLY: Ccoszen_tols => coszen_tols
USE cable_other_constants_mod,  ONLY : Crad_thresh => rad_thresh
USE cable_other_constants_mod, ONLY: clai_thresh => lai_thresh
USE cable_other_constants_mod, ONLY: cgauss_w => gauss_w
USE cable_math_constants_mod,  ONLY: cpi => pi
USE cable_math_constants_mod,  ONLY: cpi180 => pi180
USE cable_climate_type_mod, ONLY : climate_cbl

   USE cable_data_module, ONLY : icbm_type, point2constants 
   
   !ptrs to local constants 
   TYPE( icbm_type ) :: C
   ! CABLE model variables
   TYPE (air_type),       INTENT(INOUT) :: air
   TYPE (bgc_pool_type),  INTENT(INOUT) :: bgc
   TYPE (canopy_type),    INTENT(INOUT) :: canopy
   TYPE (met_type),       INTENT(INOUT) :: met
   TYPE (balances_type),  INTENT(INOUT) :: bal
   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (sum_flux_type),  INTENT(INOUT) :: sum_flux
    
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil 
   TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg  

   REAL, INTENT(IN)               :: dels ! time setp size (s)
    
   INTEGER :: k,kk,j  

CHARACTER(LEN=*), PARAMETER :: subr_name = "cbl_model_driver"
LOGICAL :: jls_standalone= .TRUE.
LOGICAL :: jls_radiation= .FALSE.
LOGICAL :: cbl_standalone = .FALSE.    

!co-efficients usoughout init_radiation ` called from _albedo as well
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: xk(mp,nrb)

   ! assign local ptrs to constants defined in cable_data_module
   CALL point2constants(C)    

      cable_runtime%um_radiation = .FALSE.
      
      IF( cable_runtime%um_explicit ) THEN
        CALL ruff_resist( veg, rough, ssnow, canopy, veg%vlai, veg%hc, canopy%vlaiw )
         met%tk = met%tk + C%grav/C%capp*(rough%zref_tq + 0.9*rough%z0m)
      ENDIF
      
      CALL define_air (met, air)
   
   CALL init_radiation(met,rad,veg, canopy) ! need to be called at every dt
!d1!CALL init_radiation( met, rad, veg, canopy,                                     &
!d1!                     rad%extkb, rad%extkd,                                     &
!d1!                     !ExtCoeff_beam, ExtCoeff_dif,                             &
!d1!                     rad%extkbm, rad%extkdm, Rad%Fbeam,                        &
!d1!                     !EffExtCoeff_beam, EffExtCoeff_dif, RadFbeam,             &
!d1!                     c1, rhoch, xk,                                            &
!d1!                     mp,nrb,                                                   &
!d1!                     Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180,         &
!d1!                     cbl_standalone, jls_standalone, jls_radiation,            &
!d1!                     subr_name,                                                &
!d1!                     veg_mask, sunlit_mask, sunlit_veg_mask,                   &
!d1!                     veg%Xfang, veg%taul, veg%refl,                            &
!d1!                     !VegXfang, VegTaul, VegRefl                               &
!d1!                     met%coszen, int(met%DoY), met%fsd,                        &
!d1!                     !coszen, metDoY, SW_down,                                 &
!d1!                     canopy%vlaiw  ) !reducedLAIdue2snow 
 
      IF( cable_runtime%um_explicit ) THEN

         CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)
!d1!CALL surface_albedo(ssnow, veg, met, rad, soil, canopy, &
!d1! ssnow%AlbSoilsn, soil%AlbSoil,                                 &
!d1!             !AlbSnow, AlbSoil,              
!d1!             mp, nrb,                                                       &
!d1!             jls_radiation,                                                 &
!d1!             veg_mask, sunlit_mask, sunlit_veg_mask,                        &  
!d1!             Ccoszen_tols, cgauss_w,                                        & 
!d1!             veg%iveg, soil%isoilm, veg%refl, veg%taul,                    & 
!d1!             !surface_type, VegRefl, VegTaul,
!d1!             met%tk, met%coszen, canopy%vlaiw,                              &
!d1!             !metTk, coszen, reducedLAIdue2snow,
!d1!             ssnow%snowd, ssnow%osnowd, ssnow%isflag,                       & 
!d1!             !SnowDepth, SnowODepth, SnowFlag_3L, 
!d1!             ssnow%ssdnn, ssnow%tgg(:,1), ssnow%tggsn(:,1), ssnow%snage,                      & 
!d1!             !SnowDensity, SoilTemp, SnowAge, 
!d1!             xk, c1, rhoch,                                                 & 
!d1!             rad%fbeam, rad%albedo,                                         &
!d1!             !RadFbeam, RadAlbedo,
!d1!             rad%extkd, rad%extkb,                                          & 
!d1!             !ExtCoeff_dif, ExtCoeff_beam,
!d1!             rad%extkdm, rad%extkbm,                                        & 
!d1!             !EffExtCoeff_dif, EffExtCoeff_beam,                
!d1!             rad%rhocdf, rad%rhocbm,                                        &
!d1!             !CanopyRefl_dif,CanopyRefl_beam,
!d1!             rad%cexpkdm, rad%cexpkbm,                                      & 
!d1!             !CanopyTransmit_dif, CanopyTransmit_beam, 
!d1!             rad%reffdf, rad%reffbm                                        &
!d1!           ) !EffSurfRefl_dif, EffSurfRefl_beam 

      ENDIF
   
   ! Calculate canopy variables:
   CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy, climate_cbl, sunlit_veg_mask, canopy%vlaiw )

   ssnow%otss_0 = ssnow%otss
   ssnow%otss = ssnow%tss
   ! RML moved out of following IF after discussion with Eva
   ssnow%owetfac = ssnow%wetfac

   IF( cable_runtime%um ) THEN
      
     IF( cable_runtime%um_implicit ) THEN
         CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
      ENDIF

   ELSE
      call soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
   ENDIF

   ssnow%deltss = ssnow%tss-ssnow%otss
   ! correction required for energy balance in online simulations
   IF( cable_runtime%um ) THEN
   
      canopy%fhs = canopy%fhs + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
      
      canopy%fhs_cor = canopy%fhs_cor + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
      
      canopy%fh = canopy%fhv + canopy%fhs

   canopy%fes = canopy%fes + ( ssnow%tss-ssnow%otss ) *                        &
                ( ssnow%dfe_ddq * ssnow%ddq_dtg )
                !( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )
   
   canopy%fes_cor = canopy%fes_cor + ( ssnow%tss-ssnow%otss ) *                &
                    ( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )

   ENDIF

   ! need to adjust fe after soilsnow
   canopy%fev  = canopy%fevc + canopy%fevw
  
   ! Calculate total latent heat flux:
   canopy%fe = canopy%fev + canopy%fes

   ! Calculate net radiation absorbed by soil + veg
   canopy%rnet = canopy%fns + canopy%fnv

   ! Calculate radiative/skin temperature:
   rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 +                             &
              rad%transd * ssnow%tss**4 )**0.25

END SUBROUTINE cbm

END MODULE cable_cbm_module


