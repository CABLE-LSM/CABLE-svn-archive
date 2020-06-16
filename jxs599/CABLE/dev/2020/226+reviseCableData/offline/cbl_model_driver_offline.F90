!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Calls CABLE routines including define_air, surface_albedo,
!          define_canopy, soilsnow, carbon
!          Note that cbm is called once per timestep in the offline case but
!          twice per timestep in the ACCESS case. Not all parts of cbm
!          are executed in each of the ACCESS calls.
!
! Called from: cable_driver for offline version
!              cable_explicit_driver, cable_implicit_driver for ACCESS
!
! Contact: Yingping.Wang@csiro.au
!
! History: Calling sequence changes for ACCESS compared to v1.4b
!
!          REV_CORR package of fixes for the sensitivity/correction terms
!
! ==============================================================================

MODULE cable_cbm_module

  IMPLICIT NONE

  PRIVATE
  PUBLIC cbm

CONTAINS

  SUBROUTINE cbm( ktau,dels, air, bgc, canopy, met,                                &
       bal, rad, rough, soil,                                      &
       ssnow, sum_flux, veg, climate )

    USE cable_common_module
    USE cable_carbon_module
    USE cbl_soil_snow_main_module, ONLY : soil_snow
    USE cable_def_types_mod
    USE cable_roughness_module, ONLY : ruff_resist
    USE cbl_init_radiation_module, ONLY : init_radiation
    USE cable_air_module, ONLY : define_air
    USE casadimension,     ONLY : icycle ! used in casa_cnp
    USE cable_data_module, ONLY : icbm_type, point2constants
    !mrd561
    USE cable_gw_hydro_module, ONLY : sli_hydrology,&
         soil_snow_gw
    USE cable_canopy_module, ONLY : define_canopy
    !USE cable_albedo_module, ONLY : surface_albedo
    USE cbl_albedo_mod, ONLY : albedo
    USE sli_main_mod, ONLY : sli_main
!data !jhan:pass these
USE cable_other_constants_mod, ONLY : CLAI_THRESH => lai_thresh
USE cable_other_constants_mod,  ONLY : Ccoszen_tols => coszen_tols
USE cable_other_constants_mod, ONLY : CGAUSS_W => gauss_w
USE cable_math_constants_mod, ONLY : CPI => pi
USE cable_math_constants_mod, ONLY : CPI180 => pi180
use cbl_masks_mod, ONLY :  fveg_mask,  fsunlit_mask,  fsunlit_veg_mask

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
    TYPE (climate_type), INTENT(IN)      :: climate

    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg

    REAL, INTENT(IN)               :: dels ! time setp size (s)
    INTEGER, INTENT(IN) :: ktau
    INTEGER :: k,kk,j
    LOGICAL, SAVE :: first_call = .TRUE.
character(len=*), parameter :: subr_name = "cbm"
LOGICAL :: cbl_standalone= .true.
LOGICAL :: jls_standalone= .false.
LOGICAL :: jls_radiation= .false.

!make local to rad_driver and also again in cbl_model_driver
!CABLE variables to keep for all CABLE pathways across the timestep 
real :: reducedLAIdue2snow(mp)

!masks
logical :: veg_mask(mp),  sunlit_mask(mp),  sunlit_veg_mask(mp) 
logical :: asunlit_veg_mask(mp) 
logical :: asunlit_mask(mp) 
!co-efficients usoughout init_radiation ` called from _albedo as well
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: xk(mp,nrb)

!iFor testing
ICYCLE = 0
cable_user%soil_struc="default"

  
! assign local ptrs to constants defined in cable_data_module
CALL point2constants(C)

IF ( cbl_standalone .OR. cable_runtime%um_explicit ) & 
  CALL ruff_resist(veg, rough, ssnow, canopy, veg%vlai, veg%hc, canopy%vlaiw)

!jhan: this call to define air may be redundant
CALL define_air (met, air)

!veg_mask =  canopy%vlaiw > .001
call fveg_mask( mp, Clai_thresh, canopy%vlaiw )
call fsunlit_mask( mp, Ccoszen_tols, met%coszen )
call fsunlit_veg_mask( mp )

CALL init_radiation( rad%extkb, rad%extkd,                                     &
                     !ExtCoeff_beam, ExtCoeff_dif,
                     rad%extkbm, rad%extkdm, Rad%Fbeam,                        &
                     !EffExtCoeff_beam, EffExtCoeff_dif, RadFbeam,
                     c1, rhoch, xk,                                            &
                     mp,nrb,                                                   &
                     Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180,         &
                     cbl_standalone, jls_standalone, jls_radiation,            &
                     subr_name,                                                &
                     veg_mask, sunlit_mask, sunlit_veg_mask,                   &
                     veg%Xfang, veg%taul, veg%refl,                            &
                     !VegXfang, VegTaul, VegRefl
                     met%coszen, int(met%DoY), met%fsd,                        &
                     !coszen, metDoY, SW_down,
                     canopy%vlaiw                                              &
                   ) !reducedLAIdue2snow 
 
IF ( cbl_standalone .OR. cable_runtime%um_explicit ) & 
  call Albedo( ssnow%AlbSoilsn, soil%AlbSoil,                                &
               !AlbSnow, AlbSoil,              
               mp, nrb,                                                      &
               jls_radiation ,                                               &
               veg_mask, sunlit_mask, sunlit_veg_mask,                       &  
               Ccoszen_tols, CGAUSS_W,                                       & 
               veg%iveg, veg%refl, veg%taul,                                 & 
               !surface_type, VegRefl, VegTaul,
               met%tk, met%coszen, canopy%vlaiw,                             &
               !metTk, coszen, reducedLAIdue2snow,
               ssnow%snowd, ssnow%osnowd, ssnow%isflag,                      & 
               !SnowDepth, SnowODepth, SnowFlag_3L, 
               ssnow%ssdnn, ssnow%tgg(:,1), ssnow%snage,                     & 
               !SnowDensity, SoilTemp, SnowAge, 
               xk, c1, rhoch,                                                & 
               rad%fbeam, rad%albedo,                                        &
               !RadFbeam, RadAlbedo,
               rad%extkd, rad%extkb,                                         & 
               !ExtCoeff_dif, ExtCoeff_beam,
               rad%extkdm, rad%extkbm,                                       & 
               !EffExtCoeff_dif, EffExtCoeff_beam,                
               rad%rhocdf, rad%rhocbm,                                       &
               !CanopyRefl_dif,CanopyRefl_beam,
               rad%cexpkdm, rad%cexpkbm,                                     & 
               !CanopyTransmit_dif, CanopyTransmit_beam, 
               rad%reffdf, rad%reffbm                                        &
             ) !EffSurfRefl_dif, EffSurfRefl_beam 


    ! Calculate canopy variables:

    !! vh_js !!
    !CABLE_LSM:check
    IF( cable_runtime%um .AND. first_call ) THEN
       ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)
       ssnow%otss = ssnow%tss
       first_call = .FALSE.
    ENDIF
    ssnow%otss_0 = ssnow%otss  ! vh should be before call to canopy?
    ssnow%otss = ssnow%tss

    CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate,sunlit_veg_mask, canopy%vlaiw)
    ! RML moved out of following IF after discussion with Eva
    ssnow%owetfac = ssnow%wetfac

    IF( cable_runtime%um ) THEN

       IF( cable_runtime%um_implicit ) THEN
          IF (cable_user%gw_model) THEN
             CALL soil_snow_gw(dels, soil, ssnow, canopy, met, bal,veg)
          ELSE
             CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
          ENDIF
       ENDIF

    ELSE
       IF(cable_user%SOIL_STRUC=='default') THEN
          IF (cable_user%gw_model) THEN
             CALL soil_snow_gw(dels, soil, ssnow, canopy, met, bal,veg)
          ELSE
             CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
          ENDIF
       ELSEIF (cable_user%SOIL_STRUC=='sli') THEN

          IF (cable_user%test_new_gw) &
               CALL sli_hydrology(dels,ssnow,soil,veg,canopy)

          CALL sli_main(ktau,dels,veg,soil,ssnow,met,canopy,air,rad,0)
       ENDIF
    ENDIF


    ssnow%deltss = ssnow%tss-ssnow%otss
    ! correction required for energy balance in online simulations
    ! REV_CORR - multiple changes to address %cls bugs and revised correction
    ! terms.  Also - do not apply correction terms if using SLI
    ! SSEB package will move these calculations to within soilsnow
    IF( cable_runtime%um .AND. cable_user%SOIL_STRUC=='default') THEN

       canopy%fhs = canopy%fhs + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
       canopy%fhs_cor = canopy%fhs_cor + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
       canopy%fh = canopy%fhv + canopy%fhs

       !canopy%fes = canopy%fes + ( ssnow%tss-ssnow%otss ) *                    &
       !          ( ssnow%dfe_ddq * ssnow%ddq_dtg )
       !          !( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )
       !
       !Ticket 137 - remove double couting of %cls
       !canopy%fes_cor = canopy%fes_cor + ( ssnow%tss-ssnow%otss ) *            &
       !                      ( ssnow%dfe_ddq * ssnow%ddq_dtg )
       !               ( ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg )

       !INH rewritten in terms of %dfe_dtg - NB factor %cls above was a bug
       canopy%fes = canopy%fes + ( ssnow%tss-ssnow%otss ) * ssnow%dfe_dtg

       !INH NB factor %cls in %fes_cor above was a bug - see Ticket #135 #137
       canopy%fes_cor = canopy%fes_cor + (ssnow%tss-ssnow%otss) * ssnow%dfe_dtg
       !canopy%fes_cor = canopy%fes_cor + ssnow%cls*(ssnow%tss-ssnow%otss) &
       !       * ssnow%dfe_dtg

       IF (cable_user%L_REV_CORR) THEN
          !INH need to add on corrections to all terms in the soil energy balance
          canopy%fns_cor = canopy%fns_cor + (ssnow%tss-ssnow%otss)*ssnow%dfn_dtg

          !NB %fns_cor also added onto out%Rnet and out%LWnet in cable_output and
          !cable_checks as the correction term needs to pass through the
          !canopy in entirity not be partially absorbed and %fns not used there
          !(as would be the case if rad%flws were changed)
          canopy%fns = canopy%fns + ( ssnow%tss-ssnow%otss )*ssnow%dfn_dtg

          canopy%ga_cor = canopy%ga_cor + ( ssnow%tss-ssnow%otss )*canopy%dgdtg
          canopy%ga = canopy%ga + ( ssnow%tss-ssnow%otss )*canopy%dgdtg

          !assign all the correction to %fes to %fess - none to %fesp
          canopy%fess = canopy%fess + ( ssnow%tss-ssnow%otss ) * ssnow%dfe_dtg

       ENDIF
    ENDIF


    ! need to adjust fe after soilsnow
    canopy%fev  = canopy%fevc + canopy%fevw

    ! Calculate total latent heat flux:
    canopy%fe = canopy%fev + canopy%fes

    ! Calculate net radiation absorbed by soil + veg
    canopy%rnet = canopy%fns + canopy%fnv

    ! Calculate radiative/skin temperature:
    IF (cable_runtime%um) THEN
       !Jan 2018: UM assumes a single emissivity for the surface in the radiation scheme
       !To accommodate this a single value of is 1. is assumed in ACCESS
       ! any leaf/soil emissivity /=1 must be incorporated into rad%trad.
       ! check that emissivities (pft and nvg) set = 1 within the UM i/o configuration
       ! CM2 - further adapted to pass the correction term onto %trad correctly
       rad%trad = ( ( 1.-rad%transd ) * C%emleaf * canopy%tv**4 +                      &
            rad%transd * C%emsoil * ssnow%otss**4 + canopy%fns_cor/C%sboltz )**0.25
    ELSE
       rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 +                             &
            rad%transd * ssnow%tss**4 )**0.25
    ENDIF

    ! rml 17/1/11 move all plant resp and soil resp calculations here
    ! from canopy. in UM only call on implicit step.
    ! put old and new soil resp calculations into soilcarb subroutine
    ! make new plantcarb subroutine
    IF (.NOT.cable_runtime%um_explicit .AND. icycle == 0) THEN

       !calculate canopy%frp
       CALL plantcarb(veg,bgc,met,canopy)

       !calculate canopy%frs
       CALL soilcarb(soil, ssnow, veg, bgc, met, canopy)

       CALL carbon_pl(dels, soil, ssnow, veg, canopy, bgc)

       canopy%fnpp = -1.0* canopy%fpn - canopy%frp
       canopy%fnee = canopy%fpn + canopy%frs + canopy%frp

    ENDIF


  END SUBROUTINE cbm

END MODULE cable_cbm_module
