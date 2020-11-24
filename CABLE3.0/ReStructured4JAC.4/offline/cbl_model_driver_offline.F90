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
! physical constants
USE cable_phys_constants_mod, ONLY : CGRAV  => GRAV
USE cable_phys_constants_mod, ONLY : CCAPP   => CAPP
USE cable_phys_constants_mod, ONLY : CEMLEAF => EMLEAF
USE cable_phys_constants_mod, ONLY : CEMSOIL => EMSOIL
USE cable_phys_constants_mod, ONLY : CSBOLTZ => SBOLTZ
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
use cbl_masks_mod, ONLY :  veg_mask,  sunlit_mask,  sunlit_veg_mask

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
character(len=*), parameter :: subr_name = "cbm"
LOGICAL :: cbl_standalone= .true.
LOGICAL :: jls_standalone= .false.
LOGICAL :: jls_radiation= .false.

!make local to rad_driver and also again in cbl_model_driver
!CABLE variables to keep for all CABLE pathways across the timestep 
real :: reducedLAIdue2snow(mp)

!masks
!logical :: veg_mask(mp),  sunlit_mask(mp),  sunlit_veg_mask(mp) 
!logical :: sunlit_veg_mask(mp) 
!logical :: sunlit_mask(mp) 
!co-efficients usoughout init_radiation ` called from _albedo as well
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: xk(mp,nrb)

!iFor testing
ICYCLE = 0
cable_user%soil_struc="default"
soil%cnsd_vec = 0.
CALL ruff_resist( veg, rough, ssnow, canopy, veg%vlai, veg%hc, canopy%vlaiw )

!jhan: this call to define air may be redundant
CALL define_air (met, air)

call fveg_mask( veg_mask, mp, Clai_thresh, canopy%vlaiw )
call fsunlit_mask( sunlit_mask, mp, Ccoszen_tols, met%coszen )
!call fsunlit_mask( sunlit_mask, mp, Ccoszen_tols,( met%fsd(:,1)+met%fsd(:,2) ) )
call fsunlit_veg_mask( sunlit_veg_mask, mp )

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
 
  call Albedo( ssnow%AlbSoilsn, soil%AlbSoil,                                &
               !AlbSnow, AlbSoil,              
               mp, nrb,                                                      &
               jls_radiation,                                                 &
               veg_mask, sunlit_mask, sunlit_veg_mask,                       &  
               Ccoszen_tols, CGAUSS_W,                                       & 
               veg%iveg, veg%refl, veg%taul,                                 & 
               !surface_type, VegRefl, VegTaul,
               met%tk, met%coszen, canopy%vlaiw,                             &
               !metTk, coszen, reducedLAIdue2snow,
               ssnow%snowd, ssnow%osnowd, ssnow%isflag,                      & 
               !SnowDepth, SnowODepth, SnowFlag_3L, 
               ssnow%ssdnn, ssnow%tgg(:,1), ssnow%tggsn(:,1), ssnow%snage,                     & 
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

!Otherwise this is only initialized for CABLE-offline IO
IF ( ktau_gl == 1) THEN
  ssnow%tss=(1 - ssnow%isflag) * ssnow%tgg(:,1) + ssnow%isflag * ssnow%tggsn(:,1)
       ssnow%otss = ssnow%tss
ENDIF
ssnow%otss_0 = ssnow%otss  ! vh should be before call to canopy?
ssnow%otss = ssnow%tss

  CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate, sunlit_veg_mask,  canopy%vlaiw)
    
ssnow%owetfac = ssnow%wetfac

     CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)
      
ssnow%deltss = ssnow%tss - ssnow%otss

! need to adjust fe after soilsnow
canopy%fev  = canopy%fevc + canopy%fevw

! Calculate total latent heat flux:
canopy%fe = canopy%fev + canopy%fes

! Calculate net radiation absorbed by soil + veg
canopy%rnet = canopy%fns + canopy%fnv

    ! Calculate radiative/skin temperature:
       rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 +                             &
            rad%transd * ssnow%tss**4 )**0.25

    !calculate canopy%frp
    CALL plantcarb(veg,bgc,met,canopy)

    !calculate canopy%frs
    CALL soilcarb(soil, ssnow, veg, bgc, met, canopy)

    CALL carbon_pl(dels, soil, ssnow, veg, canopy, bgc)

    canopy%fnpp = -1.0* canopy%fpn - canopy%frp
    canopy%fnee = canopy%fpn + canopy%frs + canopy%frp





  END SUBROUTINE cbm

END MODULE cable_cbm_module
