MODULE cbl_model_driver_mod

  IMPLICIT NONE

  PRIVATE
  PUBLIC cbl_model_driver

CONTAINS

SUBROUTINE cbl_model_driver( mp,nrb, land_pts, npft, ktau,dels, air,   &
                    bgc, canopy, met,                                &
                    bal, rad, rough, soil,                           &
                    ssnow, sum_flux, veg, z0surf_min,       &
                    LAI_pft, HGT_pft, RmetDoY  )

!subrs
use cable_wide_mod, ONLY : allocate_cable_wide     
USE cbl_albedo_mod, ONLY : albedo
USE cbl_hruff_mod,          ONLY : HgtAboveSnow
USE cbl_LAI_eff_mod,        ONLY : LAI_eff
USE radiation_albedo_mod, ONLY : allocate_rad_albedo
use cbl_masks_mod, ONLY :  fveg_mask,  fsunlit_mask,  fsunlit_veg_mask
!data
USE cable_other_constants_mod,  ONLY : Ccoszen_tols => coszen_tols
use cable_wide_mod,             ONLY : reducedLAIdue2snow
USE radiation_albedo_mod,       ONLY : AlbSnow, coszen
USE radiation_albedo_mod,       ONLY : ExtCoeff_beam, ExtCoeff_dif, &
                                       EffExtCoeff_beam, EffExtCoeff_dif, &
                                       CanopyRefl_dif, CanopyRefl_beam, &
                                       CanopyTransmit_dif, CanopyTransmit_beam, &
                                       EffSurfRefl_dif, EffSurfRefl_beam, &
                                       coszen,        &
                                       VegXfang, VegTaul, VegRefl, c1, rhoch, &
                                       SW_down, RadFbeam, xk, RadAlbedo 

!jhan:pass these
use cbl_masks_mod, ONLY :  veg_mask,  sunlit_mask,  sunlit_veg_mask
USE cable_other_constants_mod, ONLY : CLAI_THRESH => lai_thresh

USE cable_other_constants_mod, ONLY : CGAUSS_W => gauss_w
USE cable_math_constants_mod, ONLY : CPI => pi
USE cable_math_constants_mod, ONLY : CPI180 => pi180

   USE cable_common_module
   USE cable_carbon_module
   USE cable_soil_snow_module,  ONLY : soil_snow
USE cable_def_types_mod,        ONLY : air_type, bgc_pool_type, &
                                    canopy_type, met_type, balances_type, &
                                    radiation_type, roughness_type,  &
                                    soil_snow_type, sum_flux_type, &
                                    climate_type, soil_parameter_type, &
                                    veg_parameter_type

   USE cable_roughness_module, only : ruff_resist
   USE cbl_init_radiation_module, only : init_radiation

   USE cable_air_module, only : define_air
   
   USE cable_data_module, ONLY : icbm_type, point2constants
   
   USE cable_canopy_module, only : define_canopy
                                   
integer :: mp
integer :: nrb
integer :: land_pts, npft
REAL :: AlbSoil(mp,nrb)
REAL :: MetTk(mp)
real :: SnowDepth(mp) 
real :: SnowDensity(mp)
REAL :: SnowODepth(mp)
REAL :: SoilTemp(mp)
REAL :: SnowAge(mp)
integer:: SnowFlag_3L(mp)
integer:: surface_type(mp)

 
real :: z0surf_min
REAL  :: hgt_pft(mp)
REAL  :: lai_pft(mp) 
real :: RmetDoY(mp)          !Day of the Year [formerly met%doy]
integer :: metDoY(mp)          !Day of the Year [formerly met%doy]
  
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
   TYPE (climate_type) :: climate

   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type),  INTENT(INOUT)    :: veg

   REAL, INTENT(IN)               :: dels ! time setp size (s)
   INTEGER, INTENT(IN) :: ktau
   logical, save :: first_call = .true.

character(len=*), parameter :: subr_name = "cbl_model_driver"
LOGICAL :: jls_standalone= .true.
LOGICAL :: jls_radiation= .false.

  ! assign local ptrs to constants defined in cable_data_module
  CALL point2constants(C)

IF( cable_runtime%um_explicit ) CALL ruff_resist( veg, rough, ssnow, canopy, &
                                                    LAI_pft, HGT_pft )
      
! define variables common cross CABLE (or at least >2 pathways)
IF(first_call) call allocate_cable_wide( mp )
! define variables common to rad/albedo pathway
IF(first_call) call allocate_rad_albedo( mp, nrb )

CALL define_air (met, air)

!we set args in "cbl_init_radiation_args.inc" - just in case here
SW_down = met%fsd
!we can grab these from JaC later
VegXfang = Veg%Xfang
VegTaul = Veg%Taul
VegRefl = Veg%Refl
reducedLAIdue2snow = canopy%Vlaiw

!Define logical masks according to vegetation cover and sunlight. this is also 
!done in radiation pathway but that could be out of step with explicit call
call fveg_mask( mp, reducedLAIdue2snow)
call fsunlit_mask( mp, nrb, Ccoszen_tols, met%coszen )
call fsunlit_veg_mask( mp, veg_mask, sunlit_mask )

metDoy = int(RmetDoy)
CALL init_radiation( &
#include                  "cbl_init_radiation_args.inc"                  
 )
  
rad%extkd = ExtCoeff_dif
rad%extkb = ExtCoeff_beam

AlbSoil        =     soil%AlbSoil
surface_type   =     veg%iveg
metTk          =     met%Tk
coszen         =     met%coszen
SnowDepth      =     ssnow%snowd
SnowODepth     =     ssnow%osnowd
SnowFlag_3L    =     ssnow%isflag
SnowDensity    =     ssnow%ssdnn
SoilTemp       =     ssnow%tgg(:,1)
SnowAge        =     ssnow%snage

IF( cable_runtime%um_explicit ) then
  call Albedo(  &
#             include "cbl_albedo_args.inc" 
             )

rad%extkbm = EffExtCoeff_beam
rad%extkdm = EffExtCoeff_dif
rad%cexpkdm = CanopyTransmit_dif 
rad%cexpkbm = CanopyTransmit_beam

rad%reffdf = EffSurfRefl_dif
rad%reffbm = EffSurfRefl_beam

Rad%Albedo = RadAlbedo
Rad%Fbeam  = RadFbeam 
endif
!need to do this for the rest of CABLE

   !CABLE_LSM:check
  IF( first_call ) then
    ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1) 
    ssnow%otss = ssnow%tss
    first_call = .false.
  endif
  ssnow%otss_0 = ssnow%otss  ! vh should be before call to canopy?
  ssnow%otss = ssnow%tss

  CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate)

  ssnow%owetfac = ssnow%wetfac

  IF( cable_runtime%um_implicit ) CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)

  ssnow%deltss = ssnow%tss-ssnow%otss

  canopy%fhs = canopy%fhs + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
  canopy%fhs_cor = canopy%fhs_cor + ( ssnow%tss-ssnow%otss ) * ssnow%dfh_dtg
  canopy%fh = canopy%fhv + canopy%fhs

  !INH rewritten in terms of %dfe_dtg - NB factor %cls above was a bug
  canopy%fes = canopy%fes + ( ssnow%tss-ssnow%otss ) * ssnow%dfe_dtg

  !INH NB factor %cls in %fes_cor above was a bug - see Ticket #135 #137
  canopy%fes_cor = canopy%fes_cor + (ssnow%tss-ssnow%otss) * ssnow%dfe_dtg
   
  !INH need to add on corrections to all terms in the soil energy balance
  canopy%fns_cor = canopy%fns_cor + (ssnow%tss-ssnow%otss)*ssnow%dfn_dtg

  canopy%fns = canopy%fns + ( ssnow%tss-ssnow%otss )*ssnow%dfn_dtg

  canopy%ga_cor = canopy%ga_cor + ( ssnow%tss-ssnow%otss )*canopy%dgdtg
  canopy%ga = canopy%ga + ( ssnow%tss-ssnow%otss )*canopy%dgdtg

  !assign all the correction to %fes to %fess - none to %fesp
  canopy%fess = canopy%fess + ( ssnow%tss-ssnow%otss ) * ssnow%dfe_dtg

  ! need to adjust fe after soilsnow
  canopy%fev  = canopy%fevc + canopy%fevw

  ! Calculate total latent heat flux:
  canopy%fe = canopy%fev + canopy%fes

  ! Calculate net radiation absorbed by soil + veg
  canopy%rnet = canopy%fns + canopy%fnv

  ! CM2 - further adapted to pass the correction term onto %trad correctly
  rad%trad = ( ( 1.-rad%transd ) * C%emleaf * canopy%tv**4 +                      &
         rad%transd * C%emsoil * ssnow%otss**4 + canopy%fns_cor/C%sboltz )**0.25

  CALL plantcarb(veg,bgc,met,canopy)

  !calculate canopy%frs
  CALL soilcarb(soil, ssnow, veg, bgc, met, canopy)

  CALL carbon_pl(dels, soil, ssnow, veg, canopy, bgc)

  canopy%fnpp = -1.0* canopy%fpn - canopy%frp
  canopy%fnee = canopy%fpn + canopy%frs + canopy%frp

END SUBROUTINE cbl_model_driver

END MODULE cbl_model_driver_mod


