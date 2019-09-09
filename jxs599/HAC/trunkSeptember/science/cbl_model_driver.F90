MODULE cbl_model_driver_mod

  IMPLICIT NONE

  PRIVATE
  PUBLIC cbl_model_driver

CONTAINS

SUBROUTINE cbl_model_driver( mp,nrb, land_pts, npft, ktau,dels, air,   &
                    bgc, canopy, met,                                &
                    bal, rad, rough, soil,                           &
                    ssnow, sum_flux, veg, z0surf_min,       &
                    LAI_pft, HGT_pft, RmetDoY, reducedLAIdue2snow )

!subrs
USE cbl_albedo_mod, ONLY : albedo
USE cbl_hruff_mod,          ONLY : HgtAboveSnow
USE cbl_LAI_eff_mod,        ONLY : LAI_eff
USE radiation_albedo_mod, ONLY : allocate_rad_albedo
use cbl_masks_mod, ONLY :  fveg_mask,  fsunlit_mask,  fsunlit_veg_mask
!data
USE cable_other_constants_mod,  ONLY : Ccoszen_tols => coszen_tols

!jhan:pass these
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
                                   
USE legacy_support_mod, ONLY : legacy_support

integer :: mp
integer :: nrb
integer :: land_pts, npft
!REAL :: AlbSoil(mp,nrb)
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

!make local to rad_driver and also again in cbl_model_driver
!CABLE variables to keep for all CABLE pathways across the timestep 
real :: reducedLAIdue2snow(mp)
!masks
logical :: veg_mask(mp),  sunlit_mask(mp),  sunlit_veg_mask(mp) 

!co-efficients usoughout init_radiation ` called from _albedo as well
!REAL :: c1(mp,nrb)
!REAL :: rhoch(mp,nrb)
!REAL :: xk(mp,nrb)
!REAL :: CanopyRefl_dif(mp,nrb)
!REAL :: CanopyRefl_beam(mp,nrb)
! Albedos
!-------------------------------------------------------------------------------
REAL :: AlbSoil(mp,nrb)             !Bare Soil Albedo - parametrized (soil%albsoil)
REAL :: AlbSnow(mp,nrb)             !Ground Albedo given a snow coverage (ssnow%albsoilsn)
REAL :: RadAlbedo(mp,nrb)           !Total albedo given RadFbeam (rad%albedo)
REAL :: EffSurfRefl_dif(mp,nrb)     !Effective Surface Relectance as seen by atmosphere [Diffuse SW]  (rad%reffdf)
REAL :: EffSurfRefl_beam(mp,nrb)    !Effective Surface Relectance as seen by atmosphere [Direct Beam SW] (rad%reffbm)
!-------------------------------------------------------------------------------
REAL :: RadFbeam(mp,nrb)            !Computed Beam Fraction given total SW (rad%fbeam)

!common radiation scalings [computed in albedo() ]
!-------------------------------------------------------------------------------
REAL :: xk(mp,nrb)
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
!-------------------------------------------------------------------------------

!Variables shared primarily between radiation and albedo and possibly elsewhere
!-------------------------------------------------------------------------------
!Extinction co-efficients compued in init_radiation()
REAL :: ExtCoeff_beam(mp)           !"raw" Extinction co-efficient for Direct Beam component of SW radiation (rad%extkb)
REAL :: ExtCoeff_dif(mp)            !"raw"Extinction co-efficient for Diffuse component of SW radiation (rad%extkd)
REAL :: EffExtCoeff_beam(mp,nrb)    !Effective Extinction co-efficient for Direct Beam component of SW radiation (rad%extkbm)
REAL :: EffExtCoeff_dif(mp,nrb)     !Effective Extinction co-efficient for Diffuse component of SW radiation (rad%extkdm)

!Canopy reflectance/transmitance compued in albedo() 
REAL :: CanopyRefl_dif(mp,nrb)      !Canopy reflectance (rad%cexpkdm) 
REAL :: CanopyRefl_beam(mp,nrb)     !Canopy reflectance (rad%cexpkbm)   
REAL :: CanopyTransmit_dif(mp,nrb)  !Canopy Transmitance (rad%rhodf   
REAL :: CanopyTransmit_beam(mp,nrb) !Canopy Transmitance (rad%rhobm)    
!-------------------------------------------------------------------------------



  ! assign local ptrs to constants defined in cable_data_module
  CALL point2constants(C)

IF( cable_runtime%um_explicit ) CALL ruff_resist( veg, rough, ssnow, canopy, &
                                                    LAI_pft, HGT_pft,        & 
                                                    reducedLAIdue2snow )
      
CALL define_air (met, air)

!Define logical masks according to vegetation cover and sunlight. this is also 
!done in radiation pathway but that could be out of step with explicit call
!reducedLAIdue2snow = canopy%Vlaiw
canopy%Vlaiw = reducedLAIdue2snow
call fveg_mask( veg_mask,mp, Clai_thresh, reducedLAIdue2snow)
call fsunlit_mask( sunlit_mask, mp, Ccoszen_tols, met%coszen )
call fsunlit_veg_mask( sunlit_veg_mask, mp,  veg_mask, sunlit_mask )

metDoy = int(RmetDoy)
CALL init_radiation(   &
mp,                    &  
nrb,                   &
Clai_thresh,           &
Ccoszen_tols,          &
jls_standalone,        &
jls_radiation ,        &
veg_mask,              &
sunlit_mask,           &
sunlit_veg_mask,       &
reducedLAIdue2snow,    &
met%coszen,            &!coszen,                &
!get
rad%extkb,              &   !ExtCoeff_beam,         &
rad%extkd,              &   !ExtCoeff_dif,          &
rad%extkbm,            & ! EffExtCoeff_beam
rad%extkdm,            &! = EffExtCoeff_dif
!params
Veg%Xfang,             &
Veg%Taul,              &
Veg%Refl,              &
!
c1,                    &
rhoch,                 &
metDoY,                &
met%fsd,               &
Rad%Fbeam,             &
xk,                    &
!constants
CGauss_w,              &
Cpi,                   &
Cpi180,                &
subr_name              &
 )
  
!rad%extkd = ExtCoeff_dif
!rad%extkb = ExtCoeff_beam
!rad%extkbm = EffExtCoeff_beam
!rad%extkdm = EffExtCoeff_dif

!AlbSoil        =     soil%AlbSoil
!surface_type   =     veg%iveg
!metTk          =     met%Tk
!coszen         =     met%coszen
!SnowDepth      =     ssnow%snowd
!SnowODepth     =     ssnow%osnowd
!SnowFlag_3L    =     ssnow%isflag
!SnowDensity    =     ssnow%ssdnn
!SoilTemp       =     ssnow%tgg(:,1)
!SnowAge        =     ssnow%snage

IF( cable_runtime%um_explicit ) then
  call Albedo(        &
ssnow%AlbSoilsn,      &!AlbSnow,              & 
soil%AlbSoil,         &!AlbSoil,              & 
mp,                   &  
nrb,                  &
jls_radiation ,       &
veg_mask,             & 
sunlit_mask,          & 
sunlit_veg_mask,      & 
Ccoszen_tols,         &
CGAUSS_W,             & 
veg%iveg,             & !   surface_type,         &
met%tk,               & !  metTk,                & 
met%coszen,           &!  coszen,               & 
canopy%vlaiw,         &!  reducedLAIdue2snow,          &
ssnow%snowd,          &!  SnowDepth,            &  
ssnow%osnowd,         &!  SnowODepth,           & 
ssnow%isflag,         & !  SnowFlag_3L,          & 
ssnow%ssdnn,          & !   SnowDensity,          & 
ssnow%tgg(:,1),       &   !   SoilTemp,             & 
ssnow%snage,          & !   SnowAge,              &
xk,                   &  
c1,                   &  
rhoch,                &
Rad%Fbeam,            & 
Rad%Albedo,           &
rad%extkd,             &!ExtCoeff_beam,         &
rad%extkb,             & !ExtCoeff_dif,          &
rad%extkdm,           & ! EffExtCoeff_beam
rad%extkbm,           & ! = EffExtCoeff_dif
  CanopyRefl_dif,    & 
  CanopyRefl_beam,   &
rad%cexpkdm,          & ! = CanopyTransmit_dif 
rad%cexpkbm,          & ! = CanopyTransmit_beam
rad%reffdf,           &! = EffSurfRefl_dif
rad%reffbm            &! = EffSurfRefl_beam
             )
!!!rad%extkbm = EffExtCoeff_beam
!!!rad%extkdm = EffExtCoeff_dif
!rad%cexpkdm = CanopyTransmit_dif 
!rad%cexpkbm = CanopyTransmit_beam
!rad%reffdf = EffSurfRefl_dif
!rad%reffbm = EffSurfRefl_beam
!Rad%Albedo = RadAlbedo
!Rad%Fbeam  = RadFbeam 
endif
!need to do this for the rest of CABLE
call legacy_support( mp, nrb, ExtCoeff_beam, ExtCoeff_dif, &
EffExtCoeff_beam, EffExtCoeff_dif, &
CanopyRefl_dif, CanopyRefl_beam, &
CanopyTransmit_dif, CanopyTransmit_beam,&
RadFbeam, AlbSnow, &
RadAlbedo,AlbSoil, &
EffSurfRefl_dif, EffSurfRefl_beam, &
canopy, rad, ssnow )


   !CABLE_LSM:check
  IF( first_call ) then
    ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1) 
    ssnow%otss = ssnow%tss
    first_call = .false.
  endif
  ssnow%otss_0 = ssnow%otss  ! vh should be before call to canopy?
  ssnow%otss = ssnow%tss

  CALL define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate, sunlit_veg_mask, reducedLAIdue2snow )

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


