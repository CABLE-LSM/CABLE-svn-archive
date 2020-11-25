MODULE cbl_model_driver_mod

IMPLICIT NONE

PRIVATE
PUBLIC cbl_model_driver

CONTAINS

SUBROUTINE cbl_model_driver( explicit_path, mp,nrb, land_pts, npft, ktau,dels, air,          &
                    bgc, canopy, met,                                         &
                    bal, rad, rough, soil,                                    &
                    ssnow, sum_flux, veg, z0surf_min,                         &
                    LAI_pft, HGT_pft, RmetDoY, reducedLAIdue2snow )

!subrs
USE cbl_albedo_mod, ONLY: albedo
USE cbl_hruff_mod,          ONLY: HgtAboveSnow
USE cbl_LAI_eff_mod,        ONLY: LAI_eff
USE cbl_masks_mod, ONLY: fveg_mask,  fsunlit_mask,  fsunlit_veg_mask
USE cbl_masks_mod, ONLY: veg_mask,  sunlit_mask,  sunlit_veg_mask

USE cbl_init_Loobos_mod, ONLY : initialize_Lveg
USE cbl_init_Loobos_mod, ONLY : initialize_Lsoil

!data
USE cable_other_constants_mod,  ONLY: Ccoszen_tols => coszen_tols

!jhan:pass these
USE cable_other_constants_mod, ONLY: clai_thresh => lai_thresh

USE cable_other_constants_mod, ONLY: cgauss_w => gauss_w
USE cable_math_constants_mod, ONLY: cpi => pi
USE cable_math_constants_mod, ONLY: cpi180 => pi180

USE cable_common_module
USE cable_carbon_module
USE cable_air_type_mod,       ONLY: air_type
USE cable_met_type_mod,       ONLY: met_type
USE cable_radiation_type_mod, ONLY: radiation_type
USE cable_roughness_type_mod, ONLY: roughness_type
USE cable_canopy_type_mod,    ONLY: canopy_type
USE cable_soil_snow_type_mod, ONLY: soil_snow_type
USE cable_bgc_pool_type_mod,  ONLY: bgc_pool_type
USE cable_balances_type_mod,  ONLY: balances_type
USE cable_sum_flux_type_mod,  ONLY: sum_flux_type
USE cable_params_mod,         ONLY: veg_parameter_type
USE cable_params_mod,         ONLY: soil_parameter_type
USE cable_def_types_mod,      ONLY: climate_type

USE cbl_soil_snow_main_module,  ONLY: soil_snow
USE cable_roughness_module, ONLY: ruff_resist
USE cbl_init_radiation_module, ONLY: init_radiation

USE cable_air_module, ONLY: define_air
   
USE cable_phys_constants_mod, ONLY : csboltz => sboltz
USE cable_phys_constants_mod, ONLY : cemleaf => emleaf
USE cable_phys_constants_mod, ONLY : cemsoil => emsoil
   
USE cable_canopy_module, ONLY: define_canopy
USE cable_canopy_module_explicit, ONLY: define_canopy_explicit
                                   
  !diag 
  USE cable_fprint_module, ONLY : cable_fprintf
  USE cable_Pyfprint_module, ONLY : cable_Pyfprintf
  USE cable_fFile_module, ONLY : fprintf_dir_root, fprintf_dir, L_cable_fprint,&
                                 L_cable_Pyfprint, unique_subdir, qprint
  USE cable_fprint_type_mod, ONLY : fprintx 

LOGICAL :: explicit_path
INTEGER :: mp
INTEGER :: nrb
INTEGER :: land_pts, npft
REAL :: AlbSoil(mp,nrb)
REAL :: MetTk(mp)
REAL :: SnowDepth(mp) 
REAL :: SnowDensity(mp)
REAL :: SnowODepth(mp)
REAL :: SoilTemp(mp)
REAL :: SnowAge(mp)
INTEGER:: SnowFlag_3L(mp)
INTEGER:: surface_type(mp)

REAL :: z0surf_min
REAL  :: hgt_pft(mp)
REAL  :: lai_pft(mp) 
REAL :: RmetDoY(mp)          !Day of the Year [formerly met%doy]
INTEGER :: metDoY(mp)          !Day of the Year [formerly met%doy]
  
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
LOGICAL, SAVE :: first_call = .TRUE.

CHARACTER(LEN=*), PARAMETER :: subr_name = "cbl_model_driver"
LOGICAL :: jls_standalone= .TRUE.
LOGICAL :: jls_radiation= .FALSE.
LOGICAL :: cbl_standalone = .FALSE.    

!make local to rad_driver and also again in cbl_model_driver
!CABLE variables to keep for all CABLE pathways across the timestep 
REAL :: reducedLAIdue2snow(mp)

!co-efficients usoughout init_radiation ` called from _albedo as well
REAL :: c1(mp,nrb)
REAL :: rhoch(mp,nrb)
REAL :: xk(mp,nrb)
REAL :: CanopyRefl_dif(mp,nrb)
REAL :: CanopyRefl_beam(mp,nrb)

integer :: cntile
real :: diagar(1)

!commented use of fprint to include in any file that uses fprint-ing
!==============================================================================
!use cable_Pyfprint_module
!==============================================================================
  character(len=30) :: vname
  character(len=70) :: dir
  
  integer :: dimx, dimy
   
  INTEGER, SAVE ::  cDiag00=0, cDiag1=0, cDiag2=0, cDiag3=0, cDiag4=0,         &
   cDiag5=0, cDiag6=0, cDiag7=0, cDiag8=0, cDiag9=0, cDiag10=0,    cDiag11=0,  &
   cDiag12=0, cDiag13=0, cDiag14=0, cDiag15=0, cDiag16=0, cDiag17=0, cDiag18=0,&
   cDiag19=0,cDiag20=0, cDiag21=0, cDiag22=0, cDiag23=0, cDiag24=0, cDiag25=0, &
   cDiag26=0, cDiag27=0, cDiag28=0, cDiag29=0, cDiag30=0, cDiag31=0, cDiag32=0,&
   cDiag33=0, cDiag34=0, cDiag35=0, cDiag36=0, cDiag37=0, cDiag38=0, cDiag39=0,&
   cDiag40=0, cDiag41=0, cDiag42=0, cDiag43=0, cDiag44=0, cDiag45=0, cDiag46=0,&
   cDiag47=0, cDiag48=0, cDiag49=0, cDiag50=0, cDiag51=0, cDiag52=0, cDiag53=0,&
   cDiag54=0, cDiag55=0, cDiag56=0, cDiag57=0, cDiag58=0, cDiag59=0, cDiag60=0,&
   cDiag61=0, cDiag62=0, cDiag63=0, cDiag64=0, cDiag65=0, cDiag66=0, cDiag67=0,&
   cDiag68=0, cDiag69=0

  logical :: L_fprint_HW = .true.
  logical :: L_fprint= .false. ! default
 
cntile =1

IF ( first_call ) THEN 
  call initialize_Lveg( veg )
  call initialize_Lsoil( soil )
  met%ca = 3.999999e-4 
  ALLOCATE( fprintx% fprint1(mp) )
  ALLOCATE( fprintx% fprint2(mp) )
  ALLOCATE( fprintx% fprint3(mp) )
  ALLOCATE( fprintx% fprint4(mp) )
  ALLOCATE( fprintx% fprint5(mp) )
  ALLOCATE( fprintx% fprint6(mp) )
  ALLOCATE( fprintx% fprint7(mp) )
  ALLOCATE( fprintx% fprint8(mp) )
  ALLOCATE( fprintx% fprint9(mp) )
EndIf
      
metDoy = INT(RmetDoy)
!iFor testing
!ICYCLE = 0
cable_user%soil_struc="default"

CALL ruff_resist( veg, rough, ssnow, canopy, veg%vlai, veg%hc, canopy%vlaiw )
!H!CALL ruff_resist( veg, rough, ssnow, canopy, LAI_pft, HGT_pft, reducedLAIdue2snow )
reducedLAIdue2snow = canopy%vlaiw
!jhan: this call to define air may be redundant
CALL define_air (met, air)

call fveg_mask( veg_mask, mp, Clai_thresh, canopy%vlaiw )
call fsunlit_mask( sunlit_mask, mp, Ccoszen_tols, met%coszen )
!call fsunlit_mask( sunlit_mask, mp, Ccoszen_tols,( met%fsd(:,1)+met%fsd(:,2) ) )
call fsunlit_veg_mask( sunlit_veg_mask, mp )

CALL init_radiation( rad%extkb, rad%extkd,                                    &
                     !ExtCoeff_beam, ExtCoeff_dif,
                     rad%extkbm, rad%extkdm, Rad%Fbeam,                       &
                     !EffExtCoeff_beam, EffExtCoeff_dif, RadFbeam,
                     c1, rhoch, xk,                                           &
                     mp,nrb,                                                  &
                     Clai_thresh, Ccoszen_tols, CGauss_w, Cpi, Cpi180,        &
                     cbl_standalone, jls_standalone, jls_radiation,           &
                     subr_name,                                               &
                     veg_mask, sunlit_mask, sunlit_veg_mask,                  &
                     veg%Xfang, veg%taul, veg%refl,                           &
                     !VegXfang, VegTaul, VegRefl
                     met%coszen, INT(met%DoY), met%fsd,                       &
                     !coszen, metDoY, SW_down,
                     canopy%vlaiw                                              &
                   ) !reducedLAIdue2snow 
 
!!IF ( cable_runtime%um_explicit )                                              &
  CALL Albedo( ssnow%AlbSoilsn, soil%AlbSoil,                                 &
               !AlbSnow, AlbSoil,              
               mp, nrb,                                                       &
               jls_radiation,                                                 &
               veg_mask, sunlit_mask, sunlit_veg_mask,                        &  
               Ccoszen_tols, cgauss_w,                                        & 
               veg%iveg, veg%refl, veg%taul,                                  & 
               !surface_type, VegRefl, VegTaul,
               met%tk, met%coszen, canopy%vlaiw,                              &
               !metTk, coszen, reducedLAIdue2snow,
               ssnow%snowd, ssnow%osnowd, ssnow%isflag,                       & 
               !SnowDepth, SnowODepth, SnowFlag_3L, 
               ssnow%ssdnn, ssnow%tgg(:,1), ssnow%tggsn(:,1), ssnow%snage,                      & 
               !SnowDensity, SoilTemp, SnowAge, 
               xk, c1, rhoch,                                                 & 
               rad%fbeam, rad%albedo,                                         &
               !RadFbeam, RadAlbedo,
               rad%extkd, rad%extkb,                                          & 
               !ExtCoeff_dif, ExtCoeff_beam,
               rad%extkdm, rad%extkbm,                                        & 
               !EffExtCoeff_dif, EffExtCoeff_beam,                
               rad%rhocdf, rad%rhocbm,                                        &
               !CanopyRefl_dif,CanopyRefl_beam,
               rad%cexpkdm, rad%cexpkbm,                                      & 
               !CanopyTransmit_dif, CanopyTransmit_beam, 
               rad%reffdf, rad%reffbm                                        &
             ) !EffSurfRefl_dif, EffSurfRefl_beam 

!CABLE_LSM:check
!IF ( explicit_path ) THEN
IF ( ktau_gl == 1) THEN
  ssnow%tss=(1 - ssnow%isflag) * ssnow%tgg(:,1) + ssnow%isflag * ssnow%tggsn(:,1) 
  ssnow%otss = ssnow%tss
  first_call = .FALSE.
ENDIF
ssnow%otss_0 = ssnow%otss  ! vh should be before call to canopy?
ssnow%otss = ssnow%tss

IF ( explicit_path ) THEN
  CALL define_canopy_explicit(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate, sunlit_veg_mask,  canopy%vlaiw)
ELSE
  CALL define_canopy         ( bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate, &
  sunlit_veg_mask, reducedLAIdue2snow )
END IF

ssnow%owetfac = ssnow%wetfac

IF ( cable_runtime%um_implicit ) &
  CALL soil_snow(dels, soil, ssnow, canopy, met, bal,veg)

ssnow%deltss = ssnow%tss - ssnow%otss

! need to adjust fe after soilsnow
canopy%fev  = canopy%fevc + canopy%fevw

! Calculate total latent heat flux:
canopy%fe = canopy%fev + canopy%fes

! Calculate net radiation absorbed by soil + veg
canopy%rnet = canopy%fns + canopy%fnv

! Calculate radiative/skin temperature:
rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 + rad%transd * ssnow%tss**4 )**0.25
! Calculate radiative/skin temperature:
       rad%trad = ( ( 1.-rad%transd ) * canopy%tv**4 +                             &
            rad%transd * ssnow%tss**4 )**0.25

    !H!IF (.NOT.cable_runtime%um_explicit .AND. icycle == 0) THEN

    !H!   !calculate canopy%frp
    !H!   CALL plantcarb(veg,bgc,met,canopy)

    !H!   !calculate canopy%frs
    !H!   CALL soilcarb(soil, ssnow, veg, bgc, met, canopy)

    !H!   CALL carbon_pl(dels, soil, ssnow, veg, canopy, bgc)

    !H!   canopy%fnpp = -1.0* canopy%fpn - canopy%frp
    !H!   canopy%fnee = canopy%fpn + canopy%frs + canopy%frp

    !H!ENDIF


fprintf_dir="/home/599/jxs599/TestHAC5.7_HAC_impl/"
IF ( cable_runtime%um_implicit ) then 
!H!#forcing
!H!vname='met%tk'; dimx=1 
!H!diagar(1) = met%tk(cntile)
!H!call cable_Pyfprintf( cDiag16, vname, diagar, dimx, .true.)
!H!
!H!vname='met%fsd1'; dimx=1 
!H!diagar(1) = met%fsd(cntile,1)
!H!call cable_Pyfprintf( cDiag17, vname, diagar, dimx, .true.)
!H!
!H!vname='met%fsd2'; dimx=1 
!H!diagar(1) = met%fsd(cntile,2)
!H!call cable_Pyfprintf( cDiag18, vname, diagar, dimx, .true.)
!H!
!H!vname='met%fld'; dimx=1 
!H!diagar(1) = met%fld(cntile)
!H!call cable_Pyfprintf( cDiag19, vname, diagar, dimx, .true.)
!H!
!H!vname='met%precip'; dimx=1 
!H!diagar(1) = met%precip(cntile)
!H!call cable_Pyfprintf( cDiag20, vname, diagar, dimx, .true.)
!H!
!H!vname='met%precip_sn'; dimx=1 
!H!diagar(1) = met%precip_sn(cntile)
!H!call cable_Pyfprintf( cDiag21, vname, diagar, dimx, .true.)
!H!
!H!vname='met%ua'; dimx=1 
!H!diagar(1) = met%ua(cntile)
!H!call cable_Pyfprintf( cDiag22, vname, diagar, dimx, .true.)
!H!
!H!vname='met%pmb'; dimx=1 
!H!diagar(1) = met%pmb(cntile)
!H!call cable_Pyfprintf( cDiag23, vname, diagar, dimx, .true.)
!H!
!H!vname='met%qv'; dimx=1 
!H!diagar(1) = met%qv(cntile)
!H!call cable_Pyfprintf( cDiag24, vname, diagar, dimx, .true.)

vname='canopy%fhv'; dimx=1 
diagar(1) = canopy%fhv(cntile)
call cable_Pyfprintf( cDiag1, vname, diagar, dimx, .true.)

vname='canopy%fhs'; dimx=1 
diagar(1) = canopy%fhs(cntile)
call cable_Pyfprintf( cDiag2, vname, diagar, dimx, .true.)

vname='canopy%fes'; dimx=1 
diagar(1) = canopy%fes(cntile)
call cable_Pyfprintf( cDiag3, vname, diagar, dimx, .true.)

vname='canopy%fev'; dimx=1 
diagar(1) = canopy%fev(cntile)
call cable_Pyfprintf( cDiag4, vname, diagar, dimx, .true.)

vname='canopy%fess'; dimx=1 
diagar(1) = canopy%fess(cntile)
call cable_Pyfprintf( cDiag5, vname, diagar, dimx, .true.)

vname='canopy%fesp'; dimx=1 
diagar(1) = canopy%fesp(cntile)
call cable_Pyfprintf( cDiag6, vname, diagar, dimx, .true.)

vname='ssnow%wetfac'; dimx=1 
diagar(1) = ssnow%wetfac(cntile)
call cable_Pyfprintf( cDiag6, vname, diagar, dimx, .true.)

vname='ssnow%potev'; dimx=1 
diagar(1) = ssnow%potev(cntile)
call cable_Pyfprintf( cDiag8, vname, diagar, dimx, .true.)

vname='ssnow%cls'; dimx=1 
diagar(1) = ssnow%cls(cntile)
call cable_Pyfprintf( cDiag9, vname, diagar, dimx, .true.)

vname='ssnow%qstss'; dimx=1 
diagar(1) = ssnow%qstss(cntile)
call cable_Pyfprintf( cDiag10, vname, diagar, dimx, .true.)

!for soilsnow_main?
!
!vname='mwbice'; dimx=1 
!diagar(1) = fprintx%fprint3(cntile)
!call cable_Pyfprintf( cDiag29, vname, diagar, dimx, .true.)
!
!vname='mwblf'; dimx=1 
!diagar(1) = fprintx%fprint4(cntile)
!call cable_Pyfprintf( cDiag31, vname, diagar, dimx, .true.)
!
!vname='mwbfice'; dimx=1 
!diagar(1) = fprintx%fprint5(cntile)
!call cable_Pyfprintf( cDiag32, vname, diagar, dimx, .true.)
!
!vname='mwb'; dimx=1 
!diagar(1) = fprintx%fprint6(cntile)
!call cable_Pyfprintf( cDiag33, vname, diagar, dimx, .true.)
!??
!vname='mgammzz'; dimx=1 
!diagar(1) = fprintx%fprint2(cntile)
!call cable_Pyfprintf( cDiag23, vname, diagar, dimx, .true.)

endif


!!IF ( cable_runtime%um_explicit ) then 
!!fprintf_dir="/home/599/jxs599/TestHAC5.7_HAC_expl/"
!!vname='canopy%fhv'; dimx=1 
!!diagar(1) = canopy%fhv(cntile)
!!call cable_Pyfprintf( cDiag1, vname, diagar, dimx, .true.)
!!vname='canopy%fhs'; dimx=1 
!!diagar(1) = canopy%fhs(cntile)
!!call cable_Pyfprintf( cDiag2, vname, diagar, dimx, .true.)
!!vname='canopy%fes'; dimx=1 
!!diagar(1) = canopy%fes(cntile)
!!call cable_Pyfprintf( cDiag3, vname, diagar, dimx, .true.)
!!vname='canopy%fev'; dimx=1 
!!diagar(1) = canopy%fev(cntile)
!!call cable_Pyfprintf( cDiag4, vname, diagar, dimx, .true.)
!!vname='met%tk'; dimx=1 
!!diagar(1) = met%tk(cntile)
!!call cable_Pyfprintf( cDiag5, vname, diagar, dimx, .true.)
!!endif



END SUBROUTINE cbl_model_driver

END MODULE cbl_model_driver_mod


