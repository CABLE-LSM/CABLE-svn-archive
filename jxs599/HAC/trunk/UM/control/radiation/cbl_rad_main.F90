module cable_rad_main_mod
  
contains

SUBROUTINE cable_rad_main( &
#                          include "cbl_rad_main_args.inc"
                         ) 
 
!subrs
USE cable_rad_driv_mod,   ONLY : cable_rad_driver
USE cable_rad_unpack_mod, ONLY : cable_rad_unpack
USE cable_common_module, ONLY : cable_runtime
!subrs:HaC1.3 re-distrbute
USE cable_um_init_subrs_mod, ONLY : um2cable_lp !packing
USE cable_um_init_subrs_mod, ONLY : um2cable_rr !packing
USE cable_um_init_subrs_mod, ONLY : init_veg_pars_fr_vegin
use cbl_LAI_canopy_height_mod, ONLY  : limit_HGT_LAI
!subrs:HaC1.3 - new subroutines
use cable_wide_mod,          ONLY : allocate_cable_wide
USE radiation_albedo_mod,    ONLY : allocate_rad_albedo
USE cable_pft_params_mod,    ONLY : cable_pft_params
USE cable_soil_params_mod,    ONLY : cable_soil_params
use cbl_allocate_types_mod,  ONLY : alloc_cbl_types
use cbl_masks_mod,              ONLY : fveg_mask, fsunlit_mask
use cbl_masks_mod,              ONLY : fsunlit_veg_mask

!data: constants                                        
USE cable_other_constants_mod,  ONLY : z0surf_min
USE cable_other_constants_mod,  ONLY : CLAI_THRESH => lai_thresh
USE cable_other_constants_mod,  ONLY : Ccoszen_tols => coszen_tols
USE cable_other_constants_mod,  ONLY : CGAUSS_W => gauss_w
USE cable_math_constants_mod,   ONLY : CPI => pi
USE cable_math_constants_mod,   ONLY : CPI180 => pi180

!data:HaC1.3
USE cbl_masks_mod,              ONLY : L_tile_pts
use cbl_masks_mod,              ONLY : veg_mask,  sunlit_mask,  sunlit_veg_mask 
USE cable_wide_mod,             ONLY : LAI_pft_cbl, HGT_pft_cbl
use cable_wide_mod,             ONLY : SnowDepth, SnowDensity
use cable_wide_mod,             ONLY : reducedLAIdue2snow
use cable_wide_mod,             ONLY : HeightAboveSnow
USE radiation_albedo_mod,       ONLY : ExtCoeff_beam, ExtCoeff_dif, &
                                       EffExtCoeff_beam, EffExtCoeff_dif, &
                                       CanopyRefl_dif, CanopyRefl_beam, &
                                       CanopyTransmit_dif, CanopyTransmit_beam,&
                                       coszen,        &
                                       VegXfang, VegTaul, VegRefl, c1, rhoch, &
                                       SW_down, RadFbeam, xk, AlbSnow, &
                                       AlbSoil, RadAlbedo, &
                                       EffSurfRefl_dif, EffSurfRefl_beam
                                        
!we need only 4 of these - can at least veg(soil)in%get through JaC 
use cbl_allocate_types_mod, ONLY :  veg, soil

implicit none

!JaC:todo:***Hack: get from jules
integer, parameter :: mp =1
integer, parameter :: nrb =3

integer :: metDoY(mp)          !local dummy Day of the Year [formerly met%doy]
!--- IN ARGS FROM surf_couple_radiation() ------------------------------------
#include "cbl_rad_main_decs.inc"

real :: surf_down_sw_NIR(row_length, rows)
real :: surf_down_sw_VIS(row_length, rows)

REAL :: MetTk(mp) 
REAL :: SnowODepth(mp)
REAL :: SoilTemp(mp)
REAL :: SnowAge(mp)
integer:: SnowFlag_3L(mp)
integer:: surface_type(mp) 

!--- declare vars local to CABLE -------------------------------------------- 
!packed in pack
!integer :: isnow_flg_cable(land_pts, ntiles)
LOGICAL, save :: um_online = .FALSE. 
LOGICAL, save :: jls_standalone = .FALSE. 
LOGICAL, save :: jls_radiation = .FALSE. !um_radiation = jls_radiation
  
!--- declare local vars to subr -------------------------------------------- 
character(len=*), parameter :: subr_name = "cable_rad_main"
logical, save :: first_call = .true.
integer :: i, j, n
LOGICAL :: skip =.TRUE. 

!--- initialize runtime switches !JaCP: 
jls_standalone = .TRUE.  !jls_standalone can be passed as TRUE from control.F90
jls_radiation = .TRUE.

metDoY = 0

IF(first_call) call alloc_cbl_types (mp)

!******************************************************************************
!***** JaC will deliver these to us: ******************************************
!******************************************************************************
! get veg(soil) parameters to enable calling rad on first timestep
IF(first_call ) call cable_pft_params()
IF(first_call)  CALL init_veg_pars_fr_vegin( dzsoil ) 
IF(first_call ) call cable_soil_params()
  
! define variables common cross CABLE (or at least >2 pathways)
IF(first_call) call allocate_cable_wide( mp )
! define variables common to rad/albedo pathway
IF(first_call) call allocate_rad_albedo( mp, nrb )

!------------------------------------------------------------------------------
!Pack forcing and conditions for CABLE
!------------------------------------------------------------------------------
!JaCP:we can grab these veg parameters from JaC later
VegXfang = Veg%Xfang
VegTaul = Veg%Taul
VegRefl = Veg%Refl

!convert to integer from D1 where everything is a float 
!isnow_flg_cable = int(snow_flag_cable)

!pack surface type: L_tile_pts mask set FROM surf_couple_rad
surface_type = PACK( tile_frac, L_tile_pts )

!--- CABLE soil albedo forcing only needs to be done on first call
!At present only single value is used for each land point 
if (first_call) then
  albsoil(:,2) =0.; albsoil(:,3) =0.
  CALL um2cable_lp( soil_alb, soil_alb, albsoil(:,1), soil%isoilm, skip )
endif

!Pack zenith this timestep
CALL um2cable_rr( cosine_zenith_angle, coszen)

!Pack SW
if( jls_standalone )  then !Total SW is forced for jls_standalone 
  !Jhan:Assume VIS:NIR is 1:1
  surf_down_sw_NIR =  0.5 * surf_down_sw(:,:,1)
  CALL um2cable_rr( surf_down_sw_NIR, SW_down(:,1) )
  
  surf_down_sw_VIS =  0.5 * surf_down_sw(:,:,1)
  CALL um2cable_rr( surf_down_sw_VIS, SW_down(:,2) )
endif

!define logical masks that are used throughout
call fveg_mask( mp, reducedLAIdue2snow )
call fsunlit_mask( mp, nrb, Ccoszen_tols, coszen )
call fsunlit_veg_mask( mp, veg_mask, sunlit_mask )

!Store Snow Depth from previous timestep. Treat differently on 1st timestep 
SnowODepth = SnowDepth 
SnowDepth = PACK( SNOW_TILE, L_TILE_PTS )
IF(first_call) SnowODepth = SnowDepth 

SnowDensity = PACK( SNOW_avg_rho_cable, L_TILE_PTS )

!Treat snow depth across 3Layers? Typecasts from Real to integer
SnowFlag_3L = PACK( SNOW_flag_cable, L_TILE_PTS )

!Surface skin/ top layer Snow  temperature
SoilTemp =   PACK( SOIL_temp_cable(:,:,1), L_TILE_PTS )

! limit IN height, LAI  and initialize existing cable % types
call limit_HGT_LAI( LAI_pft_cbl, HGT_pft_cbl, mp, land_pts, ntiles, &
                    tile_pts, tile_index, tile_frac, L_tile_pts, &
                    LAI_pft_um, HGT_pft_um )

!------------------------------------------------------------------------------
! Call CABLE_rad_driver to run specific and necessary components of CABLE 
!------------------------------------------------------------------------------
call cable_rad_driver(  &
#                     include "cbl_rad_driver_args.inc"
                     )

! Unpack variables (CABLE computed albedos) to JULES 
!------------------------------------------------------------------------------
call cable_rad_unpack(  &
#                      include "cbl_rad_unpack_args.inc"
                     )

!flick switches before leaving  
jls_radiation= .FALSE.
first_call = .false.

return

End subroutine cable_rad_main

End module cable_rad_main_mod

