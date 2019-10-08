module cable_rad_main_mod
  
contains

SUBROUTINE cable_rad_main( &
!corresponding name (if differs) of varaible on other side of call/subroutine shown in "[]" 

!Variables to be calculated and returned by CABLE
!------------------------------------------------------------------------------
land_albedo,   & ! GridBoxMean albedo per rad band (row_length,rows,4) [land_albedo_ij]
alb_surft,     & ! albedo per rad band per tile (land_pts, nsurft, 4) [alb_tile] 
!------------------------------------------------------------------------------

!Mostly model dimensions and associated
row_length,          & !grid cell x
rows,                & !grid cell y
land_pts,            & !grid cell land points on the x,y grid
nsurft,              & !grid cell number of surface types [nsurft] 
!sm_levels,           & !grid cell number of soil levels 
npft,                & !grid cell number of PFTs 
tile_pts,            & !Number of land points per PFT [surft_pts] 
tile_index,          & !Index of land point in (land_pts) array[surft_index] 
land_index,          & !Index of land points in (x,y) array - see  corresponding *decs.inc
!------------------------------------------------------------------------------

!Surface descriptions generally parametrized
!------------------------------------------------------------------------------
!dzsoil,              & !soil thicknesses in each layer  
tile_frac,           & !fraction of each surface type per land point [frac_surft] 
LAI_pft_um,          & !Leaf area index. [LAI_pft]
HGT_pft_um,          & !Canopy height [canht_pft]
soil_alb,            & !(albsoil)Snow-free, bare soil albedo [albsoil_soilt(:,1)]
!------------------------------------------------------------------------------

!Variables passed from JULES/UM
!------------------------------------------------------------------------------
snow_tile,           & !snow depth equivalent (in water?) [snow_surft]
                       !This is the total snow depth per tile. CABLE also has depth per layer
cosine_zenith_angle, & ! cosine_zenith_angle [cosz_ij]
!------------------------------------------------------------------------------
  !The args below are passed from control() level as they usually do not exist
  !in the JULES rasiation pathway -------------------------------------------------
!Mostly model dimensions and associated!---------------------------------------
sm_levels,           & !grid cell number of soil levels 

!Surface descriptions generally parametrized!----------------------------------
dzsoil,              & !soil thicknesses in each layer  

!CABLE dimensions !------------------------------------------------------------
mp_cbl,             &!# CABLE vars assume vector of length mp(=Active patch)
msn_cable,          &!# of snow layers. at present=3 
nrb_cbl,            &!# of rad. bands(confused b/n VIS/NIR, dir/dif. wrongly=3
L_tile_pts_cbl,     &!Logical mask. TRUE where tile frac > 0. else = FALSE

!introduced prognostics. tiled soil on 6 layers. tiled snow on 3 layers etc!---
SoilTemp_CABLE,           &!soil temperature (IN for rad.)
SnowTemp_CABLE,           &!snow temperature (IN for rad.) REDUNDANT
ThreeLayerSnowFlag_CABLE, &!flag signalling 3 layer treatment (binary) IN only
OneLyrSnowDensity_CABLE,  &

!constants!--------------------------------------------------------------------
z0surf_min_cbl,           &
lai_thresh_cbl,           &
coszen_tols_cbl,          &
gauss_w_cbl,              &
pi_cbl,                   &
pi180_cbl,                &

!Vegetation parameters!--------------------------------------------------------
VegXfang,                 &
VegTaul,                  &
VegRefl                   &
) 
 
!subrs
USE cable_rad_driv_mod,   ONLY : cable_rad_driver
USE cable_rad_unpack_mod, ONLY : cable_rad_unpack
USE cable_common_module,  ONLY : cable_runtime
USE cable_pack_mod,       ONLY : cable_pack_lp !packing
USE cable_pack_mod,       ONLY : cable_pack_rr !packing

!subrs:HaC1.3 
use cbl_LAI_canopy_height_mod, ONLY  : limit_HGT_LAI

!subrs:HaC1.3 - new subroutines
USE radiation_albedo_mod,    ONLY : allocate_rad_albedo
USE radiation_albedo_mod,    ONLY : deallocate_rad_albedo

use cbl_allocate_types_mod,  ONLY : alloc_cbl_types

!we need only 4 of these - can at least veg(soil)in%get through JaC 
use cbl_allocate_types_mod, ONLY :  veg, soil

implicit none

!--- IN ARGS FROM surf_couple_radiation() ------------------------------------
integer :: row_length                       !grid cell x
integer :: rows                             !grid cell y
integer :: land_pts                         !grid cell land points on the x,y grid
integer :: nsurft                           !grid cell number of surface types 
integer :: npft                             !grid cell number of PFTs 
!!integer :: sm_levels                        !grid cell number of soil levels 
!!real :: dzsoil(sm_levels)                   !soil thicknesses in each layer  

!Variables to be calculated and returned by CABLE
real :: land_albedo(row_length,rows,4)      
real :: alb_surft(Land_pts,nsurft,4)        

real :: tile_frac(land_pts,nsurft)          !fraction of each surface type per land point 
integer :: tile_pts(nsurft)                 !Number of land points per PFT 
integer :: tile_index(land_pts,nsurft)      !Index of land point in (land_pts) array
integer :: land_index(land_pts)             !Index of land points in (x,y) array - see below
real :: LAI_pft_um(land_pts, npft)          !Leaf area index.
real :: HGT_pft_um(land_pts, npft)          !Canopy height

real :: snow_tile(land_pts,nsurft)          ! snow depth equivalent (in water?)
real :: soil_alb(land_pts)                  !(albsoil)Snow-free, bare soil albedo
real :: cosine_zenith_angle(row_length,rows)             !cosine_zenith_angle          

!--- IN ARGS FROM control() ------------------------------------
integer :: sm_levels
real :: dzsoil(sm_levels)
integer :: mp_cbl
integer :: msn_cable
integer :: nrb_cbl
LOGICAL :: L_tile_pts_cbl(land_pts,nsurft)
real :: SoilTemp_CABLE(land_pts, nsurft, sm_levels )
real :: SnowTemp_CABLE(land_pts, nsurft, msn_cable )
real :: ThreeLayerSnowFlag_CABLE(land_pts, nsurft )
real :: OneLyrSnowDensity_CABLE(land_pts, nsurft )
!constants
real :: z0surf_min_cbl                  !the minimum roughness of bare soil
real :: lai_thresh_cbl                 !The minimum LAI below which a "cell" is considred NOT vegetated
real :: coszen_tols_cbl                !threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: gauss_w_cbl(nrb_cbl)               !Gaussian integration weights
real :: pi_cbl                         !PI - describing the ratio of circumference to diameter
real :: pi180_cbl                      !PI in radians
real :: VegXfang(mp_cbl)
real :: VegTaul(mp_cbl, nrb_cbl)
real :: VegRefl(mp_cbl, nrb_cbl)
!real :: SnowDepth(mp_cbl)             !Formerly: ssnow%snowd 
!real :: SnowDensity(mp_cbl)           !Formerly: ssnow%ssdnn 
!
!real :: SnowODepth(mp_cbl)             !Formerly: ssnow%osnowd
! 
!!computed from UM HT(LAI)_PFT passed in explicit call - need at rad call
!real :: LAI_pft_cbl(mp_cbl)           !Formerly: ~veg%vlai
!real :: HGT_pft_cbl(mp_cbl)           !Formerly:  ~veg%hc 
!
!!can compute from  z0surf_min, HGT_pft_cbl, SnowDepth, SnowDensity )
!real :: HeightAboveSnow(mp_cbl)       !Formerly: rough%hruff


integer :: metDoY(mp_cbl)          !local dummy Day of the Year [formerly met%doy]
!!below decs to rename IN args **for now**
integer :: mp
integer :: nrb 
LOGICAL :: L_tile_pts(land_pts,nsurft)                                               
real :: soil_temp_cable( land_pts, nsurft, sm_levels )   !Mean snow density
real :: snow_temp_cable(land_pts, nsurft, msn_cable)             !Snow temperature (3 layer)
real :: snow_flag_cable(land_pts,nsurft)                 !surface type fraction 
real :: snow_avg_rho_cable(LAND_PTS,nsurft)              !snow cover in the ice tile.

!constants
real :: Ccoszen_tols            !threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: Cgauss_w(nrb_cbl)               !Gaussian integration weights
real :: Clai_thresh                 !The minimum LAI below which a "cell" is considred NOT vegetated
real :: Cpi                         !PI - describing the ratio of circumference to diameter
real :: Cpi180                      !PI in radians
real :: z0surf_min                  !the minimum roughness of bare soil

!Convoluted mapping using land_index(array) to get back to the row_length,rows co-ordinate
! J = ( LAND_INDEX(L)-1 ) / ROW_LENGTH + 1
! I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
! FTL_1(I,J) = where ftl_1(row_length,rows) 
!-------------------------------------------------------------------------------
REAL, ALLOCATABLE :: ExtCoeff_beam(:)
REAL, ALLOCATABLE :: ExtCoeff_dif(:)
REAL, ALLOCATABLE :: EffExtCoeff_beam(:,:)
REAL, ALLOCATABLE :: EffExtCoeff_dif(:,:)

REAL, ALLOCATABLE :: CanopyTransmit_dif(:,:)  
REAL, ALLOCATABLE :: CanopyTransmit_beam(:,:)  
REAL, ALLOCATABLE :: CanopyRefl_dif(:,:)  
REAL, ALLOCATABLE :: CanopyRefl_beam(:,:)  

REAL, ALLOCATABLE :: EffSurfRefl_dif(:,:)  
REAL, ALLOCATABLE :: EffSurfRefl_beam(:,:)  

!REAL, ALLOCATABLE :: reducedLAIdue2snow(:)
REAL, ALLOCATABLE :: coszen(:)
!REAL, ALLOCATABLE :: VegXfang(:)
!REAL, ALLOCATABLE :: VegTaul(:,:)
!REAL, ALLOCATABLE :: VegRefl(:,:)
!REAL, ALLOCATABLE :: metDoY(:)
!REAL, ALLOCATABLE :: SW_down(:,:)  
REAL, ALLOCATABLE :: RadFbeam(:,:)  
REAL, ALLOCATABLE :: RadAlbedo(:,:)  
REAL, ALLOCATABLE :: AlbSnow(:,:)  
REAL, ALLOCATABLE :: AlbSoil(:,:)  

!co-efficients usoughout init_radiation ` called from _albedo as well
REAL, ALLOCATABLE :: c1(:,:)
REAL, ALLOCATABLE :: rhoch(:,:)
REAL, ALLOCATABLE :: xk(:,:)     ! extinct. coef.for beam rad. and black leaves
!-------------------------------------------------------------------------------
!already passed at rad call 
real :: SnowDepth(mp_cbl)             !Formerly: ssnow%snowd 
real :: SnowDensity(mp_cbl)           !Formerly: ssnow%ssdnn 
real :: SnowODepth(mp_cbl)             !Formerly: ssnow%osnowd
!computed from UM HT(LAI)_PFT passed in explicit call - need at rad call
real :: LAI_pft_cbl(mp_cbl)           !Formerly: ~veg%vlai
real :: HGT_pft_cbl(mp_cbl)           !Formerly:  ~veg%hc 
!can compute from  z0surf_min, HGT_pft_cbl, SnowDepth, SnowDensity )
real :: HeightAboveSnow(mp_cbl)       !Formerly: rough%hruff

REAL :: MetTk(mp_cbl) 
!REAL :: SnowODepth(mp_cbl)
REAL :: SoilTemp(mp_cbl)
REAL :: SnowAge(mp_cbl)
integer:: SnowFlag_3L(mp_cbl)
integer:: surface_type(mp_cbl) 

!--- declare vars local to CABLE -------------------------------------------- 
!packed in pack
!integer :: isnow_flg_cable(land_pts, nsurft)
LOGICAL, save :: um_online = .FALSE. 
LOGICAL, save :: jls_standalone = .FALSE. 
LOGICAL, save :: jls_radiation = .FALSE. !um_radiation = jls_radiation
  
!make local to rad_driver and also again in cbl_model_driver
!CABLE variables to keep for all CABLE pathways across the timestep 
real :: reducedLAIdue2snow(mp_cbl)
!masks
logical :: veg_mask(mp_cbl),  sunlit_mask(mp_cbl),  sunlit_veg_mask(mp_cbl) 

!--- declare local vars to subr -------------------------------------------- 
character(len=*), parameter :: subr_name = "cable_rad_main"
logical, save :: first_call = .true.
integer :: i, j, n
LOGICAL :: skip =.TRUE. 
real :: surf_down_sw_NIR(row_length,rows)                  !1of2-band ShortWave forcing
real :: surf_down_sw_VIS(row_length,rows)                  !2of2-band ShortWave forcing
real :: SW_down(mp_cbl,2)

!******************************************************************************
!******************************************************************************
!******************************************************************************
!******************************************************************************
!******************************************************************************
!******************************************************************************

!--- initialize runtime switches !JaCP: 
jls_standalone = .TRUE.  !jls_standalone can be passed as TRUE from control.F90
jls_radiation = .TRUE.

metDoY = 0
!re-name locally 
mp = mp_cbl
nrb = nrb_cbl 
L_tile_pts = L_tile_pts_cbl
soil_temp_cable = SoilTemp_CABLE
snow_temp_cable = SnowTemp_CABLE
snow_flag_cable = ThreeLayerSnowFlag_CABLE
snow_avg_rho_cable = OneLyrSnowDensity_CABLE

z0surf_min   = z0surf_min_cbl
Ccoszen_tols = coszen_tols_cbl 
Cgauss_w     = gauss_w_cbl
Clai_thresh  = lai_thresh_cbl
Cpi          = pi_cbl
Cpi180       = pi180_cbl

IF(first_call) call alloc_cbl_types (mp)

! allocate variables common to rad/albedo pathway
call allocate_rad_albedo( mp, nrb, ExtCoeff_beam, ExtCoeff_dif, &
EffExtCoeff_beam, EffExtCoeff_dif, &
CanopyRefl_dif, CanopyRefl_beam, &
CanopyTransmit_dif, CanopyTransmit_beam,&
coszen,        &
!VegXfang, VegTaul, VegRefl, 
c1, rhoch, &
RadFbeam, xk, AlbSnow, &
RadAlbedo,AlbSoil, &
EffSurfRefl_dif, EffSurfRefl_beam &
)
!------------------------------------------------------------------------------
!Pack forcing and conditions for CABLE
!------------------------------------------------------------------------------
!JaCP:we can grab these veg parameters from JaC later
!VegXfang = Veg%Xfang
!VegTaul = Veg%Taul
!VegRefl = Veg%Refl

!convert to integer from D1 where everything is a float 
!isnow_flg_cable = int(snow_flag_cable)

!pack surface type: L_tile_pts mask set FROM surf_couple_rad
surface_type = PACK( tile_frac, L_tile_pts )

!At present only single value is used for each land point 
albsoil(:,2) =0.; albsoil(:,3) =0.
CALL cable_pack_lp( soil_alb, soil_alb, albsoil(:,1), soil%isoilm, skip )

!Pack zenith this timestep
CALL cable_pack_rr( cosine_zenith_angle, coszen)

!Pack SW
if( jls_standalone )  then !Total SW is forced for jls_standalone 
  !Jhan:dummies. on rad call dont actually need
  surf_down_sw_NIR =  0.0 
  CALL cable_pack_rr( surf_down_sw_NIR, SW_down(:,1) )
  surf_down_sw_VIS =  0.0
  CALL cable_pack_rr( surf_down_sw_VIS, SW_down(:,2) )
endif

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
call limit_HGT_LAI( LAI_pft_cbl, HGT_pft_cbl, mp, land_pts, nsurft, &
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
!flush cable
! DEallocate variables common to rad/albedo pathway
call deallocate_rad_albedo( mp, nrb, ExtCoeff_beam, ExtCoeff_dif, &
EffExtCoeff_beam, EffExtCoeff_dif, &
CanopyRefl_dif, CanopyRefl_beam, &
CanopyTransmit_dif, CanopyTransmit_beam,&
coszen,        &
!VegXfang, VegTaul, VegRefl, 
c1, rhoch, &
RadFbeam, xk, AlbSnow, &
RadAlbedo,AlbSoil, &
EffSurfRefl_dif, EffSurfRefl_beam &
)
!flick switches before leaving  
jls_radiation= .FALSE.
first_call = .false.

return

End subroutine cable_rad_main

End module cable_rad_main_mod

