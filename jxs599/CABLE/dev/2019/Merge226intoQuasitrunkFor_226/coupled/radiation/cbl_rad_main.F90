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
SurfaceType,        & 
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
integer:: SurfaceType(mp_cbl) 
real :: VegXfang(mp_cbl)
real :: VegTaul(mp_cbl, nrb_cbl)
real :: VegRefl(mp_cbl, nrb_cbl)

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
!are these local to CABLE and can be flushed every timestep
REAL :: ExtCoeff_beam(mp_cbl)
REAL :: ExtCoeff_dif(mp_cbl)
REAL :: EffExtCoeff_beam(mp_cbl, nrb_cbl)
REAL :: EffExtCoeff_dif(mp_cbl, nrb_cbl)

REAL :: CanopyTransmit_dif(mp_cbl, nrb_cbl)
REAL :: CanopyTransmit_beam(mp_cbl, nrb_cbl)
REAL :: CanopyRefl_dif(mp_cbl, nrb_cbl)
REAL :: CanopyRefl_beam(mp_cbl, nrb_cbl)

REAL :: EffSurfRefl_dif(mp_cbl, nrb_cbl)
REAL :: EffSurfRefl_beam(mp_cbl, nrb_cbl)

REAL :: coszen(mp_cbl)
REAL :: RadFbeam(mp_cbl, nrb_cbl)
REAL :: RadAlbedo(mp_cbl, nrb_cbl)
REAL :: AlbSnow(mp_cbl, nrb_cbl)
REAL :: AlbSoil(mp_cbl, nrb_cbl)

!co-efficients usoughout init_radiation ` called from _albedo as well
REAL :: c1(mp_cbl, nrb_cbl)
REAL :: rhoch(mp_cbl, nrb_cbl)
REAL :: xk(mp_cbl, nrb_cbl)
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
REAL :: SoilTemp(mp_cbl)
REAL :: SnowAge(mp_cbl)
integer:: SnowFlag_3L(mp_cbl)

!--- declare vars local to CABLE -------------------------------------------- 
!packed in pack
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
integer,save :: iradcall =0
logical, save :: albflip=.FALSE.
real, save :: ialb_surft

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

!******************************************************************************
!from rose-app.conf!!canht_ft_io= !!lai_io=
LAI_pft_um(1, :) = (/4.0,5.0,0.0,0.0,0.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/) 
HGT_pft_um(1,:) = (/16.38,19.01,0.0,0.0,0.0,0.79,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)

!--- initialize runtime switches !JaCP: 
jls_standalone = .TRUE.  !jls_standalone can be passed as TRUE from control.F90
jls_radiation = .TRUE.
!--- initialize/zero each timestep 
ExtCoeff_beam(:) = 0.0
ExtCoeff_dif(:) = 0.0
EffExtCoeff_beam(:,:) = 0.0
EffExtCoeff_dif(:,:) = 0.0
CanopyTransmit_dif(:,:) = 0.0
CanopyTransmit_beam(:,:) = 0.0
CanopyRefl_dif(:,:) = 0.0
CanopyRefl_beam(:,:) = 0.0
EffSurfRefl_dif(:,:) = 0.0
EffSurfRefl_beam(:,:) = 0.0
!coszen(:) = 0.0
RadFbeam(:,:) = 0.0
RadAlbedo(:,:) = 0.0
AlbSnow(:,:) = 0.0
AlbSoil(:,:) = 0.0
c1(:,:) = 0.0
rhoch(:,:) = 0.0
xk(:,:) = 0.0

metDoY = 0
IF(first_call) call alloc_cbl_types (mp)

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
                    LAI_pft_um, HGT_pft_um, CLAI_thresh )

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

iradcall = iradcall +1
!do i=1, Land_pts
!do j=1, nsurft
!  if(.NOT. albflip) then
!    if( alb_surft(i,j,1) > 0.0 ) then
!      albflip=.TRUE.
!      ialb_surft = alb_surft(i,j,1)
!!      iland_pts= Land_pts
!!      insurft = nsurft
!    endif
!  endif
!enddo
!enddo

if(iradcall>4104) then
!if(albflip) then
!  do i=1, Land_pts
!  do j=1, nsurft
!    if( alb_surft(Land_pts,nsurft,1) <=  ialb_surft) then
      print *, "iradcall ", iradcall
!    endif
!  enddo
!  enddo
endif

!flick switches before leaving  
jls_radiation= .FALSE.
first_call = .false.

return

End subroutine cable_rad_main

End module cable_rad_main_mod

