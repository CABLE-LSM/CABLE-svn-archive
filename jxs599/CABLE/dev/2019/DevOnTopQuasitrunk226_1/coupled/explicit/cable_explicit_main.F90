module cable_explicit_main_mod
  
contains

SUBROUTINE cable_explicit_main(                                                &
            timestep_width,                     &
            cycleno, numcycles,                                                &
            land_pts, ntiles,                                &
            npft, sm_levels,                                                   &
            land_index, tile_frac, tile_pts, tile_index,                       &
            snow_tile, lw_down, & 
            surf_down_sw, &
            tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, canopy_tile,           &
            Fland, CO2_MMR, sthu, canht_ft, lai_ft ,                           &
            FTL_TILE, FQW_TILE,                    &
            TSTAR_TILE, U_S, U_S_STD_TILE, CD_TILE, CH_TILE,                   &
            RADNET_TILE, FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,              &
            RECIP_L_MO_TILE, EPOT_TILE, doy &
            )
  !subrs called 
  USE cable_explicit_driv_mod, ONLY : cable_explicit_driver
  USE cable_expl_unpack_mod, ONLY : cable_expl_unpack

  !diag 
  USE cable_fprint_module, ONLY : cable_fprintf
  USE cable_Pyfprint_module, ONLY : cable_Pyfprintf
  USE cable_fFile_module, ONLY : fprintf_dir_root, fprintf_dir, L_cable_fprint,&
                                 L_cable_Pyfprint, unique_subdir

!data !H!
USE cable_other_constants_mod, ONLY : z0surf_min

  !processor number, timestep number / width, endstep
  USE cable_common_module, ONLY : knode_gl, ktau_gl, kwidth_gl, kend_gl
  USE cable_common_module, ONLY : cable_runtime
  USE cable_data_module, ONLY : cable
!jhan:this looks like I was testing standalone

!JULES5.3
!C!USE parallel_mod,   ONLY : mype => task_id 
USE timestep_mod,   ONLY : rtimestep => timestep

USE atm_fields_bounds_mod, ONLY : tdims

USE jules_soil_mod, ONLY: dzsoil

use model_grid_mod, ONLY : latitude, longitude

USE p_s_parms,      ONLY : bexp   => bexp_soilt,                               &
                           hcon   => hcon_soilt,                               &
                           satcon => satcon_soilt,                             &
                           sathh  => sathh_soilt,                              &
                           smvcst => smvcst_soilt,                             & 
                           smvcwt => smvcwt_soilt,                             &
                           smvccl => smvccl_soilt,                             &
                           albsoil =>albsoil_soilt,                            &
                           cosine_zenith_angle =>  cosz_ij 

!Imports for driving and flux variables
USE forcing, ONLY : ls_rain => ls_rain_ij, &
                    ls_snow => ls_snow_ij

!cable progs are set here
USE allocate_cable_progs_mod, ONLY :               &
  SoilTemp        =>  SoilTemp_CABLE,             &
  SoilMoisture    =>  SoilMoisture_CABLE,         &
  FrozenSoilFrac  => FrozenSoilFrac_CABLE,     &
  SnowDepth   => SnowDepth_CABLE,               &
  SnowMass    => SnowMass_CABLE,                &
  SnowTemp    => SnowTemp_CABLE,                &
  SnowDensity => SnowDensity_CABLE,             &
  SnowAge     => SnowAge_CABLE,                 &
  ThreeLayerSnowFlag  => ThreeLayerSnowFlag_CABLE,      &
  OneLyrSnowDensity   => OneLyrSnowDensity_CABLE

  USE cbl_allocate_types_mod, ONLY : air, bgc, canopy,      &
                                met, bal, rad, rough, soil, ssnow, sum_flux,  &
                                veg
!data
USE cable_types_mod, ONLY : L_tile_pts
 
  implicit none
 
!JaC:todo:***Hack: get from jules
integer, parameter :: mp =1
integer, parameter :: nrb=3

  !___ re-decl input args
  integer :: endstep, cycleno, numcycles
  real :: timestep_width
  INTEGER :: row_length, rows

  INTEGER ::                                                      & 
    land_pts,         & ! # of land points being processed
    ntiles,           & ! # of tiles 
    npft,             & ! # of plant functional types
    sm_levels          ! # of soil layers 

# if defined(UM_JULES)
  REAL,  DIMENSION(tdims%i_end,tdims%j_end) :: latitude,  longitude
# endif   
  INTEGER, DIMENSION(land_pts) :: land_index  ! index of land points processed

  INTEGER,  DIMENSION(ntiles) :: tile_pts ! # of land points on each tile

  INTEGER,  DIMENSION(land_pts, ntiles) :: tile_index ! index of tile points 

  REAL, DIMENSION(land_pts, ntiles) :: tile_frac
   
    REAL,  DIMENSION(land_pts, ntiles) :: snow_tile
   
  REAL,  DIMENSION(tdims%i_end,tdims%j_end) ::                                         &
  lw_down!,               &

  REAL, DIMENSION(tdims%i_end, tdims%j_end, 4) ::                         &
    surf_down_sw 
  
  REAL,  DIMENSION(tdims%i_end,tdims%j_end) ::                             &
    tl_1,       &
    qw_1,       &  
    vshr_land,  &
    pstar,      &
    z1_tq,      &
    z1_uv

  REAL, DIMENSION(land_pts, ntiles) :: canopy_tile
  
  REAL, DIMENSION(land_pts) :: fland
  
  REAL :: co2_mmr

  REAL, DIMENSION(land_pts, sm_levels) :: sthu 
   
  REAL, DIMENSION(land_pts, npft) :: canht_ft, lai_ft 
  
  REAL :: sin_theta_latitude(tdims%i_end,tdims%j_end) 
    
  !___return miscelaneous 
  REAL, DIMENSION(land_pts,ntiles) :: &
    FTL_TILE,    & ! Surface FTL for land tiles     
    FQW_TILE,    & ! Surface FQW for land tiles     
    TSTAR_TILE,  & ! radiative temperature of surface 
    U_S,         & ! Surface friction velocity (m/s)
    U_S_STD_TILE,& ! Surface friction velocity
    CD_TILE,     & ! Drag coefficient
    CH_TILE,     &  ! Transfer coefficient for heat & moisture
    RADNET_TILE, & ! Surface net radiation
    FRACA,       & ! Fraction of surface moisture
    RESFS,       & ! Combined soil, stomatal & aerodynamic resistance
                   ! factor for fraction (1-FRACA) of snow-free land tiles
    RESFT,       & ! Total resistance factor.
                   ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
                   ! 1 for snow.    
    Z0H_TILE,    &
    Z0M_TILE,    &
    RECIP_L_MO_TILE,& ! Reciprocal of the Monin-Obukhov length for tiles (m^-1)
    EPOT_TILE

  integer :: timestep_number  
  integer :: doy

  !___ local vars
  integer,  DIMENSION(land_pts, ntiles) :: iThreeLayerSnowFlag
  logical, save :: first_call = .true.
  real :: radians_degrees

  ! std template args 
  character(len=*), parameter :: subr_name = "cable_explicit_main"

  row_length = tdims%i_end
  rows = tdims%j_end
  
  !--- initialize cable_runtime% switches 
  cable_runtime%um =          .true.
  cable_runtime%um_explicit = .TRUE.
   
  ! initialize processor number, timestep width & number, endstep 
  ! UM overwrites these defaults. Not yet done in StandAlone 
  if( first_call ) then
    knode_gl =0; kwidth_gl = 1200.; kend_gl=-1
# if defined(UM_JULES)
    !knode_gl  = mype 
    kwidth_gl = int(timestep_width)
    kend_gl   = endstep   
# endif

    !--- Convert lat/long to degrees
    radians_degrees = 180.0 / ( 4.0*atan(1.0) ) ! 180 / PI
    latitude  = latitude * radians_degrees
    longitude  = longitude * radians_degrees
  endif

  timestep_number = int(rtimestep)
  ktau_gl   = timestep_number
  
  !----------------------------------------------------------------------------
  !--- CALL _driver to run specific and necessary components of CABLE with IN -
  !--- args PACKED to force CABLE
  !----------------------------------------------------------------------------
  iThreeLayerSnowFlag= int(ThreeLayerSnowFlag)

  call cable_explicit_driver(  &
#                                include "cbl_explicit_driver_args.inc"
                            ) 
  
  !----------------------------------------------------------------------------
  !--- CALL _unpack to unpack variables from CABLE back to UM format to return
  !----------------------------------------------------------------------------
   call cable_expl_unpack( latitude, longitude, FTL_TILE, FQW_TILE, TSTAR_TILE, &
                           U_S, U_S_STD_TILE, &
                           CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts,             &
                           ssnow%snowd, ssnow%cls, air%rlam, air%rho,          &
                           canopy%fe, canopy%fh, canopy%us, canopy%cdtq,       &
                           canopy%fwet, canopy%wetfac_cs, canopy%rnet,         &
                           canopy%zetar, canopy%epot, met%ua, rad%trad,        &
                           rad%transd, rough%z0m, rough%zref_tq, &
                           canopy%fes, canopy%fev )

  cable_runtime%um_explicit = .FALSE.
  first_call = .false.        

return

End subroutine cable_explicit_main
  
End module cable_explicit_main_mod

