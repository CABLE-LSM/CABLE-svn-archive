! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************

MODULE surf_couple_radiation_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: surf_couple_radiation

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SURF_COUPLE_RADIATION_MOD'

CONTAINS

!===============================================================================
! Public subroutine
!===============================================================================
SUBROUTINE surf_couple_radiation(                                             &
  !Fluxes INTENT(IN)
  tstar,                                                                      &
  !Misc INTENT(IN)
  ws10m, chloro,                                                              &
  n_band, max_n_swbands, wavelength_short, wavelength_long,                   &
  !Misc INTENT(OUT)
  sea_ice_albedo,                                                             &
  !Fluxes INTENT(OUT)
  alb_surft, land_albedo_ij,                                                  &
  !INTENT(IN)
  pond_frac_cat_sicat, pond_depth_cat_sicat,                                  &
  !(ancil_info mod)
  nsurft, land_pts, land_index, surft_pts, surft_index,                       &
  row_length, rows, ice_fract_ij, ice_fract_ncat_sicat, frac_surft,           &
  !(p_s_parms mod)
  cosz_ij, albobs_sw_gb, albobs_vis_gb, albobs_nir_gb,                        &
  z0_surft, albsoil_soilt,                                                    &
  !(coastal mod)
  flandg, tstar_sice_sicat,                                                   &
  !(prognostics mod)
  snow_mass_ij, snow_mass_sea_sicat, di_ncat_sicat, lai_pft, canht_pft,       &
  rgrain_surft, snow_surft, soot_ij, tstar_surft, ho2r2_orog_gb,              &
  !INTENT(OUT)
  albobs_sc_ij, open_sea_albedo,                                              &
  !CABLE_LSM:-----------------------------------------------------------------
  !Mostly model dimensions and associated!------------------------------------
  sm_levels,           & !grid cell number of soil levels 
  !Surface descriptions generally parametrized!-------------------------------
  dzsoil,              & !soil thicknesses in each layer  
  !CABLE dimensions !---------------------------------------------------------
  mp_cbl,             &!# CABLE vars assume vector of length mp(=Active patch)
  msn_cbl,          &!# of snow layers. at present=3 
  nrb_cbl,            &!# of rad. bands(confused b/n VIS/NIR, dir/dif. wrongly=3
  L_tile_pts_cbl,     &!Logical mask. TRUE where tile frac > 0. else = FALSE
  !introduced prognostics. tiled soil on 6 layers. tiled snow on 3 layers etc!
  SoilTemp_CABLE,           &!soil temperature (IN for rad.)
  SnowTemp_CABLE,           &!snow temperature (IN for rad.) REDUNDANT
  ThreeLayerSnowFlag_CABLE, &!flag signalling 3 layer treatment (binary) IN only
  OneLyrSnowDensity_CABLE,  &
  !constants!-----------------------------------------------------------------
  z0surf_min_cbl,           &
  lai_thresh_cbl,           &
  coszen_tols_cbl,          &
  gauss_w_cbl,              &
  pi_cbl,                   &
  pi180_cbl,                &
  !Vegetation parameters!-----------------------------------------------------
  SurfaceTypeID_cbl,        & 
  VegXfang,                 &
  VegTaul,                  &
  VegRefl                   &
)

!Module imports
USE ftsa_mod,        ONLY: ftsa
USE tile_albedo_mod, ONLY: tile_albedo

!Common modules
USE missing_data_mod,         ONLY:                                           &
  rmdi
USE ereport_mod,              ONLY:                                           &
  ereport

USE ancil_info,               ONLY:                                           &
  nsoilt

USE jules_sea_seaice_mod,     ONLY:                                           &
  nice, nice_use,                                                             &
  alpham, alphac, alphab, dtice, dt_bare, dalb_bare_wet,                      &
  pen_rad_frac, sw_beta,                                                      &
  albicev_cice, albicei_cice, albsnowv_cice, albsnowi_cice,                   &
  albpondv_cice, albpondi_cice,                                               &
  ahmax, dalb_mlt_cice, dalb_mlts_v_cice, dalb_mlts_i_cice,                   &
  dt_bare_cice, dt_snow_cice, pen_rad_frac_cice, sw_beta_cice,                &
  snowpatch
USE jules_surface_types_mod,  ONLY:                                           &
  ntype, npft

!Potential troublemaker
USE theta_field_sizes,        ONLY:                                           &
  t_i_length, t_j_length

!for testing LSM switch
USE jules_print_mgr,          ONLY: jules_message, jules_print

USE lsm_switch_mod,         ONLY:                                             &
  lsm_id, jules, cable

!CABLE_LSM:
USE cable_rad_main_mod, ONLY : cable_rad_main
!Dr Hook
USE parkind1,                 ONLY:                                           &
  jprb, jpim

USE yomhook,                  ONLY:                                           &
  lhook, dr_hook


IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Coupling routine between the UM or JULES system code and land surface
!   radiation science routines. Calls the appropriate LSM-specific code.
!
!   Some variables exist in modules only in JULES, others only in the UM
!   Options in order of preference
!   -UM and JULES share the same module names.
!   -UM and JULES have different module names and USE statements go on an ifdef
!   -The UM flavour of the variable does not live in a module. Pass in using a
!    an ifdef'ed argument list
!
!   If there are lots of ifs, then we could cosider splitting lsm_couple into
!   jules_couple and lsm_couple
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------


! Subroutine arguments

! Dimensioning variables
INTEGER, INTENT(IN) ::                                                        &
  n_band,                                                                     &
  max_n_swbands

  !UM-only args: INTENT(IN)
  !(ancil_info mod)
INTEGER, INTENT(IN)::                                                         &
  nsurft, land_pts, land_index(land_pts), surft_pts(nsurft),                  &
  surft_index(land_pts,nsurft), row_length, rows

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  ice_fract_ij(row_length, rows),                                             &
  ice_fract_ncat_sicat(row_length, rows,nice_use),                            &
  frac_surft(land_pts, ntype)

!(p_s_parms mod)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  cosz_ij(row_length, rows),                                                  &
  albobs_sw_gb(land_pts),                                                     &
  albobs_vis_gb(land_pts),                                                    &
  albobs_nir_gb(land_pts),                                                    &
  z0_surft(land_pts,nsurft),                                                  &
  albsoil_soilt(land_pts,nsoilt)

!(orog mod)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  ho2r2_orog_gb(land_pts)

!(coastal mod)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
flandg(row_length, rows),                                                     &
tstar_sice_sicat(row_length, rows)

!(prognostics mod)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  snow_mass_ij(row_length, rows),                                             &
  snow_mass_sea_sicat(row_length, rows, nice_use),                            &
  di_ncat_sicat(row_length, rows, nice_use),                                  &
  canht_pft(land_pts, npft),                                                  &
  lai_pft(land_pts, npft),                                                    &
  rgrain_surft(land_pts, nsurft),                                             &
  snow_surft(land_pts, nsurft),                                               &
  soot_ij(row_length, rows),                                                  &
  tstar_surft(land_pts, nsurft)

!UM-only args: INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  albobs_sc_ij(t_i_length,t_j_length,ntype,2),                                &
    !albedo scaling factors to obs
  open_sea_albedo(row_length,rows,2,max_n_swbands)
    !Surface albedo for Open Sea (direct and diffuse components, for each
    !band, with zeros for safety where no value applies)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  pond_frac_cat_sicat(row_length, rows, nice),                                &
!     Meltpond fraction on sea ice categories
     pond_depth_cat_sicat(row_length, rows, nice)
!     Meltpond depth on sea ice categories (m)

!Fluxes INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  tstar(row_length,rows)            !Surface temperature

!Misc INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  ws10m(row_length,rows),                                                     &
                                    !10m wind speed
  chloro(row_length,rows)           !nr surface chlorophyll content

REAL(KIND=real_jlslsm), INTENT(IN)    ::                                      &
  wavelength_short(n_band),                                                   &
  wavelength_long(n_band)

!Misc INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  sea_ice_albedo(row_length,rows,4)   !Surface Albedo for sea ice
                                      ! (*,1) - direct beam visible
                                      ! (*,2) - diffuse visible
                                      ! (*,3) - direct beam near-ir
                                      ! (*,4) - diffuse near-ir

!Fluxes INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  alb_surft(land_pts,nsurft,4),                                               &
                                      !Albedos for surface tiles.
                                      ! (*,*,1) - Direct beam visible
                                      ! (*,*,2) - Diffuse visible
                                      ! (*,*,3) - Direct beam near-IR
                                      ! (*,*,4) - Diffuse near-IR
  land_albedo_ij(t_i_length,t_j_length,4) !GBM albedos.
!CABLE_LSM: passed from control
integer :: sm_levels
real :: dzsoil(sm_levels)
integer :: mp_cbl
integer :: msn_cbl
integer :: nrb_cbl
LOGICAL :: L_tile_pts_cbl(land_pts,nsurft)
real :: SoilTemp_CABLE(land_pts, nsurft, sm_levels )
real :: SnowTemp_CABLE(land_pts, nsurft, msn_cbl )
real :: ThreeLayerSnowFlag_CABLE(land_pts, nsurft )
real :: OneLyrSnowDensity_CABLE(land_pts, nsurft )
!constants
real :: z0surf_min_cbl                  !the minimum roughness of bare soil
real :: lai_thresh_cbl                 !The minimum LAI below which a "cell" is considred NOT vegetated
real :: coszen_tols_cbl                !threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: gauss_w_cbl(nrb_cbl)               !Gaussian integration weights
real :: pi_cbl                         !PI - describing the ratio of circumference to diameter
real :: pi180_cbl                      !PI in radians
integer :: SurfaceTypeID_cbl(mp_cbl)
real :: VegXfang(mp_cbl)
real :: VegTaul(mp_cbl, nrb_cbl)
real :: VegRefl(mp_cbl, nrb_cbl)

!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------

!Land point only versions of ij variables
REAL(KIND=real_jlslsm) :: soot_gb(land_pts)
REAL(KIND=real_jlslsm) :: cosz_gb(land_pts)

!Counters
INTEGER :: i,j,l,band
INTEGER :: errorstatus

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SURF_COUPLE_RADIATION'

!-----------------------------------------------------------------------------
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE( lsm_id )
CASE ( jules )
  CALL ftsa (                                                                 &
    !INTENT(IN)
    !input fields
    flandg, ice_fract_ij, tstar, tstar_sice_sicat,                            &
    cosz_ij, ws10m, chloro,                                                   &
    snow_mass_sea_sicat, di_ncat_sicat,                                       &
    pond_frac_cat_sicat, pond_depth_cat_sicat,                                &
    !max and min sea ice albedo specifications
    alpham, alphac, alphab, dtice,                                            &
    dt_bare, dalb_bare_wet, pen_rad_frac, sw_beta,                            &
    ! parameters for CICE multi-band albedo scheme:
    albicev_cice, albicei_cice, albsnowv_cice, albsnowi_cice,                 &
    albpondv_cice, albpondi_cice,                                             &
    ahmax, dalb_mlt_cice, dalb_mlts_v_cice, dalb_mlts_i_cice,                 &
    dt_bare_cice, dt_snow_cice,                                               &
    pen_rad_frac_cice, sw_beta_cice, snowpatch,                               &
    !size and control variables
    row_length * rows, max_n_swbands,                                         &
    n_band, nice, nice_use,                                                   &
    !spectral boundaries
    wavelength_short,                                                         &
    wavelength_long,                                                          &
    !INTENT(OUT)
    !output arguments
    sea_ice_albedo,                                                           &
    open_sea_albedo)

  !Compress gridded variables to land point only
  DO l = 1, land_pts
    j = (land_index(l) - 1) / row_length + 1
    i = land_index(l) - (j-1) * row_length
    soot_gb(l)    = soot_ij(i,j)
    cosz_gb(l)    = cosz_ij(i,j)
  END DO

  CALL tile_albedo (                                                          &
    !INTENT(IN)
    t_i_length * t_j_length,                                                  &
    land_pts, nsurft,                                                         &
    land_index, surft_pts, surft_index,                                       &
    albsoil_soilt, albobs_sw_gb, albobs_vis_gb, albobs_nir_gb,                &
    cosz_gb, soot_gb, ho2r2_orog_gb,                                          &
    lai_pft, canht_pft,                                                       &
    rgrain_surft, snow_surft, tstar_surft, z0_surft, frac_surft,              &
    !INTENT(OUT)
    alb_surft,albobs_sc_ij,land_albedo_ij)

CASE ( cable )
#if defined(UM_JULES)
  errorstatus = 101
  CALL ereport('surf_couple_radiation', errorstatus,                          &
               'CABLE not yet implemented')
#else
  ! for testing LSM

  ! initialise all INTENT(OUT) fields for now until CABLE is implemented
  sea_ice_albedo(:,:,:) = 0.0
  alb_surft(:,:,:) = 0.0
  land_albedo_ij(:,:,:) = 0.0
  CALL ftsa (                                                                 &
    !INTENT(IN)
    !input fields
    flandg, ice_fract_ij, tstar, tstar_sice_sicat,                            &
    cosz_ij, ws10m, chloro,                                                   &
    snow_mass_sea_sicat, di_ncat_sicat,                                       &
    pond_frac_cat_sicat, pond_depth_cat_sicat,                                &
    !max and min sea ice albedo specifications
    alpham, alphac, alphab, dtice,                                            &
    dt_bare, dalb_bare_wet, pen_rad_frac, sw_beta,                            &
    ! parameters for CICE multi-band albedo scheme:
    albicev_cice, albicei_cice, albsnowv_cice, albsnowi_cice,                 &
    albpondv_cice, albpondi_cice,                                             &
    ahmax, dalb_mlt_cice, dalb_mlts_v_cice, dalb_mlts_i_cice,                 &
    dt_bare_cice, dt_snow_cice,                                               &
    pen_rad_frac_cice, sw_beta_cice, snowpatch,                               &
    !size and control variables
    row_length * rows, max_n_swbands,                                         &
    n_band, nice, nice_use,                                                   &
    !spectral boundaries
    wavelength_short,                                                         &
    wavelength_long,                                                          &
    !INTENT(OUT)
    !output arguments
    sea_ice_albedo,                                                           &
    open_sea_albedo)

  !Compress gridded variables to land point only
  DO l = 1, land_pts
    j = (land_index(l) - 1) / row_length + 1
    i = land_index(l) - (j-1) * row_length
    soot_gb(l)    = soot_ij(i,j)
    cosz_gb(l)    = cosz_ij(i,j)
  END DO
call cable_rad_main( & 
!corresponding name (if differs) of varaible on other side of call/subroutine shown in "[]" 

  !Variables to be calculated and returned by CABLE
  land_albedo_ij,   & ! GridBoxMean albedo per rad band (row_length,rows,4) [land_albedo]
  alb_surft,        & ! albedo per rad band per tile (land_pts, ntiles, 4) [alb_tile] 
  
  !Mostly model dimensions and associated
  row_length,          & !grid cell x
  rows,                & !grid cell y
  land_pts,            & !grid cell land points on the x,y grid
  nsurft,              & !grid cell number of surface types [ntiles] 
  !sm_levels,           & !grid cell number of soil levels 
  npft,                & !grid cell number of PFTs 
  surft_pts,           & !Number of land points per PFT [surft_pts] 
  surft_index,         & !Index of land point in (land_pts) array[surft_index] 
  land_index,          & !Index of land points in (x,y) array - see  corresponding *decs.inc
  
  !Surface descriptions generally parametrized
  !dzsoil,              & !soil thicknesses in each layer  
  frac_surft,          & !fraction of each surface type per land point [tile_frac] 
  LAI_pft,             & !Leaf area index. [LAI_pft_um]
  canht_pft,           & !Canopy height [HGT_pft_um]
  albsoil_soilt(:,1),  & !(albsoil)Snow-free, bare soil albedo [soil_alb]
  
  !Variables passed from control() level 
  snow_surft,          & ! snow depth equivalent (in water?) [snow_tile]
                         !This is the total snow depth per tile. CABLE also has depth per layer
  cosz_ij,             & ! cosine_zenith_angle [cosine_zenith_angle]  
  !The args below are passed from control() level as they usually do not exist
  !in the JULES rasiation pathway -------------------------------------------------
  !Mostly model dimensions and associated!------------------------------------
  sm_levels,           & !grid cell number of soil levels 
  !Surface descriptions generally parametrized!-------------------------------
  dzsoil,              & !soil thicknesses in each layer  
  !CABLE dimensions !---------------------------------------------------------
  mp_cbl,             &!# CABLE vars assume vector of length mp(=Active patch)
  msn_cbl,            &!# of snow layers. at present=3 
  nrb_cbl,            &!# of rad. bands(confused b/n VIS/NIR, dir/dif. wrongly=3
  L_tile_pts_cbl,     &!Logical mask. TRUE where tile frac > 0. else = FALSE
  !introduced prognostics. tiled soil on 6 layers. tiled snow on 3 layers etc!
  SoilTemp_CABLE,           &!soil temperature (IN for rad.)
  SnowTemp_CABLE,           &!snow temperature (IN for rad.) REDUNDANT
  ThreeLayerSnowFlag_CABLE, &!flag signalling 3 layer treatment (binary) IN only
  OneLyrSnowDensity_CABLE,  &
  !constants!-----------------------------------------------------------------
  z0surf_min_cbl,           &
  lai_thresh_cbl,           &
  coszen_tols_cbl,          &
  gauss_w_cbl,              &
  pi_cbl,                   &
  pi180_cbl,                &
  !Vegetation parameters!-----------------------------------------------------
  SurfaceTypeID_cbl,        & 
  VegXfang,                 &
  VegTaul,                  &
  VegRefl                   &
)
#endif

CASE DEFAULT
  errorstatus = 101
  CALL ereport('surf_couple_radiation', errorstatus,                          &
               'Unrecognised surface scheme')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE surf_couple_radiation
END MODULE surf_couple_radiation_mod

