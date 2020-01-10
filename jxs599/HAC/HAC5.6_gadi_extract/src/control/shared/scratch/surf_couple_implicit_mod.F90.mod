! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************

MODULE surf_couple_implicit_mod

USE sf_impl2_mod, ONLY: sf_impl2
USE sf_impl2_cbl_mod, ONLY: sf_impl2_cbl
USE sice_htf_mod, ONLY: sice_htf

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: surf_couple_implicit

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SURF_COUPLE_IMPLICIT_MOD'

CONTAINS

!===============================================================================
! Public subroutine
!===============================================================================
SUBROUTINE surf_couple_implicit(                                              &
  !Important switch
  l_correct,                                                                  &
  !Forcing INTENT(IN)
  pstar, lw_down, qw_1, tl_1, u_1, v_1, u_0, v_0,                             &
  !Fluxes INTENT(IN)
  sw_surft, emis_surft,                                                       &
  !Misc INTENT(IN) Many of these simply come out of explicit and into here.
  rhokm_u, rhokm_v, gamma1, gamma2, alpha1, alpha1_sea, alpha1_sice,          &
  ashtf, ashtf_sea, ashtf_surft, du, dv, fraca, resfs, resft, rhokh,          &
  rhokh_surft, rhokh_sice, rhokh_sea,                                         &
  z0hssi, z0mssi, z0h_surft, z0m_surft, chr1p5m,                              &
  chr1p5m_sice, canhc_surft, flake, tile_frac, wt_ext_surft,                  &
  cdr10m_u, cdr10m_v, r_gamma,                                                &
  !Diagnostics, INTENT(INOUT)
  sf_diag,                                                                    &
  !Fluxes INTENT(INOUT)
  fqw_ice, ftl_ice, fqw_surft, fqw_1, ftl_1, ftl_surft,                       &
  !Misc INTENT(INOUT)
  epot_surft, dtstar_surft, dtstar_sea, dtstar_sice, radnet_sice, olr,        &
  !Fluxes INTENT(OUT)
  tstar, le_surft, radnet_surft, e_sea, h_sea, taux_1, tauy_1, ecan_surft, ei,&
  esoil_soilt, ext_soilt, snowmelt, melt_surft,                               &
  ecan, ei_surft, esoil_surft, sea_ice_htf, surf_ht_flux, surf_htf_surft,     &
  !Misc INTENT(OUT)
  error,                                                                      &
  !UM-only arguments
  !JULES ancil_info module
  nsurft, land_pts, land_index, surft_index, surft_pts, ice_fract_ij,         &
  sstfrz_ij, ice_fract_ncat_sicat, z1_tq_ij,                                  &
  !JULES prognostics module
  canopy_surft, smc_soilt, k_sice_sicat, t_soil_soilt, ti_sicat, snow_surft,  &
  di_ncat_sicat, tstar_surft,                                                 &
  !JULES coastal module
  fland, flandg, tstar_sea_ij, tstar_sice_sicat, tstar_ssi_ij,                &
  taux_land_ij, tauy_land_ij, taux_ssi_ij, tauy_ssi_ij,                       &
  surf_ht_flux_land_ij, surf_ht_flux_sice_sicat, tstar_land_ij, tstar_sice_ij,&
  !JULES u_v_grid module
  dtrdz_charney_grid_1_ij,                                                    &
  !JULES switches module
  !==l_mr_physics in UM
  l_co2_interactive, lq_mix_bl,                                               &
  !JULES aero module
  co2_3d_ij,                                                                  &
  !Arguments without a JULES module
  ctctq1,dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,flandg_u,flandg_v, &
  rho1, f3_at_p, uStarGBM,tscrndcl_ssi,tscrndcl_surft,tStbTrans,              &
  taux_land_star,tauy_land_star,taux_ssi_star,                                &
  tauy_ssi_star,ei_sice,rhokh_mix, ti_gb,                                     &
!CABLE_LSM: 
  surf_down_sw )

!Module Imports

!Common modules
USE ereport_mod,           ONLY: ereport

USE jules_soil_mod,        ONLY: sm_levels
USE jules_sea_seaice_mod,  ONLY: nice_use, nice, l_sice_multilayers
USE sf_diags_mod,          ONLY: strnewsfdiag
USE ancil_info,            ONLY: nsoilt
!Potential troublemakers
USE atm_fields_bounds_mod, ONLY: tdims,   udims,   vdims,   pdims,            &
                                 tdims_s, udims_s, vdims_s, pdims_s

! Module switches name between UM and JULES-standalone
#if defined(UM_JULES)
USE rad_input_mod, ONLY: co2_mmr
#else
USE aero, ONLY:                                                               &
  co2_mmr
#endif

USE lsm_switch_mod, ONLY:                                                     &
  lsm_id, jules, cable

!CABLE_LSM:USE additional existing vars and nml defined lsm_id
!UM modules
#if defined(UM_JULES)
USE atmos_physics2_alloc_mod, only : gc
#endif
 
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

! for testing LSM switch
USE jules_print_mgr,          ONLY: jules_message, jules_print

IMPLICIT NONE


!-----------------------------------------------------------------------------
! Description:
!   Coupling routine between the UM or JULES system code and land surface
!   implicit science routines. Calls the appropriate LSM-specific code.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!Arguments

!UM-only arguments
!JULES ancil_info module
INTEGER, INTENT(IN) ::                                                        &
  nsurft, land_pts,                                                           &
  land_index(land_pts),                                                       &
  surft_index(nsurft),                                                        &
  surft_pts(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  ice_fract_ij(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
  ice_fract_ncat_sicat(pdims%i_start:pdims%i_end,                             &
                       pdims%j_start:pdims%j_end,nice),                       &
  sstfrz_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
  z1_tq_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!JULES prognostics module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  canopy_surft(land_pts,nsurft),                                              &
  smc_soilt(land_pts,nsoilt),                                                 &
  k_sice_sicat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice),     &
  t_soil_soilt(land_pts,nsoilt,sm_levels),                                    &
  snow_surft(land_pts,nsurft),                                                &
  di_ncat_sicat(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice)
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  ti_sicat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice),         &
  tstar_surft(land_pts,nsurft)
!JULES coastal module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  fland(land_pts),                                                            &
  flandg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  tstar_sice_sicat(tdims%i_start:tdims%i_end,                                 &
    tdims%j_start:tdims%j_end,nice_use),                                      &
  tstar_sea_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
  tstar_ssi_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
  taux_land_ij(udims%i_start:udims%i_end,udims%j_start:udims%j_end),          &
  tauy_land_ij(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),          &
  taux_ssi_ij(udims%i_start:udims%i_end,udims%j_start:udims%j_end),           &
  tauy_ssi_ij(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  surf_ht_flux_land_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),  &
  surf_ht_flux_sice_sicat(tdims%i_start:tdims%i_end,                          &
    tdims%j_start:tdims%j_end,nice),                                          &
  tstar_land_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
  tstar_sice_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!JULES u_v_grid module
REAL(KIND=real_jlslsm), INTENT(IN)  ::                                        &
  dtrdz_charney_grid_1_ij(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
!JULES switches module
LOGICAL, INTENT(IN) ::                                                        &
  l_co2_interactive, lq_mix_bl !==l_mr_physics
!JULES aero module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  co2_3d_ij(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end)

!Arguments without a JULES module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  ctctq1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
  dqw1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
  dtl1_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
  du_star1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),      &
  dv_star1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end),      &
  cq_cm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end),             &
    ! Coefficient in U tri-diagonal implicit matrix
  cq_cm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),             &
    ! Coefficient in V tri-diagonal implicit matrix
  flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),              &
    !Land frac (on U-grid, with 1st and last rows undefined or, at present,
    !set to "missing data")
  flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),              &
    !Land frac (on V-grid, with 1st and last rows undefined or, at present,
    !set to "missing data")
  rho1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Density on lowest level
  f3_at_p(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Coriolis parameter
  ustargbm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
    ! BM surface friction velocity
  tscrndcl_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    !Decoupled screen-level temperature over sea or sea-ice
  tscrndcl_surft(land_pts,nsurft),                                            &
    !Decoupled screen-level temperature over land tiles
  tstbtrans(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Time since the transition
  taux_land_star(udims%i_start:udims%i_end,udims%j_start:udims%j_end),        &
  taux_ssi_star(udims%i_start:udims%i_end,udims%j_start:udims%j_end),         &
  tauy_land_star(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),        &
  tauy_ssi_star(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  ei_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Sea ice sublimation
  rhokh_mix(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Exchange coeffs for moisture.
  ti_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

! Important switch
LOGICAL, INTENT(IN) ::                                                        &
  l_correct                             ! flag used by the new BL solver

!Forcing INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
    !Surface pressure (Pascals).
  lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface downward LW radiation (W/m2).
  qw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Total water content
  tl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Ice/liquid water temperature
  u_1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),           &
    !W'ly wind component (m/s)
  v_1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end),           &
    !S'ly wind component (m/s)
  u_0(udims%i_start:udims%i_end,udims%j_start:udims%j_end),                   &
    !W'ly component of surface current (m/s).
  v_0(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
    !S'ly component of surface current (m/s).

!Fluxes INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  sw_surft(land_pts,nsurft),                                                  &
    !Surface net SW radiation on land tiles (W/m2).
  emis_surft(land_pts,nsurft)
    !Emissivity for land tiles

!Misc INTENT(IN) Many of these simply come out of explicit and come
!back into here. Need a module for them.

!Arrays
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),               &
    !Exchange coefficients for momentum (on U-grid, with 1st and last rows
    !undefined or, at present, set to "missing data")
  rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),               &
    !Exchange coefficients for momentum (on V-grid, with 1st and last rows
    !undefined or, at present, set to "missing data")
  gamma1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
    !weights for new BL solver
  gamma2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
  alpha1(land_pts,nsurft),                                                    &
    !Mean gradient of saturated specific humidity with respect to temperature
    !between the bottom model layer and tile surfaces
  alpha1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
    !ALPHA1 for sea.
  alpha1_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !ALPHA1 for sea-ice.
  ashtf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),        &
    !Adjusted SEB coefficient for sea-ice
  ashtf_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Adjusted SEB coefficient for sea
  ashtf_surft(land_pts,nsurft),                                               &
    !Adjusted SEB coefficient for land tiles.
  du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),            &
    !Level 1 increment to u wind field
  dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end),            &
    !Level 1 increment to v wind field
  fraca(land_pts,nsurft),                                                     &
    !Fraction of surface moisture flux with only aerodynamic resistance for
    !snow-free land tiles.
  resfs(land_pts,nsurft),                                                     &
    !Combined soil, stomatal and aerodynamic resistance factor for fraction
    !(1-FRACA) of snow-free land tiles.
  resft(land_pts,nsurft),                                                     &
    !Total resistance factor. FRACA+(1-FRACA)*RESFS for snow-free land, 1 for
    !snow.
  rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Grid-box surface exchange coefficients (not used for JULES)
  rhokh_surft(land_pts,nsurft),                                               &
    !Surface exchange coefficients for land tiles
  rhokh_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),   &
    !Surface exchange coefficients for sea sea-ice
  rhokh_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Surface exchange coefficients for sea
  z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
  z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    !Roughness lengths over sea (m)
  z0h_surft(land_pts,nsurft),                                                 &
    !Tile roughness lengths for heat and moisture (m).
  z0m_surft(land_pts,nsurft),                                                 &
    !Tile roughness lengths for momentum.
  chr1p5m(land_pts,nsurft),                                                   &
    !Ratio of coefffs for calculation of 1.5m temp for land tiles.
  chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    !CHR1P5M for sea and sea-ice (leads ignored).
  canhc_surft(land_pts,nsurft),                                               &
    !Areal heat capacity of canopy for land tiles (J/K/m2).
  flake(land_pts,nsurft),                                                     &
    !Lake fraction.
  tile_frac(land_pts,nsurft),                                                 &
    !Tile fractions including snow cover in the ice tile.
  wt_ext_surft(land_pts,sm_levels,nsurft),                                    &
    !Fraction of evapotranspiration extracted from each soil layer by each tile.
  cdr10m_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),              &
    !Ratio of CD's reqd for calculation of 10 m wind. On U-grid; comments as
    !per RHOKM.
  cdr10m_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
    !Ratio of CD's reqd for calculation of 10 m wind. On V-grid; comments as

! Implicit weighting for B.L.
REAL(KIND=real_jlslsm), INTENT(IN)  :: r_gamma

!diagnostic array
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag

!Fluxes INTENT(INOUT)
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  fqw_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Surface FQW for sea-ice
  ftl_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Surface FTL for sea-ice
  fqw_surft(land_pts,nsurft),                                                 &
    !Surface FQW for land tiles
  ftl_surft(land_pts,nsurft),                                                 &
    !Surface FTL for land tiles
  fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Moisture flux between layers (kg per square metre per sec) FQW(,1) is
    !total water flux from surface, 'E'.
  ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !FTL(,K) contains net turbulent sensible heat flux into layer K from below;
    !so FTL(,1) is the surface sensible heat, H.(W/m2)

!Misc INTENT(INOUT)
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  epot_surft(land_pts,nsurft),                                                &
    !surface tile potential evaporation
  dtstar_surft(land_pts,nsurft),                                              &
    !Change in TSTAR over timestep for land tiles
  dtstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
    !Change is TSTAR over timestep for open sea
  dtstar_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Change is TSTAR over timestep for sea-ice
  radnet_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Surface net radiation on sea-ice (W/m2)
  olr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    ! IN    TOA - surface upward LW on last radiation timestep
    ! OUT   Corrected TOA outward LW

!Fluxes INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !GBM surface temperature (K).
  le_surft(land_pts,nsurft),                                                  &
    !Surface latent heat flux for land tiles
  radnet_surft(land_pts,nsurft),                                              &
    !Surface net radiation on land tiles (W/m2)
  e_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Evaporation from sea times leads fraction. Zero over land.
    !(kg per square metre per sec).
  h_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Surface sensible heat flux over sea times leads fraction (W/m2)
  taux_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end),                &
    !W'ly component of surface wind stress (N/sq m). (On UV-grid with first
    !and last rows undefined or, at present, set to missing data
  tauy_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),                &
    !S'ly component of surface wind stress (N/sq m).  On UV-grid; comments as
    !per TAUX
  ecan_surft(land_pts,nsurft),                                                &
    !ECAN for snow-free land tiles
  ei(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                    &
    !Sublimation from lying snow or sea-ice (kg/m2/s).
  ei_surft(land_pts,nsurft),                                                  &
    !EI for land tiles.
  esoil_soilt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nsoilt),    &
    !Surface evapotranspiration from soil moisture store (kg/m2/s).
  ext_soilt(land_pts,nsoilt,sm_levels),                                       &
    !Extraction of water from each soil layer (kg/m2/s).
  snowmelt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
    !Snowmelt (kg/m2/s).
  melt_surft(land_pts,nsurft),                                                &
    !Snowmelt on land tiles (kg/m2/s)
  ecan(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Gridbox mean evaporation from canopy/surface store (kg/m2/s). Zero over sea
  esoil_surft(land_pts,nsurft),                                               &
    !ESOIL for snow-free land tiles
  sea_ice_htf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice),      &
    !Heat flux through sea-ice (W/m2, positive downwards). (Not used for JULES)
  surf_ht_flux(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    !Net downward heat flux at surface over land and sea-ice fraction of gridbox
    !(W/m2).
  surf_htf_surft(land_pts,nsurft)
    !Net downward surface heat flux on tiles (W/m2)

!Misc INTENT(OUT)
INTEGER, INTENT(OUT) ::                                                       &
  error          !0 - AOK; 1 to 7  - bad grid definition detected


!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
INTEGER :: i,j,n

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SURF_COUPLE_IMPLICIT'
!CABLE_LSM:
!passed from control
real  :: surf_down_sw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, 4)
!CABLE_LSM:additional existing vars
!CM2!REAL, dimension(land_pts) ::                                     &
!CM2!  gs
!CM2!
!CM2!REAL ::                                                     &
!CM2!   gpp(land_pts),                                                    &
!CM2!                                ! IN Gross primary productivity
!CM2!                                !    (kg C/m2/s).
!CM2!   gpp_ft(land_pts,nsurft),                                            &
!CM2!                                ! IN Gross primary productivity
!CM2!                                !    on PFTs (kg C/m2/s).
!CM2!   npp(land_pts),                                                    &
!CM2!                                ! IN Net primary productivity
!CM2!                                !    (kg C/m2/s).
!CM2!   npp_ft(land_pts,nsurft),                                            &
!CM2!                                ! IN Net primary productivity
!CM2!                                !    on PFTs (kg C/m2/s).
!CM2!   resp_p(land_pts),                                                 &
!CM2!                                ! IN Plant respiration (kg C/m2/s).
!CM2!   resp_p_ft(land_pts,nsurft),                                         &
!CM2!                                  !IN Plant respiration on PFTs
!CM2!                                !     (kg C/m2/s).
!CM2!   resp_s(land_pts,dim_cs1),                                         &
!CM2!                                   ! IN Soil respiration (kg C/m2/s).
!CM2!   resp_s_tot(dim_cs2),                                                 &
!CM2!                                   ! IN Total soil resp'n (kg C/m2/s).
!CM2!   g_leaf(land_pts,nsurft)
!CM2!                                ! IN Leaf turnover rate (/360days).
!CM2! 
!CABLE_LSM: End


!-----------------------------------------------------------------------------
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
SELECT CASE( lsm_id )
CASE ( jules )

  CALL sf_impl2 (                                                             &
  ! IN values defining field dimensions and subset to be processed :
          land_pts,land_index,nice,nice_use,nsurft,surft_index,surft_pts,     &
          sm_levels,canhc_surft,canopy_surft,flake,smc_soilt,tile_frac,       &
          wt_ext_surft, fland,flandg,lq_mix_bl,                               &
  ! IN sea/sea-ice data :
          ice_fract_ij,ice_fract_ncat_sicat,k_sice_sicat,u_0,v_0,             &
  ! IN everything not covered so far :
          pstar,lw_down,sw_surft,                                             &
          t_soil_soilt,qw_1,tl_1,u_1,v_1,rhokm_u,rhokm_v,r_gamma,             &
          gamma1,gamma2,alpha1,alpha1_sea,alpha1_sice,                        &
          ashtf,ashtf_sea,ashtf_surft,                                        &
          dtrdz_charney_grid_1_ij,du,dv,                                      &
          fraca,resfs,resft,                                                  &
          rhokh,rhokh_surft,rhokh_sice,rhokh_sea,z1_tq_ij,                    &
          z0hssi,z0mssi,z0h_surft,z0m_surft,                                  &
          cdr10m_u,cdr10m_v,                                                  &
          chr1p5m,chr1p5m_sice,ctctq1,                                        &
          dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,                &
          l_correct,flandg_u,flandg_v,                                        &
          emis_surft,ti_sicat,snow_surft,                                     &
  ! IN variables used to calculate cooling at the screen level
          l_co2_interactive, co2_mmr, co2_3d_ij,rho1, f3_at_p, ustargbm,      &
  ! INOUT data :
          epot_surft,fqw_ice,ftl_ice,dtstar_surft,dtstar_sea,dtstar_sice,     &
          tstar_sice_sicat,tstar_ssi_ij,tstar_surft,tstar_sea_ij,radnet_sice, &
          fqw_surft, fqw_1,ftl_1,ftl_surft,olr,taux_land_ij,taux_ssi_ij,      &
          tauy_land_ij,tauy_ssi_ij,                                           &
          tscrndcl_ssi,tscrndcl_surft,tStbTrans,sf_diag,                      &
          taux_land_star,tauy_land_star,taux_ssi_star, tauy_ssi_star,         &
  ! OUT Diagnostic not requiring STASH flags :
          ecan,ei_surft,esoil_surft,sea_ice_htf,surf_ht_flux,                 &
          surf_ht_flux_land_ij,surf_ht_flux_sice_sicat,surf_htf_surft,        &
  ! OUT data required elsewhere in UM system :
          tstar,tstar_land_ij,tstar_sice_ij,le_surft,radnet_surft,e_sea,h_sea,&
          taux_1,tauy_1,ecan_surft,ei,ei_sice,esoil_soilt,ext_soilt,snowmelt, &
          melt_surft,                                                         &
          rhokh_mix,error                                                     &
        )

          !Only call on the second pass through implicit
  IF ( l_correct) THEN
    IF ( .NOT. l_sice_multilayers) THEN
      !--------------------------------------------------------------------
      ! Update sea-ice surface layer temperature,
      ! if not coupled to multilayer sea ice model.
      !--------------------------------------------------------------------
      CALL sice_htf(                                                          &
        !IN fields
        flandg,nice,                                                          &
        di_ncat_sicat,ice_fract_ij,ice_fract_ncat_sicat,                      &
        surf_ht_flux_sice_sicat,sstfrz_ij,                                    &
        !INOUT fields
        ti_sicat,sf_diag,                                                     &
        !OUT fields
        ti_gb,sea_ice_htf                                                     &
        )
    ELSE
      sea_ice_htf(:,:,:) = 0.0   ! for safety
    END IF  !l_sice_multilayers

    ! Convert sea and sea-ice fluxes to be fraction of grid-box
    ! (as required by sea and sea-ice modellers)
    DO n = 1, nice
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          surf_ht_flux_sice_sicat(i,j,n) =                                    &
            ice_fract_ncat_sicat(i,j,n) * surf_ht_flux_sice_sicat(i,j,n)
          sea_ice_htf(i,j,n)       =                                          &
            ice_fract_ncat_sicat(i,j,n) * sea_ice_htf(i,j,n)
        END DO !i
      END DO !j
    END DO  !n

    IF (sf_diag%simlt) THEN
      DO n = 1, nice
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            sf_diag%sice_mlt_htf(i,j,n) = ice_fract_ncat_sicat(i,j,n) *       &
                                        sf_diag%sice_mlt_htf(i,j,n)
          END DO
        END DO
      END DO
    END IF

  END IF !l_correct

CASE ( cable )
#if defined(UM_JULES)
  error = 101
  CALL ereport('surf_couple_implicit', error, 'CABLE not yet implemented')
#else

  ! initialise all INTENT(OUT) for now until CABLE is implemented
  !H!tstar(:,:) = 0.0
  !H!le_surft(:,:) = 0.0
  !H!radnet_surft(:,:) = 0.0
  !H!e_sea(:,:) = 0.0
  !H!h_sea(:,:) = 0.0
  !H!taux_1(:,:) = 0.0
  !H!tauy_1(:,:) = 0.0
  !H!ecan_surft(:,:) = 0.0
  !H!ei(:,:) = 0.0
  !H!esoil_soilt(:,:,:) = 0.0
  !H!ext_soilt(:,:,:) = 0.0
  !H!snowmelt(:,:) = 0.0
  !H!melt_surft(:,:) = 0.0
  !H!ecan(:,:) = 0.0
  !H!ei_surft(:,:) = 0.0
  !H!esoil_surft(:,:) = 0.0
  !H!sea_ice_htf(:,:,:) = 0.0
  !H!surf_ht_flux(:,:) = 0.0
  !H!surf_htf_surft(:,:) = 0.0
  !H!error = 0
#endif

  CALL sf_impl2_cbl (                                                        &
  ! IN values defining field dimensions and subset to be processed :
          land_pts,land_index,nice,nice_use,nsurft,surft_index,surft_pts,     &
          sm_levels,canhc_surft,canopy_surft,flake,smc_soilt,tile_frac,       &
          wt_ext_surft, fland,flandg,lq_mix_bl,                               &
  ! IN sea/sea-ice data :
          ice_fract_ij,ice_fract_ncat_sicat,k_sice_sicat,u_0,v_0,             &
  ! IN everything not covered so far :
          pstar,lw_down,sw_surft,                                             &
          t_soil_soilt,qw_1,tl_1,u_1,v_1,rhokm_u,rhokm_v,r_gamma,             &
          gamma1,gamma2,alpha1,alpha1_sea,alpha1_sice,                        &
          ashtf,ashtf_sea,ashtf_surft,                                        &
          dtrdz_charney_grid_1_ij,du,dv,                                      &
          fraca,resfs,resft,                                                  &
          rhokh,rhokh_surft,rhokh_sice,rhokh_sea,z1_tq_ij,                    &
          z0hssi,z0mssi,z0h_surft,z0m_surft,                                  &
          cdr10m_u,cdr10m_v,                                                  &
          chr1p5m,chr1p5m_sice,ctctq1,                                        &
          dqw1_1,dtl1_1,du_star1,dv_star1,cq_cm_u_1,cq_cm_v_1,                &
          l_correct,flandg_u,flandg_v,                                        &
          emis_surft,ti_sicat,snow_surft,                                     &
  ! IN variables used to calculate cooling at the screen level
          l_co2_interactive, co2_mmr, co2_3d_ij,rho1, f3_at_p, ustargbm,      &
  ! INOUT data :
          epot_surft,fqw_ice,ftl_ice,dtstar_surft,dtstar_sea,dtstar_sice,     &
          tstar_sice_sicat,tstar_ssi_ij,tstar_surft,tstar_sea_ij,radnet_sice, &
          fqw_surft, fqw_1,ftl_1,ftl_surft,olr,taux_land_ij,taux_ssi_ij,      &
          tauy_land_ij,tauy_ssi_ij,                                           &
          tscrndcl_ssi,tscrndcl_surft,tStbTrans,sf_diag,                      &
          taux_land_star,tauy_land_star,taux_ssi_star, tauy_ssi_star,         &
  ! OUT Diagnostic not requiring STASH flags :
          ecan,ei_surft,esoil_surft,sea_ice_htf,surf_ht_flux,                 &
          surf_ht_flux_land_ij,surf_ht_flux_sice_sicat,surf_htf_surft,        &
  ! OUT data required elsewhere in UM system :
          tstar,tstar_land_ij,tstar_sice_ij,le_surft,radnet_surft,e_sea,h_sea,&
          taux_1,tauy_1,ecan_surft,ei,ei_sice,esoil_soilt,ext_soilt,snowmelt, &
          melt_surft,                                                         &
          rhokh_mix,error                                                     &
        )

          !Only call on the second pass through implicit
  IF ( l_correct) THEN
    IF ( .NOT. l_sice_multilayers) THEN
      !--------------------------------------------------------------------
      ! Update sea-ice surface layer temperature,
      ! if not coupled to multilayer sea ice model.
      !--------------------------------------------------------------------
      CALL sice_htf(                                                          &
        !IN fields
        flandg,nice,                                                          &
        di_ncat_sicat,ice_fract_ij,ice_fract_ncat_sicat,                      &
        surf_ht_flux_sice_sicat,sstfrz_ij,                                    &
        !INOUT fields
        ti_sicat,sf_diag,                                                     &
        !OUT fields
        ti_gb,sea_ice_htf                                                     &
        )
    ELSE
      sea_ice_htf(:,:,:) = 0.0   ! for safety
    END IF  !l_sice_multilayers

    ! Convert sea and sea-ice fluxes to be fraction of grid-box
    ! (as required by sea and sea-ice modellers)
    DO n = 1, nice
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          surf_ht_flux_sice_sicat(i,j,n) =                                    &
            ice_fract_ncat_sicat(i,j,n) * surf_ht_flux_sice_sicat(i,j,n)
          sea_ice_htf(i,j,n)       =                                          &
            ice_fract_ncat_sicat(i,j,n) * sea_ice_htf(i,j,n)
        END DO !i
      END DO !j
    END DO  !n

    IF (sf_diag%simlt) THEN
      DO n = 1, nice
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            sf_diag%sice_mlt_htf(i,j,n) = ice_fract_ncat_sicat(i,j,n) *       &
                                        sf_diag%sice_mlt_htf(i,j,n)
          END DO
        END DO
      END DO
    END IF

  END IF !l_correct


CASE DEFAULT
  error = 101
  CALL ereport('surf_couple_implicit', error, 'Unrecognised surface scheme')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE surf_couple_implicit

END MODULE surf_couple_implicit_mod

