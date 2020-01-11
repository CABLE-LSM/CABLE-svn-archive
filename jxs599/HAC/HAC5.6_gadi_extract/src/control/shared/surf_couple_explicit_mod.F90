! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************

MODULE surf_couple_explicit_mod

USE sf_expl_l_mod, ONLY: sf_expl_l

USE um_types, ONLY: real_jlslsm
IMPLICIT NONE

PRIVATE
PUBLIC :: surf_couple_explicit

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SURF_COUPLE_EXPLICIT_MOD'

CONTAINS

!===============================================================================
! Public subroutine
!===============================================================================
SUBROUTINE surf_couple_explicit(                                              &
    !Arguments used by JULES-standalone
    !Misc INTENT(IN) DONE
    bq_1, bt_1, zh, photosynth_act_rad,                                       &
    curr_year, curr_day_number, curr_hour, curr_minute, curr_second,          &
    !Forcing INTENT(IN) DONE
    qw, tl, pstar, lw_down,                                                   &
    !Fluxes INTENT(IN) DONE
    sw_surft, tstar,                                                          &
    !Diagnostics, INTENT(INOUT)
    sf_diag,                                                                  &
    !Fluxes INTENT(OUT) DONE
    fqw_1,ftl_1,ftl_surft, fqw_surft, fqw_ice, ftl_ice, fsmc, emis_surft,     &
    !Misc INTENT(OUT)
    radnet_sice, rhokm_1, rhokm_land, rhokm_ssi,                              &
    !Out of explicit and into implicit only INTENT(OUT)
    cdr10m,                                                                   &
    alpha1, alpha1_sea, alpha1_sice, ashtf_prime, ashtf_prime_sea,            &
    ashtf_prime_surft, epot_surft,                                            &
    fraca, resfs, resft, rhokh, rhokh_surft, rhokh_sice, rhokh_sea,           &
    dtstar_surft, dtstar_sea, dtstar_sice,                                    &
    z0hssi, z0h_surft, z0mssi, z0m_surft, chr1p5m, chr1p5m_sice, canhc_surft, &
    wt_ext_surft, flake,                                                      &
    !Out of explicit and into extra only INTENT(OUT)
    hcons_soilt,                                                              &
    !Out of explicit and into implicit and extra INTENT(OUT)
    tile_frac,                                                                &
    !Additional arguments for the UM-----------------------------------------
    !JULES prognostics module
    canopy_surft, snow_surft, k_sice_sicat, cs_pool_gb_um, canht_pft, lai_pft,&
    t_soil_soilt, tsurf_elev_surft, ti_sicat, ti_cat_sicat,                   &
    tstar_surft, z0msea_ij, smc_soilt, gc_surft, gs_gb,                       &
    !JULES ancil_info module
    land_pts, z1_uv_ij, z1_tq_ij, land_index, nsurft, ice_fract_ncat_sicat,   &
    frac_surft, surft_index, surft_pts,                                       &
    ! IN input data from the wave model
    charnock_w,                                                               &
    !JULES coastal module
    fland, flandg, tstar_sea_ij, tstar_sice_sicat, vshr_land_ij, vshr_ssi_ij, &
    !JULES aero module
    co2_3d_ij, rho_aresist_ij, aresist_ij, resist_b_ij, rho_aresist_surft,    &
    aresist_surft, resist_b_surft, r_b_dust_ij, cd_std_dust_ij, u_s_std_surft,&
    !JULES trifctl module
    asteps_since_triffid, g_leaf_acc_pft, npp_acc_pft, resp_w_acc_pft,        &
    resp_s_acc_gb_um, gpp_gb, npp_gb, resp_p_gb, g_leaf_pft, gpp_pft, npp_pft,&
    resp_p_pft, resp_s_gb_um, resp_w_pft,                                     &
    !JULES p_s_parms module
    catch_surft, catch_snow_surft, hcon_soilt, smvccl_soilt, smvcst_soilt,    &
    smvcwt_soilt, sthf_soilt, sthu_soilt, z0_surft,                           &
    z0h_bare_surft, z0m_soil_gb, albsoil_soilt, cosz_ij, soil_clay_ij,        &
    !JULES orog module
    ho2r2_orog_gb, sil_orog_land_gb, h_blend_orog_ij, z0m_eff_ij,             &
    !JULES u_v_grid module
    u_1_p_ij, v_1_p_ij, u_0_p_ij, v_0_p_ij,                                   &
    !JULES switches module
    l_spec_z0,                                                                &
    !JULES c_elevate module
    z_land_ij,                                                                &
    !Not in a JULES module
    numcycles, cycleno, z1_uv_top, z1_tq_top, ddmfx,                          &
    l_aero_classic, z0m_scm, z0h_scm, recip_l_mo_sea, rib,                    &
    flandfac, fseafac, fb_surf, u_s, t1_sd, q1_sd, rhostar,                   &
    vshr, resp_s_tot_soilt, emis_soil,                                        &
    !CABLE_LSM: bring 4-band  downward SW back to control (or toplevel) 
    surf_down_sw )

!Module imports

!Common modules
USE ereport_mod,              ONLY:                                           &
  ereport
USE jules_sea_seaice_mod,     ONLY:                                           &
  nice, nice_use
USE jules_soil_mod,           ONLY:                                           &
  sm_levels
USE jules_soil_biogeochem_mod, ONLY:                                          &
! imported scalar parameters
    soil_model_rothc,                                                         &
! imported scalar variables (IN)
    soil_bgc_model
USE jules_vegetation_mod,     ONLY:                                           &
  l_phenol
USE jules_surface_types_mod,  ONLY:                                           &
  npft, ntype
USE sf_diags_mod, ONLY: strnewsfdiag
USE ozone_vars,               ONLY:                                           &
  o3_gb
USE atm_fields_bounds_mod,    ONLY:                                           &
  pdims_s, pdims, tdims
USE dust_param,               ONLY:                                           &
  ndiv
USE trif_vars_mod,            ONLY:                                           &
  resp_l_pft, resp_r_pft, n_leaf_pft, n_stem_pft, n_root_pft, lai_bal_pft
USE ancil_info,               ONLY:                                           &
  nsoilt

!Modules that change name between JULES and UM
#if defined(UM_JULES)
USE atm_step_local,         ONLY:                                             &
  land_pts_trif, npft_trif,dim_cs1, dim_cs2, co2_dim_len,co2_dim_row
USE carbon_options_mod,     ONLY:                                             &
  l_co2_interactive
USE dust_parameters_mod,    ONLY:                                             &
  l_dust, l_dust_diag
!  USE s_main_force,           ONLY: & !A SCM module currently doesn't work
!    l_spec_z0
USE rad_input_mod,          ONLY:                                             &
  co2_mmr
USE gen_phys_inputs_mod,    ONLY:                                             &
  lq_mix_bl => l_mr_physics
USE ancil_info,          ONLY:                                                &
  dim_cslayer
#else
USE ancil_info,             ONLY:                                             &
  land_pts_trif, npft_trif,dim_cs1, dim_cs2, co2_dim_len,co2_dim_row,         &
  dim_cslayer
USE switches,               ONLY:                                             &
  l_co2_interactive, l_dust, lq_mix_bl
USE aero,                   ONLY:                                             &
  co2_mmr

USE trifctl,                ONLY:                                             &
  resp_s_acc_soilt, resp_s_soilt
USE prognostics,            ONLY:                                             &
  cs_pool_soilt
USE p_s_parms,              ONLY:                                             &
  clay_soilt
#endif

! reads the lsm switch set in a namelist
USE lsm_switch_mod,        ONLY: lsm_id, jules, cable

! for testing lsm switch
USE jules_print_mgr,          ONLY: jules_message, jules_print

!Dr Hook
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

!CABLE_LSM:pass additional existing vars from divergence point here (lsm_id)
use sf_expl_l_cbl_mod, ONLY : sf_expl_l_cbl

!CABLE_LSM:pass additional these vars from control 
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
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Coupling routine between the UM or JULES system code and land surface
!   explicit science routines. Calls the appropriate LSM-specific code.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Subroutine arguments

!UM-only arguments
!JULES ancil_info module
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
  land_index(land_pts),                                                       &
  nsurft
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  z1_uv_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
  z1_tq_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
  ice_fract_ncat_sicat(tdims%i_start:tdims%i_end,                             &
                       tdims%j_start:tdims%j_end,nice_use),                   &
  frac_surft(land_pts,ntype),                                                 &
  charnock_w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
INTEGER, INTENT(OUT)  ::                                                      &
  surft_index(land_pts,nsurft),                                               &
  surft_pts(nsurft)

!JULES prognostics module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  canopy_surft(land_pts,nsurft),                                              &
  snow_surft(land_pts,nsurft),                                                &
  k_sice_sicat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use), &
  cs_pool_gb_um(land_pts,dim_cs1),                                            &
  canht_pft(land_pts,npft),                                                   &
  lai_pft(land_pts,npft),                                                     &
  t_soil_soilt(land_pts,nsoilt,sm_levels),                                    &
  tsurf_elev_surft(land_pts,nsurft),                                          &
  ti_sicat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
  ti_cat_sicat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice),     &
  tstar_surft(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  z0msea_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
  gs_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  smc_soilt(land_pts,nsoilt),                                                 &
  gc_surft(land_pts,nsurft)
!JULES coastal module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  fland(land_pts),                                                            &
  flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
  tstar_sea_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
  tstar_sice_sicat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  vshr_land_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
  vshr_ssi_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!JULES aero module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  co2_3d_ij(co2_dim_len,co2_dim_row)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  rho_aresist_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
  aresist_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
  resist_b_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
  rho_aresist_surft(land_pts,nsurft),                                         &
  aresist_surft(land_pts,nsurft),                                             &
  resist_b_surft(land_pts,nsurft),                                            &
  r_b_dust_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,ndiv),      &
  cd_std_dust_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
  u_s_std_surft(land_pts,nsurft)
!JULES trifctl module
INTEGER, INTENT(IN)::                                                         &
  asteps_since_triffid
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  g_leaf_acc_pft(land_pts,npft),                                              &
  npp_acc_pft(land_pts_trif,npft_trif),                                       &
  resp_w_acc_pft(land_pts_trif,npft_trif),                                    &
  resp_s_acc_gb_um(land_pts_trif,dim_cs1)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  gpp_gb(land_pts),                                                           &
  npp_gb(land_pts),                                                           &
  resp_p_gb(land_pts),                                                        &
  g_leaf_pft(land_pts,npft),                                                  &
  gpp_pft(land_pts,npft),                                                     &
  npp_pft(land_pts,npft),                                                     &
  resp_p_pft(land_pts,npft),                                                  &
  resp_s_gb_um(land_pts,dim_cs1),                                             &
  resp_w_pft(land_pts,npft)
!JULES p_s_parms module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  catch_surft(land_pts,nsurft),                                               &
  catch_snow_surft(land_pts,nsurft),                                          &
  hcon_soilt(land_pts,nsoilt),                                                &
  smvccl_soilt(land_pts,nsoilt,sm_levels),                                    &
  smvcst_soilt(land_pts,nsoilt,sm_levels),                                    &
  smvcwt_soilt(land_pts,nsoilt,sm_levels),                                    &
  sthf_soilt(land_pts,nsoilt,sm_levels),                                      &
  sthu_soilt(land_pts,nsoilt,sm_levels),                                      &
  z0_surft(land_pts,nsurft),                                                  &
  z0h_bare_surft(land_pts,nsurft),                                            &
  z0m_soil_gb(land_pts),                                                      &
  albsoil_soilt(land_pts,nsoilt),                                             &
  cosz_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
  soil_clay_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!JULES orog module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  ho2r2_orog_gb(land_pts),                                                    &
  sil_orog_land_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  z0m_eff_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
  h_blend_orog_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!JULES u_v_grid module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  u_1_p_ij(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  v_1_p_ij(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  u_0_p_ij(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  v_0_p_ij(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
!JULES switches module
LOGICAL, INTENT(IN) ::                                                        &
  l_spec_z0
!JULES c_elevate module
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  z_land_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
!Not in a JULES module
INTEGER, INTENT(IN) :: numcycles, cycleno
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  !These variables are INTENT(IN) to sf_expl, but not used with the
  !current configuration of standalone JULES (initialised to 0 below)
  z1_uv_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest uv-layer
  z1_tq_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    ! Height of top of lowest Tq-layer
  ddmfx(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    ! Convective downdraught mass-flux at cloud base
LOGICAL, INTENT(IN) :: l_aero_classic
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  !These variables are required for prescribed roughness lengths in
  !SCM mode in UM - not used standalone
  z0m_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Fixed Sea-surface roughness length for momentum (m).(SCM)
  z0h_scm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !Fixed Sea-surface roughness length for heat (m). (SCM)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  recip_l_mo_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
    !Reciprocal of the surface Obukhov  length at sea points. (m-1).
  rib(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !Mean bulk Richardson number for lowest layer.
  flandfac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),      &
  fseafac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
  fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface flux buoyancy over density (m^2/s^3)
  u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
    !Surface friction velocity (m/s)
  t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Standard deviation of turbulent fluctuations of layer 1 temp; used in
    !initiating convection.
  q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Standard deviation of turbulent flux of layer 1 humidity; used in
    !initiating convection.
  rhostar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
    !Surface air density
  vshr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !Magnitude of surface-to-lowest atm level wind shear (m per s).
  resp_s_tot_soilt(dim_cs2,nsoilt),                                           &
    !Total soil respiration (kg C/m2/s).
  emis_soil(land_pts)
    !Emissivity of underlying soil

!Misc INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  bq_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !A buoyancy parameter (beta q tilde).
  bt_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                  &
    !A buoyancy parameter (beta T tilde).
  zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                    &
    !Height above surface of top of boundary layer (metres).
  photosynth_act_rad(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !Net downward shortwave radiation in band 1 (w/m2).
INTEGER, INTENT(IN) ::                                                        &
  curr_year, curr_day_number, curr_hour, curr_minute, curr_second

!Forcing INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                    &
    !Total water content
  tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                    &
    !Ice/liquid water temperature
  pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                 &
    !Surface pressure (Pascals).
  lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !Surface downward LW radiation (W/m2).

!Fluxes INTENT(IN)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  sw_surft(land_pts,nsurft),                                                  &
    !Surface net SW radiation on land tiles (W/m2).
  tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
    !GBM surface temperature (K).

!diagnostic array
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag

!Fluxes INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  fqw_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !Moisture flux between layers (kg per square metre per sec).
  ftl_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                 &
    !FTL(,K) contains net turbulent sensible heat flux into layer K from below;
    !so FTL(,1) is the surface sensible heat, H.(W/m2)
  ftl_surft(land_pts,nsurft),                                                 &
    !Surface FTL for land tiles
  fqw_surft(land_pts,nsurft),                                                 &
    !Surface FQW for land tiles
  fqw_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Surface FQW for sea-ice
  ftl_ice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),      &
    !Surface FTL for sea-ice
  fsmc(land_pts,npft),                                                        &
    !Moisture availability factor.
  emis_surft(land_pts,nsurft)
    !Emissivity for land tiles

!Misc INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  radnet_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Surface net radiation on sea-ice (W/m2)
  rhokm_1(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),       &
    !Exchange coefficients for momentum on P-grid
  rhokm_land(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),    &
  rhokm_ssi(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

!Out of explicit and into implicit only INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
!   Interpolation coefficient for diagnosis of 10 m winds
  alpha1(land_pts,nsurft),                                                    &
    !Mean gradient of saturated specific humidity with respect to temperature
    !between the bottom model layer and tile surfaces
  alpha1_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !ALPHA1 for sea-ice.
  alpha1_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
    !ALPHA1 for sea.
  ashtf_prime(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Coefficient to calculate surface heat flux into sea-ice.
  ashtf_prime_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
    !Coefficient to calculate surface heat flux into sea.
  ashtf_prime_surft(land_pts,nsurft),                                         &
    !Coefficient to calculate surface heat flux into land tiles.
  epot_surft(land_pts,nsurft),                                                &
    !Local EPOT for land tiles.
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
    !Grid-box surface exchange coefficients
  rhokh_surft(land_pts,nsurft),                                               &
    !Surface exchange coefficients for land tiles
  rhokh_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),   &
    !Surface exchange coefficients for sea-ice
  rhokh_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
    !Surface exchange coefficients for sea
  dtstar_surft(land_pts,nsurft),                                              &
    !Change in TSTAR over timestep for land tiles
  dtstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
    !Change is TSTAR over timestep for open sea
  dtstar_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),  &
    !Change is TSTAR over timestep for sea-ice
  z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    !Roughness length for heat and moisture over sea (m).
  z0h_surft(land_pts,nsurft),                                                 &
    !Tile roughness lengths for heat and moisture (m).
  z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
    !Roughness length for momentum over sea (m).
  z0m_surft(land_pts,nsurft),                                                 &
    !Tile roughness lengths for momentum.
  chr1p5m(land_pts,nsurft),                                                   &
    !Ratio of coefffs for calculation of 1.5m temp for land tiles.
  chr1p5m_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
    !CHR1P5M for sea and sea-ice (leads ignored).
  canhc_surft(land_pts,nsurft),                                               &
    !Areal heat capacity of canopy for land tiles (J/K/m2).
  wt_ext_surft(land_pts,sm_levels,nsurft),                                    &
    !Fraction of evapotranspiration which is extracted from each soil layer
    !by each tile.
  flake(land_pts,nsurft)
    !Lake fraction.

!Out of explicit and into extra only INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  hcons_soilt(land_pts,nsoilt)
    !Soil thermal conductivity including water and ice

!Out of explicit and into implicit and extra INTENT(OUT)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  tile_frac(land_pts,nsurft)
    !Tile fractions including snow cover in the ice tile.

#if defined(UM_JULES)
REAL(KIND=real_jlslsm) ::                                                     &
  clay_soilt(land_pts,nsoilt,dim_cslayer)
    !Clay content of underlying soil
REAL(KIND=real_jlslsm) ::                                                     &
  cs_pool_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1),                         &
  resp_s_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1),                          &
  resp_s_acc_soilt(land_pts_trif,nsoilt,dim_cslayer,dim_cs1)
#else
LOGICAL :: l_dust_diag
  !In standalone, this switch essentially does the same job as l_dust
#endif

INTEGER ::                                                                    &
   i,j,l,n         !Various counters
INTEGER :: errorstatus

!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
! Temp until the model switching is implemented for coupled mode
! Having a parameter until then should hopefully help the compiler eliminate
! dead code

!passed from control
real :: surf_down_sw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,4)

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SURF_COUPLE_EXPLICIT'

!-----------------------------------------------------------------------------
!End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Change 2d UM clay to 1d jules clay content for soil respiration.
! Soil tiling not currently in the UM, so broadcast ij value to all tiles.
! Multi-layer clay not currently in UM so set all layers to same value.
#if defined(UM_JULES)
IF ( soil_bgc_model == soil_model_rothc ) THEN
  DO l = 1, land_pts
    j = (land_index(l) - 1) / pdims%i_end + 1
    i = land_index(l) - (j-1) * pdims%i_end
    DO n = 1, dim_cslayer
      clay_soilt(l,1,n) = soil_clay_ij(i,j)
    END DO
  END DO
END IF

DO l = 1, land_pts
  DO n = 1, dim_cs1
    cs_pool_soilt(l,1,1,n) = cs_pool_gb_um(l,n)
    resp_s_soilt(l,1,1,n)  = 0.0
    !Can't copy from resp_s_gb_um as it's INTENT(OUT) so set to zero
  END DO
END DO
DO l = 1, land_pts_trif
  DO n = 1, dim_cs1
    resp_s_acc_soilt(l,1,1,n) = resp_s_acc_gb_um(l,n)
  END DO
END DO
#else
l_dust_diag = l_dust
#endif

SELECT CASE( lsm_id )
CASE ( jules )
  CALL sf_expl_l (                                                            &
    !IN date-related values
    curr_year, curr_day_number, curr_hour, curr_minute, curr_second,          &
    !IN values defining field dimensions and subset to be processed :
    land_pts, nice, nice_use,                                                 &
    !IN  parameters for iterative SISL scheme
    numcycles, cycleno,                                                       &
    !IN parameters required from boundary-layer scheme :
    bq_1,bt_1,z1_uv_ij,z1_uv_top,z1_tq_ij,z1_tq_top,qw,tl,                    &
    !IN soil/vegetation/land surface data :
    land_index,nsurft,sm_levels,canopy_surft,catch_surft,catch_snow_surft,    &
    hcon_soilt,ho2r2_orog_gb, flandg,                                         &
    snow_surft,sil_orog_land_gb,smvccl_soilt,smvcst_soilt,smvcwt_soilt,       &
    sthf_soilt, sthu_soilt,z0_surft, z0h_bare_surft, z0m_soil_gb,             &
    !IN sea/sea-ice data :
    ice_fract_ncat_sicat,k_sice_sicat,                                        &
    !IN input data from the wave model
    charnock_w,                                                               &
    !IN everything not covered so far :
    pstar,lw_down,sw_surft,zh,ddmfx,                                          &
    co2_mmr,co2_3d_ij,l_co2_interactive,l_phenol,                             &
    asteps_since_triffid,cs_pool_soilt,frac_surft,canht_pft,                  &
    photosynth_act_rad, lai_pft,                                              &
    lq_mix_bl,t_soil_soilt,tsurf_elev_surft,                                  &
    ti_sicat,ti_cat_sicat,tstar,tstar_sea_ij,                                 &
    tstar_sice_sicat,                                                         &
    tstar_surft,z_land_ij,albsoil_soilt,cosz_ij,                              &
    l_aero_classic,l_dust,l_dust_diag,clay_soilt,o3_gb,                       &
    !IN idealised and SCM things
    l_spec_z0, z0m_scm, z0h_scm,                                              &
    !IN variables for message passing
    u_1_p_ij, v_1_p_ij, u_0_p_ij, v_0_p_ij,                                   &
    !INOUT diagnostics
    sf_diag,                                                                  &
    !INOUT data :
    z0msea_ij,gs_gb,g_leaf_acc_pft,npp_acc_pft,resp_w_acc_pft,                &
    resp_s_acc_soilt,                                                         &
    !OUT Diagnostic not requiring STASH flags :
    recip_l_mo_sea,fqw_1,ftl_1,ftl_surft,                                     &
    radnet_sice,rhokm_1,rib,                                                  &
    !OUT variables for message passing
    flandfac, fseafac, rhokm_land, rhokm_ssi,                                 &
    cdr10m,                                                                   &
    !OUT data required for tracer mixing :
    rho_aresist_ij,aresist_ij,resist_b_ij,                                    &
    rho_aresist_surft,aresist_surft,resist_b_surft,                           &
    !OUT data required for mineral dust scheme
    r_b_dust_ij,cd_std_dust_ij,u_s_std_surft,                                 &
    !OUT data required elsewhere in UM system :
    fb_surf,u_s,t1_sd,q1_sd,                                                  &
    !OUT data required elsewhere in boundary layer or surface code
    alpha1,alpha1_sea,alpha1_sice,ashtf_prime,ashtf_prime_sea,                &
    ashtf_prime_surft,fqw_surft,                                              &
    epot_surft,fqw_ice,ftl_ice,fraca,rhostar,resfs,resft,                     &
    rhokh,rhokh_surft,rhokh_sice,rhokh_sea,                                   &
    dtstar_surft,dtstar_sea,dtstar_sice,                                      &
    h_blend_orog_ij,z0hssi,z0h_surft, z0mssi,z0m_surft,                       &
    z0m_eff_ij,chr1p5m,chr1p5m_sice,smc_soilt,hcons_soilt,vshr,               &
    vshr_land_ij,                                                             &
    vshr_ssi_ij, gpp_gb,npp_gb,resp_p_gb,g_leaf_pft,gpp_pft,npp_pft,          &
    resp_p_pft,resp_s_soilt,resp_s_tot_soilt,resp_l_pft,resp_r_pft,           &
    resp_w_pft,                                                               &
    n_leaf_pft,n_root_pft,n_stem_pft,lai_bal_pft,                             &
    gc_surft,canhc_surft,wt_ext_surft,flake,                                  &
    surft_index,surft_pts,tile_frac,fsmc,emis_surft,emis_soil                 &
    )

#if defined(UM_JULES)
  DO l = 1, land_pts
    DO n = 1, dim_cs1
      resp_s_gb_um(l,n) = resp_s_soilt(l,1,1,n)
    END DO
  END DO
  DO l = 1, land_pts_trif
    DO n = 1, dim_cs1
      resp_s_acc_gb_um(l,n) = resp_s_acc_soilt(l,1,1,n)
    END DO
  END DO
#endif

CASE ( cable )
#if defined(UM_JULES)
  errorstatus = 101
  CALL ereport('surf_couple_explicit', errorstatus,                           &
               'CABLE not yet implemented')
#else
  ! for testing LSM switch
  CALL jules_print('surf_couple_explicit', jules_message)
  
! initialise all INTENT(OUT) for now until CABLE is implementedCABLE_LSM: initialise else crashes model: havent traced yet
  dtstar_surft(:,:) = 0.0

  CALL sf_expl_l_cbl (                                                            &
    !IN date-related values
    curr_year, curr_day_number, curr_hour, curr_minute, curr_second,          &
    !IN values defining field dimensions and subset to be processed :
    land_pts, nice, nice_use,                                                 &
    !IN  parameters for iterative SISL scheme
    numcycles, cycleno,                                                       &
    !IN parameters required from boundary-layer scheme :
    bq_1,bt_1,z1_uv_ij,z1_uv_top,z1_tq_ij,z1_tq_top,qw,tl,                    &
    !IN soil/vegetation/land surface data :
    land_index,nsurft,sm_levels,canopy_surft,catch_surft,catch_snow_surft,    &
    hcon_soilt,ho2r2_orog_gb, flandg,                                         &
    snow_surft,sil_orog_land_gb,smvccl_soilt,smvcst_soilt,smvcwt_soilt,       &
    sthf_soilt, sthu_soilt,z0_surft, z0h_bare_surft, z0m_soil_gb,             &
    !IN sea/sea-ice data :
    ice_fract_ncat_sicat,k_sice_sicat,                                        &
    !IN input data from the wave model
    charnock_w,                                                               &
    !IN everything not covered so far :
    pstar,lw_down,sw_surft,zh,ddmfx,                                          &
    co2_mmr,co2_3d_ij,l_co2_interactive,l_phenol,                             &
    asteps_since_triffid,cs_pool_soilt,frac_surft,canht_pft,                  &
    photosynth_act_rad, lai_pft,                                              &
    lq_mix_bl,t_soil_soilt,tsurf_elev_surft,                                  &
    ti_sicat,ti_cat_sicat,tstar,tstar_sea_ij,                                 &
    tstar_sice_sicat,                                                         &
    tstar_surft,z_land_ij,albsoil_soilt,cosz_ij,                              &
    l_aero_classic,l_dust,l_dust_diag,clay_soilt,o3_gb,                       &
    !IN idealised and SCM things
    l_spec_z0, z0m_scm, z0h_scm,                                              &
    !IN variables for message passing
    u_1_p_ij, v_1_p_ij, u_0_p_ij, v_0_p_ij,                                   &
    !INOUT diagnostics
    sf_diag,                                                                  &
    !INOUT data :
    z0msea_ij,gs_gb,g_leaf_acc_pft,npp_acc_pft,resp_w_acc_pft,                &
    resp_s_acc_soilt,                                                         &
    !OUT Diagnostic not requiring STASH flags :
    recip_l_mo_sea,fqw_1,ftl_1,ftl_surft,                                     &
    radnet_sice,rhokm_1,rib,                                                  &
    !OUT variables for message passing
    flandfac, fseafac, rhokm_land, rhokm_ssi,                                 &
    cdr10m,                                                                   &
    !OUT data required for tracer mixing :
    rho_aresist_ij,aresist_ij,resist_b_ij,                                    &
    rho_aresist_surft,aresist_surft,resist_b_surft,                           &
    !OUT data required for mineral dust scheme
    r_b_dust_ij,cd_std_dust_ij,u_s_std_surft,                                 &
    !OUT data required elsewhere in UM system :
    fb_surf,u_s,t1_sd,q1_sd,                                                  &
    !OUT data required elsewhere in boundary layer or surface code
    alpha1,alpha1_sea,alpha1_sice,ashtf_prime,ashtf_prime_sea,                &
    ashtf_prime_surft,fqw_surft,                                              &
    epot_surft,fqw_ice,ftl_ice,fraca,rhostar,resfs,resft,                     &
    rhokh,rhokh_surft,rhokh_sice,rhokh_sea,                                   &
    dtstar_surft,dtstar_sea,dtstar_sice,                                      &
    h_blend_orog_ij,z0hssi,z0h_surft, z0mssi,z0m_surft,                       &
    z0m_eff_ij,chr1p5m,chr1p5m_sice,smc_soilt,hcons_soilt,vshr,               &
    vshr_land_ij,                                                             &
    vshr_ssi_ij, gpp_gb,npp_gb,resp_p_gb,g_leaf_pft,gpp_pft,npp_pft,          &
    resp_p_pft,resp_s_soilt,resp_s_tot_soilt,resp_l_pft,resp_r_pft,           &
    resp_w_pft,                                                               &
    n_leaf_pft,n_root_pft,n_stem_pft,lai_bal_pft,                             &
    gc_surft,canhc_surft,wt_ext_surft,flake,                                  &
    surft_index,surft_pts,tile_frac,fsmc,emis_surft,emis_soil,                &
    !CABLE_LSM:delete
    fland, surf_down_sw )
  

#endif

CASE DEFAULT
  errorstatus = 101
  CALL ereport('surf_couple_explicit', errorstatus,                           &
               'Unrecognised surface scheme')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE surf_couple_explicit

END MODULE surf_couple_explicit_mod

