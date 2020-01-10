MODULE surf_couple_extra_mod
! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

PUBLIC :: surf_couple_extra

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SURF_COUPLE_EXTRA_MOD'

CONTAINS

!===============================================================================
!Note that the commented intents may not be correct. The declarations section
!will be more correct

!Routine Layout **PLEASE READ**

!-Header (please follow ordering guidance)
!-Type correction for timestep
!-SELECT for LSM
!--JULES
!---Redimensioning of arrays (Note 1)
!---Compression to land points (Note 2)
!---Science calls (Note 3)
!---Expansion to full grid
!---UM Diagnostics calls
!--CABLE
!---Details to follow
!-END

!Note 1- This is for any variable that might have differnt numbers of
!        dimensions in standalone vs UM mode. These generally arise due to the
!        differing science payloads.

!Note 2- This is where all gridded variables (with _ij suffix) are compressed
!        onto land points. This allows efficient use of OpenMP

!Note 3- New additions here should ideally consist ONLY of calls to other
!        routines. Please encapsulate longer code within a subroutine.

!===============================================================================
SUBROUTINE surf_couple_extra(                                                 &
   ! Arguments used by JULES-standalone
   !Driving data and associated INTENT(IN)
   ls_rain_ij, con_rain_ij, ls_snow_ij, con_snow_ij, tl_1_ij, lw_down,        &
   qw_1_ij, u_1_ij, v_1_ij,                                                   &
   pstar_ij,                                                                  &
   !Fluxes INTENT(IN)
   ei_surft, surf_htf_surft, ecan_surft, ext_soilt, sw_surft,                 &
   !Misc INTENT(IN)
   a_step, smlt, tile_frac, hcons_soilt,                                      &
   !Fluxes INTENT(INOUT)
   melt_surft,                                                                &
   !Fluxes INTENT(OUT)
   snomlt_surf_htf, snowmelt_ij, snomlt_sub_htf, sub_surf_roff, surf_roff,    &
   tot_tfall, snowmelt_gb, rrun, rflow, snow_soil_htf,                        &
   !Arguments for the UM-----------------------------------------
   !IN
   land_pts, row_length, rows, river_row_length, river_rows, land_index,      &
   ls_graup_ij,                                                               &
   cca_2d, nsurft, surft_pts, surft_index,                                    &
   tstar_surft,                                                               &
   lice_pts, lice_index, soil_pts, soil_index,                                &
   stf_sub_surf_roff,                                                         &
   fexp_soilt, gamtot_soilt, ti_mean_soilt, ti_sig_soilt,                     &
   cs_ch4_soilt, flash_rate_ancil, pop_den_ancil,                             &
   a_fsat_soilt, c_fsat_soilt, a_fwet_soilt, c_fwet_soilt,                    &
   ntype, fqw_surft,                                                          &
   delta_lambda, delta_phi, xx_cos_theta_latitude,                            &
   aocpl_row_length, aocpl_p_rows, xpa, xua, xva, ypa, yua, yva,              &
   g_p_field, g_r_field, n_proc, global_row_length, global_rows,              &
   global_river_row_length, global_river_rows, flandg,                        &
   trivdir, trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land,       &
   frac_agr_gb, soil_clay_ij, resp_s_soilt, npp_gb,                           &
   z0m_soil_gb,                                                               &
   !INOUT
   a_steps_since_riv,  t_soil_soilt, tsurf_elev_surft,                        &
   rgrain_surft, snow_grnd_surft, snow_surft,                                 &
   smcl_soilt, sthf_soilt, sthu_soilt, canopy_surft, fsat_soilt, fwetl_soilt, &
   zw_soilt, sthzw_soilt,                                                     &
   snow_mass_ij,  ls_rainfrac_gb,                                             &
   substore, surfstore, flowin, bflowin,                                      &
   tot_surf_runoff, tot_sub_runoff, acc_lake_evap, twatstor,                  &
   asteps_since_triffid, g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft,    &
   resp_s_acc_gb_um,                                                          &
   resp_w_acc_pft, cs_pool_gb_um, frac_surft, lai_pft, canht_pft,             &
   catch_snow_surft, catch_surft, infil_surft,                                &
   inlandout_atm_gb,                                                          &
   !OUT
   dhf_surf_minus_soil,                                                       &
   canopy_gb,  smc_soilt,                                                     &
   z0_surft, z0h_bare_surft,                                                  &
   land_sea_mask                                                              &
   )

!Module imports

!Import interfaces to subroutines called
USE hydrol_mod,               ONLY: hydrol
USE hydrol_cbl_mod,           ONLY: hydrol_cbl
USE snow_mod,                 ONLY: snow
USE jules_rivers_mod,         ONLY: l_rivers, l_inland

! Code which isn't currently suitable for building into LFRic
#if !defined(LFRIC)
USE river_control_mod,        ONLY: river_control
USE inferno_mod,              ONLY: calc_soil_carbon_pools
USE inferno_io_mod,           ONLY: inferno_io
USE irrig_dmd_mod,            ONLY: irrig_dmd
USE veg_control_mod,          ONLY: veg_control

!ifdef required to manage the different science payloads available for
!coupled and standalone use
#if defined(UM_JULES)
USE diagnostics_riv_mod,      ONLY: diagnostics_riv
#else
USE crop_mod,                 ONLY: crop
USE fire_timestep_mod,        ONLY: fire_timestep
USE metstats_timestep_mod,    ONLY: metstats_timestep
USE gridbox_mean_mod,         ONLY: soiltiles_to_gbm, surftiles_to_gbm
USE fao_evapotranspiration,   ONLY: fao_ref_evapotranspiration
USE rivers_route_mod,         ONLY: adjust_routestore
USE zenith_mod,               ONLY: photoperiod
! End of standalone Jules code
#endif

! End of code excluded from LFRic builds
#else
! Code only used for LFRic
USE missing_data_mod, ONLY: rmdi
#endif

!Variables- The rules here are:
!Each module is only USED once, using an ifdef to control the variables if req.
!Whole modules are grouped into a single ifdef/else
!Alphabetical order

!Modules common to JULES and UM

USE ancil_info,               ONLY:                                           &
#if !defined(UM_JULES)
                                    dim_cs1, land_pts_trif, npft_trif,        &
#endif
                                    dim_cslayer, nsoilt

USE atm_fields_bounds_mod,    ONLY: tdims, tdims_s, pdims, pdims_s

USE crop_vars_mod,            ONLY:                                           &
#if !defined(UM_JULES)
  phot, dphotdt, rootc_cpft, harvc_cpft, reservec_cpft, croplai_cpft,         &
  cropcanht_cpft, nday_crop,                                                  &
#endif
  plant_n_gb, frac_irr_soilt, sthu_irr_soilt, dvi_cpft

USE fire_vars,                ONLY: pop_den, flash_rate

USE jules_hydrology_mod,      ONLY: l_hydrology, l_pdm, l_top, l_var_rainfrac

USE jules_soil_biogeochem_mod, ONLY:                                          &
  soil_model_1pool, soil_model_ecosse, soil_model_rothc, soil_bgc_model

USE jules_soil_mod,           ONLY: sm_levels, l_soil_sat_down, confrac

USE jules_surface_types_mod,  ONLY: npft, ncpft, nnpft, soil

USE jules_vegetation_mod,     ONLY:                                           &
  l_crop, l_triffid, l_trif_eq, l_phenol, phenol_period, triffid_period,      &
  l_irrig_dmd, irr_crop, l_irrig_limit, l_inferno, ignition_method,           &
  l_fao_ref_evapotranspiration, frac_min

USE prognostics,              ONLY:                                           &
#if defined(UM_JULES)
   wood_prod_fast_gb, wood_prod_med_gb, wood_prod_slow_gb,                    &
#endif
  rgrainl_surft, rho_snow_grnd_surft, sice_surft, sliq_surft,                 &
  tsnow_surft, ds_surft, snowdepth_surft, nsnow_surft, cs_pool_soilt,         &
  rho_snow_surft

USE p_s_parms,                ONLY:                                           &
  bexp_soilt, sathh_soilt, hcap_soilt, hcon_soilt, satcon_soilt,              &
  smvccl_soilt, smvcwt_soilt, smvcst_soilt, clay_soilt

USE sf_diags_mod,             ONLY: sf_diag

USE theta_field_sizes,        ONLY: t_i_length, t_j_length

USE top_pdm,                  ONLY:                                           &
  qbase_soilt, qbase_zw_soilt, fch4_wetl_soilt,                               &
  fch4_wetl_cs_soilt, fch4_wetl_npp_soilt, fch4_wetl_resps_soilt

USE trif_vars_mod, ONLY:                                                      &
  frac_past_gb, fao_et0

USE trifctl,                  ONLY:                                           &
  c_veg_pft, cv_gb, lit_c_pft, lit_c_mn_gb, g_leaf_day_pft, g_leaf_phen_pft,  &
  lai_phen_pft, g_leaf_dr_out_pft, npp_dr_out_pft, resp_w_dr_out_pft, npp_pft

!Modules specific to the UM
#if defined(UM_JULES)
USE atm_fields_real_mod,      ONLY: disturb_veg_prev

USE atm_step_local,           ONLY:                                           &
  land_pts_trif, npft_trif, dim_cs1, dim_cs2,                                 &
  STASHwork19, STASHwork8, STASHwork26

USE model_domain_mod,         ONLY: model_type, mt_single_column

USE stash_array_mod,          ONLY: sf

USE timestep_mod,             ONLY: timestep

USE um_parcore,               ONLY: mype

#else
!Modules specific to JULES
USE conversions_mod,          ONLY: isec_per_day, rsec_per_day

USE fire_mod,                 ONLY: fire_prog, fire_diag, l_fire

USE fluxes,                   ONLY: surf_ht_flux_ij

USE forcing,                  ONLY: sw_down_ij, lw_down_ij

USE metstats_mod,             ONLY: metstats_prog, l_metstats

USE model_time_mod,           ONLY: timestep_len, current_time
#endif


!Technical modules
USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook
USE ereport_mod,              ONLY: ereport
USE jules_print_mgr,          ONLY: jules_message, jules_print
USE lsm_switch_mod,           ONLY: lsm_id, jules, cable

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Private subroutine for accessing the JULES science routines that are called
!   after the implicit code
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!Subroutine Arguments, ordered by intent and type

INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
  row_length,                                                                 &
  rows,                                                                       &
  river_row_length,                                                           &
  river_rows,                                                                 &
  land_index(land_pts),                                                       &
  nsurft,                                                                     &
  a_step,                                                                     &
  surft_pts(nsurft),                                                          &
  surft_index(land_pts,nsurft),                                               &
  lice_pts,                                                                   &
  lice_index(land_pts),                                                       &
  soil_pts,                                                                   &
  soil_index(land_pts),                                                       &
  ntype,                                                                      &
  aocpl_row_length,                                                           &
  aocpl_p_rows,                                                               &
  g_p_field,                                                                  &
  g_r_field,                                                                  &
  n_proc,                                                                     &
  global_row_length,                                                          &
  global_rows,                                                                &
  global_river_row_length,                                                    &
  global_river_rows

LOGICAL, INTENT(IN) ::                                                        &
  smlt,                                                                       &
  stf_sub_surf_roff,                                                          &
  land_sea_mask(row_length, rows)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  !Driving data
  ls_rain_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
  con_rain_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
  ls_snow_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
  con_snow_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
  pstar_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
  ls_graup_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
  tl_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
  lw_down(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
  qw_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
  u_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
  v_1_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &

  !Fluxes
  ei_surft(land_pts,nsurft),                                                  &
       !Sublimation of snow (kg/m2/s)
  surf_htf_surft(land_pts,nsurft),                                            &
       !Surface heat flux (W/m2)
  ecan_surft(land_pts,nsurft),                                                &
       !Canopy evaporation from land tiles (kg/m2/s).
  ext_soilt(land_pts,nsoilt,sm_levels),                                       &
       !Extraction of water from each soil layer (kg/m2/s).
  tile_frac(land_pts,nsurft),                                                 &
  sw_surft(land_pts,nsurft),                                                  &
         !Surface net SW radiation on land tiles (W/m2)
  cca_2d(row_length,rows),                                                    &
  tstar_surft(land_pts,nsurft),                                               &
  fexp_soilt(land_pts,nsoilt),                                                &
  gamtot_soilt(land_pts,nsoilt),                                              &
  ti_mean_soilt(land_pts,nsoilt),                                             &
  ti_sig_soilt(land_pts,nsoilt),                                              &
  npp_gb(land_pts),                                                           &
  a_fsat_soilt(land_pts,nsoilt),                                              &
  c_fsat_soilt(land_pts,nsoilt),                                              &
  a_fwet_soilt(land_pts,nsoilt),                                              &
  c_fwet_soilt(land_pts,nsoilt),                                              &
  fqw_surft(land_pts,nsurft),                                                 &
  frac_agr_gb(land_pts),                                                      &
  soil_clay_ij(row_length,rows),                                              &
  z0m_soil_gb(land_pts),                                                      &
  flash_rate_ancil(row_length,rows),                                          &
  pop_den_ancil(row_length,rows),                                             &

  !River routing
  delta_lambda,                                                               &
  delta_phi,                                                                  &
  xx_cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,                        &
                        tdims_s%j_start:tdims_s%j_end),                       &
  xpa(aocpl_row_length+1),                                                    &
  xua(0:aocpl_row_length),                                                    &
  xva(aocpl_row_length+1),                                                    &
  ypa(aocpl_p_rows),                                                          &
  yua(aocpl_p_rows),                                                          &
  yva(0:aocpl_p_rows),                                                        &
  flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),        &
  trivdir(river_row_length, river_rows),                                      &
  trivseq(river_row_length, river_rows),                                      &
  r_area(row_length, rows),                                                   &
  slope(row_length, rows),                                                    &
  flowobs1(row_length, rows),                                                 &
  r_inext(row_length, rows),                                                  &
  r_jnext(row_length, rows),                                                  &
  r_land(row_length, rows)

INTEGER, INTENT(INOUT)  ::                                                    &
  a_steps_since_riv,                                                          &
  asteps_since_triffid

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  melt_surft(land_pts,nsurft),                                                &
        !Surface or canopy snowmelt rate (kg/m2/s)
        !On output, this is the total melt rate for the tile
        !(i.e. sum of  melt on canopy and ground).
  snomlt_sub_htf(land_pts),                                                   &
        !Sub-canopy snowmelt heat flux (W/m2)
  hcons_soilt(land_pts,nsoilt),                                               &
  dhf_surf_minus_soil(land_pts),                                              &
        ! Heat flux difference across the FLake snowpack (W/m2)
  resp_s_soilt(land_pts,nsoilt,1,dim_cs1),                                    &
  ls_rainfrac_gb(land_pts),                                                   &
  catch_snow_surft(land_pts,nsurft),                                          &
  t_soil_soilt(land_pts,nsoilt,sm_levels),                                    &
  tsurf_elev_surft(land_pts,nsurft),                                          &
  rgrain_surft(land_pts,nsurft),                                              &
  snow_grnd_surft(land_pts,nsurft),                                           &
  snow_surft(land_pts,nsurft),                                                &
  smcl_soilt(land_pts,nsoilt,sm_levels),                                      &
  sthf_soilt(land_pts,nsoilt,sm_levels),                                      &
  infil_surft(land_pts,nsurft),                                               &
  catch_surft(land_pts,nsurft),                                               &
  sthu_soilt(land_pts,nsoilt,sm_levels),                                      &
  canopy_surft(land_pts,nsurft),                                              &
  fsat_soilt(land_pts,nsoilt),                                                &
  fwetl_soilt(land_pts,nsoilt),                                               &
  zw_soilt(land_pts,nsoilt),                                                  &
  sthzw_soilt(land_pts,nsoilt),                                               &
  substore(row_length, rows),                                                 &
  surfstore(row_length, rows),                                                &
  flowin(row_length, rows),                                                   &
  bflowin(row_length, rows),                                                  &
  tot_surf_runoff(land_pts),                                                  &
  tot_sub_runoff(land_pts),                                                   &
  acc_lake_evap(row_length,rows),                                             &
  twatstor(river_row_length, river_rows),                                     &
  g_leaf_acc_pft(land_pts,npft),                                              &
  g_leaf_phen_acc_pft(land_pts,npft),                                         &
  npp_acc_pft(land_pts_trif,npft_trif),                                       &
  resp_s_acc_gb_um(land_pts_trif,dim_cs1),                                    &
  resp_w_acc_pft(land_pts_trif,npft_trif),                                    &
  cs_ch4_soilt(land_pts,nsoilt),                                              &
         ! soil carbon used in wetland CH4 emissions model if TRIFFID
         ! is switched off
  cs_pool_gb_um(land_pts,dim_cs1),                                            &
  frac_surft(land_pts,ntype),                                                 &
  lai_pft(land_pts,npft),                                                     &
  canht_pft(land_pts,npft),                                                   &
  inlandout_atm_gb(land_pts)

REAL(KIND=real_jlslsm), INTENT(OUT)::                                         &
  sub_surf_roff(land_pts),                                                    &
        !Sub-surface runoff (kg/m2/s).
  surf_roff(land_pts),                                                        &
        !Surface runoff (kg/m2/s).
  tot_tfall(land_pts),                                                        &
        !Total throughfall (kg/m2/s).
  snomlt_surf_htf(row_length,rows),                                           &
        !Gridbox snowmelt heat flux (W/m2)
  snow_soil_htf(land_pts,nsurft),                                             &
        !Tiled snow->soil heat flux (W/m2)
  snowmelt_gb(land_pts),                                                      &
  snowmelt_ij(row_length,rows),                                               &
  snow_mass_ij(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
  rrun(land_pts),                                                             &
         !Surface runoff after river routing (kg/m2/s)
  rflow(land_pts),                                                            &
         !River runoff (kg/m2/s)
  canopy_gb(land_pts),                                                        &
  smc_soilt(land_pts,nsoilt),                                                 &
  z0_surft(land_pts,nsurft),                                                  &
  z0h_bare_surft(land_pts,nsurft)

!==============================================================================
!Local variables

LOGICAL ::                                                                    &
  trip_call

! Need INVERT_OCEAN (hence N->S configuration of ATMOS grid) for river routing
LOGICAL, PARAMETER :: invert_ocean = .FALSE.

INTEGER ::                                                                    &
   i,j,l,n,m,                                                                 &
        !Various counters
  crop_call
        !indicates whether crop model is to be called

INTEGER, PARAMETER :: crop_period = 1
        ! Crop code hard wired to run daily : crop_period = 1

REAL(KIND=real_jlslsm) :: timestep_real ! Model timestep (s) in REAL

REAL(KIND=real_jlslsm) ::                                                     &
  !Gridbox versions of forcing data (ls_rainfrac_gb is an argument)
  con_rainfrac_gb(land_pts),                                                  &
  ls_rain_gb(land_pts),                                                       &
  con_rain_gb(land_pts),                                                      &
  ls_snow_gb(land_pts),                                                       &
  ls_graup_gb(land_pts),                                                      &
  con_snow_gb(land_pts),                                                      &
  pstar_gb(land_pts),                                                         &
  tl_1_gb(land_pts),                                                          &
  qw_1_gb(land_pts),                                                          &
  u_1_gb(land_pts),                                                           &
  v_1_gb(land_pts),                                                           &

  !Passed between different science schemes (various)
  qbase_l_soilt(land_pts,nsoilt,sm_levels+1),                                 &
        ! Base flow from each soil layer (kg m-2 s-1).
  w_flux_soilt(land_pts,nsoilt,0:sm_levels),                                  &
        ! Fluxes of water between layers (kg m-2 s-1).
  surf_ht_flux_ld(land_pts),                                                  &

  !Passed between fire & INFERNO routines only
  smc_gb(land_pts),                                                           &
         !To allow GBM soil moisture to be passed down to fire
  c_soil_dpm_gb(land_pts),                                                    &
         ! Gridbox soil C in the Decomposable Plant Material pool (kg m-2).
  c_soil_rpm_gb(land_pts),                                                    &
         ! Gridbox soil C in the Resistant Plant Material pool (kg m-2).

  !Passed between evapotranspiration routines only
  trad(land_pts),                                                             &
         ! gridbox effective radiative temperature (assuming emissivity=1)

  !Passed from snow or hydrol to diagnostics_hyd only
  snow_mass_gb(land_pts),                                                     &
  dun_roff_soilt(land_pts,nsoilt),                                            &
  drain_soilt(land_pts,nsoilt),                                               &

  !Passed between veg routines and diagnostics_veg only
  resp_s_dr_out_gb_um(land_pts,dim_cs1+1),                                    &

  !Passed between river routing and diagnostics_riv only
  riverout(row_length, rows),                                                 &
  riverout_rgrid(river_row_length, river_rows),                               &
  box_outflow(river_row_length, river_rows),                                  &
  box_inflow(river_row_length, river_rows),                                   &
  inlandout_riv(river_row_length,river_rows)

INTEGER                      :: errorstatus
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'SURF_COUPLE_EXTRA'

!Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!===============================================================================
!End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Care needed as the timestep is used in arithmetic to control the
!timing of when various science schemes are called.
#if defined(UM_JULES)
timestep_real =  timestep
#else
timestep_real = REAL(timestep_len)
#endif

!CABLE_LSM: implement switching based on lsm_id
SELECT CASE( lsm_id )

  !=============================================================================
CASE ( jules )

  !=============================================================================
  ! Redimensioning of arrays

#if defined(UM_JULES)
  !Dimensionality of variables differ betweeen UM and standalone, so copy across

  DO l = 1,land_pts
    DO n = 1,dim_cs1
      cs_pool_soilt(l,1,1,n) = cs_pool_gb_um(l,n)
    END DO
  END DO

  ! Change 2d to 1d soil clay content for soil respiration- dimensionalities
  ! differ between UM and standalone.
  ! Soil tiling not currently in the UM, so broadcast ij value to all tiles.
  ! Multi-layer clay not currently in UM so set all layers to same value.
  IF ( soil_bgc_model == soil_model_rothc ) THEN
    m = 1
    DO l = 1, land_pts
      j = (land_index(l) - 1) / row_length + 1
      i = land_index(l) - (j-1) * row_length
      DO n = 1, dim_cslayer
        clay_soilt(l,m,n) = soil_clay_ij(i,j)
      END DO
    END DO
  END IF

  !Compress pop_den and flash rate fields to land points- in the UM the ancil
  !is 2D
  IF (l_inferno) THEN
    DO l = 1, land_pts
      j = (land_index(l) - 1) / row_length + 1
      i = land_index(l) - (j-1) * row_length
      pop_den(l)    = pop_den_ancil(i,j)
      flash_rate(l) = flash_rate_ancil(i,j)
    END DO
  END IF
#endif

  !Encompassing IF statement for land_pts > 0. We exit/re-enter this IF about
  !half way down to allow for river routing

  IF (l_hydrology .AND. land_pts > 0 ) THEN

    !===========================================================================
    ! Compression to land points

    DO l = 1, land_pts
      j = (land_index(l) - 1) / row_length + 1
      i = land_index(l) - (j-1) * row_length
      ls_rain_gb(l)    = ls_rain_ij(i,j)
      con_rain_gb(l)   = con_rain_ij(i,j)
      con_snow_gb(l)   = con_snow_ij(i,j)
      ls_snow_gb(l)    = ls_snow_ij(i,j)
      ls_graup_gb(l)   = ls_graup_ij(i,j)
      pstar_gb(l)      = pstar_ij(i,j)
      tl_1_gb(l)       = tl_1_ij(i,j)
      qw_1_gb(l)       = qw_1_ij(i,j)
      u_1_gb(l)        = u_1_ij(i,j)
      v_1_gb(l)        = v_1_ij(i,j)
    END DO

    !Pass jules the modelled rain fractions
    !In standalone mode, both rainfracs are passed in as zeroed arrays
    IF (l_var_rainfrac) THEN
      DO l = 1, land_pts
        j = (land_index(l) - 1) / row_length + 1
        i = land_index(l) - (j-1) * row_length

        con_rainfrac_gb(l) = MIN(cca_2d(i,j),0.5) !As in diagnostics_conv

        !provide some safety checking for convective rain with no CCA
        IF (con_rain_gb(l) > 0.0 ) THEN
          IF (con_rainfrac_gb(l) == 0.0) THEN
            con_rainfrac_gb(l) = confrac
            !and for very small CCA amounts
          ELSE IF (con_rainfrac_gb(l) < 0.01) THEN
            con_rainfrac_gb(l) = 0.01
          END IF
        END IF

        !provide some safety checking for ls rain with no rainfrac
        IF (ls_rain_gb(l) > 0.0) THEN
          IF (ls_rainfrac_gb(l) == 0.0) THEN
            ls_rainfrac_gb(l) = 0.5
            !and for very small rainfrac amounts
          ELSE IF (ls_rainfrac_gb(l) < 0.01) THEN
            ls_rainfrac_gb(l) = 0.01
          END IF
        END IF
      END DO

    ELSE !use original default values
      DO l = 1, land_pts
        con_rainfrac_gb(l) = confrac
        ls_rainfrac_gb(l)  = 1.0
      END DO
    END IF !l_var_rainfrac

    !===========================================================================
    ! Science calls

    !Snow (standalone and UM)
    CALL snow ( land_pts,timestep_real,smlt,nsurft,surft_pts,                 &
                surft_index,catch_snow_surft,con_snow_gb,con_rain_gb,         &
                tile_frac,ls_snow_gb,ls_graup_gb,ls_rain_gb,                  &
                ei_surft,hcap_soilt(:,:,1),hcons_soilt,melt_surft,            &
                smcl_soilt(:,:,1),sthf_soilt(:,:,1),surf_htf_surft,           &
                t_soil_soilt(:,:,1),tsurf_elev_surft,                         &
                tstar_surft,smvcst_soilt(:,:,1),rgrain_surft,rgrainl_surft,   &
                rho_snow_grnd_surft,                                          &
                sice_surft,sliq_surft,snow_grnd_surft,snow_surft,             &
                snowdepth_surft, tsnow_surft,nsnow_surft,ds_surft,            &
                snomlt_surf_htf,snow_mass_gb,rho_snow_surft,snomlt_sub_htf,   &
                snowmelt_gb,snow_soil_htf,surf_ht_flux_ld,sf_diag,            &
                dhf_surf_minus_soil )

    !Hydrology (standalone and UM)
    CALL hydrol (                                                             &
      lice_pts,lice_index,soil_pts,soil_index, nsnow_surft,                   &
      land_pts,sm_levels,bexp_soilt,catch_surft,con_rain_gb,                  &
      ecan_surft,ext_soilt,hcap_soilt,hcon_soilt,ls_rain_gb,                  &
      con_rainfrac_gb, ls_rainfrac_gb,                                        &
      satcon_soilt,sathh_soilt,snowdepth_surft, snow_soil_htf,                &
      surf_ht_flux_ld,timestep_real,                                          &
      smvcst_soilt,smvcwt_soilt,canopy_surft,                                 &
      stf_sub_surf_roff,smcl_soilt,sthf_soilt,sthu_soilt,                     &
      t_soil_soilt,tsurf_elev_surft,canopy_gb,smc_soilt,snowmelt_gb,          &
      sub_surf_roff,surf_roff,tot_tfall,                                      &
      ! add new inland basin variable
      inlandout_atm_gb,l_inland,                                              &
      ! Additional variables for MOSES II
      nsurft,surft_pts,surft_index,                                           &
      infil_surft, melt_surft,tile_frac,                                      &
      ! Additional variables required for large-scale hydrology:
      l_top,l_pdm,fexp_soilt,ti_mean_soilt,cs_ch4_soilt,cs_pool_soilt,        &
      dun_roff_soilt,drain_soilt,fsat_soilt,fwetl_soilt,qbase_soilt,          &
      qbase_l_soilt, qbase_zw_soilt, w_flux_soilt,                            &
      zw_soilt,sthzw_soilt,a_fsat_soilt,c_fsat_soilt,a_fwet_soilt,            &
      c_fwet_soilt,                                                           &
      resp_s_soilt,npp_gb,fch4_wetl_soilt,                                    &
      fch4_wetl_cs_soilt,fch4_wetl_npp_soilt,fch4_wetl_resps_soilt,           &
      dim_cs1,l_soil_sat_down,l_triffid,asteps_since_triffid)

  END IF ! ( l_hydrology .AND. land_pts /= 0 )

! Code not yet ported to LFRic
#if !defined(LFRIC)

  !Here we need to exit the land_pts IF to allow river routing to be called
  !on all MPI ranks. This is because it is possible that a rank with
  !land_pts = 0 still has a river routing point.

  ! River Routing (standalone and UM, but differing options)
  IF ( l_rivers ) THEN
#if defined(UM_JULES)
    CALL river_control(                                                       &
      !LOGICAL, INTENT(IN)
      invert_ocean,                                                           &
      !INTEGER, INTENT(IN)
      n_proc, land_pts, row_length, rows, river_row_length, river_rows,       &
      land_index, ntype, aocpl_row_length, aocpl_p_rows, g_p_field,           &
      g_r_field, mype, global_row_length, global_rows,                        &
      global_river_row_length,                                                &
      global_river_rows, nsurft,                                              &
      !REAL, INTENT(IN)
      fqw_surft, delta_lambda, delta_phi, xx_cos_theta_latitude,              &
      xpa, xua, xva, ypa, yua, yva, flandg, trivdir,                          &
      trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land,             &
      smvcst_soilt, smvcwt_soilt,                                             &
      surf_roff, sub_surf_roff, frac_surft,                                   &
      !INTEGER, INTENT(INOUT)
      a_steps_since_riv,                                                      &
      !REAL, INTENT(INOUT)
      substore, surfstore, flowin, bflowin,                                   &
      tot_surf_runoff, tot_sub_runoff, acc_lake_evap, twatstor,               &
      smcl_soilt, sthu_soilt,                                                 &
      !LOGICAL, INTENT(OUT)
      trip_call,                                                              &
      !REAL, INTENT(OUT)
      inlandout_atm_gb, inlandout_riv, riverout,                              &
      riverout_rgrid, box_outflow, box_inflow, rflow, rrun                    &
      )
#else
    CALL river_control( land_pts,sub_surf_roff,surf_roff,rflow,rrun)
#endif
  END IF ! l_rivers (ATMOS)

  !Irrigation (standalone and UM for some options)
  IF ( l_irrig_dmd ) THEN
    CALL irrig_dmd(land_pts, sm_levels, frac_irr_soilt,                       &
                   a_step, plant_n_gb,                                        &
                   sthf_soilt, smvccl_soilt, smvcst_soilt, smvcwt_soilt,      &
                   sthzw_soilt, sthu_irr_soilt, sthu_soilt,                   &
                   smcl_soilt, irr_crop, dvi_cpft)

#if !defined(UM_JULES)
    ! Apply irrigation to planting dates as specified, apply limitations to
    ! irrigation based on the water in rivers/groundwater
    ! irr_crop == 1 is not advisable for UM use. It requires detailed
    ! information about planting dates to be calculated and there is a large
    ! amount of technical debt which needs to be sorted before going into the UM
    IF ( irr_crop == 1 ) THEN
      CALL calc_crop_date(land_index, land_pts, t_i_length, t_j_length,       &
                          nsurft,                                             &
                          frac_surft, sw_surft, tstar_surft, lw_down, tl_1_ij,&
                          con_rain_ij, ls_rain_ij, con_snow_ij, ls_snow_ij,   &
                          plant_n_gb, nday_crop)
    END IF

    ! Limitation code requires l_irrig_dmd = TRUE, l_top = TRUE, l_rivers = TRUE
    ! and i_river_UM = rivers_trip (3) . The technical parts of this code are
    ! designed only to run with standalone TRIP routing code. UM TRIP/RFM code is
    ! deprecated.
    IF ( l_irrig_limit ) THEN
      CALL adjust_routestore()
    END IF
#endif

  ELSE
    ! Set sthu_irr_soilt to 0.0 in case it is still reported
    ! It would be better to not allocate this array when it is not being used.
    ! This led to a memory leak in the D1 array (UM).
    sthu_irr_soilt(:,:,:) = 0.0
  END IF ! l_irrig_dmd

  !Restart the IF for land_pts > 0 now that river routing is done with
  IF (land_pts > 0) THEN

    !Crops (standalone only)
#if !defined(UM_JULES)
    IF ( l_crop ) THEN
      crop_call = MOD ( REAL(a_step),                                         &
                        REAL(crop_period) * rsec_per_day / timestep_real )

      DO n = 1,ncpft
        DO l = 1,land_pts
          npp_acc_pft(l,nnpft + n) = npp_acc_pft(l,nnpft + n)                 &
                                + (npp_pft(l,nnpft + n) * timestep_real)
        END DO
      END DO

      CALL photoperiod(t_i_length * t_j_length, phot, dphotdt)

      CALL crop(t_i_length * t_j_length, land_pts, land_index, a_step,        &
                crop_call, sm_levels, frac_surft, phot, dphotdt,              &
                sf_diag%t1p5m_surft, t_soil_soilt, sthu_soilt, smvccl_soilt,  &
                smvcst_soilt, npp_acc_pft,                                    &
                canht_pft, lai_pft, dvi_cpft, rootc_cpft, harvc_cpft,         &
                reservec_cpft, croplai_cpft, cropcanht_cpft,                  &
                catch_surft, z0_surft)
    END IF  ! l_crop

    !Metstats (standalone only)
        !Beware- does not currently account for graupel
    IF ( l_metstats ) THEN
      CALL metstats_timestep(tl_1_gb, qw_1_gb, u_1_gb, v_1_gb, ls_rain_gb,    &
                            con_rain_gb, ls_snow_gb, con_snow_gb, pstar_gb,   &
                            metstats_prog,                                    &
                            !Vars that should be USED but can't due to
                            !UM/standalone differences
                            current_time%time, timestep_real,land_pts)
    END IF

    !Fire (standalone only)
    IF ( l_fire ) THEN
      !Calculate the gridbox mean soil moisture
      smc_gb = soiltiles_to_gbm(smc_soilt)
      CALL fire_timestep(metstats_prog, smc_gb, fire_prog, fire_diag,         &
                        !Vars that should be USED but can't due to
                        !UM/standalone differences
                        current_time%time, current_time%month, timestep_real, &
                        land_pts)
    END IF
#endif

    !INFERNO (standalone and UM)
    IF ( l_inferno ) THEN
      CALL calc_soil_carbon_pools(land_pts, soil_pts, soil_index, dim_cs1,    &
                                  cs_pool_soilt,                              &
                                  c_soil_dpm_gb, c_soil_rpm_gb)

      CALL inferno_io( sf_diag%t1p5m_surft, sf_diag%q1p5m_surft, pstar_gb,    &
                      sthu_soilt, sm_levels,                                  &
                      frac_surft, c_soil_dpm_gb, c_soil_rpm_gb, canht_pft,    &
                      ls_rain_gb, con_rain_gb,                                &
                      pop_den, flash_rate,                                    &
                      land_pts, ignition_method,                              &
                      nsurft, asteps_since_triffid)
    END IF  !  l_inferno

    !Vegetation (standalone and UM- differing functionality)
    IF (l_phenol .OR. l_triffid) THEN
      CALL veg_control(                                                       &
        land_pts, nsurft, dim_cs1,                                            &
        a_step, asteps_since_triffid,                                         &
        land_pts_trif, npft_trif,                                             &
        phenol_period, triffid_period,                                        &
        l_phenol, l_triffid, l_trif_eq,                                       &
        timestep_real, frac_agr_gb, frac_past_gb, satcon_soilt,               &
        g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft,                     &
        resp_s_acc_gb_um, resp_w_acc_pft,                                     &
        cs_pool_gb_um, frac_surft, lai_pft, clay_soilt, z0m_soil_gb,          &
        canht_pft,                                                            &
        catch_snow_surft, catch_surft, infil_surft, z0_surft, z0h_bare_surft, &
        c_veg_pft, cv_gb, lit_c_pft, lit_c_mn_gb, g_leaf_day_pft,             &
        g_leaf_phen_pft,                                                      &
        lai_phen_pft, g_leaf_dr_out_pft, npp_dr_out_pft, resp_w_dr_out_pft,   &
        resp_s_dr_out_gb_um, qbase_l_soilt, sthf_soilt, sthu_soilt,           &
        w_flux_soilt, t_soil_soilt,cs_pool_soilt)
    END IF

#if !defined(UM_JULES)
    !Reference evapotranspiration (standalone only)
    IF (l_fao_ref_evapotranspiration) THEN
      trad = ( surftiles_to_gbm(tstar_surft**4) )**0.25
      CALL fao_ref_evapotranspiration(soil_pts, soil_index,                   &
        land_pts, land_index, sf_diag%t1p5m,                                  &
        sw_down_ij, lw_down_ij, surf_ht_flux_ij, sf_diag%u10m,                &
        sf_diag%v10m, sf_diag%q1p5m, pstar_ij, trad, fao_et0)
    END IF
! End of standalone Jules code
#endif

! End of code excluded from LFRic builds
#else
  ! These variables are intent out and so must be initialised for LFRic to
  ! compile. Set to rmdi in the hope this will show errors if they are used!
  rrun = rmdi
  rflow = rmdi
  z0_surft = rmdi
  z0h_bare_surft = rmdi

  IF (land_pts > 0) THEN

! End of bespoke LFRic code
#endif

    !End of science calls
    !===========================================================================
    !Expansion to full grid

    DO l = 1, land_pts
      j=(land_index(l) - 1) / row_length + 1
      i = land_index(l) - (j-1) * row_length
      snow_mass_ij(i,j) = snow_mass_gb(l)
      snowmelt_ij(i,j)  = snowmelt_gb(l)
    END DO

  END IF ! land_pts > 0

  !=============================================================================
  !UM Diagnostics calls

    !For the UM, soil tiling has not been implemented, ie nsoilt = 1, so we can
    !hard-code _soilt variables with index 1 using m = 1
#if defined(UM_JULES) && !defined(LFRIC)
  IF (model_type /= mt_single_column) THEN
    m = 1
    IF (l_hydrology .AND. sf(0,8) ) THEN
      ! DEPENDS ON: diagnostics_hyd
      CALL diagnostics_hyd(                                                   &
        row_length, rows,                                                     &
        land_pts, sm_levels,                                                  &
        land_index,inlandout_atm_gb,                                          &
        smc_soilt(:,m), surf_roff, sub_surf_roff,                             &
        snow_mass_gb, snowmelt_gb,                                            &
        canopy_gb,t_soil_soilt(:,m,:),                                        &
        tsurf_elev_surft,snow_soil_htf,                                       &
        smcl_soilt(:,m,:),                                                    &
        nsurft, snomlt_surf_htf, sthu_soilt(:,m,:), sthf_soilt(:,m,:),        &
        tot_tfall, snow_surft, melt_surft,                                    &
        rgrain_surft, land_sea_mask,                                          &
        dun_roff_soilt(:,m), drain_soilt(:,m), qbase_soilt(:,m),              &
        qbase_zw_soilt(:,m), fch4_wetl_soilt(:,m),fch4_wetl_cs_soilt(:,m),    &
        fch4_wetl_npp_soilt(:,m),fch4_wetl_resps_soilt(:,m),                  &
        fexp_soilt(:,m),gamtot_soilt(:,m),ti_mean_soilt(:,m),                 &
        ti_sig_soilt(:,m),                                                    &
        fsat_soilt(:,m),fwetl_soilt(:,m),zw_soilt(:,m),sthzw_soilt(:,m),      &
        timestep_real,                                                        &
        STASHwork8,                                                           &
        sf_diag                                                               &
        )
    END IF

    IF ( l_rivers .AND. trip_call .AND. sf(0,26) ) THEN
      CALL diagnostics_riv(                                                   &
        row_length, rows,                                                     &
        river_row_length, river_rows,                                         &
        riverout,                                                             &
        riverout_rgrid,                                                       &
        box_outflow, box_inflow,                                              &
        twatstor,inlandout_riv,                                               &
        STASHwork26                                                           &
        )
    END IF

    IF (sf(0,19)) THEN
      ! DEPENDS ON: diagnostics_veg
      CALL diagnostics_veg(                                                   &
        row_length, rows,                                                     &
        dim_cs1,                                                              &
        land_pts,                                                             &
        land_index,                                                           &
        ntype,npft,                                                           &
        c_veg_pft,cv_gb,g_leaf_phen_pft,                                      &
        lit_c_pft,lit_c_mn_gb,g_leaf_day_pft,                                 &
        lai_phen_pft,g_leaf_dr_out_pft,npp_dr_out_pft,                        &
        resp_w_dr_out_pft,resp_s_dr_out_gb_um,frac_agr_gb,disturb_veg_prev,   &
        wood_prod_fast_gb, wood_prod_med_gb, wood_prod_slow_gb,               &
        frac_surft,lai_pft,canht_pft,cs_pool_gb_um,                           &
        STASHwork19                                                           &
        )
    END IF
  END IF
#endif

  !=============================================================================
CASE ( cable )
#if defined(UM_JULES)
  errorstatus = 101
  CALL ereport('surf_couple_extra', errorstatus,                              &
               'CABLE not yet implemented')
#else

  ! initialise all INTENT(OUT) for now until CABLE is implemented
  melt_surft(:,:) = 0.0
  snomlt_surf_htf(:,:) = 0.0
  snowmelt_ij(:,:) = 0.0
  snomlt_sub_htf(:) = 0.0
  sub_surf_roff(:) = 0.0
  surf_roff(:) = 0.0
  tot_tfall(:) = 0.0
  snowmelt_gb(:) = 0.0
  rrun(:) = 0.0
  rflow(:) = 0.0
  snow_soil_htf(:,:) = 0.0
!CABLE_LSM: copied from JULES case
  !=============================================================================
  ! Redimensioning of arrays

#if defined(UM_JULES)
  !Dimensionality of variables differ betweeen UM and standalone, so copy across

  DO l = 1,land_pts
    DO n = 1,dim_cs1
      cs_pool_soilt(l,1,1,n) = cs_pool_gb_um(l,n)
    END DO
  END DO

  ! Change 2d to 1d soil clay content for soil respiration- dimensionalities
  ! differ between UM and standalone.
  ! Soil tiling not currently in the UM, so broadcast ij value to all tiles.
  ! Multi-layer clay not currently in UM so set all layers to same value.
  IF ( soil_bgc_model == soil_model_rothc ) THEN
    m = 1
    DO l = 1, land_pts
      j = (land_index(l) - 1) / row_length + 1
      i = land_index(l) - (j-1) * row_length
      DO n = 1, dim_cslayer
        clay_soilt(l,m,n) = soil_clay_ij(i,j)
      END DO
    END DO
  END IF

  !Compress pop_den and flash rate fields to land points- in the UM the ancil
  !is 2D
  IF (l_inferno) THEN
    DO l = 1, land_pts
      j = (land_index(l) - 1) / row_length + 1
      i = land_index(l) - (j-1) * row_length
      pop_den(l)    = pop_den_ancil(i,j)
      flash_rate(l) = flash_rate_ancil(i,j)
    END DO
  END IF
#endif

  !Encompassing IF statement for land_pts > 0. We exit/re-enter this IF about
  !half way down to allow for river routing

  IF (l_hydrology .AND. land_pts > 0 ) THEN

    !===========================================================================
    ! Compression to land points

    DO l = 1, land_pts
      j = (land_index(l) - 1) / row_length + 1
      i = land_index(l) - (j-1) * row_length
      ls_rain_gb(l)    = ls_rain_ij(i,j)
      con_rain_gb(l)   = con_rain_ij(i,j)
      con_snow_gb(l)   = con_snow_ij(i,j)
      ls_snow_gb(l)    = ls_snow_ij(i,j)
      ls_graup_gb(l)   = ls_graup_ij(i,j)
      pstar_gb(l)      = pstar_ij(i,j)
      tl_1_gb(l)       = tl_1_ij(i,j)
      qw_1_gb(l)       = qw_1_ij(i,j)
      u_1_gb(l)        = u_1_ij(i,j)
      v_1_gb(l)        = v_1_ij(i,j)
    END DO

    !Pass jules the modelled rain fractions
    !In standalone mode, both rainfracs are passed in as zeroed arrays
    IF (l_var_rainfrac) THEN
      DO l = 1, land_pts
        j = (land_index(l) - 1) / row_length + 1
        i = land_index(l) - (j-1) * row_length

        con_rainfrac_gb(l) = MIN(cca_2d(i,j),0.5) !As in diagnostics_conv

        !provide some safety checking for convective rain with no CCA
        IF (con_rain_gb(l) > 0.0 ) THEN
          IF (con_rainfrac_gb(l) == 0.0) THEN
            con_rainfrac_gb(l) = confrac
            !and for very small CCA amounts
          ELSE IF (con_rainfrac_gb(l) < 0.01) THEN
            con_rainfrac_gb(l) = 0.01
          END IF
        END IF

        !provide some safety checking for ls rain with no rainfrac
        IF (ls_rain_gb(l) > 0.0) THEN
          IF (ls_rainfrac_gb(l) == 0.0) THEN
            ls_rainfrac_gb(l) = 0.5
            !and for very small rainfrac amounts
          ELSE IF (ls_rainfrac_gb(l) < 0.01) THEN
            ls_rainfrac_gb(l) = 0.01
          END IF
        END IF
      END DO

    ELSE !use original default values
      DO l = 1, land_pts
        con_rainfrac_gb(l) = confrac
        ls_rainfrac_gb(l)  = 1.0
      END DO
    END IF !l_var_rainfrac

    !===========================================================================
    ! Science calls

    !Snow (standalone and UM)
!HaC!!    CALL snow ( land_pts,timestep_real,smlt,nsurft,surft_pts,                 &
!HaC!!                surft_index,catch_snow_surft,con_snow_gb,con_rain_gb,         &
!HaC!!                tile_frac,ls_snow_gb,ls_graup_gb,ls_rain_gb,                  &
!HaC!!                ei_surft,hcap_soilt(:,:,1),hcons_soilt,melt_surft,            &
!HaC!!                smcl_soilt(:,:,1),sthf_soilt(:,:,1),surf_htf_surft,           &
!HaC!!                t_soil_soilt(:,:,1),tsurf_elev_surft,                         &
!HaC!!                tstar_surft,smvcst_soilt(:,:,1),rgrain_surft,rgrainl_surft,   &
!HaC!!                rho_snow_grnd_surft,                                          &
!HaC!!                sice_surft,sliq_surft,snow_grnd_surft,snow_surft,             &
!HaC!!                snowdepth_surft, tsnow_surft,nsnow_surft,ds_surft,            &
!HaC!!                snomlt_surf_htf,snow_mass_gb,rho_snow_surft,snomlt_sub_htf,   &
!HaC!!                snowmelt_gb,snow_soil_htf,surf_ht_flux_ld,sf_diag,            &
!HaC!!                dhf_surf_minus_soil )
!HaC!!
    !Hydrology (standalone and UM)
    CALL hydrol_cbl (                                                             &
      lice_pts,lice_index,soil_pts,soil_index, nsnow_surft,                   &
      land_pts,sm_levels,bexp_soilt,catch_surft,con_rain_gb,                  &
      ecan_surft,ext_soilt,hcap_soilt,hcon_soilt,ls_rain_gb,                  &
      con_rainfrac_gb, ls_rainfrac_gb,                                        &
      satcon_soilt,sathh_soilt,snowdepth_surft, snow_soil_htf,                &
      surf_ht_flux_ld,timestep_real,                                          &
      smvcst_soilt,smvcwt_soilt,canopy_surft,                                 &
      stf_sub_surf_roff,smcl_soilt,sthf_soilt,sthu_soilt,                     &
      t_soil_soilt,tsurf_elev_surft,canopy_gb,smc_soilt,snowmelt_gb,          &
      sub_surf_roff,surf_roff,tot_tfall,                                      &
      ! add new inland basin variable
      inlandout_atm_gb,l_inland,                                              &
      ! Additional variables for MOSES II
      nsurft,surft_pts,surft_index,                                           &
      infil_surft, melt_surft,tile_frac,                                      &
      ! Additional variables required for large-scale hydrology:
      l_top,l_pdm,fexp_soilt,ti_mean_soilt,cs_ch4_soilt,cs_pool_soilt,        &
      dun_roff_soilt,drain_soilt,fsat_soilt,fwetl_soilt,qbase_soilt,          &
      qbase_l_soilt, qbase_zw_soilt, w_flux_soilt,                            &
      zw_soilt,sthzw_soilt,a_fsat_soilt,c_fsat_soilt,a_fwet_soilt,            &
      c_fwet_soilt,                                                           &
      resp_s_soilt,npp_gb,fch4_wetl_soilt,                                    &
      fch4_wetl_cs_soilt,fch4_wetl_npp_soilt,fch4_wetl_resps_soilt,           &
      dim_cs1,l_soil_sat_down,l_triffid,asteps_since_triffid,                 &
      !CABLE_LSM:additional existing vars
      snow_surft, snow_mass_gb ) 

  END IF ! ( l_hydrology .AND. land_pts /= 0 )

  !Here we need to exit the land_pts IF to allow river routing to be called
  !on all MPI ranks. This is because it is possible that a rank with
  !land_pts = 0 still has a river routing point.

  ! River Routing (standalone and UM, but differing options)
  IF ( l_rivers ) THEN
#if defined(UM_JULES)
    CALL river_control(                                                       &
      !LOGICAL, INTENT(IN)
      invert_ocean,                                                           &
      !INTEGER, INTENT(IN)
      n_proc, land_pts, row_length, rows, river_row_length, river_rows,       &
      land_index, ntype, aocpl_row_length, aocpl_p_rows, g_p_field,           &
      g_r_field, mype, global_row_length, global_rows,                        &
      global_river_row_length,                                                &
      global_river_rows, halo_i, halo_j, model_levels, nsurft,                &
      !REAL, INTENT(IN)
      fqw_surft, delta_lambda, delta_phi, xx_cos_theta_latitude,              &
      xpa, xua, xva, ypa, yua, yva, flandg, trivdir,                          &
      trivseq, r_area, slope, flowobs1, r_inext, r_jnext, r_land,             &
      smvcst_soilt, smvcwt_soilt,                                             &
      surf_roff, sub_surf_roff, frac_surft,                                   &
      !INTEGER, INTENT(INOUT)
      a_steps_since_riv,                                                      &
      !REAL, INTENT(INOUT)
      substore, surfstore, flowin, bflowin,                                   &
      tot_surf_runoff, tot_sub_runoff, acc_lake_evap, twatstor,               &
      smcl_soilt, sthu_soilt,                                                 &
      !LOGICAL, INTENT(OUT)
      trip_call,                                                              &
      !REAL, INTENT(OUT)
      inlandout_atm_gb, inlandout_riv, riverout,                              &
      riverout_rgrid, box_outflow, box_inflow, rflow, rrun                    &
      )
#else
    CALL river_control( land_pts,sub_surf_roff,surf_roff,rflow,rrun)
#endif
  END IF ! l_rivers (ATMOS)

  !Irrigation (standalone and UM for some options)
  IF ( l_irrig_dmd ) THEN
    CALL irrig_dmd(land_pts, sm_levels, frac_irr_soilt,                       &
                   a_step, plant_n_gb,                                        &
                   sthf_soilt, smvccl_soilt, smvcst_soilt, smvcwt_soilt,      &
                   sthzw_soilt, sthu_irr_soilt, sthu_soilt,                   &
                   smcl_soilt, irr_crop, dvi_cpft)

#if !defined(UM_JULES)
    ! Apply irrigation to planting dates as specified, apply limitations to 
    ! irrigation based on the water in rivers/groundwater
    ! irr_crop == 1 is not advisable for UM use. It requires detailed 
    ! information about planting dates to be calculated and there is a large 
    ! amount of technical debt which needs to be sorted before going into the UM
    IF ( irr_crop == 1 ) THEN
      CALL calc_crop_date(land_index, land_pts, t_i_length, t_j_length,       &
                          nsurft,                                             &
                          frac_surft, sw_surft, tstar_surft, lw_down, tl_1_ij,&
                          con_rain_ij, ls_rain_ij, con_snow_ij, ls_snow_ij,   &
                          plant_n_gb, nday_crop)
    END IF

    ! Limitation code requires l_irrig_dmd = TRUE, l_top = TRUE, l_rivers = TRUE
    ! and i_river_UM = rivers_trip (3) . The technical parts of this code are
    ! designed only to run with standalone TRIP routing code. UM TRIP/RFM code is
    ! deprecated.
    IF ( l_irrig_limit ) THEN
      CALL adjust_routestore()
    END IF
#endif

  ELSE
    ! Set sthu_irr_soilt to 0.0 in case it is still reported
    sthu_irr_soilt(:,:,:) = 0.0
  END IF ! l_irrig_dmd

  !Restart the IF for land_pts > 0 now that river routing is done with
  IF (land_pts > 0) THEN

    !Crops (standalone only)
#if !defined(UM_JULES)
    IF ( l_crop ) THEN
      crop_call = MOD ( REAL(a_step),                                         &
                        REAL(crop_period) * rsec_per_day / timestep_real )

      DO n = 1,ncpft
        DO l = 1,land_pts
          npp_acc_pft(l,nnpft + n) = npp_acc_pft(l,nnpft + n)                 &
                                + (npp_pft(l,nnpft + n) * timestep_real)
        END DO
      END DO

      CALL photoperiod(t_i_length * t_j_length, phot, dphotdt)

      CALL crop(t_i_length * t_j_length, land_pts, land_index, a_step,        &
                crop_call, sm_levels, frac_surft, phot, dphotdt,              &
                sf_diag%t1p5m_surft, t_soil_soilt, sthu_soilt, smvccl_soilt,  &
                smvcst_soilt, npp_acc_pft,                                    &
                canht_pft, lai_pft, dvi_cpft, rootc_cpft, harvc_cpft,         &
                reservec_cpft, croplai_cpft, cropcanht_cpft,                  &
                catch_surft, z0_surft)
    END IF  ! l_crop

    !Metstats (standalone only)
        !Beware- does not currently account for graupel
    IF ( l_metstats ) THEN
      CALL metstats_timestep(tl_1_gb, qw_1_gb, u_1_gb, v_1_gb, ls_rain_gb,    &
                            con_rain_gb, ls_snow_gb, con_snow_gb, pstar_gb,   &
                            metstats_prog,                                    &
                            !Vars that should be USED but can't due to
                            !UM/standalone differences
                            current_time%time, timestep_real,land_pts)
    END IF

    !Fire (standalone only)
    IF ( l_fire ) THEN
      !Calculate the gridbox mean soil moisture
      smc_gb = soiltiles_to_gbm(smc_soilt)
      CALL fire_timestep(metstats_prog, smc_gb, fire_prog, fire_diag,         &
                        !Vars that should be USED but can't due to
                        !UM/standalone differences
                        current_time%time, current_time%month, timestep_real, &
                        land_pts)
    END IF
#endif

    !INFERNO (standalone and UM)
    IF ( l_inferno ) THEN
      CALL calc_soil_carbon_pools(land_pts, soil_pts, soil_index, dim_cs1,    &
                                  cs_pool_soilt,                              &
                                  c_soil_dpm_gb, c_soil_rpm_gb)

      CALL inferno_io( sf_diag%t1p5m_surft, sf_diag%q1p5m_surft, pstar_gb,    &
                      sthu_soilt, sm_levels,                                  &
                      frac_surft, c_soil_dpm_gb, c_soil_rpm_gb, canht_pft,    &
                      ls_rain_gb, con_rain_gb,                                &
                      pop_den, flash_rate,                                    &
                      land_pts, ignition_method,                              &
                      nsurft, asteps_since_triffid)
    END IF  !  l_inferno

    !Vegetation (standalone and UM- differing functionality)
    IF (l_phenol .OR. l_triffid) THEN
      CALL veg_control(                                                       &
        land_pts, nsurft, dim_cs1,                                            &
        a_step, asteps_since_triffid,                                         &
        land_pts_trif, npft_trif,                                             &
        phenol_period, triffid_period,                                        &
        l_phenol, l_triffid, l_trif_eq,                                       &
        timestep_real, frac_agr_gb, frac_past_gb, satcon_soilt,               &
        g_leaf_acc_pft, g_leaf_phen_acc_pft, npp_acc_pft,                     &
        resp_s_acc_gb_um, resp_w_acc_pft,                                     &
        cs_pool_gb_um, frac_surft, lai_pft, clay_soilt, z0m_soil_gb,          &
        canht_pft,                                                            &
        catch_snow_surft, catch_surft, infil_surft, z0_surft, z0h_bare_surft, &
        c_veg_pft, cv_gb, lit_c_pft, lit_c_mn_gb, g_leaf_day_pft,             &
        g_leaf_phen_pft,                                                      &
        lai_phen_pft, g_leaf_dr_out_pft, npp_dr_out_pft, resp_w_dr_out_pft,   &
        resp_s_dr_out_gb_um, qbase_l_soilt, sthf_soilt, sthu_soilt,           &
        w_flux_soilt, t_soil_soilt,cs_pool_soilt)
    END IF

#if !defined(UM_JULES)
    !Reference evapotranspiration (standalone only)
    IF (l_fao_ref_evapotranspiration) THEN
      trad = ( surftiles_to_gbm(tstar_surft**4) )**0.25
      CALL fao_ref_evapotranspiration(soil_pts, soil_index,                   &
        land_pts, land_index, sf_diag%t1p5m,                                  &
        sw_down_ij, lw_down_ij, surf_ht_flux_ij, sf_diag%u10m,                &
        sf_diag%v10m, sf_diag%q1p5m, pstar_ij, trad, fao_et0)
    END IF
#endif

    !End of science calls
    !===========================================================================
    !Expansion to full grid

    DO l = 1, land_pts
      j=(land_index(l) - 1) / row_length + 1
      i = land_index(l) - (j-1) * row_length
      snow_mass_ij(i,j) = snow_mass_gb(l)
      snowmelt_ij(i,j)  = snowmelt_gb(l)
    END DO

  END IF ! land_pts > 0

  !=============================================================================
  !UM Diagnostics calls

    !For the UM, soil tiling has not been implemented, ie nsoilt = 1, so we can
    !hard-code _soilt variables with index 1 using m = 1
#if defined(UM_JULES)
  IF (model_type /= mt_single_column) THEN
    m = 1
    IF (l_hydrology .AND. sf(0,8) ) THEN
      ! DEPENDS ON: diagnostics_hyd
      CALL diagnostics_hyd(                                                   &
        row_length, rows, model_levels,                                       &
        n_rows, global_row_length, global_rows,                               &
        halo_i, halo_j, offx, offy, mype,                                     &
        n_proc, n_procx, n_procy,                                             &
        g_rows, g_row_length,                                                 &
        at_extremity,                                                         &
        land_pts, sm_levels,                                                  &
        !Put inland basin outflow in call to diagriv
        land_index,inlandout_atm_gb,                                          &
        smc_soilt(:,m), surf_roff, sub_surf_roff,                             &
        snow_mass_gb, snowmelt_gb,                                            &
        canopy_gb,t_soil_soilt(:,m,:),                                        &
        tsurf_elev_surft,snow_soil_htf,                                       &
        smcl_soilt(:,m,:),                                                    &
        nsurft, snomlt_surf_htf, sthu_soilt(:,m,:), sthf_soilt(:,m,:),        &
        tot_tfall, snow_surft, melt_surft,                                    &
        rgrain_surft, land_sea_mask,                                          &
        dun_roff_soilt(:,m), drain_soilt(:,m), qbase_soilt(:,m),              &
        qbase_zw_soilt(:,m), fch4_wetl_soilt(:,m),fch4_wetl_cs_soilt(:,m),    &
        fch4_wetl_npp_soilt(:,m),fch4_wetl_resps_soilt(:,m),                  &
        fexp_soilt(:,m),gamtot_soilt(:,m),ti_mean_soilt(:,m),                 &
        ti_sig_soilt(:,m),                                                    &
        fsat_soilt(:,m),fwetl_soilt(:,m),zw_soilt(:,m),sthzw_soilt(:,m),      &
        timestep_real,                                                        &
        STASHwork8,                                                           &
        sf_diag                                                               &
        )
    END IF

    IF ( l_rivers .AND. trip_call .AND. sf(0,26) ) THEN
      CALL diagnostics_riv(                                                   &
        row_length, rows,                                                     &
        river_row_length, river_rows,                                         &
        at_extremity,                                                         &
        at_extremity,                                                         &
        riverout,                                                             &
        riverout_rgrid,                                                       &
        box_outflow, box_inflow,                                              &
        !Put inland basin outflow in call to diagriv
        twatstor,inlandout_riv,                                               &
        STASHwork26                                                           &
        )
    END IF

    IF (sf(0,19)) THEN
      ! DEPENDS ON: diagnostics_veg
      CALL diagnostics_veg(                                                   &
        row_length, rows, n_rows,                                             &
        global_row_length, global_rows,                                       &
        dim_cs1, dim_cs2,                                                     &
        halo_i, halo_j, offx, offy, mype,                                     &
        n_proc, n_procx, n_procy,                                             &
        g_rows, g_row_length,                                                 &
        at_extremity,                                                         &
        land_pts,                                                             &
        land_index,                                                           &
        ntype,npft,                                                           &
        c_veg_pft,cv_gb,g_leaf_phen_pft,                                      &
        lit_c_pft,lit_c_mn_gb,g_leaf_day_pft,                                 &
        lai_phen_pft,g_leaf_dr_out_pft,npp_dr_out_pft,                        &
        resp_w_dr_out_pft,resp_s_dr_out_gb_um,frac_agr_gb,disturb_veg_prev,   &
        wood_prod_fast_gb, wood_prod_med_gb, wood_prod_slow_gb,               &
        frac_surft,lai_pft,canht_pft,cs_pool_gb_um,                           &
        STASHwork19                                                           &
        )
    END IF
  END IF
#endif


#endif
CASE DEFAULT
  errorstatus = 101
  CALL ereport('surf_couple_extra', errorstatus,                              &
               'Unrecognised surface scheme')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE surf_couple_extra
END MODULE surf_couple_extra_mod
