MODULE allocate_jules_arrays_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOCATE_JULES_ARRAYS_MOD'

CONTAINS
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine ALLOCATE_JULES_ARRAYS
!
! Description: Routine that allocates memory to the JULES arrays
! This assume that the values in the jules_surface_types module have been set
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: Please refer to ModuleLeaders.txt and UM file CodeOwners.txt
!   This file belongs in section: Land

#if defined(UM_JULES)
SUBROUTINE allocate_jules_arrays(                                             &
  land_pts,nsurft,sm_levels,nice,nice_use                                     &
  )
#else
SUBROUTINE allocate_jules_arrays()

!Replacements for the argument list
USE ancil_info,               ONLY:                                           &
  land_pts, nsurft
USE jules_soil_mod,           ONLY:                                           &
  sm_levels
USE jules_rivers_mod,         ONLY:                                           &
   tot_surf_runoff_gb, tot_sub_runoff_gb, acc_lake_evap_gb
USE jules_sea_seaice_mod,     ONLY:                                           &
  nice, nice_use
#endif

!Common Non-science modules
USE parkind1,                 ONLY:                                           &
  jprb, jpim
USE yomhook,                  ONLY:                                           &
  lhook, dr_hook
USE jules_print_mgr,          ONLY:                                           &
  jules_message, jules_print, PrNorm
USE ereport_mod,              ONLY:                                           &
  ereport
USE ancil_info,               ONLY:                                           &
  dim_cslayer

!USE statements that apply to both UM and JULES
USE ancil_info,               ONLY:                                           &
  ssi_index, sea_index, sice_index, sice_pts_ncat, fssi_ij,                   &
  sea_frac, sice_frac, sice_frac_ncat, sice_index_ncat,                       &
  l_lice_point, l_lice_surft,                                                 &
  l_soil_point, nsoilt, soilt_pts, soilt_index, frac_soilt, land_pts_trif,    &
  npft_trif, dim_cs1

USE bvoc_vars,                ONLY:                                           &
  isoprene_gb, isoprene_pft, terpene_gb , terpene_pft,                        &
  methanol_gb, methanol_pft, acetone_gb, acetone_pft

USE c_elevate,                ONLY:                                           &
  surf_hgt_surft, lw_down_elevcorr_surft

USE crop_vars_mod,            ONLY:                                           &
  sthu_irr_soilt, frac_irr_all, frac_irr_soilt, frac_irr_old_soilt,           &
  frac_irr_surft,                                                             &
  plant_n_gb, smc_irr_soilt, wt_ext_irr_surft, gs_irr_surft, dvimax_gb,       &
  gc_irr_surft, resfs_irr_surft, ext_irr_soilt, wt_ext_irr_gb, fsmc_irr_gb,   &
  irrdaysdiag_gb, icntmax_gb, tl_1_day_av_gb, tl_1_day_av_use_gb,             &
  rn_1_day_av_gb, rn_1_day_av_use_gb, ndpy, nyav, prec_1_day_av_gb,           &
  prec_1_day_av_use_gb, irrig_water_gb, irrtiles

USE c_z0h_z0m,                ONLY:                                           &
  z0h_z0m, z0h_z0m_classic

USE fire_vars,                ONLY:                                           &
  burnt_area, burnt_area_ft,                                                  &
  emitted_carbon, emitted_carbon_ft, emitted_carbon_DPM, emitted_carbon_RPM,  &
  fire_em_CO2, fire_em_CO2_ft, fire_em_CO2_DPM, fire_em_CO2_RPM,              &
  fire_em_CO, fire_em_CO_ft, fire_em_CO_DPM, fire_em_CO_RPM,                  &
  fire_em_CH4, fire_em_CH4_ft, fire_em_CH4_DPM, fire_em_CH4_RPM,              &
  fire_em_NOx, fire_em_NOx_ft, fire_em_NOx_DPM, fire_em_NOx_RPM,              &
  fire_em_SO2, fire_em_SO2_ft, fire_em_SO2_DPM, fire_em_SO2_RPM,              &
  fire_em_OC,  fire_em_OC_ft, fire_em_OC_DPM,  fire_em_OC_RPM,                &
  fire_em_BC, fire_em_BC_ft, fire_em_BC_DPM, fire_em_BC_RPM,                  &
  pop_den, flash_rate, flammability_ft

USE fluxes,                   ONLY:                                           &
  anthrop_heat_surft, surf_ht_store_surft,                                    &
  sw_sicat, sw_rts_sicat, swup_rts_sicat, swdn_rts_sicat, alb_sicat,          &
  sw_sea, sw_rts_sea

!These should be spirited away to other modules, allowing it to be retired
USE jules_mod,                ONLY:                                           &
  snowdep_surft, albobs_scaling_surft

USE jules_radiation_mod,      ONLY:                                           &
  l_albedo_obs, l_spec_albedo

USE jules_snow_mod,           ONLY:                                           &
  nsmax, cansnowtile

USE jules_internal,           ONLY:                                           &
  unload_backgrnd_pft

USE jules_soil_biogeochem_mod, ONLY:                                          &
  soil_model_ecosse, soil_model_rothc, soil_bgc_model, l_layeredc

USE jules_soil_mod,           ONLY:                                           &
  l_bedrock, ns_deep

USE jules_surface_mod,        ONLY:                                           &
  diff_frac, l_urban2t

USE jules_surface_types_mod,  ONLY:                                           &
  npft, nnvg, ntype

USE jules_vegetation_mod,     ONLY:                                           &
  l_triffid, l_phenol, irr_crop, l_nitrogen, l_irrig_dmd, l_use_pft_psi

USE nvegparm,                 ONLY:                                           &
  albsnc_nvg, albsnf_nvgu, albsnf_nvg, albsnf_nvgl, catch_nvg, emis_nvg,      &
  gs_nvg, infil_nvg, z0_nvg, ch_nvg, vf_nvg

USE ozone_vars,               ONLY:                                           &
  o3_gb, flux_o3_pft, fo3_pft

USE pftparm,                  ONLY:                                           &
  act_jmax, act_vcmax, albsnc_max, albsnc_min, albsnf_maxu, albsnf_max,       &
  albsnf_maxl, alpha, alpha_elec, alniru, alnir, alnirl, alparu, alpar,       &
  alparl, a_wl, a_ws, b_wl, catch0, c3, dcatch_dlai, deact_jmax, deact_vcmax, &
  dgl_dm, dgl_dt, dqcrit, ds_jmax, ds_vcmax, dz0v_dh, z0v, emis_pft,          &
  eta_sl, fd, fsmc_of, f0, g1_stomata, glmin, g_leaf_0, infil_f, jv25_ratio,  &
  kext, kpar, lai_alb_lim, neff,                                              &
  nl0, nr_nl, ns_nl, omegau, omega, omegal, omniru, omnir, omnirl, orient,    &
  r_grow, rootd_ft, sigl, tleaf_of, tlow,  tupp, lma, nmass, vsl, vint, kn,   &
  knl, q10_leaf, fl_o3_ct, dfp_dcuo, ci_st, gpp_st, ief, tef, mef, aef,       &
  fef_co2, fef_co, fef_ch4, fef_nox, fef_so2, fef_oc, fef_bc, ccleaf_min,     &
  ccleaf_max, ccwood_min, ccwood_max, avg_ba, fire_mort,                      &
  nsw, nr, hw_sw, psi_close, psi_open, fsmc_p0, fsmc_mod, can_struct_a

USE prognostics,              ONLY:                                           &
  lai_pft, canht_pft, smcl_soilt, t_soil_soilt, tsurf_elev_surft,             &
  rgrain_surft, snow_surft, soot_ij, t_soil_soilt_acc,                        &
  tstar_surft, canopy_surft, canopy_gb, cs_pool_soilt, ti_sicat, z0msea_ij,   &
  gs_gb, gc_surft, smc_soilt, di_ncat_sicat, k_sice_sicat, snow_grnd_surft,   &
  snow_mass_ij, snow_mass_sea_sicat, nsnow_surft, rho_snow_grnd_surft,        &
  snowdepth_surft, ds_surft, rgrainl_surft, sice_surft, sliq_surft,           &
  tsnow_surft, wood_prod_fast_gb, wood_prod_med_gb, wood_prod_slow_gb,        &
  frac_agr_prev_gb, frac_past_prev_gb,                                        &
  !This is required in the UM dump file but not for standalone.
  !It may be redundant information and therefore could be removed
  rho_snow_surft, tsoil_deep_gb, ns_pool_gb, n_inorg_gb, n_inorg_soilt_lyrs,  &
  n_inorg_avail_pft, triffid_co2_gb

USE p_s_parms,                ONLY:                                           &
  bexp_soilt, sathh_soilt, hcap_soilt, hcon_soilt, satcon_soilt,              &
  smvccl_soilt, smvcwt_soilt, smvcst_soilt, clay_soilt

USE switches_urban,           ONLY:                                           &
  l_moruses

USE theta_field_sizes,        ONLY:                                           &
  t_i_length, t_j_length

USE top_pdm,                  ONLY:                                           &
  a_fsat_soilt, a_fwet_soilt, c_fsat_soilt, c_fwet_soilt, drain_soilt,        &
  dun_roff_soilt, fexp_soilt, fsat_soilt, fch4_wetl_soilt,                    &
  fch4_wetl_cs_soilt, fch4_wetl_npp_soilt, fch4_wetl_resps_soilt,             &
  fwetl_soilt, gamtot_soilt, qbase_soilt, qbase_zw_soilt, sthzw_soilt,        &
  ti_mean_soilt, ti_sig_soilt, zw_soilt, inlandout_atm_gb,                    &
  fch4_wetl_acc_soilt

USE trif,                     ONLY:                                           &
  crop, g_area, g_grow, g_root, g_wood, lai_max, lai_min, alloc_fast,         &
  alloc_med, alloc_slow, dpm_rpm_ratio, retran_r, retran_l

USE trifctl,                  ONLY:                                           &
  lai_phen_pft, c_veg_pft, cv_gb, g_leaf_day_pft, g_leaf_dr_out_pft,          &
  lit_c_pft, lit_c_mn_gb, npp_dr_out_pft, resp_w_dr_out_pft,                  &
  resp_s_dr_out_gb, frac_agr_gb, g_leaf_acc_pft, npp_acc_pft,                 &
  resp_w_acc_pft, resp_s_acc_soilt, g_leaf_phen_acc_pft, gpp_gb, npp_gb,      &
  resp_p_gb, g_leaf_pft, g_leaf_phen_pft, gpp_pft, npp_pft, resp_p_pft,       &
  resp_s_soilt, resp_w_pft

USE trif_vars_mod,            ONLY:                                           &
  wp_fast_in_gb, wp_med_in_gb, wp_slow_in_gb, wp_fast_out_gb, wp_med_out_gb,  &
  wp_slow_out_gb, lit_c_orig_pft, lit_c_ag_pft, n_leaf_pft, n_root_pft,       &
  n_stem_pft, resp_r_pft, resp_l_pft, lai_bal_pft, pc_s_pft, resp_s_diag_gb,  &
  resp_s_pot_diag_gb, minl_n_gb, minl_n_pot_gb, immob_n_gb, immob_n_pot_gb,   &
  fn_gb, leafC_pft, rootC_pft, stemC_pft, woodC_pft, droot_pft, dleaf_pft,    &
  dwood_pft, n_uptake_extract, n_uptake_pft, n_demand_gb,                     &
  n_uptake_growth_pft,                                                        &
  n_demand_growth_pft, n_demand_spread_pft, n_uptake_spread_pft,              &
  n_demand_lit_pft, n_uptake_gb, n_demand_pft, exudates_pft, exudates_gb,     &
  dcveg_pft, dnveg_pft, dcveg_gb, dnveg_gb, n_veg_pft, n_veg_gb, n_loss_gb,   &
  dpm_ratio_gb, root_litC_pft, leaf_litC_pft, wood_litC_pft, root_litN_pft,   &
  leaf_litN_pft, wood_litN_pft, litterC_pft, litterN_pft, lit_n_pft,          &
  lit_n_t_gb, deposition_N_gb, n_fix_add, n_fix_gb, n_fix_pft,                &
  n_leach_soilt, n_gas_gb,                                                    &
  fapar_diag_pft, apar_diag_pft, apar_diag_gb, fao_et0,                       &
  cnsrv_carbon_veg2_gb, cnsrv_carbon_triffid_gb,                              &
  cnsrv_veg_triffid_gb, cnsrv_soil_triffid_gb, cnsrv_prod_triffid_gb,         &
  cnsrv_nitrogen_triffid_gb,                                                  &
  cnsrv_vegN_triffid_gb, cnsrv_soilN_triffid_gb, cnsrv_N_inorg_triffid_gb,    &
  root_abandon_pft, harvest_pft, harvest_gb, root_abandon_n_pft,              &
  harvest_n_pft, harvest_n_gb, n_fertiliser_add, n_fertiliser_pft,            &
  n_fertiliser_gb,                                                            &
  frac_past_gb, lit_n_orig_pft, lit_n_ag_pft, npp_n_gb, npp_n,                &
  lit_c_fire_pft, lit_c_nofire_pft, lit_n_fire_pft, lit_n_nofire_pft,         &
  veg_c_fire_emission_gb, veg_c_fire_emission_pft, g_burn_pft,                &
  n_leaf_trif_pft, n_leaf_alloc_trif_pft, n_leaf_labile_trif_pft,             &
  n_root_trif_pft, n_stem_trif_pft, lit_n_ag_pft_diag, n_luc,                 &
  lit_n_pft_diag, leafC_gbm, woodC_gbm, rootC_gbm, root_abandon_gb,           &
  root_abandon_n_gb, resp_s_to_atmos_gb, g_burn_gb, burnt_carbon_dpm,         &
  burnt_carbon_rpm, g_burn_pft_acc,                                           &
  gpp_gb_out, gpp_pft_out, gpp_gb_acc, gpp_pft_acc,                           &
  resp_p_actual_gb, resp_p_actual_pft, n_leach_gb_acc

USE urban_param,              ONLY:                                           &
   hgt_gb, hwr_gb, wrr_gb, disp_gb, ztm_gb, albwl_gb, albrd_gb, emisw_gb,     &
   emisr_gb

!Model-dependent USE statements
#if defined(UM_JULES)
  !Used for getting t_i_length and t_jlength and for FLAKE allocations
USE atm_fields_bounds_mod,    ONLY:                                           &
  tdims

!For FLAKE
USE jules_surface_mod,        ONLY:                                           &
  l_flake_model

!For FLAKE
USE lake_mod,                 ONLY:                                           &
  surf_ht_flux_lake_ij, surf_ht_flux_lk_gb, sw_down_gb, coriolis_param_gb,    &
  u_s_lake_gb, lake_depth_gb, lake_fetch_gb, lake_albedo_gb, lake_t_snow_gb,  &
  lake_t_ice_gb, lake_t_mean_gb, lake_t_mxl_gb, lake_shape_factor_gb,         &
  lake_h_snow_gb, lake_h_ice_gb, lake_h_mxl_gb, lake_t_sfc_gb, ts1_lake_gb,   &
  nusselt_gb, g_dt_gb

USE pftparm,                  ONLY:                                           &
  !dust_veg_scj is only used in the boundary layer scheme?
  dust_veg_scj

#else
USE update_mod, ONLY: sthuf_soilt

USE aero,                     ONLY:                                           &
  co2_3d_ij, rho_cd_modv1_ij, rho_aresist_ij, aresist_ij, resist_b_ij,        &
  rho_aresist_surft, aresist_surft, resist_b_surft, r_b_dust_ij,              &
  cd_std_dust_ij, u_s_std_surft

USE ancil_info,               ONLY:                                           &
  land_pts, surft_index, soil_index, lice_index, ice_fract_ij,                &
  ice_fract_ncat_sicat,                                                       &
  ti_cat_sicat, pond_frac_cat_sicat, pond_depth_cat_sicat, sstfrz_ij,         &
  z1_uv_ij, z1_tq_ij, nsurft, frac_surft, surft_pts

USE boundary_layer_var_mod,   ONLY:                                           &
  dzl, zh

USE c_elevate,                ONLY:                                           &
  z_land_ij

USE coastal,                  ONLY:                                           &
  tstar_land_ij, tstar_sea_ij, tstar_sice_ij, tstar_sice_sicat,               &
  tstar_ssi_ij, taux_land_ij, taux_ssi_ij, tauy_land_ij, tauy_ssi_ij,         &
  vshr_land_ij, vshr_ssi_ij, surf_ht_flux_land_ij, surf_ht_flux_sice_sicat,   &
  taux_land_star, tauy_land_star, taux_ssi_star, tauy_ssi_star

USE cropparm,                 ONLY:                                           &
  beta1, beta2, beta3, r_gamma, delta, remob, cfrac_s, cfrac_r, cfrac_l,      &
  allo1, allo2, t_bse, t_opt, t_max, tt_emr, crit_pp, pp_sens, rt_dir,        &
  alpha1, alpha2, alpha3, mu, nu, yield_frac, initial_carbon,                 &
  initial_c_dvi, sen_dvi, t_mort

USE crop_vars_mod,            ONLY:                                           &
  phot, dphotdt, dvi_cpft, croprootc_cpft => rootc_cpft, harvc_cpft,          &
  reservec_cpft, yield_diag_cpft, nonyield_diag_cpft, leafc_diag_cpft,        &
  stemc_diag_cpft, croplai_cpft, cropcanht_cpft, sow_date_cpft,               &
  tt_veg_cpft, tt_rep_cpft, latestharv_date_cpft,                             &
  harvest_trigger_cpft, harvest_counter_cpft

USE dust_parameters_mod,      ONLY:                                           &
  ndiv

USE fluxes,                   ONLY:                                           &
  alb_surft, e_sea_ij, ecan_ij, ecan_surft, ei_ij, ei_surft, esoil_ij_soilt,  &
  esoil_surft, ext_soilt, fqw_surft, fqw_sicat, fsmc_pft,                     &
  ftl_sicat, ftl_surft, h_sea_ij, hf_snow_melt_gb, land_albedo_ij,            &
  le_surft, melt_surft, sea_ice_htf_sicat, snomlt_sub_htf_gb, snow_melt_gb,   &
  snowmelt_ij, sub_surf_roff_gb, surf_ht_flux_ij, surf_htf_surft,             &
  surf_roff_gb, radnet_surft, tot_tfall_gb, tstar_ij,                         &
  sw_surft, emis_surft, alb_sicat, rflow_gb, rrun_gb,                         &
  snow_soil_htf, z0m_surft, z0h_surft

USE gridmean_fluxes,          ONLY:                                           &
  fqw_1_ij, ftl_1_ij, taux_1_ij, tauy_1_ij

USE forcing,                  ONLY:                                           &
  qw_1_ij, tl_1_ij, u_0_ij, v_0_ij, u_1_ij, v_1_ij, pstar_ij, ls_rain_ij,     &
  con_rain_ij, ls_snow_ij, con_snow_ij, sw_down_ij, lw_down_ij, diff_rad_ij,  &
  diurnal_temperature_range_ij

USE jules_deposition_mod,     ONLY:                                           &
! imported parameters
  dry_dep_ukca_jules,                                                         &
! imported variables
  dry_dep_model, l_deposition, l_deposition_flux, ndry_dep_species,           &
  tracer_field

USE deposition_species_mod,   ONLY:                                           &
  dd_ice_coeff, dep_species_name, diffusion_coeff, diffusion_corr,            &
  h2dd_c, h2dd_m, h2dd_q, r_tundra, rsurf_std

USE jules_sea_seaice_mod,     ONLY:                                           &
  nice, nice_use

USE jules_soil_mod,           ONLY:                                           &
  sm_levels

USE jules_surface_types_mod,  ONLY:                                           &
  ncpft

USE jules_vegetation_mod,     ONLY:                                           &
  l_crop

USE max_dimensions,           ONLY:                                           &
  npft_max

USE nlsizes_namelist_mod,     ONLY:                                           &
  bl_levels

USE orog,                     ONLY:                                           &
  sil_orog_land_gb, ho2r2_orog_gb, h_blend_orog_ij, z0m_eff_ij

USE p_s_parms,                ONLY:                                           &
  albsoil_soilt, albobs_sw_gb, albobs_vis_gb, albobs_nir_gb, catch_surft,     &
  catch_snow_surft, cosz_ij, infil_surft, z0_surft, z0h_bare_surft,           &
  z0m_soil_gb, sthu_soilt, sthf_soilt, sthu_min_soilt,                        &
  soil_ph_soilt, v_close_pft, v_open_pft

USE soil_ecosse_vars_mod,     ONLY:                                           &
  co2_soil_gb, deposition_n_driver, n2_denitrif_gb, n2o_denitrif_gb,          &
  n2o_nitrif_gb,  n2o_partial_nitrif_gb, n2o_soil_gb, n_amm_soilt,            &
  n_denitrification_gb, n_leach_amm_gb, n_leach_nit_gb,                       &
  n_nit_soilt, n_nitrification_gb, no_soil_gb,                                &
  plant_input_c_gb, plant_input_n_gb, qbase_l_driver, sthf_driver,            &
  sthu_driver, tsoil_driver, wflux_driver, soil_c_add, soil_n_add

USE pdm_vars,                 ONLY:                                           &
  slope_gb

USE u_v_grid,                 ONLY:                                           &
  u_0_p_ij, v_0_p_ij, u_1_p_ij, v_1_p_ij, dtrdz_charney_grid_1_ij

#endif

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Allocates the model arrays using sizes determined during initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!
! Modules are ordered alphabetically. Please respect it if you add something
! new.
!
! The allocation statements are also grouped into 'JULES', 'UM' and 'Common'
! according to whether they are needed for just JULES, the UM or both.
!
!-----------------------------------------------------------------------------

#if defined(UM_JULES)
! Input variables for dimensioning
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
                       ! Number of land points
  nsurft,                                                                     &
                       ! Number of surface tiles
  sm_levels,                                                                  &
                       ! Number of soil layers
  nice,                                                                       &
                       ! Number of sea ice categories
  nice_use
                       ! Number of sea ice cats used in radiation and
                       !  explicit part of surface exchange
#endif

!-----------------------------------------------------------------------
! Local variables for error trapping
!-----------------------------------------------------------------------
INTEGER ::                                                                    &
  error        = 0,                                                           &
                       ! Variable for trapping the error from each
                       ! individual call to allocate
  error_sum    = 0,                                                           &

                       ! Variable to track the sum of all errors
                       ! resulting from calls to allocate. Hence we
                       ! know that everything was successful if and
                       ! only if this is zero at the end

  temp_size,                                                                  &
  temp_tiles,                                                                 &
  temp_layers,                                                                &
                       ! For storing the size of array to allocate for variables
                       ! that are sometimes set to size 1.
                       ! Removes some duplicate allocate statements
  errcode 
                       ! Variable to use in error report

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_JULES_ARRAYS'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
errcode = 101

#if defined(UM_JULES)
!DO WE NEED THIS? Only existed in the UM version.
!-----------------------------------------------------------------------
! Compute length of theta field in i and j directions.
! (The module keeps these values for future use, this is
!  earliest place in the code that they are needed.)
!-----------------------------------------------------------------------
t_i_length = tdims%i_end - tdims%i_start + 1
t_j_length = tdims%j_end - tdims%j_start + 1

!We need to set the sizes of a few dimensions. These are also exist in 
!atm_step_local and set by the same logic in atm_step_phys_init. However,
!This happens after the call to this routine.
!Pending a better approach to JULES memory management we'll repeat the logic 
!here

IF (l_triffid) THEN
  land_pts_trif = land_pts
  npft_trif = npft
  dim_cs1 = 4
ELSE
  land_pts_trif = 1
  npft_trif = 1
  dim_cs1 = 1
END IF

#endif

!------------------------------------------------------------------------------
!Begin the array allocation. We will follow the order that the USE statements
!have been done.
!------------------------------------------------------------------------------

!USE statements that apply to both UM and JULES

!  ====ancil_info module common====
ALLOCATE( ssi_index(t_i_length * t_j_length), stat = error )
error_sum = error
ALLOCATE( sea_index(t_i_length * t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( sice_index(t_i_length * t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( sice_pts_ncat(nice), stat = error )
error_sum = error_sum + error
ALLOCATE( fssi_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( sea_frac(t_i_length * t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( sice_frac(t_i_length * t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( sice_frac_ncat(t_i_length * t_j_length,nice), stat = error )
error_sum = error_sum + error
ALLOCATE( sice_index_ncat(t_i_length * t_j_length,nice), stat = error )
error_sum = error_sum + error
ALLOCATE( l_lice_point(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( l_soil_point(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( soilt_pts(nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( soilt_index(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( frac_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  ssi_index(:)         = 0
  sea_index(:)         = 0
  sice_index(:)        = 0
  sice_pts_ncat(:)     = 0
  fssi_ij(:,:)         = 0.0
  sea_frac(:)          = 0.0
  sice_frac(:)         = 0.0
  sice_frac_ncat(:,:)  = 0.0
  sice_index_ncat(:,:) = 0
  l_lice_point(:)      = .FALSE.
  l_soil_point(:)      = .FALSE.
  soilt_pts(:)         = 0
  soilt_index(:,:)     = 0
  frac_soilt(:,:)      = 0.0
ELSE
  CALL ereport ("Something in ancil_info module common causing error.",       &
                errcode, "please check allocate_jules_arrays")
END IF
ALLOCATE( l_lice_surft(ntype), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  l_lice_surft(:)      = .FALSE.
END IF

!  ====bvoc_vars module common====
! BVOC diagnostics
ALLOCATE( isoprene_gb(land_pts), stat = error )
error_sum = error
ALLOCATE( isoprene_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( terpene_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( terpene_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( methanol_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( methanol_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( acetone_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( acetone_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  isoprene_gb(:)    = 0.0
  isoprene_pft(:,:) = 0.0
  terpene_gb(:)     = 0.0
  terpene_pft(:,:)  = 0.0
  methanol_gb(:)    = 0.0
  methanol_pft(:,:) = 0.0
  acetone_gb(:)     = 0.0
  acetone_pft(:,:)  = 0.0
ELSE
  CALL ereport ("Something in bvoc_vars module common is causing error.",     &
                errcode, "please check allocate_jules_arrays")
END IF

!  ====c_elevate module common====
! Height above mean grid-box
ALLOCATE( surf_hgt_surft(land_pts,nsurft), stat = error )
error_sum = error
ALLOCATE( lw_down_elevcorr_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  surf_hgt_surft(:,:)         = 0.0
  lw_down_elevcorr_surft(:,:) = 0.0
ELSE
  CALL ereport ("Something in c_elevate module common is causing error.",     &
                errcode, "please check allocate_jules_arrays")
END IF

!  ====crop_vars_mod module common====
! Irrigation variables
ALLOCATE( sthu_irr_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error
ALLOCATE( frac_irr_all(land_pts,1), stat = error )
error_sum = error_sum + error
ALLOCATE( frac_irr_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( frac_irr_old_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( frac_irr_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( plant_n_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( smc_irr_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( wt_ext_irr_surft(land_pts,sm_levels,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( gs_irr_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( dvimax_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( gc_irr_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( resfs_irr_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( ext_irr_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE( wt_ext_irr_gb(land_pts,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE( fsmc_irr_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( irrDaysDiag_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( irrig_water_gb(land_pts), stat = error )
error_sum = error_sum + error

IF ( error_sum == 0 ) THEN
  sthu_irr_soilt(:,:,:)   = 0.0
  frac_irr_all(:,:)       = 0.0
  frac_irr_soilt(:,:)     = 0.0
  frac_irr_old_soilt(:,:) = 0.0
  frac_irr_surft(:,:)     = 0.0
  plant_n_gb(:)           = 0
  smc_irr_soilt(:,:)      = 0.0
  wt_ext_irr_surft(:,:,:) = 0.0
  gs_irr_surft(:,:)       = 0.0
  dvimax_gb(:)            = 0.0
  gc_irr_surft(:,:)       = 0.0
  resfs_irr_surft(:,:)    = 0.0
  ext_irr_soilt(:,:,:)    = 0.0
  wt_ext_irr_gb(:,:)      = 0.0
  fsmc_irr_gb(:)          = 0.0
  irrDaysDiag_gb(:)       = 0.0
  irrig_water_gb(:)       = 0.0
ELSE
  CALL ereport ("Something in crop_vars_mod module common causing error.",    &
                errcode, "please check allocate_jules_arrays")
END IF

!-----------------------------------------------------------------------
! Allocate space for inferno diagnostic variables
!-----------------------------------------------------------------------
ALLOCATE( burnt_area(land_pts), stat = error )
error_sum = error
ALLOCATE( burnt_area_ft(land_pts,npft), stat = error )
error_sum = error_sum + error

ALLOCATE( emitted_carbon(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( emitted_carbon_ft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( emitted_carbon_DPM(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( emitted_carbon_RPM(land_pts), stat = error )
error_sum = error_sum + error

ALLOCATE( fire_em_CO2(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_CO2_ft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_CO2_DPM(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_CO2_RPM(land_pts), stat = error )
error_sum = error_sum + error

ALLOCATE( fire_em_CO(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_CO_ft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_CO_DPM(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_CO_RPM(land_pts), stat = error )
error_sum = error_sum + error

ALLOCATE( fire_em_CH4(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_CH4_ft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_CH4_DPM(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_CH4_RPM(land_pts), stat = error )
error_sum = error_sum + error

ALLOCATE( fire_em_NOx(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_NOx_ft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_NOx_DPM(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_NOx_RPM(land_pts), stat = error )
error_sum = error_sum + error

ALLOCATE( fire_em_SO2(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_SO2_ft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_SO2_DPM(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_SO2_RPM(land_pts), stat = error )
error_sum = error_sum + error

ALLOCATE( fire_em_OC(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_OC_ft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_OC_DPM(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_OC_RPM(land_pts), stat = error )
error_sum = error_sum + error

ALLOCATE( fire_em_BC(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_BC_ft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_BC_DPM(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_em_BC_RPM(land_pts), stat = error )
error_sum = error_sum + error

ALLOCATE( pop_den(land_pts), stat = error )
error_sum = error_sum + error

ALLOCATE( flash_rate(land_pts), stat = error )
error_sum = error_sum + error

ALLOCATE( flammability_ft(land_pts,npft), stat = error )
error_sum = error_sum + error

IF ( error_sum == 0 ) THEN
  burnt_area(:)          = 0.0
  burnt_area_ft(:,:)     = 0.0
  emitted_carbon(:)      = 0.0
  emitted_carbon_ft(:,:) = 0.0
  emitted_carbon_DPM(:)  = 0.0
  emitted_carbon_RPM(:)  = 0.0
  fire_em_CO2(:)         = 0.0
  fire_em_CO2_ft(:,:)    = 0.0
  fire_em_CO2_DPM(:)     = 0.0
  fire_em_CO2_RPM(:)     = 0.0
  fire_em_CO(:)          = 0.0
  fire_em_CO_ft(:,:)     = 0.0
  fire_em_CO_DPM(:)      = 0.0
  fire_em_CO_RPM(:)      = 0.0
  fire_em_CH4(:)         = 0.0
  fire_em_CH4_ft(:,:)    = 0.0
  fire_em_CH4_DPM(:)     = 0.0
  fire_em_CH4_RPM(:)     = 0.0
  fire_em_NOx(:)         = 0.0
  fire_em_NOx_ft(:,:)    = 0.0
  fire_em_NOx_DPM(:)     = 0.0
  fire_em_NOx_RPM(:)     = 0.0
  fire_em_SO2(:)         = 0.0
  fire_em_SO2_ft(:,:)    = 0.0
  fire_em_SO2_DPM(:)     = 0.0
  fire_em_SO2_RPM(:)     = 0.0
  fire_em_OC(:)          = 0.0
  fire_em_OC_ft(:,:)     = 0.0
  fire_em_OC_DPM(:)      = 0.0
  fire_em_OC_RPM(:)      = 0.0
  fire_em_BC(:)          = 0.0
  fire_em_BC_ft(:,:)     = 0.0
  fire_em_BC_DPM(:)      = 0.0
  fire_em_BC_RPM(:)      = 0.0
  pop_den(:)             = 0.0
  flash_rate(:)          = 0.0
  flammability_ft(:,:)   = 0.0
ELSE
  CALL ereport ("Something in inferno diagnostic variables causing error.",   &
                 errcode, "please check allocate_jules_arrays")
END IF

! vars for Doell and Siebert crop calendar
IF ( irr_crop == 1 ) THEN
  ALLOCATE( icntmax_gb(land_pts), stat = error )
  error_sum = error
  ALLOCATE( tl_1_day_av_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( tl_1_day_av_use_gb(land_pts,ndpy,nyav), stat = error )
  error_sum = error_sum + error
  ALLOCATE( prec_1_day_av_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( prec_1_day_av_use_gb(land_pts,ndpy,nyav), stat = error )
  error_sum = error_sum + error
  ALLOCATE( rn_1_day_av_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( rn_1_day_av_use_gb(land_pts,ndpy,nyav), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    icntmax_gb(:)                = 0
    tl_1_day_av_gb(:)            = 0.0
    tl_1_day_av_use_gb(:,:,:)    = 0.0
    prec_1_day_av_gb(:)          = 0.0
    prec_1_day_av_use_gb(:,:,:)  = 0.0
    rn_1_day_av_gb(:)            = 0.0
    rn_1_day_av_use_gb(:,:,:)    = 0.0
  ELSE
    CALL ereport ("Something in vars for Doell and Siebert crop calendar "//  &
                  "causing error.", errcode, "please check "  //              &
                  "allocate_jules_arrays")
  END IF
END IF

!  ====c_z0h_z0m module, common
! Surface type variables.
ALLOCATE( z0h_z0m(ntype), stat = error )
error_sum = error
ALLOCATE( z0h_z0m_classic(ntype), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  z0h_z0m(:)         = 0.0
  z0h_z0m_classic(:) = 0.0
ELSE
  CALL ereport ("Something in c_z0h_z0m module, common causing error.",       &
                 errcode, "please check allocate_jules_arrays")
END IF

!  ====fluxes module, common
ALLOCATE( surf_ht_store_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( anthrop_heat_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( sw_rts_sicat(t_i_length * t_j_length, nice_use), stat = error )
error_sum = error_sum + error
ALLOCATE( swup_rts_sicat(t_i_length * t_j_length, nice_use), stat = error )
error_sum = error_sum + error
ALLOCATE( swdn_rts_sicat(t_i_length * t_j_length, nice_use), stat = error )
error_sum = error_sum + error
ALLOCATE( sw_sicat(t_i_length * t_j_length, nice_use), stat = error )
error_sum = error_sum + error
ALLOCATE( alb_sicat(t_i_length * t_j_length, nice_use, 4), stat = error )
error_sum = error_sum + error
ALLOCATE( sw_rts_sea(t_i_length * t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( sw_sea(t_i_length * t_j_length), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  surf_ht_store_surft(:,:) = 0.0
  anthrop_heat_surft(:,:)  = 0.0
  sw_rts_sicat(:,:)        = 0.0
  swup_rts_sicat(:,:)      = 0.0
  swdn_rts_sicat(:,:)      = 0.0
  sw_sicat(:,:)            = 0.0
  alb_sicat(:,:,:)         = 0.0
  sw_rts_sea(:)            = 0.0
  sw_sea(:)                = 0.0
ELSE
  CALL ereport ("Something in fluxes module, common causing error.",          &
                 errcode, "please check allocate_jules_arrays")
END IF

!  ====jules_mod module common====
ALLOCATE( snowdep_surft(land_pts,nsurft), stat = error )
error_sum = error
IF ( error_sum == 0 ) THEN
  snowdep_surft(:,:)          = 0.0
ELSE
  CALL ereport ("Something in jules_mod module common causing error.",        &
                 errcode, "please check allocate_jules_arrays")
END IF

IF ( l_albedo_obs ) THEN
  IF ( l_spec_albedo ) THEN
    ALLOCATE( albobs_scaling_surft(land_pts,ntype,2), stat = error )
  ELSE
    ALLOCATE( albobs_scaling_surft(land_pts,ntype,1), stat = error )
  END IF
  error_sum = error
  IF ( error_sum == 0 ) THEN
    albobs_scaling_surft(:,:,:) = 0.0
  ELSE
    CALL ereport ("Something in jules_mod module common causing error 2.",    &
                 errcode, "please check allocate_jules_arrays")
  END IF
END IF

!========jules_radiation_mod module========
!Nothing to allocate

!========jules_snow_mod module========
IF (ANY(cansnowtile(1:npft)) .EQV. .TRUE.) THEN
  ALLOCATE( unload_backgrnd_pft(land_pts,npft), stat = error )
  error_sum = error
ELSE
  ALLOCATE( unload_backgrnd_pft(1,1), stat = error )
  error_sum = error
END IF
IF ( error_sum == 0 ) THEN
  unload_backgrnd_pft(:,:) = 0.0
ELSE
  CALL ereport ("Something in jules_snow_mod module common causing " //       &
                "error.", errcode, "please check allocate_jules_arrays")
END IF


!  ====jules_surface_mod module common====
ALLOCATE( diff_frac(t_i_length * t_j_length), stat = error )
error_sum = error
IF ( error_sum == 0 ) THEN
  diff_frac(:) = 0.0
ELSE
  CALL ereport ("Something in jules_surfce_mod module common causing " //     &
                "error.", errcode, "please check allocate_jules_arrays")
END IF

!========jules_surface_types_mod module========
!Nothing to allocate

!  ====jules_vegetation_mod module UM====
!Nothing to allocate

!  ====nvegparm module common====
! Non-veg surface type variables.
ALLOCATE( albsnc_nvg(nnvg), stat = error )
error_sum = error
ALLOCATE( albsnf_nvgu(nnvg), stat = error )
error_sum = error_sum + error
ALLOCATE( albsnf_nvg(nnvg), stat = error )
error_sum = error_sum + error
ALLOCATE( albsnf_nvgl(nnvg), stat = error )
error_sum = error_sum + error
ALLOCATE( catch_nvg(nnvg), stat = error )
error_sum = error_sum + error
ALLOCATE( emis_nvg(nnvg), stat = error )
error_sum = error_sum + error
ALLOCATE( gs_nvg(nnvg), stat = error )
error_sum = error_sum + error
ALLOCATE( infil_nvg(nnvg), stat = error )
error_sum = error_sum + error
ALLOCATE( z0_nvg(nnvg), stat = error )
error_sum = error_sum + error
ALLOCATE( ch_nvg(nnvg), stat = error )
error_sum = error_sum + error
ALLOCATE( vf_nvg(nnvg), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  albsnc_nvg(:)  = 0.0
  albsnf_nvgu(:) = 0.0
  albsnf_nvg(:)  = 0.0
  albsnf_nvgl(:) = 0.0
  catch_nvg(:)   = 0.0
  emis_nvg(:)    = 0.0
  gs_nvg(:)      = 0.0
  infil_nvg(:)   = 0.0
  z0_nvg(:)      = 0.0
  ch_nvg(:)      = 0.0
  vf_nvg(:)      = 0.0
ELSE
  CALL ereport ("Something in nvegparm module common causing error.",         &
                 errcode, "please check allocate_jules_arrays")
END IF

!  ====ozone_vars module common====
! Ozone variables
ALLOCATE( o3_gb(land_pts), stat = error )
error_sum = error
ALLOCATE( flux_o3_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fo3_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  o3_gb(:)         = 0.0
  flux_o3_pft(:,:) = 0.0
  fo3_pft(:,:)     = 0.0
ELSE
  CALL ereport ("Something in ozone_vars module common causing error.",       &
                 errcode, "please check allocate_jules_arrays")
END IF

!  ====pftparm module common====
! Veg surface type variables.
ALLOCATE( act_jmax(npft), stat = error )
error_sum = error
ALLOCATE( act_vcmax(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( albsnc_max(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( albsnc_min(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( albsnf_maxu(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( albsnf_max(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( albsnf_maxl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( alpha_elec(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( alpha(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( alniru(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( alnir(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( alnirl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( alparu(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( alpar(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( alparl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( a_wl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( a_ws(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( b_wl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( catch0(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( c3(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( dcatch_dlai(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( deact_jmax(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( deact_vcmax(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( dgl_dm(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( dgl_dt(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( dqcrit(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( ds_jmax(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( ds_vcmax(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( dz0v_dh(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( z0v(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( emis_pft(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( eta_sl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fd(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fsmc_of(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( f0(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( g1_stomata(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( glmin(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( g_leaf_0(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( infil_f(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( jv25_ratio(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( kext(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( kpar(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( lai_alb_lim(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( neff(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( nl0(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( nr_nl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( ns_nl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( nsw(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( nr(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( hw_sw(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( can_struct_a(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( omegau(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( omega(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( omegal(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( omniru(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( omnir(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( omnirl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( orient(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( r_grow(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( rootd_ft(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( psi_close(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( psi_open(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fsmc_p0(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fsmc_mod(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( sigl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( tleaf_of(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( tlow(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( tupp(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( lma(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( nmass(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( vsl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( vint(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( kn(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( knl(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( q10_leaf(npft), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  act_jmax(:)     = 0.0
  act_vcmax(:)    = 0.0
  deact_jmax(:)   = 0.0
  deact_vcmax(:)  = 0.0
  albsnc_max(:)   = 0.0
  albsnc_min(:)   = 0.0
  albsnf_maxu(:)  = 0.0
  albsnf_max(:)   = 0.0
  albsnf_maxl(:)  = 0.0
  alpha_elec(:)   = 0.0
  alpha(:)        = 0.0
  alniru(:)       = 0.0
  alnir(:)        = 0.0
  alnirl(:)       = 0.0
  alparu(:)       = 0.0
  alpar(:)        = 0.0
  alparl(:)       = 0.0
  a_wl(:)         = 0.0
  a_ws(:)         = 0.0
  b_wl(:)         = 0.0
  catch0(:)       = 0.0
  c3(:)           = 0
  dcatch_dlai(:)  = 0.0
  dgl_dm(:)       = 0.0
  dgl_dt(:)       = 0.0
  dqcrit(:)       = 0.0
  ds_jmax(:)      = 0.0
  ds_vcmax(:)     = 0.0
  dz0v_dh(:)      = 0.0
  z0v(:)          = 0.0
  emis_pft(:)     = 0.0
  eta_sl(:)       = 0.0
  fd(:)           = 0.0
  fsmc_of(:)      = 0.0
  f0(:)           = 0.0
  g1_stomata(:)   = 0.0
  glmin(:)        = 0.0
  g_leaf_0(:)     = 0.0
  infil_f(:)      = 0.0
  jv25_ratio(:)   = 0.0
  kext(:)         = 0.0
  kpar(:)         = 0.0
  lai_alb_lim(:)  = 0.0
  neff(:)         = 0.0
  nl0(:)          = 0.0
  nr_nl(:)        = 0.0
  ns_nl(:)        = 0.0
  nsw(:)          = 0.0
  nr(:)           = 0.0
  hw_sw(:)        = 0.0
  can_struct_a(:) = 0.0
  omegau(:)       = 0.0
  omega(:)        = 0.0
  omegal(:)       = 0.0
  omniru(:)       = 0.0
  omnir(:)        = 0.0
  omnirl(:)       = 0.0
  orient(:)       = 0
  r_grow(:)       = 0.0
  rootd_ft(:)     = 0.0
  psi_close(:)    = 0.0
  psi_open(:)     = 0.0
  fsmc_p0(:)      = 0.0
  fsmc_mod(:)     = 0
  sigl(:)         = 0.0
  tleaf_of(:)     = 0.0
  tlow(:)         = 0.0
  tupp(:)         = 0.0
  lma(:)          = 0.0
  nmass(:)        = 0.0
  vsl(:)          = 0.0
  vint(:)         = 0.0
  kn(:)           = 0.0
  knl(:)          = 0.0
  q10_leaf(:)     = 0.0
ELSE
  CALL ereport ("Something in pftparm module common causing error.",          &
                 errcode, "please check allocate_jules_arrays")
END IF

! Ozone damage parameters
ALLOCATE( fl_o3_ct(npft), stat = error )
error_sum = error
ALLOCATE( dfp_dcuo(npft), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  fl_o3_ct(:) = 0.0
  dfp_dcuo(:) = 0.0
ELSE
  CALL ereport ("Something in Ozone damage parameters causing error.",        &
                 errcode, "please check allocate_jules_arrays")
END IF

! BVOC emission parameters
ALLOCATE( ci_st(npft), stat = error )
error_sum = error
ALLOCATE( gpp_st(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( ief(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( tef(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( mef(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( aef(npft), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  ci_st(:)  = 0.0
  gpp_st(:) = 0.0
  ief(:)    = 0.0
  tef(:)    = 0.0
  mef(:)    = 0.0
  aef(:)    = 0.0
ELSE
  CALL ereport ("Something in BVOC emission parameters causing error.",       &
                errcode, "please check allocate_jules_arrays")
END IF

! INFERNO combustion parameters
ALLOCATE( ccleaf_min(npft), stat = error )
error_sum = error
ALLOCATE( ccleaf_max(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( ccwood_min(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( ccwood_max(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( avg_ba(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fire_mort(npft), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  ccleaf_min(:) = 0.0
  ccleaf_max(:) = 0.0
  ccwood_min(:) = 0.0
  ccwood_max(:) = 0.0
  avg_ba(:)     = 0.0
  fire_mort(:)  = 0.0
ELSE
  CALL ereport ("Something in INFERNO combustion parameters causing error.",  &
                 errcode, "please check allocate_jules_arrays")
END IF

! INFERNO emission parameters
ALLOCATE( fef_co2(npft), stat = error )
error_sum = error
ALLOCATE( fef_co(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fef_ch4(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fef_nox(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fef_so2(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fef_oc(npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fef_bc(npft), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  fef_co2(:) = 0.0
  fef_co(:)  = 0.0
  fef_ch4(:) = 0.0
  fef_nox(:) = 0.0
  fef_so2(:) = 0.0
  fef_oc(:)  = 0.0
  fef_bc(:)  = 0.0
ELSE
  CALL ereport ("Something in INFERNO emission parameters causing error.",    &
                 errcode, "please check allocate_jules_arrays")
END IF

!  ====prognostics module common====
ALLOCATE( nsnow_surft(land_pts,nsurft), stat = error )
error_sum = error
ALLOCATE( rho_snow_grnd_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( snowdepth_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( cs_pool_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1), stat = error )
error_sum = error_sum + error
IF (l_layeredc) THEN
  ALLOCATE( t_soil_soilt_acc(land_pts,nsoilt,sm_levels), stat = error )
  error_sum = error_sum + error
END IF
IF ( error_sum == 0 ) THEN
  nsnow_surft(:,:)         = 0
  rho_snow_grnd_surft(:,:) = 0.0
  snowdepth_surft(:,:)     = 0.0
  IF (l_layeredc) THEN
    t_soil_soilt_acc(:,:,:)  = 0.0
  END IF
  cs_pool_soilt(:,:,:,:)     = 0.0
ELSE
  CALL ereport ("Something in prognostics module common causing error.",      &
                 errcode, "please check allocate_jules_arrays")
END IF

! Nitrogen scheme variables
ALLOCATE( ns_pool_gb(land_pts,dim_cslayer,dim_cs1), stat = error )
error_sum = error
ALLOCATE( n_inorg_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( n_inorg_soilt_lyrs(land_pts,nsoilt,dim_cslayer), stat = error )
error_sum = error_sum + error
ALLOCATE( n_inorg_avail_pft(land_pts,npft,dim_cslayer), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  ns_pool_gb(:,:,:)         = 0.001
  n_inorg_soilt_lyrs(:,:,:) = 0.0
  n_inorg_gb(:)             = 0.0
  n_inorg_avail_pft(:,:,:)  = 0.0
ELSE
  CALL ereport ("Something in Nitrogen scheme variables causing error.",      &
                 errcode, "please check allocate_jules_arrays")
END IF

!For multilayer snow variables, only allocate to full size if the scheme is
!being used, ie nsmax > 0
IF (nsmax > 0) THEN
  temp_size   = land_pts
  temp_tiles  = nsurft
  temp_layers = nsmax
ELSE
  temp_size   = 1
  temp_tiles  = 1
  temp_layers = 1
END IF

ALLOCATE( ds_surft(temp_size,temp_tiles,temp_layers),      stat = error )
error_sum = error
ALLOCATE( rgrainl_surft(temp_size,temp_tiles,temp_layers), stat = error )
error_sum = error_sum + error
ALLOCATE( sice_surft(temp_size,temp_tiles,temp_layers),    stat = error )
error_sum = error_sum + error
ALLOCATE( sliq_surft(temp_size,temp_tiles,temp_layers),    stat = error )
error_sum = error_sum + error
ALLOCATE( tsnow_surft(temp_size,temp_tiles,temp_layers),   stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  ds_surft(:,:,:)      = 0.0
  rgrainl_surft(:,:,:) = 0.0
  sice_surft(:,:,:)    = 0.0
  sliq_surft(:,:,:)    = 0.0
  tsnow_surft(:,:,:)   = 0.0
ELSE
  CALL ereport ("Something in multilayer snow variables causing error.",      &
                 errcode, "please check allocate_jules_arrays")
END IF

!See comment in prognostics module
ALLOCATE( rho_snow_surft(land_pts,nsurft,nsmax), stat = error )
error_sum = error
IF ( error_sum == 0 ) THEN
  rho_snow_surft(:,:,:) = 0.0
ELSE
  CALL ereport ("Something in See comment in prognostics module causing " //  &
                "error.", errcode, "please check allocate_jules_arrays")
END IF

! Allocate WP Pools
IF ( l_triffid .OR. l_phenol ) THEN
  ALLOCATE( wood_prod_fast_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( wood_prod_med_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( wood_prod_slow_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( frac_agr_prev_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( frac_past_prev_gb(land_pts), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    wood_prod_fast_gb(:) = 0.0
    wood_prod_med_gb(:)  = 0.0
    wood_prod_slow_gb(:) = 0.0
    frac_agr_prev_gb(:)  = 0.0
    frac_past_prev_gb(:) = 0.0
  ELSE
    CALL ereport ("Something in Allocate WP Pools causing error.",            &
                 errcode, "please check allocate_jules_arrays")
  END IF
END IF

! If TRIFFID is switched on, the TRIFFID-derived CO2 flux variable, 
! which is passed to the atmosphere in UM simulations with interactive CO2,
! needs to be allocated. Not testing for L_CO2_INTERACTIVE so that the
! diagnostic of triffid_co2_gb can be output with or without interactive CO2.   
IF ( l_triffid ) THEN
  ALLOCATE( triffid_co2_gb(land_pts), stat = error )
  error_sum = error

  IF ( error_sum == 0 ) THEN
    triffid_co2_gb(:)    = 0.0
  ELSE
    CALL ereport ("Something in allocate triffid triffid_co2_gb  " //         &
                  "causing error.", errcode, "please check "  //              &
                  "allocate_jules_arrays")
  END IF
END IF

! Only allocate the bedrock tsoil_deep_gb if bedrock is being used
IF ( l_bedrock ) THEN
  ALLOCATE( tsoil_deep_gb(land_pts,ns_deep), stat = error )
  error_sum = error
  IF ( error_sum == 0 ) THEN
    tsoil_deep_gb(:,:) = 0.0
  ELSE
    CALL ereport ("Something in allocate the bedrock tsoil_deep_gb  " //      &
                  "causing error.", errcode, "please check " //               &
                  "allocate_jules_arrays")
  END IF
END IF

!  ====p_s_parms module common====
ALLOCATE(bexp_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error
ALLOCATE(sathh_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE(hcap_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE(hcon_soilt(land_pts,nsoilt,0:sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE(satcon_soilt(land_pts,nsoilt,0:sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE(smvccl_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE(smvcwt_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE(smvcst_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE( clay_soilt(land_pts,nsoilt,dim_cslayer), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  bexp_soilt(:,:,:)   = 0.0
  sathh_soilt(:,:,:)  = 0.0
  hcap_soilt(:,:,:)   = 0.0
  hcon_soilt(:,:,:)   = 0.0
  satcon_soilt(:,:,:) = 0.0
  smvccl_soilt(:,:,:) = 0.0
  smvcwt_soilt(:,:,:) = 0.0
  smvcst_soilt(:,:,:) = 0.0
  clay_soilt(:,:,:)     = 0.0
ELSE
  CALL ereport ("Something in p_s_parms module common causing error.",        &
                 errcode, "please check allocate_jules_arrays")
END IF

!  ====switches_urban module common====
!Nothing to allocate

!========theta_field_sizes module========
!Nothing to allocate

!  ====top_pdm module JULES====
!       TOPMODEL and PDM variables
ALLOCATE( a_fsat_soilt(land_pts,nsoilt), stat = error )
error_sum = error
ALLOCATE( a_fwet_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( c_fsat_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( c_fwet_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( drain_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( dun_roff_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( fch4_wetl_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( fch4_wetl_acc_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( fch4_wetl_cs_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( fch4_wetl_npp_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( fch4_wetl_resps_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( fexp_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( fsat_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( fwetl_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( gamtot_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( qbase_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( qbase_zw_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( sthzw_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( ti_mean_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( ti_sig_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( zw_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( inlandout_atm_gb(land_pts), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  a_fsat_soilt(:,:)          = 0.0
  a_fwet_soilt(:,:)          = 0.0
  c_fsat_soilt(:,:)          = 0.0
  c_fwet_soilt(:,:)          = 0.0
  drain_soilt(:,:)           = 0.0
  dun_roff_soilt(:,:)        = 0.0
  fch4_wetl_soilt(:,:)       = 0.0
  fch4_wetl_acc_soilt(:,:)   = 0.0
  fch4_wetl_cs_soilt(:,:)    = 0.0
  fch4_wetl_resps_soilt(:,:) = 0.0
  fch4_wetl_npp_soilt(:,:)   = 0.0
  fexp_soilt(:,:)            = 0.0
  fsat_soilt(:,:)            = 0.0
  fwetl_soilt(:,:)           = 0.0
  gamtot_soilt(:,:)          = 0.0
  qbase_soilt(:,:)           = 0.0
  qbase_zw_soilt(:,:)        = 0.0
  sthzw_soilt(:,:)           = 0.0
  ti_mean_soilt(:,:)         = 0.0
  ti_sig_soilt(:,:)          = 0.0
  zw_soilt(:,:)              = 0.0
  inlandout_atm_gb(:)        = 0.0
ELSE
  CALL ereport ("Something in top_pdm module JULES " //                       &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF

!  ====trif module common====
! TRIFFID variables - only needed if TRIFFID and/or phenology is selected.
IF ( l_triffid .OR. l_phenol ) THEN
  ALLOCATE( crop(npft), stat = error )
  error_sum = error
  ALLOCATE( g_area(npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( g_grow(npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( g_root(npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( g_wood(npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( lai_max(npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( lai_min(npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( alloc_fast(npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( alloc_med(npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( alloc_slow(npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( dpm_rpm_ratio(npft), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    crop(:)          = 0
    g_area(:)        = 0.0
    g_grow(:)        = 0.0
    g_root(:)        = 0.0
    g_wood(:)        = 0.0
    lai_max(:)       = 0.0
    lai_min(:)       = 0.0
    alloc_fast(:)    = 0.0
    alloc_med(:)     = 0.0
    alloc_slow(:)    = 0.0
    dpm_rpm_ratio(:) = 0.0
  ELSE
    CALL ereport ("Something in trif module common causing error", errcode,   &
                  "please check allocate_jules_arrays")
  END IF
END IF

!  ====trifctl module common====
! Triffid variables
ALLOCATE( g_leaf_acc_pft(land_pts,npft), stat = error )
error_sum = error
ALLOCATE( npp_acc_pft(land_pts_TRIF,npft_trif), stat = error )
error_sum = error_sum + error
ALLOCATE( resp_w_acc_pft(land_pts_TRIF,npft_trif), stat = error )
error_sum = error_sum + error
ALLOCATE( resp_s_acc_soilt(land_pts_TRIF,nsoilt,dim_cslayer,dim_cs1),         &
                           stat = error )
error_sum = error_sum + error
ALLOCATE( g_leaf_phen_acc_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( gpp_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( npp_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( resp_p_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( g_leaf_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( g_leaf_phen_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( gpp_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( npp_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( resp_p_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( resp_s_soilt(land_pts,nsoilt,dim_cslayer,dim_cs1), stat = error )
error_sum = error_sum + error
ALLOCATE( resp_w_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( lai_phen_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( c_veg_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( cv_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( g_leaf_day_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( g_leaf_dr_out_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( lit_c_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( lit_c_mn_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( npp_dr_out_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( resp_w_dr_out_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( resp_s_dr_out_gb(land_pts,dim_cslayer,5), stat = error )
error_sum = error_sum + error
ALLOCATE( frac_agr_gb(land_pts), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  g_leaf_acc_pft(:,:)       = 0.0
  npp_acc_pft(:,:)          = 0.0
  resp_w_acc_pft(:,:)       = 0.0
  resp_s_acc_soilt(:,:,:,:) = 0.0
  g_leaf_phen_acc_pft(:,:)  = 0.0
  gpp_gb(:)                 = 0.0
  npp_gb(:)                 = 0.0
  resp_p_gb(:)              = 0.0
  g_leaf_pft(:,:)           = 0.0
  g_leaf_phen_pft(:,:)      = 0.0
  gpp_pft(:,:)              = 0.0
  npp_pft(:,:)              = 0.0
  resp_p_pft(:,:)           = 0.0
  resp_s_soilt(:,:,:,:)     = 0.0
  resp_w_pft(:,:)           = 0.0
  lai_phen_pft(:,:)         = 0.0
  c_veg_pft(:,:)            = 0.0
  cv_gb(:)                  = 0.0
  g_leaf_day_pft(:,:)       = 0.0
  g_leaf_dr_out_pft(:,:)    = 0.0
  lit_c_pft(:,:)            = 0.0
  lit_c_mn_gb(:)            = 0.0
  npp_dr_out_pft(:,:)       = 0.0
  resp_w_dr_out_pft(:,:)    = 0.0
  resp_s_dr_out_gb(:,:,:)   = 0.0
  frac_agr_gb(:)            = 0.0
ELSE
  CALL ereport ("Something in trifctl module JULES " //                       &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF

!  ====trif_vars_mod module common====
ALLOCATE( resp_l_pft(land_pts,npft), stat = error )
error_sum = error
ALLOCATE( resp_r_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( n_leaf_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( n_root_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( n_stem_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( lai_bal_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( pc_s_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( fapar_diag_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( apar_diag_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( apar_diag_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( fao_et0(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( frac_past_gb(land_pts), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  resp_l_pft(:,:)     = 0.0
  resp_r_pft(:,:)     = 0.0
  n_leaf_pft(:,:)     = 0.0
  n_root_pft(:,:)     = 0.0
  n_stem_pft(:,:)     = 0.0
  lai_bal_pft(:,:)    = 0.0
  pc_s_pft(:,:)       = 0.0
  fapar_diag_pft(:,:) = 0.0
  apar_diag_pft(:,:)  = 0.0
  apar_diag_gb(:)     = 0.0
  fao_et0(:)          = 0.0
  frac_past_gb(:)     = 0.0
ELSE
  CALL ereport ("Something in trif_vars_mod common 1 causing error.",         &
                 errcode, "please check allocate_jules_arrays")
END IF

IF ( l_triffid .OR. l_phenol ) THEN
  ALLOCATE( wp_fast_in_gb(land_pts), stat = error )
  error_sum = error
  ALLOCATE( wp_med_in_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( wp_slow_in_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( wp_fast_out_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( wp_med_out_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( wp_slow_out_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( lit_c_orig_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( lit_c_ag_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( lit_c_fire_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( lit_c_nofire_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( lit_n_fire_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( lit_n_nofire_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( veg_c_fire_emission_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( veg_c_fire_emission_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error

  IF ( error_sum == 0 ) THEN
    wp_fast_in_gb(:)             = 0.0
    wp_med_in_gb(:)              = 0.0
    wp_slow_in_gb(:)             = 0.0
    wp_fast_out_gb(:)            = 0.0
    wp_med_out_gb(:)             = 0.0
    wp_slow_out_gb(:)            = 0.0
    lit_c_orig_pft(:,:)          = 0.0
    lit_c_ag_pft(:,:)            = 0.0
    lit_c_fire_pft(:,:)          = 0.0
    lit_c_nofire_pft(:,:)        = 0.0
    lit_n_fire_pft(:,:)          = 0.0
    lit_n_nofire_pft(:,:)        = 0.0
    veg_c_fire_emission_gb(:)    = 0.0
    veg_c_fire_emission_pft(:,:) = 0.0
  ELSE
    CALL ereport ("Something in trif_vars_mod common 2 causing error.",       &
                   errcode, "please check allocate_jules_arrays")
  END IF

  ALLOCATE( cnsrv_carbon_veg2_gb(land_pts), stat = error )
  error_sum = error
  ALLOCATE( cnsrv_carbon_triffid_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( cnsrv_veg_triffid_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( cnsrv_soil_triffid_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( cnsrv_prod_triffid_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( resp_s_to_atmos_gb(land_pts,dim_cslayer), stat = error )
  error_sum = error_sum + error
  ALLOCATE( root_abandon_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( root_abandon_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( cnsrv_nitrogen_triffid_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( cnsrv_vegN_triffid_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( cnsrv_soilN_triffid_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( cnsrv_N_inorg_triffid_gb(land_pts), stat = error )
  error_sum = error_sum + error    
  ALLOCATE( harvest_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( harvest_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( root_abandon_n_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( root_abandon_n_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( harvest_n_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( harvest_n_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_fertiliser_add(land_pts,nsoilt,dim_cslayer), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_fertiliser_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_fertiliser_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( g_burn_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_leaf_trif_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_leaf_alloc_trif_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_leaf_labile_trif_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_stem_trif_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_root_trif_pft(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( lit_n_ag_pft_diag(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( lit_n_pft_diag(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_luc(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( g_burn_pft_acc(land_pts,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( g_burn_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( burnt_carbon_dpm(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( burnt_carbon_rpm(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE ( gpp_gb_out(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( gpp_pft_out(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( gpp_gb_acc(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( gpp_pft_acc(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( resp_p_actual_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( resp_p_actual_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_leach_gb_acc(land_pts), stat = error)
  error_sum = error_sum + error

  IF ( error_sum == 0 ) THEN
    cnsrv_carbon_veg2_gb(:)      = 0.0
    cnsrv_carbon_triffid_gb(:)   = 0.0
    cnsrv_veg_triffid_gb(:)      = 0.0
    cnsrv_soil_triffid_gb(:)     = 0.0
    cnsrv_prod_triffid_gb(:)     = 0.0
    resp_s_to_atmos_gb(:,:)      = 0.0
    root_abandon_pft(:,:)        = 0.0
    root_abandon_gb(:)           = 0.0
    cnsrv_nitrogen_triffid_gb(:) = 0.0
    cnsrv_vegN_triffid_gb(:)     = 0.0
    cnsrv_soilN_triffid_gb(:)    = 0.0
    cnsrv_N_inorg_triffid_gb(:)  = 0.0
    harvest_pft(:,:)             = 0.0
    harvest_gb(:)                = 0.0
    root_abandon_n_pft(:,:)      = 0.0
    root_abandon_n_gb(:)         = 0.0
    harvest_n_pft(:,:)           = 0.0
    harvest_n_gb(:)              = 0.0
    n_fertiliser_add(:,:,:)      = 0.0
    n_fertiliser_pft(:,:)        = 0.0
    n_fertiliser_gb(:)           = 0.0
    g_burn_pft(:,:)              = 0.0
    n_leaf_trif_pft(:,:)         = 0.0
    n_leaf_alloc_trif_pft(:,:)   = 0.0
    n_leaf_labile_trif_pft(:,:)  = 0.0
    n_stem_trif_pft(:,:)         = 0.0
    n_root_trif_pft(:,:)         = 0.0
    lit_n_ag_pft_diag(:,:)       = 0.0
    lit_n_pft_diag(:,:)          = 0.0
    n_luc(:)                     = 0.0
    g_burn_pft_acc(:,:)          = 0.0   
    g_burn_gb(:)                 = 0.0
    burnt_carbon_dpm             = 0.0
    burnt_carbon_rpm             = 0.0      
    gpp_gb_out(:)                = 0.0
    gpp_pft_out(:,:)             = 0.0
    gpp_gb_acc(:)                = 0.0
    gpp_pft_acc(:,:)             = 0.0
    resp_p_actual_gb(:)          = 0.0
    resp_p_actual_pft(:,:)       = 0.0
    n_leach_gb_acc(:)            = 0.0
  ELSE
    CALL ereport ("Something in trif_vars_mod common 3 causing error.",       &
                   errcode, "please check allocate_jules_arrays")
  END IF

  ! Variables added for #7 (nitrogen scheme)
  ALLOCATE( resp_s_diag_gb(land_pts,dim_cslayer,dim_cs1+1), stat = error )
  error_sum = error
  ALLOCATE( resp_s_pot_diag_gb(land_pts,dim_cslayer,dim_cs1+1), stat = error )
  error_sum = error_sum + error
  ALLOCATE ( n_veg_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_veg_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( dnveg_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( dcveg_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( dnveg_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( dcveg_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_demand_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_uptake_extract(land_pts,nsoilt,dim_cslayer), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_uptake_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( deposition_N_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_uptake_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_demand_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( dleaf_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( droot_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( dwood_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_demand_lit_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_demand_spread_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_demand_growth_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_uptake_growth_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_uptake_spread_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_fix_add(land_pts,nsoilt,dim_cslayer), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_fix_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_fix_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_leach_soilt(land_pts,nsoilt), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_gas_gb(land_pts,dim_cslayer), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( exudates_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( exudates_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( npp_n_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( npp_n(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( root_litC_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( leaf_litC_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( wood_litC_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( root_litN_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( leaf_litN_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( wood_litN_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( lit_N_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( lit_n_t_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE( minl_n_pot_gb(land_pts,dim_cslayer,dim_cs1+1), stat = error )
  error_sum = error_sum + error
  ALLOCATE( immob_n_gb(land_pts,dim_cslayer,dim_cs1+1), stat = error )
  error_sum = error_sum + error
  ALLOCATE( immob_n_pot_gb(land_pts,dim_cslayer,dim_cs1+1), stat = error )
  error_sum = error_sum + error
  ALLOCATE( fn_gb(land_pts,dim_cslayer), stat = error )
  error_sum = error_sum + error
  ALLOCATE( minl_n_gb(land_pts,dim_cslayer,dim_cs1+1), stat = error )
  error_sum = error_sum + error
  ALLOCATE ( dpm_ratio_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( n_loss_gb(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( leafC_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( rootC_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( woodC_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( leafC_gbm(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( rootC_gbm(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( woodC_gbm(land_pts), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( litterC_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( stemC_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( litterN_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( retran_r(npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( retran_l(npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( lit_n_orig_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  ALLOCATE ( lit_n_ag_pft(land_pts,npft), stat = error)
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    resp_s_diag_gb(:,:,:)        = 0.0
    resp_s_pot_diag_gb(:,:,:)    = 0.0
    n_veg_pft(:,:)               = 0.0
    n_veg_gb(:)                  = 0.0
    dnveg_pft(:,:)               = 0.0
    dcveg_pft(:,:)               = 0.0
    dnveg_gb(:)                  = 0.0
    dcveg_gb(:)                  = 0.0
    n_demand_gb(:)               = 0.0
    n_uptake_extract(:,:,:)      = 0.0
    n_uptake_gb(:)               = 0.0
    deposition_n_gb(:)           = 0.0
    n_uptake_pft(:,:)            = 0.0
    n_demand_pft(:,:)            = 0.0
    dleaf_pft(:,:)               = 0.0
    droot_pft(:,:)               = 0.0
    dwood_pft(:,:)               = 0.0
    n_demand_lit_pft(:,:)        = 0.0
    n_demand_spread_pft(:,:)     = 0.0
    n_demand_growth_pft(:,:)     = 0.0
    n_uptake_growth_pft(:,:)     = 0.0
    n_uptake_spread_pft(:,:)     = 0.0
    n_fix_add(:,:,:)             = 0.0
    n_fix_gb(:)                  = 0.0
    n_fix_pft(:,:)               = 0.0
    n_leach_soilt(:,:)           = 0.0
    n_gas_gb(:,:)                = 0.0
    exudates_pft(:,:)            = 0.0
    exudates_gb(:)               = 0.0
    npp_n_gb(:)                  = 0.0
    npp_n(:,:)                   = 0.0
    root_litC_pft(:,:)           = 0.0
    leaf_litC_pft(:,:)           = 0.0
    wood_litC_pft(:,:)           = 0.0
    root_litN_pft(:,:)           = 0.0
    leaf_litN_pft(:,:)           = 0.0
    wood_litN_pft(:,:)           = 0.0
    lit_N_pft(:,:)               = 0.0
    lit_n_t_gb(:)                = 0.0
    minl_n_pot_gb(:,:,:)         = 0.0
    immob_n_gb(:,:,:)            = 0.0
    immob_n_pot_gb(:,:,:)        = 0.0
    fn_gb(:,:)                   = 0.0
    minl_n_gb(:,:,:)             = 0.0
    dpm_ratio_gb(:)              = 0.0
    n_loss_gb(:)                 = 0.0
    leafC_pft(:,:)               = 0.0
    rootC_pft(:,:)               = 0.0
    woodC_pft(:,:)               = 0.0
    leafC_gbm(:)                 = 0.0
    rootC_gbm(:)                 = 0.0
    woodC_gbm(:)                 = 0.0
    litterC_pft(:,:)             = 0.0
    stemC_pft(:,:)               = 0.0
    litterN_pft(:,:)             = 0.0
    retran_r(:)                  = 0.0
    retran_l(:)                  = 0.0
    lit_n_orig_pft(:,:)          = 0.0
    lit_n_ag_pft(:,:)            = 0.0

  ELSE
    CALL ereport ("Something in Variables added for #7 (nitrogen scheme) "//  &
                 "causing error.", errcode, "please check " //                &
                 "allocate_jules_arrays")
  END IF
END IF

!  ====urban_param module UM====
IF ( l_urban2t .OR. l_moruses ) THEN
  CALL jules_print(                                                           &
      'allocate_jules_arrays',                                                &
      'Allocating URBAN-2T / MORUSES arrays',                                 &
      level = PrNorm)
  ALLOCATE( wrr_gb(land_pts), stat = error )
  error_sum = error
  ALLOCATE( hgt_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( hwr_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( disp_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( ztm_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( albwl_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( albrd_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( emisw_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( emisr_gb(land_pts), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    wrr_gb(:)   = 0.0
    hgt_gb(:)   = 0.0
    hwr_gb(:)   = 0.0
    albwl_gb(:) = 0.0
    albrd_gb(:) = 0.0
    emisw_gb(:) = 0.0
    emisr_gb(:) = 0.0
    ztm_gb(:)   = 0.0
    disp_gb(:)  = 0.0
  ELSE
    CALL ereport ("Something in urban_param module UM causing error.",        &
                   errcode, "please check allocate_jules_arrays")
  END IF
END IF

!Model-dependent USE statements

#if defined(UM_JULES)
!#############################################################################
! Section for UM_JULES.
!#############################################################################
!  ====atm_fields_bounds_mod module UM
!Nothing to allocate

!  ====jules_surface_mod module UM====
!Nothing to allocate

!  ====lake_mod module UM====
IF ( l_flake_model ) THEN
  ALLOCATE( surf_ht_flux_lake_ij(tdims%i_start:tdims%i_end                    &
                             ,tdims%j_start:tdims%j_end), stat = error )
  error_sum = error
  temp_size = land_pts
ELSE
  ALLOCATE( surf_ht_flux_lake_ij(1,1), stat = error )
  error_sum = error
  temp_size = 1
END IF
IF ( error_sum == 0 ) THEN
  surf_ht_flux_lake_ij(:,:) = 0.0
ELSE
  CALL ereport ("Something in lake_mod module UM 1 causing error.",           &
               errcode, "please check allocate_jules_arrays")
END IF

ALLOCATE( surf_ht_flux_lk_gb(temp_size), stat = error )
error_sum = error
ALLOCATE( sw_down_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( coriolis_param_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( u_s_lake_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_depth_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_fetch_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_albedo_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_t_snow_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_t_ice_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_t_mean_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_t_mxl_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_shape_factor_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_h_snow_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_h_ice_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_h_mxl_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( lake_t_sfc_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( ts1_lake_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( nusselt_gb(temp_size), stat = error )
error_sum = error_sum + error
ALLOCATE( g_dt_gb(temp_size), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  surf_ht_flux_lk_gb(:)   = 0.0
  sw_down_gb(:)           = 0.0
  coriolis_param_gb(:)    = 0.0
  u_s_lake_gb(:)          = 0.0
  lake_depth_gb(:)        = 0.0
  lake_fetch_gb(:)        = 0.0
  lake_albedo_gb(:)       = 0.0
  lake_t_snow_gb(:)       = 0.0
  lake_t_ice_gb(:)        = 0.0
  lake_t_mean_gb(:)       = 0.0
  lake_t_mxl_gb(:)        = 0.0
  lake_shape_factor_gb(:) = 0.0
  lake_h_snow_gb(:)       = 0.0
  lake_h_ice_gb(:)        = 0.0
  lake_h_mxl_gb(:)        = 0.0
  lake_t_sfc_gb(:)        = 0.0
  ts1_lake_gb(:)          = 0.0
  nusselt_gb(:)           = 0.0
  g_dt_gb(:)              = 0.0
ELSE
  CALL ereport ("Something in lake_mod module UM 3 causing error.",           &
               errcode, "please check allocate_jules_arrays")
END IF

!  ====pftparm module UM====
ALLOCATE( dust_veg_scj(npft), stat = error )
error_sum = error
IF ( error_sum == 0 ) THEN
  dust_veg_scj(:) = 0.0
ELSE
  CALL ereport ("Something in pftparm module UM causing error.",              &
               errcode, "please check allocate_jules_arrays")
END IF
#else
!#############################################################################
! Section for standalone JULES.
!#############################################################################
!  ====aero module JULES====
ALLOCATE( co2_3d_ij(t_i_length,t_j_length), stat = error )
error_sum = error
ALLOCATE( rho_cd_modv1_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( rho_aresist_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( aresist_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( resist_b_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( rho_aresist_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( aresist_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( resist_b_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( r_b_dust_ij(t_i_length,t_j_length,ndiv), stat = error )
error_sum = error_sum + error
ALLOCATE( cd_std_dust_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( u_s_std_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  co2_3d_ij(:,:)         = 0.0
  rho_cd_modv1_ij(:,:)   = 0.0
  rho_aresist_ij(:,:)    = 0.0
  aresist_ij(:,:)        = 0.0
  resist_b_ij(:,:)       = 0.0
  rho_aresist_surft(:,:) = 0.0
  aresist_surft(:,:)     = 0.0
  resist_b_surft(:,:)    = 0.0
  r_b_dust_ij(:,:,:)     = 0.0
  cd_std_dust_ij(:,:)    = 0.0
  u_s_std_surft(:,:)     = 0.0
ELSE
  CALL ereport ("Something in aero module JULES causing error.",              &
               errcode, "please check allocate_jules_arrays")
END IF

!  ====ancil_info module JULES====
ALLOCATE( surft_index(land_pts,ntype), stat = error )
error_sum = error
ALLOCATE( soil_index(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( lice_index(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( ice_fract_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( ice_fract_ncat_sicat(t_i_length,t_j_length,nice), stat = error )
error_sum = error_sum + error
ALLOCATE( ti_cat_sicat(t_i_length,t_j_length,nice), stat = error )
error_sum = error_sum + error
ALLOCATE( pond_frac_cat_sicat(t_i_length,t_j_length,nice), stat = error )
error_sum = error_sum + error
ALLOCATE( pond_depth_cat_sicat(t_i_length,t_j_length,nice), stat = error )
error_sum = error_sum + error
ALLOCATE( sstfrz_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( z1_uv_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( z1_tq_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( frac_surft(land_pts,ntype), stat = error )
error_sum = error_sum + error
ALLOCATE( surft_pts(ntype), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  surft_index(:,:)            = 0
  soil_index(:)               = 0
  lice_index(:)               = 0
  ice_fract_ij(:,:)           = 0.0
  ice_fract_ncat_sicat(:,:,:) = 0.0
  ti_cat_sicat(:,:,:)         = 0.0
  pond_frac_cat_sicat(:,:,:)  = 0.0
  pond_depth_cat_sicat(:,:,:) = 0.0
  sstfrz_ij(:,:)              = 0.0
  z1_uv_ij(:,:)               = 0.0
  z1_tq_ij(:,:)               = 0.0
  frac_surft(:,:)             = 0.0
  surft_pts(:)                = 0
ELSE
  CALL ereport ("Something in ancil_info module JULES causing error.",        &
               errcode, "please check allocate_jules_arrays")
END IF

!  ====boundary_layer_vars_mod module JULES====
ALLOCATE( zh(t_i_length,t_j_length), stat = error )
error_sum = error
IF ( l_deposition ) THEN
  ALLOCATE( dzl(t_i_length,t_j_length,bl_levels), stat = error )
  error_sum = error_sum + error
END IF
IF ( error_sum == 0 ) THEN
  zh(:,:) = 0.0
  IF ( l_deposition ) dzl(:,:,:) = 0.0
ELSE
  CALL ereport ("Something in boundary_layer_vars_mod module JULES " //       &
                "causing error.",                                             &
               errcode, "please check allocate_jules_arrays")
END IF

!  ====c_elevate module JULES====
! Land height
ALLOCATE( z_land_ij(t_i_length,t_j_length), stat = error )
error_sum = error
IF ( error_sum == 0 ) THEN
  z_land_ij(:,:) = 0.0
ELSE
  CALL ereport ("Something in c_elevate module JULES causing error.",         &
               errcode, "please check allocate_jules_arrays")
END IF

!  ====coastal module JULES====
! Coastal tiling variables
ALLOCATE( tstar_land_ij(t_i_length,t_j_length), stat = error )
error_sum = error
ALLOCATE( tstar_sea_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( tstar_sice_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( tstar_sice_sicat(t_i_length,t_j_length,nice_use), stat = error )
error_sum = error_sum + error
ALLOCATE( tstar_ssi_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( taux_land_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( taux_ssi_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( tauy_land_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( tauy_ssi_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( vshr_land_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( vshr_ssi_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( surf_ht_flux_land_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( surf_ht_flux_sice_sicat(t_i_length,t_j_length,nice), stat = error )
error_sum = error_sum + error
ALLOCATE( taux_land_star(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( tauy_land_star(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( taux_ssi_star(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( tauy_ssi_star(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  tstar_land_ij(:,:)             = 0.0
  tstar_sea_ij(:,:)              = 0.0
  tstar_sice_ij(:,:)             = 0.0
  tstar_sice_sicat(:,:,:)        = 0.0
  tstar_ssi_ij(:,:)              = 0.0
  taux_land_ij(:,:)              = 0.0
  taux_ssi_ij(:,:)               = 0.0
  tauy_land_ij(:,:)              = 0.0
  tauy_ssi_ij(:,:)               = 0.0
  vshr_land_ij(:,:)              = 0.0
  vshr_ssi_ij(:,:)               = 0.0
  surf_ht_flux_land_ij(:,:)      = 0.0
  surf_ht_flux_sice_sicat(:,:,:) = 0.0
  taux_land_star(:,:)            = 0.0
  tauy_land_star(:,:)            = 0.0
  taux_ssi_star(:,:)             = 0.0
  tauy_ssi_star(:,:)             = 0.0
ELSE
  CALL ereport ("Something in coastal module JULES causing error.",           &
               errcode, "please check allocate_jules_arrays")
END IF

!  ====cropparm module JULES====
!  ====crop_vars_mod module JULES====
!  Done together as they share the l_crop IF
IF ( l_crop ) THEN
  ! Crop variables
  ALLOCATE( t_bse(ncpft), stat = error )
  error_sum = error
  ALLOCATE( t_opt(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( t_max(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( tt_emr(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( crit_pp(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( pp_sens(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( rt_dir(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( alpha1(ncpft), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    t_bse(:)   = 0.0
    t_opt(:)   = 0.0
    t_max(:)   = 0.0
    tt_emr(:)  = 0.0
    crit_pp(:) = 0.0
    pp_sens(:) = 0.0
    rt_dir(:)  = 0.0
    alpha1(:)  = 0.0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_crop IF 1 causing error.", errcode, "please check " //   &
                  "allocate_jules_arrays")
  END IF

  ALLOCATE( alpha2(ncpft), stat = error )
  error_sum = error
  ALLOCATE( alpha3(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( beta1(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( beta2(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( beta3(ncpft), stat = error )
  error_sum = error
  IF ( error_sum == 0 ) THEN
    alpha2(:) = 0.0
    alpha3(:) = 0.0
    beta1(:)  = 0.0
    beta2(:)  = 0.0
    beta3(:)  = 0.0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_crop IF 2 causing error.", errcode, "please check " //   &
                  "allocate_jules_arrays")
  END IF

  ALLOCATE( r_gamma(ncpft), stat = error )
  error_sum = error
  ALLOCATE( delta(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( remob(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( cfrac_s(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( cfrac_r(ncpft), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    r_gamma(:) = 0.0
    delta(:)   = 0.0
    remob(:)   = 0.0
    cfrac_s(:) = 0.0
    cfrac_r(:) = 0.0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_crop IF 3 causing error.", errcode, "please check " //   &
                  "allocate_jules_arrays")
  END IF

  ALLOCATE( cfrac_l(ncpft), stat = error )
  error_sum = error
  ALLOCATE( allo1(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( allo2(ncpft), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    cfrac_l(:) = 0.0
    allo1(:)   = 0.0
    allo2(:)   = 0.0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_crop IF 4 causing error.", errcode, "please check " //   &
                  "allocate_jules_arrays")
  END IF

  ALLOCATE( mu(ncpft), stat = error )
  error_sum = error
  ALLOCATE( nu(ncpft), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    mu(:) = 0.0
    nu(:) = 0.0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_crop IF 5 causing error.", errcode, "please check " //   &
                  "allocate_jules_arrays")
  END IF

  ALLOCATE( yield_frac(ncpft), stat = error )
  error_sum = error
  ALLOCATE( initial_carbon(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( initial_c_dvi(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( sen_dvi(ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( t_mort(ncpft), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    yield_frac(:)     = 0.0
    initial_carbon(:) = 0.0
    initial_c_dvi(:)  = 0.0
    sen_dvi(:)        = 0.0
    t_mort(:)         = 0.0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_crop IF 6 causing error.", errcode, "please check " //   &
                  "allocate_jules_arrays")
  END IF

  ALLOCATE( dvi_cpft(land_pts,ncpft), stat = error )
  error_sum = error
  ALLOCATE( croprootc_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( harvc_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( reservec_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    dvi_cpft(:,:)       = 0.0
    croprootc_cpft(:,:) = 0.0
    harvc_cpft(:,:)     = 0.0
    reservec_cpft(:,:)  = 0.0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_crop IF 7 causing error.", errcode, "please check " //   &
                  "allocate_jules_arrays")
  END IF

  ALLOCATE( yield_diag_cpft(land_pts,ncpft), stat = error )
  error_sum = error
  ALLOCATE( nonyield_diag_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( harvest_trigger_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( harvest_counter_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( leafc_diag_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( stemc_diag_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    yield_diag_cpft(:,:)      = 0.0
    nonyield_diag_cpft(:,:)   = 0.0
    harvest_trigger_cpft(:,:) = 0
    harvest_counter_cpft(:,:) = 0
    stemc_diag_cpft(:,:)      = 0.0
    leafc_diag_cpft(:,:)      = 0.0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_crop IF 8 causing error.", errcode, "please check " //   &
                  "allocate_jules_arrays")
  END IF

  ALLOCATE( croplai_cpft(land_pts,ncpft), stat = error )
  error_sum = error
  ALLOCATE( cropcanht_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    croplai_cpft(:,:)   = 0.0
    cropcanht_cpft(:,:) = 0.0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_crop IF 9 causing error.", errcode, "please check " //   &
                  "allocate_jules_arrays")
  END IF

  ALLOCATE( sow_date_cpft(land_pts,ncpft), stat = error )
  error_sum = error
  ALLOCATE( tt_veg_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( tt_rep_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( latestharv_date_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error

  IF ( error_sum == 0 ) THEN
    sow_date_cpft(:,:)        = 0.0
    tt_veg_cpft(:,:)          = 0.0
    tt_rep_cpft(:,:)          = 0.0
    latestharv_date_cpft(:,:) = 0.0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_crop IF 10 causing error.", errcode, "please check " //  &
                  "allocate_jules_arrays")
  END IF

ELSE
  ALLOCATE( dvi_cpft(land_pts,ncpft), stat = error )
  error_sum = error
  ALLOCATE( croprootc_cpft(land_pts,ncpft), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    dvi_cpft(:,:)       = 0.0
    croprootc_cpft(:,:) = 0.0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_crop ELSE 1 causing error.", errcode, "please " //       &
                  "check allocate_jules_arrays")
  END IF
END IF

ALLOCATE( phot(t_i_length * t_j_length), stat = error )
error_sum = error
ALLOCATE( dphotdt(t_i_length * t_j_length), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  phot(:)    = 0.0
  dphotdt(:) = 0.0
ELSE
  CALL ereport ("Something in crop_vars_mod module JULES " //                 &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF

IF ( l_irrig_dmd ) THEN
  ALLOCATE(irrtiles(npft_max), stat = error )
  error_sum = error
  IF ( error_sum == 0 ) THEN
    irrtiles(:) = 0
  ELSE
    CALL ereport ("Something in crop_vars_mod module JULES " //               &
                  "l_irrig_dmd IF causing error.", errcode, "please " //      &
                  "check allocate_jules_arrays")
  END IF
END IF

!  ====diag_swchs module JULES====
!Nothing to allocate

!  ====dust_parameters_mod module JULES====
!Nothing to allocate

!  ====fluxes module JULES====
! Runoff components.
ALLOCATE( sub_surf_roff_gb(land_pts), stat = error )
error_sum = error
ALLOCATE( surf_roff_gb(land_pts), stat = error )
error_sum = error_sum + error
! (Forcing) fluxes
ALLOCATE( alb_surft(land_pts,nsurft,4) )
error_sum = error_sum + error
ALLOCATE( tstar_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( e_sea_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( fsmc_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( ftl_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( le_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( h_sea_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( fqw_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( fqw_sicat(t_i_length,t_j_length,nice_use), stat = error )
error_sum = error_sum + error
ALLOCATE( ftl_sicat(t_i_length,t_j_length,nice_use), stat = error )
error_sum = error_sum + error
ALLOCATE( ecan_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( esoil_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( sea_ice_htf_sicat(t_i_length,t_j_length,nice), stat = error )
error_sum = error_sum + error
ALLOCATE( surf_ht_flux_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( surf_htf_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( snow_soil_htf(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( land_albedo_ij(t_i_length,t_j_length,4) )
error_sum = error_sum + error
ALLOCATE( ei_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( ei_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( ecan_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( esoil_ij_soilt(t_i_length,t_j_length,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( ext_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE( snowmelt_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( hf_snow_melt_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( radnet_surft(land_pts,nsurft) )
error_sum = error_sum + error
ALLOCATE( sw_surft(land_pts,nsurft) )
error_sum = error_sum + error
ALLOCATE( emis_surft(land_pts,nsurft) )
error_sum = error_sum + error
ALLOCATE( snow_melt_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( snomlt_sub_htf_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( tot_tfall_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( melt_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( z0m_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( z0h_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  sub_surf_roff_gb(:)          = 0.0
  surf_roff_gb(:)              = 0.0
  alb_surft(:,:,:)             = 0.0
  tstar_ij(:,:)                = 0.0
  e_sea_ij(:,:)                = 0.0
  fsmc_pft(:,:)                = 0.0
  ftl_surft(:,:)               = 0.0
  le_surft(:,:)                = 0.0
  h_sea_ij(:,:)                = 0.0
  fqw_surft(:,:)               = 0.0
  fqw_sicat(:,:,:)             = 0.0
  ftl_sicat(:,:,:)             = 0.0
  ecan_ij(:,:)                 = 0.0
  esoil_surft(:,:)             = 273.15
  sea_ice_htf_sicat(:,:,:)     = 0.0
  surf_ht_flux_ij(:,:)         = 0.0
  surf_htf_surft(:,:)          = 0.0
  snow_soil_htf(:,:)           = 0.0
  land_albedo_ij(:,:,:)        = 0.0
  ei_ij(:,:)                   = 0.0
  ei_surft(:,:)                = 0.0
  ecan_surft(:,:)              = 0.0
  esoil_ij_soilt(:,:,:)        = 0.0
  ext_soilt(:,:,:)             = 0.0
  snowmelt_ij(:,:)             = 0.0
  hf_snow_melt_gb(:)           = 0.0
  radnet_surft(:,:)            = 0.0
  sw_surft(:,:)                = 0.0
  emis_surft(:,:)              = 0.0
  snow_melt_gb(:)              = 0.0
  snomlt_sub_htf_gb(:)         = 0.0
  tot_tfall_gb(:)              = 0.0
  melt_surft(:,:)              = 0.0
  z0m_surft(:,:)               = 0.0
  z0h_surft(:,:)               = 0.0
ELSE
  CALL ereport ("Something in fluxes module JULES " //                        &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF

!  ====gridmean_fluxes module JULES====
ALLOCATE( fqw_1_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( ftl_1_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( taux_1_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( tauy_1_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  fqw_1_ij(:,:)                = 0.0
  ftl_1_ij(:,:)                = 0.0
  taux_1_ij(:,:)               = 0.0
  tauy_1_ij(:,:)               = 0.0
ELSE
  CALL ereport ("Something in gridmean_fluxes module JULES " //               &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF


ALLOCATE( sthuf_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error
IF ( error_sum == 0 ) THEN
  sthuf_soilt(:,:,:)  = 0.0
ELSE
  CALL ereport ("Something in fluxes module 2 JULES " //                      &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF

!  ====forcing module JULES====
! Forcing variables
ALLOCATE( qw_1_ij(t_i_length,t_j_length), stat = error )
error_sum = error
ALLOCATE( tl_1_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( u_0_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( v_0_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( u_1_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( v_1_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( pstar_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( ls_rain_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( con_rain_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( ls_snow_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( con_snow_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( sw_down_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( lw_down_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( diff_rad_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( diurnal_temperature_range_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  qw_1_ij(:,:)                       = 0.0
  tl_1_ij(:,:)                       = 0.0
  u_0_ij(:,:)                        = 0.0
  v_0_ij(:,:)                        = 0.0
  u_1_ij(:,:)                        = 0.0
  v_1_ij(:,:)                        = 0.0
  pstar_ij(:,:)                      = 0.0
  ls_rain_ij(:,:)                    = 0.0
  con_rain_ij(:,:)                   = 0.0
  ls_snow_ij(:,:)                    = 0.0
  con_snow_ij(:,:)                   = 0.0
  sw_down_ij(:,:)                    = 0.0
  lw_down_ij(:,:)                    = 0.0
  diff_rad_ij(:,:)                   = 0.0
  diurnal_temperature_range_ij(:,:)  = 0.0
ELSE
  CALL ereport ("Something in forcing module JULES " //                       &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF

!  ====deposition_mod module JULES====
! Atmospheric deposition variables.
IF ( l_deposition ) THEN

  ! Model-dependent variables.
  SELECT CASE ( dry_dep_model )
  CASE ( dry_dep_ukca_jules )
    ! Deposition variables.
    IF ( l_deposition_flux ) THEN
      ALLOCATE( tracer_field(land_pts,ndry_dep_species), stat = error )
      error_sum = error
      IF ( error_sum == 0 ) THEN
        tracer_field(:,:) = 0.0
      ELSE
        CALL ereport ("Something in deposition_mod module JULES " //          &
                      "causing error.",                                       &
                      errcode, "please check allocate_jules_arrays")
      END IF
    END IF  !  l_deposition_flux
  CASE DEFAULT
    CALL ereport ("Invalid value for dry_dep_model. ", errcode,               &
                   "Please check allocate_jules_arrays (deposition_mod).")
  END SELECT
END IF  !  l_deposition

!  ====deposition_species_mod module JULES====
! Atmospheric deposition parameters.
IF ( l_deposition ) THEN
  ALLOCATE( dep_species_name(ndry_dep_species), stat = error )
  error_sum = error

  ! Model-dependent variables.
  SELECT CASE ( dry_dep_model )
  CASE ( dry_dep_ukca_jules )
    ! Deposition parameters.
    ALLOCATE( dd_ice_coeff(ndry_dep_species,3), stat = error )
    error_sum = error_sum + error
    ALLOCATE( diffusion_coeff(ndry_dep_species), stat = error )
    error_sum = error_sum + error
    ALLOCATE( diffusion_corr(ndry_dep_species), stat = error )
    error_sum = error_sum + error
    ALLOCATE( r_tundra(ndry_dep_species), stat = error )
    error_sum = error_sum + error
    ALLOCATE( rsurf_std(ntype,ndry_dep_species), stat = error )
    error_sum = error_sum + error
    ! Always allocate space for H2 parameters, even if not used.
    ALLOCATE( h2dd_c(ntype), stat = error )
    error_sum = error_sum + error
    ALLOCATE( h2dd_m(ntype), stat = error )
    error_sum = error_sum + error
    ALLOCATE( h2dd_q(ntype), stat = error )
    error_sum = error_sum + error
    IF ( error_sum == 0 ) THEN
      dd_ice_coeff(:,:)   = 0.0
      dep_species_name(:) = ""
      diffusion_coeff(:)  = 0.0
      diffusion_corr(:)   = 0.0
      r_tundra(:)         = 0.0
      rsurf_std(:,:)      = 0.0
      h2dd_c(:)           = 0.0
      h2dd_m(:)           = 0.0
      h2dd_q(:)           = 0.0
    ELSE
      CALL ereport ("Something in deposition_species_mod module JULES " //    &
                    "causing error.",                                         &
                    errcode, "please check allocate_jules_arrays")
    END IF

  CASE DEFAULT
    CALL ereport ("Invalid value for dry_dep_model. ", errcode,               &
                   "Please check allocate_jules_arrays " //                   &
                   "(deposition_species_mod).")

  END SELECT
END IF  !  l_deposition

!  ====jules_sea_seaice_mod module JULES====
!Nothing to allocate

!  ====jules_soil_mod module JULES====
!Nothing to allocate

!  ====jules_surface_types_mod module JULES====
!Nothing to allocate

!  ====jules_vegetation_mod module JULES====
!Nothing to allocate

!  ====orog module JULES====
! Orographic roughness variables
ALLOCATE( sil_orog_land_gb(land_pts), stat = error )
error_sum = error
ALLOCATE( ho2r2_orog_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( h_blend_orog_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( z0m_eff_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  sil_orog_land_gb(:)  = 0.0
  ho2r2_orog_gb(:)  = 0.0
  h_blend_orog_ij(:,:)  = 0.0
  z0m_eff_ij(:,:)  = 0.0
ELSE
  CALL ereport ("Something in orog module JULES " //                          &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF

!  ====prognostics module JULES====
ALLOCATE( lai_pft(land_pts,npft), stat = error )
error_sum = error
ALLOCATE( canht_pft(land_pts,npft), stat = error )
error_sum = error_sum + error
ALLOCATE( smcl_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE( t_soil_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE( tsurf_elev_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( rgrain_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( snow_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( soot_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( tstar_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( canopy_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( canopy_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( ti_sicat(t_i_length,t_j_length,nice), stat = error )
error_sum = error_sum + error
ALLOCATE( z0msea_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( gs_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( gc_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( smc_soilt(land_pts,nsoilt), stat = error )
error_sum = error_sum + error
ALLOCATE( di_ncat_sicat(t_i_length,t_j_length,nice), stat = error )
error_sum = error_sum + error
ALLOCATE( k_sice_sicat(t_i_length,t_j_length,nice), stat = error )
error_sum = error_sum + error
ALLOCATE( snow_grnd_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( snow_mass_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( snow_mass_sea_sicat(t_i_length,t_j_length,nice_use), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  lai_pft(:,:)               = 0.0
  canht_pft(:,:)             = 0.0
  smcl_soilt(:,:,:)          = 0.0
  t_soil_soilt(:,:,:)        = 0.0
  tsurf_elev_surft(:,:)      = 0.0
  rgrain_surft(:,:)          = 0.0
  snow_surft(:,:)            = 0.0
  soot_ij(:,:)               = 0.0
  tstar_surft(:,:)           = 0.0
  canopy_surft(:,:)          = 0.0
  canopy_gb(:)               = 0.0
  ti_sicat(:,:,:)            = 0.0
  z0msea_ij(:,:)             = 0.0
  gs_gb(:)                   = 0.0
  gc_surft(:,:)              = 0.0
  smc_soilt(:,:)             = 0.0
  di_ncat_sicat(:,:,:)       = 0.0
  k_sice_sicat(:,:,:)        = 0.0
  snow_grnd_surft(:,:)       = 0.0
  snow_mass_ij(:,:)          = 0.0
  snow_mass_sea_sicat(:,:,:) = 0.0
ELSE
  CALL ereport ("Something in prognostics module JULES " //                   &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF

!  ====p_s_parms module JULES====
! Plant and soil parameters
ALLOCATE( albsoil_soilt(land_pts,nsoilt), stat = error )
error_sum = error
ALLOCATE( albobs_sw_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( albobs_vis_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( albobs_nir_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( catch_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( catch_snow_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( cosz_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( infil_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( z0_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( z0h_bare_surft(land_pts,nsurft), stat = error )
error_sum = error_sum + error
ALLOCATE( z0m_soil_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( sthu_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE( sthf_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE( sthu_min_soilt(land_pts,nsoilt,sm_levels), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  albsoil_soilt(:,:)    = 0.0
  albobs_sw_gb(:)       = 0.0
  albobs_vis_gb(:)      = 0.0
  albobs_nir_gb(:)      = 0.0
  catch_surft(:,:)      = 0.0
  catch_snow_surft(:,:) = 0.0
  cosz_ij(:,:)          = 0.0
  infil_surft(:,:)      = 0.0
  z0_surft(:,:)         = 0.0
  z0h_bare_surft(:,:)   = 0.0
  z0m_soil_gb(:)        = 0.0
  sthu_soilt(:,:,:)     = 0.0
  sthf_soilt(:,:,:)     = 0.0
  sthu_min_soilt(:,:,:) = 0.0
ELSE
  CALL ereport ("Something in p_s_parms module JULES " //                     &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF
! Soil ancillaries on soil carbon layers. Only used with ECOSSE.
IF ( soil_bgc_model == soil_model_ecosse ) THEN
  ALLOCATE( soil_ph_soilt(land_pts,nsoilt,dim_cslayer), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    soil_ph_soilt(:,:,:)   = 0.0
  ELSE
    CALL ereport ( RoutineName, errcode,                                      &
                   "Error related to soil ancillaries on soil carbon " //     &
                   "layers. Only used with ECOSSE." )
  END IF
END IF

!Initialise to zero to prevent junk data existing in tiles with zero
!occupancy
!z0_surft(:,:) = 0.0 !moved to in the loop 
  
IF ( l_use_pft_psi ) THEN
  ALLOCATE( v_close_pft(land_pts,sm_levels,npft), stat = error )
  error_sum = error_sum + error
  ALLOCATE( v_open_pft(land_pts,sm_levels,npft), stat = error )
  error_sum = error
  IF ( error_sum == 0 ) THEN
    v_close_pft(:,:,:) = 0.0
    v_open_pft(:,:,:) = 0.0
  ELSE
    CALL ereport ( "Something in l_use_pft_psi causing an error.", errcode,   &
                   " please check allocate_jules_arrays.")
  END IF
END IF

!  ====soil_ecosse_vars module JULES====
IF ( soil_bgc_model == soil_model_ecosse ) THEN
  ! Note that at present we allocate all ECOSSE variables even if this
  ! is a carbon-only run that does not require the nitrogen variables.
  ! ECOSSE prognostics.
  ALLOCATE( n_amm_soilt(land_pts,nsoilt,dim_cslayer), stat = error )
  error_sum = error
  ALLOCATE( n_nit_soilt(land_pts,nsoilt,dim_cslayer), stat = error )
  error_sum = error_sum + error
  ! ECOSSE time-averaged drivers.
  ALLOCATE( deposition_n_driver(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( qbase_l_driver(land_pts,nsoilt,sm_levels), stat = error )
  error_sum = error_sum + error
  ALLOCATE( sthf_driver(land_pts,nsoilt,sm_levels), stat = error )
  error_sum = error_sum + error
  ALLOCATE( sthu_driver(land_pts,nsoilt,sm_levels), stat = error )
  error_sum = error_sum + error
  ALLOCATE( tsoil_driver(land_pts,nsoilt,sm_levels), stat = error )
  error_sum = error_sum + error
  ALLOCATE( wflux_driver(land_pts,nsoilt,0:sm_levels), stat = error )
  error_sum = error_sum + error
  ! Diagnostics.
  ALLOCATE( co2_soil_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_denitrification_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_leach_amm_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_leach_nit_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n_nitrification_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( no_soil_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n2o_soil_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n2o_denitrif_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n2o_nitrif_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n2o_partial_nitrif_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( n2_denitrif_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( plant_input_c_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( plant_input_n_gb(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( soil_c_add(land_pts), stat = error )
  error_sum = error_sum + error
  ALLOCATE( soil_n_add(land_pts), stat = error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) THEN
    n_amm_soilt(:,:,:)     = 0.0
    n_nit_soilt(:,:,:)     = 0.0
    deposition_n_driver(:) = 0.0
    qbase_l_driver(:,:,:)  = 0.0
    sthf_driver(:,:,:)     = 0.0
    sthu_driver(:,:,:)     = 0.0
    tsoil_driver(:,:,:)    = 0.0
    wflux_driver(:,:,:)    = 0.0
    co2_soil_gb(:)         = 0.0
    n_denitrification_gb(:)  = 0.0
    n_leach_amm_gb(:)        = 0.0
    n_leach_nit_gb(:)        = 0.0
    n_nitrification_gb(:)    = 0.0
    n2o_nitrif_gb(:)         = 0.0
    n2o_partial_nitrif_gb(:) = 0.0
    n2o_denitrif_gb(:)       = 0.0
    n2_denitrif_gb(:)        = 0.0
    no_soil_gb(:)            = 0.0
    n2o_soil_gb(:)           = 0.0
    plant_input_c_gb(:)    = 0.0
    plant_input_n_gb(:)    = 0.0
    soil_c_add(:)          = 0.0
    soil_n_add(:)          = 0.0
  ELSE
    CALL ereport ( RoutineName, errcode,                                      &
                   "Error related to ECOSSE variables." )
  END IF
END IF

!  ====pdm_vars module JULES====
!       PDM variables
ALLOCATE( slope_gb(land_pts), stat = error )
error_sum = error_sum + error

!  ====u_v_grid module JULES====
! Grid-change variables
ALLOCATE( u_0_p_ij(t_i_length,t_j_length), stat = error )
error_sum = error
ALLOCATE( v_0_p_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( u_1_p_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( v_1_p_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
ALLOCATE( dtrdz_charney_grid_1_ij(t_i_length,t_j_length), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  u_0_p_ij(:,:)    = 0.0
  v_0_p_ij(:,:)    = 0.0
  u_1_p_ij(:,:)    = 0.0
  v_1_p_ij(:,:)    = 0.0
  dtrdz_charney_grid_1_ij(:,:) = 0.0
ELSE
  CALL ereport ("Something in u_v_grid module JULES " //                      &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF

!  ====jules_rivers module JULES====
ALLOCATE( tot_surf_runoff_gb(land_pts), stat = error )
error_sum = error
ALLOCATE( tot_sub_runoff_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( acc_lake_evap_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( rflow_gb(land_pts), stat = error )
error_sum = error_sum + error
ALLOCATE( rrun_gb(land_pts), stat = error )
error_sum = error_sum + error
IF ( error_sum == 0 ) THEN
  tot_surf_runoff_gb(:) = 0.0
  tot_sub_runoff_gb(:)  = 0.0
  acc_lake_evap_gb(:)   = 0.0
  rflow_gb(:)           = 0.0
  rrun_gb(:)            = 0.0
ELSE
  CALL ereport ("Something in jules rivers module JULES " //                  &
                "causing error.", errcode, "please check " //                 &
                "allocate_jules_arrays")
END IF

#endif

!-----------------------------------------------------------------------
! Write out an error if there was one
!-----------------------------------------------------------------------
! Check for error.
IF ( error_sum /= 0 )                                                         &
  CALL ereport("allocate_jules_arrays", errcode,                              &
               "Error allocating JULES model arrays")


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE allocate_jules_arrays

END MODULE allocate_jules_arrays_mod
