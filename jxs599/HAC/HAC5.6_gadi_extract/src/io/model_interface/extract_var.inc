#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

FUNCTION extract_var(var_id) RESULT(cube)

  !Use subroutines
USE gridbox_mean_mod, ONLY:                                                   &
  surftiles_to_gbm, soiltiles_to_gbm

!Use variables
!Science variables
USE water_constants_mod, ONLY:                                                &
  rho_water, tm

USE csigma, ONLY:                                                             &
  sbcon

USE jules_surface_types_mod, ONLY:                                            &
  ntype, npft, lake, ncpft

USE theta_field_sizes, ONLY:                                                  &
  t_i_length, t_j_length

USE jules_surface_mod, ONLY:                                                  &
  l_aggregate, i_aggregate_opt, l_point_data

USE jules_radiation_mod, ONLY:                                                &
  wght_alb

USE jules_vegetation_mod, ONLY:                                               &
  l_vegdrag_surft

USE c_z0h_z0m, ONLY:                                                          &
  z0h_z0m

USE model_grid_mod, ONLY:                                                     &
  latitude, longitude

USE ancil_info, ONLY:                                                         &
  dim_cslayer, frac_surft, land_pts, land_index, lice_index, nsurft,          &
  soil_pts, soil_index, surft_index, surft_pts, nsoilt, frac_soilt

USE p_s_parms, ONLY:                                                          &
  bexp_soilt, clay_soilt, sathh_soilt, satcon_soilt, smvcst_soilt,            &
  smvccl_soilt, smvcwt_soilt, hcap_soilt, hcon_soilt,                         &
  z0_surft, z0h_bare_surft, sthu_soilt, sthf_soilt,                           &
  catch_surft, cosz_ij, sthu_min_soilt, soil_ph_soilt

USE soil_ecosse_vars_mod, ONLY:                                               &
  co2_soil_gb, n2_denitrif_gb, n2o_denitrif_gb, n2o_nitrif_gb,                &
  n2o_partial_nitrif_gb, n2o_soil_gb, n_amm_soilt, n_denitrification_gb,      &
  n_leach_amm_gb, n_leach_nit_gb, n_nitrification_gb, n_nit_soilt,            &
  no_soil_gb, plant_input_c_gb, plant_input_n_gb

USE fluxes, ONLY:                                                             &
  alb_surft, fqw_surft, hf_snow_melt_gb, sub_surf_roff_gb,                    &
  surf_roff_gb, snomlt_sub_htf_gb, snow_melt_gb,                              &
  tot_tfall_gb, ecan_ij, ei_ij, esoil_ij_soilt,                               &
  land_albedo_ij, surf_ht_flux_ij,                                            &
  melt_surft, anthrop_heat_surft, emis_surft,                                 &
  ext_soilt, fsmc_pft, ftl_surft,le_surft, radnet_surft,                      &
  surf_ht_store_surft, surf_htf_surft, ecan_surft,                            &
  ei_surft, esoil_surft, rflow_gb, rrun_gb,                                   &
  snow_soil_htf, sw_surft, z0m_surft, z0h_surft

USE gridmean_fluxes, ONLY:                                                    &
  fqw_1_ij,ftl_1_ij,taux_1_ij,tauy_1_ij

USE prognostics, ONLY:                                                        &
  cs_pool_soilt, t_soil_soilt, gs_gb, sice_surft,                             &
  sliq_surft, snow_surft, snow_grnd_surft,                                    &
  rgrainl_surft, tsnow_surft, canht_pft, canopy_surft,                        &
  gc_surft, lai_pft, rgrain_surft,                                            &
  rho_snow_grnd_surft, smcl_soilt,                                            &
  snowdepth_surft, tstar_surft, nsnow_surft,                                  &
  wood_prod_fast_gb, wood_prod_med_gb,                                        &
  wood_prod_slow_gb, ds_surft, tsoil_deep_gb,                                 &
  n_inorg_soilt_lyrs, ns_pool_gb, tsurf_elev_surft

USE crop_vars_mod, ONLY:                                                      &
  sow_date_cpft, tt_veg_cpft, tt_rep_cpft, latestharv_date_cpft,              &
  croprootc => rootc_cpft,                                                    &
  harvc_cpft, reservec_cpft, dvi_cpft,                                        &
  croplai_cpft, cropcanht_cpft,                                               &
  yield_diag_cpft, stemc_diag_cpft,                                           &
  leafc_diag_cpft, sthu_irr_soilt,                                            &
  frac_irr_soilt, irrDaysDiag_gb, irrig_water_gb,                             &
  harvest_trigger_cpft, harvest_counter_cpft

USE trifctl, ONLY:                                                            &
  cv_gb, gpp_gb, lit_c_mn_gb, npp_gb, resp_p_gb,                              &
  resp_s_soilt, resp_s_dr_out_gb, c_veg_pft, g_leaf_pft,                      &
  g_leaf_day_pft, g_leaf_dr_out_pft, g_leaf_phen_pft,                         &
  gpp_pft, lit_c_pft, npp_dr_out_pft, npp_pft, resp_p_pft,                    &
  resp_w_dr_out_pft, resp_w_pft, lai_phen_pft, frac_agr_gb

USE trif_vars_mod, ONLY:                                                      &
  lit_c_ag_pft, lit_c_orig_pft, WP_fast_in_gb, WP_slow_in_gb,                 &
  WP_med_in_gb, WP_fast_out_gb, WP_med_out_gb, WP_slow_out_gb,                &
  pc_s_pft, lai_bal_pft, n_leaf_pft, n_root_pft, n_stem_pft, resp_r_pft,      &
  resp_l_pft, leafc_pft, rootc_pft, stemc_pft, woodc_pft,                     &
  droot_pft, dleaf_pft, dwood_pft, root_litc_pft, leaf_litc_pft,              &
  wood_litc_pft, root_litn_pft, leaf_litn_pft, wood_litn_pft,                 &
  litterc_pft, lit_n_pft, lit_n_t_gb, n_leach_soilt, n_fix_gb,                &
  n_fix_pft, n_gas_gb,                                                        &
  n_uptake_pft, n_uptake_gb, n_demand_pft, n_demand_gb,                       &
  exudates_pft, exudates_gb, littern_pft, n_uptake_growth_pft,                &
  n_demand_growth_pft, n_demand_spread_pft, n_demand_lit_pft,                 &
  n_veg_gb, n_veg_pft, dcveg_pft, dcveg_gb, dnveg_pft,                        &
  dnveg_gb, minl_n_gb, minl_n_pot_gb, immob_n_gb,                             &
  immob_n_pot_gb, fn_gb, resp_s_diag_gb, resp_s_pot_diag_gb,                  &
  deposition_n_gb, g_burn_pft, n_loss_gb, dpm_ratio_gb,                       &
  fapar_diag_pft, apar_diag_pft, apar_diag_gb, fao_et0,                       &
  n_fertiliser_pft, n_fertiliser_gb, lit_n_ag_pft, lit_n_ag_pft,              &
  root_abandon_pft, root_abandon_gb, harvest_pft, harvest_gb,                 &
  harvest_n_pft, harvest_n_gb, root_abandon_n_pft, root_abandon_n_gb,         &
  lit_n_orig_pft, frac_past_gb,                                               &
  cnsrv_carbon_veg2_gb, cnsrv_carbon_triffid_gb,                              &
  cnsrv_veg_triffid_gb, cnsrv_soil_triffid_gb, cnsrv_prod_triffid_gb,         &
  cnsrv_nitrogen_triffid_gb, cnsrv_vegN_triffid_gb, cnsrv_soilN_triffid_gb,   &
  cnsrv_n_inorg_triffid_gb, resp_s_to_atmos_gb,                               &
  npp_n, npp_n_gb, lit_c_fire_pft, burnt_carbon_dpm, lit_n_fire_pft,          &
  burnt_carbon_rpm, veg_c_fire_emission_gb, veg_c_fire_emission_pft,          &
  g_burn_gb

USE veg_param, ONLY:                                                          &
  secs_per_360days

USE jules_deposition_mod, ONLY:                                               &
  tracer_field

USE jules_soil_biogeochem_mod, ONLY: soil_bgc_model, soil_model_ecosse,       &
  soil_model_rothc

USE jules_soil_ecosse_mod, ONLY: dt_soilc

USE jules_soil_mod, ONLY:                                                     &
  dzsoil, zsmc, sm_levels, ns_deep

USE top_pdm, ONLY:                                                            &
  drain_soilt, fch4_wetl_soilt,                                               &
  fch4_wetl_cs_soilt, fch4_wetl_npp_soilt, fch4_wetl_resps_soilt,             &
  fsat_soilt, fwetl_soilt, qbase_soilt, qbase_zw_soilt,                       &
  dun_roff_soilt, sthzw_soilt, zw_soilt

USE jules_mod, ONLY:                                                          &
  snowdep_surft

USE forcing, ONLY:                                                            &
  con_rain_ij, con_snow_ij, ls_rain_ij, ls_snow_ij,                           &
  lw_down_ij, pstar_ij, qw_1_ij, sw_down_ij, tl_1_ij,                         &
  u_1_ij, v_1_ij, diurnal_temperature_range_ij

USE sf_diags_mod, ONLY:                                                       &
  sf_diag

USE jules_snow_mod, ONLY:                                                     &
  canSnowTile

USE jules_surface_mod, ONLY:                                                  &
  diff_frac

USE jules_rivers_mod, ONLY:                                                   &
  rivers_xgrid, rivers_ygrid

USE overbank_inundation_mod, ONLY:                                            &
  frac_fplain_lp

USE ozone_vars, ONLY:                                                         &
  flux_o3_pft, fo3_pft

USE bvoc_vars,                ONLY:                                           &
  isoprene_gb, isoprene_pft, terpene_gb , terpene_pft,                        &
  methanol_gb, methanol_pft, acetone_gb, acetone_pft

USE aero, ONLY:                                                               &
  co2_mmr

USE fire_mod,  ONLY:                                                          &
  fire_prog, fire_diag

USE fire_vars,  ONLY:                                                         &
  flammability_ft, burnt_area, burnt_area_ft,                                 &
  emitted_carbon, emitted_carbon_ft, emitted_carbon_DPM, emitted_carbon_RPM,  &
  fire_em_CO2, fire_em_CO2_ft, fire_em_CO2_DPM, fire_em_CO2_RPM,              &
  fire_em_CO, fire_em_CO_ft, fire_em_CO_DPM, fire_em_CO_RPM,                  &
  fire_em_CH4, fire_em_CH4_ft, fire_em_CH4_DPM, fire_em_CH4_RPM,              &
  fire_em_NOx, fire_em_NOx_ft, fire_em_NOx_DPM, fire_em_NOx_RPM,              &
  fire_em_SO2, fire_em_SO2_ft, fire_em_SO2_DPM, fire_em_SO2_RPM,              &
  fire_em_OC,  fire_em_OC_ft, fire_em_OC_DPM,  fire_em_OC_RPM,                &
  fire_em_BC, fire_em_BC_ft, fire_em_BC_DPM, fire_em_BC_RPM

USE imogen_progs, ONLY:                                                       &
  d_land_atmos_co2, d_ocean_atmos, c_emiss_out

!Others
USE io_constants, ONLY:                                                       &
  mdi
USE string_utils_mod, ONLY:                                                   &
  to_string
USE logging_mod, ONLY:                                                        &
  log_fatal
USE data_cube_mod, ONLY:                                                      &
  data_cube, cube_from_array, cube_safe_copy,                                 &
  cube_free, operator (*), operator (-)

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an identifer for a model variable, returns the data currently
!   associated with that variable as a cube
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
INTEGER, INTENT(IN) :: var_id
                     ! Identifies the variable to extract data from


! Return
TYPE(data_cube) :: cube  ! The extracted data


! Work variables
TYPE(data_cube) :: cube_land  ! Workspace cube for land data
                              ! This is required so that it can be deallocated
                              ! to avoid memory leaks

TYPE(data_cube) :: upward_cube, emis_cube, downward_cube, lw_down_cube
                              ! Work cubes for use in calculation of lw_net
                              ! Required to avoid memory leaks

REAL :: workspace_land(land_pts)  ! Used as a space for calculations
                                  ! before creating a cube

REAL :: workspace_surft(land_pts, nsurft)
                                ! Used in calculation of tile variables

REAL :: workspace_levs(land_pts, sm_levels)
                                ! Used in calculation of gridbox mean (GBM)
                                ! on soil levels

REAL :: workspace_cpft(land_pts, ncpft)
                                ! Used in calculation of crop tile variables

REAL :: workspace_cs(land_pts, dim_cslayer)
                                ! Used in calculation of soil C layer
                                ! variables.

REAL :: ones(land_pts, sm_levels)  ! An array full of ones to aid with
                                   ! calculation of first frozen/unfrozen
                                   ! layer
INTEGER :: layer(land_pts)  ! Used in calculation of first frozen/unfrozen
                            ! soil layer

REAL :: sum_frac(land_pts)  ! Used in calculation of fsmc_gb - the sum
                            ! of frac_surft over all pft tiles

REAL :: dz  ! Used in calculation of available soil moisture
            ! Size of current soil layer
REAL :: ztop  ! Used in calculation of available soil moisture
              ! Depth to top of current soil layer

INTEGER :: i,j,k,l,n,m  ! Index variables

REAL :: co2_tmp(land_pts) ! Used to out put the 1-d co2_mmr variable

!-----------------------------------------------------------------------------

! Initialise data to missing data value
workspace_land(:)     = mdi
workspace_surft(:,:)  = mdi
workspace_cpft(:,:)   = mdi
workspace_levs(:,:)   = mdi
workspace_cs(:,:)     = mdi
ones(:,:)             = 1.0

! We use the string identifier to search for CASEs in the SELECT, rather than
! being sensitive to the implementation of integer variable ids
SELECT CASE ( get_string_identifier(var_id) )
CASE ( 'latitude' )
  cube = cube_from_array(latitude)

CASE ( 'longitude' )
  cube = cube_from_array(longitude)

CASE ( 'frac' )
  cube_land = cube_from_array(frac_surft)
  cube = map_from_land(cube_land)

CASE ( 'b' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(bexp_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for b: " //                   &
                   "calculation is unphysical. Use b_soilt instead")
  END IF

CASE ( 'b_soilt' )
  cube_land = cube_from_array(bexp_soilt)
  cube = map_from_land(cube_land)

CASE ( 'clay' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, dim_cslayer
      workspace_levs(:,k) = soiltiles_to_gbm(clay_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:dim_cslayer))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for clay: " //                &
                   "calculation is unphysical. Use clay_soilt instead")
  END IF

CASE ( 'clay_soilt' )
  cube_land = cube_from_array(clay_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sathh' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(sathh_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sathh: " //               &
                   "calculation is unphysical. Use sathh_soilt instead")
  END IF

CASE ( 'sathh_soilt' )
  cube_land = cube_from_array(sathh_soilt)
  cube = map_from_land(cube_land)

CASE ( 'satcon' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(satcon_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for satcon: " //              &
                   "calculation is unphysical. Use satcon_soilt instead")
  END IF

CASE ( 'satcon_soilt' )
  cube_land = cube_from_array(satcon_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sm_sat' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(smvcst_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sm_sat: " //              &
                   "calculation is unphysical. Use sm_sat_soilt instead")
  END IF

CASE ( 'sm_sat_soilt' )
  cube_land = cube_from_array(smvcst_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sm_crit' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(smvccl_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sm_crit: " //             &
                   "calculation is unphysical. Use sm_crit_soilt instead")
  END IF

CASE ( 'sm_crit_soilt' )
  cube_land = cube_from_array(smvccl_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sm_wilt' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(smvcwt_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sm_wilt: " //             &
                   "calculation is unphysical. Use sm_wilt_soilt instead")
  END IF

CASE ( 'sm_wilt_soilt' )
  cube_land = cube_from_array(smvcwt_soilt)
  cube = map_from_land(cube_land)

CASE ( 'hcap' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(hcap_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for hcap: " //                &
                   "calculation is unphysical. Use hcap_soilt instead")
  END IF

CASE ( 'hcap_soilt' )
  cube_land = cube_from_array(hcap_soilt)
  cube = map_from_land(cube_land)

CASE ( 'hcon' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(hcon_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for hcon: " //                &
                   "calculation is unphysical. Use hcon_soilt instead")
  END IF

CASE ( 'hcon_soilt' )
  cube_land = cube_from_array(hcon_soilt)
  cube = map_from_land(cube_land)

  !-----------------------------------------------------------------------------
  ! Soil ancillaries on soil carbon layers.
  !-----------------------------------------------------------------------------
CASE ( 'soil_ph' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, dim_cslayer
      workspace_cs(:,k) = soiltiles_to_gbm(soil_ph_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_cs(:,1:dim_cslayer))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal( "extract_var",                                            &
                    "Cannot output GB mean value for soil_ph: "     //        &
                    "calculation is unphysical. Use soil_ph_soilt " //        &
                    "instead" )
  END IF

CASE ( 'soil_ph_soilt' )
  cube_land = cube_from_array(soil_ph_soilt)
  cube = map_from_land(cube_land)

CASE ( 'albedo_land' )
  ! Calculate the albedo as used in subroutine control when calculating the net
  ! shortwave on tiles
  ! Here we take the average of diffuse albedos in VIS and NIR
  cube_land = cube_from_array(                                                &
    surftiles_to_gbm( (wght_alb(1) * alb_surft(:,:,1) +                       &
                   wght_alb(2) * alb_surft(:,:,2) +                           &
                   wght_alb(3) * alb_surft(:,:,3) +                           &
                   wght_alb(4) * alb_surft(:,:,4)) )                          &
  )
  cube = map_from_land(cube_land)

CASE ( 'canopy_gb' )
  !     Don't use the canopy_gb variable, as this is calculated part-way through
  !     a timestep.
  cube_land = cube_from_array(surftiles_to_gbm( canopy_surft ))
  cube = map_from_land(cube_land)

CASE ( 'cs_gb' )
  IF ( nsoilt == 1) THEN
    !Case for a single soil tile
    cube_land = cube_from_array(SUM(sum(cs_pool_soilt(:,1,:,:), 3), 2))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output soil_CN_gb for soil tiled runs")
  END IF

CASE ( 'cv' )
  cube_land = cube_from_array(cv_gb)
  cube = map_from_land(cube_land)

CASE ( 'depth_frozen' )
  IF ( nsoilt == 1) THEN
    ! Get frozen depth from surface
    ! Start by assuming that there are no frozen layers anywhere,
    !and so frozen depth is 0
    workspace_land(:) = 0.0

    ! Get the first unfrozen layer for every land point using MINLOC
    ! We give an array of ones as the array to take mins from,
    ! so that we get the minimum index at which the mask is true
    layer(:) = MINLOC(ones, 2, t_soil_soilt(:,1,:) >= tm)

    ! If the layer found above is 0 for any point, that means
    !no unfrozen layers were found and hence the whole column is frozen
    WHERE ( layer(:) == 0 )
      workspace_land(:) = SUM(dzsoil(:))
    END WHERE

    ! At points where the layer found above is 1, that means
    ! no frozen layers so we can ignore those points
    ! Otherwise, interpolate to estimate depth of zero degC isotherm
    ! Do all points with the same unfrozen layer at once
    DO i = 2,sm_levels
      WHERE ( layer(:) == i )
        workspace_land(:) = SUM(dzSoil(1:i-1)) + dzSoil(i) *                  &
                          (tm - t_soil_soilt(:,1,i-1)) /                      &
                          (t_soil_soilt(:,1,i) - t_soil_soilt(:,1,i-1))
      END WHERE
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for depth_frozen: " //        &
                   "soil tiled diagnostic not yet implemented")
  END IF

CASE ( 'depth_frozen_sthf' )
  IF ( nsoilt == 1) THEN
    ! Get frozen depth from surface in terms on frozen/unfrozen soil moisture
    ! Method similar in construction to depth_frozen
    workspace_land(:) = 0.0
    layer(:) = MINLOC(ones, 2, sthf_soilt(:,1,:) == 0.0)
    WHERE ( layer(:) == 0 )
      workspace_land(:) = SUM(dzsoil(:)) - dzsoil(sm_levels) *                &
            (1 - sthf_soilt(:,1,sm_levels) /                                  &
            (sthf_soilt(:,1,sm_levels) + sthu_soilt(:,1,sm_levels) -          &
            sthu_min_soilt(:,1,sm_levels)))
      ! The bottom layer can be below freezing temp but only partially frozen.
    END WHERE
    DO i = 2,sm_levels
      WHERE ( layer(:) == i )
        workspace_land(:) = SUM(dzSoil(1:i-1)) - dzsoil(i-1) *                &
                (1 - sthf_soilt(:,1,i-1) / (sthf_soilt(:,1,i-1) +             &
                sthu_soilt(:,1,i-1) - sthu_min_soilt(:,1,i-1)))
      END WHERE
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for"                   //     &
                   "depth_frozen_sthf: "                               //     &
                   "soil tiled diagnostic not yet implemented")
  END IF

CASE ( 'depth_unfrozen' )
  IF ( nsoilt == 1) THEN
    ! Get unfrozen depth from surface
    ! See depth_frozen for method description, but swapping frozen for unfrozen
    ! and vica-versa
    workspace_land(:) = 0.0
    layer(:) = MINLOC(ones, 2, t_soil_soilt(:,1,:) < tm)
    WHERE ( layer(:) == 0 )
      workspace_land(:) = SUM(dzsoil(:))
    END WHERE
    DO i = 2,sm_levels
      WHERE ( layer(:) == i )
        workspace_land(:) = SUM(dzSoil(1:i-1)) + dzSoil(i) *                  &
                          (tm - t_soil_soilt(:,1,i-1)) /                      &
                          (t_soil_soilt(:,1,i) - t_soil_soilt(:,1,i-1))
      END WHERE
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for depth_unfrozen: "  //     &
                   "soil tiled diagnostic not yet implemented")
  END IF

CASE ( 'depth_unfrozen_sthf' )
  IF ( nsoilt == 1) THEN
    ! Get unfrozen depth from surface in terms of frozen/unfrozen
    ! soil moisture
    ! See depth_frozen for method, but swapping frozen for unfrozen
    ! & vica-versa
    workspace_land(:) = 0.0
    layer(:) = MINLOC(ones, 2, sthf_soilt(:,1,:) > 0.0)
    WHERE ( layer(:) == 0 )
      workspace_land(:) = SUM(dzsoil(:))
    END WHERE
    DO i = 1,sm_levels           ! HERE WE CAN NO LONGER IGNORE i=1.
      WHERE ( layer(:) == i )
        workspace_land(:) = SUM(dzSoil(1:i)) - dzSoil(i) *                    &
                            sthf_soilt(:,1,i) / (sthf_soilt(:,1,i) +          &
                            sthu_soilt(:,1,i) - sthu_min_soilt(:,1,i))
      END WHERE
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for"                   //     &
                   "depth_unfrozen_sthf: "                             //     &
                   "soil tiled diagnostic not yet implemented")
  END IF

CASE ( 'drain' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(drain_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for drain: "           //     &
                   "calculation is unphysical. Use drain_soilt instead")
  END IF

CASE ( 'drain_soilt' )
  cube_land = cube_from_array(drain_soilt)
  cube = map_from_land(cube_land)

CASE ( 'elake' )
  workspace_land(:) = 0.0
  IF ( .NOT. l_aggregate .AND. lake > 0 )                                     &
    workspace_land(:) = fqw_surft(:,lake) * frac_surft(:,lake)
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'emis_gb' )
  cube_land = cube_from_array(surftiles_to_gbm(emis_surft))
  cube = map_from_land(cube_land)

CASE ( 'fch4_wetl' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(fch4_wetl_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fch4_wetl: "       //     &
                   "Use fch4_wetl_soilt instead")
  END IF

CASE ( 'fch4_wetl_soilt' )
  cube_land = cube_from_array(fch4_wetl_soilt)
  cube = map_from_land(cube_land)

CASE ( 'fch4_wetl_cs' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(fch4_wetl_cs_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fch4_wetl_cs: "    //     &
                   "Use fch4_wetl_cs_soilt instead")
  END IF

CASE ( 'fch4_wetl_cs_soilt' )
  cube_land = cube_from_array(fch4_wetl_cs_soilt)
  cube = map_from_land(cube_land)

CASE ( 'fch4_wetl_npp' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(fch4_wetl_npp_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fch4_wetl_npp: "   //     &
                   "Use fch4_wetl_npp_soilt instead")
  END IF

CASE ( 'fch4_wetl_npp_soilt' )
  cube_land = cube_from_array(fch4_wetl_npp_soilt)
  cube = map_from_land(cube_land)

CASE ( 'fch4_wetl_resps' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(fch4_wetl_resps_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fch4_wetl_resps: " //     &
                   "Use fch4_wetl_resps_soilt instead")
  END IF

CASE ( 'fch4_wetl_resps_soilt' )
  cube_land = cube_from_array(fch4_wetl_resps_soilt)
  cube = map_from_land(cube_land)

CASE ( 'fsat' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(fsat_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fsat: "            //     &
                   "Use fsat_soilt instead")
  END IF

CASE ( 'fsat_soilt' )
  cube_land = cube_from_array(fsat_soilt)
  cube = map_from_land(cube_land)

CASE ( 'fsmc_gb' )
  ! Calculate gridbox mean over PFTs.
  ! Calculate the weighted sum over pfts
  workspace_land(:) = SUM(fsmc_pft * frac_surft(:,1:npft), 2)
  sum_frac(:) = SUM(frac_surft(:,1:npft), 2)
  ! Normalise to the vegetation fraction
  WHERE ( sum_frac > EPSILON(1.0) )
    workspace_land(:) = workspace_land(:) / sum_frac(:)
  ELSEWHERE
    ! Where there is no veg, set an impossible value
    workspace_land(:) = mdi
  END WHERE
  ! Copy the result onto the full grid
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'fwetl' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(fwetl_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for fwetl: "           //     &
                   "Use fwetl_soilt instead")
  END IF

CASE ( 'fwetl_soilt' )
  cube_land = cube_from_array(fwetl_soilt)
  cube = map_from_land(cube_land)

CASE ( 'gpp_gb' )
  cube_land = cube_from_array(gpp_gb)
  cube = map_from_land(cube_land)

CASE ( 'lai_gb' )
  ! n.b. lai_gb is calculated here to make sure it is
  ! consistent with lai_pft after lai_pft is updated
  workspace_land(:) = SUM(lai_pft * frac_surft(:,1:npft), 2)
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'et_stom_gb' )
  cube = cube_from_array(sf_diag%et_stom_ij)

CASE ( 'et_stom' )
  cube_land = cube_from_array(sf_diag%et_stom_surft)
  cube = map_from_land(cube_land)

CASE ( 'fprf' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(sf_diag%fprf)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output fprf with tiled soil"           )
  END IF

CASE ( 'fsth' )
  IF ( nsoilt == 1) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = sf_diag%fsth(:,k)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output fsth with tiled soil"           )
  END IF

CASE ( 'ftemp' )
  IF ( nsoilt == 1) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = sf_diag%ftemp(:,k)
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output ftemp with tiled soil"           )
  END IF


CASE ( 'gs' )
  cube_land = cube_from_array(gs_gb)
  cube = map_from_land(cube_land)

CASE ( 'hf_snow_melt' )
  cube_land = cube_from_array(hf_snow_melt_gb)
  cube = map_from_land(cube_land)

CASE ( 'land_index' )
  cube_land = cube_from_array(REAL(land_index))
  cube = map_from_land(cube_land)

CASE ( 'lice_index' )
  cube_land = cube_from_array(REAL(lice_index))
  cube = map_from_land(cube_land)

CASE ( 'lit_c_mean' )
  cube_land = cube_from_array(lit_c_mn_gb)
  cube = map_from_land(cube_land)

CASE ( 'lw_net' )
  ! Calculate gridbox mean upwards longwave
  cube_land = cube_from_array(                                                &
    sbcon * surftiles_to_gbm(emis_surft * tstar_surft**4)                     &
  )
  upward_cube = map_from_land(cube_land)
  CALL cube_free(cube_land)

  ! Calculate gridbox mean downward longwave
  cube_land = cube_from_array(surftiles_to_gbm(emis_surft))
  emis_cube = map_from_land(cube_land)
  lw_down_cube = cube_from_array(lw_down_ij)
  downward_cube = emis_cube * lw_down_cube

  ! Now calculate the net flux
  cube = downward_cube - upward_cube

  ! Free work cubes
  CALL cube_free(upward_cube)
  CALL cube_free(emis_cube)
  CALL cube_free(lw_down_cube)
  CALL cube_free(downward_cube)

CASE ( 'lw_up' )
  cube_land = cube_from_array(                                                &
    sbcon * surftiles_to_gbm(emis_surft * tstar_surft**4)                     &
  )
  cube = map_from_land(cube_land)

CASE ( 'npp_gb' )
  cube_land = cube_from_array(npp_gb)
  cube = map_from_land(cube_land)

CASE ( 'qbase' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(qbase_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for qbase: "           //     &
                   "Use qbase_soilt instead")
  END IF

CASE ( 'qbase_soilt' )
  cube_land = cube_from_array(qbase_soilt)
  cube = map_from_land(cube_land)

CASE ( 'qbase_zw' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(qbase_zw_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for qbase_zw: "        //     &
                   "Use qbase_zw_soilt instead")
  END IF

CASE ( 'qbase_zw_soilt' )
  cube_land = cube_from_array(qbase_zw_soilt)
  cube = map_from_land(cube_land)

CASE ( 'rad_net' )
  cube_land = cube_from_array(surftiles_to_gbm(radnet_surft))
  cube = map_from_land(cube_land)

CASE ( 'resp_p_gb' )
  cube_land = cube_from_array(resp_p_gb)
  cube = map_from_land(cube_land)

CASE ( 'resp_s_gb' )
  IF ( nsoilt == 1) THEN
    !Case for a single soil tile
    cube_land = cube_from_array(SUM(sum(resp_s_soilt(:,1,:,:), 3), 2))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output resp_s_gb for soil tiled runs")
  END IF

CASE ( 'resp_s_dr_out' )
  ! HACK: We only output the total respiration for now
  cube_land = cube_from_array(SUM(resp_s_dr_out_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'resp_s_diag' )
  cube_land = cube_from_array(resp_s_diag_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'resp_s_pot_diag' )
  cube_land = cube_from_array(resp_s_pot_diag_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'immob_n' )
  cube_land = cube_from_array(immob_n_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'immob_n_pot' )
  cube_land = cube_from_array(immob_n_pot_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'minl_n' )
  cube_land = cube_from_array(minl_n_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'minl_n_pot' )
  cube_land = cube_from_array(minl_n_pot_gb(:,:,1:4))
  cube = map_from_land(cube_land)

CASE ( 'ns' )
  cube_land = cube_from_array(ns_pool_gb(:,:,:))
  cube = map_from_land(cube_land)

CASE ( 'fN' )
  cube_land = cube_from_array(fn_gb(:,:))
  cube = map_from_land(cube_land)

CASE ( 'soil_CN' )
  IF ( nsoilt == 1) THEN
    !Case for a single soil tile
    cube_land = cube_from_array(cs_pool_soilt(:,1,:,:) / ns_pool_gb(:,:,:))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output soil_CN for soil tiled runs")
  END IF

CASE ( 'resp_s_diag_gb' )
  cube_land = cube_from_array(SUM(resp_s_diag_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'resp_s_pot_diag_gb' )
  cube_land = cube_from_array(SUM(resp_s_pot_diag_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'immob_n_gb' )
  cube_land = cube_from_array(SUM(immob_n_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'immob_n_pot_gb' )
  cube_land = cube_from_array(SUM(immob_n_pot_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'minl_n_gb' )
  cube_land = cube_from_array(SUM(minl_n_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'minl_n_pot_gb' )
  cube_land = cube_from_array(SUM(minl_n_pot_gb(:,:,5), 2))
  cube = map_from_land(cube_land)

CASE ( 'ns_gb' )
  cube_land = cube_from_array(SUM(sum(ns_pool_gb(:,:,1:4),dim = 3), 2))
  cube = map_from_land(cube_land)

CASE ( 'soil_CN_gb' )
  IF ( nsoilt == 1) THEN
    cube_land =                                                               &
      cube_from_array(SUM(sum(cs_pool_soilt(:,1,:,1:4),dim = 3), 2) /         &
      SUM(sum(ns_pool_gb(:,:,1:4),dim = 3), 2))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output soil_CN_gb for soil tiled runs")
  END IF

CASE ( 'runoff' )
  cube_land = cube_from_array(sub_surf_roff_gb(:) + surf_roff_gb(:))
  cube = map_from_land(cube_land)

CASE ( 'sat_excess_roff' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(dun_roff_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sat_excess_roff: " //     &
                   "Use sat_excess_roff_soilt instead")
  END IF

CASE ( 'sat_excess_roff_soilt' )
  cube_land = cube_from_array(dun_roff_soilt)
  cube = map_from_land(cube_land)

CASE ( 'smc_avail_top' )
  IF ( nsoilt == 1) THEN
    workspace_land(:) = 0.0
    ! We maintain the depth of the top of the current layer as we go down
    ! through the soil column
    ztop = 0.0
    DO k = 1,sm_levels
      ! If the top of this layer is below where we want to calculate to,
      !we are done
      IF ( ztop >= zsmc ) EXIT

      ! Calculate the amount of this layer that we want to take into account
      dz = dzsoil(k)
      ! If the layer contains the calculation limit, update dz
      IF ( ztop + dz > zsmc ) dz = zsmc - ztop

      ! Add on the contribution for this layer
      DO j = 1,soil_pts
        i = soil_index(j)
        workspace_land(i) = workspace_land(i) + rho_water * dz *              &
                            MAX(0.0,                                          &
                            sthu_soilt(i,1,k) * smvcst_soilt(i,1,k)           &
                            - smvcwt_soilt(i,1,k))
      END DO
      ztop = ztop + dzsoil(k)
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for"                   //     &
                   "smc_avail_top: "                                   //     &
                   "soil tiled diagnostic not yet implemented")
  END IF

CASE ( 'smc_avail_tot' )
  IF ( nsoilt == 1) THEN
    ! This is the same as smc_avail_top, but for whole column
    workspace_land(:) = 0.0
    DO k = 1,sm_levels
      DO j = 1,soil_pts
        i = soil_index(j)
        workspace_land(i) = workspace_land(i) + rho_water * dzsoil(k) *       &
                            MAX(0.0,                                          &
                            sthu_soilt(i,1,k) * smvcst_soilt(i,1,k) -         &
                            smvcwt_soilt(i,1,k))
      END DO
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for"                   //     &
                   "smc_avail_tot: "                                   //     &
                   "soil tiled diagnostic not yet implemented")
  END IF

CASE ( 'smc_tot' )
  workspace_land(:) = 0.0
  IF ( nsoilt == 1) THEN
    DO k = 1,sm_levels
      DO j = 1,soil_pts
        i = soil_index(j)
        workspace_land(i) = workspace_land(i) + rho_water * dzsoil(k) *       &
                            MAX(0.0,                                          &
                            (sthu_soilt(i,1,k) + sthf_soilt(i,1,k)) *         &
                            smvcst_soilt(i,1,k))
      END DO
    END DO
  ELSE
    !We can't take a gbm of smvcst_soilt so add up manually
    DO k = 1,sm_levels
      DO m = 1,nsoilt
        DO j = 1,soil_pts
          i = soil_index(j)
          workspace_land(i) = workspace_land(i) + (frac_soilt(i,m) *          &
                              (rho_water * dzsoil(k) *                        &
                              MAX(0.0,                                        &
                              (sthu_soilt(i,m,k) + sthf_soilt(i,m,k)) *       &
                              smvcst_soilt(i,m,k))))
        END DO !j
      END DO !m
    END DO !k
  END IF
  ! Map the result onto the full grid
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'snomlt_sub_htf' )
  cube_land = cube_from_array(snomlt_sub_htf_gb)
  cube = map_from_land(cube_land)

CASE ( 'snow_can_gb' )
  ! Only include tiles where canopy snow model is used
  cube_land = cube_from_array(surftiles_to_gbm(snow_surft, canSnowTile))
  cube = map_from_land(cube_land)

CASE ( 'snow_depth_gb' )
  cube_land = cube_from_array(surftiles_to_gbm(snowDepth_surft))
  cube = map_from_land(cube_land)

CASE ( 'snow_frac' )
  ! Sum frac over tiles with snow.
  workspace_land(:) = 0.0
  IF ( l_aggregate ) THEN
    WHERE ( snow_surft(:,1) + snow_grnd_surft(:,1) > EPSILON(1.0) )
      workspace_land(:) = 1.0
    END WHERE
  ELSE
    workspace_land(:) = SUM(frac_surft, 2, snow_surft +                       &
                            snow_grnd_surft > EPSILON(1.0))
  END IF
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'snow_grnd_gb' )
  ! Only include tiles where canopy snow model is used
  cube_land = cube_from_array(surftiles_to_gbm(snow_grnd_surft,canSnowTile))
  cube = map_from_land(cube_land)

CASE ( 'snow_ice_gb' )
  ! Calculate sum of sice along the snow layers dimension
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      i = surft_index(j,n)
      workspace_surft(i,n) = SUM(sice_surft(i,n,1:nsnow_surft(i,n)))
    END DO
  END DO
  cube_land = cube_from_array(surftiles_to_gbm(workspace_surft))
  cube = map_from_land(cube_land)

CASE ( 'snow_liq_gb' )
  ! Calculate sum of sliq along the snow layers dimension
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      i = surft_index(j,n)
      workspace_surft(i,n) = SUM(sliq_surft(i,n,1:nsnow_surft(i,n)))
    END DO
  END DO
  cube_land = cube_from_array(surftiles_to_gbm(workspace_surft))
  cube = map_from_land(cube_land)

CASE ( 'snow_melt_gb' )
  cube_land = cube_from_array(snow_melt_gb)
  cube = map_from_land(cube_land)

CASE ( 'soil_index' )
  DO l = 1,soil_pts
    workspace_land(soil_index(l)) = REAL(soil_index(l))
  END DO
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'sthzw' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(sthzw_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sthzw: "           //     &
                   "calculation is unphysical. Use sthzw_soilt instead")
  END IF

CASE ( 'sthzw_soilt' )
  cube_land = cube_from_array(sthzw_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sub_surf_roff' )
  cube_land = cube_from_array(sub_surf_roff_gb)
  cube = map_from_land(cube_land)

CASE ( 'surf_roff' )
  cube_land = cube_from_array(surf_roff_gb)
  cube = map_from_land(cube_land)

CASE ( 'swet_liq_tot' )
  IF ( nsoilt == 1) THEN
    ! Divide column unfrozen moisture content by saturated moisture content
    workspace_land(:) = 0.0
    DO j = 1,soil_pts
      i = soil_index(j)
      workspace_land(i) = SUM(dzsoil(:) * MAX(0.0,                            &
                        sthu_soilt(i,1,:) * smvcst_soilt(i,1,:)))             &
                        / SUM(dzsoil(:) * MAX(0.0, smvcst_soilt(i,1,:)))
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for"                   //     &
                   "swet_liq_tot: "                                    //     &
                   "soil tiled diagnostic not yet implemented")
  END IF

CASE ( 'swet_tot' )
  IF ( nsoilt == 1) THEN
    ! Divide column total moisture content by saturated moisture content
    workspace_land(:) = 0.0
    DO j = 1,soil_pts
      i = soil_index(j)
      workspace_land(i) = SUM(dzsoil(:) *                                     &
                            MAX(0.0,                                          &
                            (sthu_soilt(i,1,:) + sthf_soilt(i,1,:)) *         &
                            smvcst_soilt(i,1,:)))                             &
                            /SUM(dzsoil(:) * MAX(0.0, smvcst_soilt(i,1,:)))
    END DO
    ! Map the result onto the full grid
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for"                   //     &
                   "swet_liq_tot: "                                    //     &
                   "soil tiled diagnostic not yet implemented")
  END IF

CASE ( 'sw_net' )
  ! Calculate the albedo as used in subroutine control when calculating
  ! the net shortwave on tiles
  ! Here we take the average of diffuse albedos in VIS and NIR.
  DO l = 1,land_pts
    j = (land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    workspace_land(l) = ( 1.0 - (                                             &
                          wght_alb(1) * land_albedo_ij(i,j,1) +               &
                          wght_alb(2) * land_albedo_ij(i,j,2) +               &
                          wght_alb(3) * land_albedo_ij(i,j,3) +               &
                          wght_alb(4) * land_albedo_ij(i,j,4)                 &
                ) ) * sw_down_ij(i,j)

  END DO
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'NDVI_land' )
  ! Calculate NDVI on tiles
  ! Using the average of direct & diffuse albedos in VIS and NIR.
  DO l = 1,land_pts
    j = (land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    workspace_land(l) = (((land_albedo_ij(i,j,3) + land_albedo_ij(i,j,4)) -   &
                         (land_albedo_ij(i,j,1) + land_albedo_ij(i,j,2))) /   &
                        ((land_albedo_ij(i,j,3) + land_albedo_ij(i,j,4)) +    &
                         (land_albedo_ij(i,j,1) + land_albedo_ij(i,j,2))))

  END DO
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)


CASE ( 'tfall' )
  cube_land = cube_from_array(tot_tfall_gb)
  cube = map_from_land(cube_land)

CASE ( 'trad' )
  ! Assuming emissivity=1.
  cube_land = cube_from_array(( surftiles_to_gbm(tstar_surft**4) )**0.25)
  cube = map_from_land(cube_land)

CASE ( 'zw' )
  IF ( nsoilt == 1) THEN
    cube_land = cube_from_array(soiltiles_to_gbm(zw_soilt))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for zw: "              //     &
                   "Use zw_soilt instead")
  END IF

CASE ( 'zw_soilt' )
  cube_land = cube_from_array(zw_soilt)
  cube = map_from_land(cube_land)

CASE ( 'c_veg' )
  cube_land = cube_from_array(c_veg_pft)
  cube = map_from_land(cube_land)

CASE ( 'fapar' )
  cube_land = cube_from_array(fapar_diag_pft)
  cube = map_from_land(cube_land)

CASE ( 'apar' )
  cube_land = cube_from_array(apar_diag_pft)
  cube = map_from_land(cube_land)

CASE ( 'apar_gb' )
  cube_land = cube_from_array(apar_diag_gb)
  cube = map_from_land(cube_land)

CASE ( 'fao_et0' )
  cube_land = cube_from_array(fao_et0)
  cube = map_from_land(cube_land)

CASE ( 'canht' )
  cube_land = cube_from_array(canht_pft)
  cube = map_from_land(cube_land)

CASE ( 'flux_o3_stom' )
  cube_land = cube_from_array(flux_o3_pft)
  cube = map_from_land(cube_land)

CASE ( 'fsmc' )
  cube_land = cube_from_array(fsmc_pft)
  cube = map_from_land(cube_land)

CASE ( 'g_leaf' )
  cube_land = cube_from_array(g_leaf_pft)
  cube = map_from_land(cube_land)

CASE ( 'g_leaf_day' )
  cube_land = cube_from_array(g_leaf_day_pft)
  cube = map_from_land(cube_land)

CASE ( 'g_leaf_dr_out' )
  cube_land = cube_from_array(g_leaf_dr_out_pft)
  cube = map_from_land(cube_land)

CASE ( 'g_leaf_phen' )
  cube_land = cube_from_array(g_leaf_phen_pft)
  cube = map_from_land(cube_land)

CASE ( 'gpp' )
  cube_land = cube_from_array(gpp_pft)
  cube = map_from_land(cube_land)

CASE ( 'lai' )
  cube_land = cube_from_array(lai_pft)
  cube = map_from_land(cube_land)

CASE ( 'lai_phen' )
  cube_land = cube_from_array(lai_phen_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_c' )
  cube_land = cube_from_array(lit_c_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_c_ag' )
  cube_land = cube_from_array(lit_c_ag_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_c_orig' )
  cube_land = cube_from_array(lit_c_orig_pft)
  cube = map_from_land(cube_land)

CASE ( 'harvest' )
  cube_land = cube_from_array(harvest_pft)
  cube = map_from_land(cube_land)

CASE ( 'harvest_gb' )
  cube_land = cube_from_array(harvest_gb)
  cube = map_from_land(cube_land)

CASE ( 'harvest_n' )
  cube_land = cube_from_array(harvest_n_pft)
  cube = map_from_land(cube_land)

CASE ( 'harvest_n_gb' )
  cube_land = cube_from_array(harvest_n_gb)
  cube = map_from_land(cube_land)

CASE ( 'n_fertiliser' )
  cube_land = cube_from_array(n_fertiliser_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_fertiliser_gb' )
  cube_land = cube_from_array(n_fertiliser_gb)
  cube = map_from_land(cube_land)

CASE ( 'root_abandon' )
  cube_land = cube_from_array(root_abandon_pft)
  cube = map_from_land(cube_land)

CASE ( 'root_abandon_gb' )
  cube_land = cube_from_array(root_abandon_gb)
  cube = map_from_land(cube_land)

CASE ( 'root_abandon_n' )
  cube_land = cube_from_array(root_abandon_n_pft)
  cube = map_from_land(cube_land)

CASE ( 'root_abandon_n_gb' )
  cube_land = cube_from_array(root_abandon_n_gb)
  cube = map_from_land(cube_land)

CASE ( 'lit_n_ag' )
  cube_land = cube_from_array(lit_n_ag_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_n_orig' )
  cube_land = cube_from_array(lit_n_orig_pft)
  cube = map_from_land(cube_land)

CASE ( 'npp_dr_out' )
  cube_land = cube_from_array(npp_dr_out_pft)
  cube = map_from_land(cube_land)

CASE ( 'npp' )
  cube_land = cube_from_array(npp_pft)
  cube = map_from_land(cube_land)

CASE ( 'o3_exp_fac' )
  cube_land = cube_from_array(fo3_pft)
  cube = map_from_land(cube_land)

CASE ( 'co2_mmr' )
  co2_tmp(:)=co2_mmr
  cube_land = cube_from_array(co2_tmp)
  cube = map_from_land(cube_land)

CASE ( 'frac_agr' )
  cube_land = cube_from_array(frac_agr_gb)
  cube = map_from_land(cube_land)

CASE ( 'frac_past' )
  cube_land = cube_from_array(frac_past_gb)
  cube = map_from_land(cube_land)

CASE ( 'pc_s' )
  cube_land = cube_from_array(pc_s_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_leaf' )
  cube_land = cube_from_array(n_leaf_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_root' )
  cube_land = cube_from_array(n_root_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_stem' )
  cube_land = cube_from_array(n_stem_pft)
  cube = map_from_land(cube_land)

CASE ( 'lai_bal' )
  cube_land = cube_from_array(lai_bal_pft)
  cube = map_from_land(cube_land)

CASE ( 'resp_p' )
  cube_land = cube_from_array(resp_p_pft)
  cube = map_from_land(cube_land)

CASE ( 'resp_r' )
  cube_land = cube_from_array(resp_r_pft)
  cube = map_from_land(cube_land)

CASE ( 'resp_l' )
  cube_land = cube_from_array(resp_l_pft)
  cube = map_from_land(cube_land)

CASE ( 'resp_w_dr_out' )
  cube_land = cube_from_array(resp_w_dr_out_pft)
  cube = map_from_land(cube_land)

CASE ( 'resp_w' )
  cube_land = cube_from_array(resp_w_pft)
  cube = map_from_land(cube_land)

CASE ( 'cs' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(cs_pool_soilt(:,1,:,:))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output GB mean value for cs when nsoilt > 1")
  END IF

CASE ( 'resp_s_to_atmos_gb' )
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_rothc )
    cube_land = cube_from_array( SUM(resp_s_to_atmos_gb(:,:),2) )
  CASE ( soil_model_ecosse )
    ! Convert units from kg m-2 s-1 to kg m-2 (360 days)-1.
    cube_land = cube_from_array( co2_soil_gb(:) * secs_per_360days            &
                                 / dt_soilc )
  CASE DEFAULT
    CALL log_fatal("extract_var",                                             &
                   "Cannot output resps_to_atmos_gb - no code for " //        &
                   "this soil model.")
  END SELECT
  cube = map_from_land(cube_land)

CASE ( 'resp_s' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(resp_s_soilt(:,1,:,:))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output GB mean value for resp_s when nsoilt > 1")
  END IF

CASE ( 'wood_prod_fast' )
  cube_land = cube_from_array(wood_prod_fast_gb)
  cube = map_from_land(cube_land)

CASE ( 'wood_prod_med' )
  cube_land = cube_from_array(wood_prod_med_gb)
  cube = map_from_land(cube_land)

CASE ( 'wood_prod_slow' )
  cube_land = cube_from_array(wood_prod_slow_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_fast_in' )
  cube_land = cube_from_array(WP_fast_in_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_med_in' )
  cube_land = cube_from_array(WP_med_in_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_slow_in' )
  cube_land = cube_from_array(WP_slow_in_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_fast_out' )
  cube_land = cube_from_array(WP_fast_out_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_med_out' )
  cube_land = cube_from_array(WP_med_out_gb)
  cube = map_from_land(cube_land)

CASE ( 'WP_slow_out' )
  cube_land = cube_from_array(WP_slow_out_gb)
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_carbon_veg2' )
  cube_land = cube_from_array(REAL(cnsrv_carbon_veg2_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_carbon_triffid' )
  cube_land = cube_from_array(REAL(cnsrv_carbon_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_veg_triffid' )
  cube_land = cube_from_array(REAL(cnsrv_veg_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_soil_triffid' )
  cube_land = cube_from_array(REAL(cnsrv_soil_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_prod_triffid' )
  cube_land = cube_from_array(REAL(cnsrv_prod_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_nitrogen_triffid' )
  cube_land = cube_from_array(REAL(cnsrv_nitrogen_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_vegN_triffid' )
  cube_land = cube_from_array(REAL(cnsrv_vegN_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_soilN_triffid' )
  cube_land = cube_from_array(REAL(cnsrv_soilN_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'cnsrv_n_inorg_triffid' )
  cube_land = cube_from_array(REAL(cnsrv_n_inorg_triffid_gb))
  cube = map_from_land(cube_land)

CASE ( 'lit_c_fire' )
  cube_land = cube_from_array(REAL(lit_c_fire_pft))
  cube = map_from_land(cube_land)

CASE ( 'burnt_carbon_dpm' )
  cube_land = cube_from_array(REAL(burnt_carbon_dpm))
  cube = map_from_land(cube_land)

CASE ( 'lit_n_fire' )
  cube_land = cube_from_array(REAL(lit_n_fire_pft))
  cube = map_from_land(cube_land)

CASE ( 'burnt_carbon_rpm' )
  cube_land = cube_from_array(REAL(burnt_carbon_rpm))
  cube = map_from_land(cube_land)

CASE ( 'veg_c_fire_emission_gb' )
  cube_land = cube_from_array(REAL(veg_c_fire_emission_gb))
  cube = map_from_land(cube_land)

CASE ( 'veg_c_fire_emission_pft' )
  cube_land = cube_from_array(REAL(veg_c_fire_emission_pft))
  cube = map_from_land(cube_land)

CASE ( 'n_inorg' )
  cube_land = cube_from_array(n_inorg_soilt_lyrs(:,1,:))
  cube = map_from_land(cube_land)

CASE ( 'n_inorg_gb' )
  IF ( nsoilt == 1 ) THEN
    IF ( soil_bgc_model == soil_model_rothc ) THEN
      cube_land = cube_from_array( SUM( n_inorg_soilt_lyrs(:,1,:), 2 ) )
    ELSE IF (  soil_bgc_model == soil_model_ecosse ) THEN
      cube_land = cube_from_array( SUM( n_nit_soilt(:,1,:) +                  &
                                        n_amm_soilt(:,1,:), 2 ) )
    END IF
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for n_inorg_gb with " //      &
                   "nsoilt>1")
  END IF

CASE ( 'deposition_n' )
  cube_land = cube_from_array(deposition_n_gb)
  cube = map_from_land(cube_land)

CASE ( 'g_burn_pft' )
  cube_land = cube_from_array(REAL(g_burn_pft))
  cube = map_from_land(cube_land)

CASE ( 'g_burn_gb' )
  cube_land = cube_from_array(REAL(g_burn_gb))
  cube = map_from_land(cube_land)

CASE ( 'leafC' )
  cube_land = cube_from_array(leafc_pft)
  cube = map_from_land(cube_land)

CASE ( 'rootC' )
  cube_land = cube_from_array(rootc_pft)
  cube = map_from_land(cube_land)

CASE ( 'stemC' )
  cube_land = cube_from_array(stemc_pft)
  cube = map_from_land(cube_land)

CASE ( 'woodC' )
  cube_land = cube_from_array(woodc_pft)
  cube = map_from_land(cube_land)

CASE ( 'dleaf' )
  cube_land = cube_from_array(dleaf_pft)
  cube = map_from_land(cube_land)

CASE ( 'droot' )
  cube_land = cube_from_array(droot_pft)
  cube = map_from_land(cube_land)

CASE ( 'dwood' )
  cube_land = cube_from_array(dwood_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_uptake' )
  cube_land = cube_from_array(n_uptake_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_uptake_gb' )
  cube_land = cube_from_array(n_uptake_gb)
  cube = map_from_land(cube_land)

CASE ( 'leaf_litC' )
  cube_land = cube_from_array(leaf_litc_pft)
  cube = map_from_land(cube_land)

CASE ( 'root_litC' )
  cube_land = cube_from_array(root_litc_pft)
  cube = map_from_land(cube_land)

CASE ( 'wood_litC' )
  cube_land = cube_from_array(wood_litc_pft)
  cube = map_from_land(cube_land)

CASE ( 'leaf_litN' )
  cube_land = cube_from_array(leaf_litn_pft)
  cube = map_from_land(cube_land)

CASE ( 'root_litN' )
  cube_land = cube_from_array(root_litn_pft)
  cube = map_from_land(cube_land)

CASE ( 'wood_litN' )
  cube_land = cube_from_array(wood_litn_pft)
  cube = map_from_land(cube_land)

CASE ( 'litterC' )
  cube_land = cube_from_array(litterc_pft)
  cube = map_from_land(cube_land)

CASE ( 'litterN' )
  cube_land = cube_from_array(littern_pft)
  cube = map_from_land(cube_land)

CASE ( 'exudates' )
  cube_land = cube_from_array(exudates_pft)
  cube = map_from_land(cube_land)

CASE ( 'exudates_gb' )
  cube_land = cube_from_array(exudates_gb)
  cube = map_from_land(cube_land)

CASE ( 'npp_n' )
  cube_land = cube_from_array(npp_n)
  cube = map_from_land(cube_land)

CASE ( 'npp_n_gb' )
  cube_land = cube_from_array(npp_n_gb)
  cube = map_from_land(cube_land)

CASE ( 'n_demand' )
  cube_land = cube_from_array(n_demand_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_fix' )
  cube_land = cube_from_array(n_fix_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_fix_gb' )
  cube_land = cube_from_array(n_fix_gb)
  cube = map_from_land(cube_land)

CASE ( 'n_gas' )
  cube_land = cube_from_array(n_gas_gb(:,:))
  cube = map_from_land(cube_land)

CASE ( 'n_gas_gb' )
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_rothc )
    cube_land = cube_from_array(SUM(n_gas_gb(:,:),2))
  CASE ( soil_model_ecosse )
    cube_land = cube_from_array( n2_denitrif_gb + no_soil_gb                  &
                                                + n2o_soil_gb )
  END SELECT
  cube = map_from_land(cube_land)

CASE ( 'n_leach' )
  SELECT CASE ( soil_bgc_model )
  CASE ( soil_model_rothc )
    IF (nsoilt == 1) THEN
      cube_land = cube_from_array(n_leach_soilt(:,1))
    ELSE
      CALL log_fatal("extract_var",                                           &
              "Cannot output GB mean value for n_leach when nsoilt > 1")
    END IF
  CASE ( soil_model_ecosse )
    cube_land = cube_from_array( n_leach_amm_gb(:) + n_leach_nit_gb(:) )
  END SELECT
  cube = map_from_land(cube_land)

CASE ( 'n_demand_gb' )
  cube_land = cube_from_array(n_demand_gb)
  cube = map_from_land(cube_land)

CASE ( 'n_demand_growth' )
  cube_land = cube_from_array(n_demand_growth_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_uptake_growth' )
  cube_land = cube_from_array(n_uptake_growth_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_demand_lit' )
  cube_land = cube_from_array(n_demand_lit_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_demand_spread' )
  cube_land = cube_from_array(n_demand_spread_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_veg' )
  cube_land = cube_from_array(n_veg_pft)
  cube = map_from_land(cube_land)

CASE ( 'n_veg_gb' )
  cube_land = cube_from_array(n_veg_gb)
  cube = map_from_land(cube_land)

CASE ( 'n_loss' )
  cube_land = cube_from_array(n_loss_gb)
  cube = map_from_land(cube_land)

CASE ( 'dpm_ratio' )
  cube_land = cube_from_array(dpm_ratio_gb)
  cube = map_from_land(cube_land)

CASE ( 'dnveg' )
  cube_land = cube_from_array(dnveg_pft)
  cube = map_from_land(cube_land)

CASE ( 'dnveg_gb' )
  cube_land = cube_from_array(dnveg_gb)
  cube = map_from_land(cube_land)

CASE ( 'dcveg' )
  cube_land = cube_from_array(dcveg_pft)
  cube = map_from_land(cube_land)

CASE ( 'dcveg_gb' )
  cube_land = cube_from_array(dcveg_gb)
  cube = map_from_land(cube_land)

CASE ( 'lit_N' )
  cube_land = cube_from_array(lit_n_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_N_t' )
  cube_land = cube_from_array(lit_n_t_gb)
  cube = map_from_land(cube_land)

CASE ( 'lit_n' )
  cube_land = cube_from_array(lit_n_pft)
  cube = map_from_land(cube_land)

CASE ( 'lit_n_t' )
  cube_land = cube_from_array(lit_n_t_gb)
  cube = map_from_land(cube_land)

CASE ( 'con_rain' )
  cube = cube_from_array(con_rain_ij)

CASE ( 'con_snow' )
  cube = cube_from_array(con_snow_ij)

CASE ( 'cosz' )
  cube = cube_from_array(cosz_ij)

CASE ( 'diff_frac' )
  cube = cube_from_array(RESHAPE(diff_frac, (/ t_i_length, t_j_length /)))

CASE ( 'ecan_gb' )
  cube = cube_from_array(ecan_ij)

CASE ( 'ei_gb' )
  cube = cube_from_array(ei_ij)

CASE ( 'esoil_gb' )
  IF (nsoilt == 1) THEN
    cube = cube_from_array(esoil_ij_soilt(:,:,1))
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output GB mean value for esoil_gb when nsoilt > 1")
  END IF

CASE ( 'fqw_gb' )
  cube = cube_from_array(fqw_1_ij)

CASE ( 'ftl_gb' )
  cube = cube_from_array(ftl_1_ij)

CASE ( 'tau_gb' )
  cube = cube_from_array(sf_diag%tau_1)

CASE ( 'land_albedo_1' )
  cube = cube_from_array(land_albedo_ij(:,:,1))

CASE ( 'land_albedo_2' )
  cube = cube_from_array(land_albedo_ij(:,:,2))

CASE ( 'land_albedo_3' )
  cube = cube_from_array(land_albedo_ij(:,:,3))

CASE ( 'land_albedo_4' )
  cube = cube_from_array(land_albedo_ij(:,:,4))

CASE ( 'latent_heat' )
  cube = cube_from_array(sf_diag%latent_heat)

CASE ( 'ls_rain' )
  cube = cube_from_array(ls_rain_ij)

CASE ( 'ls_snow' )
  cube = cube_from_array(ls_snow_ij)

CASE ( 'lw_down' )
  cube = cube_from_array(lw_down_ij)

CASE ( 'precip' )
  cube = cube_from_array(ls_rain_ij + con_rain_ij +                           &
                         ls_snow_ij + con_snow_ij)

CASE ( 'pstar' )
  cube = cube_from_array(pstar_ij)

CASE ( 'q1p5m_gb' )
  cube = cube_from_array(sf_diag%q1p5m)

CASE ( 'qw1' )
  cube = cube_from_array(qw_1_ij)

CASE ( 'rainfall' )
  cube = cube_from_array(ls_rain_ij + con_rain_ij)

CASE ( 'snomlt_surf_htf' )
  cube = cube_from_array(sf_diag%snomlt_surf_htf)

CASE ( 'snowfall' )
  cube = cube_from_array(ls_snow_ij + con_snow_ij)

CASE ( 'snow_mass_gb' )
  !     Don't use the snow_mass_ij variable as that is calculated under control.
  cube_land = cube_from_array(surftiles_to_gbm(snow_grnd_surft + snow_surft))
  cube = map_from_land(cube_land)

CASE ( 'surf_ht_flux_gb' )
  cube = cube_from_array(surf_ht_flux_ij)

CASE ( 'sw_down' )
  cube = cube_from_array(sw_down_ij)

CASE ( 't1p5m_gb' )
  cube = cube_from_array(sf_diag%t1p5m)

CASE ( 'taux1' )
  cube = cube_from_array(taux_1_ij)

CASE ( 'tauy1' )
  cube = cube_from_array(tauy_1_ij)

CASE ( 'tl1' )
  cube = cube_from_array(tl_1_ij)

CASE ( 'tstar_gb' )
  !     Don't use tstar as that is calculated part-way through a timestep.
  cube_land = cube_from_array(surftiles_to_gbm( tstar_surft ))
  cube = map_from_land(cube_land)

CASE ( 'u1' )
  cube = cube_from_array(u_1_ij)

CASE ( 'u10m' )
  cube = cube_from_array(sf_diag%u10m)

CASE ( 'v1' )
  cube = cube_from_array(v_1_ij)

CASE ( 'v10m' )
  cube = cube_from_array(sf_diag%v10m)

CASE ( 'wind' )
  cube = cube_from_array(SQRT(u_1_ij**2 + v_1_ij**2))

CASE ( 'dt_range' )
  cube = cube_from_array(diurnal_temperature_range_ij)

CASE ( 'ext' )
  IF ( nsoilt == 1 ) THEN
    cube_land = cube_from_array(ext_soilt(:,1,:))
    cube      = map_from_land(cube_land)
  ELSE
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(ext_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube      = map_from_land(cube_land)
  END IF

CASE ( 'ext_soilt' )
  cube_land = cube_from_array(ext_soilt)
  cube = map_from_land(cube_land)

CASE ( 'smcl' )
  IF ( nsoilt == 1 ) THEN
    cube_land = cube_from_array(smcl_soilt(:,1,:))
    cube      = map_from_land(cube_land)
  ELSE
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(smcl_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube      = map_from_land(cube_land)
  END IF

CASE ( 'smcl_soilt' )
  cube_land = cube_from_array(smcl_soilt)
  cube = map_from_land(cube_land)

CASE ( 'soil_wet' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(sthu_soilt(:,:,k)) +             &
                            soiltiles_to_gbm(sthf_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for soil_wet: " //            &
                   "Diagnostic not yet implemented")
  END IF

CASE ( 'sthf' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(sthf_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sthf: " //                &
                   "calculation is unphysical. Use sthf_soilt instead")
  END IF

CASE ( 'sthf_soilt' )
  cube_land = cube_from_array(sthf_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sthu' )
  IF ( nsoilt == 1 ) THEN
    DO k = 1, sm_levels
      workspace_levs(:,k) = soiltiles_to_gbm(sthu_soilt(:,:,k))
    END DO
    cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
                   "Cannot output GB mean value for sthu: " //                &
                   "calculation is unphysical. Use sthu_soilt instead")
  END IF

CASE ( 'sthu_soilt' )
  cube_land = cube_from_array(sthu_soilt)
  cube = map_from_land(cube_land)

CASE ( 't_soil' )
  DO k = 1, sm_levels
    workspace_levs(:,k) = soiltiles_to_gbm(t_soil_soilt(:,:,k))
  END DO
  cube_land = cube_from_array(workspace_levs(:,1:sm_levels))
  cube = map_from_land(cube_land)

CASE ( 't_soil_soilt' )
  cube_land = cube_from_array(t_soil_soilt)
  cube = map_from_land(cube_land)

CASE ( 'tsoil_deep' )
  cube_land = cube_from_array(tsoil_deep_gb(:,1:ns_deep))
  cube = map_from_land(cube_land)

CASE ( 'alb_tile_1' )
  cube_land = cube_from_array(alb_surft(:,:,1))
  cube = map_from_land(cube_land)

CASE ( 'alb_tile_2' )
  cube_land = cube_from_array(alb_surft(:,:,2))
  cube = map_from_land(cube_land)

CASE ( 'alb_tile_3' )
  cube_land = cube_from_array(alb_surft(:,:,3))
  cube = map_from_land(cube_land)

CASE ( 'alb_tile_4' )
  cube_land = cube_from_array(alb_surft(:,:,4))
  cube = map_from_land(cube_land)

CASE ( 'anthrop_heat' )
  cube_land = cube_from_array(anthrop_heat_surft)
  cube = map_from_land(cube_land)

CASE ( 'canopy' )
  cube_land = cube_from_array(canopy_surft)
  cube = map_from_land(cube_land)

CASE ( 'catch' )
  cube_land = cube_from_array(catch_surft)
  cube = map_from_land(cube_land)

CASE ( 'ecan' )
  cube_land = cube_from_array(ecan_surft)
  cube = map_from_land(cube_land)

CASE ( 'ei' )
  cube_land = cube_from_array(ei_surft)
  cube = map_from_land(cube_land)

CASE ( 'emis' )
  cube_land = cube_from_array(emis_surft)
  cube = map_from_land(cube_land)

CASE ( 'esoil' )
  cube_land = cube_from_array(esoil_surft)
  cube = map_from_land(cube_land)

CASE ( 'fqw' )
  ! Note that fqw_surft does not (always) do this job!
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      l = surft_index(j,n)
      workspace_surft(l,n) = ecan_surft(l,n) + ei_surft(l,n) +                &
                             esoil_surft(l,n)
      ! Add lake evaporation
      IF ( n == lake )                                                        &
        workspace_surft(l,n) = workspace_surft(l,n) + fqw_surft(l,n)
    END DO
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'ftl' )
  cube_land = cube_from_array(ftl_surft)
  cube = map_from_land(cube_land)

CASE ( 'tau' )
  cube_land = cube_from_array(sf_diag%tau_surft)
  cube = map_from_land(cube_land)

CASE ( 'gc' )
  cube_land = cube_from_array(gc_surft)
  cube = map_from_land(cube_land)

CASE ( 'le' )
  cube_land = cube_from_array(le_surft)
  cube = map_from_land(cube_land)

CASE ( 'nsnow' )
  cube_land = cube_from_array(REAL(nsnow_surft))
  cube = map_from_land(cube_land)

CASE ( 'q1p5m' )
  cube_land = cube_from_array(sf_diag%q1p5m_surft)
  cube = map_from_land(cube_land)

CASE ( 'rad_net_tile' )
  cube_land = cube_from_array(radnet_surft)
  cube = map_from_land(cube_land)

CASE ( 'sw_surft' )
  cube_land = cube_from_array(sw_surft)
  cube = map_from_land(cube_land)

CASE ( 'lw_up_surft' )
  cube_land = cube_from_array(sf_diag%lw_up_surft)
  cube = map_from_land(cube_land)

CASE ( 'lw_down_surft' )
  cube_land = cube_from_array(sf_diag%lw_down_surft)
  cube = map_from_land(cube_land)

CASE ( 'rgrain' )
  cube_land = cube_from_array(rgrain_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_can_melt' )
  ! Only include tiles where canopy snow model is used
  DO i = 1,nsurft
    IF ( canSnowTile(i) ) workspace_surft(:,i) = melt_surft(:,i)
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_can' )
  ! Only include tiles where canopy snow model is used
  DO i = 1,nsurft
    IF ( canSnowTile(i) ) workspace_surft(:,i) = snow_surft(:,i)
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_depth' )
  cube_land = cube_from_array(snowdepth_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_grnd_rho' )
  cube_land = cube_from_array(rho_snow_grnd_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_grnd' )
  ! Only include tiles where canopy snow model is used
  DO i = 1,nsurft
    IF ( canSnowTile(i) ) workspace_surft(:,i) = snow_grnd_surft(:,i)
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_ground' )
  DO i = 1,nsurft
    IF ( canSnowTile(i) ) THEN
      workspace_surft(:,i) = snow_grnd_surft(:,i)
    ELSE
      workspace_surft(:,i) = snow_surft(:,i)
    END IF
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_ice_tile' )
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      i = surft_index(j,n)
      workspace_surft(i,n) = SUM(sice_surft(i,n,1:nsnow_surft(i,n)))
    END DO
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_liq_tile' )
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      i = surft_index(j,n)
      workspace_surft(i,n) = SUM(sliq_surft(i,n,1:nsnow_surft(i,n)))
    END DO
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_mass' )
  workspace_surft(:,:) = snow_surft(:,:)
  ! Add snow below canopy
  DO n = 1,nsurft
    IF ( canSnowTile(n) )                                                     &
      workspace_surft(:,n) = workspace_surft(:,n) + snow_grnd_surft(:,n)
  END DO
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_melt' )
  workspace_surft(:,:) = melt_surft(:,:)
  ! Add melting of snow below canopy.
  !     IF ( can_model == 4 )
  !workspace_surft(:,:) = workspace_surft(:,:) + snowGMeltDiag(:,:)
  cube_land = cube_from_array(workspace_surft)
  cube = map_from_land(cube_land)

CASE ( 'surf_ht_flux' )
  cube_land = cube_from_array(surf_htf_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_soil_htf' )
  cube_land = cube_from_array(snow_soil_htf)
  cube = map_from_land(cube_land)

CASE ( 'snice_smb_surft' )
  cube_land = cube_from_array(sf_diag%snice_smb_surft)
  cube = map_from_land(cube_land)

CASE ( 'snice_m_surft' )
  cube_land = cube_from_array(sf_diag%snice_m_surft)
  cube = map_from_land(cube_land)

CASE ( 'snice_freez_surft' )
  cube_land = cube_from_array(sf_diag%snice_freez_surft)
  cube = map_from_land(cube_land)

CASE ( 'snice_sicerate_surft' )
  cube_land = cube_from_array(sf_diag%snice_sicerate_surft)
  cube = map_from_land(cube_land)

CASE ( 'snice_sliqrate_surft' )
  cube_land = cube_from_array(sf_diag%snice_sliqrate_surft)
  cube = map_from_land(cube_land)

CASE ( 'snice_runoff_surft' )
  cube_land = cube_from_array(sf_diag%snice_runoff_surft)
  cube = map_from_land(cube_land)

CASE ( 'surf_ht_store' )
  cube_land = cube_from_array(surf_ht_store_surft)
  cube = map_from_land(cube_land)

CASE ( 't1p5m' )
  cube_land = cube_from_array(sf_diag%t1p5m_surft)
  cube = map_from_land(cube_land)

CASE ( 'tstar' )
  cube_land = cube_from_array(tstar_surft)
  cube = map_from_land(cube_land)

CASE ( 'tsurf_elev_surft' )
  cube_land = cube_from_array(tsurf_elev_surft)
  cube = map_from_land(cube_land)

CASE ( 'z0' )
  IF ( ANY(l_vegdrag_surft) ) THEN
    cube_land = cube_from_array(z0m_surft)
  ELSE
    cube_land = cube_from_array(z0_surft)
  END IF
  cube = map_from_land(cube_land)

CASE ( 'z0h' )
  !     This diagnostic is set from z0_surft unless separately
  !     aggregated.
  IF ( l_aggregate .AND. (i_aggregate_opt == 1) ) THEN
    cube_land = cube_from_array(z0h_bare_surft)
  ELSE
    IF ( ANY(l_vegdrag_surft) ) THEN
      cube_land = cube_from_array(z0h_surft)
    ELSE
      DO n = 1,nsurft
        workspace_surft(:,n) =  z0h_z0m(n) * z0_surft(:,n)
      END DO
      cube_land = cube_from_array(workspace_surft)
    END IF
  END IF
  cube = map_from_land(cube_land)

CASE ( 'tile_index' )
  cube_land = cube_from_array(REAL(surft_index))
  cube = map_from_land(cube_land)

CASE ( 'isoprene_gb' )
  cube_land = cube_from_array(isoprene_gb)
  cube = map_from_land(cube_land)

CASE ( 'isoprene' )
  cube_land = cube_from_array(isoprene_pft)
  cube = map_from_land(cube_land)

CASE ( 'terpene_gb' )
  cube_land = cube_from_array(terpene_gb)
  cube = map_from_land(cube_land)

CASE ( 'terpene' )
  cube_land = cube_from_array(terpene_pft)
  cube = map_from_land(cube_land)

CASE ( 'methanol_gb' )
  cube_land = cube_from_array(methanol_gb)
  cube = map_from_land(cube_land)

CASE ( 'methanol' )
  cube_land = cube_from_array(methanol_pft)
  cube = map_from_land(cube_land)

CASE ( 'acetone_gb' )
  cube_land = cube_from_array(acetone_gb)
  cube = map_from_land(cube_land)

CASE ( 'acetone' )
  cube_land = cube_from_array(acetone_pft)
  cube = map_from_land(cube_land)

CASE ( 'croprootc' )
  cube_land = cube_from_array(croprootc)
  cube = map_from_land(cube_land)

CASE ( 'cropharvc' )
  cube_land = cube_from_array(harvc_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropreservec' )
  cube_land = cube_from_array(reservec_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropdvi' )
  cube_land = cube_from_array(dvi_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropyield' )
  cube_land = cube_from_array(yield_diag_cpft)
  cube = map_from_land(cube_land)

CASE ( 'harvest_trigger' )
  workspace_cpft(:,:) = REAL(harvest_trigger_cpft(:,:))
  cube_land = cube_from_array(workspace_cpft)
  cube = map_from_land(cube_land)

CASE ( 'harvest_counter' )
  workspace_cpft(:,:) = REAL(harvest_counter_cpft(:,:))
  cube_land = cube_from_array(workspace_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropstemc' )
  cube_land = cube_from_array(stemc_diag_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropleafc' )
  cube_land = cube_from_array(leafc_diag_cpft)
  cube = map_from_land(cube_land)

CASE ( 'croplai' )
  cube_land = cube_from_array(croplai_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropcanht' )
  cube_land = cube_from_array(cropcanht_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropsowdate' )
  cube_land = cube_from_array(sow_date_cpft)
  cube = map_from_land(cube_land)

CASE ( 'croplatestharvdate' )
  cube_land = cube_from_array(latestharv_date_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropttveg' )
  cube_land = cube_from_array(tt_veg_cpft)
  cube = map_from_land(cube_land)

CASE ( 'cropttrep' )
  cube_land = cube_from_array(tt_rep_cpft)
  cube = map_from_land(cube_land)

CASE ( 'frac_irrig' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(frac_irr_soilt(:,1))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output frac_irrig when nsoilt > 1. " //                   &
            "Use frac_irrig_soilt instead")
  END IF

CASE ( 'frac_irrig_soilt' )
  cube_land = cube_from_array(frac_irr_soilt)
  cube = map_from_land(cube_land)

CASE ( 'sthu_irr' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(sthu_irr_soilt(:,1,:))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output sthu_irr when nsoilt > 1. " //                     &
            "Use sthu_irr_soilt instead")
  END IF

CASE ( 'sthu_irr_soilt' )
  cube_land = cube_from_array(sthu_irr_soilt)
  cube = map_from_land(cube_land)

CASE ( 'irrDaysDiag' )
  cube_land = cube_from_array(irrDaysDiag_gb)
  cube = map_from_land(cube_land)

CASE ( 'irrig_water' )
  cube_land = cube_from_array(irrig_water_gb)
  cube = map_from_land(cube_land)

CASE ( 'rgrainl' )
  cube_land = cube_from_array(rgrainl_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_ds' )
  cube_land = cube_from_array(ds_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_ice' )
  cube_land = cube_from_array(sice_surft)
  cube = map_from_land(cube_land)

CASE ( 'snow_liq' )
  cube_land = cube_from_array(sliq_surft)
  cube = map_from_land(cube_land)

CASE ( 'tsnow' )
  cube_land = cube_from_array(tsnow_surft)
  cube = map_from_land(cube_land)

CASE ( 'rflow' )
  cube_land = cube_from_array(rflow_gb)
  cube = map_from_land(cube_land)

CASE ( 'rrun' )
  cube_land = cube_from_array(rrun_gb)
  cube = map_from_land(cube_land)

CASE ( 'frac_fplain_lp' )
  cube_land = cube_from_array(frac_fplain_lp)
  cube = map_from_land(cube_land)

  !
  !-----------------------------------------------------------------------------
  ! Equivalent neutral winds
  !-----------------------------------------------------------------------------
CASE ( 'u10m_n' )
  cube = cube_from_array(sf_diag%u10m_n)

CASE ( 'v10m_n' )
  cube = cube_from_array(sf_diag%v10m_n)

CASE ( 'mu10m_n' )
  cube = cube_from_array(sf_diag%mu10m_n)

CASE ( 'mv10m_n' )
  cube = cube_from_array(sf_diag%mv10m_n)

  ! Fire module variables
CASE ('fire_mcarthur')
  cube_land = cube_from_array(fire_diag(:)%mcarthur%ffdi)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian_ffmc')
  cube_land = cube_from_array(fire_prog(:)%canadian%ffmc)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian_dmc')
  cube_land = cube_from_array(fire_prog(:)%canadian%dmc)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian_dc')
  cube_land = cube_from_array(fire_prog(:)%canadian%dc)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian_isi')
  cube_land = cube_from_array(fire_diag(:)%canadian%isi)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian_bui')
  cube_land = cube_from_array(fire_diag(:)%canadian%bui)
  cube      = map_from_land(cube_land)

CASE ('fire_canadian')
  cube_land = cube_from_array(fire_diag(:)%canadian%fwi)
  cube      = map_from_land(cube_land)

CASE ('fire_nesterov')
  cube_land = cube_from_array(fire_prog(:)%nesterov%findex)
  cube      = map_from_land(cube_land)

  ! INFERNO variables
CASE ( 'flammability' )
  cube_land = cube_from_array(flammability_ft)
  cube = map_from_land(cube_land)

CASE ( 'burnt_area_gb' )
  cube_land = cube_from_array(burnt_area)
  cube = map_from_land(cube_land)

CASE ( 'burnt_area' )
  cube_land = cube_from_array(burnt_area_ft)
  cube = map_from_land(cube_land)

CASE ( 'emitted_carbon_gb' )
  cube_land = cube_from_array(emitted_carbon)
  cube = map_from_land(cube_land)
CASE ( 'emitted_carbon' )
  cube_land = cube_from_array(emitted_carbon_ft)
  cube = map_from_land(cube_land)
CASE ( 'emitted_carbon_DPM' )
  cube_land = cube_from_array(emitted_carbon_DPM)
  cube = map_from_land(cube_land)
CASE ( 'emitted_carbon_RPM' )
  cube_land = cube_from_array(emitted_carbon_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_CO2_gb' )
  cube_land = cube_from_array(fire_em_CO2)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO2' )
  cube_land = cube_from_array(fire_em_CO2_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO2_DPM' )
  cube_land = cube_from_array(fire_em_CO2_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO2_RPM' )
  cube_land = cube_from_array(fire_em_CO2_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_CO_gb' )
  cube_land = cube_from_array(fire_em_CO)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO' )
  cube_land = cube_from_array(fire_em_CO_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO_DPM' )
  cube_land = cube_from_array(fire_em_CO_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CO_RPM' )
  cube_land = cube_from_array(fire_em_CO_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_CH4_gb' )
  cube_land = cube_from_array(fire_em_CH4)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CH4' )
  cube_land = cube_from_array(fire_em_CH4_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CH4_DPM' )
  cube_land = cube_from_array(fire_em_CH4_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_CH4_RPM' )
  cube_land = cube_from_array(fire_em_CH4_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_NOx_gb' )
  cube_land = cube_from_array(fire_em_NOx)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_NOx' )
  cube_land = cube_from_array(fire_em_NOx_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_NOx_DPM' )
  cube_land = cube_from_array(fire_em_NOx_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_NOx_RPM' )
  cube_land = cube_from_array(fire_em_NOx_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_SO2_gb' )
  cube_land = cube_from_array(fire_em_SO2)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_SO2' )
  cube_land = cube_from_array(fire_em_SO2_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_SO2_DPM' )
  cube_land = cube_from_array(fire_em_SO2_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_SO2_RPM' )
  cube_land = cube_from_array(fire_em_SO2_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_OC_gb' )
  cube_land = cube_from_array(fire_em_OC)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_OC' )
  cube_land = cube_from_array(fire_em_OC_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_OC_DPM' )
  cube_land = cube_from_array(fire_em_OC_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_OC_RPM' )
  cube_land = cube_from_array(fire_em_OC_RPM)
  cube = map_from_land(cube_land)

CASE ( 'fire_em_BC_gb' )
  cube_land = cube_from_array(fire_em_BC)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_BC' )
  cube_land = cube_from_array(fire_em_BC_ft)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_BC_DPM' )
  cube_land = cube_from_array(fire_em_BC_DPM)
  cube = map_from_land(cube_land)
CASE ( 'fire_em_BC_RPM' )
  cube_land = cube_from_array(fire_em_BC_RPM)
  cube = map_from_land(cube_land)

  !-----------------------------------------------------------------------------
  ! Variables available with any multi-pool soil C model (e.g. RothC, ECOSSE).
  !-----------------------------------------------------------------------------
  ! Soil prognostic variables and single pools from multi-pool prognostic
  ! variables. These are all coded for a single soil tile only.
  !-----------------------------------------------------------------------------
CASE ( 'c_bio' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(cs_pool_soilt(:,1,:,3))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output c_bio when nsoilt > 1")
  END IF

CASE ( 'c_dpm' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(cs_pool_soilt(:,1,:,1))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output c_dpm when nsoilt > 1")
  END IF

CASE ( 'c_hum' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(cs_pool_soilt(:,1,:,4))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output c_hum when nsoilt > 1")
  END IF

CASE ( 'c_rpm' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(cs_pool_soilt(:,1,:,2))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output c_rpm when nsoilt > 1")
  END IF

CASE ( 'n_bio' )
  cube_land = cube_from_array(ns_pool_gb(:,:,3))
  cube = map_from_land(cube_land)

CASE ( 'n_dpm' )
  cube_land = cube_from_array(ns_pool_gb(:,:,1))
  cube = map_from_land(cube_land)

CASE ( 'n_hum' )
  cube_land = cube_from_array(ns_pool_gb(:,:,4))
  cube = map_from_land(cube_land)

CASE ( 'n_rpm' )
  cube_land = cube_from_array(ns_pool_gb(:,:,2))
  cube = map_from_land(cube_land)

  ! Gridbox totals of multi-pool soil prognostic variables.
  ! All are coded for a single soil tile only.
CASE ( 'c_bio_gb' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array( SUM(cs_pool_soilt(:,1,:,3),2) )
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output c_bio_gb when nsoilt > 1")
  END IF

CASE ( 'c_dpm_gb' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array( SUM(cs_pool_soilt(:,1,:,1),2) )
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output c_dpm_gb when nsoilt > 1")
  END IF

CASE ( 'c_hum_gb' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array( SUM(cs_pool_soilt(:,1,:,4),2) )
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output c_hum_gb when nsoilt > 1")
  END IF

CASE ( 'c_rpm_gb' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array( SUM(cs_pool_soilt(:,1,:,2),2) )
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output c_rpm_gb when nsoilt > 1")
  END IF

CASE ( 'n_bio_gb' )
  cube_land = cube_from_array( SUM(ns_pool_gb(:,:,3),2) )
  cube = map_from_land(cube_land)

CASE ( 'n_dpm_gb' )
  cube_land = cube_from_array( SUM(ns_pool_gb(:,:,1),2) )
  cube = map_from_land(cube_land)

CASE ( 'n_hum_gb' )
  cube_land = cube_from_array( SUM(ns_pool_gb(:,:,4),2) )
  cube = map_from_land(cube_land)

CASE ( 'n_rpm_gb' )
  cube_land = cube_from_array( SUM(ns_pool_gb(:,:,2),2) )
  cube = map_from_land(cube_land)

  ! Gridbox totals of derived variables for multi-pool soil models.
CASE ( 'n_soil_gb' )
  ! Coded for a single soil tile only.
  IF (nsoilt == 1) THEN
    ! Sum organic pools.
    workspace_land(:) = SUM( SUM(ns_pool_gb(:,:,:),3), 2 ) 
    ! Add inorganic pools.
    IF ( soil_bgc_model == soil_model_rothc ) THEN
      workspace_land(:) = workspace_land(:) +                                 &
                          SUM( n_inorg_soilt_lyrs(:,1,:), 2 )
    ELSE IF ( soil_bgc_model == soil_model_ecosse ) THEN
      workspace_land(:) = workspace_land(:) +                                 &
                          SUM( n_amm_soilt(:,1,:) + n_nit_soilt(:,1,:), 2 )
    END IF
    !     Map the result onto the full grid.
    cube_land = cube_from_array(workspace_land)
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output n_soil_gb when nsoilt > 1")
  END IF

  ! ECOSSE prognostic variables.
  ! All are coded for a single soil tile only.
CASE ( 'n_amm' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(n_amm_soilt(:,1,:))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output n_amm when nsoilt > 1")
  END IF

CASE ( 'n_nit' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array(n_nit_soilt(:,1,:))
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output n_nit when nsoilt > 1")
  END IF

  ! Gridbox totals of ECOSSE prognostic variables.
CASE ( 'n_amm_gb' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array( SUM( SUM(n_amm_soilt,2), 2 ) )
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output n_amm_gb when nsoilt > 1")
  END IF

CASE ( 'n_nit_gb' )
  IF (nsoilt == 1) THEN
    cube_land = cube_from_array( SUM( SUM(n_nit_soilt,2), 2 ) )
    cube = map_from_land(cube_land)
  ELSE
    CALL log_fatal("extract_var",                                             &
            "Cannot output n_nit_gb when nsoilt > 1")
  END IF

  !   ECOSSE fluxes.
CASE ( 'plant_input_c_gb' )
  cube_land = cube_from_array( plant_input_c_gb(:) )
  cube = map_from_land(cube_land)

CASE ( 'plant_input_n_gb' )
  cube_land = cube_from_array( plant_input_n_gb(:) )
  cube = map_from_land(cube_land)

  ! Imogen prognostics
CASE ( 'd_land_atmos_co2' )
  workspace_land(:) = d_land_atmos_co2
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'd_ocean_atmos' )
  workspace_land(:) = d_ocean_atmos
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)

CASE ( 'c_emiss_out' )
  workspace_land(:) = c_emiss_out
  cube_land = cube_from_array(workspace_land)
  cube = map_from_land(cube_land)
  ! Imogen end

CASE ( 'n_denitrif_gb' )
  cube_land = cube_from_array( n_denitrification_gb )
  cube = map_from_land(cube_land)

CASE ( 'n_nitrif_gb' )
  cube_land = cube_from_array( n_nitrification_gb )
  cube = map_from_land(cube_land)

CASE ( 'no_soil_gb' )
  cube_land = cube_from_array( no_soil_gb )
  cube = map_from_land(cube_land)

CASE ( 'n2_denitrif_gb' )
  cube_land = cube_from_array( n2_denitrif_gb )
  cube = map_from_land(cube_land)

CASE ( 'n2o_denitrif_gb' )
  cube_land = cube_from_array( n2o_denitrif_gb )
  cube = map_from_land(cube_land)

CASE ( 'n2o_nitrif_gb' )
  cube_land = cube_from_array( n2o_nitrif_gb )
  cube = map_from_land(cube_land)

CASE ( 'n2o_part_nitrif_gb' )
  cube_land = cube_from_array( n2o_partial_nitrif_gb )
  cube = map_from_land(cube_land)

CASE ( 'n2o_soil_gb' )
  cube_land = cube_from_array( n2o_soil_gb )
  cube = map_from_land(cube_land)

  ! ECOSSE leaching diagnostics.
CASE ( 'n_leach_amm_gb' )
  cube_land = cube_from_array( n_leach_amm_gb )
  cube = map_from_land(cube_land)

CASE ( 'n_leach_nit_gb' )
  cube_land = cube_from_array( n_leach_nit_gb )
  cube = map_from_land(cube_land)

  ! Atmospheric deposition variables.
CASE ( 'tracer_field' )
  cube_land = cube_from_array( tracer_field )
  cube = map_from_land(cube_land)

CASE DEFAULT
  CALL log_fatal("extract_var",                                               &
                 "Unrecognised variable for output - '" //                    &
                 TRIM(get_string_identifier(var_id)) // "'. " //              &
                 "See docs for available variables")
END SELECT

!-----------------------------------------------------------------------------
! Free the land cube
!
! Note that is safe to call this routine even if the cube has not been
! allocated - it just won't do anything
!-----------------------------------------------------------------------------
CALL cube_free(cube_land)


RETURN

END FUNCTION extract_var
#endif
