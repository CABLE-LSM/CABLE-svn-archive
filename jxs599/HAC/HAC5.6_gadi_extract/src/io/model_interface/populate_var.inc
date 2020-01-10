#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE populate_var(var_id, cube, const_val)

  !Science variables
USE ancil_info, ONLY:                                                         &
  dim_cslayer, land_pts, frac_surft, nsurft, z1_tq_ij, nsoilt

USE boundary_layer_var_mod, ONLY:                                             &
  dzl, zh

USE model_grid_mod, ONLY:                                                     &
  latitude, longitude

USE coastal, ONLY:                                                            &
  flandg

USE c_elevate, ONLY:                                                          &
  surf_hgt_surft, z_land_land

USE p_s_parms, ONLY:                                                          &
  albsoil_soilt, albobs_sw_gb, albobs_vis_gb,                                 &
  albobs_nir_gb, bexp_soilt, sathh_soilt, satcon_soilt,                       &
  smvcst_soilt, smvccl_soilt, smvcwt_soilt, hcap_soilt, hcon_soilt,           &
  sthu_soilt, clay_soilt, soil_ph_soilt

USE soil_ecosse_vars_mod, ONLY:                                               &
  n_amm_soilt, n_nit_soilt

USE top_pdm, ONLY:                                                            &
  fexp_soilt, ti_mean_soilt, ti_sig_soilt, sthzw_soilt, zw_soilt

USE pdm_vars, ONLY:                                                           &
  slope_gb

USE trifctl, ONLY:                                                            &
  frac_agr_gb, cv_gb

USE trif_vars_mod, ONLY:                                                      &
  deposition_n_gb, frac_past_gb, g_burn_pft

USE urban_param, ONLY:                                                        &
  wrr_gb, hwr_gb, hgt_gb, ztm_gb, disp_gb, albwl_gb,                          &
  albrd_gb, emisw_gb, emisr_gb

USE prognostics, ONLY:                                                        &
  canopy_surft, cs_pool_soilt, gs_gb, snow_surft,                             &
  t_soil_soilt, tstar_surft, lai_pft, canht_pft,                              &
  rgrain_surft, rho_snow_grnd_surft, snowdepth_surft,                         &
  snow_grnd_surft, nsnow_surft, sice_surft,                                   &
  sliq_surft, tsnow_surft, rgrainl_surft, ds_surft,                           &
  wood_prod_fast_gb,wood_prod_med_gb,                                         &
  wood_prod_slow_gb, frac_agr_prev_gb, tsoil_deep_gb,                         &
  n_inorg_soilt_lyrs, ns_pool_gb, frac_past_prev_gb,                          &
  tsurf_elev_surft, l_broadcast_soilt

USE crop_vars_mod, ONLY:                                                      &
  dvi_cpft, rootc_cpft, harvc_cpft,                                           &
  reservec_cpft, croplai_cpft, cropcanht_cpft,                                &
  sow_date_cpft, tt_veg_cpft, tt_rep_cpft, latestharv_date_cpft,              &
  sthu_irr_soilt, frac_irr_all

USE forcing, ONLY:                                                            &
  pstar_ij, qw_1_ij, tl_1_ij, lw_down_ij, sw_down_ij,                         &
  diff_rad_ij, ls_rain_ij, ls_snow_ij, con_rain_ij,                           &
  con_snow_ij, u_1_ij, v_1_ij, diurnal_temperature_range_ij

USE ozone_vars, ONLY:                                                         &
  o3_gb

USE jules_deposition_mod, ONLY:                                               &
  tracer_field

USE jules_soil_mod, ONLY:                                                     &
  sm_levels, l_tile_soil, l_broadcast_ancils

USE jules_rivers_mod, ONLY:                                                   &
  rivers_dir, rivers_seq, rivers_dra,                                         &
  rivers_xgrid, rivers_ygrid,                                                 &
  rivers_lat2d, rivers_lon2d,                                                 &
  rivers_sto_rp, rfm_surfstore_rp,                                            &
  rfm_substore_rp, rfm_flowin_rp, rfm_bflowin_rp

USE overbank_inundation_mod, ONLY:                                            &
  logn_mean, logn_stdev  

USE aero, ONLY:                                                               &
  co2_mmr

USE fire_mod, ONLY:                                                           &
  fire_prog

USE metstats_mod, ONLY:                                                       &
  metstats_prog

USE fire_vars, ONLY:                                                          &
  pop_den, flash_rate

!Others
USE string_utils_mod, ONLY:                                                   &
  to_string

USE data_cube_mod, ONLY:                                                      &
  data_cube, cube_get_data, cube_free, cube_create

USE update_mod, ONLY: sthuf_soilt

!CABLE_LSM: tiles soil prognostics to be initialized from IO
USE allocate_cable_progs_mod, ONLY : SoilTemp_CABLE,  SoilMoisture_CABLE,     &
                          FrozenSoilFrac_CABLE, SnowDepth_CABLE,               & 
                          SnowMass_CABLE,  SnowTemp_CABLE,  SnowDensity_CABLE, & 
                          ThreeLayerSnowFlag_CABLE, OneLyrSnowDensity_CABLE, &
                          SnowAge_CABLE

IMPLICIT NONE

! Argument types
INTEGER, INTENT(IN) :: var_id  ! Identifies the variable to fill
TYPE(data_cube), INTENT(IN), OPTIONAL :: cube
                     ! The data to put in to the variable as a cube
REAL, INTENT(IN), OPTIONAL :: const_val
                     ! A constant value to fill all elements of the variable
                     ! with


! Work variables
TYPE(data_cube) :: cube_land  ! Workspace cube for land and river data
                              ! This is required so that it can be deallocated
                              ! to avoid memory leaks

REAL :: nsnow_real(land_pts, nsurft)  ! Real version of nsnow

REAL :: co2_tmp(land_pts)

!-----------------------------------------------------------------------------


IF ( .NOT. PRESENT(cube) .AND. .NOT. PRESENT(const_val) )                     &
  CALL log_fatal("populate_var",                                              &
                 "Neither data cube or const_val have been provided for " //  &
                 "variable '" // TRIM(get_string_identifier(var_id)) // "'")

IF ( PRESENT(cube) .AND. PRESENT(const_val) )                                 &
  CALL log_warn("populate_var",                                               &
                "data cube and const_val both provided for variable '" //     &
                TRIM(get_string_identifier(var_id)) //                        &
                "' - using data in preference")

! We use the string identifier to select a CASE from the SELECT statement to
! avoid being dependent on how integer variable ids are implemented
SELECT CASE ( get_string_identifier(var_id) )
CASE ( 'latitude' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, latitude)
  ELSE
    latitude(:,:) = const_val
  END IF

CASE ( 'longitude' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, longitude)
  ELSE
    longitude(:,:) = const_val
  END IF

CASE ( 'land_fraction' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, flandg)
  ELSE
    flandg(:,:) = const_val
  END IF

CASE ( 'surf_hgt_surft' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, surf_hgt_surft)
  ELSE
    surf_hgt_surft(:,:) = const_val
  END IF
CASE ( 'z_land_land' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, z_land_land)
  ELSE
    z_land_land(:) = const_val
  END IF

CASE ( 'frac' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, frac_surft)
  ELSE
    frac_surft(:,:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Soil properties have a corresponding 0-level version for input of constant
  ! z values. albsoil_gb, albobs_sw_gb, etc have no levels anyway
  !
  ! Each variable can be read in with _soilt if the soil tile dimension is
  ! present. Where appropriate, this can be combined with _const_z as well.
  ! There are a few possible outcomes
  !
  ! No suffix:
  !  l_tile_soil = T and l_broadcast_ancils = T
  !  --Spread around soilt dim
  !  l_tile_soil = F
  !  --Do nothing
  !
  ! _const_z suffix:
  !  l_tile_soil = T and l_broadcast_ancils = T (and l_const_z = T)
  !  --Spread around soilt and layer dims
  !  l_tile_soil = F
  !  --Spread around layer dim
  !
  ! _soilt suffix:
  ! l_tile_soil = T and l_broadcast_ancils = F
  ! --Do nothing
  !
  ! _const_z_soilt:
  ! l_tile_soil = T and l_broadcast_ancils = F (and l_const_z = T)
  ! --Spread around layer dim
  !
  !-----------------------------------------------------------------------------

      !=== Start cases for albsoil ===
CASE ( 'albsoil' )  !no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, albsoil_soilt(:,1))
      albsoil_soilt(:,:) = SPREAD(albsoil_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, albsoil_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate albsoil- failure in logic")
    END IF
  ELSE
    albsoil_soilt(:,:) = const_val
  END IF

CASE ( 'albsoil_soilt' ) !with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, albsoil_soilt(:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate albsoil_soilt- failure in logic")
    END IF
  ELSE
    albsoil_soilt(:,:) = const_val
  END IF
  !=== End cases for albsoil ===

  !=== Start cases for clay ===
CASE ( 'clay' )  ! With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt dim.
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, clay_soilt(:,1,1:dim_cslayer))
      clay_soilt(:,:,:) = SPREAD(clay_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Do nothing special.
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, clay_soilt(:,1,1:dim_cslayer))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate clay - failure in logic")
    END IF
  ELSE
    clay_soilt(:,:,:) = const_val
  END IF

CASE ( 'clay_const_z' )  ! No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, clay_soilt(:,1,1))
      clay_soilt(:,1,1:dim_cslayer) = SPREAD(clay_soilt(:,1,1), 2,            &
                                             dim_cslayer)
      clay_soilt(:,:,:)             = SPREAD(clay_soilt(:,1,:), 2,            &
                                             nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, clay_soilt(:,1,1))
      clay_soilt(:,1,1:dim_cslayer) = SPREAD(clay_soilt(:,1,1), 2,            &
                                             dim_cslayer)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate clay_const_z - failure in logic")
    END IF
  ELSE
    clay_soilt(:,:,:) = const_val
  END IF

CASE ( 'clay_soilt' )  ! With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       ! Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, clay_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate clay_soilt- failure in logic")
    END IF
  ELSE
    clay_soilt(:,:,:) = const_val
  END IF

CASE ( 'clay_const_z_soilt' ) ! No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, clay_soilt(:,1:nsoilt,1))
      clay_soilt(:,:,:) = SPREAD(clay_soilt(:,:,1), 3, dim_cslayer)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate clay_const_z_soilt- failure in logic")
    END IF
  ELSE
    clay_soilt(:,:,:) = const_val
  END IF
  !=== End cases for clay ===

  !=== Start cases for b ===
CASE ( 'b' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, bexp_soilt(:,1,1:sm_levels))
      bexp_soilt(:,:,:) = SPREAD(bexp_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, bexp_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate b- failure in logic")
    END IF
  ELSE
    bexp_soilt(:,:,:) = const_val
  END IF

CASE ( 'b_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, bexp_soilt(:,1,1))
      bexp_soilt(:,1,1:sm_levels) = SPREAD(bexp_soilt(:,1,1), 2, sm_levels)
      bexp_soilt(:,:,:)           = SPREAD(bexp_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, bexp_soilt(:,1,1))
      bexp_soilt(:,1,1:sm_levels) = SPREAD(bexp_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate b_const_z- failure in logic")
    END IF
  ELSE
    bexp_soilt(:,:,:) = const_val
  END IF

CASE ( 'b_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, bexp_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate b_soilt- failure in logic")
    END IF
  ELSE
    bexp_soilt(:,:,:) = const_val
  END IF

CASE ( 'b_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, bexp_soilt(:,1:nsoilt,1))
      bexp_soilt(:,:,:) = SPREAD(bexp_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate b_const_z_soilt- failure in logic")
    END IF
  ELSE
    bexp_soilt(:,:,:) = const_val
  END IF
  !=== End cases for b ===

  !=== Start cases for sathh ===
CASE ( 'sathh' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sathh_soilt(:,1,1:sm_levels))
      sathh_soilt(:,:,:) = SPREAD(sathh_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sathh_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sathh- failure in logic")
    END IF
  ELSE
    sathh_soilt(:,:,:) = const_val
  END IF

CASE ( 'sathh_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sathh_soilt(:,1,1))
      sathh_soilt(:,1,1:sm_levels) = SPREAD(sathh_soilt(:,1,1), 2, sm_levels)
      sathh_soilt(:,:,:)           = SPREAD(sathh_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sathh_soilt(:,1,1))
      sathh_soilt(:,1,1:sm_levels) = SPREAD(sathh_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sathh_const_z- failure in logic")
    END IF
  ELSE
    sathh_soilt(:,:,:) = const_val
  END IF

CASE ( 'sathh_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sathh_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sathh_soilt- failure in logic")
    END IF
  ELSE
    sathh_soilt(:,:,:) = const_val
  END IF

CASE ( 'sathh_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sathh_soilt(:,1:nsoilt,1))
      sathh_soilt(:,:,:) = SPREAD(sathh_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sathh_const_z_soilt- failure in logic")
    END IF
  ELSE
    sathh_soilt(:,:,:) = const_val
  END IF
  !=== End cases for sathh ===

  !=== Start cases for sathh ===
CASE ( 'satcon' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, satcon_soilt(:,1,1:sm_levels))
      satcon_soilt(:,:,:) = SPREAD(satcon_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, satcon_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate satcon- failure in logic")
    END IF
  ELSE
    satcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'satcon_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, satcon_soilt(:,1,1))
      satcon_soilt(:,1,1:sm_levels) = SPREAD(satcon_soilt(:,1,1), 2, sm_levels)
      satcon_soilt(:,:,:)           = SPREAD(satcon_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, satcon_soilt(:,1,1))
      satcon_soilt(:,1,1:sm_levels) = SPREAD(satcon_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate satcon_const_z- failure in logic")
    END IF
  ELSE
    satcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'satcon_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, satcon_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate satcon_soilt- failure in logic")
    END IF
  ELSE
    satcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'satcon_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, satcon_soilt(:,1:nsoilt,1))
      satcon_soilt(:,:,:) = SPREAD(satcon_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate satcon_const_z_soilt- failure in logic")
    END IF
  ELSE
    satcon_soilt(:,:,:) = const_val
  END IF
  !=== End cases for satcon ===

  !=== Start cases for sm_sat ===
CASE ( 'sm_sat' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcst_soilt(:,1,1:sm_levels))
      smvcst_soilt(:,:,:) = SPREAD(smvcst_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcst_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_sat- failure in logic")
    END IF
  ELSE
    smvcst_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_sat_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcst_soilt(:,1,1))
      smvcst_soilt(:,1,1:sm_levels) = SPREAD(smvcst_soilt(:,1,1), 2, sm_levels)
      smvcst_soilt(:,:,:)           = SPREAD(smvcst_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcst_soilt(:,1,1))
      smvcst_soilt(:,1,1:sm_levels) = SPREAD(smvcst_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_sat_const_z- failure in logic")
    END IF
  ELSE
    smvcst_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_sat_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcst_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_sat_soilt- failure in logic")
    END IF
  ELSE
    smvcst_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_sat_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcst_soilt(:,1:nsoilt,1))
      smvcst_soilt(:,:,:) = SPREAD(smvcst_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_sat_const_z_soilt- failure in logic")
    END IF
  ELSE
    smvcst_soilt(:,:,:) = const_val
  END IF
  !=== End cases for sm_sat ===

  !=== Start cases for sm_crit ===
CASE ( 'sm_crit' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvccl_soilt(:,1,1:sm_levels))
      smvccl_soilt(:,:,:) = SPREAD(smvccl_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvccl_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_crit- failure in logic")
    END IF
  ELSE
    smvccl_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_crit_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvccl_soilt(:,1,1))
      smvccl_soilt(:,1,1:sm_levels) = SPREAD(smvccl_soilt(:,1,1), 2, sm_levels)
      smvccl_soilt(:,:,:)           = SPREAD(smvccl_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvccl_soilt(:,1,1))
      smvccl_soilt(:,1,1:sm_levels) = SPREAD(smvccl_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_crit_const_z- failure in logic")
    END IF
  ELSE
    smvccl_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_crit_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvccl_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_crit_soilt- failure in logic")
    END IF
  ELSE
    smvccl_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_crit_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvccl_soilt(:,1:nsoilt,1))
      smvccl_soilt(:,:,:) = SPREAD(smvccl_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_crit_const_z_soilt- failure in logic")
    END IF
  ELSE
    smvccl_soilt(:,:,:) = const_val
  END IF
  !=== End cases for sm_crit ===

  !=== Start cases for sm_wilt ===
CASE ( 'sm_wilt' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcwt_soilt(:,1,1:sm_levels))
      smvcwt_soilt(:,:,:) = SPREAD(smvcwt_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcwt_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_wilt- failure in logic")
    END IF
  ELSE
    smvcwt_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_wilt_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcwt_soilt(:,1,1))
      smvcwt_soilt(:,1,1:sm_levels) = SPREAD(smvcwt_soilt(:,1,1), 2, sm_levels)
      smvcwt_soilt(:,:,:)           = SPREAD(smvcwt_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcwt_soilt(:,1,1))
      smvcwt_soilt(:,1,1:sm_levels) = SPREAD(smvcwt_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_wilt_const_z- failure in logic")
    END IF
  ELSE
    smvcwt_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_wilt_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcwt_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_wilt_soilt- failure in logic")
    END IF
  ELSE
    smvcwt_soilt(:,:,:) = const_val
  END IF

CASE ( 'sm_wilt_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, smvcwt_soilt(:,1:nsoilt,1))
      smvcwt_soilt(:,:,:) = SPREAD(smvcwt_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sm_wilt_const_z_soilt- failure in logic")
    END IF
  ELSE
    smvcwt_soilt(:,:,:) = const_val
  END IF
  !=== End cases for sm_wilt ===

  !=== Start cases for hcap ===
CASE ( 'hcap' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcap_soilt(:,1,1:sm_levels))
      hcap_soilt(:,:,:) = SPREAD(hcap_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcap_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcap- failure in logic")
    END IF
  ELSE
    hcap_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcap_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcap_soilt(:,1,1))
      hcap_soilt(:,1,1:sm_levels) = SPREAD(hcap_soilt(:,1,1), 2, sm_levels)
      hcap_soilt(:,:,:)           = SPREAD(hcap_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcap_soilt(:,1,1))
      hcap_soilt(:,1,1:sm_levels) = SPREAD(hcap_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcap_const_z- failure in logic")
    END IF
  ELSE
    hcap_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcap_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcap_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcap_soilt- failure in logic")
    END IF
  ELSE
    hcap_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcap_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcap_soilt(:,1:nsoilt,1))
      hcap_soilt(:,:,:) = SPREAD(hcap_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcap_const_z_soilt- failure in logic")
    END IF
  ELSE
    hcap_soilt(:,:,:) = const_val
  END IF
  !=== End cases for hcap ===

  !=== Start cases for hcon ===
CASE ( 'hcon' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcon_soilt(:,1,1:sm_levels))
      hcon_soilt(:,:,:) = SPREAD(hcon_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcon_soilt(:,1,1:sm_levels))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcon- failure in logic")
    END IF
  ELSE
    hcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcon_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcon_soilt(:,1,1))
      hcon_soilt(:,1,1:sm_levels) = SPREAD(hcon_soilt(:,1,1), 2, sm_levels)
      hcon_soilt(:,:,:)           = SPREAD(hcon_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcon_soilt(:,1,1))
      hcon_soilt(:,1,1:sm_levels) = SPREAD(hcon_soilt(:,1,1), 2, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcon_const_z- failure in logic")
    END IF
  ELSE
    hcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcon_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcon_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcon_soilt- failure in logic")
    END IF
  ELSE
    hcon_soilt(:,:,:) = const_val
  END IF

CASE ( 'hcon_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, hcon_soilt(:,1:nsoilt,1))
      hcon_soilt(:,:,:) = SPREAD(hcon_soilt(:,:,1), 3, sm_levels)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate hcon_const_z_soilt- failure in logic")
    END IF
  ELSE
    hcon_soilt(:,:,:) = const_val
  END IF
  !=== End cases for hcon ===

  !-----------------------------------------------------------------------------
  ! Soil ancillaries on soil carbon layers.
  !-----------------------------------------------------------------------------

      !=== Start cases for soil_ph ===
CASE ( 'soil_ph' )  !With levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, soil_ph_soilt(:,1,1:dim_cslayer))
      soil_ph_soilt(:,:,:) = SPREAD(soil_ph_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, soil_ph_soilt(:,1,1:dim_cslayer))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate soil_ph- failure in logic")
    END IF
  ELSE
    soil_ph_soilt(:,:,:) = const_val
  END IF

CASE ( 'soil_ph_const_z' ) !No levels, no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      ! Spread around soilt and layer dims
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, soil_ph_soilt(:,1,1))
      soil_ph_soilt(:,1,1:dim_cslayer) = SPREAD(soil_ph_soilt(:,1,1),         &
                                                  2, dim_cslayer)
      soil_ph_soilt(:,:,:)           = SPREAD(soil_ph_soilt(:,1,:),           &
                                                2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      ! Copy the values from the first vertical level to all others
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, soil_ph_soilt(:,1,1))
      soil_ph_soilt(:,1,1:dim_cslayer) = SPREAD(soil_ph_soilt(:,1,1),         &
                                                  2, dim_cslayer)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate soil_ph_const_z- failure in logic")
    END IF
  ELSE
    soil_ph_soilt(:,:,:) = const_val
  END IF

CASE ( 'soil_ph_soilt' ) !With levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, soil_ph_soilt(:,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate soil_ph_soilt- failure in logic")
    END IF
  ELSE
    soil_ph_soilt(:,:,:) = const_val
  END IF

CASE ( 'soil_ph_const_z_soilt' ) !No levels, with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Spread around layer dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, soil_ph_soilt(:,1:nsoilt,1))
      soil_ph_soilt(:,:,:) = SPREAD(soil_ph_soilt(:,:,1), 3,                  &
                                      dim_cslayer)
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate soil_ph_const_z_soilt- " //             &
                     "failure in logic")
    END IF
  ELSE
    soil_ph_soilt(:,:,:) = const_val
  END IF
  !=== End cases for soil_ph ===

  !-----------------------------------------------------------------------------
  ! TOPMODEL variables are just on land points
  !
  ! Each variable can be read in with _soilt if the soil tile dimension is
  ! present.
  ! There are a few possible outcomes:
  !
  ! No suffix:
  !  l_tile_soil = T and l_broadcast_ancils = T
  !  --Spread around soilt dim
  !  l_tile_soil = F
  !  --Do nothing
  !
  ! _soilt suffix:
  ! l_tile_soil = T and l_broadcast_ancils = F
  ! --Do nothing
  !-----------------------------------------------------------------------------
      !=== Start cases for fexp ===
CASE ( 'fexp' )  !no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, fexp_soilt(:,1))
      fexp_soilt(:,:) = SPREAD(fexp_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, fexp_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate fexp- failure in logic")
    END IF
  ELSE
    fexp_soilt(:,:) = const_val
  END IF

CASE ( 'fexp_soilt' ) !with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, fexp_soilt(:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate fexp_soilt- failure in logic")
    END IF
  ELSE
    fexp_soilt(:,:) = const_val
  END IF
  !=== End cases for fexp ===

  !=== Start cases for ti_mean ===
CASE ( 'ti_mean' )  !no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, ti_mean_soilt(:,1))
      ti_mean_soilt(:,:) = SPREAD(ti_mean_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, ti_mean_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate ti_mean- failure in logic")
    END IF
  ELSE
    ti_mean_soilt(:,:) = const_val
  END IF

CASE ( 'ti_mean_soilt' ) !with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, ti_mean_soilt(:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate ti_mean_soilt- failure in logic")
    END IF
  ELSE
    ti_mean_soilt(:,:) = const_val
  END IF
  !=== End cases for ti_mean ===

  !=== Start cases for ti_sig ===
CASE ( 'ti_sig' )  !no tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_ancils) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, ti_sig_soilt(:,1))
      ti_sig_soilt(:,:) = SPREAD(ti_sig_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, ti_sig_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate ti_sig- failure in logic")
    END IF
  ELSE
    ti_sig_soilt(:,:) = const_val
  END IF

CASE ( 'ti_sig_soilt' ) !with tiles
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. .NOT. l_broadcast_ancils) THEN
       !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, ti_sig_soilt(:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate ti_sig_soilt- failure in logic")
    END IF
  ELSE
    ti_sig_soilt(:,:) = const_val
  END IF
  !=== End cases for ti_sig ===

  !-----------------------------------------------------------------------------
  ! Agricultural fraction is on land points only
  !-----------------------------------------------------------------------------
CASE ( 'frac_agr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, frac_agr_gb)
  ELSE
    frac_agr_gb(:) = const_val
  END IF
CASE ( 'frac_past' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, frac_past_gb)
  ELSE
    frac_past_gb(:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Irrigation fraction is on land points only
  !-----------------------------------------------------------------------------
CASE ( 'frac_irrig' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, frac_irr_all(:,1))
  ELSE
    frac_irr_all(:,1) = const_val
  END IF

CASE ( 'frac_agr_prev' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, frac_agr_prev_gb)
  ELSE
    frac_agr_prev_gb(:) = const_val
  END IF

CASE ( 'frac_past_prev' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, frac_past_prev_gb)
  ELSE
    frac_past_prev_gb(:) = const_val
  END IF

CASE ( 'wood_prod_fast' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, wood_prod_fast_gb)
  ELSE
    wood_prod_fast_gb(:) = const_val
  END IF

CASE ( 'wood_prod_med' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, wood_prod_med_gb)
  ELSE
    wood_prod_med_gb(:) = const_val
  END IF

CASE ( 'wood_prod_slow' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, wood_prod_slow_gb)
  ELSE
    wood_prod_slow_gb(:) = const_val
  END IF

CASE ( 'n_inorg' )
  IF ( nsoilt == 1 ) THEN
    IF ( PRESENT(cube) ) THEN
      cube_land = map_to_land(cube)
      CALL cube_get_data( cube_land, n_inorg_soilt_lyrs(:,1,:) )
    ELSE
      n_inorg_soilt_lyrs(:,:,:) = const_val
    END IF
  ELSE
    CALL log_fatal("populate_var",                                            &
                   "Cannot populate n_inorg when nsoilt > 1.")
  END IF

  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, n_inorg_soilt_lyrs(:,1,:))
      n_inorg_soilt_lyrs(:,:,:) =                                             &
        SPREAD(n_inorg_soilt_lyrs(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, n_inorg_soilt_lyrs(:,1,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate n_inorg- failure in logic")
    END IF
  ELSE
    n_inorg_soilt_lyrs(:,:,:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Urban variables are on land points only
  !-----------------------------------------------------------------------------
CASE ( 'wrr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, wrr_gb)
  ELSE
    wrr_gb(:) = const_val
  END IF

CASE ( 'hwr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, hwr_gb)
  ELSE
    hwr_gb(:) = const_val
  END IF

CASE ( 'hgt' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, hgt_gb)
  ELSE
    hgt_gb(:) = const_val
  END IF

CASE ( 'ztm' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, ztm_gb)
  ELSE
    ztm_gb(:) = const_val
  END IF

CASE ( 'disp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, disp_gb)
  ELSE
    disp_gb(:) = const_val
  END IF

CASE ( 'albwl' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, albwl_gb)
  ELSE
    albwl_gb(:) = const_val
  END IF

CASE ( 'albrd' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, albrd_gb)
  ELSE
    albrd_gb(:) = const_val
  END IF

CASE ( 'emisw' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, emisw_gb)
  ELSE
    emisw_gb(:) = const_val
  END IF

CASE ( 'emisr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, emisr_gb)
  ELSE
    emisr_gb(:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Variables that are set as initial conditions
  ! These are all land points only, but with varying numbers of levels

  ! Some variables can be read in with _soilt if the soil tile dimension is
  ! present. There are a few possible outcomes
  !
  ! No suffix:
  !  l_tile_soil = T and l_broadcast_soilt = T
  !  --Spread around soilt dim
  !  l_tile_soil = F
  !  --Do nothing
  !
  ! _soilt suffix:
  ! l_tile_soil = T and l_broadcast_soilt = F
  ! --Do nothing
  !-----------------------------------------------------------------------------

CASE ( 'canopy' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, canopy_surft)
  ELSE
    canopy_surft(:,:) = const_val
  END IF

CASE ( 'cs' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, cs_pool_soilt(:,1,:,:))
      cs_pool_soilt(:,:,:,:) = SPREAD(cs_pool_soilt(:,1,:,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, cs_pool_soilt(:,1,:,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate cs- failure in logic")
    END IF
  ELSE
    cs_pool_soilt(:,:,:,:) = const_val
  END IF

CASE ( 'ns' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, ns_pool_gb)
  ELSE
    ns_pool_gb(:,:,:) = const_val
  END IF

CASE ( 'deposition_n' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, deposition_n_gb)
  ELSE
    deposition_n_gb(:) = const_val
  END IF

CASE ( 'g_burn_pft' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, g_burn_pft)
  ELSE
    g_burn_pft(:,:) = const_val
  END IF

CASE ( 'gs' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, gs_gb)
  ELSE
    gs_gb(:) = const_val
  END IF

CASE ( 'snow_tile' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, snow_surft)
  ELSE
    snow_surft(:,:) = const_val
  END IF

CASE ( 'sthuf' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sthuf_soilt(:,1,:))
      sthuf_soilt(:,:,:) = SPREAD(sthuf_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sthuf_soilt(:,1,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sthuf- failure in logic")
    END IF
  ELSE
    sthuf_soilt(:,:,:) = const_val
  END IF

CASE ( 't_soil' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, t_soil_soilt(:,1,:))
      t_soil_soilt(:,:,:) = SPREAD(t_soil_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, t_soil_soilt(:,1,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate t_soil- failure in logic")
    END IF
  ELSE
    t_soil_soilt(:,:,:) = const_val
  END IF

CASE ( 'tsoil_deep' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, tsoil_deep_gb)
  ELSE
    tsoil_deep_gb(:,:) = const_val
  END IF

CASE ( 'sthu_irr' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sthu_irr_soilt(:,1,:))
      sthu_irr_soilt(:,:,:) = SPREAD(sthu_irr_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sthu_irr_soilt(:,1,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sthu_irr- failure in logic")
    END IF
  ELSE
    sthu_irr_soilt(:,:,:) = const_val
  END IF

CASE ( 'tstar_tile' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, tstar_surft)
  ELSE
    tstar_surft(:,:) = const_val
  END IF

CASE ( 'tsurf_elev_surft' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, tsurf_elev_surft)
  ELSE
    tsurf_elev_surft(:,:) = const_val
  END IF

CASE ( 'cropdvi' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, dvi_cpft)
  ELSE
    dvi_cpft(:,:) = const_val
  END IF

CASE ( 'croprootc' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, rootc_cpft)
  ELSE
    rootc_cpft(:,:) = const_val
  END IF

CASE ( 'cropharvc' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, harvc_cpft)
  ELSE
    harvc_cpft(:,:) = const_val
  END IF

CASE ( 'cropreservec' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, reservec_cpft)
  ELSE
    reservec_cpft(:,:) = const_val
  END IF

CASE ( 'croplai' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, croplai_cpft)
  ELSE
    croplai_cpft(:,:) = const_val
  END IF

CASE ( 'cropcanht' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, cropcanht_cpft)
  ELSE
    cropcanht_cpft(:,:) = const_val
  END IF

CASE ( 'cropsowdate' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, sow_date_cpft)
  ELSE
    sow_date_cpft(:,:) = const_val
  END IF

CASE ( 'croplatestharvdate' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, latestharv_date_cpft)
  ELSE
    latestharv_date_cpft(:,:) = const_val
  END IF

CASE ( 'cropttveg' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, tt_veg_cpft)
  ELSE
    tt_veg_cpft(:,:) = const_val
  END IF

CASE ( 'cropttrep' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, tt_rep_cpft)
  ELSE
    tt_rep_cpft(:,:) = const_val
  END IF

CASE ( 'lai' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, lai_pft)
  ELSE
    lai_pft(:,:) = const_val
  END IF

CASE ( 'canht' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, canht_pft)
  ELSE
    canht_pft(:,:) = const_val
  END IF

CASE ( 'sthzw' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sthzw_soilt(:,1))
      sthzw_soilt(:,:) = SPREAD(sthzw_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, sthzw_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate sthzw- failure in logic")
    END IF
  ELSE
    sthzw_soilt(:,:) = const_val
  END IF

CASE ( 'zw' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, zw_soilt(:,1))
      zw_soilt(:,:) = SPREAD(zw_soilt(:,1), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, zw_soilt(:,1))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate zw- failure in logic")
    END IF
  ELSE
    zw_soilt(:,:) = const_val
  END IF

CASE ( 'slope' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, slope_gb)
  ELSE
    slope_gb(:) = const_val
  END IF

CASE ( 'rgrain' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, rgrain_surft)
  ELSE
    rgrain_surft(:,:) = const_val
  END IF

CASE ( 'cv' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, cv_gb)
  ELSE
    cv_gb(:) = const_val
  END IF

CASE ( 'rho_snow' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, rho_snow_grnd_surft)
  ELSE
    rho_snow_grnd_surft(:,:) = const_val
  END IF

CASE ( 'snow_depth' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, snowdepth_surft)
  ELSE
    snowdepth_surft(:,:) = const_val
  END IF

CASE ( 'snow_grnd' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, snow_grnd_surft)
  ELSE
    snow_grnd_surft(:,:) = const_val
  END IF

CASE ( 'nsnow' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, nsnow_real)
    nsnow_surft(:,:) = NINT(nsnow_real)
  ELSE
    nsnow_surft(:,:) = NINT(const_val)
  END IF

CASE ( 'snow_ds' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, ds_surft)
  ELSE
    ds_surft(:,:,:) = const_val
  END IF

CASE ( 'snow_ice' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, sice_surft)
  ELSE
    sice_surft(:,:,:) = const_val
  END IF

CASE ( 'snow_liq' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, sliq_surft)
  ELSE
    sliq_surft(:,:,:) = const_val
  END IF

CASE ( 'tsnow' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, tsnow_surft)
  ELSE
    tsnow_surft(:,:,:) = const_val
  END IF

CASE ( 'rgrainl' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, rgrainl_surft)
  ELSE
    rgrainl_surft(:,:,:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Forcing variables
  !-----------------------------------------------------------------------------
CASE ( 'pstar' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, pstar_ij)
  ELSE
    pstar_ij(:,:) = const_val
  END IF

CASE ( 'q' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, qw_1_ij)
  ELSE
    qw_1_ij(:,:) = const_val
  END IF
  qw_1_ij(:,:) = MAX(qw_1_ij, 0.0)

CASE ( 't' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, tl_1_ij)
  ELSE
    tl_1_ij(:,:) = const_val
  END IF

CASE ( 'rad_net' )
  ! Net downward radiation is stored in lw_down until it is processed
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, lw_down_ij)
  ELSE
    lw_down_ij(:,:) = const_val
  END IF

CASE ( 'lw_net' )
  ! Net LW downward radiation is stored in lw_down until it is processed
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, lw_down_ij)
  ELSE
    lw_down_ij(:,:) = const_val
  END IF

CASE ( 'sw_net' )
  ! Net SW downward radiation is stored in sw_down until it is processed
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, sw_down_ij)
  ELSE
    sw_down_ij(:,:) = const_val
  END IF

CASE ( 'lw_down' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, lw_down_ij)
  ELSE
    lw_down_ij(:,:) = const_val
  END IF

CASE ( 'sw_down' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, sw_down_ij)
  ELSE
    sw_down_ij(:,:) = const_val
  END IF

CASE ( 'diff_rad' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, diff_rad_ij)
  ELSE
    diff_rad_ij(:,:) = const_val
  END IF

CASE ( 'dt_range' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, diurnal_temperature_range_ij)
  ELSE
    diurnal_temperature_range_ij(:,:) = const_val
  END IF

CASE ( 'precip' )
  ! Store total precip as large-scale rainfall until it is partitioned
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, ls_rain_ij)
  ELSE
    ls_rain_ij(:,:) = const_val
  END IF

CASE ( 'tot_rain' )
  ! Store total rainfall as large-scale until it is partitioned
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, ls_rain_ij)
  ELSE
    ls_rain_ij(:,:) = const_val
  END IF

CASE ( 'tot_snow' )
  ! If given total snow, we assume it is all large-scale
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, ls_snow_ij)
  ELSE
    ls_snow_ij(:,:) = const_val
  END IF

CASE ( 'con_rain' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, con_rain_ij)
  ELSE
    con_rain_ij(:,:) = const_val
  END IF

CASE ( 'ls_rain' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, ls_rain_ij)
  ELSE
    ls_rain_ij(:,:) = const_val
  END IF

CASE ( 'con_snow' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, con_snow_ij)
  ELSE
    con_snow_ij(:,:) = const_val
  END IF

CASE ( 'ls_snow' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, ls_snow_ij)
  ELSE
    ls_snow_ij(:,:) = const_val
  END IF

CASE ( 'wind' )
  ! Wind speed just goes directly into u component
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, u_1_ij)
  ELSE
    u_1_ij(:,:) = const_val
  END IF

CASE ( 'u' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, u_1_ij)
  ELSE
    u_1_ij(:,:) = const_val
  END IF

CASE ( 'v' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, v_1_ij)
  ELSE
    v_1_ij(:,:) = const_val
  END IF

CASE ( 'z1_tq_in' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, z1_tq_ij)
  ELSE
    z1_tq_ij(:,:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! River routing variables
  !-----------------------------------------------------------------------------

CASE ( 'latitude_2d' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, rivers_lat2d)
  ELSE
    rivers_lat2d(:,:) = const_val
  END IF

CASE ( 'longitude_2d' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, rivers_lon2d)
  ELSE
    rivers_lon2d(:,:) = const_val
  END IF

CASE ( 'direction' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, rivers_dir)
  ELSE
    rivers_dir(:,:) = const_val
  END IF

CASE ( 'sequence' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, rivers_seq)
  ELSE
    rivers_seq(:,:) = const_val
  END IF

CASE ( 'area' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, rivers_dra)
  ELSE
    rivers_dra(:,:) = const_val
  END IF

CASE ( 'rivers_ygrid' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    cube_land%values(:) = RESHAPE(cube%values, (/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rivers_ygrid)
  ELSE
    rivers_ygrid(:) = const_val
  END IF

CASE ( 'rivers_xgrid' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    cube_land%values(:) = RESHAPE(cube%values, (/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rivers_xgrid)
  ELSE
    rivers_xgrid(:) = const_val
  END IF

CASE ( 'rivers_sto_rp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rivers_sto_rp)
  ELSE
    rivers_sto_rp(:) = const_val
  END IF

CASE ( 'rfm_surfstore_rp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rfm_surfstore_rp)
  ELSE
    rfm_surfstore_rp(:) = const_val
  END IF

CASE ( 'rfm_substore_rp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rfm_substore_rp)
  ELSE
    rfm_substore_rp(:) = const_val
  END IF

CASE ( 'rfm_flowin_rp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rfm_flowin_rp)
  ELSE
    rfm_flowin_rp(:) = const_val
  END IF

CASE ( 'rfm_bflowin_rp' )
  IF ( PRESENT(cube) ) THEN
    cube_land = cube_create((/ SIZE(cube%values) /))
    CALL cube_get_data(cube_land, rfm_bflowin_rp)
  ELSE
    rfm_bflowin_rp(:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Overbank inundation variables
  !-----------------------------------------------------------------------------

CASE ( 'logn_mean' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, logn_mean)
  ELSE
    logn_mean(:,:) = const_val
  END IF

CASE ( 'logn_stdev' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, logn_stdev)
  ELSE
    logn_stdev(:,:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Fire and metstats variables- land points only
  !-----------------------------------------------------------------------------
  ! metstats module cases- land points only
CASE ( 'temp_max_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%temp_max_00h%run)
  ELSE
    metstats_prog(:)%temp_max_00h%run = const_val
  END IF

CASE ( 'temp_ave_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%temp_ave_00h%run)
  ELSE
    metstats_prog(:)%temp_ave_00h%run = const_val
  END IF

CASE ( 'prec_tot_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%prec_tot_00h%run)
  ELSE
    metstats_prog(:)%prec_tot_00h%run = const_val
  END IF

CASE ( 'prec_tot_12h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%prec_tot_12h%run)
  ELSE
    metstats_prog(:)%prec_tot_12h%run = const_val
  END IF

CASE ( 'rhum_min_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%rhum_min_00h%run)
  ELSE
    metstats_prog(:)%rhum_min_00h%run = const_val
  END IF

CASE ( 'dewp_ave_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%dewp_ave_00h%run)
  ELSE
    metstats_prog(:)%dewp_ave_00h%run = const_val
  END IF

CASE ( 'wind_ave_00h_r' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%wind_ave_00h%run)
  ELSE
    metstats_prog(:)%wind_ave_00h%run = const_val
  END IF

CASE ( 'temp_max_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%temp_max_00h%fin)
  ELSE
    metstats_prog(:)%temp_max_00h%fin = const_val
  END IF

CASE ( 'temp_ave_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%temp_ave_00h%fin)
  ELSE
    metstats_prog(:)%temp_ave_00h%fin = const_val
  END IF

CASE ( 'temp_pnt_12h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%temp_pnt_12h%fin)
  ELSE
    metstats_prog(:)%temp_pnt_12h%fin = const_val
  END IF

CASE ( 'prec_tot_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%prec_tot_00h%fin)
  ELSE
    metstats_prog(:)%prec_tot_00h%fin = const_val
  END IF

CASE ( 'prec_tot_12h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%prec_tot_12h%fin)
  ELSE
    metstats_prog(:)%prec_tot_12h%fin = const_val
  END IF

CASE ( 'rhum_min_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%rhum_min_00h%fin)
  ELSE
    metstats_prog(:)%rhum_min_00h%fin = const_val
  END IF

CASE ( 'rhum_pnt_12h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%rhum_pnt_12h%fin)
  ELSE
    metstats_prog(:)%rhum_pnt_12h%fin = const_val
  END IF

CASE ( 'dewp_ave_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%dewp_ave_00h%fin)
  ELSE
    metstats_prog(:)%dewp_ave_00h%fin = const_val
  END IF

CASE ( 'wind_ave_00h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%wind_ave_00h%fin)
  ELSE
    metstats_prog(:)%wind_ave_00h%fin = const_val
  END IF

CASE ( 'wind_pnt_12h' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, metstats_prog(:)%wind_pnt_12h%fin)
  ELSE
    metstats_prog(:)%wind_pnt_12h%fin = const_val
  END IF

  ! Fire module variables- land points only
CASE ( 'fire_mcarthur_r_dr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%mcarthur%r_dr)
  ELSE
    fire_prog(:)%mcarthur%r_dr = const_val
  END IF

CASE ( 'fire_mcarthur_n_dr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%mcarthur%n_dr)
  ELSE
    fire_prog(:)%mcarthur%n_dr = const_val
  END IF

CASE ( 'fire_canadian_ffmc' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%canadian%ffmc)
  ELSE
    fire_prog(:)%canadian%ffmc = const_val
  END IF

CASE ( 'fire_canadian_ffmc_mois' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%canadian%ffmc_mois)
  ELSE
    fire_prog(:)%canadian%ffmc_mois = const_val
  END IF

CASE ( 'fire_canadian_dmc' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%canadian%dmc)
  ELSE
    fire_prog(:)%canadian%dmc = const_val
  END IF

CASE ( 'fire_canadian_dc' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%canadian%dc)
  ELSE
    fire_prog(:)%canadian%dc = const_val
  END IF

CASE ( 'fire_nesterov' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, fire_prog(:)%nesterov%findex)
  ELSE
    fire_prog(:)%nesterov%findex = const_val
  END IF

  !-----------------------------------------------------------------------------
  ! Other variables that might be prescribed
  !-----------------------------------------------------------------------------
CASE ( 'bl_height' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, zh)
  ELSE
    zh(:,:) = const_val
  END IF

CASE ( 'level_separation' )
  IF ( PRESENT(cube) ) THEN
    CALL cube_get_data(cube, dzl)
  ELSE
    dzl(:,:,:) = const_val
  END IF

CASE ( 'ozone' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, o3_gb)
  ELSE
    o3_gb(:) = const_val
  END IF

CASE ( 'tracer_field' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, tracer_field)
  ELSE
    tracer_field(:,:) = const_val
  END IF

CASE ( 'albobs_sw' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, albobs_sw_gb)
  ELSE
    albobs_sw_gb(:) = const_val
  END IF

CASE ( 'albobs_vis' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, albobs_vis_gb)
  ELSE
    albobs_vis_gb(:) = const_val
  END IF

CASE ( 'albobs_nir' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, albobs_nir_gb)
  ELSE
    albobs_nir_gb(:) = const_val
  END IF

CASE ( 'co2_mmr' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, co2_tmp)

    IF ((MAXVAL(co2_tmp) /= co2_tmp(1)) .OR.                                  &
        (MINVAL(co2_tmp) /= co2_tmp(1))) THEN
      CALL log_fatal("populate_var",                                          &
                     "All land points must have same CO2 concetration")
    END IF
    co2_mmr = co2_tmp(1)
  ELSE
    co2_mmr = const_val
  END IF

CASE ( 'pop_den' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, pop_den)
  ELSE
    pop_den(:) = const_val
  END IF

CASE ( 'flash_rate' )
  IF ( PRESENT(cube) ) THEN
    cube_land = map_to_land(cube)
    CALL cube_get_data(cube_land, flash_rate)
  ELSE
    flash_rate(:) = const_val
  END IF

  !-----------------------------------------------------------------------------
  !   ECOSSE variables.
  !-----------------------------------------------------------------------------
  !   ECOSSE prognostic variables.
CASE ( 'n_amm' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, n_amm_soilt(:,1,:))
      n_amm_soilt(:,:,:) = SPREAD(n_amm_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, n_amm_soilt(:,1,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate n_amm- failure in logic")
    END IF
  ELSE
    n_amm_soilt(:,:,:) = const_val
  END IF

CASE ( 'n_nit' )
  IF ( PRESENT(cube) ) THEN
    IF ( l_tile_soil .AND. l_broadcast_soilt) THEN
      !Spread around soilt dim
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, n_nit_soilt(:,1,:))
      n_nit_soilt(:,:,:) = SPREAD(n_nit_soilt(:,1,:), 2, nsoilt)
    ELSE IF ( .NOT. l_tile_soil ) THEN
      !Do nothing special
      cube_land = map_to_land(cube)
      CALL cube_get_data(cube_land, n_nit_soilt(:,1,:))
    ELSE
      CALL log_fatal("populate_var",                                          &
                     "Cannot populate n_nit- failure in logic")
    END IF
  ELSE
    n_nit_soilt(:,:,:) = const_val
  END IF

    CASE ( 'SoilTemp_CABLE' )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, SoilTemp_CABLE )
      ELSE
        SoilTemp_CABLE(:,:,:) = const_val
      END IF

    CASE ( 'FrozenSoilFrac_CABLE' )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, FrozenSoilFrac_CABLE )
      ELSE
        FrozenSoilFrac_CABLE(:,:,:) = const_val
      END IF

    CASE ( 'SoilMoisture_CABLE' )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, SoilMoisture_CABLE )
      ELSE
        SoilMoisture_CABLE(:,:,:) = const_val
      END IF

    CASE ( 'SnowDepth_CABLE' )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, SnowDepth_CABLE )
      ELSE
        SnowDepth_CABLE(:,:,:) = const_val
      END IF

    CASE ( 'SnowAge_CABLE' )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, SnowAge_CABLE )
      ELSE
        SnowAge_CABLE(:,:) = const_val
      END IF
    
    CASE ( 'SnowMass_CABLE' )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, SnowMass_CABLE )
      ELSE
        SnowMass_CABLE(:,:,:) = const_val
      END IF
    
    CASE ( 'ThreeLayerSnowFlag_CABLE' )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, ThreeLayerSnowFlag_CABLE)
      ELSE
        ThreeLayerSnowFlag_CABLE(:,:) = const_val
      END IF
    
    CASE ( 'SnowTemp_CABLE' )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, SnowTemp_CABLE )
      ELSE
        SnowTemp_CABLE(:,:,:) = const_val
      END IF
    
    CASE ( 'SnowDensity_CABLE' )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, SnowDensity_CABLE)
      ELSE
        SnowDensity_CABLE(:,:,:) = const_val
      END IF
    
    CASE ( 'OneLyrSnowDensity_CABLE' )
      IF ( PRESENT(cube) ) THEN
        cube_land = map_to_land(cube)
        CALL cube_get_data(cube_land, OneLyrSnowDensity_CABLE)
      ELSE
        OneLyrSnowDensity_CABLE(:,:) = const_val
      END IF

CASE DEFAULT
  CALL log_fatal("populate_var",                                              &
                 "Unrecognised variable for input - '" //                     &
                 TRIM(get_string_identifier(var_id)) // "'. " //              &
                 "See docs for available variables")
END SELECT


! Free the land cube
CALL cube_free(cube_land)

RETURN

END SUBROUTINE populate_var
#endif
