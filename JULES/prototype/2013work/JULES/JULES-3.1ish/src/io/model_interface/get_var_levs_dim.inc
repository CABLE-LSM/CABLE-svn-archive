#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/model_interface/get_var_levs_dim.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-09-27 09:33:34 +0100 (Thu, 27 Sep 2012) $
!!   $LastChangedRevision: 609 $
!!
!!****************************************************************************

SUBROUTINE get_var_levs_dim(var_id, dim_name_in, dim_name_out, dim_size)

  USE string_utils_mod, ONLY : to_string

  USE max_dimensions, ONLY : ntiles_max

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an identifer for a model variable, returns the constant that
!   identifies what levels dimension to use
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  INTEGER, INTENT(IN) :: var_id  ! Identifies the model variable
  CHARACTER(len=*), INTENT(OUT) :: dim_name_in  ! The dimension name to use for
                                                ! the variable in input files
  CHARACTER(len=*), INTENT(OUT) :: dim_name_out  ! The dimension name to use for
                                                 ! the variable in output files
  INTEGER, INTENT(OUT) :: dim_size  ! The dimension size to use for the variable

!-----------------------------------------------------------------------------

  SELECT CASE ( var_id )
    CASE ( var_id_latitude, var_id_longitude, var_id_land_fraction )
! No vertical dimension required, so return the empty string as dimension name
      dim_name_in  = ''
      dim_name_out = ''
      dim_size     = surface_dim_size

    CASE ( var_id_surf_hgt )
      dim_name_in  = tile_dim_name
      dim_name_out = TILE_DIM_NAME_OUT
      dim_size     = tile_dim_size

    CASE ( var_id_frac )
      dim_name_in  = type_dim_name
      dim_name_out = TYPE_DIM_NAME_OUT
      dim_size     = type_dim_size

!-----------------------------------------------------------------------------
! Soil variables
!-----------------------------------------------------------------------------
    CASE ( var_id_b, var_id_sathh, var_id_satcon, var_id_sm_sat,              &
           var_id_sm_crit, var_id_sm_wilt, var_id_hcap, var_id_hcon )
      dim_name_in  = soil_dim_name
      dim_name_out = SOIL_DIM_NAME_OUT
      dim_size     = soil_dim_size

! Soil properties have a corresponding 0-level version for input of constant
! z values. albsoil has no levels anyway
    CASE ( var_id_albsoil, var_id_b_const_z, var_id_sathh_const_z,            &
           var_id_satcon_const_z, var_id_sm_sat_const_z,                      &
           var_id_sm_crit_const_z, var_id_sm_wilt_const_z,                    &
           var_id_hcap_const_z, var_id_hcon_const_z )
! No vertical dimension required, so return the empty string as dimension name
      dim_name_in  = ''
      dim_name_out = ''
      dim_size     = surface_dim_size

!-----------------------------------------------------------------------------
! TOPMODEL variables have no vertical levels
!-----------------------------------------------------------------------------
    CASE ( var_id_fexp, var_id_ti_mean, var_id_ti_sig )
! No vertical dimension required, so return the empty string as dimension name
      dim_name_in  = ''
      dim_name_out = ''
      dim_size     = surface_dim_size

!-----------------------------------------------------------------------------
! Agricultural fraction
!-----------------------------------------------------------------------------
    CASE ( var_id_frac_agr )
! No vertical dimension required, so return the empty string as dimension name
      dim_name_in  = ''
      dim_name_out = ''
      dim_size     = surface_dim_size

!-----------------------------------------------------------------------------
! Urban parameters have no vertical levels
!-----------------------------------------------------------------------------
    CASE ( var_id_wrr, var_id_hwr, var_id_hgt, var_id_ztm, var_id_disp,       &
           var_id_albwl, var_id_albrd, var_id_emisw, var_id_emisr )
! No vertical dimension required, so return the empty string as dimension name
      dim_name_in  = ''
      dim_name_out = ''
      dim_size     = surface_dim_size

!-----------------------------------------------------------------------------
! Variables that can be set as initial conditions
!-----------------------------------------------------------------------------
    CASE ( var_id_gs, var_id_sthzw, var_id_zw, var_id_cv )
! No vertical dimension required, so return the empty string as dimension name
      dim_name_in  = ''
      dim_name_out = ''
      dim_size     = surface_dim_size

    CASE ( var_id_canht, var_id_lai )
      dim_name_in  = pft_dim_name
      dim_name_out = PFT_DIM_NAME_OUT
      dim_size     = pft_dim_size

    CASE ( var_id_cs )
      dim_name_in  = scpool_dim_name
      dim_name_out = SCPOOL_DIM_NAME_OUT
      dim_size     = scpool_dim_size

    CASE ( var_id_sthuf, var_id_t_soil )
      dim_name_in  = soil_dim_name
      dim_name_out = SOIL_DIM_NAME_OUT
      dim_size     = soil_dim_size

    CASE ( var_id_canopy, var_id_nsnow, var_id_rgrain, var_id_rho_snow,       &
           var_id_snow_tile, var_id_snow_depth, var_id_snow_grnd,             &
           var_id_tstar_tile )
      dim_name_in  = tile_dim_name
      dim_name_out = TILE_DIM_NAME_OUT
      dim_size     = tile_dim_size

! CABLE
    CASE ( var_id_isnow_flg3l, var_id_snow_rho1l, var_id_snage_tile )
      dim_name_in  = tile_dim_name
      dim_name_out = TILE_DIM_NAME_OUT
      dim_size     = tile_dim_size

!jhan:snow vars can only be set to a constant in JULES io
    CASE ( var_id_snow_rho1, var_id_snow_rho2, var_id_snow_rho3, & 
           var_id_snow_mass1, var_id_snow_mass2, var_id_snow_mass3, & 
           var_id_snow_depth1, var_id_snow_depth2, var_id_snow_depth3, & 
           var_id_snow_tmp1, var_id_snow_tmp2, var_id_snow_tmp3 )
      dim_name_in  = tile_dim_name
      dim_name_out = TILE_DIM_NAME_OUT
      dim_size     = tile_dim_size

    CASE ( var_id_sthu_tile1, &
           var_id_sthu_tile2, & 
           var_id_sthu_tile3, & 
           var_id_sthu_tile4, & 
           var_id_sthu_tile5, & 
           var_id_sthu_tile6 )
      dim_name_in  = tile_dim_name
      dim_name_out = TILE_DIM_NAME_OUT
      dim_size     = tile_dim_size

    CASE ( var_id_sthf_tile1, &
           var_id_sthf_tile2, & 
           var_id_sthf_tile3, & 
           var_id_sthf_tile4, & 
           var_id_sthf_tile5, & 
           var_id_sthf_tile6 )
      dim_name_in  = tile_dim_name
      dim_name_out = TILE_DIM_NAME_OUT
      dim_size     = tile_dim_size


    CASE ( var_id_smcl_tile1, &
           var_id_smcl_tile2, & 
           var_id_smcl_tile3, & 
           var_id_smcl_tile4, & 
           var_id_smcl_tile5, & 
           var_id_smcl_tile6 )
      dim_name_in  = tile_dim_name
      dim_name_out = TILE_DIM_NAME_OUT
      dim_size     = tile_dim_size


    CASE ( var_id_tsoil_tile1, &
           var_id_tsoil_tile2, & 
           var_id_tsoil_tile3, & 
           var_id_tsoil_tile4, & 
           var_id_tsoil_tile5, & 
           var_id_tsoil_tile6 )
      dim_name_in  = tile_dim_name
      dim_name_out = TILE_DIM_NAME_OUT
      dim_size     = tile_dim_size


!-----------------------------------------------------------------------------
! Driving variables have no vertical levels
!-----------------------------------------------------------------------------
    CASE ( var_id_pstar, var_id_q, var_id_t, var_id_rad_net, var_id_lw_net,   &
           var_id_sw_net, var_id_lw_down, var_id_sw_down, var_id_diff_rad,    &
           var_id_precip, var_id_tot_rain, var_id_tot_snow, var_id_con_rain,  &
           var_id_ls_rain, var_id_con_snow, var_id_ls_snow, var_id_wind,      &
           var_id_u, var_id_v )
! No vertical dimension required, so return the empty string as dimension name
      dim_name_in  = ''
      dim_name_out = ''
      dim_size     = surface_dim_size

!-----------------------------------------------------------------------------
! Ozone has no vertical levels
!-----------------------------------------------------------------------------
    CASE ( var_id_ozone )
! No vertical dimension required, so return the empty string as dimension name
      dim_name_in  = ''
      dim_name_out = ''
      dim_size     = surface_dim_size

!-----------------------------------------------------------------------------
! Output variables that have not been previously defined
!-----------------------------------------------------------------------------
    CASE ( var_id_albedo_land, var_id_canopy_gb, var_id_cs_gb,                &
           var_id_depth_frozen, var_id_depth_unfrozen, var_id_drain,          &
           var_id_elake, var_id_emis_gb, var_id_fch4_wetl, var_id_fsat,       &
           var_id_fsmc_gb, var_id_fwetl, var_id_gpp_gb, var_id_hf_snow_melt,  &
           var_id_land_index, var_id_lice_index, var_id_lit_c_mean,           &
           var_id_lw_up, var_id_npp_gb, var_id_qbase, var_id_qbase_zw,        &
           var_id_resp_p_gb, var_id_resp_s_gb, var_id_resp_s_dr_out,          &
           var_id_runoff, var_id_sat_excess_roff, var_id_smc_avail_top,       &
           var_id_smc_avail_tot, var_id_smc_tot, var_id_snomlt_sub_htf,       &
           var_id_snow_can_gb, var_id_snow_depth_gb, var_id_snow_frac,        &
           var_id_snow_frac_alb, var_id_snow_grnd_gb, var_id_snow_ice_gb,     &
           var_id_snow_liq_gb, var_id_snow_melt_gb, var_id_soil_index,        &
           var_id_sub_surf_roff, var_id_surf_roff, var_id_swet_liq_tot,       &
           var_id_swet_tot, var_id_tfall, var_id_trad )
! No vertical dimension required, so return the empty string as dimension name
      dim_name_in  = ''
      dim_name_out = ''
      dim_size     = surface_dim_size

    CASE ( var_id_c_veg, var_id_flux_o3_stom, var_id_fsmc, var_id_g_leaf,     &
           var_id_g_leaf_day, var_id_g_leaf_dr_out, var_id_g_leaf_phen,       &
           var_id_gpp, var_id_lai_phen, var_id_lit_c, var_id_npp_dr_out,      &
           var_id_npp, var_id_o3_exp_fac, var_id_resp_p,                      &
           var_id_resp_w_dr_out, var_id_resp_w,                               &
           var_id_isoprene, var_id_terpene, var_id_methanol, var_id_acetone)
      dim_name_in  = pft_dim_name
      dim_name_out = PFT_DIM_NAME_OUT
      dim_size     = pft_dim_size

    CASE ( var_id_resp_s )
      dim_name_in  = scpool_dim_name
      dim_name_out = SCPOOL_DIM_NAME_OUT
      dim_size     = scpool_dim_size

    CASE ( var_id_cosz, var_id_diff_frac, var_id_ecan_gb, var_id_ei_gb,       &
           var_id_esoil_gb, var_id_fqw_gb, var_id_ftl_gb,                     &
           var_id_land_albedo_1, var_id_land_albedo_2, var_id_land_albedo_3,  &
           var_id_land_albedo_4, var_id_latent_heat, var_id_q1p5m_gb,         &
           var_id_qw1, var_id_rainfall, var_id_snomlt_surf_htf,               &
           var_id_snowfall, var_id_snow_mass_gb, var_id_surf_ht_flux_gb,      &
           var_id_t1p5m_gb, var_id_taux1, var_id_tauy1, var_id_tl1,           &
           var_id_tstar_gb, var_id_u1, var_id_u10m, var_id_v1, var_id_v10m,   &
           var_id_isoprene_gb, var_id_terpene_gb, var_id_methanol_gb,         &
           var_id_acetone_gb )
! No vertical dimension required, so return the empty string as dimension name
      dim_name_in  = ''
      dim_name_out = ''
      dim_size     = surface_dim_size

    CASE ( var_id_ext, var_id_smcl, var_id_soil_wet, var_id_sthf,             &
           var_id_sthu )
      dim_name_in  = soil_dim_name
      dim_name_out = SOIL_DIM_NAME_OUT
      dim_size     = soil_dim_size

    CASE ( var_id_alb_tile_1, var_id_alb_tile_2, var_id_alb_tile_3,           &
           var_id_alb_tile_4, var_id_anthrop_heat, var_id_catch, var_id_ecan, &
           var_id_ei, var_id_emis, var_id_esoil, var_id_fqw, var_id_ftl,      &
           var_id_gc, var_id_le, var_id_q1p5m, var_id_rad_net_tile,           &
           var_id_snow_can_melt, var_id_snow_can, var_id_snow_grnd_rho,       &
           var_id_snow_ground, var_id_snow_ice_tile, var_id_snow_liq_tile,    &
           var_id_snow_mass, var_id_snow_melt, var_id_surf_ht_flux,           &
           var_id_surf_ht_store, var_id_t1p5m, var_id_tstar, var_id_z0 )
      dim_name_in  = tile_dim_name
      dim_name_out = TILE_DIM_NAME_OUT
      dim_size     = tile_dim_size

    CASE ( var_id_tile_index )
      dim_name_in  = type_dim_name
      dim_name_out = TYPE_DIM_NAME_OUT
      dim_size     = type_dim_size

!-----------------------------------------------------------------------------
! Snow layer variables can be identified using ranges of ids
!-----------------------------------------------------------------------------
    CASE ( var_id_rgrainl:(var_id_rgrainl + ntiles_max - 1) )
      dim_name_in  = snow_dim_name
      dim_name_out = SNOW_DIM_NAME_OUT
      dim_size     = snow_dim_size

    CASE ( var_id_snow_ds:(var_id_snow_ds + ntiles_max - 1) )
      dim_name_in  = snow_dim_name
      dim_name_out = SNOW_DIM_NAME_OUT
      dim_size     = snow_dim_size

    CASE ( var_id_snow_ice:(var_id_snow_ice + ntiles_max - 1) )
      dim_name_in  = snow_dim_name
      dim_name_out = SNOW_DIM_NAME_OUT
      dim_size     = snow_dim_size

    CASE ( var_id_snow_liq:(var_id_snow_liq + ntiles_max - 1) )
      dim_name_in  = snow_dim_name
      dim_name_out = SNOW_DIM_NAME_OUT
      dim_size     = snow_dim_size

    CASE ( var_id_tsnow:(var_id_tsnow + ntiles_max - 1) )
      dim_name_in  = snow_dim_name
      dim_name_out = SNOW_DIM_NAME_OUT
      dim_size     = snow_dim_size

    CASE DEFAULT
      CALL log_fatal("get_var_levs_dim",                                      &
                     "Unrecognised variable - '" //                           &
                     TRIM(get_string_identifier(var_id)) // "'. " //          &
                     "See docs for available variables")
  END SELECT

  RETURN

END SUBROUTINE get_var_levs_dim
#endif
