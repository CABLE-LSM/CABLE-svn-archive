#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/model_interface/get_string_identifier.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-09-27 09:33:34 +0100 (Thu, 27 Sep 2012) $
!!   $LastChangedRevision: 609 $
!!
!!****************************************************************************

FUNCTION get_string_identifier(var_id) RESULT(identifier)

  USE string_utils_mod, ONLY : to_string

  USE max_dimensions, ONLY : ntiles_max

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an integer id for a model variable, returns the string identifier
!   for that variable
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  INTEGER, INTENT(IN) :: var_id  ! The integer id for the variable

! Return type
  CHARACTER(len=IDENTIFIER_LEN) :: identifier  ! The string identifier for
                                               ! the model variable


!-----------------------------------------------------------------------------


  SELECT CASE ( var_id )
    CASE ( var_id_latitude )
      identifier = 'latitude'

    CASE ( var_id_longitude )
      identifier = 'longitude'

    CASE ( var_id_land_fraction )
      identifier = 'land_fraction'

    CASE ( var_id_surf_hgt )
      identifier = 'surf_hgt'

    CASE ( var_id_frac )
      identifier = 'frac'

    CASE ( var_id_albsoil )
      identifier = 'albsoil'

    CASE ( var_id_b )
      identifier = 'b'

    CASE ( var_id_b_const_z )
      identifier = 'b_const_z'

    CASE ( var_id_sathh )
      identifier = 'sathh'

    CASE ( var_id_sathh_const_z )
      identifier = 'sathh_const_z'

    CASE ( var_id_satcon )
      identifier = 'satcon'

    CASE ( var_id_satcon_const_z )
      identifier = 'satcon_const_z'

    CASE ( var_id_sm_sat )
      identifier = 'sm_sat'

    CASE ( var_id_sm_sat_const_z )
      identifier = 'sm_sat_const_z'

    CASE ( var_id_sm_crit )
      identifier = 'sm_crit'

    CASE ( var_id_sm_crit_const_z )
      identifier = 'sm_crit_const_z'

    CASE ( var_id_sm_wilt )
      identifier = 'sm_wilt'

    CASE ( var_id_sm_wilt_const_z )
      identifier = 'sm_wilt_const_z'

    CASE ( var_id_hcap )
      identifier = 'hcap'

    CASE ( var_id_hcap_const_z )
      identifier = 'hcap_const_z'

    CASE ( var_id_hcon )
      identifier = 'hcon'

    CASE ( var_id_hcon_const_z )
      identifier = 'hcon_const_z'

    CASE ( var_id_fexp )
      identifier = 'fexp'

    CASE ( var_id_ti_mean )
      identifier = 'ti_mean'

    CASE ( var_id_ti_sig )
      identifier = 'ti_sig'

    CASE ( var_id_frac_agr )
      identifier = 'frac_agr'

    CASE ( var_id_wrr )
      identifier = 'wrr'

    CASE ( var_id_hwr )
      identifier = 'hwr'

    CASE ( var_id_hgt )
      identifier = 'hgt'

    CASE ( var_id_ztm )
      identifier = 'ztm'

    CASE ( var_id_disp )
      identifier = 'disp'

    CASE ( var_id_albwl )
      identifier = 'albwl'

    CASE ( var_id_albrd )
      identifier = 'albrd'

    CASE ( var_id_emisw )
      identifier = 'emisw'

    CASE ( var_id_emisr )
      identifier = 'emisr'

    CASE ( var_id_canopy )
      identifier = 'canopy'

    CASE ( var_id_cs )
      identifier = 'cs'

    CASE ( var_id_gs )
      identifier = 'gs'

    CASE ( var_id_snow_tile )
      identifier = 'snow_tile'

    CASE ( var_id_sthuf )
      identifier = 'sthuf'

    CASE ( var_id_t_soil )
      identifier = 't_soil'

    CASE ( var_id_tstar_tile )
      identifier = 'tstar_tile'

    CASE ( var_id_lai )
      identifier = 'lai'

    CASE ( var_id_canht )
      identifier = 'canht'

    CASE ( var_id_sthzw )
      identifier = 'sthzw'

    CASE ( var_id_zw )
      identifier = 'zw'

    CASE ( var_id_rgrain )
      identifier = 'rgrain'

    CASE ( var_id_cv )
      identifier = 'cv'

    CASE ( var_id_rho_snow )
      identifier = 'rho_snow'

    CASE ( var_id_snow_depth )
      identifier = 'snow_depth'

    CASE ( var_id_snow_grnd )
      identifier = 'snow_grnd'

    CASE ( var_id_nsnow )
      identifier = 'nsnow'

    CASE ( var_id_pstar )
      identifier = 'pstar'

    CASE ( var_id_q )
      identifier = 'q'

    CASE ( var_id_t )
      identifier = 't'

    CASE ( var_id_rad_net )
      identifier = 'rad_net'

    CASE ( var_id_lw_net )
      identifier = 'lw_net'

    CASE ( var_id_sw_net )
      identifier = 'sw_net'

    CASE ( var_id_lw_down )
      identifier = 'lw_down'

    CASE ( var_id_sw_down )
      identifier = 'sw_down'

    CASE ( var_id_diff_rad )
      identifier = 'diff_rad'

    CASE ( var_id_precip )
      identifier = 'precip'

    CASE ( var_id_tot_rain )
      identifier = 'tot_rain'

    CASE ( var_id_tot_snow )
      identifier = 'tot_snow'

    CASE ( var_id_con_rain )
      identifier = 'con_rain'

    CASE ( var_id_ls_rain )
      identifier = 'ls_rain'

    CASE ( var_id_con_snow )
      identifier = 'con_snow'

    CASE ( var_id_ls_snow )
      identifier = 'ls_snow'

    CASE ( var_id_wind )
      identifier = 'wind'

    CASE ( var_id_u )
      identifier = 'u'

    CASE ( var_id_v )
      identifier = 'v'

    CASE ( var_id_albedo_land )
      identifier = 'albedo_land'

    CASE ( var_id_canopy_gb )
      identifier = 'canopy_gb'

    CASE ( var_id_cs_gb )
      identifier = 'cs_gb'

    CASE ( var_id_depth_frozen )
      identifier = 'depth_frozen'

    CASE ( var_id_depth_unfrozen )
      identifier = 'depth_unfrozen'

    CASE ( var_id_drain )
      identifier = 'drain'

    CASE ( var_id_elake )
      identifier = 'elake'

    CASE ( var_id_emis_gb )
      identifier = 'emis_gb'

    CASE ( var_id_fch4_wetl )
      identifier = 'fch4_wetl'

    CASE ( var_id_fsat )
      identifier = 'fsat'

    CASE ( var_id_fsmc_gb )
      identifier = 'fsmc_gb'

    CASE ( var_id_fwetl )
      identifier = 'fwetl'

    CASE ( var_id_gpp_gb )
      identifier = 'gpp_gb'

    CASE ( var_id_hf_snow_melt )
      identifier = 'hf_snow_melt'

    CASE ( var_id_land_index )
      identifier = 'land_index'

    CASE ( var_id_lice_index )
      identifier = 'lice_index'

    CASE ( var_id_lit_c_mean )
      identifier = 'lit_c_mean'

    CASE ( var_id_lw_up )
      identifier = 'lw_up'

    CASE ( var_id_npp_gb )
      identifier = 'npp_gb'

    CASE ( var_id_qbase )
      identifier = 'qbase'

    CASE ( var_id_qbase_zw )
      identifier = 'qbase_zw'

    CASE ( var_id_resp_p_gb )
      identifier = 'resp_p_gb'

    CASE ( var_id_resp_s_gb )
      identifier = 'resp_s_gb'

    CASE ( var_id_resp_s_dr_out )
      identifier = 'resp_s_dr_out'

    CASE ( var_id_runoff )
      identifier = 'runoff'

    CASE ( var_id_sat_excess_roff )
      identifier = 'sat_excess_roff'

    CASE ( var_id_smc_avail_top )
      identifier = 'smc_avail_top'

    CASE ( var_id_smc_avail_tot )
      identifier = 'smc_avail_tot'

    CASE ( var_id_smc_tot )
      identifier = 'smc_tot'

    CASE ( var_id_snomlt_sub_htf )
      identifier = 'snomlt_sub_htf'

    CASE ( var_id_snow_can_gb )
      identifier = 'snow_can_gb'

    CASE ( var_id_snow_depth_gb )
      identifier = 'snow_depth_gb'

    CASE ( var_id_snow_frac )
      identifier = 'snow_frac'

    CASE ( var_id_snow_frac_alb )
      identifier = 'snow_frac_alb'

    CASE ( var_id_snow_grnd_gb )
      identifier = 'snow_grnd_gb'

    CASE ( var_id_snow_ice_gb )
      identifier = 'snow_ice_gb'

    CASE ( var_id_snow_liq_gb )
      identifier = 'snow_liq_gb'

    CASE ( var_id_snow_melt_gb )
      identifier = 'snow_melt_gb'

    CASE ( var_id_soil_index )
      identifier = 'soil_index'

    CASE ( var_id_sub_surf_roff )
      identifier = 'sub_surf_roff'

    CASE ( var_id_surf_roff )
      identifier = 'surf_roff'

    CASE ( var_id_swet_liq_tot )
      identifier = 'swet_liq_tot'

    CASE ( var_id_swet_tot )
      identifier = 'swet_tot'

    CASE ( var_id_tfall )
      identifier = 'tfall'

    CASE ( var_id_trad )
      identifier = 'trad'

    CASE ( var_id_c_veg )
      identifier = 'c_veg'

    CASE ( var_id_flux_o3_stom )
      identifier = 'flux_o3_stom'

    CASE ( var_id_fsmc )
      identifier = 'fsmc'

    CASE ( var_id_g_leaf )
      identifier = 'g_leaf'

    CASE ( var_id_g_leaf_day )
      identifier = 'g_leaf_day'

    CASE ( var_id_g_leaf_dr_out )
      identifier = 'g_leaf_dr_out'

    CASE ( var_id_g_leaf_phen )
      identifier = 'g_leaf_phen'

    CASE ( var_id_gpp )
      identifier = 'gpp'

    CASE ( var_id_lai_phen )
      identifier = 'lai_phen'

    CASE ( var_id_lit_c )
      identifier = 'lit_c'

    CASE ( var_id_npp_dr_out )
      identifier = 'npp_dr_out'

    CASE ( var_id_npp )
      identifier = 'npp'

    CASE ( var_id_o3_exp_fac )
      identifier = 'o3_exp_fac'

    CASE ( var_id_resp_p )
      identifier = 'resp_p'

    CASE ( var_id_resp_w_dr_out )
      identifier = 'resp_w_dr_out'

    CASE ( var_id_resp_w )
      identifier = 'resp_w'

    CASE ( var_id_resp_s )
      identifier = 'resp_s'

    CASE ( var_id_cosz )
      identifier = 'cosz'

    CASE ( var_id_diff_frac )
      identifier = 'diff_frac'

    CASE ( var_id_ecan_gb )
      identifier = 'ecan_gb'

    CASE ( var_id_ei_gb )
      identifier = 'ei_gb'

    CASE ( var_id_esoil_gb )
      identifier = 'esoil_gb'

    CASE ( var_id_fqw_gb )
      identifier = 'fqw_gb'

    CASE ( var_id_ftl_gb )
      identifier = 'ftl_gb'

    CASE ( var_id_land_albedo_1 )
      identifier = 'land_albedo_1'

    CASE ( var_id_land_albedo_2 )
      identifier = 'land_albedo_2'

    CASE ( var_id_land_albedo_3 )
      identifier = 'land_albedo_3'

    CASE ( var_id_land_albedo_4 )
      identifier = 'land_albedo_4'

    CASE ( var_id_latent_heat )
      identifier = 'latent_heat'

    CASE ( var_id_q1p5m_gb )
      identifier = 'q1p5m_gb'

    CASE ( var_id_qw1 )
      identifier = 'qw1'

    CASE ( var_id_rainfall )
      identifier = 'rainfall'

    CASE ( var_id_snomlt_surf_htf )
      identifier = 'snomlt_surf_htf'

    CASE ( var_id_snowfall )
      identifier = 'snowfall'

    CASE ( var_id_snow_mass_gb )
      identifier = 'snow_mass_gb'

    CASE ( var_id_surf_ht_flux_gb )
      identifier = 'surf_ht_flux_gb'

    CASE ( var_id_t1p5m_gb )
      identifier = 't1p5m_gb'

    CASE ( var_id_taux1 )
      identifier = 'taux1'

    CASE ( var_id_tauy1 )
      identifier = 'tauy1'

    CASE ( var_id_tl1 )
      identifier = 'tl1'

    CASE ( var_id_tstar_gb )
      identifier = 'tstar_gb'

    CASE ( var_id_u1 )
      identifier = 'u1'

    CASE ( var_id_u10m )
      identifier = 'u10m'

    CASE ( var_id_v1 )
      identifier = 'v1'

    CASE ( var_id_v10m )
      identifier = 'v10m'

    CASE ( var_id_ext )
      identifier = 'ext'

    CASE ( var_id_smcl )
      identifier = 'smcl'

    CASE ( var_id_soil_wet )
      identifier = 'soil_wet'

    CASE ( var_id_sthf )
      identifier = 'sthf'

    CASE ( var_id_sthu )
      identifier = 'sthu'

    CASE ( var_id_alb_tile_1 )
      identifier = 'alb_tile_1'

    CASE ( var_id_alb_tile_2 )
      identifier = 'alb_tile_2'

    CASE ( var_id_alb_tile_3 )
      identifier = 'alb_tile_3'

    CASE ( var_id_alb_tile_4 )
      identifier = 'alb_tile_4'

    CASE ( var_id_anthrop_heat )
      identifier = 'anthrop_heat'

    CASE ( var_id_catch )
      identifier = 'catch'

    CASE ( var_id_ecan )
      identifier = 'ecan'

    CASE ( var_id_ei )
      identifier = 'ei'

    CASE ( var_id_emis )
      identifier = 'emis'

    CASE ( var_id_esoil )
      identifier = 'esoil'

    CASE ( var_id_fqw )
      identifier = 'fqw'

    CASE ( var_id_ftl )
      identifier = 'ftl'

    CASE ( var_id_gc )
      identifier = 'gc'

    CASE ( var_id_le )
      identifier = 'le'

    CASE ( var_id_q1p5m )
      identifier = 'q1p5m'

    CASE ( var_id_rad_net_tile )
      identifier = 'rad_net_tile'

    CASE ( var_id_snow_can_melt )
      identifier = 'snow_can_melt'

    CASE ( var_id_snow_can )
      identifier = 'snow_can'

    CASE ( var_id_snow_grnd_rho )
      identifier = 'snow_grnd_rho'

    CASE ( var_id_snow_ground )
      identifier = 'snow_ground'

    CASE ( var_id_snow_ice_tile )
      identifier = 'snow_ice_tile'

    CASE ( var_id_snow_liq_tile )
      identifier = 'snow_liq_tile'

    CASE ( var_id_snow_mass )
      identifier = 'snow_mass'

    CASE ( var_id_snow_melt )
      identifier = 'snow_melt'

    CASE ( var_id_surf_ht_flux )
      identifier = 'surf_ht_flux'

    CASE ( var_id_surf_ht_store )
      identifier = 'surf_ht_store'

    CASE ( var_id_t1p5m )
      identifier = 't1p5m'

    CASE ( var_id_tstar )
      identifier = 'tstar'

    CASE ( var_id_z0 )
      identifier = 'z0'

    CASE ( var_id_tile_index )
      identifier = 'tile_index'

    CASE ( var_id_ozone )
      identifier = 'ozone'

    CASE ( var_id_isoprene_gb )
      identifier = 'isoprene_gb'

    CASE ( var_id_isoprene )
      identifier = 'isoprene'

    CASE ( var_id_terpene_gb )
      identifier = 'terpene_gb'

    CASE ( var_id_terpene )
      identifier = 'terpene'

    CASE ( var_id_methanol_gb )
      identifier = 'methanol_gb'

    CASE ( var_id_methanol )
      identifier = 'methanol'

    CASE ( var_id_acetone_gb )
      identifier = 'acetone_gb'

    CASE ( var_id_acetone )
      identifier = 'acetone'
! CABLE
    CASE ( var_id_isnow_flg3l)
      identifier = 'isnow_flg3l'
    
    CASE ( var_id_snow_rho1l)
      identifier = 'snow_rho1l'

    CASE ( var_id_snage_tile)
      identifier = 'snage_tile'

    CASE ( var_id_snow_rho1 )
      identifier = 'snow_rho1'

    CASE ( var_id_snow_rho2 )
      identifier = 'snow_rho2'

    CASE ( var_id_snow_rho3 )
      identifier = 'snow_rho3'
    
    CASE ( var_id_snow_mass1 )
      identifier = 'snow_mass1'

    CASE ( var_id_snow_mass2 )
      identifier = 'snow_mass2'

    CASE ( var_id_snow_mass3 )
      identifier = 'snow_mass3'

    CASE ( var_id_snow_depth1 )
      identifier = 'snow_depth1'

    CASE ( var_id_snow_depth2 )
      identifier = 'snow_depth2'

    CASE ( var_id_snow_depth3 )
      identifier = 'snow_depth3'

    CASE ( var_id_snow_tmp1 )
      identifier = 'snow_tmp1'

    CASE ( var_id_snow_tmp2 )
      identifier = 'snow_tmp2'

    CASE ( var_id_snow_tmp3 )
      identifier = 'snow_tmp3'

    CASE ( var_id_sthu_tile1)
      identifier = 'sthu_tile1'

    CASE ( var_id_sthu_tile2)
      identifier = 'sthu_tile2'

    CASE ( var_id_sthu_tile3)
      identifier = 'sthu_tile3'

    CASE ( var_id_sthu_tile4)
      identifier = 'sthu_tile4'

    CASE ( var_id_sthu_tile5)
      identifier = 'sthu_tile5'

    CASE ( var_id_sthu_tile6)
      identifier = 'sthu_tile6'

    CASE ( var_id_sthf_tile1)
      identifier = 'sthf_tile1'

    CASE ( var_id_sthf_tile2)
      identifier = 'sthf_tile2'

    CASE ( var_id_sthf_tile3)
      identifier = 'sthf_tile3'

    CASE ( var_id_sthf_tile4)
      identifier = 'sthf_tile4'

    CASE ( var_id_sthf_tile5)
      identifier = 'sthf_tile5'

    CASE ( var_id_sthf_tile6)
      identifier = 'sthf_tile6'

    CASE ( var_id_smcl_tile1)
      identifier = 'smcl_tile1'

    CASE ( var_id_smcl_tile2)
      identifier = 'smcl_tile2'

    CASE ( var_id_smcl_tile3)
      identifier = 'smcl_tile3'

    CASE ( var_id_smcl_tile4)
      identifier = 'smcl_tile4'

    CASE ( var_id_smcl_tile5)
      identifier = 'smcl_tile5'

    CASE ( var_id_smcl_tile6)
      identifier = 'smcl_tile6'

    CASE ( var_id_tsoil_tile1)
      identifier = 'tsoil_tile1'

    CASE ( var_id_tsoil_tile2)
      identifier = 'tsoil_tile2'

    CASE ( var_id_tsoil_tile3)
      identifier = 'tsoil_tile3'

    CASE ( var_id_tsoil_tile4)
      identifier = 'tsoil_tile4'

    CASE ( var_id_tsoil_tile5)
      identifier = 'tsoil_tile5'

    CASE ( var_id_tsoil_tile6)
      identifier = 'tsoil_tile6'



!-----------------------------------------------------------------------------
! Snow layer variables can be identified using ranges of ids
!-----------------------------------------------------------------------------
    CASE ( var_id_rgrainl:(var_id_rgrainl + ntiles_max - 1) )
      identifier = 'rgrainl_' // TRIM(to_string(var_id - var_id_rgrainl + 1))

    CASE ( var_id_snow_ds:(var_id_snow_ds + ntiles_max - 1) )
      identifier = 'snow_ds_' // TRIM(to_string(var_id - var_id_snow_ds + 1))

    CASE ( var_id_snow_ice:(var_id_snow_ice + ntiles_max - 1) )
      identifier = 'snow_ice_' // TRIM(to_string(var_id - var_id_snow_ice + 1))

    CASE ( var_id_snow_liq:(var_id_snow_liq + ntiles_max - 1) )
      identifier = 'snow_liq_' // TRIM(to_string(var_id - var_id_snow_liq + 1))

    CASE ( var_id_tsnow:(var_id_tsnow + ntiles_max - 1) )
      identifier = 'tsnow_' // TRIM(to_string(var_id - var_id_tsnow + 1))
! CABLE	
!    CASE ( var_id_isnow_flg3l:(var_id_isnow_flg3l + ntiles_max - 1) )
!      identifier = 'isnow_flg3l_' // TRIM(to_string(var_id - var_id_isnow_flg3l + 1))

    CASE DEFAULT
      CALL log_fatal("get_string_identifier",                                 &
                     "Unrecognised variable id - '" // TRIM(to_string(var_id)))
  END SELECT

  RETURN

END FUNCTION get_string_identifier
#endif
