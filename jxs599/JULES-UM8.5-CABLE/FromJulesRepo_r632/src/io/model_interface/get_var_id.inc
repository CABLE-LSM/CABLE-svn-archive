#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION get_var_id(identifier) RESULT(var_id)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given a string identifier for a model variable, returns the integer id for
!   that variable
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  CHARACTER(len=*), INTENT(IN) :: identifier  ! Identifies the model variable

! Return type
  INTEGER :: var_id  ! The integer id for the variable

! Work variables


!-----------------------------------------------------------------------------


  SELECT CASE ( identifier )
    CASE ( 'latitude' )
      var_id = var_id_latitude

    CASE ( 'longitude' )
      var_id = var_id_longitude

    CASE ( 'land_fraction' )
      var_id = var_id_land_fraction

    CASE ( 'surf_hgt' )
      var_id = var_id_surf_hgt

    CASE ( 'frac' )
      var_id = var_id_frac

    CASE ( 'albsoil' )
      var_id = var_id_albsoil

    CASE ( 'albobs_sw' )
      var_id = var_id_albobs_sw

    CASE ( 'albobs_vis' )
      var_id = var_id_albobs_vis

    CASE ( 'albobs_nir' )
      var_id = var_id_albobs_nir

    CASE ( 'b' )
      var_id = var_id_b

    CASE ( 'b_const_z' )
      var_id = var_id_b_const_z

    CASE ( 'sathh' )
      var_id = var_id_sathh

    CASE ( 'sathh_const_z' )
      var_id = var_id_sathh_const_z

    CASE ( 'satcon' )
      var_id = var_id_satcon

    CASE ( 'satcon_const_z' )
      var_id = var_id_satcon_const_z

    CASE ( 'sm_sat' )
      var_id = var_id_sm_sat

    CASE ( 'sm_sat_const_z' )
      var_id = var_id_sm_sat_const_z

    CASE ( 'sm_crit' )
      var_id = var_id_sm_crit

    CASE ( 'sm_crit_const_z' )
      var_id = var_id_sm_crit_const_z

    CASE ( 'sm_wilt' )
      var_id = var_id_sm_wilt

    CASE ( 'sm_wilt_const_z' )
      var_id = var_id_sm_wilt_const_z

    CASE ( 'hcap' )
      var_id = var_id_hcap

    CASE ( 'hcap_const_z' )
      var_id = var_id_hcap_const_z

    CASE ( 'hcon' )
      var_id = var_id_hcon

    CASE ( 'hcon_const_z' )
      var_id = var_id_hcon_const_z

    CASE ( 'fexp' )
      var_id = var_id_fexp

    CASE ( 'ti_mean' )
      var_id = var_id_ti_mean

    CASE ( 'ti_sig' )
      var_id = var_id_ti_sig

    CASE ( 'frac_agr' )
      var_id = var_id_frac_agr

    CASE ( 'wrr' )
      var_id = var_id_wrr

    CASE ( 'hwr' )
      var_id = var_id_hwr

    CASE ( 'hgt' )
      var_id = var_id_hgt

    CASE ( 'ztm' )
      var_id = var_id_ztm

    CASE ( 'disp' )
      var_id = var_id_disp

    CASE ( 'albwl' )
      var_id = var_id_albwl

    CASE ( 'albrd' )
      var_id = var_id_albrd

    CASE ( 'emisw' )
      var_id = var_id_emisw

    CASE ( 'emisr' )
      var_id = var_id_emisr

    CASE ( 'canopy' )
      var_id = var_id_canopy

    CASE ( 'cs' )
      var_id = var_id_cs

    CASE ( 'gs' )
      var_id = var_id_gs

    CASE ( 'snow_tile' )
      var_id = var_id_snow_tile

    CASE ( 'sthuf' )
      var_id = var_id_sthuf

    CASE ( 't_soil' )
      var_id = var_id_t_soil

    CASE ( 'tstar_tile' )
      var_id = var_id_tstar_tile

    CASE ( 'lai' )
      var_id = var_id_lai

    CASE ( 'canht' )
      var_id = var_id_canht

    CASE ( 'sthzw' )
      var_id = var_id_sthzw

    CASE ( 'zw' )
      var_id = var_id_zw

    CASE ( 'rgrain' )
      var_id = var_id_rgrain

    CASE ( 'cv' )
      var_id = var_id_cv

    CASE ( 'rho_snow' )
      var_id = var_id_rho_snow

    CASE ( 'snow_depth' )
      var_id = var_id_snow_depth

    CASE ( 'snow_grnd' )
      var_id = var_id_snow_grnd

    CASE ( 'nsnow' )
      var_id = var_id_nsnow

    CASE ( 'snow_ds' )
      var_id = var_id_snow_ds

    CASE ( 'snow_ice' )
      var_id = var_id_snow_ice

    CASE ( 'snow_liq' )
      var_id = var_id_snow_liq

    CASE ( 'tsnow' )
      var_id = var_id_tsnow

    CASE ( 'rgrainl' )
      var_id = var_id_rgrainl

    CASE ( 'pstar' )
      var_id = var_id_pstar

    CASE ( 'q' )
      var_id = var_id_q

    CASE ( 't' )
      var_id = var_id_t

    CASE ( 'rad_net' )
      var_id = var_id_rad_net

    CASE ( 'lw_net' )
      var_id = var_id_lw_net

    CASE ( 'sw_net' )
      var_id = var_id_sw_net

    CASE ( 'lw_down' )
      var_id = var_id_lw_down

    CASE ( 'sw_down' )
      var_id = var_id_sw_down

    CASE ( 'diff_rad' )
      var_id = var_id_diff_rad

    CASE ( 'precip' )
      var_id = var_id_precip

    CASE ( 'tot_rain' )
      var_id = var_id_tot_rain

    CASE ( 'tot_snow' )
      var_id = var_id_tot_snow

    CASE ( 'con_rain' )
      var_id = var_id_con_rain

    CASE ( 'ls_rain' )
      var_id = var_id_ls_rain

    CASE ( 'con_snow' )
      var_id = var_id_con_snow

    CASE ( 'ls_snow' )
      var_id = var_id_ls_snow

    CASE ( 'wind' )
      var_id = var_id_wind

    CASE ( 'u' )
      var_id = var_id_u

    CASE ( 'v' )
      var_id = var_id_v

    CASE ( 'albedo_land' )
      var_id = var_id_albedo_land

    CASE ( 'canopy_gb' )
      var_id = var_id_canopy_gb

    CASE ( 'cs_gb' )
      var_id = var_id_cs_gb

    CASE ( 'depth_frozen' )
      var_id = var_id_depth_frozen

    CASE ( 'depth_unfrozen' )
      var_id = var_id_depth_unfrozen

    CASE ( 'drain' )
      var_id = var_id_drain

    CASE ( 'elake' )
      var_id = var_id_elake

    CASE ( 'emis_gb' )
      var_id = var_id_emis_gb

    CASE ( 'fch4_wetl' )
      var_id = var_id_fch4_wetl

    CASE ( 'fsat' )
      var_id = var_id_fsat

    CASE ( 'fsmc_gb' )
      var_id = var_id_fsmc_gb

    CASE ( 'fwetl' )
      var_id = var_id_fwetl

    CASE ( 'gpp_gb' )
      var_id = var_id_gpp_gb

    CASE ( 'hf_snow_melt' )
      var_id = var_id_hf_snow_melt

    CASE ( 'land_index' )
      var_id = var_id_land_index

    CASE ( 'lice_index' )
      var_id = var_id_lice_index

    CASE ( 'lit_c_mean' )
      var_id = var_id_lit_c_mean

    CASE ( 'lw_up' )
      var_id = var_id_lw_up

    CASE ( 'npp_gb' )
      var_id = var_id_npp_gb

    CASE ( 'qbase' )
      var_id = var_id_qbase

    CASE ( 'qbase_zw' )
      var_id = var_id_qbase_zw

    CASE ( 'resp_p_gb' )
      var_id = var_id_resp_p_gb

    CASE ( 'resp_s_gb' )
      var_id = var_id_resp_s_gb

    CASE ( 'resp_s_dr_out' )
      var_id = var_id_resp_s_dr_out

    CASE ( 'runoff' )
      var_id = var_id_runoff

    CASE ( 'sat_excess_roff' )
      var_id = var_id_sat_excess_roff

    CASE ( 'smc_avail_top' )
      var_id = var_id_smc_avail_top

    CASE ( 'smc_avail_tot' )
      var_id = var_id_smc_avail_tot

    CASE ( 'smc_tot' )
      var_id = var_id_smc_tot

    CASE ( 'snomlt_sub_htf' )
      var_id = var_id_snomlt_sub_htf

    CASE ( 'snow_can_gb' )
      var_id = var_id_snow_can_gb

    CASE ( 'snow_depth_gb' )
      var_id = var_id_snow_depth_gb

    CASE ( 'snow_frac' )
      var_id = var_id_snow_frac

    CASE ( 'snow_grnd_gb' )
      var_id = var_id_snow_grnd_gb

    CASE ( 'snow_ice_gb' )
      var_id = var_id_snow_ice_gb

    CASE ( 'snow_liq_gb' )
      var_id = var_id_snow_liq_gb

    CASE ( 'snow_melt_gb' )
      var_id = var_id_snow_melt_gb

    CASE ( 'soil_index' )
      var_id = var_id_soil_index

    CASE ( 'sub_surf_roff' )
      var_id = var_id_sub_surf_roff

    CASE ( 'surf_roff' )
      var_id = var_id_surf_roff

    CASE ( 'swet_liq_tot' )
      var_id = var_id_swet_liq_tot

    CASE ( 'swet_tot' )
      var_id = var_id_swet_tot

    CASE ( 'tfall' )
      var_id = var_id_tfall

    CASE ( 'trad' )
      var_id = var_id_trad

    CASE ( 'c_veg' )
      var_id = var_id_c_veg

    CASE ( 'flux_o3_stom' )
      var_id = var_id_flux_o3_stom

    CASE ( 'fsmc' )
      var_id = var_id_fsmc

    CASE ( 'g_leaf' )
      var_id = var_id_g_leaf

    CASE ( 'g_leaf_day' )
      var_id = var_id_g_leaf_day

    CASE ( 'g_leaf_dr_out' )
      var_id = var_id_g_leaf_dr_out

    CASE ( 'g_leaf_phen' )
      var_id = var_id_g_leaf_phen

    CASE ( 'gpp' )
      var_id = var_id_gpp

    CASE ( 'lai_phen' )
      var_id = var_id_lai_phen

    CASE ( 'lit_c' )
      var_id = var_id_lit_c

    CASE ( 'npp_dr_out' )
      var_id = var_id_npp_dr_out

    CASE ( 'npp' )
      var_id = var_id_npp

    CASE ( 'o3_exp_fac' )
      var_id = var_id_o3_exp_fac

    CASE ( 'resp_p' )
      var_id = var_id_resp_p

    CASE ( 'resp_w_dr_out' )
      var_id = var_id_resp_w_dr_out

    CASE ( 'resp_w' )
      var_id = var_id_resp_w

    CASE ( 'resp_s' )
      var_id = var_id_resp_s

    CASE ( 'cosz' )
      var_id = var_id_cosz

    CASE ( 'diff_frac' )
      var_id = var_id_diff_frac

    CASE ( 'ecan_gb' )
      var_id = var_id_ecan_gb

    CASE ( 'ei_gb' )
      var_id = var_id_ei_gb

    CASE ( 'esoil_gb' )
      var_id = var_id_esoil_gb

    CASE ( 'fqw_gb' )
      var_id = var_id_fqw_gb

    CASE ( 'ftl_gb' )
      var_id = var_id_ftl_gb

    CASE ( 'land_albedo_1' )
      var_id = var_id_land_albedo_1

    CASE ( 'land_albedo_2' )
      var_id = var_id_land_albedo_2

    CASE ( 'land_albedo_3' )
      var_id = var_id_land_albedo_3

    CASE ( 'land_albedo_4' )
      var_id = var_id_land_albedo_4

    CASE ( 'latent_heat' )
      var_id = var_id_latent_heat

    CASE ( 'q1p5m_gb' )
      var_id = var_id_q1p5m_gb

    CASE ( 'qw1' )
      var_id = var_id_qw1

    CASE ( 'rainfall' )
      var_id = var_id_rainfall

    CASE ( 'snomlt_surf_htf' )
      var_id = var_id_snomlt_surf_htf

    CASE ( 'snowfall' )
      var_id = var_id_snowfall

    CASE ( 'snow_mass_gb' )
      var_id = var_id_snow_mass_gb

    CASE ( 'surf_ht_flux_gb' )
      var_id = var_id_surf_ht_flux_gb

    CASE ( 't1p5m_gb' )
      var_id = var_id_t1p5m_gb

    CASE ( 'taux1' )
      var_id = var_id_taux1

    CASE ( 'tauy1' )
      var_id = var_id_tauy1

    CASE ( 'tl1' )
      var_id = var_id_tl1

    CASE ( 'tstar_gb' )
      var_id = var_id_tstar_gb

    CASE ( 'u1' )
      var_id = var_id_u1

    CASE ( 'u10m' )
      var_id = var_id_u10m

    CASE ( 'v1' )
      var_id = var_id_v1

    CASE ( 'v10m' )
      var_id = var_id_v10m

    CASE ( 'ext' )
      var_id = var_id_ext

    CASE ( 'smcl' )
      var_id = var_id_smcl

    CASE ( 'soil_wet' )
      var_id = var_id_soil_wet

    CASE ( 'sthf' )
      var_id = var_id_sthf

    CASE ( 'sthu' )
      var_id = var_id_sthu

    CASE ( 'alb_tile_1' )
      var_id = var_id_alb_tile_1

    CASE ( 'alb_tile_2' )
      var_id = var_id_alb_tile_2

    CASE ( 'alb_tile_3' )
      var_id = var_id_alb_tile_3

    CASE ( 'alb_tile_4' )
      var_id = var_id_alb_tile_4

    CASE ( 'anthrop_heat' )
      var_id = var_id_anthrop_heat

    CASE ( 'catch' )
      var_id = var_id_catch

    CASE ( 'ecan' )
      var_id = var_id_ecan

    CASE ( 'ei' )
      var_id = var_id_ei

    CASE ( 'emis' )
      var_id = var_id_emis

    CASE ( 'esoil' )
      var_id = var_id_esoil

    CASE ( 'fqw' )
      var_id = var_id_fqw

    CASE ( 'ftl' )
      var_id = var_id_ftl

    CASE ( 'gc' )
      var_id = var_id_gc

    CASE ( 'le' )
      var_id = var_id_le

    CASE ( 'q1p5m' )
      var_id = var_id_q1p5m

    CASE ( 'rad_net_tile' )
      var_id = var_id_rad_net_tile

    CASE ( 'snow_can_melt' )
      var_id = var_id_snow_can_melt

    CASE ( 'snow_can' )
      var_id = var_id_snow_can

    CASE ( 'snow_grnd_rho' )
      var_id = var_id_snow_grnd_rho

    CASE ( 'snow_ground' )
      var_id = var_id_snow_ground

    CASE ( 'snow_ice_tile' )
      var_id = var_id_snow_ice_tile

    CASE ( 'snow_liq_tile' )
      var_id = var_id_snow_liq_tile

    CASE ( 'snow_mass' )
      var_id = var_id_snow_mass

    CASE ( 'snow_melt' )
      var_id = var_id_snow_melt

    CASE ( 'surf_ht_flux' )
      var_id = var_id_surf_ht_flux

    CASE ( 'surf_ht_store' )
      var_id = var_id_surf_ht_store

    CASE ( 't1p5m' )
      var_id = var_id_t1p5m

    CASE ( 'tstar' )
      var_id = var_id_tstar

    CASE ( 'z0' )
      var_id = var_id_z0

    CASE ( 'tile_index' )
      var_id = var_id_tile_index

    CASE ( 'ozone' )
      var_id = var_id_ozone

    CASE ( 'isoprene_gb' )
      var_id = var_id_isoprene_gb

    CASE ( 'isoprene' )
      var_id = var_id_isoprene

    CASE ( 'terpene_gb' )
      var_id = var_id_terpene_gb

    CASE ( 'terpene' )
      var_id = var_id_terpene

    CASE ( 'methanol_gb' )
      var_id = var_id_methanol_gb

    CASE ( 'methanol' )
      var_id = var_id_methanol

    CASE ( 'acetone_gb' )
      var_id = var_id_acetone_gb

    CASE ( 'acetone' )
      var_id = var_id_acetone

    CASE DEFAULT
      CALL log_fatal("get_var_id",                                            &
                     "Unrecognised variable identifier - '" //                &
                     TRIM(identifier) // "'. See docs for available variables")
  END SELECT

  RETURN

END FUNCTION get_var_id
#endif
