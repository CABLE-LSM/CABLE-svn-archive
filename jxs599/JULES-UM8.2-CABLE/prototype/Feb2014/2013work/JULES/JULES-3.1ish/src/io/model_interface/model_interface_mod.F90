#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/model_interface/model_interface_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-09-27 09:33:34 +0100 (Thu, 27 Sep 2012) $
!!   $LastChangedRevision: 609 $
!!
!!****************************************************************************

MODULE model_interface_mod

  USE io_constants, ONLY : MAX_SDF_NAME_LEN

  USE max_dimensions, ONLY : ntiles_max

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

!CABLE
  USE cable_vars_mod!, ONLY : isnow_flg3l, snow_rho1l, snage_tile

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This module provides the interface between the IO routines and the
!   model variables
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
! Variable ids
  INTEGER, PARAMETER ::                                                       &
    var_id_latitude        = 1,                                               &
    var_id_longitude       = var_id_latitude + 1,                             &
    var_id_land_fraction   = var_id_longitude + 1,                            &
    var_id_surf_hgt        = var_id_land_fraction + 1,                        &
    var_id_frac            = var_id_surf_hgt + 1,                             &
    var_id_albsoil         = var_id_frac + 1,                                 &
    var_id_b               = var_id_albsoil + 1,                              &
    var_id_b_const_z       = var_id_b + 1,                                    &
    var_id_sathh           = var_id_b_const_z + 1,                            &
    var_id_sathh_const_z   = var_id_sathh + 1,                                &
    var_id_satcon          = var_id_sathh_const_z + 1,                        &
    var_id_satcon_const_z  = var_id_satcon + 1,                               &
    var_id_sm_sat          = var_id_satcon_const_z + 1,                       &
    var_id_sm_sat_const_z  = var_id_sm_sat + 1,                               &
    var_id_sm_crit         = var_id_sm_sat_const_z + 1,                       &
    var_id_sm_crit_const_z = var_id_sm_crit + 1,                              &
    var_id_sm_wilt         = var_id_sm_crit_const_z + 1,                      &
    var_id_sm_wilt_const_z = var_id_sm_wilt + 1,                              &
    var_id_hcap            = var_id_sm_wilt_const_z + 1,                      &
    var_id_hcap_const_z    = var_id_hcap + 1,                                 &
    var_id_hcon            = var_id_hcap_const_z + 1,                         &
    var_id_hcon_const_z    = var_id_hcon + 1,                                 &
    var_id_fexp            = var_id_hcon_const_z + 1,                         &
    var_id_ti_mean         = var_id_fexp + 1,                                 &
    var_id_ti_sig          = var_id_ti_mean + 1,                              &
    var_id_frac_agr        = var_id_ti_sig + 1,                               &
    var_id_wrr             = var_id_frac_agr + 1,                             &
    var_id_hwr             = var_id_wrr + 1,                                  &
    var_id_hgt             = var_id_hwr + 1,                                  &
    var_id_ztm             = var_id_hgt + 1,                                  &
    var_id_disp            = var_id_ztm + 1,                                  &
    var_id_albwl           = var_id_disp + 1,                                 &
    var_id_albrd           = var_id_albwl + 1,                                &
    var_id_emisw           = var_id_albrd + 1,                                &
    var_id_emisr           = var_id_emisw + 1,                                &
    var_id_canopy          = var_id_emisr + 1,                                &
    var_id_cs              = var_id_canopy + 1,                               &
    var_id_gs              = var_id_cs + 1,                                   &
    var_id_snow_tile       = var_id_gs + 1,                                   &
    var_id_sthuf           = var_id_snow_tile + 1,                            &
    var_id_t_soil          = var_id_sthuf + 1,                                &
    var_id_tstar_tile      = var_id_t_soil + 1,                               &
    var_id_lai             = var_id_tstar_tile + 1,                           &
    var_id_canht           = var_id_lai + 1,                                  &
    var_id_sthzw           = var_id_canht + 1,                                &
    var_id_zw              = var_id_sthzw + 1,                                &
    var_id_rgrain          = var_id_zw + 1,                                   &
    var_id_cv              = var_id_rgrain + 1,                               &
    var_id_rho_snow        = var_id_cv + 1,                                   &
    var_id_snow_depth      = var_id_rho_snow + 1,                             &
    var_id_snow_grnd       = var_id_snow_depth + 1,                           &
    var_id_nsnow           = var_id_snow_grnd + 1,                            &
! Snow variables need a range of ids to be allocated to them to accomodate
! the tiles
    var_id_snow_ds         = var_id_nsnow + 1,                                &
    var_id_snow_ice        = var_id_snow_ds + ntiles_max,                     &
    var_id_snow_liq        = var_id_snow_ice + ntiles_max,                    &
    var_id_tsnow           = var_id_snow_liq + ntiles_max,                    &
    var_id_rgrainl         = var_id_tsnow + ntiles_max,                       &
    var_id_pstar           = var_id_rgrainl + ntiles_max,                     &
    var_id_q               = var_id_pstar + 1,                                &
    var_id_t               = var_id_q + 1,                                    &
    var_id_rad_net         = var_id_t + 1,                                    &
    var_id_lw_net          = var_id_rad_net + 1,                              &
    var_id_sw_net          = var_id_lw_net + 1,                               &
    var_id_lw_down         = var_id_sw_net + 1,                               &
    var_id_sw_down         = var_id_lw_down + 1,                              &
    var_id_diff_rad        = var_id_sw_down + 1,                              &
    var_id_precip          = var_id_diff_rad + 1,                             &
    var_id_tot_rain        = var_id_precip + 1,                               &
    var_id_tot_snow        = var_id_tot_rain + 1,                             &
    var_id_con_rain        = var_id_tot_snow + 1,                             &
    var_id_ls_rain         = var_id_con_rain + 1,                             &
    var_id_con_snow        = var_id_ls_rain + 1,                              &
    var_id_ls_snow         = var_id_con_snow + 1,                             &
    var_id_wind            = var_id_ls_snow + 1,                              &
    var_id_u               = var_id_wind + 1,                                 &
    var_id_v               = var_id_u + 1,                                    &
    var_id_albedo_land     = var_id_v + 1,                                    &
    var_id_canopy_gb       = var_id_albedo_land + 1,                          &
    var_id_cs_gb           = var_id_canopy_gb + 1,                            &
    var_id_depth_frozen    = var_id_cs_gb + 1,                                &
    var_id_depth_unfrozen  = var_id_depth_frozen + 1,                         &
    var_id_drain           = var_id_depth_unfrozen + 1,                       &
    var_id_elake           = var_id_drain + 1,                                &
    var_id_emis_gb         = var_id_elake + 1,                                &
    var_id_fch4_wetl       = var_id_emis_gb + 1,                              &
    var_id_fsat            = var_id_fch4_wetl + 1,                            &
    var_id_fsmc_gb         = var_id_fsat + 1,                                 &
    var_id_fwetl           = var_id_fsmc_gb + 1,                              &
    var_id_gpp_gb          = var_id_fwetl + 1,                                &
    var_id_hf_snow_melt    = var_id_gpp_gb + 1,                               &
    var_id_land_index      = var_id_hf_snow_melt + 1,                         &
    var_id_lice_index      = var_id_land_index + 1,                           &
    var_id_lit_c_mean      = var_id_lice_index + 1,                           &
    var_id_lw_up           = var_id_lit_c_mean + 1,                           &
    var_id_npp_gb          = var_id_lw_up + 1,                                &
    var_id_qbase           = var_id_npp_gb + 1,                               &
    var_id_qbase_zw        = var_id_qbase + 1,                                &
    var_id_resp_p_gb       = var_id_qbase_zw + 1,                             &
    var_id_resp_s_gb       = var_id_resp_p_gb + 1,                            &
    var_id_resp_s_dr_out   = var_id_resp_s_gb + 1,                            &
    var_id_runoff          = var_id_resp_s_dr_out + 1,                        &
    var_id_sat_excess_roff = var_id_runoff + 1,                               &
    var_id_smc_avail_top   = var_id_sat_excess_roff + 1,                      &
    var_id_smc_avail_tot   = var_id_smc_avail_top + 1,                        &
    var_id_smc_tot         = var_id_smc_avail_tot + 1,                        &
    var_id_snomlt_sub_htf  = var_id_smc_tot + 1,                              &
    var_id_snow_can_gb     = var_id_snomlt_sub_htf + 1,                       &
    var_id_snow_depth_gb   = var_id_snow_can_gb + 1,                          &
    var_id_snow_frac       = var_id_snow_depth_gb + 1,                        &
    var_id_snow_frac_alb   = var_id_snow_frac + 1,                            &
    var_id_snow_grnd_gb    = var_id_snow_frac_alb + 1,                        &
    var_id_snow_ice_gb     = var_id_snow_grnd_gb + 1,                         &
    var_id_snow_liq_gb     = var_id_snow_ice_gb + 1,                          &
    var_id_snow_melt_gb    = var_id_snow_liq_gb + 1,                          &
    var_id_soil_index      = var_id_snow_melt_gb + 1,                         &
    var_id_sub_surf_roff   = var_id_soil_index + 1,                           &
    var_id_surf_roff       = var_id_sub_surf_roff + 1,                        &
    var_id_swet_liq_tot    = var_id_surf_roff + 1,                            &
    var_id_swet_tot        = var_id_swet_liq_tot + 1,                         &
    var_id_tfall           = var_id_swet_tot + 1,                             &
    var_id_trad            = var_id_tfall + 1,                                &
    var_id_c_veg           = var_id_trad + 1,                                 &
    var_id_flux_o3_stom    = var_id_c_veg + 1,                                &
    var_id_fsmc            = var_id_flux_o3_stom + 1,                         &
    var_id_g_leaf          = var_id_fsmc + 1,                                 &
    var_id_g_leaf_day      = var_id_g_leaf + 1,                               &
    var_id_g_leaf_dr_out   = var_id_g_leaf_day + 1,                           &
    var_id_g_leaf_phen     = var_id_g_leaf_dr_out + 1,                        &
    var_id_gpp             = var_id_g_leaf_phen + 1,                          &
    var_id_lai_phen        = var_id_gpp + 1,                                  &
    var_id_lit_c           = var_id_lai_phen + 1,                             &
    var_id_npp_dr_out      = var_id_lit_c + 1,                                &
    var_id_npp             = var_id_npp_dr_out + 1,                           &
    var_id_o3_exp_fac      = var_id_npp + 1,                                  &
    var_id_resp_p          = var_id_o3_exp_fac + 1,                           &
    var_id_resp_w_dr_out   = var_id_resp_p + 1,                               &
    var_id_resp_w          = var_id_resp_w_dr_out + 1,                        &
    var_id_resp_s          = var_id_resp_w + 1,                               &
    var_id_cosz            = var_id_resp_s + 1,                               &
    var_id_diff_frac       = var_id_cosz + 1,                                 &
    var_id_ecan_gb         = var_id_diff_frac + 1,                            &
    var_id_ei_gb           = var_id_ecan_gb + 1,                              &
    var_id_esoil_gb        = var_id_ei_gb + 1,                                &
    var_id_fqw_gb          = var_id_esoil_gb + 1,                             &
    var_id_ftl_gb          = var_id_fqw_gb + 1,                               &
    var_id_land_albedo_1   = var_id_ftl_gb + 1,                               &
    var_id_land_albedo_2   = var_id_land_albedo_1 + 1,                        &
    var_id_land_albedo_3   = var_id_land_albedo_2 + 1,                        &
    var_id_land_albedo_4   = var_id_land_albedo_3 + 1,                        &
    var_id_latent_heat     = var_id_land_albedo_4 + 1,                        &
    var_id_q1p5m_gb        = var_id_latent_heat + 1,                          &
    var_id_qw1             = var_id_q1p5m_gb + 1,                             &
    var_id_rainfall        = var_id_qw1 + 1,                                  &
    var_id_snomlt_surf_htf = var_id_rainfall + 1,                             &
    var_id_snowfall        = var_id_snomlt_surf_htf + 1,                      &
    var_id_snow_mass_gb    = var_id_snowfall + 1,                             &
    var_id_surf_ht_flux_gb = var_id_snow_mass_gb + 1,                         &
    var_id_t1p5m_gb        = var_id_surf_ht_flux_gb + 1,                      &
    var_id_taux1           = var_id_t1p5m_gb + 1,                             &
    var_id_tauy1           = var_id_taux1 + 1,                                &
    var_id_tl1             = var_id_tauy1 + 1,                                &
    var_id_tstar_gb        = var_id_tl1 + 1,                                  &
    var_id_u1              = var_id_tstar_gb + 1,                             &
    var_id_u10m            = var_id_u1 + 1,                                   &
    var_id_v1              = var_id_u10m + 1,                                 &
    var_id_v10m            = var_id_v1 + 1,                                   &
    var_id_ext             = var_id_v10m + 1,                                 &
    var_id_smcl            = var_id_ext + 1,                                  &
    var_id_soil_wet        = var_id_smcl + 1,                                 &
    var_id_sthf            = var_id_soil_wet + 1,                             &
    var_id_sthu            = var_id_sthf + 1,                                 &
    var_id_alb_tile_1      = var_id_sthu + 1,                                 &
    var_id_alb_tile_2      = var_id_alb_tile_1 + 1,                           &
    var_id_alb_tile_3      = var_id_alb_tile_2 + 1,                           &
    var_id_alb_tile_4      = var_id_alb_tile_3 + 1,                           &
    var_id_anthrop_heat    = var_id_alb_tile_4 + 1,                           &
    var_id_catch           = var_id_anthrop_heat + 1,                         &
    var_id_ecan            = var_id_catch + 1,                                &
    var_id_ei              = var_id_ecan + 1,                                 &
    var_id_emis            = var_id_ei + 1,                                   &
    var_id_esoil           = var_id_emis + 1,                                 &
    var_id_fqw             = var_id_esoil + 1,                                &
    var_id_ftl             = var_id_fqw + 1,                                  &
    var_id_gc              = var_id_ftl + 1,                                  &
    var_id_le              = var_id_gc + 1,                                   &
    var_id_q1p5m           = var_id_le + 1,                                   &
    var_id_rad_net_tile    = var_id_q1p5m + 1,                                &
    var_id_snow_can_melt   = var_id_rad_net_tile + 1,                         &
    var_id_snow_can        = var_id_snow_can_melt + 1,                        &
    var_id_snow_grnd_rho   = var_id_snow_can + 1,                             &
    var_id_snow_ground     = var_id_snow_grnd_rho + 1,                        &
    var_id_snow_ice_tile   = var_id_snow_ground + 1,                          &
    var_id_snow_liq_tile   = var_id_snow_ice_tile + 1,                        &
    var_id_snow_mass       = var_id_snow_liq_tile + 1,                        &
    var_id_snow_melt       = var_id_snow_mass + 1,                            &
    var_id_surf_ht_flux    = var_id_snow_melt + 1,                            &
    var_id_surf_ht_store   = var_id_surf_ht_flux + 1,                         &
    var_id_t1p5m           = var_id_surf_ht_store + 1,                        &
    var_id_tstar           = var_id_t1p5m + 1,                                &
    var_id_z0              = var_id_tstar + 1,                                &
    var_id_tile_index      = var_id_z0 + 1,                                   &
    var_id_ozone           = var_id_tile_index + 1,                           &
    var_id_isoprene_gb     = var_id_ozone + 1,                                &
    var_id_isoprene        = var_id_isoprene_gb + 1,                          &
    var_id_terpene_gb      = var_id_isoprene + 1,                             &
    var_id_terpene         = var_id_terpene_gb + 1,                           &
    var_id_methanol_gb     = var_id_terpene + 1,                              &
    var_id_methanol        = var_id_methanol_gb + 1,                          &
    var_id_acetone_gb      = var_id_methanol + 1,                             &
    var_id_acetone         = var_id_acetone_gb + 1, &
! CABLE    
    !snow layer vars only a CONST
    !var_id_isnow_flg3l= var_id_acetone + ntiles_max 
    var_id_isnow_flg3l= var_id_acetone + 1,& 
    var_id_snow_rho1l = var_id_isnow_flg3l + 1,& 
    var_id_snage_tile = var_id_snow_rho1l + 1,& 
    var_id_snow_rho1 = var_id_snage_tile + 1, & 
    var_id_snow_rho2 = var_id_snow_rho1 + 1, & 
    var_id_snow_rho3 = var_id_snow_rho2 + 1, & 
    var_id_snow_mass1 = var_id_snow_rho3 + 1, & 
    var_id_snow_mass2 = var_id_snow_mass1 + 1, & 
    var_id_snow_mass3 = var_id_snow_mass2 + 1, & 
    var_id_snow_depth1 = var_id_snow_mass3 + 1, & 
    var_id_snow_depth2 = var_id_snow_depth1 + 1, & 
    var_id_snow_depth3 = var_id_snow_depth2 + 1, & 
    var_id_snow_tmp1 = var_id_snow_depth3 + 1, & 
    var_id_snow_tmp2 = var_id_snow_tmp1 + 1, & 
    var_id_snow_tmp3 = var_id_snow_tmp2 + 1,                                  & 
    var_id_sthu_tile1 = var_id_snow_tmp3 + 1,                                 &
    var_id_sthu_tile2 = var_id_sthu_tile1 + 1,                                & 
    var_id_sthu_tile3 = var_id_sthu_tile2 + 1,                                & 
    var_id_sthu_tile4 = var_id_sthu_tile3 + 1,                                & 
    var_id_sthu_tile5 = var_id_sthu_tile4 + 1,                                & 
    var_id_sthu_tile6 = var_id_sthu_tile5 + 1,                                & 
    var_id_sthf_tile1 = var_id_sthu_tile6 + 1,                                & 
    var_id_sthf_tile2 = var_id_sthf_tile1 + 1,                                & 
    var_id_sthf_tile3 = var_id_sthf_tile2 + 1,                                & 
    var_id_sthf_tile4 = var_id_sthf_tile3 + 1,                                & 
    var_id_sthf_tile5 = var_id_sthf_tile4 + 1,                                & 
    var_id_sthf_tile6 = var_id_sthf_tile5 + 1,                                & 
    var_id_smcl_tile1 = var_id_sthf_tile6 + 1,                                & 
    var_id_smcl_tile2 = var_id_smcl_tile1 + 1,                                & 
    var_id_smcl_tile3 = var_id_smcl_tile2 + 1,                                & 
    var_id_smcl_tile4 = var_id_smcl_tile3 + 1,                                & 
    var_id_smcl_tile5 = var_id_smcl_tile4 + 1,                                & 
    var_id_smcl_tile6 = var_id_smcl_tile5 + 1,                                & 
    var_id_tsoil_tile1 = var_id_smcl_tile6 + 1,                               & 
    var_id_tsoil_tile2 = var_id_tsoil_tile1 + 1,                              & 
    var_id_tsoil_tile3 = var_id_tsoil_tile2 + 1,                              & 
    var_id_tsoil_tile4 = var_id_tsoil_tile3 + 1,                              & 
    var_id_tsoil_tile5 = var_id_tsoil_tile4 + 1,                              & 
    var_id_tsoil_tile6 = var_id_tsoil_tile5 + 1



  INTEGER, PARAMETER :: IDENTIFIER_LEN = 20  ! Length for identifiers of model
                                             ! variables

! The name of each levels dimension in output files
  CHARACTER(len=MAX_SDF_NAME_LEN), PARAMETER ::                               &
    PFT_DIM_NAME_OUT    = 'pft',                                              &
    NVG_DIM_NAME_OUT    = 'nvg',                                              &
    TYPE_DIM_NAME_OUT   = 'type',                                             &
    TILE_DIM_NAME_OUT   = 'tile',                                             &
    SNOW_DIM_NAME_OUT   = 'snow',                                             &
    SOIL_DIM_NAME_OUT   = 'soil',                                             &
    SCPOOL_DIM_NAME_OUT = 'scpool'


!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
! The name of each levels dimension in input files
  CHARACTER(len=MAX_SDF_NAME_LEN) ::                                          &
    pft_dim_name    = 'pft',                                                  &
    nvg_dim_name    = 'nvg',                                                  &
    type_dim_name   = 'type',                                                 &
    tile_dim_name   = 'tile',                                                 &
    snow_dim_name   = 'snow',                                                 &
    soil_dim_name   = 'soil',                                                 &
    scpool_dim_name = 'scpool'

! The size of each possible levels dimension
  INTEGER ::                                                                  &
    surface_dim_size = 1,                                                     &
        ! 'Surface' variables have a levels dimension of size 1
    pft_dim_size     = 5,                                                     &
    nvg_dim_size     = 4,                                                     &
    type_dim_size    = 9,                                                     &
    tile_dim_size    = 9,                                                     &
    snow_dim_size    = 0,                                                     &
    soil_dim_size    = 4,                                                     &
    scpool_dim_size  = 1


CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "get_var_id.inc"
#include "get_string_identifier.inc"
#include "get_var_levs_dim.inc"
#include "get_var_attrs.inc"
#include "extract_var.inc"
#include "populate_var.inc"

#include "map_to_land.inc"
#include "map_from_land.inc"
#include "tiles_to_gbm.inc"

END MODULE model_interface_mod
#endif
