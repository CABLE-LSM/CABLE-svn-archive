#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE get_var_attrs(var_id, int_attributes, real_attributes, char_attributes)

  USE io_constants, ONLY : MAX_ATTR_VAR, MDI

  USE dictionary_mod, ONLY : dict, dict_create, dict_set

  USE string_utils_mod, ONLY : to_string

  USE switches, ONLY : l_triffid

  USE soil_param, ONLY : zsmc

  USE model_grid_mod, ONLY : latitude, longitude

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an identifer for a model variable, returns the attributes that
!   should be defined on output files
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  INTEGER, INTENT(IN) :: var_id  ! Identifies the model variable
  TYPE(dict), INTENT(OUT) :: int_attributes  ! Attributes with integer values
  TYPE(dict), INTENT(OUT) :: real_attributes  ! Attributes with real values
  TYPE(dict), INTENT(OUT) :: char_attributes  ! Attributes with character values


!-----------------------------------------------------------------------------

! Create the dictionary objects
  int_attributes = dict_create(MAX_ATTR_VAR, INT(1))
  real_attributes = dict_create(MAX_ATTR_VAR, REAL(1.0))
  char_attributes = dict_create(MAX_ATTR_VAR, "char")

! All variables except lat and lon will have a missing data attribute and
! a coordinate attribute
! lat and lon do not need these since they are interpreted as coordinate variables
! as per the CF convention
  IF ( var_id /= var_id_latitude .AND. var_id /= var_id_longitude ) THEN
    CALL dict_set(real_attributes, "missing_value", MDI)
    CALL dict_set(real_attributes, "_FillValue", MDI)
    CALL dict_set(char_attributes, "coordinates", "latitude longitude")
  END IF

  SELECT CASE ( var_id )
    CASE ( var_id_latitude )
      CALL dict_set(char_attributes, "standard_name", "latitude")
      CALL dict_set(char_attributes, "long_name", "Gridbox latitude")
      CALL dict_set(char_attributes, "units", "degrees")
      CALL dict_set(real_attributes, "valid_min", MINVAL(latitude))
      CALL dict_set(real_attributes, "valid_max", MAXVAL(latitude))

    CASE ( var_id_longitude )
      CALL dict_set(char_attributes, "standard_name", "longitude")
      CALL dict_set(char_attributes, "long_name", "Gridbox longitude")
      CALL dict_set(char_attributes, "units", "degrees")
      CALL dict_set(real_attributes, "valid_min", MINVAL(longitude))
      CALL dict_set(real_attributes, "valid_max", MAXVAL(longitude))

    CASE ( var_id_frac )
      CALL dict_set(char_attributes, "long_name", "Fractional cover of each surface type")
      CALL dict_set(char_attributes, "units", "1")
      CALL dict_set(real_attributes, "valid_min", 0.0)
      CALL dict_set(real_attributes, "valid_max", 1.0)

    CASE ( var_id_b )
      CALL dict_set(char_attributes, "long_name", "Gridbox Brooks-Corey exponent for each soil layer")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_sathh )
      CALL dict_set(char_attributes, "long_name", "Gridbox saturated soil water pressure for each soil layer")
      CALL dict_set(char_attributes, "units", "m")

    CASE ( var_id_satcon )
      CALL dict_set(char_attributes, "long_name", "Gridbox saturated hydraulic conductivity for each soil layer")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_sm_sat )
      CALL dict_set(char_attributes, "long_name", "Gridbox volumetric moisture content at saturation for each soil layer")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_sm_crit )
      CALL dict_set(char_attributes, "long_name", "Gridbox volumetric moisture content at critical point for each soil layer")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_sm_wilt )
      CALL dict_set(char_attributes, "long_name", "Gridbox  volumetric moisture content at wilting point for each soil layer")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_hcap )
      CALL dict_set(char_attributes, "long_name", "Gridbox dry soil heat capacity for each soil layer")
      CALL dict_set(char_attributes, "units", "J K-1 m-3")

    CASE ( var_id_hcon )
      CALL dict_set(char_attributes, "long_name", "Gridbox dry soil thermal conductivity for each soil layer")
      CALL dict_set(char_attributes, "units", "W m-1 K-1")

    CASE ( var_id_albedo_land )
      CALL dict_set(char_attributes, "long_name", "Gridbox albedo (as used for net shortwave calculation)")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_canopy_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox canopy water content")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_cs_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox soil carbon (total)")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_cv )
      CALL dict_set(char_attributes, "long_name", "Gridbox mean vegetation carbon at end of model timestep.")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_depth_frozen )
      CALL dict_set(char_attributes, "long_name", "Gridbox depth of frozen ground at surface")
      CALL dict_set(char_attributes, "units", "m")

    CASE ( var_id_depth_unfrozen )
      CALL dict_set(char_attributes, "long_name", "Gridbox depth of unfrozen ground at surface")
      CALL dict_set(char_attributes, "units", "m")

    CASE ( var_id_drain )
      CALL dict_set(char_attributes, "long_name", "Drainage from bottom (nshyd) soil layer")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_elake )
      CALL dict_set(char_attributes, "long_name", "Gridbox mean evaporation from lakes")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_emis_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox emissivity")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_fch4_wetl )
      CALL dict_set(char_attributes, "long_name", "Scaled methane flux from wetland fraction")
      CALL dict_set(char_attributes, "units", "10^-9 kg m-2 s-1")

    CASE ( var_id_fsat )
      CALL dict_set(char_attributes, "long_name", "Surface saturated fraction")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_fsmc_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox soil moisture availability factor (beta)")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_fwetl )
      CALL dict_set(char_attributes, "long_name", "Wetland fraction at end of model timestep.")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_gpp_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox gross primary productivity")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_gs )
      CALL dict_set(char_attributes, "long_name", "Gridbox surface conductance to evaporation")
      CALL dict_set(char_attributes, "units", "m s-1")

    CASE ( var_id_hf_snow_melt )
      CALL dict_set(char_attributes, "long_name", "Gridbox snowmelt heat flux")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_land_index )
      CALL dict_set(char_attributes, "long_name", "Index (gridbox number) of land points")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_lice_index )
      CALL dict_set(char_attributes, "long_name", "Index (gridbox number) of land ice points")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_lit_c_mean )
      CALL dict_set(char_attributes, "long_name", "Gridbox mean carbon litter")
      CALL dict_set(char_attributes, "units", "kg m-2 per 360days")

    CASE ( var_id_lw_net )
      CALL dict_set(char_attributes, "long_name", "Gridbox net downward longwave radiation at surface")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_lw_up )
      CALL dict_set(char_attributes, "long_name", "Gridbox surface upward LW radiation of land points")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_npp_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox net primary productivity")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_qbase )
      CALL dict_set(char_attributes, "long_name", "Baseflow (lateral subsurface runoff)")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_qbase_zw )
      CALL dict_set(char_attributes, "long_name", "Baseflow from deep layer")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_rad_net )
      CALL dict_set(char_attributes, "long_name", "Surface net radiation of land points")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_resp_p_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox plant respiration")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_resp_s_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox soil respiration (total)")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_resp_s_dr_out )
      CALL dict_set(char_attributes, "long_name", "Gridbox mean soil respiration for driving TRIFFID")
      CALL dict_set(char_attributes, "units", "kg  m-2 per 360days")

    CASE ( var_id_runoff )
      CALL dict_set(char_attributes, "long_name", "Gridbox runoff rate")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_sat_excess_roff )
      CALL dict_set(char_attributes, "long_name", 'Saturation excess surface ("Dunne") runoff')
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_smc_avail_top )
      CALL dict_set(char_attributes, "long_name",                             &
                    "Gridbox available moisture in top " //                   &
                    TRIM(to_string(zsmc)) // "m of soil")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_smc_avail_tot )
      CALL dict_set(char_attributes, "long_name", "Gridbox available moisture in soil column")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_smc_tot )
      CALL dict_set(char_attributes, "long_name", "Gridbox total soil moisture in column")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snomlt_sub_htf )
      CALL dict_set(char_attributes, "long_name", "Gridbox sub-canopy snowmelt heat flux")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_snow_can_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox snow on canopy")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snow_depth_gb )
      CALL dict_set(char_attributes, "long_name", "Snow depth (on ground)")
      CALL dict_set(char_attributes, "units", "m")

    CASE ( var_id_snow_frac )
      CALL dict_set(char_attributes, "long_name", "Gridbox snow-covered fraction of land points")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_snow_grnd_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox snow below canopy (snow_grnd)")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snow_ice_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox ice content of snow on ground")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snow_liq_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox liquid content of snow on ground")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snow_melt_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox rate of snowmelt")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_soil_index )
      CALL dict_set(char_attributes, "long_name", "Index (gridbox number) of soil points")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_sthzw )
      CALL dict_set(char_attributes, "long_name", "Soil wetness in deep (water table/TOPMODEL) layer")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_sub_surf_roff )
      CALL dict_set(char_attributes, "long_name", "Gridbox sub-surface runoff")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_surf_roff )
      CALL dict_set(char_attributes, "long_name", "Gridbox surface runoff")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

!    CASE ( var_id_surf_roff_inf )
!      CALL dict_set(char_attributes, "long_name", "Gridbox infiltration excess surface runoff")
!      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_swet_liq_tot )
      CALL dict_set(char_attributes, "long_name", "Gridbox unfrozen soil moisture as fraction of saturation")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_swet_tot )
      CALL dict_set(char_attributes, "long_name", "Gridbox soil moisture as fraction of saturation")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_sw_net )
      CALL dict_set(char_attributes, "long_name", "Gridbox net downward shortwave radiation at surface")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_tfall )
      CALL dict_set(char_attributes, "long_name", "Gridbox throughfall")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_trad )
      CALL dict_set(char_attributes, "long_name", "Gridbox effective radiative temperature (assuming emissivity=1)")
      CALL dict_set(char_attributes, "units", "K")

!    CASE ( var_id_w_flux_sfc )
!      CALL dict_set(char_attributes, "long_name", "Gridbox downwards moisture flux at soil surface")
!      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_zw )
      CALL dict_set(char_attributes, "long_name", "Gridbox mean depth to water table")
      CALL dict_set(char_attributes, "units", "m")

    CASE ( var_id_c_veg )
      CALL dict_set(char_attributes, "long_name", "PFT total carbon content of the vegetation at the end of model timestep.")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_canht )
      CALL dict_set(char_attributes, "long_name", "PFT canopy height")
      CALL dict_set(char_attributes, "units", "m")

!    CASE ( var_id_ci )
!      CALL dict_set(char_attributes, "long_name", "PFT internal CO2 pressure")
!      CALL dict_set(char_attributes, "units", "Pa")

    CASE ( var_id_flux_o3_stom )
      CALL dict_set(char_attributes, "long_name", "Flux of O3 to stomata")
      CALL dict_set(char_attributes, "units", "mol m-2 s-1")

    CASE ( var_id_fsmc )
      CALL dict_set(char_attributes, "long_name", "PFT soil moisture availability factor (beta)")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_g_leaf )
      CALL dict_set(char_attributes, "long_name", "PFT leaf turnover rate")
      CALL dict_set(char_attributes, "units", "per 360days")

    CASE ( var_id_g_leaf_day )
      CALL dict_set(char_attributes, "long_name", "PFT mean leaf turnover rate for input to PHENOL")
      CALL dict_set(char_attributes, "units", "per 360days")

    CASE ( var_id_g_leaf_dr_out )
      CALL dict_set(char_attributes, "long_name", "PFT mean leaf turnover rate for driving TRIFFID")
      CALL dict_set(char_attributes, "units", "per 360days")

    CASE ( var_id_g_leaf_phen )
      CALL dict_set(char_attributes, "long_name", "PFT mean leaf turnover rate over phenology period")
      CALL dict_set(char_attributes, "units", "per 360days")

!    CASE ( var_id_g_stom )
!      CALL dict_set(char_attributes, "long_name", "PFT bulk (canopy) stomatal conductance for water vapour")
!      CALL dict_set(char_attributes, "units", "m s-1")

    CASE ( var_id_gpp )
      CALL dict_set(char_attributes, "long_name", "PFT gross primary productivity")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_lai )
      CALL dict_set(char_attributes, "long_name", "PFT leaf area index")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_lai_phen )
      CALL dict_set(char_attributes, "long_name", "PFT leaf area index after phenology")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_lit_c )
      CALL dict_set(char_attributes, "long_name", "PFT carbon Litter")
      CALL dict_set(char_attributes, "units", "kg  m-2 per 360days")

    CASE ( var_id_npp_dr_out )
      CALL dict_set(char_attributes, "long_name", "PFT mean NPP for driving TRIFFID")
      CALL dict_set(char_attributes, "units", "kg  m-2 per 360days")

    CASE ( var_id_npp )
      CALL dict_set(char_attributes, "long_name", "PFT net primary productivity")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_o3_exp_fac )
      CALL dict_set(char_attributes, "long_name", "Ozone exposure factor")
      CALL dict_set(char_attributes, "units", "1")

!    CASE ( var_id_rdc )
!      CALL dict_set(char_attributes, "long_name", "Canopy dark respiration, without soil water dependence")
!      CALL dict_set(char_attributes, "units", "molCO2 m-2 s-1")

    CASE ( var_id_resp_p )
      CALL dict_set(char_attributes, "long_name", "PFT plant respiration")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_resp_w_dr_out )
      CALL dict_set(char_attributes, "long_name", "PFT mean wood respiration for driving TRIFFID")
      CALL dict_set(char_attributes, "units", "kg C m-2 per 360days")

    CASE ( var_id_resp_w )
      CALL dict_set(char_attributes, "long_name", "PFT wood maintenance respiration")
      CALL dict_set(char_attributes, "units", "kg C m-2 s-1")

    CASE ( var_id_cs )
      IF ( l_triffid ) THEN
        CALL dict_set(char_attributes, "long_name", "Gridbox soil carbon in each pool (DPM,RPM,bio,hum)")
      ELSE
        CALL dict_set(char_attributes, "long_name", "Gridbox soil carbon (single pool)")
      END IF
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_resp_s )
      IF ( l_triffid ) THEN
        CALL dict_set(char_attributes, "long_name", "Gridbox soil respiration from each pool (DPM,RPM,bio,hum)")
      ELSE
        CALL dict_set(char_attributes, "long_name", "Gridbox soil respiration (single pool)")
      END IF
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_con_rain )
      CALL dict_set(char_attributes, "long_name", "Gridbox convective rainfall")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_con_snow )
      CALL dict_set(char_attributes, "long_name", "Gridbox convective snowfall")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_cosz )
      CALL dict_set(char_attributes, "long_name", "Cosine of the zenith angle")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_diff_frac )
      CALL dict_set(char_attributes, "long_name", "Gridbox fraction of radiation that is diffuse")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_ecan_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox mean evaporation from canopy/surface store")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_ei_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox sublimation from lying snow or sea-ice")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_esoil_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox surface evapotranspiration from soil moisture store")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_fqw_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox moisture flux from surface")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_ftl_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox surface sensible heat flux")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_land_albedo_1 )
      CALL dict_set(char_attributes, "long_name", "Gridbox albedo for waveband 1")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_land_albedo_2 )
      CALL dict_set(char_attributes, "long_name", "Gridbox albedo for waveband 2")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_land_albedo_3 )
      CALL dict_set(char_attributes, "long_name", "Gridbox albedo for waveband 3")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_land_albedo_4 )
      CALL dict_set(char_attributes, "long_name", "Gridbox albedo for waveband 4")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_latent_heat )
      CALL dict_set(char_attributes, "long_name", "Gridbox surface latent heat flux")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_ls_rain )
      CALL dict_set(char_attributes, "long_name", "Gridbox large-scale rainfall")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_ls_snow )
      CALL dict_set(char_attributes, "long_name", "Gridbox large-scale snowfall")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_lw_down )
      CALL dict_set(char_attributes, "long_name", "Gridbox surface downward LW radiation")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_precip )
      CALL dict_set(char_attributes, "long_name", "Gridbox precipitation rate")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_pstar )
      CALL dict_set(char_attributes, "long_name", "Gridbox surface pressure")
      CALL dict_set(char_attributes, "units", "Pa")

    CASE ( var_id_q1p5m_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox specific humidity at 1.5m height")
      CALL dict_set(char_attributes, "units", "kg kg-1")

    CASE ( var_id_qw1 )
      CALL dict_set(char_attributes, "long_name", "Gridbox specific humidity (total water content)")
      CALL dict_set(char_attributes, "units", "kg kg-1")

    CASE ( var_id_rainfall )
      CALL dict_set(char_attributes, "long_name", "Gridbox rainfall rate")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_snomlt_surf_htf )
      CALL dict_set(char_attributes, "long_name", "Gridbox heat flux used for surface melting of snow")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_snowfall )
      CALL dict_set(char_attributes, "long_name", "Gridbox snowfall rate")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_snow_mass_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox snowmass (on land this is lying_snow=snow_tile+snow_grnd)")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_surf_ht_flux_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox net downward heat flux at surface over land and sea-ice fraction")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_sw_down )
      CALL dict_set(char_attributes, "long_name", "Gridbox surface downward SW radiation")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_t1p5m_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox temperature at 1.5m height")
      CALL dict_set(char_attributes, "units", "K")

    CASE ( var_id_taux1 )
      CALL dict_set(char_attributes, "long_name", "Gridbox westerly component of surface wind stress")
      CALL dict_set(char_attributes, "units", "N m-2")

    CASE ( var_id_tauy1 )
      CALL dict_set(char_attributes, "long_name", "Gridbox southerly component of surface wind stress")
      CALL dict_set(char_attributes, "units", "N m-2")

    CASE ( var_id_tl1 )
      CALL dict_set(char_attributes, "long_name", "Gridbox ice/liquid water temperature")
      CALL dict_set(char_attributes, "units", "K")

    CASE ( var_id_tstar_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox surface temperature")
      CALL dict_set(char_attributes, "units", "K")

    CASE ( var_id_u1 )
      CALL dict_set(char_attributes, "long_name", "Gridbox westerly wind component")
      CALL dict_set(char_attributes, "units", "m s-1")

    CASE ( var_id_u10m )
      CALL dict_set(char_attributes, "long_name", "Gridbox westerly wind component at 10 m height")
      CALL dict_set(char_attributes, "units", "m s-1")

    CASE ( var_id_v1 )
      CALL dict_set(char_attributes, "long_name", "Gridbox southerly wind component")
      CALL dict_set(char_attributes, "units", "m s-1")

    CASE ( var_id_v10m )
      CALL dict_set(char_attributes, "long_name", "Gridbox southerly wind component at 10 m height")
      CALL dict_set(char_attributes, "units", "m s-1")

    CASE ( var_id_wind )
      CALL dict_set(char_attributes, "long_name", "Gridbox wind speed")
      CALL dict_set(char_attributes, "units", "m s-1")

    CASE ( var_id_ext )
      CALL dict_set(char_attributes, "long_name", "Gridbox extraction of water from each soil layer")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_smcl )
      CALL dict_set(char_attributes, "long_name", "Gridbox moisture content of each soil layer")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_soil_wet )
      CALL dict_set(char_attributes, "long_name",                             &
                                     "Gridbox total moisture content of each soil layer, as fraction of saturation")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_sthf )
      CALL dict_set(char_attributes, "long_name",                             &
                                     "Gridbox frozen moisture content of each soil layer as a fraction of saturation")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_sthu )
      CALL dict_set(char_attributes, "long_name",                             &
                    "Gridbox unfrozen moisture content of each soil layer as a fraction of saturation")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_t_soil )
      CALL dict_set(char_attributes, "long_name", "Gridbox sub-surface temperature of each layer")
      CALL dict_set(char_attributes, "units", "K")

!    CASE ( var_id_w_flux )
!      CALL dict_set(char_attributes, "long_name", "Downwards moisture flux at bottom of each soil layer")
!      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_alb_tile_1 )
      CALL dict_set(char_attributes, "long_name", "Tile land albedo, waveband 1")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_alb_tile_2 )
      CALL dict_set(char_attributes, "long_name", "Tile land albedo, waveband 2")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_alb_tile_3 )
      CALL dict_set(char_attributes, "long_name", "Tile land albedo, waveband 3")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_alb_tile_4 )
      CALL dict_set(char_attributes, "long_name", "Tile land albedo, waveband 4")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_anthrop_heat )
      CALL dict_set(char_attributes, "long_name", "Anthropogenic heat flux for each tile")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_canopy )
      CALL dict_set(char_attributes, "long_name", "Tile surface/canopy water for snow-free land tiles")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_catch )
      CALL dict_set(char_attributes, "long_name", "Tile surface/canopy water capacity of snow-free land tiles")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_ecan )
      CALL dict_set(char_attributes, "long_name", "Tile evaporation from canopy/surface store for snow-free land tiles")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_ei )
      CALL dict_set(char_attributes, "long_name", "Tile sublimation from lying snow for land tiles")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_emis )
      CALL dict_set(char_attributes, "long_name", "Tile emissivity")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_esoil )
      CALL dict_set(char_attributes, "long_name",                             &
                                     "Tile surface evapotranspiration from soil moisture store for snow-free land tiles")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_fqw )
      CALL dict_set(char_attributes, "long_name", "Tile surface moisture flux for land tiles")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_ftl )
      CALL dict_set(char_attributes, "long_name", "Tile surface sensible heat flux for land tiles")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_gc )
      CALL dict_set(char_attributes, "long_name", "Tile surface conductance to evaporation for land tiles")
      CALL dict_set(char_attributes, "units", "m s-1")

    CASE ( var_id_le )
      CALL dict_set(char_attributes, "long_name", "Tile surface latent heat flux for land tiles")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_nsnow )
      CALL dict_set(char_attributes, "long_name", "Number of snow layers")
      CALL dict_set(char_attributes, "units", "1")

    CASE ( var_id_q1p5m )
      CALL dict_set(char_attributes, "long_name", "Tile specific humidity at 1.5m over land tiles")
      CALL dict_set(char_attributes, "units", "kg kg-1")

    CASE ( var_id_rad_net_tile )
      CALL dict_set(char_attributes, "long_name", "Tile surface net radiation")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_rgrain )
      CALL dict_set(char_attributes, "long_name", "Tile snow surface grain size")
      CALL dict_set(char_attributes, "units", "microns")

    CASE ( var_id_snow_can_melt )
      CALL dict_set(char_attributes, "long_name", "Tile melt of canopy snow")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_snow_can )
      CALL dict_set(char_attributes, "long_name", "Tile canopy snow")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snow_depth )
      CALL dict_set(char_attributes, "long_name", "Tile snow depth (on ground)")
      CALL dict_set(char_attributes, "units", "m")

!    CASE ( var_id_snow_grnd_can_melt )
!      CALL dict_set(char_attributes, "long_name", "Tile melt of snow under canopy")
!      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_snow_grnd_rho )
      CALL dict_set(char_attributes, "long_name", "Tile density of snow on ground")
      CALL dict_set(char_attributes, "units", "kg m-3")

    CASE ( var_id_snow_grnd )
      CALL dict_set(char_attributes, "long_name", "Tile snow on ground below canopy (snow_grnd)")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snow_ground )
      CALL dict_set(char_attributes, "long_name", "Tile snow on ground (snow_tile or snow_grnd)")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snow_ice_tile )
      CALL dict_set(char_attributes, "long_name", "Tile total ice content in snow on ground")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snow_liq_tile )
      CALL dict_set(char_attributes, "long_name", "Tile total liquid content in snow on ground")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snow_mass )
      CALL dict_set(char_attributes, "long_name", "Tile lying snow (snow_tile+snow_grnd)")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snow_melt )
      CALL dict_set(char_attributes, "long_name", "Tile snow melt (melt_tile)")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_surf_ht_flux )
      CALL dict_set(char_attributes, "long_name", "Downward heat flux for each tile")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_surf_ht_store )
      CALL dict_set(char_attributes, "long_name", "C*(dT/dt) for each tile")
      CALL dict_set(char_attributes, "units", "W m-2")

    CASE ( var_id_t1p5m )
      CALL dict_set(char_attributes, "long_name", "Tile temperature at 1.5m over land tiles")
      CALL dict_set(char_attributes, "units", "K")

    CASE ( var_id_tstar )
      CALL dict_set(char_attributes, "long_name", "Tile surface temperature")
      CALL dict_set(char_attributes, "units", "K")

    CASE ( var_id_z0 )
      CALL dict_set(char_attributes, "long_name", "Tile surface roughness")
      CALL dict_set(char_attributes, "units", "m")

    CASE ( var_id_tile_index )
      CALL dict_set(char_attributes, "long_name", "Index (gridbox number) of land points with each surface type")
      CALL dict_set(char_attributes, "units", "m")

    CASE ( var_id_isoprene_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox mean isoprene emission flux")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_isoprene )
      CALL dict_set(char_attributes, "long_name", "PFT isoprene emission flux")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_terpene_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox mean (mono-)terpene emission flux")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_terpene )
      CALL dict_set(char_attributes, "long_name", "PFT (mono-)terpene emission flux")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_methanol_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox mean methanol emission flux")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_methanol )
      CALL dict_set(char_attributes, "long_name", "PFT methanol emission flux")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_acetone_gb )
      CALL dict_set(char_attributes, "long_name", "Gridbox mean acetone emission flux")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_acetone )
      CALL dict_set(char_attributes, "long_name", "PFT acetone emission flux")
      CALL dict_set(char_attributes, "units", "kg m-2 s-1")

    CASE ( var_id_rgrainl )
      CALL dict_set(char_attributes, "long_name", "Grain size in snow layers for each tile")
      CALL dict_set(char_attributes, "units", "microns")

    CASE ( var_id_snow_ds )
      CALL dict_set(char_attributes, "long_name", "Depth of snow layers for each tile")
      CALL dict_set(char_attributes, "units", "m")

    CASE ( var_id_snow_ice )
      CALL dict_set(char_attributes, "long_name", "Ice mass in snow layers for each tile")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_snow_liq )
      CALL dict_set(char_attributes, "long_name", "Liquid mass in snow layers for each tile")
      CALL dict_set(char_attributes, "units", "kg m-2")

    CASE ( var_id_tsnow )
      CALL dict_set(char_attributes, "long_name", "Temperature of snow layers for each tile")
      CALL dict_set(char_attributes, "units", "K")

    CASE DEFAULT
      CALL log_fatal("get_var_attrs",                                         &
                     "Unrecognised variable for output - '" //                &
                     TRIM(get_string_identifier(var_id)) // "'. " //          &
                     "See docs for available variables")
  END SELECT

  RETURN

END SUBROUTINE get_var_attrs
#endif
