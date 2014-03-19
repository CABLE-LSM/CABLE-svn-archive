#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!-----------------------------------------------------------------------------
! DATA statement for variable metadata
!-----------------------------------------------------------------------------


DATA metadata/                                                                &
!-----------------------------------------------------------------------------
! Metadata for latitude
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'latitude',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox latitude",                                                       &
! Units
    "degrees"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for longitude
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'longitude',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox longitude",                                                      &
! Units
    "degrees"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_fraction
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_fraction',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for surf_hgt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'surf_hgt',                                                               &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for frac
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'frac',                                                                   &
! Variable type
    VAR_TYPE_TYPE,                                                            &
! Long name
    "Fractional cover of each surface type",                                  &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albsoil
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albsoil',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albobs_sw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albobs_sw',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albobs_vis
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albobs_vis',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albobs_nir
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albobs_nir',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for b
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'b',                                                                      &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox Brooks-Corey exponent for each soil layer",                      &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for b_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'b_const_z',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sathh
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sathh',                                                                  &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox saturated soil water pressure for each soil layer",              &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sathh_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sathh_const_z',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for satcon
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'satcon',                                                                 &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox saturated hydraulic conductivity for each soil layer",           &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for satcon_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'satcon_const_z',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_sat
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_sat',                                                                 &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox volumetric moisture content at saturation for each soil layer",  &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_sat_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_sat_const_z',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_crit
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_crit',                                                                &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox volumetric moisture content at critical point for each soil layer", &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_crit_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_crit_const_z',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_wilt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_wilt',                                                                &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox  volumetric moisture content at wilting point for each soil layer", &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sm_wilt_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sm_wilt_const_z',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hcap
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hcap',                                                                   &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox dry soil heat capacity for each soil layer",                     &
! Units
    "J K-1 m-3"                                                               &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hcap_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hcap_const_z',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hcon
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hcon',                                                                   &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox dry soil thermal conductivity for each soil layer",              &
! Units
    "W m-1 K-1"                                                               &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hcon_const_z
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hcon_const_z',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fexp
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fexp',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ti_mean
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ti_mean',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ti_sig
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ti_sig',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for frac_agr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'frac_agr',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for wrr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'wrr',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hwr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hwr',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hgt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hgt',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ztm
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ztm',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for disp
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'disp',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albwl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albwl',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albrd
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albrd',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emisw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emisw',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emisr
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emisr',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for canopy
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'canopy',                                                                 &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile surface/canopy water for snow-free land tiles",                     &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cs
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cs',                                                                     &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name is dynamic - set in get_var_attrs
    "",                                                                       &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for gs
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'gs',                                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface conductance to evaporation",                             &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_tile
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_tile',                                                              &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sthuf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sthuf',                                                                  &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for t_soil
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    't_soil',                                                                 &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox sub-surface temperature of each layer",                          &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tstar_tile
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tstar_tile',                                                             &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lai
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lai',                                                                    &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT leaf area index",                                                    &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for canht
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'canht',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT canopy height",                                                      &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sthzw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sthzw',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Soil wetness in deep (water table/TOPMODEL) layer",                      &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for zw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'zw',                                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean depth to water table",                                      &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rgrain
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rgrain',                                                                 &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile snow surface grain size",                                           &
! Units
    "microns"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cv
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cv',                                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean vegetation carbon at end of model timestep.",               &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rho_snow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rho_snow',                                                               &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_depth
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_depth',                                                             &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile snow depth (on ground)",                                            &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_grnd
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_grnd',                                                              &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile snow on ground below canopy (snow_grnd)",                           &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for nsnow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'nsnow',                                                                  &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
   "Number of snow layers",                                                  &
! Units
   "1"                                                                       &
 ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_ds
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_ds',                                                                &
! Variable type
    VAR_TYPE_SNOW,                                                            &
! Long name
    "Depth of snow layers for each tile",                                     &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_ice
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_ice',                                                               &
! Variable type
    VAR_TYPE_SNOW,                                                            &
! Long name
    "Ice mass in snow layers for each tile",                                  &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_liq
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_liq',                                                               &
! Variable type
    VAR_TYPE_SNOW,                                                            &
! Long name
    "Liquid mass in snow layers for each tile",                               &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tsnow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tsnow',                                                                  &
! Variable type
    VAR_TYPE_SNOW,                                                            &
! Long name
    "Temperature of snow layers for each tile",                               &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rgrainl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rgrainl',                                                                &
! Variable type
    VAR_TYPE_SNOW,                                                            &
! Long name
    "Grain size in snow layers for each tile",                                &
! Units
    "microns"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for pstar
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'pstar',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface pressure",                                               &
! Units
    "Pa"                                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for q
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'q',                                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for t
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    't',                                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rad_net
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rad_net',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Surface net radiation of land points",                                   &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lw_net
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lw_net',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox net downward longwave radiation at surface",                     &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sw_net
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sw_net',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox net downward shortwave radiation at surface",                    &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lw_down
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lw_down',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface downward LW radiation",                                  &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sw_down
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sw_down',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface downward SW radiation",                                  &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for diff_rad
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'diff_rad',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for precip
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'precip',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox precipitation rate",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tot_rain
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tot_rain',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tot_snow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tot_snow',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for con_rain
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'con_rain',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox convective rainfall",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ls_rain
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ls_rain',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox large-scale rainfall",                                           &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for con_snow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'con_snow',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox convective snowfall",                                            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ls_snow
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ls_snow',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox large-scale snowfall",                                           &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for wind
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'wind',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox wind speed",                                                     &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for u
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'u',                                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for v
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'v',                                                                      &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for albedo_land
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'albedo_land',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox albedo (as used for net shortwave calculation)",                 &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for canopy_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'canopy_gb',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox canopy water content",                                           &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cs_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cs_gb',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox soil carbon (total)",                                            &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for depth_frozen
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'depth_frozen',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox depth of frozen ground at surface",                              &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for depth_unfrozen
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'depth_unfrozen',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox depth of unfrozen ground at surface",                            &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for drain
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'drain',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Drainage from bottom (nshyd) soil layer",                                &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for elake
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'elake',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean evaporation from lakes",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emis_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emis_gb',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox emissivity",                                                     &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fch4_wetl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fch4_wetl',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Scaled methane flux from wetland fraction",                              &
! Units
    "10^-9 kg m-2 s-1"                                                        &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fsat
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fsat',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Surface saturated fraction",                                             &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fsmc_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fsmc_gb',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox soil moisture availability factor (beta)",                       &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fwetl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fwetl',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Wetland fraction at end of model timestep.",                             &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for gpp_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'gpp_gb',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox gross primary productivity",                                     &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for hf_snow_melt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'hf_snow_melt',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snowmelt heat flux",                                             &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_index
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_index',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Index (gridbox number) of land points",                                  &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lice_index
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lice_index',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Index (gridbox number) of land ice points",                              &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lit_c_mean
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lit_c_mean',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean carbon litter",                                             &
! Units
    "kg m-2 per 360days"                                                      &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lw_up
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lw_up',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface upward LW radiation of land points",                     &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for npp_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'npp_gb',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox net primary productivity",                                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for qbase
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'qbase',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Baseflow (lateral subsurface runoff)",                                   &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for qbase_zw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'qbase_zw',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Baseflow from deep layer",                                               &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_p_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_p_gb',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox plant respiration",                                              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_s_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_s_gb',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox soil respiration (total)",                                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_s_dr_out
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_s_dr_out',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean soil respiration for driving TRIFFID",                      &
! Units
    "kg  m-2 per 360days"                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for runoff
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'runoff',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox runoff rate",                                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sat_excess_roff
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sat_excess_roff',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Saturation excess surface ('Dunne') runoff",                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for smc_avail_top
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'smc_avail_top',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name is dynamic - set in get_var_attrs
    "",                                                                       &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for smc_avail_tot
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'smc_avail_tot',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox available moisture in soil column",                              &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for smc_tot
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'smc_tot',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox total soil moisture in column",                                  &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snomlt_sub_htf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snomlt_sub_htf',                                                         &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox sub-canopy snowmelt heat flux",                                  &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_can_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_can_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snow on canopy",                                                 &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_depth_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_depth_gb',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Snow depth (on ground)",                                                 &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_frac
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_frac',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snow-covered fraction of land points",                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_grnd_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_grnd_gb',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snow below canopy (snow_grnd)",                                  &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_ice_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_ice_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox ice content of snow on ground",                                  &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_liq_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_liq_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox liquid content of snow on ground",                               &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_melt_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_melt_gb',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox rate of snowmelt",                                               &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for soil_index
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'soil_index',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Index (gridbox number) of soil points",                                  &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sub_surf_roff
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sub_surf_roff',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox sub-surface runoff",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for surf_roff
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'surf_roff',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface runoff",                                                 &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for swet_liq_tot
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'swet_liq_tot',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox unfrozen soil moisture as fraction of saturation",               &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for swet_tot
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'swet_tot',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox soil moisture as fraction of saturation",                        &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tfall
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tfall',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox throughfall",                                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for trad
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'trad',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox effective radiative temperature (assuming emissivity=1)",        &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for c_veg
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'c_veg',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT total carbon content of the vegetation at the end of model timestep.", &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for flux_o3_stom
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'flux_o3_stom',                                                           &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Flux of O3 to stomata",                                                  &
! Units
    "mol m-2 s-1"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fsmc
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fsmc',                                                                   &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT soil moisture availability factor (beta)",                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for g_leaf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'g_leaf',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT leaf turnover rate",                                                 &
! Units
    "per 360days"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for g_leaf_day
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'g_leaf_day',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT mean leaf turnover rate for input to PHENOL",                        &
! Units
    "per 360days"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for g_leaf_dr_out
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'g_leaf_dr_out',                                                          &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT mean leaf turnover rate for driving TRIFFID",                        &
! Units
    "per 360days"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for g_leaf_phen
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'g_leaf_phen',                                                            &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT mean leaf turnover rate over phenology period",                      &
! Units
    "per 360days"                                                             &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for gpp
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'gpp',                                                                    &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT gross primary productivity",                                         &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lai_phen
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lai_phen',                                                               &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT leaf area index after phenology",                                    &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for lit_c
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'lit_c',                                                                  &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT carbon Litter",                                                      &
! Units
    "kg  m-2 per 360days"                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for npp_dr_out
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'npp_dr_out',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT mean NPP for driving TRIFFID",                                       &
! Units
    "kg  m-2 per 360days"                                                     &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for npp
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'npp',                                                                    &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT net primary productivity",                                           &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for o3_exp_fac
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'o3_exp_fac',                                                             &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "Ozone exposure factor",                                                  &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_p
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_p',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT plant respiration",                                                  &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_w_dr_out
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_w_dr_out',                                                          &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT mean wood respiration for driving TRIFFID",                          &
! Units
    "kg C m-2 per 360days"                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_w
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_w',                                                                 &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT wood maintenance respiration",                                       &
! Units
    "kg C m-2 s-1"                                                            &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for resp_s
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'resp_s',                                                                 &
! Variable type
    VAR_TYPE_SCPOOL,                                                          &
! Long name is dynamic - set in get_var_attrs
    "",                                                                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for cosz
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'cosz',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Cosine of the zenith angle",                                             &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for diff_frac
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'diff_frac',                                                              &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox fraction of radiation that is diffuse",                          &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ecan_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ecan_gb',                                                                &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean evaporation from canopy/surface store",                     &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ei_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ei_gb',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox sublimation from lying snow or sea-ice",                         &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for esoil_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'esoil_gb',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface evapotranspiration from soil moisture store",            &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fqw_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fqw_gb',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox moisture flux from surface",                                     &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ftl_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ftl_gb',                                                                 &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface sensible heat flux",                                     &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_albedo_1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_albedo_1',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox albedo for waveband 1",                                          &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_albedo_2
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_albedo_2',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox albedo for waveband 2",                                          &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_albedo_3
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_albedo_3',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox albedo for waveband 3",                                          &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for land_albedo_4
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'land_albedo_4',                                                          &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox albedo for waveband 4",                                          &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for latent_heat
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'latent_heat',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface latent heat flux",                                       &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for q1p5m_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'q1p5m_gb',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox specific humidity at 1.5m height",                               &
! Units
    "kg kg-1"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for qw1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'qw1',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox specific humidity (total water content)",                        &
! Units
    "kg kg-1"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rainfall
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rainfall',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox rainfall rate",                                                  &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snomlt_surf_htf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snomlt_surf_htf',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox heat flux used for surface melting of snow",                     &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snowfall
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snowfall',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snowfall rate",                                                  &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_mass_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_mass_gb',                                                           &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox snowmass (on land this is lying_snow=snow_tile+snow_grnd)",      &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for surf_ht_flux_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'surf_ht_flux_gb',                                                        &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox net downward heat flux at surface over land and sea-ice fraction", &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for t1p5m_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    't1p5m_gb',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox temperature at 1.5m height",                                     &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for taux1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'taux1',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox westerly component of surface wind stress",                      &
! Units
    "N m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tauy1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tauy1',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox southerly component of surface wind stress",                     &
! Units
    "N m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tl1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tl1',                                                                    &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox ice/liquid water temperature",                                   &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tstar_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tstar_gb',                                                               &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox surface temperature",                                            &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for u1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'u1',                                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox westerly wind component",                                        &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for u10m
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'u10m',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox westerly wind component at 10 m height",                         &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for v1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'v1',                                                                     &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox southerly wind component",                                       &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for v10m
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'v10m',                                                                   &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox southerly wind component at 10 m height",                        &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ext
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ext',                                                                    &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox extraction of water from each soil layer",                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for smcl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'smcl',                                                                   &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox moisture content of each soil layer",                            &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for soil_wet
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'soil_wet',                                                               &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox total moisture content of each soil layer, as fraction of saturation", &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sthf
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sthf',                                                                   &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox frozen moisture content of each soil layer as a fraction of saturation", &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sthu
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sthu',                                                                   &
! Variable type
    VAR_TYPE_SOIL,                                                            &
! Long name
    "Gridbox unfrozen moisture content of each soil layer as a fraction of saturation", &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for alb_tile_1
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'alb_tile_1',                                                             &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile land albedo, waveband 1",                                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for alb_tile_2
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'alb_tile_2',                                                             &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile land albedo, waveband 2",                                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for alb_tile_3
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'alb_tile_3',                                                             &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile land albedo, waveband 3",                                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for alb_tile_4
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'alb_tile_4',                                                             &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile land albedo, waveband 4",                                           &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for anthrop_heat
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'anthrop_heat',                                                           &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Anthropogenic heat flux for each tile",                                  &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for catch
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'catch',                                                                  &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile surface/canopy water capacity of snow-free land tiles",             &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ecan
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ecan',                                                                   &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile evaporation from canopy/surface store for snow-free land tiles",    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ei
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ei',                                                                     &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile sublimation from lying snow for land tiles",                        &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for emis
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'emis',                                                                   &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile emissivity",                                                        &
! Units
    "1"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for esoil
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'esoil',                                                                  &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile surface evapotranspiration from soil moisture store for snow-free land tiles", &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for fqw
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'fqw',                                                                    &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile surface moisture flux for land tiles",                              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ftl
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ftl',                                                                    &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile surface sensible heat flux for land tiles",                         &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for gc
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'gc',                                                                     &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile surface conductance to evaporation for land tiles",                 &
! Units
    "m s-1"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for le
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'le',                                                                     &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile surface latent heat flux for land tiles",                           &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for q1p5m
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'q1p5m',                                                                  &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile specific humidity at 1.5m over land tiles",                         &
! Units
    "kg kg-1"                                                                 &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for rad_net_tile
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'rad_net_tile',                                                           &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile surface net radiation",                                             &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_can_melt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_can_melt',                                                          &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile melt of canopy snow",                                               &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_can
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_can',                                                               &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile canopy snow",                                                       &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_grnd_rho
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_grnd_rho',                                                          &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile density of snow on ground",                                         &
! Units
    "kg m-3"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_ground
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_ground',                                                            &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile snow on ground (snow_tile or snow_grnd)",                           &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_ice_tile
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_ice_tile',                                                          &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile total ice content in snow on ground",                               &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_liq_tile
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_liq_tile',                                                          &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile total liquid content in snow on ground",                            &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_mass
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_mass',                                                              &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile lying snow (snow_tile+snow_grnd)",                                  &
! Units
    "kg m-2"                                                                  &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_melt
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_melt',                                                              &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile snow melt (melt_tile)",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for surf_ht_flux
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'surf_ht_flux',                                                           &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Downward heat flux for each tile",                                       &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for surf_ht_store
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'surf_ht_store',                                                          &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "C*(dT/dt) for each tile",                                                &
! Units
    "W m-2"                                                                   &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for t1p5m
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    't1p5m',                                                                  &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile temperature at 1.5m over land tiles",                               &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tstar
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tstar',                                                                  &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile surface temperature",                                               &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for z0
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'z0',                                                                     &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile surface roughness",                                                 &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for z0h
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'z0h',                                                                    &
! Variable type
    VAR_TYPE_TILE,                                                            &
! Long name
    "Tile surface thermal roughness",                                         &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for tile_index
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'tile_index',                                                             &
! Variable type
    VAR_TYPE_TYPE,                                                            &
! Long name
    "Index (gridbox number) of land points with each surface type",           &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for ozone
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'ozone',                                                                  &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Variable is not available for output, so give dummy values for long name and units
    "", ""                                                                    &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for isoprene_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'isoprene_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean isoprene emission flux",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for isoprene
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'isoprene',                                                               &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT isoprene emission flux",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for terpene_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'terpene_gb',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean (mono-)terpene emission flux",                              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for terpene
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'terpene',                                                                &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT (mono-)terpene emission flux",                                       &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for methanol_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'methanol_gb',                                                            &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean methanol emission flux",                                    &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for methanol
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'methanol',                                                               &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT methanol emission flux",                                             &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for acetone_gb
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'acetone_gb',                                                             &
! Variable type
    VAR_TYPE_SURFACE,                                                         &
! Long name
    "Gridbox mean acetone emission flux",                                     &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for acetone
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'acetone',                                                                &
! Variable type
    VAR_TYPE_PFT,                                                             &
! Long name
    "PFT acetone emission flux",                                              &
! Units
    "kg m-2 s-1"                                                              &
  ),                                                                          &
!jhan: revise units here
!-----------------------------------------------------------------------------
! Metadata for soil_temp_CABLE 
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'soil_temp_CABLE',                                                              &
! Variable type
    VAR_TYPE_CABLE_SOIL,                                                      &
! Long name
    "PFT soil temperature in layers",                                         &
! Units
    "K"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for sthf_CABLE 
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'sthf_CABLE',                                                              &
! Variable type
    VAR_TYPE_CABLE_SOIL,                                                      &
! Long name
    "PFT frozen soil in layers",                                         &
! Units
    "fraction"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for smcl_CABLE 
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'smcl_CABLE',                                                              &
! Variable type
    VAR_TYPE_CABLE_SOIL,                                                      &
! Long name
    "PFT soil moisture in layers",                                         &
! Units
    ""                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_depth_CABLE 
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_depth_CABLE',                                                             &
! Variable type
    VAR_TYPE_CABLE_SNOW,                                                      &
! Long name
    "PFT snow depth in layers",                                               &
! Units
    ""                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_mass_CABLE 
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_mass_CABLE',                                                         &
! Variable type
    VAR_TYPE_CABLE_SNOW,                                                      &
! Long name
    "PFT snow mass in layers",                                      &
! Units
    ""                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_rho_CABLE 
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_rho_CABLE',                                                         &
! Variable type
    VAR_TYPE_CABLE_SNOW,                                                      &
! Long name
    "PFT snow density in layers",                                             &
! Units
    ""                                                                        &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_temp_CABLE 
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_temp_CABLE',                                                        &
! Variable type
    VAR_TYPE_CABLE_SNOW,                                                      &
! Long name
    "PFT snow temperature in layers",                                         &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_flg3l_CABLE 
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_flg3l_CABLE',                                                       &
! Variable type
    VAR_TYPE_CABLE1l_SNOW,                                                      &
! Long name
    "3 layer snow flag for tile",                                             &
! Units
    "m"                                                                       &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_age_CABLE 
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_age_CABLE',                                                         &
! Variable type
    VAR_TYPE_CABLE1l_SNOW,                                                      &
! Long name
    "snow age for tile",                                                      &
! Units
    ""                                                                        &
  ),                                                                          &
!-----------------------------------------------------------------------------
! Metadata for snow_rho1l_CABLE 
!-----------------------------------------------------------------------------
  var_metadata(                                                               &
! String identifier
    'snow_rho1l_CABLE',                                                       &
! Variable type
    VAR_TYPE_CABLE1l_SNOW,                                                      &
! Long name
    "PFT snow density in aggregate 1 layer",                                  &
! Units
    ""                                                                        &
  )                                                                           &
/
#endif
