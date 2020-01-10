#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE model_interface_mod

USE io_constants, ONLY: max_sdf_name_len, max_attr_val_len

USE logging_mod, ONLY: log_info, log_debug, log_warn, log_error, log_fatal

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This module provides the interface between the IO routines and the
!   model variables
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
! Length for identifiers of model variables
INTEGER, PARAMETER :: identifier_len = 24

! The name of each levels dimension in output files
CHARACTER(LEN=max_sdf_name_len), PARAMETER ::                                 &
  bl_level_dim_name_out = 'bllevel',                                          &
  pft_dim_name_out     = 'pft',                                               &
  cpft_dim_name_out    = 'cpft',                                              &
  nvg_dim_name_out     = 'nvg',                                               &
  type_dim_name_out    = 'type',                                              &
  tile_dim_name_out    = 'tile',                                              &
  soilt_dim_name_out   = 'soilt',                                             &
  snow_dim_name_out    = 'snow',                                              &
  soil_dim_name_out    = 'soil',                                              &
  scpool_dim_name_out  = 'scpool',                                            &
  sclayer_dim_name_out = 'sclayer',                                           &
  tracer_dim_name_out  = 'tracer',                                            &
  bedrock_dim_name_out = 'bedrock'


!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
! The name of each levels dimension in input files
CHARACTER(LEN=max_sdf_name_len) ::                                            &
  bl_level_dim_name = 'bllevel',                                              &
  pft_dim_name     = 'pft',                                                   &
  cpft_dim_name    = 'cpft',                                                  &
  nvg_dim_name     = 'nvg',                                                   &
  type_dim_name    = 'type',                                                  &
  tile_dim_name    = 'tile',                                                  &
  soilt_dim_name   = 'soilt',                                                 &
  snow_dim_name    = 'snow',                                                  &
  soil_dim_name    = 'soil',                                                  &
  scpool_dim_name  = 'scpool',                                                &
  sclayer_dim_name = 'sclayer',                                               &
  tracer_dim_name  = 'tracer',                                                &
  bedrock_dim_name = 'bedrock'

! The size of each possible levels dimension
INTEGER ::                                                                    &
  bl_level_dim_size = 0,                                                      &
  pft_dim_size     = 5,                                                       &
  cpft_dim_size    = 0,                                                       &
  nvg_dim_size     = 4,                                                       &
  type_dim_size    = 9,                                                       &
  tile_dim_size    = 9,                                                       &
  soilt_dim_size   = 9,                                                       &
  snow_dim_size    = 0,                                                       &
  soil_dim_size    = 4,                                                       &
  scpool_dim_size  = 1,                                                       &
  sclayer_dim_size = 1,                                                       &
  tracer_dim_size  = 0,                                                       &
  bedrock_dim_size = 100


!-----------------------------------------------------------------------------
! Information about the actual variables available for output
!-----------------------------------------------------------------------------
! Constants for the different 'types' of variable
INTEGER, PARAMETER ::                                                         &
  var_type_bl_level       = 0,                                                &
  var_type_surface        = var_type_bl_level+1,                              &
  var_type_pft            = var_type_surface+1,                               &
  var_type_cpft           = var_type_pft+1,                                   &
  var_type_nvg            = var_type_cpft+1,                                  &
  var_type_type           = var_type_nvg+1,                                   &
  var_type_surft          = var_type_type+1,                                  &
  var_type_soilt          = var_type_surft+1,                                 &
  var_type_soilt_soil     = var_type_soilt+1,                                 &
  var_type_snow           = var_type_soilt_soil+1,                            &
  var_type_soil           = var_type_snow+1,                                  &
  var_type_scpool         = var_type_soil+1,                                  &
  var_type_bedrock        = var_type_scpool+1,                                &
  var_type_rp             = var_type_bedrock+1,                               &
  var_type_sclayer        = var_type_rp+1,                                    &
  var_type_soilt_sclayer  = var_type_sclayer+1,                               &
  var_type_tracer         = var_type_soilt_sclayer + 1,                       &
  var_type_cable_soil     = var_type_tracer + 1,                              &
  var_type_cable_snow     = var_type_cable_soil + 1,                          & 
  var_type_cable1L_snow   = var_type_cable_snow +1

! Derived type to contain metadata about model variables
TYPE var_metadata

  CHARACTER(LEN=identifier_len) :: identifier
                            ! The string identifier of the variable

  INTEGER :: var_type  ! The type of the variable - must be one of the above

  CHARACTER(LEN=max_attr_val_len) :: long_name
                            ! The value of the long_name attribute
  CHARACTER(LEN=max_attr_val_len) :: units
                            ! The value of the units attribute

END TYPE var_metadata

! Array holding the metadata for all model variables that we can use for input
! or output !CABLE_LSM adds 10 prognostics for tiled soil/snow prognostics
INTEGER, PARAMETER :: n_vars = 558

TYPE(var_metadata) :: metadata(n_vars)

! Include the metadata DATA statement
#include "variable_metadata.inc"


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
PRIVATE
PUBLIC                                                                        &
! Parameters
    identifier_len, bl_level_dim_name_out,                                    &
    pft_dim_name_out, cpft_dim_name_out, nvg_dim_name_out, type_dim_name_out, &
    tile_dim_name_out, soilt_dim_name_out, snow_dim_name_out,                 &
    soil_dim_name_out, scpool_dim_name_out, bedrock_dim_name_out,             &
    sclayer_dim_name_out, tracer_dim_name_out,                                &
! Variables
    bl_level_dim_name,                                                        &
    pft_dim_name, cpft_dim_name, nvg_dim_name, type_dim_name, tile_dim_name,  &
    soilt_dim_name, snow_dim_name, soil_dim_name, scpool_dim_name,            &
    bedrock_dim_name, sclayer_dim_name, tracer_dim_name,                      &
    bl_level_dim_size,                                                        &
    pft_dim_size, cpft_dim_size, nvg_dim_size,                                &
    type_dim_size, tile_dim_size, soilt_dim_size, snow_dim_size,              &
    soil_dim_size, scpool_dim_size, bedrock_dim_size, sclayer_dim_size,       &
    tracer_dim_size,                                                          &
! Routines for changing between string and integer identifiers
    get_var_id, get_string_identifier,                                        &
! Variable inquiry routines
    get_var_levs_dims, get_var_attrs,                                         &
! Routines for getting and setting values
    extract_var, populate_var


CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "get_var_id.inc"
#include "get_string_identifier.inc"
#include "get_var_levs_dims.inc"
#include "get_var_attrs.inc"
#include "extract_var.inc"
#include "populate_var.inc"

#include "map_to_land.inc"
#include "map_from_land.inc"

END MODULE model_interface_mod
#endif
