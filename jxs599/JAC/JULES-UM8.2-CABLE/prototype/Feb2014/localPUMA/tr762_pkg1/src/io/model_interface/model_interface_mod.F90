#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE model_interface_mod

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_ATTR_VAL_LEN

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

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
! Length for identifiers of model variables
  INTEGER, PARAMETER :: IDENTIFIER_LEN = 20

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
    pft_dim_size     = 5,                                                     &
    nvg_dim_size     = 4,                                                     &
    type_dim_size    = 9,                                                     &
    tile_dim_size    = 9,                                                     &
    snow_dim_size    = 0,                                                     &
    soil_dim_size    = 4,                                                     &
    scpool_dim_size  = 1


!-----------------------------------------------------------------------------
! Information about the actual variables available for output
!-----------------------------------------------------------------------------
! Constants for the different 'types' of variable
  INTEGER, PARAMETER ::                                                       &
    VAR_TYPE_SURFACE = 0,                                                     &
    VAR_TYPE_PFT     = 1,                                                     &
    VAR_TYPE_NVG     = 2,                                                     &
    VAR_TYPE_TYPE    = 3,                                                     &
    VAR_TYPE_TILE    = 4,                                                     &
    VAR_TYPE_SNOW    = 5,                                                     &
    VAR_TYPE_SOIL    = 6,                                                     &
    VAR_TYPE_SCPOOL  = 7,                                                     &         
    VAR_TYPE_CABLE_SOIL    = 8,                                               &
    VAR_TYPE_CABLE_SNOW    = 9,                                               &  
    VAR_TYPE_CABLE1l_SNOW    = 10

! Derived type to contain metadata about model variables
  TYPE var_metadata

    CHARACTER(len=IDENTIFIER_LEN) :: identifier
                              ! The string identifier of the variable

    INTEGER :: var_type  ! The type of the variable - must be one of the above

    CHARACTER(len=MAX_ATTR_VAL_LEN) :: long_name
                              ! The value of the long_name attribute
    CHARACTER(len=MAX_ATTR_VAL_LEN) :: units
                              ! The value of the units attribute

  END TYPE var_metadata

! Array holding the metadata for all model variables that we can use for input
! or output
  !INTEGER, PARAMETER :: N_VARS = 213
  INTEGER, PARAMETER :: N_VARS = 223
  TYPE(var_metadata) :: metadata(N_VARS)

! Include the metadata DATA statement
#include "variable_metadata.inc"


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
  PRIVATE
  PUBLIC                                                                      &
! Parameters
    IDENTIFIER_LEN,                                                           &
    PFT_DIM_NAME_OUT, NVG_DIM_NAME_OUT, TYPE_DIM_NAME_OUT, TILE_DIM_NAME_OUT, &
    SNOW_DIM_NAME_OUT, SOIL_DIM_NAME_OUT, SCPOOL_DIM_NAME_OUT,                &
! Variables
    pft_dim_name, nvg_dim_name, type_dim_name, tile_dim_name, snow_dim_name,  &
    soil_dim_name, scpool_dim_name,                                           &
    pft_dim_size, nvg_dim_size, type_dim_size, tile_dim_size, snow_dim_size,  &
    soil_dim_size, scpool_dim_size,                                           &
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
#include "tiles_to_gbm.inc"

END MODULE model_interface_mod
#endif
