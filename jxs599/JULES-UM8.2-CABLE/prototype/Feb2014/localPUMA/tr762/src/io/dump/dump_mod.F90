#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE dump_mod

  USE io_constants, ONLY : FORMAT_LEN, FORMAT_ASCII, FORMAT_NCDF,             &
                           MAX_SDF_NAME_LEN

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
  INTEGER, PARAMETER :: MAX_DIM_DUMP = 11
  INTEGER, PARAMETER :: MAX_VAR_DUMP = 28

  CHARACTER(len=MAX_SDF_NAME_LEN), PARAMETER ::                               &
    LAND_DIM_NAME    = "land",                                                &
    PFT_DIM_NAME     = "pft",                                                 &
    SC_POOL_DIM_NAME = "scpool",                                              &
    SNOW_DIM_NAME    = "snow",                                                &
    SOIL_DIM_NAME    = "soil",                                                &
    TILE_DIM_NAME    = "tile",                                                &
    TYPE_DIM_NAME    = "type",                                                &
    SCALAR_DIM_NAME  = "scalar",                                              &
    NOLEVS_DIM_NAME  = "olevs",                                               &
    NFARRAY_DIM_NAME = "nfarray",                                             &
    SEED_DIM_NAME    = "seed"

#if defined(NCDF_DUMMY)
  CHARACTER(len=FORMAT_LEN), PARAMETER :: DUMP_FORMAT = FORMAT_ASCII
#else
  CHARACTER(len=FORMAT_LEN), PARAMETER :: DUMP_FORMAT = FORMAT_NCDF
#endif

!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
  PRIVATE
  PUBLIC MAX_VAR_DUMP, required_vars_for_configuration, read_dump, write_dump


CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "required_vars_for_configuration.inc"
#include "read_dump.inc"
#include "write_dump.inc"

END MODULE dump_mod
#endif
