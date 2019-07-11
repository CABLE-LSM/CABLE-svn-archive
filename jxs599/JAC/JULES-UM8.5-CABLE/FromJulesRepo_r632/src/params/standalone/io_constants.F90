#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE io_constants

  IMPLICIT NONE

! File formats
  INTEGER, PARAMETER :: FORMAT_LEN = 3
  CHARACTER(len=FORMAT_LEN), PARAMETER ::                                     &
    FORMAT_ASCII = 'asc',                                                     &
      ! Indicates an ASCII file
    FORMAT_NCDF = 'nc'
      ! Indicates a NetCDF file


! Modes for opening files
  INTEGER, PARAMETER ::                                                       &
    MODE_READ  = 1,                                                           &
    MODE_WRITE = 2


! 'Special' units
  INTEGER, PARAMETER ::                                                       &
    UNIT_STDIN       = 5,                                                     &
    UNIT_STDOUT      = 6,                                                     &
    NAMELIST_UNIT    = 1,                                                     &
    POINTS_FILE_UNIT = 2,                                                     &
    FILE_LIST_UNIT   = 3,                                                     &
    IMOGEN_UNIT      = 99


! Constant to specify that an attribute is global (passed in place of var_id)
  INTEGER, PARAMETER :: ATTR_GLOBAL = -1


  REAL, PARAMETER :: MDI = -1e20  ! Missing data indicator for output files


! Various maximums for quantities required in IO
  INTEGER, PARAMETER ::                                                       &
    MAX_FILE_NAME_LEN = 500,                                                  &
    MAX_SDF_NAME_LEN = 50,                                                    &
    MAX_ATTR_VAL_LEN = 50,                                                    &
    MAX_DIM_FILE = 40,                                                        &
    MAX_VAR_FILE = 50,                                                        &
    MAX_ATTR_FILE = 500,                                                      &
    MAX_DIM_VAR = 7,                                                          &
    MAX_ATTR_VAR = 10

END MODULE io_constants
#endif
