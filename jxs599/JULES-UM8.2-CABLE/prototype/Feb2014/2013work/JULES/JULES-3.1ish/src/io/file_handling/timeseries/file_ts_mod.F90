#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE file_ts_mod

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, MAX_SDF_NAME_LEN,               &
                           MAX_ATTR_VAL_LEN, MAX_DIM_FILE, MAX_VAR_FILE,      &
                           MAX_ATTR_FILE

  USE datetime_mod, ONLY : datetime,                                          &
! Also import the comparison operators on the datetime type
                           OPERATOR(.eq.), OPERATOR(.ne.), OPERATOR(.lt.),    &
                           OPERATOR(.gt.), OPERATOR(.le.), OPERATOR(.ge.)

  USE string_utils_mod, ONLY : to_string

  USE file_gridded_mod, ONLY : file_gridded

  USE dictionary_mod, ONLY : dict

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------------
  CHARACTER(len=MAX_SDF_NAME_LEN) :: TIME_INDEX_VAR_NAME = 'time'

!-----------------------------------------------------------------------------
! Type definitions
!-----------------------------------------------------------------------------
! Internal type containing information about the grid
  TYPE grid_ts

    LOGICAL :: is_1d  ! Indicates if the grid is 1d or 2d

    CHARACTER(len=MAX_SDF_NAME_LEN) :: grid_name  ! ONLY USED IF is_1d=T
                                                  ! The name of the single
                                                  ! grid dimension
    INTEGER :: grid_npoints = 0  ! ONLY USED IF is_1d=T
                                 ! The size of the single grid dimension

    CHARACTER(len=MAX_SDF_NAME_LEN) :: grid_x_name, grid_y_name
                                 ! ONLY USED IF is_1d=F
                                 ! The names of the x and y dimensions of
                                 ! the grid respectively
    INTEGER :: grid_nx, grid_ny  ! ONLY USED IF is_1d=F
                                 ! The sizes of the x and y dimensions,
                                 ! respectively

  END TYPE grid_ts


! Internal type containing information about a dimension
  TYPE dim_ts

    INTEGER :: current_id  ! The id of the dimension in the currently open file
    CHARACTER(len=MAX_SDF_NAME_LEN) :: name ! The name of the dimension
    INTEGER :: length  ! The length of the dimension

  END TYPE dim_ts


! Internal type containing information about a variable
  TYPE var_ts

    INTEGER :: current_id  ! The id of the variable in the currently open file

    CHARACTER(len=MAX_SDF_NAME_LEN) :: name  ! The name of the variable in the
                                             ! file(s)

    LOGICAL :: has_levels  ! T - variable has a vertical levels dimension
                           !     in the file(s)
                           ! F - variable does not have a vertical levels
                           !     dimension

    INTEGER :: levels_dim  ! The id of the vertical levels dimension (as an
                           ! index in the dims array on the file_ts object)

    LOGICAL :: use_time   ! Indicates whether the variable uses the time dimension

! Dictionaries to hold attribute values
    TYPE(dict) :: attrs_real  ! The values for real valued attributes
                              ! Indexed by attribute name
    TYPE(dict) :: attrs_int  ! The values for integer valued attributes
                             ! Indexed by attribute name
    TYPE(dict) :: attrs_char  ! The values for char valued attributes
                              ! Indexed by attribute name

  END TYPE var_ts


  TYPE file_ts
!-----------------------------------------------------------------------------
! This type and associated functions and subroutines are intended to provide
! access to timeseries data in a unified way, whether that data is contained
! in a single file or multiple files
!-----------------------------------------------------------------------------
    INTEGER :: mode ! The mode in which to open files
                    ! This can be one of MODE_READ or MODE_WRITE


!-----------------------------------------------------------------------------
! Properties defining the characteristics of the data
!-----------------------------------------------------------------------------
    TYPE(datetime) :: data_start ! The date and time of the first data
    TYPE(datetime) :: data_end   ! The date and time of the last data
    INTEGER :: data_period ! The period of the data in the files
                           ! (in seconds or a 'special' period)
                           ! NOTE: that this means all data in a file must
                           ! have the same period
    LOGICAL :: is_climatology = .FALSE. ! Indicates whether the data represents
                                        ! a climatology.
                                        ! In this case, the same data are
                                        ! reused for every year


!-----------------------------------------------------------------------------
! Properties containing information about the dimensions, variables and
! attributes contained in each file - this will need to be redefined for
! each file that is opened
!-----------------------------------------------------------------------------
    LOGICAL :: define_mode = .TRUE.  ! Indicates whether object is in
                                     ! 'define mode', i.e. if we can define
                                     ! dimensions and variables

    TYPE(grid_ts) :: grid  ! The grid that variables in the file(s) are on

    LOGICAL :: has_time_dim = .FALSE.  ! Indicates whether the time dimension
                                       ! has been defined
    TYPE(dim_ts) :: time_dim ! The time dimension for the file(s)

    INTEGER :: ndims = 0  ! The number of non-grid, non-time dimensions
                          ! that have been defined
    TYPE(dim_ts) :: dims(MAX_DIM_FILE) ! The non-grid, non-time dimensions
                                       ! in the file(s)

    INTEGER :: nvars = 0  ! The number of variables that have been defined
    TYPE(var_ts) :: vars(MAX_VAR_FILE) ! The variables in the file(s)

! Dictionaries to hold global attribute values
    TYPE(dict) :: attrs_real  ! The values for real valued attributes
                                  ! Indexed by attribute name
    TYPE(dict) :: attrs_int  ! The values for integer valued attributes
                                 ! Indexed by attribute name
    TYPE(dict) :: attrs_char  ! The values for char valued attributes
                                  ! Indexed by attribute name

    INTEGER :: time_index_var_id  ! The variable id for the time index
                                  ! in the currently open file - note that
                                  ! since the time index does not have a grid,
                                  ! we have to use the file_handle underlying
                                  ! the file_gridded object to manipulate
                                  ! it, hence this is a variable id in the
                                  ! file_handle object, i.e. file_ts%open_file%fh


!-----------------------------------------------------------------------------
! Properties used to determine what file to open
!-----------------------------------------------------------------------------
    LOGICAL :: use_time_template

! With time templating
    CHARACTER(len=MAX_FILE_NAME_LEN) :: template
    INTEGER :: tpl_period

! With a file list
    INTEGER :: nfiles ! The number of files
    CHARACTER(len=MAX_FILE_NAME_LEN), POINTER :: files(:) => NULL()
                       ! List of file names
    TYPE(datetime), POINTER :: file_times(:) => NULL()
                       ! Time of first data for each file


!-----------------------------------------------------------------------------
! Properties related to the currently open file
!-----------------------------------------------------------------------------
    LOGICAL :: has_open_file = .FALSE.
                              ! T - the file in open_file represents an
                              !     actual open file
                              ! F - the file object in open_file is uninitialised

    INTEGER :: open_file_index ! The index in the files/file_times arrays
                               ! of the currently open file

    TYPE(file_gridded) :: open_file ! The currently open file

    TYPE(datetime) :: next_file_start ! The start time of the next file

    TYPE(datetime) :: current_datetime ! The datetime that this file_ts
                                       ! object is currently pointing to
                                       ! NOTE: it is possible for this datetime
                                       ! to contain a different value to
                                       ! the current model datetime - it is
                                       ! maintained to ensure the correct file
                                       ! is kept open

  END TYPE file_ts


! Overloads for file_ts_def_attr
  INTERFACE file_ts_def_attr
    MODULE PROCEDURE file_ts_def_attr_real, file_ts_def_attr_int,             &
                     file_ts_def_attr_char
  END INTERFACE file_ts_def_attr

CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "file_ts_open.inc"
#include "file_ts_def_grid.inc"
#include "file_ts_def_dim.inc"
#include "file_ts_def_time_dim.inc"
#include "file_ts_def_var.inc"
#include "file_ts_def_attr.inc"
#include "file_ts_enddef.inc"
#include "file_ts_seek_to_datetime.inc"
#include "file_ts_advance.inc"
#include "file_ts_read_var.inc"
#include "file_ts_write_var.inc"
#include "file_ts_close.inc"

! Utility routines
#include "file_ts_internal_open_file.inc"

END MODULE file_ts_mod
#endif
