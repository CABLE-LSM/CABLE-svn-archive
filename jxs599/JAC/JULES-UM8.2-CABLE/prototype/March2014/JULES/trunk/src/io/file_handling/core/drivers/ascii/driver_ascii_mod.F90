#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE driver_ascii_mod

  USE io_constants, ONLY : MAX_DIM_FILE, MAX_VAR_FILE, MAX_ATTR_FILE,         &
                           MAX_FILE_NAME_LEN, MAX_SDF_NAME_LEN,               &
                           MAX_ATTR_VAL_LEN, MAX_DIM_VAR

  USE string_utils_mod, ONLY : to_string

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
  CHARACTER(len=15) :: EXTENSIONS(3)       ! The file name extensions recognised
  DATA EXTENSIONS / 'asc', 'txt', 'dat' /  ! for ASCII files

  INTEGER, PARAMETER :: MAX_UNIT_NUMBER = 80

  CHARACTER(len=1) :: COMMENT_CHARS(2)
  DATA COMMENT_CHARS / '!', '#' /


!-----------------------------------------------------------------------------
! Type definitions
!-----------------------------------------------------------------------------
  TYPE file_ascii
!-----------------------------------------------------------------------------
! This type contains information required to process an ASCII file
!-----------------------------------------------------------------------------
    INTEGER :: unit  ! The unit number that the file is open on

    CHARACTER(len=MAX_FILE_NAME_LEN) :: name  ! The name of the file
                                              ! This exists purely for
                                              ! debugging purposes

    INTEGER :: mode ! The mode the file is open in
                    ! This can be one of MODE_READ or MODE_WRITE


    LOGICAL :: define_mode = .TRUE.  ! Determines if file is in define mode

    INTEGER :: ndims = 0  ! The number of dimensions in the file
    CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_names(MAX_DIM_FILE)
                                        ! The names of the dimensions
    INTEGER :: dim_sizes(MAX_DIM_FILE)  ! The sizes of the dimensions

    LOGICAL :: has_record_dim = .FALSE.  ! Indicates if the file has a
                                         ! record dimension
    CHARACTER(len=MAX_SDF_NAME_LEN) :: record_dim_name
                                         ! The name of the record dimension

    INTEGER :: nvars = 0  ! The number of variables in the file
    CHARACTER(len=MAX_SDF_NAME_LEN) :: var_names(MAX_VAR_FILE)
                          ! The names of the variables
    INTEGER :: var_sizes(MAX_VAR_FILE)  ! The sizes of the variables
    INTEGER :: var_offsets(MAX_VAR_FILE)  ! The offset of each variable
                                          ! in the buffer
    INTEGER :: var_ndims(MAX_VAR_FILE)  ! The number of dimensions that
                                        ! each variable has
    INTEGER :: var_dims(MAX_VAR_FILE, MAX_DIM_VAR)
                                        ! The dimension ids for each variable

    INTEGER :: nattrs = 0
    INTEGER :: attr_var_ids(MAX_ATTR_FILE) ! The id of the variable that each
                                           ! attribute is associated with
    CHARACTER(len=MAX_SDF_NAME_LEN+MAX_ATTR_VAL_LEN+3) :: attr_values(MAX_ATTR_FILE)
                                           ! The representation of each attribute
                                           ! that will be printed to file
                                           ! This is of the form
                                           !   name = value

    INTEGER :: record_len  ! The record length to use when reading lines from
                           ! this file
    REAL, POINTER :: buffer(:)  ! Stores data for all variables for the
                                ! current record (i.e. line)

    LOGICAL :: buffer_is_dirty = .FALSE.  ! Only used in write mode
                                          ! Indicates if the buffer has been written
                                          ! to

  END TYPE file_ascii


!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
  CHARACTER(len=20) :: out_format_str  ! The format string to use to write
                                       ! output. This is calculated once all
                                       ! variables have been defined


!-----------------------------------------------------------------------------
! Interface definitions for the overloads of read/write var for variables
! of different ranks
!-----------------------------------------------------------------------------
  INTERFACE file_ascii_read_var
    MODULE PROCEDURE file_ascii_read_var_scalar, file_ascii_read_var_1d,      &
                     file_ascii_read_var_2d,     file_ascii_read_var_3d,      &
                     file_ascii_read_var_4d,     file_ascii_read_var_5d,      &
                     file_ascii_read_var_6d,     file_ascii_read_var_7d
  END INTERFACE file_ascii_read_var

  INTERFACE file_ascii_write_var
    MODULE PROCEDURE file_ascii_write_var_scalar, file_ascii_write_var_1d,    &
                     file_ascii_write_var_2d,     file_ascii_write_var_3d,    &
                     file_ascii_write_var_4d,     file_ascii_write_var_5d,    &
                     file_ascii_write_var_6d,     file_ascii_write_var_7d
  END INTERFACE file_ascii_write_var


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
  PRIVATE
  PUBLIC EXTENSIONS,                                                          &
! Types
         file_ascii,                                                          &
! Routines for opening and closing
         file_ascii_open, file_ascii_close,                                   &
! Routines for definitions
         file_ascii_def_dim, file_ascii_def_record_dim, file_ascii_def_var,   &
         file_ascii_def_attr_real, file_ascii_def_attr_int,                   &
         file_ascii_def_attr_char, file_ascii_enddef,                         &
! Routines for seeking
         file_ascii_seek, file_ascii_advance,                                 &
! Routines for reading and writing
         file_ascii_read_var, file_ascii_write_var,                           &
! Routines for inquiry
         file_ascii_introspect, file_ascii_inquire_dim, file_ascii_inquire_var


CONTAINS

! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "file_ascii_open.inc"
#include "file_ascii_def_dim.inc"
#include "file_ascii_def_record_dim.inc"
#include "file_ascii_def_var.inc"
#include "file_ascii_def_attr.inc"
#include "file_ascii_enddef.inc"
#include "file_ascii_seek.inc"
#include "file_ascii_advance.inc"
#include "file_ascii_read_var.inc"
#include "file_ascii_write_var.inc"
#include "file_ascii_close.inc"

! Inquiry routines
#include "file_ascii_introspect.inc"
#include "file_ascii_inquire_dim.inc"
#include "file_ascii_inquire_var.inc"

! Other utility routines
#include "file_ascii_fill_buffer.inc"
#include "file_ascii_flush_buffer.inc"

END MODULE driver_ascii_mod
#endif
