#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


MODULE driver_ncdf_mod

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, MAX_VAR_FILE

  USE string_utils_mod, ONLY : to_string

! Import this globally as it is universally used
  USE netcdf, ONLY : NF90_NOERR

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE


!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
  CHARACTER(len=15) :: EXTENSIONS(2)  ! The file name extensions recognised
  DATA EXTENSIONS / 'nc', 'cdf' /     ! for NetCDF files

!-----------------------------------------------------------------------------
! Type definitions
!-----------------------------------------------------------------------------
  TYPE var_ncdf

    INTEGER :: id  ! The NetCDF id of the variable
    INTEGER :: ndims  ! The number of NON-RECORD dimensions that the
                       ! variable has
    LOGICAL :: is_record  ! Does the variable use the record dimension?

  END TYPE var_ncdf


  TYPE file_ncdf
!-----------------------------------------------------------------------------
! This type contains information required to process a NetCDF file
!-----------------------------------------------------------------------------
    INTEGER :: id  ! NetCDF id of file

    CHARACTER(len=MAX_FILE_NAME_LEN) :: name  ! The name of the file
                                              ! This exists purely for
                                              ! debugging purposes

    INTEGER :: mode ! The mode the file is open in
                    ! This can be one of MODE_READ or MODE_WRITE

    INTEGER :: record_dim = -1  ! The id of the record dimension
                                ! -1 means no record dimension

    INTEGER :: current_record = 1  ! By default, we read from the 1st record

    INTEGER :: nvars = 0  ! The number of variables in the file
    TYPE(var_ncdf) :: vars(MAX_VAR_FILE)  ! The variables in the file

  END TYPE file_ncdf


!-----------------------------------------------------------------------------
! Interface definitions for the overloads of read/write var for variables
! of different ranks
!-----------------------------------------------------------------------------
  INTERFACE file_ncdf_read_var
    MODULE PROCEDURE file_ncdf_read_var_scalar, file_ncdf_read_var_1d,        &
                     file_ncdf_read_var_2d,     file_ncdf_read_var_3d,        &
                     file_ncdf_read_var_4d,     file_ncdf_read_var_5d,        &
                     file_ncdf_read_var_6d,     file_ncdf_read_var_7d
  END INTERFACE file_ncdf_read_var

  INTERFACE file_ncdf_write_var
    MODULE PROCEDURE file_ncdf_write_var_scalar, file_ncdf_write_var_1d,      &
                     file_ncdf_write_var_2d,     file_ncdf_write_var_3d,      &
                     file_ncdf_write_var_4d,     file_ncdf_write_var_5d,      &
                     file_ncdf_write_var_6d,     file_ncdf_write_var_7d
  END INTERFACE file_ncdf_write_var


!-----------------------------------------------------------------------------
! Visibility declarations
!-----------------------------------------------------------------------------
  PRIVATE
  PUBLIC EXTENSIONS,                                                          &
! Types
         file_ncdf,                                                           &
! Routines for opening and closing
         file_ncdf_open, file_ncdf_close,                                     &
! Routines for definitions
         file_ncdf_def_dim, file_ncdf_def_record_dim, file_ncdf_def_var,      &
         file_ncdf_def_attr_real, file_ncdf_def_attr_int,                     &
         file_ncdf_def_attr_char, file_ncdf_enddef,                           &
! Routines for seeking
         file_ncdf_seek, file_ncdf_advance,                                   &
! Routines for reading and writing
         file_ncdf_read_var, file_ncdf_write_var,                             &
! Routines for inquiry
         file_ncdf_introspect, file_ncdf_inquire_dim, file_ncdf_inquire_var


CONTAINS

! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "file_ncdf_open.inc"
#include "file_ncdf_def_dim.inc"
#include "file_ncdf_def_record_dim.inc"
#include "file_ncdf_def_var.inc"
#include "file_ncdf_def_attr.inc"
#include "file_ncdf_enddef.inc"
#include "file_ncdf_seek.inc"
#include "file_ncdf_advance.inc"
#include "file_ncdf_read_var.inc"
#include "file_ncdf_write_var.inc"
#include "file_ncdf_close.inc"

! Inquiry routines
#include "file_ncdf_introspect.inc"
#include "file_ncdf_inquire_dim.inc"
#include "file_ncdf_inquire_var.inc"

#include "log_fatal_ncdf.inc"

END MODULE driver_ncdf_mod
#endif
