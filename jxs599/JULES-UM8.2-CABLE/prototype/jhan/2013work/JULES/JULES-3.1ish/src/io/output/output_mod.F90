#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/output/output_mod.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

MODULE output_mod

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN,               &
                           FORMAT_LEN, FORMAT_ASCII, FORMAT_NCDF

  USE string_utils_mod, ONLY : to_string

  USE datetime_mod, ONLY : datetime,                                          &
! Also import the comparison operators on the datetime type
                           OPERATOR(.eq.), OPERATOR(.ne.), OPERATOR(.lt.),    &
                           OPERATOR(.gt.), OPERATOR(.le.), OPERATOR(.ge.)

  USE file_ts_mod, ONLY : file_ts

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
  INTEGER, PARAMETER :: NPROFILES_MAX = 50

! Constants for types of output available
  CHARACTER(len=1), PARAMETER ::                                              &
    OUTPUT_SNAPSHOT = 'S',                                                    &
        ! Indicates output should be a snapshot every output period
    OUTPUT_ACCUM    = 'A',                                                    &
        ! Indicates output should be an accumulation over the output period
    OUTPUT_MEAN     = 'M'
        ! Indicates output should be a mean over the output period

! Dimension names for the output grid
! USED ONLY FOR 1D GRID
  CHARACTER(len=MAX_SDF_NAME_LEN), PARAMETER :: GRID_DIM_NAME = "points"
      ! The dimension name of the single grid dimension

! USED ONLY FOR 2D GRID
  CHARACTER(len=MAX_SDF_NAME_LEN), PARAMETER ::                               &
    X_DIM_NAME = "x",                                                         &
    Y_DIM_NAME = "y"

  CHARACTER(len=MAX_SDF_NAME_LEN) :: TIME_DIM_NAME = "tstep"
                         ! The dimension name for the time dimension

#if defined(NCDF_DUMMY)
  CHARACTER(len=FORMAT_LEN), PARAMETER :: OUTPUT_FORMAT = FORMAT_ASCII
#else
  CHARACTER(len=FORMAT_LEN), PARAMETER :: OUTPUT_FORMAT = FORMAT_NCDF
#endif
                                           ! The file format to use for output


!-----------------------------------------------------------------------------
! Module types
!-----------------------------------------------------------------------------
  TYPE output_profile

    CHARACTER(len=MAX_SDF_NAME_LEN) :: profile_name  ! The name of the profile

    LOGICAL :: has_open_file = .FALSE.  ! Indicates if the profile is currently
                                        ! open
    TYPE(file_ts) :: fh  ! Handle to the underlying file

    LOGICAL :: output_spinup  ! T - this profile is outputting during spinup
                              ! F - this profile is not outputting during spinup

    LOGICAL :: output_main_run  ! T - this profile is outputting (part of) the
                                !     main run
                                ! F - this profile is not outputting any of
                                !     the main run
    TYPE(datetime) :: output_start, output_end
                              ! The start and end times for output during the
                              ! main run

    INTEGER :: file_period    ! The period to use for files for this profile
    INTEGER :: output_period  ! The output period for this profile
    INTEGER :: sample_period  ! The sampling period for this profile

    TYPE(datetime) :: next_sample_time  ! The time that the data will next be
                                        ! sampled
    TYPE(datetime) :: current_output_time  ! The time that the current output
                                           ! period started
    TYPE(datetime) :: next_output_time  ! The time that the next output period
                                        ! starts

    INTEGER :: samples_in_period  ! The number of times that data has been
                                  ! sampled in this output period

    INTEGER :: nfields = 0  ! The number of fields in the profile
    TYPE(output_field), POINTER :: fields(:) => NULL()  ! The fields in the profile

  END TYPE output_profile


  TYPE output_field

    INTEGER :: var_id  ! The variable id of the model variable that
                       ! this field outputs - see model_interface_mod

    CHARACTER(len=MAX_SDF_NAME_LEN) :: output_name
                       ! The name of the variable in output files

    CHARACTER(len=1) :: type  ! The type of output to use

    INTEGER :: file_id  ! The id of the variable in the underlying file

    REAL, POINTER :: data(:,:,:) => NULL() ! The data collected for output

  END TYPE output_field


!-----------------------------------------------------------------------------
! Module variables
!-----------------------------------------------------------------------------
  CHARACTER(len=MAX_FILE_NAME_LEN) :: output_dir = "./output"
                                           ! The directory to use for output

  CHARACTER(len=MAX_SDF_NAME_LEN) :: run_id = ""  ! Identifier for the run

  INTEGER :: nx = 1, ny = 1
      ! The number of points in the x and y dimensions of the output grid
      ! respectively
      ! For a 1d grid, ny will be 1 and nx is the number of points

  INTEGER :: nprofiles = 0  ! The number of output profiles requested
  TYPE(output_profile), SAVE :: profiles(NPROFILES_MAX)  ! The requested profiles


CONTAINS


! Fortran INCLUDE statements would be preferred, but (at least) the pgf90
! compiler objects to their use when the included file contains pre-processor
! directives. At present, such directives are used to exclude files from
! the UM build, so are required. This may change in the future.
#include "register_output_profile.inc"
#include "sample_data.inc"
#include "output_data.inc"
#include "close_all.inc"

! Internal procedures
#include "internal_open_output_file.inc"
#include "internal_define_vars.inc"

END MODULE output_mod
#endif
