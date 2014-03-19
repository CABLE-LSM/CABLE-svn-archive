#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/core/drivers/ncdf/file_ncdf_open.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

FUNCTION file_ncdf_open(name, mode) RESULT(file)

  USE io_constants, ONLY : MODE_READ, MODE_WRITE

  USE netcdf, ONLY :                                                          &
! Constants
    NF90_NOWRITE, NF90_CLOBBER,                                               &
! Procedures
    nf90_open, nf90_create

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a NetCDF file and returns a file_ncdf object representing it
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  CHARACTER(len=*), INTENT(IN) :: name
                              ! The name of the file
  INTEGER, INTENT(IN) :: mode ! The mode to open the file
                              ! One of MODE_READ or MODE_WRITE

! Return type
  TYPE(file_ncdf) :: file


! Work variables
  INTEGER :: ncid ! The NetCDF id of the opened file
  INTEGER :: error ! Error code for any errors that occur

!-----------------------------------------------------------------------------


  SELECT CASE ( mode )
    CASE ( MODE_READ )
      CALL log_info("file_ncdf_open",                                         &
                    "Opening file " // TRIM(name) // " for reading")
! Open file for reading only - file must already exist
      error = nf90_open(name, NF90_NOWRITE, ncid)

    CASE ( MODE_WRITE )
      CALL log_info("file_ncdf_open",                                         &
                    "Opening file " // TRIM(name) // " for writing")
! Create an empty file for (reading and) writing - if a file with the
! same name exists, overwrite it
      error = nf90_create(name, NF90_CLOBBER, ncid)

    CASE DEFAULT
! Read and write are the only supported modes
      CALL log_fatal("file_ncdf_open",                                        &
                     "Unsupported mode - " // TRIM(to_string(mode)))

  END SELECT

  IF ( error /= NF90_NOERR )                                                  &
    CALL log_fatal_ncdf("file_ncdf_open",                                     &
                        "Error opening file " // TRIM(name), error)


! Initialise the file_ncdf object
  file%name = name
  file%mode = mode
  file%id   = ncid

  RETURN

END FUNCTION file_ncdf_open
#endif
