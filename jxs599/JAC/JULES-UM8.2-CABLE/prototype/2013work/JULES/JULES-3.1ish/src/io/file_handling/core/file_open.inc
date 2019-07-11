#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/core/file_open.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

FUNCTION file_open(name, mode) RESULT(file)

  USE driver_ascii_mod, ONLY : EXTENSIONS_ASCII => EXTENSIONS, file_ascii_open
  USE driver_ncdf_mod, ONLY : EXTENSIONS_NCDF => EXTENSIONS, file_ncdf_open

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a file and returns a file_handle object representing it
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
  TYPE(file_handle) :: file


! Work variables
  INTEGER :: last_dot  ! The index of the last . (dot) in the file name
                       ! Used to extract the file extension
  CHARACTER(len=15) :: extension  ! The extension for the file name


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Get the extension for the file name (i.e. everything after the last dot)
!-----------------------------------------------------------------------------
  last_dot = INDEX(name, '.', BACK=.TRUE.)
  extension = name(last_dot+1:)

! Select a driver based on the file extension
  IF ( ANY(EXTENSIONS_ASCII == extension) ) THEN
! Assign the correct driver
    file%driver = DRIVER_ASCII
! Initialise the ASCII representation of the file
    file%ascii = file_ascii_open(name, mode)

  ELSE IF ( ANY(EXTENSIONS_NCDF == extension) ) THEN
! Assign the correct driver
    file%driver = DRIVER_NCDF
! Initialise the NetCDF representation of the file
    file%ncdf = file_ncdf_open(name, mode)

  ELSE
! File type not recognised
    CALL log_fatal("file_open",                                               &
                   "Unrecognised file type - see docs for supported types")
  END IF

  RETURN

END FUNCTION file_open
#endif
