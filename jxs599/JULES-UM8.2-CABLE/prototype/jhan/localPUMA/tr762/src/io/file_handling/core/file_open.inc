#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_open(name, mode, comm, info) RESULT(file)

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
  INTEGER, OPTIONAL, INTENT(IN) :: comm  ! MPI communicator to use for parallel
                                         ! I/O
                                         ! If not given, serial I/O is used
  INTEGER, OPTIONAL, INTENT(IN) :: info  ! MPI info object to use for parallel
                                         ! I/O
                                         ! If not given, serial I/O is used

! Return type
  TYPE(file_handle) :: file


! Work variables
  INTEGER :: last_dot  ! The index of the last . (dot) in the file name
                       ! Used to extract the file extension
  CHARACTER(len=15) :: extension  ! The extension for the file name


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the MPI variables are either specified together or not at all
!-----------------------------------------------------------------------------
  IF ( PRESENT(comm) .NEQV. PRESENT(info) )                                   &
    CALL log_fatal("file_open",                                               &
                   "Only one of comm and info is present - either give a " // &
                   "value for both MPI variables for parallel access or " //  &
                   "omit both for serial access")

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
    IF ( PRESENT(comm) ) THEN
      file%ascii = file_ascii_open(name, mode, comm, info)
    ELSE
      file%ascii = file_ascii_open(name, mode)
    END IF

  ELSE IF ( ANY(EXTENSIONS_NCDF == extension) ) THEN
! Assign the correct driver
    file%driver = DRIVER_NCDF
! Initialise the NetCDF representation of the file
    IF ( PRESENT(comm) ) THEN
      file%ncdf = file_ncdf_open(name, mode, comm, info)
    ELSE
      file%ncdf = file_ncdf_open(name, mode)
    END IF

  ELSE
! File type not recognised
    CALL log_fatal("file_open",                                               &
                   "Unrecognised file extension for file " // TRIM(name) //   &
                   " - see docs for supported file types")
  END IF

  RETURN

END FUNCTION file_open
#endif
