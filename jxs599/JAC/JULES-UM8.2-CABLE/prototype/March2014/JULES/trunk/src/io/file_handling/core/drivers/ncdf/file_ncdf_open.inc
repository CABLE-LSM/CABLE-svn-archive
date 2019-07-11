#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ncdf_open(name, mode, comm, info) RESULT(file)

  USE io_constants, ONLY : MODE_READ, MODE_WRITE

  USE netcdf, ONLY :                                                          &
! Constants
    NF90_NOWRITE, NF90_CLOBBER, NF90_NETCDF4, NF90_MPIIO,                     &
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
  INTEGER, OPTIONAL, INTENT(IN) :: comm  ! MPI communicator to use for parallel
                                         ! I/O
                                         ! If not given, serial I/O is used
  INTEGER, OPTIONAL, INTENT(IN) :: info  ! MPI info object to use for parallel
                                         ! I/O
                                         ! If not given, serial I/O is used

!-----------------------------------------------------------------------------
! We know from the conditions imposed by file_open that comm and info are
! either both present or both not present
!-----------------------------------------------------------------------------

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
! We don't need to specify comm or info to get parallel access here - any
! number of parallel readers are possible anyway...
      error = nf90_open(name, NF90_NOWRITE, ncid)

    CASE ( MODE_WRITE )
      CALL log_info("file_ncdf_open",                                         &
                    "Opening file " // TRIM(name) // " for writing")
! Create an empty file for (reading and) writing - if a file with the
! same name exists, overwrite it
      IF ( PRESENT(comm) ) THEN
! Don't try to specify MPI variables to NetCDF routines if we are using the
! dummy MPI library
#if defined(MPI_DUMMY)
        error = nf90_create(name, NF90_CLOBBER, ncid)
#else
! If we want to use parallel I/O, we have to force NetCDF-4 files to be created
        error = nf90_create(                                                  &
          name, IOR(NF90_CLOBBER, IOR(NF90_NETCDF4, NF90_MPIIO)), ncid,       &
          COMM=comm, INFO=info                                                &
        )
#endif
      ELSE
        error = nf90_create(name, NF90_CLOBBER, ncid)
      END IF

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
