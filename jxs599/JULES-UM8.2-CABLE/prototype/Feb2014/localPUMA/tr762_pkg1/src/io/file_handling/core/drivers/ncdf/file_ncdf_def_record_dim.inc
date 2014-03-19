#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_ncdf_def_record_dim(file, dim_name) RESULT(dim_id)

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MODE_READ, MODE_WRITE

  USE netcdf, ONLY : nf90_inq_dimid, nf90_def_dim, NF90_UNLIMITED

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines the record dimension on the given file, returning the dimension id
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(INOUT) :: file
                                  ! The file to define the dimension on
  CHARACTER(len=*), INTENT(IN) :: dim_name
                                  ! The name of the dimension

! Return type
  INTEGER :: dim_id               ! The dimension id

! Work variables
  INTEGER :: error  ! The current error code (if any)


!-----------------------------------------------------------------------------

! We do different things depending on whether we are open for reading
! or writing
  SELECT CASE ( file%mode )
    CASE ( MODE_READ )
! In read mode, we just get the id of the named dimension
! We don't actually care whether it is a 'proper' unlimited dimension or not,
! since we won't be writing to it
      error = nf90_inq_dimid(file%id, dim_name, dim_id)
      IF ( error /= NF90_NOERR )                                              &
        CALL log_fatal_ncdf("file_ncdf_def_record_dim",                       &
                            "Error getting dimension id", error)

    CASE ( MODE_WRITE )
! In write mode, we define an unlimited dimension
      error = nf90_def_dim(file%id, dim_name, NF90_UNLIMITED, dim_id)
      IF ( error /= NF90_NOERR )                                              &
        CALL log_fatal_ncdf("file_ncdf_def_record_dim",                       &
                            "Error defining unlimited dimension", error)

! No default case as we already know that MODE_READ and MODE_WRITE are
! the only options
  END SELECT

! Store the id of the record dimension for future use
  file%record_dim = dim_id

  RETURN

END FUNCTION file_ncdf_def_record_dim
#endif
