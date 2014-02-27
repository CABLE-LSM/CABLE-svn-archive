#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/core/drivers/ncdf/file_ncdf_def_dim.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

FUNCTION file_ncdf_def_dim(file, dim_name, dim_len) RESULT(dim_id)

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MODE_READ, MODE_WRITE

  USE netcdf, ONLY : nf90_inq_dimid, nf90_inquire_dimension, nf90_def_dim

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines a dimension on the given file, returning the dimension id
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
  INTEGER, INTENT(IN) :: dim_len  ! The length of the dimension

! Return type
  INTEGER :: dim_id               ! The dimension id

! Work variables
  INTEGER :: error  ! The current error code (if any)
  INTEGER :: len_in_file  ! The length of the dimension in the file (only
                          ! used in read mode to validate dimension size)


!-----------------------------------------------------------------------------

! We do different things depending on whether we are open for reading
! or writing
  SELECT CASE ( file%mode )
    CASE ( MODE_READ )
! In read mode, we just get the id of the dimension and validate it's size
      error = nf90_inq_dimid(file%id, dim_name, dim_id)
      IF ( error /= NF90_NOERR )                                              &
        CALL log_fatal_ncdf("file_ncdf_def_dim",                              &
                            "Error getting dimension id", error)

      error = nf90_inquire_dimension(file%id, dim_id, len=len_in_file)
      IF ( error /= NF90_NOERR )                                              &
        CALL log_fatal_ncdf("file_ncdf_def_dim",                              &
                            "Error getting dimension size", error)
      IF ( dim_len /= len_in_file )                                           &
        CALL log_fatal("file_ncdf_def_dim",                                   &
                       "Dimension length mismatch - length in file: " //      &
                       TRIM(to_string(len_in_file)) // ", expected: " //      &
                       TRIM(to_string(dim_len)))

    CASE ( MODE_WRITE )
! In write mode, we define the dimension
      error = nf90_def_dim(file%id, dim_name, dim_len, dim_id)
      IF ( error /= NF90_NOERR )                                              &
        CALL log_fatal_ncdf("file_ncdf_def_dim",                              &
                            "Error defining dimension", error)

! No default case as we already know that MODE_READ and MODE_WRITE are
! the only options
  END SELECT

  RETURN

END FUNCTION file_ncdf_def_dim
#endif
