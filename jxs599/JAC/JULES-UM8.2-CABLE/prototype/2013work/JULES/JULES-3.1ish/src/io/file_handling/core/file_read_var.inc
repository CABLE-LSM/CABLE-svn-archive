#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/core/file_read_var.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE file_read_var_scalar(file, var_id, value, start)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE driver_ascii_mod, ONLY : file_ascii_read_var_scalar
  USE driver_ncdf_mod, ONLY : file_ncdf_read_var_scalar

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a scalar value from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_handle), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, INTENT(OUT) :: value  ! The value read from file

! Optional arguments
  INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                             ! (one value for each NON-RECORD
                                             ! dimension of the variable)


! Work variables
  INTEGER :: local_start(MAX_DIM_VAR)

!-----------------------------------------------------------------------------
  local_start(:) = 1
  IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

  SELECT CASE ( file%driver )
    CASE ( DRIVER_ASCII )
      CALL file_ascii_read_var_scalar(file%ascii, var_id, value, local_start)

    CASE ( DRIVER_NCDF )
      CALL file_ncdf_read_var_scalar(file%ncdf, var_id, value, local_start)

    CASE DEFAULT
      CALL log_fatal("file_read_var_scalar",                                  &
                     "Unrecognised driver - see docs for available drivers")
  END SELECT

  RETURN

END SUBROUTINE file_read_var_scalar


SUBROUTINE file_read_var_1d(file, var_id, values, start, count)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE driver_ascii_mod, ONLY : file_ascii_read_var_1d
  USE driver_ncdf_mod, ONLY : file_ncdf_read_var_1d

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 1d array from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_handle), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, DIMENSION(:), INTENT(OUT) :: values  ! The values read from file

! Optional arguments
  INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                             ! (one value for each NON-RECORD
                                             ! dimension of the variable)
  INTEGER, INTENT(IN), OPTIONAL :: count(:)  ! The number of points to read
                                             ! in each dimension of the variable


! Work variables
  INTEGER :: local_start(MAX_DIM_VAR)
  INTEGER :: local_count(MAX_DIM_VAR)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
  local_start(:) = 1
  IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default count is the size of the array to fill, with the rest of the
! dimensions being one
  local_count(:) = 1
  local_count(1) = SIZE(values)
  IF ( PRESENT(count) ) local_count(:SIZE(count)) = count(:)

  SELECT CASE ( file%driver )
    CASE ( DRIVER_ASCII )
      CALL file_ascii_read_var_1d(                                            &
        file%ascii, var_id, values, local_start, local_count                  &
      )

    CASE ( DRIVER_NCDF )
      CALL file_ncdf_read_var_1d(                                             &
        file%ncdf, var_id, values, local_start, local_count                   &
      )

    CASE DEFAULT
      CALL log_fatal("file_read_var_1d",                                      &
                     "Unrecognised driver - see docs for available drivers")
  END SELECT

  RETURN

END SUBROUTINE file_read_var_1d


SUBROUTINE file_read_var_2d(file, var_id, values, start, count)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE driver_ascii_mod, ONLY : file_ascii_read_var_2d
  USE driver_ncdf_mod, ONLY : file_ncdf_read_var_2d

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 2d array from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_handle), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, DIMENSION(:,:), INTENT(OUT) :: values  ! The values read from file

! Optional arguments
  INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                             ! (one value for each NON-RECORD
                                             ! dimension of the variable)
  INTEGER, INTENT(IN), OPTIONAL :: count(:)  ! The number of points to read
                                             ! in each dimension of the variable


! Work variables
  INTEGER :: local_start(MAX_DIM_VAR)
  INTEGER :: local_count(MAX_DIM_VAR)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
  local_start(:) = 1
  IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default count is the size of the array to fill from the first
! available dimensions
  local_count(:) = 1
  local_count(1:2) = SHAPE(values)
  IF ( PRESENT(count) ) local_count(:SIZE(count)) = count(:)

  SELECT CASE ( file%driver )
    CASE ( DRIVER_ASCII )
      CALL file_ascii_read_var_2d(                                            &
        file%ascii, var_id, values, local_start, local_count                  &
      )

    CASE ( DRIVER_NCDF )
      CALL file_ncdf_read_var_2d(                                             &
        file%ncdf, var_id, values, local_start, local_count                   &
      )

    CASE DEFAULT
      CALL log_fatal("file_read_var_2d",                                      &
                     "Unrecognised driver - see docs for available drivers")
  END SELECT

  RETURN

END SUBROUTINE file_read_var_2d


SUBROUTINE file_read_var_3d(file, var_id, values, start, count)

  USE io_constants, ONLY : MAX_DIM_VAR

  USE driver_ascii_mod, ONLY : file_ascii_read_var_3d
  USE driver_ncdf_mod, ONLY : file_ncdf_read_var_3d

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads a 3d array from a variable in the given file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_handle), INTENT(INOUT) :: file  ! The file to read from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, DIMENSION(:,:,:), INTENT(OUT) :: values  ! The values read from file

! Optional arguments
  INTEGER, INTENT(IN), OPTIONAL :: start(:)  ! The point to start reading from
                                             ! (one value for each NON-RECORD
                                             ! dimension of the variable)
  INTEGER, INTENT(IN), OPTIONAL :: count(:)  ! The number of points to read
                                             ! in each dimension of the variable


! Work variables
  INTEGER :: local_start(MAX_DIM_VAR)
  INTEGER :: local_count(MAX_DIM_VAR)

!-----------------------------------------------------------------------------

! The default is to start at the 1st element of each dimension
  local_start(:) = 1
  IF ( PRESENT(start) ) local_start(:SIZE(start)) = start(:)

! The default count is the size of the array to fill from the first
! available dimensions
  local_count(:) = 1
  local_count(1:3) = SHAPE(values)
  IF ( PRESENT(count) ) local_count(:SIZE(count)) = count(:)

  SELECT CASE ( file%driver )
    CASE ( DRIVER_ASCII )
      CALL file_ascii_read_var_3d(                                            &
        file%ascii, var_id, values, local_start, local_count                  &
      )

    CASE ( DRIVER_NCDF )
      CALL file_ncdf_read_var_3d(                                             &
        file%ncdf, var_id, values, local_start, local_count                   &
      )

    CASE DEFAULT
      CALL log_fatal("file_read_var_3d",                                      &
                     "Unrecognised driver - see docs for available drivers")
  END SELECT

  RETURN

END SUBROUTINE file_read_var_3d
#endif