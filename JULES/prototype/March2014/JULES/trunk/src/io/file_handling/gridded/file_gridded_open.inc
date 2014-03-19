#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


FUNCTION file_gridded_open(name, mode, comm, info) RESULT(file)

  USE file_mod, ONLY : file_open

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Opens a gridded file and returns a file_gridded object representing it
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
  TYPE(file_gridded) :: file


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the MPI variables are either specified together or not at all
!-----------------------------------------------------------------------------
  IF ( PRESENT(comm) .NEQV. PRESENT(info) )                                   &
    CALL log_fatal("file_gridded_open",                                       &
                   "Only one of comm and info is present - either give a " // &
                   "value for both MPI variables for parallel access or " //  &
                   "omit both for serial access")

! All we need to do is open the underlying file
  IF ( PRESENT(comm) ) THEN
    file%fh = file_open(name, mode, comm, info)
  ELSE
    file%fh = file_open(name, mode)
  END IF

  RETURN

END FUNCTION file_gridded_open
#endif
