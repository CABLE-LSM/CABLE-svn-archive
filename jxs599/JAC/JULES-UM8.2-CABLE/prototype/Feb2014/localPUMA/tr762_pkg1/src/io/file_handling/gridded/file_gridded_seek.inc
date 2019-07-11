#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_gridded_seek(file, record)

  USE file_mod, ONLY : file_seek

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Seeks the file to before a particular record (i.e. the next time a
!   values are read from the file using file_read_var, they will be read from
!   the requested record)
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_gridded), INTENT(INOUT) :: file  ! The file to seek
  INTEGER, INTENT(IN) :: record              ! The record number to seek to


!-----------------------------------------------------------------------------

! We just need to defer to the underlying file
  CALL file_seek(file%fh, record)

  RETURN

END SUBROUTINE file_gridded_seek
#endif
