#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/core/drivers/ascii/file_ascii_seek.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE file_ascii_seek(file, record)

  USE io_constants, ONLY : MODE_READ

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Seeks the file to a particular record (i.e. the next time
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
  TYPE(file_ascii), INTENT(INOUT) :: file  ! The file to seek
  INTEGER, INTENT(IN) :: record            ! The record number to seek to


! Work variables
  INTEGER :: error  ! Error indicator

  INTEGER :: i  ! Loop counter

!-----------------------------------------------------------------------------

! We can't seek the file if it's in define mode
  IF ( file%define_mode )                                                     &
    CALL log_fatal("file_ascii_seek",                                         &
                   "Cannot advance - file is still in define mode")

! Arbitrary seeks are only allowed in read mode
  IF ( file%mode /= MODE_READ )                                               &
    CALL log_fatal("file_ascii_seek",                                         &
                   "Arbitrary seeks are only allowed in read mode")

! We can only seek if there is a record dimension
  IF ( .NOT. file%has_record_dim )                                            &
    CALL log_fatal("file_ascii_seek",                                         &
                   "Cannot advance file - no record dimension defined")

!-----------------------------------------------------------------------------
! Actually seek the file
! We do this by rewinding the file to the start, then using fill_buffer to
! skip comment lines
!-----------------------------------------------------------------------------
  REWIND(file%unit, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("file_ascii_seek",                                         &
                   "Error rewinding file " //                                 &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  DO i = 1,record
    CALL file_ascii_fill_buffer(file)
  END DO

! By the time we get to here, the buffer is filled with the data from record
! #record, which is what we want...!

  RETURN

END SUBROUTINE file_ascii_seek
#endif
