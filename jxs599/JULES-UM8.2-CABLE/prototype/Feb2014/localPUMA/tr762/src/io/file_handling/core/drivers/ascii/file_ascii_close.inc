#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_close(file)

  USE io_constants, ONLY : MODE_WRITE

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Closes and frees any resources consumed by the given file object
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ascii), INTENT(INOUT) :: file ! The file to close


! Work variables
  INTEGER :: error ! Error code for any errors that occur

!-----------------------------------------------------------------------------

  CALL log_info("file_ascii_close",                                           &
                "Closing file " // TRIM(file%name))

! Flush the buffer to disk only if it has uncommitted writes
  IF ( file%mode == MODE_WRITE .AND. file%buffer_is_dirty )                   &
    CALL file_ascii_flush_buffer(file)

! Deallocate the memory associated with the buffer
  DEALLOCATE(file%buffer)
  NULLIFY(file%buffer)

  CLOSE(file%unit, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("file_ascii_close",                                        &
                   "Error closing file " // TRIM(file%name) //                &
                   " (IOSTAT=" // TRIM(to_string(error)) // ")")

  RETURN

END SUBROUTINE file_ascii_close
#endif
