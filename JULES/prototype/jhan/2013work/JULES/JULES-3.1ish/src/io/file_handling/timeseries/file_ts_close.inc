#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/timeseries/file_ts_close.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE file_ts_close(file)

  USE file_gridded_mod, ONLY : file_gridded_close

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
  TYPE(file_ts), INTENT(INOUT) :: file ! The file to close


!-----------------------------------------------------------------------------

! Make sure the currently open file is closed
  CALL file_gridded_close(file%open_file)

! Free the memory associated with the files and file_times arrays if it was
! allocated (i.e. if not using time templating)
  IF ( .NOT. file%use_time_template ) THEN
    DEALLOCATE(file%files)
    DEALLOCATE(file%file_times)
  END IF

  RETURN

END SUBROUTINE file_ts_close
#endif
