#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/input/time_varying/time_varying_close_all.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE close_all()

  USE file_ts_mod, ONLY : file_ts_close

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Closes all the time varying input files and frees all resources
!   consumed by them
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: i,j  ! Loop counters


!-----------------------------------------------------------------------------

  DO i = 1,nfiles
! Deallocate all the data arrays for the input_fields
    DO j = 1,files(i)%nfields
      DEALLOCATE(files(i)%fields(j)%data)
      NULLIFY(files(i)%fields(j)%data)
    END DO

    DEALLOCATE(files(i)%fields)
    NULLIFY(files(i)%fields)

    CALL file_ts_close(files(i)%fh)
  END DO


END SUBROUTINE close_all
#endif
