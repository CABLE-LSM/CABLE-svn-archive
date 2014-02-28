#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/input/time_varying/interpolation/interpolate.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

FUNCTION interpolate(data, interp_flag, tsteps_in_period, tstep)              &
                                                     RESULT(interpolated_data)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Returns data interpolated to the current time
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  REAL, POINTER, INTENT(IN) :: data(:,:,:,:)
                      ! The data to interpolate
                      ! Dimensions are x, y, levels, time
                      ! Marked as a pointer to preserve array indexing
  CHARACTER(len=*), INTENT(IN) :: interp_flag
                      ! The interpolation scheme to use
  INTEGER, INTENT(IN) :: tsteps_in_period(-1:1)
                      ! The number of timesteps in the interpolated data
                      ! for each data timestep
                      ! I.e. tsteps_in_period(-1) is the number of interpolated
                      !      (i.e. model) timesteps between data(-1) and data(0)
                      !      tsteps_in_period(0) is the number of interpolated
                      !      timesteps between data(0) and data(1)
                      ! etc...
  INTEGER, INTENT(IN) :: tstep
                      ! The interpolated (i.e. model) timestep that we want data for
                      ! This is a number from 0 to (tsteps_in_period(0) - 1)


! Return type
  REAL :: interpolated_data(SIZE(data,1), SIZE(data,2), SIZE(data,3))
                      ! The interpolated data for the current time

! Work variables
  REAL :: weights(-1:2)  ! The weights for each time

! Used in CASE ( INTERP_AVE_BACKWARD ) and CASE ( INTERP_AVE_FORWARD )
  REAL :: n1, n2  ! The number of timesteps in the first and second intervals

! Used in CASE ( INTERP_AVE_CENTRED )
  REAL :: n  ! The number of timesteps in the interval

! Used in all INTERP_AVE_* cases
  REAL :: weighted_data(SIZE(data,1), SIZE(data,2), SIZE(data,3))
  REAL :: numer(SIZE(data,1), SIZE(data,2), SIZE(data,3))
  REAL :: denom(SIZE(data,1), SIZE(data,2), SIZE(data,3))

  REAL :: t  ! Real version of tstep


!-----------------------------------------------------------------------------

! Reject values that we can't interpolate for - it only makes sense to give
! values for times between current and next data
  IF ( tstep < 0 .OR. tsteps_in_period(0) <= tstep )                          &
    CALL log_fatal("interpolate",                                             &
                   "tstep to interpolate to is invalid for given data " //    &
                   "periods")


  SELECT CASE ( interp_flag )
    CASE ( INTERP_AVE_BACKWARD )
!-----------------------------------------------------------------------------
! Inputs are backward time averages, i.e. time average ending at given time
!
! This is a slightly modified version of the scheme in JULES v3.0. It has
! been modified to attempt to provide sensible values for varying length
! data timesteps (e.g. monthly data) - the previous scheme was limited to
! fixed size data timesteps. This extension should be considered experimental.
!
! In the case of fixed length timesteps, the new scheme reduces to exactly
! the scheme in JULES v3.0, and so can be safely used in this case.
!-----------------------------------------------------------------------------
      n1 = REAL(tsteps_in_period(0))
      n2 = REAL(tsteps_in_period(1))
      t  = REAL(tstep)

      weights(1) = 1.0 - ( ABS(t + (n2/n1)*t - n1 + 1.0) / (n1 + n2) )
      weights(0) = MAX(1.0 - ( (t + (n2/n1)*t + n1 + 1.0) / (n1 + n2) ), 0.0)
      weights(2) = MAX(1.0 - ((n1 + 2.0*n2 - t - (n2/n1)*t - 1.0) / (n1 + n2)), 0.0)

      weighted_data(:,:,:) = ( weights(0) * data(:,:,:,0) +                   &
                               weights(1) * data(:,:,:,1) +                   &
                               weights(2) * data(:,:,:,2) )
      denom(:,:,:) = 0.5 * (data(:,:,:,0) + data(:,:,:,2)) + 3.0 * data(:,:,:,1)
      numer(:,:,:) = 4.0 * data(:,:,:,1)

      WHERE ( denom > EPSILON(denom) )
        interpolated_data = weighted_data * numer / denom
      ELSEWHERE
        interpolated_data = 0.0
      END WHERE


    CASE ( INTERP_AVE_CENTRED )
!-----------------------------------------------------------------------------
! Inputs are centred time averages, i.e. time average centred on given time
!-----------------------------------------------------------------------------
      n = REAL(tsteps_in_period(0))
      t = REAL(tstep)

      weights(0)  = 1.0 - ( (2.0*t + 1.0) / (2.0*n) )
      weights(1)  = 1.0 - weights(0)
      weighted_data(:,:,:) = ( weights(0) * data(:,:,:,0)  +                  &
                               weights(1) * data(:,:,:,1) )

      IF ( t > n/2 ) THEN
        denom(:,:,:) = 0.5 * (data(:,:,:,0) + data(:,:,:,2))                  &
                     + 3.0 * data(:,:,:,1)
        numer(:,:,:) = 4.0 * data(:,:,:,1)
      ELSE
        denom(:,:,:) = 0.5 * (data(:,:,:,-1) + data(:,:,:,1))                 &
                     + 3.0 * data(:,:,:,0)
        numer(:,:,:) = 4.0 * data(:,:,:,0)
      ENDIF

      WHERE ( denom > EPSILON(denom) )
        interpolated_data = weighted_data * numer / denom
      ELSEWHERE
        interpolated_data = 0.0
      END WHERE


    CASE ( INTERP_AVE_FORWARD )
!-----------------------------------------------------------------------------
! Inputs are forward time averages, i.e. time average starting at given time
!
! This is a slightly modified version of the scheme in JULES v3.0. It has
! been modified to attempt to provide sensible values for varying length
! data timesteps (e.g. monthly data) - the previous scheme was limited to
! fixed size data timesteps. This extension should be considered experimental.
!
! In the case of fixed length timesteps, the new scheme reduces to exactly
! the scheme in JULES v3.0, and so can be safely used in this case.
!-----------------------------------------------------------------------------
      n1 = REAL(tsteps_in_period(-1))
      n2 = REAL(tsteps_in_period(0))
      t  = REAL(tstep)

      weights(0)  = 1.0 - ( ABS(t + (n1/n2)*t - n1 + 1) / (n1 + n2) )
      weights(-1) = MAX(1.0 - ( (t + (n1/n2)*t + n1 + 1) / (n1 + n2) ), 0.0)
      weights(1)  = MAX(1.0 - ((n1 + 2*n2 - t - (n1/n2)*t - 1) / (n1 + n2)), 0.0)

      weighted_data(:,:,:) = ( weights(-1) * data(:,:,:,-1) +                 &
                               weights(0)  * data(:,:,:,0)  +                 &
                               weights(1)  * data(:,:,:,1) )
      denom(:,:,:) = 0.5 * (data(:,:,:,-1) + data(:,:,:,1)) + 3.0 * data(:,:,:,0)
      numer(:,:,:) = 4.0 * data(:,:,:,0)

      WHERE ( denom > EPSILON(denom) )
        interpolated_data = weighted_data * numer / denom
      ELSEWHERE
        interpolated_data = 0.0
      END WHERE


    CASE ( INTERP_INSTANT )
!-----------------------------------------------------------------------------
! Inputs are instantaneous values at given times - interpolation is linear
! between times
!-----------------------------------------------------------------------------
! Calculate the weights based on how far through the time period we are
      weights(0) = REAL(tsteps_in_period(0) - tstep) / REAL(tsteps_in_period(0))
      weights(1) = 1.0 - weights(0)

! Calculate the data to return
      interpolated_data(:,:,:) = weights(0) * data(:,:,:,0) +                 &
                                 weights(1) * data(:,:,:,1)


    CASE ( NO_INTERP_END )
!-----------------------------------------------------------------------------
! No interpolation, (current) value is valid over time interval ending at
! given time
!-----------------------------------------------------------------------------
      interpolated_data(:,:,:) = data(:,:,:,1)


    CASE ( NO_INTERP_CENTRED )
!-----------------------------------------------------------------------------
! No interpolation, (current) value is valid over time interval centred on
! given time
!-----------------------------------------------------------------------------
      IF ( tstep < (tsteps_in_period(0) / 2) ) THEN
! If we are closer to current data than next data, use the current data
        interpolated_data(:,:,:) = data(:,:,:,0)
      ELSE
! Otherwise use the next data
        interpolated_data(:,:,:) = data(:,:,:,1)
      END IF


    CASE ( NO_INTERP_START )
!-----------------------------------------------------------------------------
! No interpolation, (current) value is valid over time interval starting at
! given time
!-----------------------------------------------------------------------------
      interpolated_data(:,:,:) = data(:,:,:,0)


    CASE DEFAULT
      CALL log_fatal("interpolate",                                           &
                     "Unrecognised interpolation flag - " // interp_flag)
  END SELECT

  RETURN

END FUNCTION interpolate
#endif
