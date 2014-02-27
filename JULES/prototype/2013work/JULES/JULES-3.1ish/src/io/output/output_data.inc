#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/output/output_data.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE output_data()

  USE datetime_mod, ONLY : datetime_advance

  USE model_time_mod, ONLY : current_time, is_spinup, spinup_cycle,           &
                             timestep_len

  USE file_ts_mod, ONLY : file_ts_write_var, file_ts_advance

  USE datetime_mod, ONLY : datetime_to_string

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   For each profile, decide if we need to output data at this timestep, and
!   output it if we do
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  TYPE(datetime) :: next_time  ! The time at the next timestep

  INTEGER :: data_shape(3)  ! The shape to give the data array
  REAL, ALLOCATABLE :: data(:,:,:)  ! Array to temporarily store the data for
                                    ! the current timestep

  INTEGER :: i,j  ! Loop counter

!-----------------------------------------------------------------------------


  next_time = datetime_advance(current_time, timestep_len)

  DO i = 1,nprofiles
!-----------------------------------------------------------------------------
! First check if the profile is currently active
!-----------------------------------------------------------------------------
    IF ( is_spinup .AND. .NOT. profiles(i)%output_spinup )                    &
! If we are in spinup but profile is not outputting, we don't need to do anything
      CYCLE

    IF ( .NOT. is_spinup .AND. ( .NOT. profiles(i)%output_main_run .OR.       &
         ( current_time < profiles(i)%output_start .OR.                       &
           profiles(i)%output_end <= current_time ) ) )                       &
! If we are in the main run but either not outputting any of the main run or
! the time is out of the range of times the profile is outputting, we don't
! need to do anything
      CYCLE

! We already know that we have a suitable file open due to the check in
! sample_data

!-----------------------------------------------------------------------------
! We need to output data on this timestep if it is the last timestep in an
! output period, i.e. the next timestep is in the next output period
!-----------------------------------------------------------------------------
    IF ( next_time >= profiles(i)%next_output_time ) THEN
! For each field in the profile, output the data
      DO j = 1,profiles(i)%nfields
! Get the data that we want to output
        data_shape(:) = SHAPE(profiles(i)%fields(j)%data)
        ALLOCATE(data(data_shape(1),data_shape(2),data_shape(3)))

        SELECT CASE ( profiles(i)%fields(j)%type )
          CASE ( OUTPUT_SNAPSHOT )
! In the case of a snapshot output, just output the current data
            data(:,:,:) = profiles(i)%fields(j)%data(:,:,:)

          CASE ( OUTPUT_ACCUM )
! In the case of an accumulation, we output the accumulated value multiplied
! by the number of timesteps in a sample period
            data(:,:,:) = profiles(i)%fields(j)%data(:,:,:) *                 &
                          REAL(profiles(i)%sample_period) / REAL(timestep_len)

          CASE ( OUTPUT_MEAN )
! In the case of a mean, divide the accumulated value by the number of samples
! taken
            data(:,:,:) = profiles(i)%fields(j)%data(:,:,:) /                 &
                                                 profiles(i)%samples_in_period
! Only in the case of means do we zero the field's data property and start
! gathering again
            profiles(i)%fields(j)%data(:,:,:) = 0.0

          CASE DEFAULT
            CALL log_fatal("sample_data",                                     &
                           "Unrecognised output type - " //                   &
                            profiles(i)%fields(j)%type)
        END SELECT

! Actually write the data to file
! For now, we always write a full grid - hence write_subgrid=F
        CALL file_ts_write_var(                                               &
          profiles(i)%fh, profiles(i)%fields(j)%file_id, data, .FALSE.        &
        )

        DEALLOCATE(data)
      END DO

!-----------------------------------------------------------------------------
! Update tracking variables and, if this is not the last output for the
! profile, advance the file
!-----------------------------------------------------------------------------
      profiles(i)%samples_in_period = 0

      profiles(i)%current_output_time = profiles(i)%next_output_time
      profiles(i)%next_output_time = datetime_advance(                        &
        profiles(i)%next_output_time, profiles(i)%output_period               &
      )

! Check if this is the last output time for the currently open output file
! If it is not, then advance the file
      IF ( profiles(i)%next_output_time <= profiles(i)%fh%data_end )          &
        CALL file_ts_advance(profiles(i)%fh)

    END IF

! Nothing more to do!!
  END DO

  RETURN

END SUBROUTINE output_data
#endif
