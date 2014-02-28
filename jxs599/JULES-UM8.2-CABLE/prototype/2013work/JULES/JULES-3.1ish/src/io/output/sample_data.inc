#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/output/sample_data.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE sample_data()

  USE datetime_mod, ONLY : datetime_advance

  USE model_time_mod, ONLY : current_time, is_spinup, spinup_cycle, spinup_start

  USE model_interface_mod, ONLY : extract_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   For each profile, decide if we need to sample data at this timestep, and
!   sample it if we do
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: data_shape(3)  ! The shape to give the data array
  REAL, ALLOCATABLE :: data(:,:,:)  ! Array to temporarily store the data for
                                    ! the current timestep

  INTEGER :: i,j  ! Loop counter

!-----------------------------------------------------------------------------

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

!-----------------------------------------------------------------------------
! Check if we need to open a new file
! This is the case if:
!   - The profile is outputting during spinup and we are on the first
!     timestep of a spinup cycle
!   - The profile is outputting part of the main run and we are in the main
!     run at the time when output starts
! We do this check in two stages for readability
!-----------------------------------------------------------------------------
    IF ( is_spinup .AND. profiles(i)%output_spinup .AND.                      &
         current_time == spinup_start )                                       &
      CALL internal_open_output_file(profiles(i))

    IF ( .NOT. is_spinup .AND. profiles(i)%output_main_run .AND.              &
         current_time == profiles(i)%output_start )                           &
      CALL internal_open_output_file(profiles(i))

!-----------------------------------------------------------------------------
! Sample data if required
!-----------------------------------------------------------------------------
    IF ( profiles(i)%next_sample_time <= current_time ) THEN
! For each field in the profile, sample the data
      DO j = 1,profiles(i)%nfields
        data_shape(:) = SHAPE(profiles(i)%fields(j)%data)
        ALLOCATE(data(data_shape(1),data_shape(2),data_shape(3)))

! First, extract the data on the model grid from the model variable
        CALL extract_var(profiles(i)%fields(j)%var_id, data)

        SELECT CASE ( profiles(i)%fields(j)%type )
          CASE ( OUTPUT_SNAPSHOT )
! In the case of a snapshot output, just update the data to the current value
! of the variable
            profiles(i)%fields(j)%data(:,:,:) = data(:,:,:)

          CASE ( OUTPUT_ACCUM, OUTPUT_MEAN )
! In the case of accumulation or mean, update the accumulation
            profiles(i)%fields(j)%data(:,:,:) =                               &
                               profiles(i)%fields(j)%data(:,:,:) + data(:,:,:)

          CASE DEFAULT
            CALL log_fatal("sample_data",                                     &
                           "Unrecognised output type - " //                   &
                            profiles(i)%fields(j)%type)
        END SELECT

        DEALLOCATE(data)
      END DO

      profiles(i)%samples_in_period = profiles(i)%samples_in_period + 1

      profiles(i)%next_sample_time = datetime_advance(                        &
        profiles(i)%next_sample_time, profiles(i)%sample_period               &
      )
    END IF

! Nothing more to do!!
  END DO

  RETURN

END SUBROUTINE sample_data
#endif
