#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE spinup_init()

  USE ancil_info, ONLY : land_pts, sm_levels

  USE prognostics, ONLY : smcl, t_soil

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Allocates variables required for spinup and stores initial state
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: i  ! Loop counter


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Allocate space to store the state of the spinup variables for comparison
!-----------------------------------------------------------------------------
  DO i = 1,nvars
    SELECT CASE ( spinup_vars(i)%identifier )
      CASE ( 'smcl' )
        ALLOCATE(spinup_vars(i)%data(land_pts,sm_levels))
        spinup_vars(i)%data(:,:) = smcl(:,:)

      CASE ( 't_soil' )
        ALLOCATE(spinup_vars(i)%data(land_pts,sm_levels))
        spinup_vars(i)%data(:,:) = t_soil(:,:)

      CASE DEFAULT
        CALL log_fatal("spinup_init",                                         &
                       "Unrecognised variable for spinup - " //               &
                       TRIM(spinup_vars(i)%identifier))
    END SELECT
  END DO

  RETURN

END SUBROUTINE spinup_init
#endif
