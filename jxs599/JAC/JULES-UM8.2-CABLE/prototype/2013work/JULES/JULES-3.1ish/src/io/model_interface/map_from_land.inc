#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/model_interface/map_from_land.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE map_from_land(data_land, data_grid)

  USE ancil_info, ONLY : land_pts, land_index

  USE theta_field_sizes, ONLY : t_i_length, t_j_length

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes data on land points only and maps it onto the full model grid
!   Only sets values at the land points - other values are left untouched
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  REAL, INTENT(IN) :: data_land(:,:)  ! The land point data to map to the
                                      ! full model grid
                                      ! Dimensions are land_pts, levels
                                      ! Vertical levels are preserved in the mapping
  REAL, INTENT(OUT) :: data_grid(t_i_length,t_j_length,SIZE(data_land,2))
                                      ! The land data mapped onto the full
                                      ! model grid

! Work variables
  INTEGER :: i, j, l  ! Indexing variables


!-----------------------------------------------------------------------------


! We assume the data is on land points
! This bound can't be specified in the type signature because we want the
! 2nd dimension to be able to assume any size
  IF ( SIZE(data_land, 1) /= land_pts )                                       &
    CALL log_fatal("map_from_land",                                           &
                   "Input data can have any number of vertical levels but "// &
                   "must be on land points")

  DO l = 1,land_pts
    j = (land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j - 1) * t_i_length

    data_grid(i,j,:) = data_land(l,:)
  END DO

  RETURN

END SUBROUTINE map_from_land
#endif
