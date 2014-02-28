#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/model_interface/map_to_land.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE map_to_land(data, data_land)

  USE ancil_info, ONLY : land_pts, land_index

  USE theta_field_sizes, ONLY : t_i_length, t_j_length

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes data on the full model grid and maps it to land points only
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  REAL, INTENT(IN) :: data(:,:,:)  ! The data to map to land points
                                   ! Dimensions are x, y, vertical levels
                                   ! Vertical levels are preserved in the mapping
  REAL, INTENT(OUT) :: data_land(land_pts,SIZE(data,3))
                                   ! The data mapped onto land points

! Work variables
  INTEGER :: i, j, l  ! Indexing variables


!-----------------------------------------------------------------------------


! We assume the data is on the full model grid
! These bounds can't be specified in the type signature because we want the
! 3rd dimension to be able to assume any size
  IF ( SIZE(data, 1) /= t_i_length .OR. SIZE(data, 2) /= t_j_length )         &
    CALL log_fatal("map_to_land",                                             &
                   "Input data can have any number of vertical levels but "// &
                   "must be on the full model grid")

  DO l = 1,land_pts
    j = (land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j - 1) * t_i_length

    data_land(l,:) = data(i,j,:)
  END DO

  RETURN

END SUBROUTINE map_to_land
#endif
