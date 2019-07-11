#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/model_interface/tiles_to_gbm.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

FUNCTION tiles_to_gbm(tile_data, tile_mask) RESULT(gbm_data)

  USE ancil_info, ONLY : land_pts, ntiles, tile_index, frac, tile_pts

  USE switches, ONLY : l_aggregate

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes data on land points * ntiles and calculates a gridbox mean value
!   for each land point. If mask is given, only the tiles for which mask
!   is true are included
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  REAL, INTENT(IN) :: tile_data(land_pts,ntiles)  ! Per tile data
  LOGICAL, INTENT(IN), OPTIONAL :: tile_mask(ntiles)
                                    ! T - include tile in calculation of gbm
                                    ! F - do not include tile in calculation

! Return type
  REAL :: gbm_data(land_pts)

! Work variables
  INTEGER :: i, p, t  ! Index variables

  LOGICAL :: tile_mask_local(ntiles)  ! Local version of tile_mask that is
                                      ! always present

!-------------------------------------------------------------------------------
! Set mask.
  IF (PRESENT( tile_mask )) THEN
    tile_mask_local(:) = tile_mask(:)
  ELSE
    tile_mask_local(:) = .TRUE.
  ENDIF

! Initialise the average.
  gbm_data(:) = 0.0

  IF ( l_aggregate ) THEN
! If l_aggregate is .TRUE., then all tile variables are essentially gridbox
! means already
    gbm_data(:) = tile_data(:,1)
  ELSE
! Otherwise, we can just use frac, since ntiles=ntype
    DO t = 1,ntiles
      IF ( tile_mask_local(t) ) THEN
        DO i = 1,tile_pts(t)
          p = tile_index(i,t)
          gbm_data(p) = gbm_data(p) + frac(p,t) * tile_data(p,t)
        ENDDO
      ENDIF
    ENDDO
  END IF

  RETURN

END FUNCTION tiles_to_gbm
#endif
