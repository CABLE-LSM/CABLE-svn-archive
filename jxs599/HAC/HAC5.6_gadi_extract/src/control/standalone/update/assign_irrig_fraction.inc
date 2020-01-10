#if !defined(UM_JULES)
! *****************************COPYRIGHT***************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use
! and distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT***************************************

SUBROUTINE assign_irrig_fraction ()

USE logging_mod, ONLY: log_fatal

USE string_utils_mod, ONLY: to_string

USE ancil_info, ONLY:                                                         &
land_pts, nsurft, frac_surft, nsoilt

USE jules_surface_types_mod, ONLY: ntype

USE crop_vars_mod, ONLY:                                                      &
!  imported logicals with intent(in)
   frac_irr_all_tiles,                                                        &
!  imported arrays with intent(in)
   frac_irr_all,irrtiles,                                                     &
   frac_irr_soilt, frac_irr_surft,                                            &
   frac_irr_old_soilt,                                                        &
   nirrtile

  !USE jules_vegetation_mod, ONLY : irr_crop

IMPLICIT NONE

INTEGER ::                                                                    &
 m,l,k,n
                       ! Loop indices

INTEGER ::   alltiles(nsurft)

! Local arrays
REAL ::                                                                       &
 frac_irr_rest_soilt(land_pts,nsoilt)
                      ! Remaining irrigation fraction to assign to tiles

!-------------------------------------------------------------------------------
!   Update irrigation ancillary
!-------------------------------------------------------------------------------

! Loop over land points and assign irrigation fraction
! note frac_irr_all has shape (land_points, year) or (land_points, 1)
DO m = 1,nsoilt
  DO l = 1,land_pts
    frac_irr_old_soilt(l,m) = frac_irr_soilt(l,m)
    frac_irr_soilt(l,m)     = frac_irr_all(l,1)

    !     Set a minimum threshold for irrigation:
    IF ( frac_irr_soilt(l,m) < 1.0e-3 ) frac_irr_soilt(l,m) = 0.0
  END DO
END DO

!-------------------------------------------------------------------------------
!   Assign irrigation fraction to each tile
!   hadrd - this was originally done in physiol.f90
!-------------------------------------------------------------------------------
IF ( .NOT. frac_irr_all_tiles ) THEN
  ! assign irrigated fraction to pre-defined tiles ONLY
  ! in decreasing order of importance
  frac_irr_surft(:,:) = 0.0
  frac_irr_rest_soilt(:,:)    = frac_irr_soilt(:,:)

  ! create index of all tiles (veg+non-veg)
  DO m = 1,nsurft
    alltiles(m) = m
  END DO

  DO l = 1,land_pts
    DO k = 1,nirrtile ! loop over irrigated pfts

      !==============================================================================
      !**NOTICE REGARDING SOIL TILING**
      !
      !The following section facilitates the use of soil tiling. As implemented,
      !there are two soil tiling options:
      !
      !nsoilt == 1
      !Operate as with a single soil tile, functionally identical to JULES upto
      ! at least vn4.7 (Oct 2016)
      ! This means that a soilt variable being passed 'up' to the surface is
      ! broadcast to the surft variable (with weighting by frac if requred)
      !
      !nsoilt > 1
      !Operate with nsoilt = nsurft, with a direct mapping between them
      ! This means that a soilt variable being passed 'up' to the surface is simply
      ! copied into the surft variable
      !
      ! This will need to be refactored for other tiling approaches. This note
      ! will be replicated elsewhere in the code as required
      !
      !These comments apply until **END NOTICE REGARDING SOIL TILING**
      !==============================================================================

              !Set the current soil tile (see notice above)
      IF (nsoilt == 1) THEN
        !There is only 1 soil tile
        m = 1
      ELSE ! nsoilt == nsurft
        !Soil tiles map directly on to surface tiles
        m = k
      END IF !nsoilt

      !==============================================================================
      !**END NOTICE REGARDING SOIL TILING**
      !==============================================================================

      ! for each tile, check if the index corresponds with index in irrtiles
      DO n = 1,ntype
        IF ( alltiles(n) == irrtiles(k) ) THEN
          ! assign (remaining) irrigated fraction to this tile
          frac_irr_surft(l,n)      = MIN(frac_irr_rest_soilt(l,m),            &
                                         frac_surft(l,n) )
          frac_irr_rest_soilt(l,m) = frac_irr_rest_soilt(l,m)                 &
                                     - frac_irr_surft(l,n)
          ! check for negative remaining frac_irr_soilt
          IF (frac_irr_rest_soilt(l,m) < 0.0) THEN
            CALL log_fatal("update_ancil_irrig",                              &
              "Error in assigning irrigated fraction at point " //            &
              to_string(l) // " - " //                                        &
              "irrigated fraction: " //                                       &
              to_string(frac_irr_soilt(l,m)) // ", " //                       &
              "tile fraction: " // to_string(frac_surft(l,n)) // ", " //      &
              "irrigated fraction assigned:"// to_string(frac_irr_surft(l,n)))
          END IF
        END IF
      END DO ! ntype
    END DO ! irrtilenames
  END DO ! land_pts

  !to prevent negative fractions (this can happen due to numerical inaccuracy)
  frac_irr_surft(:,:) = MAX(frac_irr_surft(:,:),0.0)

  ! in sf_evap, frac_irr_surft is used as a multiplier to frac
  WHERE ( frac_surft > 0.0 ) frac_irr_surft = frac_irr_surft / frac_surft

ELSE ! if frac_irr_all_tiles = TRUE
  ! assign irrigation fraction to all tiles
  ! to reproduce the original results
  frac_irr_surft(:,:) = 0.0
  DO n = 1,nsurft

    !==============================================================================
    !**NOTICE REGARDING SOIL TILING**
    !
    !The following section facilitates the use of soil tiling. As implemented,
    !there are two soil tiling options:
    !
    !nsoilt == 1
    !Operate as with a single soil tile, functionally identical to JULES upto
    ! at least vn4.7 (Oct 2016)
    ! This means that a soilt variable being passed 'up' to the surface is
    ! broadcast to the surft variable (with weighting by frac if requred)
    !
    !nsoilt > 1
    !Operate with nsoilt = nsurft, with a direct mapping between them
    ! This means that a soilt variable being passed 'up' to the surface is simply
    ! copied into the surft variable
    !
    ! This will need to be refactored for other tiling approaches. This note
    ! will be replicated elsewhere in the code as required
    !
    !These comments apply until **END NOTICE REGARDING SOIL TILING**
    !==============================================================================

          !Set the current soil tile (see notice above)
    IF (nsoilt == 1) THEN
      !There is only 1 soil tile
      m = 1
    ELSE ! nsoilt == nsurft
      !Soil tiles map directly on to surface tiles
      m = n
    END IF !nsoilt

    !==============================================================================
    !**END NOTICE REGARDING SOIL TILING**
    !==============================================================================

    WHERE ( frac_irr_soilt(:,m) > EPSILON(1.0) )
      frac_irr_surft(:,n) = frac_irr_soilt(:,m)
    END WHERE
  END DO

END IF ! .not. frac_irr_all_tiles

END SUBROUTINE assign_irrig_fraction
#endif
