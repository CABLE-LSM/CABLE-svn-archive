#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/initial_conditions/total_snow_init.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE total_snow_init()

  USE ancil_info, ONLY : land_pts, nsmax, ntiles, tile_index, tile_pts

  USE prognostics, ONLY : nsnow, rgrain, snow_tile, t_soil, rgrainL,          &
                          rho_snow_grnd, sice, sliq, snow_grnd, snowDepth,    &
                          tsnow

  USE snow_param, ONLY : rho_snow_const, rho_snow_fresh, canSnowTile, ds

  USE switches, ONLY : l_spec_albedo

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Sets up the initial snow conditions from snow_tile only
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: i, j, k, n  ! Index variables

  REAL :: snow_on_ground(land_pts,ntiles)
                ! Snow considered to be on the ground (not in canopy) (kg m-2)
                ! This is all the snow

!-------------------------------------------------------------------------------
! Put all snow onto the ground and zero canopy snow.
! Currently all snow is held in snow_tile.
! For can_model=4 tiles, put snow into snow_grnd and zero snow_tile.
!-------------------------------------------------------------------------------

! Save input value.
  snow_on_ground(:,:) = snow_tile(:,:)

! Initialise stores to zero.
  snow_grnd(:,:) = 0.0
  snow_tile(:,:) = 0.0

! Initialise other variables with values that will be retained where there is
! no tile - using "sensible" values for when these are printed.
  snowDepth(:,:) = 0.0
  IF ( nsmax < 1 ) THEN
    rho_snow_grnd(:,:) = rho_snow_const
  ELSE
    rho_snow_grnd(:,:) = rho_snow_fresh
    tsnow(:,:,:) = 273.15
    ds(:,:,:) = 0.0
    IF ( l_spec_albedo ) rgrainL(:,:,:) = 0.0
  END IF

  DO n = 1,ntiles
    IF ( canSnowTile(n) ) THEN
      DO j = 1,tile_pts(n)
        i = tile_index(j,n)
        snow_grnd(i,n) = snow_on_ground(i,n)
      END DO
    ELSE
      DO j = 1,tile_pts(n)
        i = tile_index(j,n)
        snow_tile(i,n) = snow_on_ground(i,n)
      END DO
    END IF
  END DO

!-------------------------------------------------------------------------------
! Set snow density, calculate snow depth and set temperature of snow to equal
! that of soil.
!-------------------------------------------------------------------------------
  DO n=1,ntiles
    DO j=1,tile_pts(n)
      i = tile_index(j,n)
!     Use the constant (snowpack) density for nsmax=0 and if there is an
!     existing pack. If nsmax>0 and there is no pack, initialise the density
!     to the fresh snow value so that this value is used when/if a snowpack
!     next develops.
      IF ( nsmax == 0 .OR.                                                    &
           ( snow_on_ground(i,n) > EPSILON(snow_on_ground) ) ) THEN
        rho_snow_grnd(i,n) = rho_snow_const
      ELSE
        rho_snow_grnd(i,n) = rho_snow_fresh
      END IF
      snowDepth(i,n) = snow_on_ground(i,n) / rho_snow_grnd(i,n)
      IF ( nsmax > 0 ) THEN
        tsnow(i,n,:) = t_soil(i,1)
        IF ( l_spec_albedo ) rgrainl(i,n,:) = rgrain(i,n)
      END IF
    END DO
  END DO

!-------------------------------------------------------------------------------
! Calculate snow layer thicknesses.
!-------------------------------------------------------------------------------
  DO n=1,ntiles
    CALL layersnow(land_pts, tile_pts(n), tile_index(:,n), snowdepth(:,n),    &
                                                        nsnow(:,n), ds(:,n,:))
  END DO

!-------------------------------------------------------------------------------
! Set layer frozen and liquid contents.
!-------------------------------------------------------------------------------
  IF ( nsmax > 0 ) THEN
    sice(:,:,:) = 0.0
    sliq(:,:,:) = 0.0
    DO n = 1,ntiles
      DO j = 1,tile_pts(n)
        i = tile_index(j,n)
        DO k = 1,nsnow(i,n)
          sice(i,n,k) =  snow_on_ground(i,n) * ds(i,n,k) / snowdepth(i,n)
        END DO
      END DO
    END DO
  END IF

  RETURN

END SUBROUTINE total_snow_init
#endif
