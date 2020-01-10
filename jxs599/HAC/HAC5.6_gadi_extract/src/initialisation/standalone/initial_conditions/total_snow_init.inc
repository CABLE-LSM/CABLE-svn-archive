#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE total_snow_init()

USE layersnow_mod,   ONLY: layersnow

USE ancil_info, ONLY: land_pts, nsurft, surft_index, surft_pts,               &
                                         lice_index, lice_pts, nsoilt

USE prognostics, ONLY: nsnow_surft, rgrain_surft, snow_surft, t_soil_soilt,   &
                        rgrainl_surft, rho_snow_grnd_surft, sice_surft,       &
                        sliq_surft, snow_grnd_surft, snowdepth_surft,         &
                        tsnow_surft, ds_surft, tsurf_elev_surft

USE jules_snow_mod, ONLY: nsmax, rho_snow_const, rho_snow_fresh, canSnowTile

USE jules_surface_mod, ONLY: l_elev_land_ice

USE jules_radiation_mod, ONLY: l_snow_albedo, l_embedded_snow

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Sets up the initial snow conditions from snow_surft only
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER :: i, j, k, n, m  ! Index variables

REAL(KIND=real_jlslsm) :: snow_on_ground(land_pts,nsurft)
              ! Snow considered to be on the ground (not in canopy) (kg m-2)
              ! This is all the snow

!-------------------------------------------------------------------------------
! Put all snow onto the ground and zero canopy snow.
! Currently all snow is held in snow_surft.
! For can_model=4 tiles, put snow into snow_grnd_surft and zero snow_surft.
!-------------------------------------------------------------------------------

! Save input value.
snow_on_ground(:,:) = snow_surft(:,:)

! Initialise stores to zero.
snow_grnd_surft(:,:) = 0.0
snow_surft(:,:) = 0.0

! Initialise other variables with values that will be retained where there is
! no tile - using "sensible" values for when these are printed.
snowdepth_surft(:,:) = 0.0
IF ( nsmax < 1 ) THEN
  rho_snow_grnd_surft(:,:) = rho_snow_const
ELSE
  rho_snow_grnd_surft(:,:) = rho_snow_fresh
  tsnow_surft(:,:,:) = 273.15
  ds_surft(:,:,:) = 0.0
  IF ( l_snow_albedo .OR. l_embedded_snow ) rgrainl_surft(:,:,:) = 0.0
END IF

DO n = 1,nsurft
  IF ( canSnowTile(n) ) THEN
    DO j = 1,surft_pts(n)
      i = surft_index(j,n)
      snow_grnd_surft(i,n) = snow_on_ground(i,n)
    END DO
  ELSE
    DO j = 1,surft_pts(n)
      i = surft_index(j,n)
      snow_surft(i,n) = snow_on_ground(i,n)
    END DO
  END IF
END DO

!-------------------------------------------------------------------------------
! Set snow density, calculate snow depth and set temperature of snow to equal
! that of soil.
!-------------------------------------------------------------------------------

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

DO n = 1,nsurft

  !Set the current soil tile (see notice above)
  IF (nsoilt == 1) THEN
    !There is only 1 soil tile
    m = 1
  ELSE ! nsoilt == nsurft
    !Soil tiles map directly on to surface tiles
    m = n
  END IF !nsoilt

  DO j = 1,surft_pts(n)
    i = surft_index(j,n)
    !     Use the constant (snowpack) density for nsmax=0 and if there is an
    !     existing pack. If nsmax>0 and there is no pack, initialise the density
    !     to the fresh snow value so that this value is used when/if a snowpack
    !     next develops.
    IF ( nsmax == 0 .OR.                                                      &
         ( snow_on_ground(i,n) > EPSILON(snow_on_ground) ) ) THEN
      rho_snow_grnd_surft(i,n) = rho_snow_const
    ELSE
      rho_snow_grnd_surft(i,n) = rho_snow_fresh
    END IF
    snowdepth_surft(i,n) = snow_on_ground(i,n) / rho_snow_grnd_surft(i,n)
    IF ( nsmax > 0 ) THEN

      tsnow_surft(i,n,:) = t_soil_soilt(i,m,1)
      IF ( l_snow_albedo .OR. l_embedded_snow )                               &
        rgrainl_surft(i,n,:) = rgrain_surft(i,n)
    END IF
  END DO
END DO

!==============================================================================
!**END NOTICE REGARDING SOIL TILING**
!==============================================================================

IF (l_elev_land_ice) THEN
  DO j = 1,lice_pts
    i = lice_index(j)
    IF ( nsmax > 0 ) THEN
      DO n = 1,nsurft
        tsnow_surft(i,n,:) = tsurf_elev_surft(i,n)
      END DO
    END IF
  END DO
END IF

nsnow_surft(:,:) = 0
IF ( nsmax > 0 ) THEN
  !-------------------------------------------------------------------------------
  ! Calculate snow layer thicknesses.
  !-------------------------------------------------------------------------------
  DO n = 1,nsurft
    CALL layersnow(land_pts, surft_pts(n), surft_index(:,n),                  &
                   snowdepth_surft(:,n), nsnow_surft(:,n), ds_surft(:,n,:))
  END DO

  !-------------------------------------------------------------------------------
  ! Set layer frozen and liquid contents.
  !-------------------------------------------------------------------------------
  sice_surft(:,:,:) = 0.0
  sliq_surft(:,:,:) = 0.0
  DO n = 1,nsurft
    DO j = 1,surft_pts(n)
      i = surft_index(j,n)
      DO k = 1,nsnow_surft(i,n)
        sice_surft(i,n,k) =  snow_on_ground(i,n) * ds_surft(i,n,k)            &
                              / snowdepth_surft(i,n)
      END DO
    END DO
  END DO
END IF

RETURN

END SUBROUTINE total_snow_init
#endif
