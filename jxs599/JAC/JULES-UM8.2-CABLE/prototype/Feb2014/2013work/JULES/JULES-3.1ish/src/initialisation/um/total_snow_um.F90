#if defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Land

! UM version of JULES routine totalSnowInit

  SUBROUTINE total_snow_um(  land_pts      &
                            ,ntype         &
                            ,ntiles        &
                            ,frac          &
                            ,dsm_levels    &
                            ,l_spec_albedo &
                            ,rgrain        &
                            ,snow_tile     &
                            ,t_soil        &
                            ,snow_grnd     &
                            ,nsnow         &
                            ,rgrainl       &
                            ,rho_snow_grnd &
                            ,sice          &
                            ,sliq          &
                            ,snowdepth     &
                            ,ds            &
                            ,tsnow )

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     nsmax

  USE snow_param, ONLY :  &
!  imported scalars with intent(in)
     rho_snow_const  &
    ,rho_snow_fresh  &
!  imported arrays with intent(in)
    ,canSnowTile

  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook

  IMPLICIT NONE

! subroutine arguments
  INTEGER :: land_pts                       &
!
          ,ntype                            &
!
          ,ntiles                           &
!
          ,dsm_levels

  INTEGER ::                                &
    nsnow(land_pts,ntiles)

  REAL ::  frac(      land_pts,ntype)       &
!
          ,rgrain(    land_pts,ntiles)      &
!
          ,snow_tile( land_pts,ntiles)      &
!
          ,snow_grnd( land_pts,ntiles)      &
!
          ,t_soil(    land_pts,dsm_levels)  &
!
          ,snowdepth(    land_pts,ntiles)          &
!
          ,rho_snow_grnd(land_pts,ntiles)          &
!
          ,tsnow(        land_pts,ntiles,nsmax)    &
!
          ,rgrainl(      land_pts,ntiles,nsmax)    &
!
          ,sice(         land_pts,ntiles,nsmax)    &
!
          ,sliq(         land_pts,ntiles,nsmax)    &
!
          ,ds(           land_pts,ntiles,nsmax)    &
!
          ,rho_snow(     land_pts,ntiles,nsmax)


  LOGICAL :: l_spec_albedo


! Local scalar variables.
  INTEGER :: i,j,k,n         !  work/loop counters

! Local array variables.
  REAL :: snow_on_ground(land_pts,ntiles)  !  snow considered to be on the
!        ground (not in canopy) (kg m-2). This is all the snow.

  INTEGER :: tile_index(land_pts,ntype)       &
            ,tile_pts(           ntype)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('TOTAL_SNOW_UM',zhook_in,zhook_handle)

!-------------------------------------
! set up tile counting and indexing
!-------------------------------------
! DEPENDS ON: tilepts
 CALL tilepts(land_pts,frac,tile_pts,tile_index)

!-------------------------------------------------------------------------------
! Put all snow onto the ground and zero canopy snow.
! Currently all snow is NOT held in snow_tile.
! For can_model=4 tiles, put snow into snow_grnd and zero snow_tile.
!-------------------------------------------------------------------------------

! Save input value.

  DO n=1,ntiles
    IF ( canSnowTile(n) ) THEN
      DO j=1,tile_pts(n)
        i = tile_index(j,n)
        snow_on_ground(i,n) = snow_tile(i,n) + snow_grnd(i,n)
      ENDDO
    ELSE
      DO j=1,tile_pts(n)
        i = tile_index(j,n)
        snow_on_ground(i,n) = snow_tile(i,n)
      ENDDO
    ENDIF
  ENDDO

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

  DO n=1,ntiles
    IF ( canSnowTile(n) ) THEN
      DO j=1,tile_pts(n)
        i = tile_index(j,n)
        snow_grnd(i,n) = snow_on_ground(i,n)
      END DO
    ELSE
      DO j=1,tile_pts(n)
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
      IF ( nsmax < 1 .OR.                                                     &
           ( snow_on_ground(i,n) > EPSILON(snow_on_ground) ) ) THEN
        rho_snow_grnd(i,n) = rho_snow_const
      ELSE
        rho_snow_grnd(i,n) = rho_snow_fresh
      END IF
      snowDepth(i,n) = max(snow_on_ground(i,n),0.0) / rho_snow_grnd(i,n)
      IF ( nsmax > 0 ) THEN
        tsnow(i,n,:) = t_soil(i,1)
        IF ( l_spec_albedo ) rgrainl(i,n,:) = rgrain(i,n)
      END IF
    END DO
  END DO

!-------------------------------------------------------------------------------
! Calculate snow layer thicknesses.
!-------------------------------------------------------------------------------
  DO N=1,NTILES
! DEPENDS ON: layersnow
    CALL layersnow( land_pts,tile_pts(n),tile_index(:,n),snowdepth(:,n)  &
                   ,nsnow(:,n),ds(:,n,:) )
  END DO

!-------------------------------------------------------------------------------
! Set layer frozen and liquid contents.
!-------------------------------------------------------------------------------
  IF ( nsmax > 0 ) THEN
    sice(:,:,:) = 0.0
    sliq(:,:,:) = 0.0
    DO n=1,ntiles
      DO j=1,tile_pts(n)
        i = tile_index(j,n)
        DO k=1,nsnow(i,n)
          IF (snowdepth(i,n) > 0.0) THEN
            sice(i,n,k) =  snow_on_ground(i,n) * ds(i,n,k) / snowdepth(i,n)
          END IF
        END DO
      END DO
    END DO
  END IF

  IF (lhook) CALL dr_hook('TOTAL_SNOW_UM',zhook_out,zhook_handle)

  END SUBROUTINE total_snow_um
#endif
