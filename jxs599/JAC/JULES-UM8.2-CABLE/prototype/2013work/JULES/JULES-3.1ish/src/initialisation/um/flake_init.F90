#if defined(UM_JULES) && defined(UM_FLAKE)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine FLAKE_INIT ----------------------------------
!
! Description: Initialisation of FLake.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
!   Code Owner: See Unified Model Code Owner's HTML page
!   This file belongs in section: Land


SUBROUTINE flake_init(  land_field,ntiles,sm_levels,land_index  &
                       ,frac_land                               &
                       ,tstar_tile,snow_tile,deep_soil_temp )


USE dyn_coriolis_mod, ONLY : f3_at_u

  USE nstypes, ONLY : ntype

  USE theta_field_sizes, ONLY : t_i_length

  USE lake_mod
  USE flake_parameters, ONLY : h_snow_min_flk   &
                              ,h_ice_min_flk    &
                              ,h_ice_max

! lake tile number
  USE nstypes, ONLY : lake

! constant snow density for FLake (could do better?)
  USE snow_param, ONLY : rho_snow_const

! get initial lake albedo from albsnc_nnvg, albsnf_nnvg
  USE nstypes,  ONLY : npft              &
                      ,nnvg
  USE nvegparm, ONLY : albsnc_nvg        &
                      ,albsnf_nvg

  USE c_0_dg_c, ONLY : tm

  USE soil_param, ONLY : dzsoil

IMPLICIT NONE

! Subroutine arguments

!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) ::                                  &
  land_field,                                           &
                 ! Number of land points
  ntiles,                                               &
                 ! Number of surface tiles
  sm_levels      ! Number of soil layers

INTEGER, INTENT(IN) ::                                  &
  land_index(land_field)

  REAL ::  frac_land(land_field,ntype)

  REAL ::  tstar_tile(land_field,ntiles)                &
          ,snow_tile( land_field,ntiles)                &
          ,deep_soil_temp(land_field,sm_levels)

! WORK variables:
  INTEGER :: i,j,l

  REAL :: soil_mean_temp


! set the module variables to initial values
        lake_depth(       :) = lake_depth_0
        lake_fetch(       :) = lake_fetch_0
        lake_t_mean(      :) = lake_t_mean_0
        lake_t_mxl(       :) = lake_t_mxl_0
        lake_t_ice(       :) = lake_t_ice_0
        lake_h_mxl(       :) = lake_h_mxl_0
        lake_h_ice(       :) = lake_h_ice_0
        lake_shape_factor(:) = lake_shape_0

    DO l=1,land_field
      IF (frac_land(l,lake) > 0.0) THEN

! get the mean temperature of the soil column
      soil_mean_temp = SUM(deep_soil_temp(l,:)*dzsoil(:))  &
                      /SUM(dzsoil(:))

! set surface T based on T* of the lake tile
        lake_t_sfc(l) = tstar_tile(l, lake)

! set mixed-layer T based on the temperature of the 1st soil level
        lake_t_mxl(l) = MAX(deep_soil_temp(l, 1), tm + 0.10)

! set mean T to the mean temperature over the soil column,
! BUT constrain to between the mixed-layer temperature and freezing.
        lake_t_mean(l) = MIN(soil_mean_temp, lake_t_mxl(l) - 0.05)
        lake_t_mean(l) = MAX(lake_t_mean(l), tm + 0.05)

! initialise the ice thickness
        IF(soil_mean_temp             <= tm) THEN
          lake_h_ice(l) = lake_depth(l)
        ELSE IF (deep_soil_temp(l, 1) <= tm) THEN
          lake_h_ice(l) = lake_h_mxl(l)
        ELSE IF (tstar_tile(l, lake ) <= tm) THEN
          lake_h_ice(l) = 2.0 * h_ice_min_flk
        ELSE
          lake_h_ice(l) = 0.0
        END IF
! ...and bound by the Mironov & Ritter (2004) maximum
!    to avoid SQRT(-ve) in FLake_driver
        lake_h_ice(l) = MIN( lake_h_ice(l),h_ice_max )

! bound the mixed-layer depth by the available unfrozen depth
        lake_h_mxl(l) = MIN(lake_h_mxl(l),lake_depth(l)-lake_h_ice(l))
        lake_h_mxl(l) = MAX(lake_h_mxl(l),0.0)

! set the snow depth from the mass and density of snow on the tile,
! BUT depending on the presence of ice
        lake_h_snow(l) = 0.0
        IF(lake_h_ice(l) > 0.0)THEN
          lake_h_snow(l) = snow_tile(l, lake) / rho_snow_const
          IF(lake_h_snow(l) > EPSILON(0.0))THEN
            lake_h_snow(l) = MAX(lake_h_snow(l),2.0 * h_snow_min_flk)
          END IF
        END IF

! set the snow surface T based on T* of the lake tile
        lake_t_snow(l) = tstar_tile(l, lake)

! set the ice upper-boundary temperature depending on the
! presence of snow
        IF(lake_h_snow(l) <= h_snow_min_flk)THEN
          lake_t_ice(l) = tstar_tile(l, lake)
        ELSE
          lake_t_ice(l) = deep_soil_temp(l, 1)
        END IF

! put an upper bound on the ice and snow T's
        lake_t_snow(l) = MIN(lake_t_snow(l),tm - 0.05)
        lake_t_ice( l) = MIN(lake_t_ice( l),tm - 0.05)

! initialise the lake albedo depending on snow
        IF(lake_h_snow(l) <= h_snow_min_flk)THEN
          lake_albedo(l) = albsnf_nvg(npft+nnvg-lake)
        ELSE
          lake_albedo(l) = albsnc_nvg(npft+nnvg-lake)
        END IF

! initialise the Nusselt number
        nusselt(l) = 100.0

      ELSE
        lake_t_sfc( l) = 0.0
        lake_t_snow(l) = 0.0
        lake_t_ice( l) = 0.0
        lake_t_mean(l) = 0.0
        lake_t_mxl( l) = 0.0
        lake_h_ice( l) = 0.0
        lake_h_mxl( l) = 0.0
        lake_h_snow(l) = 0.0
        lake_albedo(l) = 0.0
        nusselt(    l) = 0.0
      END IF

    END DO

    DO l=1,land_field
      IF (frac_land(l,lake) > 0.0) THEN
        j=(land_index(l)-1)/t_i_length + 1
        i = land_index(l) - (j-1)*t_i_length

! set the Coriolis parameter : ABSOLUTE VALUE
!
! To get the value at theta points,
! average the adjacent values at u points.
!
        coriolis_param(l) = ABS( (f3_at_u(i,j)+f3_at_u(i-1,j))/2.0 )

      END IF
    END DO

RETURN

END SUBROUTINE flake_init
#endif
