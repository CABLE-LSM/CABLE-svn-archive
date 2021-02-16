! *****************************COPYRIGHT*************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use
! and distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT*************************************
!
! Description:
! Calculate irrigation demand for each irrigated grid box.
!-----------------------------------------------------------------------------

MODULE irrig_dmd_mod

CONTAINS

!###############################################################################

SUBROUTINE irrig_dmd ( land_pts, sm_levels, frac_irr_soilt,                   &
                       a_step, plant_n_gb,                                    &
                       sthf_soilt, smvccl_soilt, smvcst_soilt, smvcwt_soilt,  &
                       sthzw_soilt, sthu_irr_soilt, sthu_soilt,               &
                       smcl_soilt, irr_crop, dvi_cpft,                        &
                       !New arguments replacing USE statements
                       !jules_rivers_mod
                       rivers_sto_per_m2_on_landpts, rivers_adj_on_landpts,   &
                       !crop_vars_mod
                       nday_crop, irrDaysDiag_gb, irrig_water_gb, dvimax_gb)
! Description:
!   Calculates irrigation demand over unfrozen soils as the amount of
!     water needed to alleviate soil water deficit
!
! Method:
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in VEGETATION
! Code originally developed by Nic Gedney
!
! Code Description:
!   Language: Fortran 90.
!   This code is partially modified to JULES coding standards v1.
!   Rutger Dankers, July 2011
!
!-----------------------------------------------------------------------------

USE water_constants_mod, ONLY: rho_water

USE conversions_mod, ONLY: secs_in_day=>isec_per_day

USE jules_soil_mod, ONLY: dzsoil

USE jules_irrig_mod, ONLY: l_irrig_dmd, l_irrig_limit

USE time_info_mod, ONLY: l_360, l_leap, current_model_time

USE datetime_utils_mod, ONLY: day_of_year, days_in_year

USE timestep_mod, ONLY: timestep

USE ancil_info, ONLY: nsoilt

USE ereport_mod, ONLY: ereport

USE jules_hydrology_mod, ONLY: zw_max


USE jules_surface_types_mod,  ONLY: ncpft

!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Scalar arguments with intent(IN) :
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
    ! Number of land ice points.
  sm_levels,                                                                  &
    ! No. of soil moisture levels.
  a_step,                                                                     &
    ! Atmospheric timestep number
    ! necessary when called from control.f90 each timestep
  irr_crop
    ! Switch for irrigation cropping model.

! Array arguments with intent(IN) :
INTEGER, INTENT(IN) ::                                                        &
  plant_n_gb(land_pts)
    ! Best plant date for non-rice.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  frac_irr_soilt(land_pts,nsoilt),                                            &
    ! Irrigation fraction for this year.
  sthf_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Frozen soil moisture content of each layer as a fraction of saturation.
  smvccl_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Critical volumetric SMC (cubic m per cubic m of soil).
  smvcst_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Volumetric saturation point (m3/m3 of soil).
  smvcwt_soilt(land_pts,nsoilt,sm_levels),                                    &
    ! Volumetric wilting point (m3/m3 of soil).
   dvi_cpft(land_pts,ncpft)
    ! Development index for crop tiles

! Array arguments with intent(INOUT) :
REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                      &
  sthu_irr_soilt(land_pts,nsoilt,sm_levels),                                  &
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation in irrigated fraction
  sthu_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation.
  smcl_soilt(land_pts,nsoilt,sm_levels),                                      &
    ! Soil moisture content of each layer (kg/m2).
  sthzw_soilt(land_pts,nsoilt)
    ! Soil moist fraction in deep layer.

!New arguments replacing USE statements
!jules_rivers_mod
REAL(KIND=real_jlslsm), INTENT(IN) :: rivers_sto_per_m2_on_landpts(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: rivers_adj_on_landpts(land_pts)
!crop_vars_mod
INTEGER, INTENT(IN) :: nday_crop
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: irrDaysDiag_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: irrig_water_gb(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: dvimax_gb(land_pts)


! LOCAL scalars :
LOGICAL , PARAMETER :: irr_zw_all = .TRUE.
    ! Withdraw all available water from the groundwater store (T),
    ! or withdraw until the wilting point only (F).

INTEGER ::                                                                    &
  l,n,m,                                                                      &
    ! Loop counters.
  tspd,                                                                       &
    ! Timesteps per day.
  days_in_yr,                                                                 &
    ! Total number of days in current year.
  day_of_yr,                                                                  &
    ! Current day of year.
  day,                                                                        &
    ! Current day.
  month,                                                                      &
    ! Current month.
  year,                                                                       &
    ! Current year.
  errorstatus

REAL(KIND=real_jlslsm) ::                                                     &
  irrig_sum,                                                                  &
    ! Total irrigwater_levels_day over soil layers (kg/m2).
  irrig_lim,                                                                  &
    ! Total irrigwater_levels_day constrained by water in deep soil (kg/m2).
  irrig_max,                                                                  &
    ! Maximum available water for irrigation from deep soil (kg/m2).
  irrig_dif,                                                                  &
    ! Amount by which irrigwater_levels_day is decreased (kg/m2).
  sthu_tmp
    ! Dummy for sthu_irr_soilt.

LOGICAL :: l_julian_crop

! LOCAL arrays :
REAL(KIND=real_jlslsm) ::                                                     &
  irrigwater_levels_day(land_pts,sm_levels),                                  &
    ! Addition of irrigation water to each soil layer (kg/m2/day).
  smclsatzw_soilt(land_pts,nsoilt),                                           &
    ! Moisture content in deep layer at saturation (kg/m2).
  smclwiltzw_soilt(land_pts,nsoilt),                                          &
    ! Moisture content in deep layer at wilting point (kg/m2).
  smclzw_soilt(land_pts,nsoilt),                                              &
    ! Actual moisture content in deep layer (kg/m2).
  smclzw_rest_soilt(land_pts,nsoilt),                                         &
    ! Deep moisture content after extraction of irrigwater_levels_day (kg/m2).
  zdepth(0:sm_levels),                                                        &
    ! Lower soil layer boundary depth (m).
  irrig_riv(land_pts)
    ! Total irrigation water extracted from river routing storage (kg/m2).

LOGICAL :: l_irrigated(land_pts)
    ! Flag for cells being irrigated

!-----------------------------------------------------------------------------
! End of header
!-----------------------------------------------------------------------------

! dvi_cpft is a prognostic variable calculated only by the standalone crop model
! The crop model is required for this option and is not switched on in the UM
! irr_crop=0 is the only available option for UM JULES
! Calculates maximum dvi per grid cell
IF ( irr_crop == 2 ) THEN
  DO l = 1,land_pts
    dvimax_gb(l) = MAXVAL(dvi_cpft(l,:))
  END DO
END IF

CALL current_model_time(year, month, day)

days_in_yr = days_in_year(year, l_360, l_leap)
day_of_yr  = day_of_year(year, month, day, l_360, l_leap)

tspd = secs_in_day / timestep ! #timesteps per day
l_irrigated(:) = .FALSE.

! Initialise water for irrigation to be zero:
irrigwater_levels_day(:,:) = 0.0

! Initialise variables for constraining irrigation supply
IF ( l_irrig_limit ) THEN

  ! Calculate depth of soil column
  zdepth(:) = 0.0
  DO n = 1,sm_levels
    zdepth(n) = zdepth(n-1) + dzsoil(n)
  END DO

  ! Calculate water availability in deep GW store
  smclsatzw_soilt(:,:)  = rho_water * smvcst_soilt(:,:,sm_levels)             &
                          * (zw_max - zdepth(sm_levels))
  smclwiltzw_soilt(:,:) = rho_water * smvcwt_soilt(:,:,sm_levels)             &
                          * (zw_max - zdepth(sm_levels))
  smclzw_soilt(:,:)     = sthzw_soilt(:,:) * smclsatzw_soilt(:,:)
  ! initialise other variables
  smclzw_rest_soilt(:,:)   = smclzw_soilt(:,:)
  rivers_adj_on_landpts(:) = 1.0
  irrig_riv(:)             = 0.0

END IF ! l_irrig_limit

!-----------------------------------------------------------------------------
! Calculate irrigation demand
!-----------------------------------------------------------------------------

! Calculate demand for irrigated water at start of every day:
IF ( MOD(a_step, tspd) == 0 .AND. ( a_step > 1 ) ) THEN
  !Loop over soil tiles
  DO m = 1, nsoilt
    ! Loop over land points
    DO l = 1, land_pts

      ! Determine if irrigation should be applied
      l_julian_crop = .FALSE.

      SELECT CASE ( irr_crop )

      CASE ( 0 )
        ! Continuous irrigation
        ! Effectively the same as if crop season lasts the whole year.
        l_julian_crop = .TRUE.

      CASE ( 1 )
        ! original crop model based on Doell & Siebert, 2002
        ! khalla- If today is before harvest date (plant_n_gb + 150) and
        ! after plant_n_gb
        IF ( (day_of_yr - plant_n_gb(l) < nday_crop) .AND.                    &
             (day_of_yr - plant_n_gb(l) >= 0) ) THEN
          l_julian_crop = .TRUE.
        END IF

        ! khalla- If plant_n_gb is less than 150 days before end of year and
        !         today is early in year, before harvest date
        IF ( (plant_n_gb(l) > days_in_yr - nday_crop) .AND.                   &
             (day_of_yr <= plant_n_gb(l) + nday_crop - days_in_yr) ) THEN
          l_julian_crop = .TRUE.
        END IF

        IF ( plant_n_gb(l) <= 0 ) l_julian_crop = .FALSE.
                                                ! crop code not yet called

      CASE ( 2 )
        ! use JULES-crop development index
        ! for the time being, the maximum dvi across all crop tiles is used
        ! to trigger irrigation. However, it may be better to vary the
        ! irrigated fraction according to which tiles have a suitable dvi
        IF ( dvimax_gb(l) > -1.0 .AND. dvimax_gb(l) < 2.0 ) THEN
          l_julian_crop = .TRUE.
        END IF

      CASE DEFAULT
        errorstatus = 101
        CALL ereport('irrig_dmd', errorstatus, 'Invalid value for irr_crop')

      END SELECT

      ! Irrigate if crop is planted and irrigation fraction > 0
      IF ( l_julian_crop .AND. frac_irr_soilt(l,m) > 0.0 ) THEN

        ! Irrigate top two soil layers
        DO n = 1,2       !SM_LEVELS
          ! Only irrigate if there is no frozen soil:
          IF ( sthf_soilt(l,m,n) <= 0.0 ) THEN
            IF ( sthu_irr_soilt(l,m,n) <                                      &
                 (smvccl_soilt(l,m,n) / smvcst_soilt(l,m,n)) ) THEN
              sthu_tmp              = sthu_irr_soilt(l,m,n)
              sthu_irr_soilt(l,m,n) = smvccl_soilt(l,m,n)                     &
                                      / smvcst_soilt(l,m,n)

              ! Ensure that irrigated soil moisture is less than saturation:
              sthu_irr_soilt(l,m,n) = MIN(sthu_irr_soilt(l,m,n),              &
                                          1.0 - sthf_soilt(l,m,n) )

              ! Ensure that gridbox mean soil moisture is less than saturation:
              sthu_irr_soilt(l,m,n) = MIN(sthu_irr_soilt(l,m,n),              &
                                          (( 1.0 - sthf_soilt(l,m,n)          &
                                             - sthu_soilt(l,m,n))             &
                                             / frac_irr_soilt(l,m) + sthu_tmp)&
                                         )

              irrigwater_levels_day(l,n) = (sthu_irr_soilt(l,m,n) - sthu_tmp) &
                                           * rho_water * dzsoil(n)            &
                                           * smvcst_soilt(l,m,n)              &
                                           * frac_irr_soilt(l,m)

              sthu_soilt(l,m,n) = sthu_soilt(l,m,n) + frac_irr_soilt(l,m)     &
                                  * (sthu_irr_soilt(l,m,n) - sthu_tmp)

              IF ( sthu_soilt(l,m,n) + sthf_soilt(l,m,n) > 1.0 ) THEN
                errorstatus = 101
                CALL ereport('irrig_dmd', errorstatus, 'GDM super saturation')
              END IF

              ! Update gridbox moisture content
              smcl_soilt(l,m,n) = (sthu_soilt(l,m,n) + sthf_soilt(l,m,n))     &
                                  * rho_water * dzsoil(n)                     &
                                  * smvcst_soilt(l,m,n)

            END IF   ! ( sthu_irr <  (smvccl / smvcst) )
          END IF   ! ( sthf_soilt(l,m,N) <= 0.0 )

        END DO ! SM_LEVELS

        !---------------------------------------------------------------------
        ! Constrain irrigation supply if required
        !---------------------------------------------------------------------

        IF ( l_irrig_limit ) THEN

          irrig_sum = 0.0
          irrig_max = 0.0
          irrig_lim = 0.0
          irrig_dif = 0.0

          ! total irrigation water added
          ! note that irrigwater_levels_day is already multiplied by
          ! frac_irr_soilt
          ! i.e. units are kg/m2 for entire grid box
          DO n = 1,sm_levels
            irrig_sum = irrig_sum + irrigwater_levels_day(l,n)
          END DO

          IF ( irrig_sum > 0.0 ) THEN
            ! amount of water added should not be more than available
            ! from deep groundwater store
            IF ( .NOT. irr_zw_all ) THEN
              ! option (1) withdraw water until the wilting point
              irrig_max = MAX(smclzw_soilt(l,m) - smclwiltzw_soilt(l,m), 0.0)
            ELSE
              ! option (2) withdraw all available water
              irrig_max = MAX( smclzw_soilt(l,m), 0.0)
            END IF

            irrig_lim = MIN(irrig_max, irrig_sum)
            smclzw_rest_soilt(l,m) = smclzw_soilt(l,m) - irrig_lim
                                                  !min(irrig_max,irrig_sum)

            ! re-calculate soil moisture fraction in deep layer
            sthzw_soilt(l,m) = MAX(smclzw_rest_soilt(l,m)                     &
                                   /smclsatzw_soilt(l,m)                      &
                                   ,0.0 )

            ! if irrigation is constrained, try extracting from river
            ! routestore
            IF ( rivers_sto_per_m2_on_landpts(l) > 0.0) THEN
              irrig_riv(l) = MAX(irrig_sum - irrig_lim, 0.0)
              irrig_riv(l) = MIN(rivers_sto_per_m2_on_landpts(l),             &
                                 irrig_riv(l))
              rivers_adj_on_landpts(l) = (rivers_sto_per_m2_on_landpts(l)     &
                                          - irrig_riv(l))                     &
                                         / rivers_sto_per_m2_on_landpts(l)
              irrig_lim = irrig_lim + irrig_riv(l)
            END IF ! rivers_sto_per_m2_on_landpts gt 0

            ! re-calculate soil moisture
            ! if irrigation is constrained, some or all of the extra water
            ! that was added should be subtracted again
            DO n = 1,sm_levels
              ! amount of water to be subtracted
              irrig_dif = (1.0 - irrig_lim / irrig_sum ) *                    &
                            irrigwater_levels_day(l,n)
              ! update soil moisture variables
              irrigwater_levels_day(l,n) = MAX(irrigwater_levels_day(l,n)     &
                                                - irrig_dif                   &
                                               , 0.0)
              sthu_irr_soilt(l,m,n) = sthu_irr_soilt(l,m,n)                   &
                                      - irrig_dif                             &
                                        / (rho_water * dzsoil(n)              &
                                           * smvcst_soilt(l,m,n)              &
                                           * frac_irr_soilt(l,m))
              smcl_soilt(l,m,n) = smcl_soilt(l,m,n) - irrig_dif
              sthu_soilt(l,m,n) = smcl_soilt(l,m,n)                           &
                                  / (rho_water * dzsoil(n)                    &
                                     * smvcst_soilt(l,m,n))                   &
                                  - sthf_soilt(l,m,n)
            END DO ! n
          END IF ! irrigsum gt 0
        END IF ! l_irrig_limit

      END IF   ! ( L_JULIAN_CROP .AND. frac_irr_soilt(l,m) >  0.0 )

      irrig_water_gb(l) = SUM(irrigwater_levels_day(l,:)) / REAL(secs_in_day)
      IF ( irrig_water_gb(l) > EPSILON(0.0) ) l_irrigated(l) = .TRUE.

    END DO ! LAND_PTS
  END DO ! Soil tiles
END IF ! (MOD(A_STEP,TSPD)

! number of days on which irrigation is applied
WHERE ( l_irrigated ) irrDaysDiag_gb = irrDaysDiag_gb + 1

RETURN

END SUBROUTINE irrig_dmd


END MODULE irrig_dmd_mod
