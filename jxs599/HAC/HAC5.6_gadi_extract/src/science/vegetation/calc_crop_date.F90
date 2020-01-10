! *****************************COPYRIGHT**************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use
! and distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
!
! [Met Office Ref SC0237]
! ****************************COPYRIGHT***************************************
! Optimum crops (non-rice) and plant date for each irrigated grid box.

SUBROUTINE calc_crop_date(land_index, land_pts, row_length, rows,             &
                          nsurft, frac,                                       &
                          sw_surft, tstar_surft, lw_down, tl_1,               &
                          con_rain, ls_rain, con_snow, ls_snow,               &
                          plant_n, nday_crop)

! Description:
!   Calculates average temperature, precipitation and net radiation
!     over NYAV years to determine best planting date for
!     non-rice crops for use in irrigation code
!
! Method:
!
! Notes hadrd 2011-07-15:
!   This subroutine was originally included within control.f90
!     but is now put in a separate file
!   Alternatively, it may be included in a separate irrigation subroutine?
! Notes hadrd 2011-09-30:
!   Variables that need to be saved between timesteps are now in
!   module irrcrop_ctl
!     and allocated at the start (and de-allocated at the end)
!     when <irr_crop> is set to 1 in the run control file
!
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

USE conversions_mod, ONLY:                                                    &
  ! imported scalar parameters
  isec_per_day

USE csigma, ONLY:                                                             &
  ! imported scalar parameters
  sbcon ! Stefan-Boltzmann constant (W/m**2/K**4).

USE crop_vars_mod, ONLY:                                                      &
  ! imported scalars
  iyear_old,nyav,ndpy,                                                        &
  ! imported arrays
  icntmax_gb,                                                                 &
  tl_1_day_av_gb, tl_1_day_av_use_gb,                                         &
  prec_1_day_av_gb, prec_1_day_av_use_gb,                                     &
  rn_1_day_av_gb, rn_1_day_av_use_gb,                                         &
  startyr, startmon, startday, starttime

USE jules_surface_types_mod, ONLY:                                            &
  ! imported scalars
  ntype

USE time_info_mod, ONLY: l_360, l_leap, current_model_time

USE datetime_utils_mod, ONLY: day_of_year, days_in_year

USE timestep_mod, ONLY: timestep

!-----------------------------------------------------------------------------

IMPLICIT NONE

! Scalar arguments with intent(IN) :
INTEGER, INTENT(IN) :: land_pts,row_length,rows,nsurft
INTEGER, INTENT(IN) :: nday_crop

! Array arguments with intent(IN) :
INTEGER, INTENT(IN) :: land_index(land_pts)
REAL, INTENT(IN) :: frac(land_pts,ntype)        ! Fractions of surface types

REAL, INTENT(IN) :: sw_surft(land_pts,nsurft)   ! Surface net SW radiation on
                                                ! land tiles (W/m2).

REAL, INTENT(IN) :: tstar_surft(land_pts,nsurft) ! Surface tile temperatures
                                                 ! (K).

REAL, INTENT(IN) :: lw_down(row_length,rows)    ! Surface downward LW
                                                ! radiation (W/m2).
REAL, INTENT(IN) :: tl_1(row_length,rows)       ! Liquid/frozen water
                                                ! temperature for lowest
                                                ! atmospheric layer (K).
REAL, INTENT(IN) :: con_rain(row_length,rows)   ! Convective rain (kg/m2/s).
REAL, INTENT(IN) :: ls_rain(row_length,rows)    ! Large-scale rain (kg/m2/s).
REAL, INTENT(IN) :: con_snow(row_length,rows)   ! Convective snow (kg/m2/s).
REAL, INTENT(IN) :: ls_snow(row_length,rows)    ! Large-scale snow (kg/m2/s).

! Array arguments with intent(OUT) :
INTEGER, INTENT(OUT) :: plant_n(land_pts)       ! best plant date for non-rice

! LOCAL scalar variables :
INTEGER :: i,j,l,n                  ! Loop counters
INTEGER :: tspd                     ! timesteps per day

INTEGER :: iyear
INTEGER :: day_of_yr                ! current day of year (== IJULIAN)
INTEGER :: secs_in_year             ! Number of seconds in the year.

INTEGER :: year, month, day, time_of_day
LOGICAL :: datetime_ne

! LOCAL array variables :
INTEGER :: icount(ndpy,nyav)
REAL :: rn(land_pts)                ! Net downward radiation (W m-2).

!-----------------------------------------------------------------------------
! End of header
!-----------------------------------------------------------------------------

CALL current_model_time(year, month, day, time_of_day)

tspd         = isec_per_day / timestep ! #timesteps per day
day_of_yr    = day_of_year(year, month, day, l_360, l_leap)
secs_in_year = ndpy * isec_per_day

!-----------------------------------------------------------------------------
! Calculate grid box mean net radiation for use in crop irrigation:
!-----------------------------------------------------------------------------

rn(:) = 0.0 ! first set RN to zero

DO n = 1,nsurft ! loop over tiles
  DO l = 1,land_pts
    rn(l) = rn(l) + frac(l,n) * ( sw_surft(l,n) & ! net SW per tile
            - sbcon * tstar_surft(l,n)**4.0 )     ! upward LW radiation
                                                  ! per tile
  END DO
END DO

! Add gridbox-average downward LW radiation
DO l = 1,land_pts
  j = (land_index(l) - 1) / row_length + 1
  i =  land_index(l) - (j-1) * row_length
  rn(l) = rn(l) + lw_down(i,j)
END DO

!-----------------------------------------------------------------------------
! Calculate averages and reset values at end of every cycle:
!-----------------------------------------------------------------------------
! nyav = nr of years averaged for crop plant estimates
! iyear varies between 1, 2, and 3 if nyav = 3
iyear = MOD ( ABS(year - startyr), nyav)

IF ( iyear == 0 ) iyear = nyav
IF ( day_of_yr <= ndpy ) THEN
  icount(day_of_yr,iyear) = icount(day_of_yr,iyear) + 1

  DO l = 1,land_pts
    j = (land_index(l) - 1) / row_length + 1
    i = land_index(l) - (j-1) * row_length

    ! Calculate average temperature, precip and net radiaton over NYAV years.
    ! First calculate average over 1 day.
    ! hadrd - daily averages are reset after 1 day and therefore do not need
    !         the dimensions (points, ijulian, iyear)
    tl_1_day_av_gb(l)   = tl_1_day_av_gb(l) + tl_1(i,j) / REAL(tspd)
    prec_1_day_av_gb(l) = prec_1_day_av_gb(l) +                               &
                          (con_rain(i,j) + con_snow(i,j)                      &
                          + ls_rain(i,j) + ls_snow(i,j))                      &
                          / REAL(tspd)
    rn_1_day_av_gb(l)   = rn_1_day_av_gb(l) + rn(l) / REAL(tspd)

  END DO

  IF ( icount(day_of_yr,iyear) == tspd) THEN ! after 1 day

    ! Save averages for cropping irrigation code:
    DO l = 1,land_pts
      tl_1_day_av_use_gb(l,day_of_yr,iyear)   = tl_1_day_av_gb(l)
      prec_1_day_av_use_gb(l,day_of_yr,iyear) = prec_1_day_av_gb(l)
      rn_1_day_av_use_gb(l,day_of_yr,iyear)   = rn_1_day_av_gb(l)
    END DO

    ! Reset to zero:
    ! hadrd - this is done already after 1 day
    icount(day_of_yr,iyear) = 0
    DO l = 1,land_pts
      tl_1_day_av_gb(l)     = 0.0
      prec_1_day_av_gb(l)   = 0.0
      rn_1_day_av_gb(l)     = 0.0
    END DO

  END IF

  !---------------------------------------------------------------------------
  ! Call cropping irrigation code at end of every YEAR of data stored:
  !---------------------------------------------------------------------------
  datetime_ne = ( year /= startyr )   .OR.                                    &
                ( month /= startmon ) .OR.                                    &
                ( day  /= startday )  .OR.                                    &
                ( time_of_day /= starttime )

  IF ( (iyear /= iyear_old) .AND. datetime_ne ) THEN
    ! When year has changed i.e. first day of new year

    CALL opt_crop_date(land_pts, ndpy, nyav, nday_crop,                       &
                       plant_n, icntmax_gb,                                   &
                       tl_1_day_av_use_gb,                                    &
                       prec_1_day_av_use_gb,                                  &
                       rn_1_day_av_use_gb)

  END IF

END IF
iyear_old = iyear

RETURN

END SUBROUTINE calc_crop_date
