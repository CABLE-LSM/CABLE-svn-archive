! *****************************COPYRIGHT**************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use
! and distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT**************************************

! Optimum crops (non-rice) and plant date for each irrigated grid box.

SUBROUTINE opt_crop_date(land_pts, ndpy, nyav, nday_crop,                     &
                         plant_n, icntmax,                                    &
                         tl_1_day_av_use,                                     &
                         prec_1_day_av_use,                                   &
                         epot_1_day_av_use)


! Description:
!   Calculates best planting date for non-rice
!   as a function of temperature, precipitation and
!     potential evaporation criteria over NYAV year
!   Assumes potential evap = RN
!   Note routine is only called after year has completed
!
! Method:
!
! Notes RD 2011-07-12:
!   This subroutine was originally included at the end of control.f90
!     but is now put in a separate file
!   Alternatively, it may be included in a separate irrigation subroutine?
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
  zerodegc  ! ZeroDegC (K)

USE water_constants_mod, ONLY:                                                &
  ! imported scalar parameters
  lc   ! latent heat of condensation of water at 0degc (J kg-1).

USE conversions_mod, ONLY:                                                    &
  ! imported scalar parameters
  rsec_per_day

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts
INTEGER, INTENT(IN) :: ndpy           ! No. of days per year
INTEGER, INTENT(IN) :: nyav           ! No. of years averaged for crop plant
                                      ! estimates
INTEGER, INTENT(IN) :: nday_crop      ! khalla - set to 150 in irrcrop_ctl

REAL(KIND=real_jlslsm), INTENT(IN) :: tl_1_day_av_use(land_pts,ndpy,nyav)
  ! Daily mean air temperature (K).
REAL(KIND=real_jlslsm), INTENT(IN) :: prec_1_day_av_use(land_pts,ndpy,nyav)
  ! Daily mean precipitation rate (kg m-2 s-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: epot_1_day_av_use(land_pts,ndpy,nyav)
  ! Daily mean net radiation (W m-2).

INTEGER, INTENT(OUT) :: plant_n(land_pts)  ! best plant date for non-rice
INTEGER, INTENT(OUT) :: icntmax(land_pts)  ! general T criteria for non-rice

! Time periods for non-rice criteria:
INTEGER, PARAMETER :: bday_n1 = 1, bday_n2 = 21, bday_n3 = 51, eday_n3 = 110

REAL(KIND=real_jlslsm), PARAMETER ::                                          &
  tmin_n   = zerodegc +  5.0,                                                 &
  tmin_n23 = zerodegc + 15.0,                                                 &
  tmax_n23 = zerodegc + 30.0

REAL(KIND=real_jlslsm), PARAMETER :: prec_thr = 0.5

INTEGER :: iend,iend1,iend2,ibeg
INTEGER :: i,l,ii,iy,j,icc
INTEGER :: max_day  ! Max no. of days which reach general crop criteria
INTEGER :: eday_n1, eday_n2

INTEGER ::                                                                    &
  ic(land_pts,ndpy),                                                          &
    ! Counter.
  icnt(land_pts,ndpy),                                                        &
    ! General T criteria for non-rice.
  icn(land_pts,ndpy),                                                         &
    ! Criteria for non-rice WHOLE growing period.
  icn1(land_pts,ndpy),                                                        &
    ! Criteria for non-rice period 1.
  icn2(land_pts,ndpy),                                                        &
    ! Criteria for non-rice period 2.
  icn3(land_pts,ndpy) ,                                                       &
    ! Criteria for non-rice period 3.
  plant_n_tmp(land_pts),                                                      &
    ! Best plant date for non-rice.
  iicn(land_pts)
    ! Counter for start date for non-rice.

REAL(KIND=real_jlslsm) ::                                                     &
  tl_av(land_pts,ndpy),                                                       &
    ! Multi-year average daily temperature (K).
  prec_av(land_pts,ndpy),                                                     &
    ! Multi-year average daily precipitation (mm).
  prec_tmp(land_pts,ndpy),                                                    &
    ! Precipitation (mm).
  epot_av(land_pts,ndpy),                                                     &
    ! Multi-year average daily potential evaporation (mm).
  epot_tmp(land_pts,ndpy)
    ! Potential evaporation (mm).

!-----------------------------------------------------------------------------
! End of header
!-----------------------------------------------------------------------------

! First set to zero
DO l = 1,land_pts
  DO i = 1,ndpy
    tl_av(l,i)   = 0.0
    prec_av(l,i) = 0.0
    epot_av(l,i) = 0.0
    icnt(l,i)    = 0
    icn1(l,i)    = 0
    icn2(l,i)    = 0
    icn3(l,i)    = 0
    ic(l,i)      = 0
  END DO
END DO

DO l = 1,land_pts    ! loop over land points
  DO iy = 1,nyav     ! loop over 3 years
    DO i = 1,ndpy    ! loop over each day in year

      IF ( tl_1_day_av_use(l,i,iy) > 0.0) THEN  ! If data

        ic(l,i) = ic(l,i) + 1
        ! Calculate sum of daily temp, precip and net rad for this day over
        ! NYAV years.
        tl_av(l,i)   = tl_av(l,i) + tl_1_day_av_use(l,i,iy)
        prec_av(l,i) = prec_av(l,i) + prec_1_day_av_use(l,i,iy)
        epot_av(l,i) = epot_av(l,i) + epot_1_day_av_use(l,i,iy)

        ! Dont allow negative potential evaporation:
        IF ( epot_av(l,i) < 0.0 ) epot_av(l,i) = 0.0

        prec_tmp(l,i) = 0.0
        epot_tmp(l,i) = 0.0

      END IF ! TL_1_DAY_AV_USE GT 0
    END DO ! NDPY
  END DO ! NYAV
END DO ! LAND_PTS

!-----------------------------------------------------------------------------
! Calculate 10 day avg for EPOT and PREC
!-----------------------------------------------------------------------------
DO l = 1,land_pts
  DO i = 1,ndpy
    icc = 0 ! counter for calculating averages

    DO j = -4,5 ! loop from -4 to +5 days

      ii = i + j ! new index for day number
      IF ( ii < 1 )    ii = ii + ndpy
      IF ( ii > ndpy ) ii = ii - ndpy

      IF ( tl_av(l,ii) > 0.0 ) THEN ! If data
        icc = icc + 1
        prec_tmp(l,i) = prec_tmp(l,i) + prec_av(l,ii)
        epot_tmp(l,i) = epot_tmp(l,i) + epot_av(l,ii)
      END IF

    END DO ! loop from -4 to +5 days

    IF ( icc > 0 ) THEN
      ! Calculate average over ICC days
      prec_tmp(l,i) = prec_tmp(l,i) / REAL(icc)
      epot_tmp(l,i) = epot_tmp(l,i) / REAL(icc)
    END IF

  END DO ! NDPY
END DO ! LAND_PTS

! Copy average to PREC_AV and EPOT_AV
DO l = 1,land_pts
  DO i = 1,ndpy
    prec_av(l,i) = prec_tmp(l,i)
    epot_av(l,i) = epot_tmp(l,i)
  END DO
END DO

!-----------------------------------------------------------------------------
! Calculate multi-year avg for temperature, prec and epot
!-----------------------------------------------------------------------------
DO l = 1,land_pts
  DO i = 1,ndpy

    IF ( ic(l,i) > 0 ) THEN
      tl_av(l,i)   = tl_av(l,i) / REAL(ic(l,i))
      prec_av(l,i) = prec_av(l,i) / REAL(ic(l,i)) * rsec_per_day
        ! convert to mm/day
      epot_av(l,i) = epot_av(l,i) / REAL(ic(l,i)) * rsec_per_day / lc
        ! convert from W/m2 to mm/day
    END IF

  END DO
END DO

!-----------------------------------------------------------------------------
! Non-Rice
! khalla - similar code for rice exists in r1951_irrig_crop
!-----------------------------------------------------------------------------

! Temperature threshold for whole of the growing season:
DO l = 1,land_pts
  DO i = 1,ndpy
    iend  = i + nday_crop - 1
    iend1 = iend
    iend2 = 0

    ! At end of year, calculation should continue with data from start of year
    IF ( iend > ndpy ) THEN
      iend1 = ndpy
      iend2 = nday_crop - 1 + i - iend1
    END IF

    ! Counter for when general T criteria for non-rice are met
    DO ii = i,iend1
      IF ( tl_av(l,ii) > tmin_n ) icnt(l,i) = icnt(l,i) + 1
    END DO
    DO ii = 1,iend2
      IF ( tl_av(l,ii) > tmin_n ) icnt(l,i) = icnt(l,i) + 1
    END DO

    icntmax(l) = MAXVAL( icnt(l,1:ndpy) )

  END DO ! NDPY

  ! Now consider components of the growing season:
  !------------------------------------------------
  ! Period 1: Days 1-20 avg(Prec)>0.5 avg(Epot)
  eday_n1 = bday_n2 - 1 ! hadrd - BDAY_N2 set to 21 in header

  DO i = 1,ndpy
    ibeg = i + bday_n1 - 1 ! hadrd - BDAY_N1 set to 1 in header
    iend = i + eday_n1 - 1

    ! At end of year, calculation should continue with data from start of year
    IF ( ibeg > ndpy ) THEN
      ibeg = ibeg - ndpy
      iend = iend - ndpy
    END IF

    iend1 = iend
    iend2 = 0
    IF ( (iend > ndpy) .AND. (ibeg <= ndpy) ) THEN
      ! Calculation is split between remainder of year + start of year
      iend1 = ndpy
      iend2 = iend - iend1
    END IF

    DO ii = ibeg,iend1
      IF ( prec_av(l,ii) > prec_thr * epot_av(l,ii) ) THEN
        icn1(l,i) = icn1(l,i) + 1
      END IF
    END DO

    DO ii = 1,iend2
      IF ( prec_av(l,ii) > prec_thr * epot_av(l,ii) ) THEN
        icn1(l,i) = icn1(l,i) + 1
      END IF
    END DO
  END DO

  !------------------------------------------------
  ! Period 2: Days 21-50 T 18-30, avg(Prec)>0.5 avg(Epot)
  eday_n2 = bday_n3 - 1 ! hadrd - BDAY_N3 set to 51 in header

  DO i = 1,ndpy
    ibeg = i + bday_n2 - 1
    iend = i + eday_n2 - 1

    ! At end of year, calculation should continue with data from start of year
    IF ( ibeg > ndpy ) THEN
      ibeg = ibeg - ndpy
      iend = iend - ndpy
    END IF

    iend1 = iend
    iend2 = 0

    IF ( (iend > ndpy) .AND. (ibeg <= ndpy) ) THEN
      ! Calculation is split between remainder of year + start of year
      iend1 = ndpy
      iend2 = iend - iend1
    END IF

    DO ii = ibeg,iend1
      IF ( (tl_av(l,ii) >= tmin_n23) .AND. (tl_av(l,ii) <= tmax_n23) ) THEN
        icn2(l,i) = icn2(l,i) + 1
      END IF
      IF ( prec_av(l,ii) > prec_thr * epot_av(l,ii) ) THEN
        icn2(l,i) = icn2(l,i) + 1
      END IF
    END DO

    DO ii = 1,iend2
      IF ( (tl_av(l,ii) >= tmin_n23) .AND. (tl_av(l,ii) <= tmax_n23) ) THEN
        icn2(l,i) = icn2(l,i) + 1
      END IF
      IF ( prec_av(l,ii) > prec_thr * epot_av(l,ii) ) THEN
        icn2(l,i) = icn2(l,i) + 1
      END IF
    END DO
  END DO

  !------------------------------------------------
  ! Period 3: Days 51-110 T 15-30
  DO i = 1,ndpy
    ibeg = i + bday_n3 - 1
    iend = i + eday_n3 - 1

    ! At end of year, calculation should continue with data from start of year
    IF ( ibeg > ndpy ) THEN
      ibeg = ibeg - ndpy
      iend = iend - ndpy
    END IF

    iend1 = iend
    iend2 = 0

    IF ( iend > ndpy .AND. ibeg <= ndpy ) THEN
      ! Calculation is split between remainder of year + start of year
      iend1 = ndpy
      iend2 = iend - iend1
    END IF

    DO ii = ibeg,iend1
      IF ( (tl_av(l,ii) >= tmin_n23) .AND. (tl_av(l,ii) <= tmax_n23) ) THEN
        icn3(l,i) = icn3(l,i) + 1
      END IF
    END DO
    DO ii = 1,iend2
      IF ( (tl_av(l,ii) >= tmin_n23) .AND. (tl_av(l,ii) <= tmax_n23) ) THEN
        icn3(l,i) = icn3(l,i) + 1
      END IF
    END DO
  END DO

  !------------------------------------------------
  ! Combine criteria over WHOLE growing period
  DO i = 1,ndpy
    icn(l,i) = icn1(l,i) + icn2(l,i) + icn3(l,i)
  END DO

  !------------------------------------------------
  ! Now eliminate growing periods when T is not over minimum threshold
  ! (5C) for non-rice for 150 consecutive days:

  max_day        = MAXVAL(icnt(l,:)) ! max number of days that temp
                                     ! criterion is met
  plant_n_tmp(l) = 0
  iicn(l)        = -1000

  DO i = 1,ndpy
    IF ( icnt(l,i) >= max_day ) THEN ! temperature threshold
      IF ( icn(l,i) > iicn(l) ) THEN ! combined criteria
        plant_n_tmp(l) = i
        iicn(l)        = icn(l,i)
      END IF
    END IF
  END DO

  plant_n(l) = plant_n_tmp(l) ! best plant day for non-rice

END DO ! end l land points loop

END SUBROUTINE opt_crop_date
