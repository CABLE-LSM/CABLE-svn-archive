#if !defined(UM_JULES)
!**********************************************************************
! Routine to calculate the normalised solar radiation at each time and
! the time of the daily maximum temperature (GMT).
!
! Written by Peter Cox (March 1996)
!**********************************************************************
  SUBROUTINE SUNNY(                                               &
    DAYNUMBER,JDAY,POINTS,YEAR,LAT,LONG,SUN,TIME_MAX              &
  )

    USE c_pi

    USE datetime_mod, ONLY : SECS_IN_DAY

    IMPLICIT NONE

    INTEGER ::                                                    &
      DAYNUMBER,                                                  &
                ! IN Day of the year.
      JDAY,                                                       &
                ! IN Number of timesteps in the day.
      POINTS,                                                     &
                ! IN Number of spatial points.
      YEAR      ! IN Calender year.

    REAL ::                                                       &
      LAT(POINTS),                                                &
                ! IN Latitude (degrees).
      LONG(POINTS),                                               &
                ! IN Longitude (degrees).
      SUN(POINTS,JDAY),                                           &
                ! OUT Normalised solar radiation at each time.
      TIME_MAX(POINTS),                                           &
                ! OUT GMT at which temperature is maximum (hrs).
      COSDEC,                                                     &
                ! WORK COS (solar declination).
      COSLAT,                                                     &
                ! WORK COS (latitude).
      COSZ(POINTS),                                               &
                ! WORK Timestep mean COSZ.
      COSZM(POINTS),                                              &
                ! WORK Daily mean COSZ.
      LATRAD,                                                     &
                ! WORK Latitude (radians).
      LIT(POINTS),                                                &
                ! WORK Sunlit fraction of the day
      LONGRAD(POINTS),                                            &
                ! WORK Longitude (radians).
      SCS,                                                        &
                ! WORK Factor for TOA solar.
      SINDEC,                                                     &
                ! WORK SIN (solar declination).
      SINLAT(POINTS),                                             &
                ! WORK SIN (latitude).
      TANDEC,                                                     &
                ! WORK TAN (solar declination).
      TANLAT,                                                     &
                ! WORK TAN (latitude).
      TANTAN,                                                     &
                ! WORK TANDEC*TANLAT.
      OMEGA_UP,OMEGA_DOWN,                                        &
                ! WORK Solar angle of sunrise and sunset (radians).
      TIME_UP,TIME_DOWN,                                          &
                ! WORK GMT of sunrise and sunset (hrs).
      TIMESTEP,                                                   &
                ! WORK Timestep (s).
      TIME

    INTEGER :: I,J    ! WORK Loop counter.


    CALL SOLPOS (DAYNUMBER, YEAR, SINDEC, SCS)

    DO I=1,POINTS
      LATRAD     = PI_OVER_180*LAT(I)
      LONGRAD(I) = PI_OVER_180*LONG(I)
      SINLAT(I)  = SIN(LATRAD)
      COSZM(I)   = 0.0
    ENDDO

    COSDEC   = SQRT(1 - SINDEC**2)
    TANDEC   = SINDEC / COSDEC
    TIMESTEP = REAL(SECS_IN_DAY) / REAL(JDAY)

!----------------------------------------------------------------------
! Calculate the COSZ at each time
!----------------------------------------------------------------------
    DO J=1,JDAY

      TIME = (J-1) * TIMESTEP
      CALL SOLANG(                                                &
        SINDEC,TIME,TIMESTEP,SINLAT,LONGRAD,POINTS,LIT,COSZ       &
      )

      DO I=1,POINTS
        SUN(I,J) = COSZ(I) * LIT(I)
        COSZM(I) = COSZM(I) + SUN(I,J) / REAL(JDAY)
      ENDDO

    ENDDO

!----------------------------------------------------------------------
! Calculate the normalised solar radiation
!----------------------------------------------------------------------
    DO J=1,JDAY
      DO I=1,POINTS

        IF (COSZM(I) > EPSILON(1.0)) THEN
          SUN(I,J) = SUN(I,J) / COSZM(I)
        ELSE
          SUN(I,J)=0.0
        ENDIF

      ENDDO
    ENDDO

!----------------------------------------------------------------------
! Calculate the time of maximum temperature. Assume this occurs 0.15
! of the daylength after local noon (guess !).
!----------------------------------------------------------------------
    DO I=1,POINTS

      COSLAT = SQRT(1 - SINLAT(I)**2)
      TANLAT = SINLAT(I) / COSLAT
      TANTAN = TANLAT * TANDEC

      IF (ABS(TANTAN) <= 1.0) THEN      ! Sun sets and rises

        OMEGA_UP   = -ACOS(-TANTAN)
        TIME_UP    = 0.5 * 24.0                                   &
                   * ((OMEGA_UP - LONGRAD(I)) / PI + 1)
        OMEGA_DOWN = ACOS(-TANTAN)
        TIME_DOWN  = 0.5 * 24.0                                   &
                   * ((OMEGA_DOWN - LONGRAD(I)) / PI + 1)

      ELSE                             ! Perpertual day or night
        TIME_UP   = 0.0
        TIME_DOWN = 0.0
      ENDIF

      TIME_MAX(I) = 0.5 * (TIME_UP + TIME_DOWN)                   &
                  + 0.15 * (TIME_DOWN - TIME_UP)

    ENDDO

    RETURN

  END SUBROUTINE SUNNY
#endif
