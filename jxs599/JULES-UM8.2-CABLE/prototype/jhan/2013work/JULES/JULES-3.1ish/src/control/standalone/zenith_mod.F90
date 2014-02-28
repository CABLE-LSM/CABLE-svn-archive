#if !defined(UM_JULES)
MODULE ZENITH_MOD

  IMPLICIT NONE
! This module contains subroutine zenith and some supporting time/calendar code.

CONTAINS

!###############################################################################
!###############################################################################

  SUBROUTINE ZENITH(NPOINTS,COSZ)

! Subroutine to calculate the cosine of the zenith angle

    USE datetime_mod, ONLY : SECS_IN_DAY, days_in_year, day_of_year

    USE model_time_mod, ONLY : current_time

    USE c_pi, ONLY : pi,pi_over_180

    USE theta_field_sizes, ONLY : t_i_length

    USE latlon_mod, ONLY : latitude, longitude

    IMPLICIT NONE

    INTEGER                   :: NPOINTS  ! Number of points
    REAL, DIMENSION(NPOINTS)  :: COSZ     ! Cosine of zenith angle

    REAL    ::  HH     ! Hour angle of sun relative to solar noon (radians)
    REAL    ::  SD     ! Solar declination angle (radians)

    INTEGER ::  DAY    ! Current day
    INTEGER ::  MONTH  ! Current month
    INTEGER ::  YEAR   ! Current year
    INTEGER ::  time   ! Current time

    INTEGER :: I,J,L   ! Loop counters


!-----------------------------------------------------------------------------


! Split the current date into its parts for ease
    year  = current_time%year
    month = current_time%month
    day   = current_time%day
    time  = current_time%time

    DO L = 1,NPOINTS
      J = (L - 1) / t_i_length + 1
      I = L - (J - 1) * t_i_length

      HH = PI * (2.0 * REAL(time) / REAL(SECS_IN_DAY) + LONGITUDE(I,J) / 180.0 - 1.0)

      SD = -23.4 * PI_OVER_180 *                                              &
           COS(2.0 * PI * (day_of_year(year, month, day) + 10) / days_in_year(year))

      COSZ(L) = SIN(PI_OVER_180 * LATITUDE(I,J)) * SIN(SD)                    &
              + COS(PI_OVER_180 * LATITUDE(I,J)) * COS(SD) * COS(HH)

      IF ( COSZ(L) < 0.0 ) COSZ(L) = 0.0

    ENDDO

    RETURN

  END SUBROUTINE ZENITH

!###############################################################################
!###############################################################################

END MODULE zenith_mod
#endif
