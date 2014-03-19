#if !defined(UM_JULES)
!----------------------------------------------------------------------
! This routine calculates sub-daily variability
! C. Huntingford (April 2001) - based on earlier version by P. Cox
!----------------------------------------------------------------------
  SUBROUTINE DAY_CALC(                                            &
    POINTSM,STEP_DAY,DAY_MON,SW_L,PRECIP_L,T_L,DTEMP_DAY_L,LW_L,  &
    PSTAR_L,WIND_L,RH15M_L,SW_SD,T_SD,LW_SD,CONV_RAIN_SD,         &
    LS_RAIN_SD,LS_SNOW_SD,PSTAR_SD,WIND_SD,QHUM_SD,MONTH,IDAY,    &
    LAT,LONG,SEC_DAY,NSDMAX,SEED_RAIN                             &
  )

    USE c_pi, ONLY : pi

    IMPLICIT NONE

    INTEGER ::                                                    &
      DAY_MON,                                                    &
              !IN "Number of days in a month" (day)
      POINTSM,                                                    &
              !IN Maximum number of points in grid.
      STEP_DAY,                                                   &
              !IN WORK Calculated number of timesteps per day
      SEC_DAY,                                                    &
              !IN Number of seconds in day (sec)
      NSDMAX,                                                     &
              !IN Maximum possible number of sub-daily timesteps.
      SEED_RAIN(4),                                               &
              !IN Seeding numbers required to disaggregate
              !   the rainfall
      N_TALLY !WORK Counting up number of precipitation periods.

!-----------------------------------------------------------------------
! Single day arrays of incoming variables
!-----------------------------------------------------------------------
    REAL ::                                                       &
      SW_L(POINTSM),                                              &
              ! IN Daily values for downward shortwave
              !    radiation (W/m2)
      PRECIP_L(POINTSM),                                          &
              ! In Daily values of rain+snow (mm/day)
      T_L(POINTSM),                                               &
              ! IN Daily values of temperature (K)
      DTEMP_DAY_L(POINTSM),                                       &
              ! IN Daily values of diurnal temperature range (K)
      LW_L(POINTSM),                                              &
              ! IN Daily values of surface longwave
              !    radiation (W/m**2).
      PSTAR_L(POINTSM),                                           &
              ! IN Daily values of pressure at level 1 (Pa).
      WIND_L(POINTSM),                                            &
              ! IN Daily values of wind speed (m/s)
      RH15M_L(POINTSM)
              ! Humidity (%)

!-----------------------------------------------------------------------
! Single day (global) arrays of variables above, but for up to hourly
! periods. NOTE: DTEMP_DAY_SD does not exist as T_SD combines T_L
! and DTEMP_DAY_L
!-----------------------------------------------------------------------

    REAL ::                                                       &
      SW_SD(POINTSM,NSDMAX),                                      &
              ! OUT Sub-daily values for downward
              !     shortwave radiation (W/m2)
      T_SD(POINTSM,NSDMAX),                                       &
              ! OUT Sub-daily values of temperature (K)
      LW_SD(POINTSM,NSDMAX),                                      &
              ! OUT Sub-daily values of surface
              !     longwave radiation (W/m**2).
      PSTAR_SD(POINTSM,NSDMAX),                                   &
              ! OUT Sub-daily values of pressure at level 1 (Pa).
      WIND_SD(POINTSM,NSDMAX),                                    &
              ! OUT Sub-daily values of wind speed (m/s)
      QHUM_SD(POINTSM,NSDMAX),                                    &
              ! OUT Sub_daily humidity deficit calculated
              !     from RH15M (kg/kg)
      CONV_RAIN_SD(POINTSM,NSDMAX),                               &
              ! OUT Sub-daily convective rain (mm/day)
      LS_RAIN_SD(POINTSM,NSDMAX),                                 &
              ! OUT Sub-daily large scale rain (mm/day)
      LS_SNOW_SD(POINTSM,NSDMAX)
              ! OUT Sub-daily large scale snow (mm/day)

    INTEGER ::                                                    &
      N_EVENT(POINTSM,NSDMAX),                                    &
              ! 1: if rains/snows during timestep period
              ! 0: otherwise
      N_EVENT_LOCAL(NSDMAX)
              ! As N_EVENT, but for each gridpoint

    REAL ::                                                       &
      PREC_LOC(NSDMAX)
              ! Temporary value of rainfall for each gridbox.

    REAL ::                                                       &
      T_SD_LOCAL(POINTSM),                                        &
              ! WORK Intermediate calculation of temperature (K)
      PSTAR_SD_LOCAL(POINTSM),                                    &
              ! WORK Intermediate calculation of pressure (Pa)
      QS_SD_LOCAL(POINTSM)
              ! WORK Saturated humidity deficit associated with
              !      T_SD_LOCAL and PSTAR_SD_LOCAL

    INTEGER ::                                                    &
      ISTEP,                                                      &
              ! WORK Loop over sub-daily periods
      MONTH,                                                      &
              ! IN Month of interest
      IDAY,                                                       &
              ! IN Day number since beginning of month
      DAYNUMBER,                                                  &
              ! WORK Daynumber since beginning of year
      L       ! WORK Loop over land points

    REAL ::                                                       &
      SUN(POINTSM,NSDMAX),                                        &
              ! WORK Normalised solar radiation
      TIME_MAX(POINTSM),                                          &
              ! WORK GMT at which temperature is maximum (hours)
      LAT(POINTSM),                                               &
              ! IN Latitude (degrees)
      LONG(POINTSM),                                              &
              ! IN Longitude (degrees)
      TIME_DAY,                                                   &
              ! WORK Time of day (hours)
      TIMESTEP,                                                   &
              ! WORK Timestep (seconds)
      RANDOM_NUM_SD
              ! WORK Random number associated with rai

    REAL ::                                                       &
      TEMP_CONV,                                                  &
              ! WORK Temperature above which rainfall is convective (K)
      TEMP_SNOW
              ! WORK Temperature below which snow occurs (K)
    PARAMETER(TEMP_CONV = 293.15)
    PARAMETER(TEMP_SNOW = 275.15)

    REAL ::                                                       &
      INIT_HOUR_CONV_RAIN,                                        &
              !WORK Start of convective rain event (hour)
      INIT_HOUR_LS_RAIN,                                          &
              !WORK Start of large scale rain event (hour)
      INIT_HOUR_LS_SNOW,                                          &
              !WORK Start of large scale snow event (hour)
      DUR_CONV_RAIN,                                              &
              !WORK Start of convective rain event (hour)
      DUR_LS_RAIN,                                                &
              !WORK Duration of large scale rain event (hour)
      DUR_LS_SNOW,                                                &
              !WORK Start of large scale snow event (hour)
      END_HOUR_CONV_RAIN,                                         &
              !WORK End of convective rain event (hour)
      END_HOUR_LS_RAIN,                                           &
              !WORK End of large scale rain event (hour)
      END_HOUR_LS_SNOW,                                           &
              !WORK End of large scale snow event (hour)
      HOUREVENT,                                                  &
              !WORK Local variable giving hours during diurnal
              !     period for checking if precip. event occurs
      MAX_PRECIP_RATE,                                            &
              !WORK Maximum allowed precip. rate allowed within
              !     each sub-daily timestep (mm/day). (This only
              !     applies when STEP_DAY >= 2)
      PERIOD_LEN
              !WORK Length of period (hr)


!-----------------------------------------------------------------------
! Calculate the maximum precipitation rate. It is noted that 58 mm/day
! over 8 timesteps, and where all fell within a single 3 hour period
! caused numerical issues for MOSES. This corresponded to a rate of
! 464 mm/day during the 3-hour period. Hence, place a limit of 350
! mm/day.
!-----------------------------------------------------------------------
    PARAMETER(MAX_PRECIP_RATE = 350.0)

!-----------------------------------------------------------------------
! First check whether sub-daily calculations are required.
!-----------------------------------------------------------------------

!      DUR_CONV_RAIN = 2.0
!      DUR_LS_RAIN = 5.0
!      DUR_LS_SNOW = 5.0

    DUR_CONV_RAIN = 6.0
    DUR_LS_RAIN = 1.0
    DUR_LS_SNOW = 1.0

    PERIOD_LEN = 24.0 / FLOAT(STEP_DAY)

!-----------------------------------------------------------------------
! First check whether sub-daily calculations are required.
!-----------------------------------------------------------------------
    IF(STEP_DAY >= 2) THEN

!-----------------------------------------------------------------------
! Ensure that the durations are at least as long as a time period for
! the model to prevent solution "falling through gaps"
!-----------------------------------------------------------------------
      IF(DUR_CONV_RAIN <= PERIOD_LEN)                             &
        DUR_CONV_RAIN = PERIOD_LEN+1.0E-6
      IF(DUR_LS_RAIN <= PERIOD_LEN)                               &
        DUR_LS_RAIN = PERIOD_LEN+1.0E-6
      IF(DUR_LS_SNOW <= PERIOD_LEN)                               &
        DUR_LS_SNOW = PERIOD_LEN+1.0E-6

      TIMESTEP = FLOAT(SEC_DAY) / FLOAT(STEP_DAY)

!-----------------------------------------------------------------------
! Calculate the diurnal cycle in the SW radiation
!-----------------------------------------------------------------------
      DAYNUMBER = INT((MONTH-1.0) * REAL(DAY_MON)) + IDAY
      CALL SUNNY(                                                 &
        DAYNUMBER,STEP_DAY,POINTSM,1990,LAT,LONG,SUN,TIME_MAX     &
      )

!-----------------------------------------------------------------------
! Loop over timesteps
!-----------------------------------------------------------------------
      DO ISTEP=1,STEP_DAY

!-----------------------------------------------------------------------
! Calculate timestep values of the driving data.
!-----------------------------------------------------------------------
        TIME_DAY = (REAL(ISTEP) - 0.5) * TIMESTEP

        DO L=1,POINTSM
          T_SD(L,ISTEP) = T_L(L) + 0.5 * DTEMP_DAY_L(L) *         &
                   COS(2*PI*(TIME_DAY-3600.0*TIME_MAX(L))/SEC_DAY)
          LW_SD(L,ISTEP) = LW_L(L)                                &
                         * (4.0 * T_SD(L,ISTEP) / T_L(L) - 3.0)
          SW_SD(L,ISTEP) = SW_L(L) * SUN(L,ISTEP)
        ENDDO

!-----------------------------------------------------------------------
! Calculate timestep values of the driving data that is not split up
! into diurnal behaviour.
!-----------------------------------------------------------------------
        DO L=1,POINTSM
          PSTAR_SD(L,ISTEP) = PSTAR_L(L)
          WIND_SD(L,ISTEP)  = WIND_L(L)
        ENDDO

!-----------------------------------------------------------------------
! Check that humidity value is not greater than QSAT (but otherwise, QHU
! is not split up into diurnal behaviour).
!-----------------------------------------------------------------------
        DO L=1,POINTSM
          PSTAR_SD_LOCAL(L) = PSTAR_SD(L,ISTEP)
          T_SD_LOCAL(L)     = T_SD(L,ISTEP)
        ENDDO

        CALL QSAT(QS_SD_LOCAL,T_SD_LOCAL,PSTAR_SD_LOCAL,POINTSM)

        DO L=1,POINTSM
          QHUM_SD(L,ISTEP) = 0.01 * RH15M_L(L) * QS_SD_LOCAL(L)
        ENDDO

      ENDDO   ! End of timestep loop within the individual days.

!-----------------------------------------------------------------------
! Calculate daily rainfall disaggregation
!-----------------------------------------------------------------------
      DO L=1,POINTSM

!-----------------------------------------------------------------------
! Precipitation is split into four components,
! these being large scale rain, convective rain, large scale snow,
! convective snow. Call random number generator for different
! durations.
!-----------------------------------------------------------------------
        CALL RNDM(RANDOM_NUM_SD,SEED_RAIN)

!-----------------------------------------------------------------------
! Calculate type of precipitation. The decision is based purely up
! mean daily temperature, T_L. The cutoffs are:
!
! Convective scale rain (duration CONV_RAIN_DUR): T_L > 20.0oC
! Large scale rain (duration LS_RAIN_DUR) : 20.0oC > T_L > 2oC
! Large scale snow (duration LS_SNOW_DUR) : T_L < 2oC
! Convective snow - IGNORED
!
!-----------------------------------------------------------------------
! Initialise arrays
        DO ISTEP=1,STEP_DAY
          CONV_RAIN_SD(L,ISTEP) = 0.0
          LS_RAIN_SD(L,ISTEP)   = 0.0
          LS_SNOW_SD(L,ISTEP)   = 0.0
          N_EVENT(L,ISTEP)      = 0
          N_EVENT_LOCAL(ISTEP)  = 0
          PREC_LOC(ISTEP)       = 0.0
        ENDDO

! Calculate rainfall disaggregation.

! Start with convective rain
! (temperatures based upon mean daily temperature)

! First check if warm enough for convective rain
        IF(T_L(L) >= TEMP_CONV) THEN

          INIT_HOUR_CONV_RAIN = RANDOM_NUM_SD*(24.0-DUR_CONV_RAIN)
          END_HOUR_CONV_RAIN = INIT_HOUR_CONV_RAIN+DUR_CONV_RAIN

          N_TALLY = 0

          DO ISTEP=1,STEP_DAY
            HOUREVENT = (REAL(ISTEP)-0.5)*PERIOD_LEN
            IF(HOUREVENT >= INIT_HOUR_CONV_RAIN .AND.             &
               HOUREVENT < END_HOUR_CONV_RAIN) THEN
              N_EVENT(L,ISTEP) = 1
              N_TALLY = N_TALLY + 1
            ENDIF
          ENDDO

          DO ISTEP=1,STEP_DAY
            IF(N_EVENT(L,ISTEP).EQ.1) THEN      !Rains on this day
              CONV_RAIN_SD(L,ISTEP) =                             &
                        (REAL(STEP_DAY)/REAL(N_TALLY))*PRECIP_L(L)
              PREC_LOC(ISTEP) = CONV_RAIN_SD(L,ISTEP)
              N_EVENT_LOCAL(ISTEP) = N_EVENT(L,ISTEP)
            ENDIF
          ENDDO

! Check that no convective rain periods
! exceed MAX_PRECIP_RATE, or if so,
! then redistribute. The variable that is redistributed is local
! variable PREC_LOC - CONV_RAIN_SD is then set to this after the
! call to REDIS.

          CALL REDIS(                                             &
            NSDMAX,STEP_DAY,MAX_PRECIP_RATE,PREC_LOC,             &
            N_EVENT_LOCAL,N_TALLY                                 &
          )

          DO ISTEP=1,STEP_DAY
            CONV_RAIN_SD(L,ISTEP) = PREC_LOC(ISTEP)
          ENDDO

! Now look at large scale rainfall components
        ELSE IF(T_L(L) < TEMP_CONV .AND. T_L(L) >= TEMP_SNOW) THEN

          INIT_HOUR_LS_RAIN = RANDOM_NUM_SD*(24.0-DUR_LS_RAIN)
          END_HOUR_LS_RAIN = INIT_HOUR_LS_RAIN+DUR_LS_RAIN

          N_TALLY = 0

          DO ISTEP=1,STEP_DAY
            HOUREVENT = (REAL(ISTEP)-0.5)*PERIOD_LEN
            IF(HOUREVENT >= INIT_HOUR_LS_RAIN .AND.               &
               HOUREVENT < END_HOUR_LS_RAIN) THEN
              N_EVENT(L,ISTEP) = 1
              N_TALLY = N_TALLY + 1
            ENDIF
          ENDDO

          DO ISTEP=1,STEP_DAY
            IF(N_EVENT(L,ISTEP) == 1) THEN      !Rains on this day
              LS_RAIN_SD(L,ISTEP) =                               &
                        (REAL(STEP_DAY)/REAL(N_TALLY))*PRECIP_L(L)
              PREC_LOC(ISTEP) = LS_RAIN_SD(L,ISTEP)
              N_EVENT_LOCAL(ISTEP) = N_EVENT(L,ISTEP)
            ENDIF
          ENDDO

! Check that no large scale rain periods
! exceed MAX_PRECIP_RATE, or if so,
! then redistribute.

          CALL REDIS(                                             &
            NSDMAX,STEP_DAY,MAX_PRECIP_RATE,PREC_LOC,             &
            N_EVENT_LOCAL,N_TALLY                                 &
          )

          DO ISTEP=1,STEP_DAY
            LS_RAIN_SD(L,ISTEP) = PREC_LOC(ISTEP)
          ENDDO

! Now look at large scale snow components
        ELSE

          INIT_HOUR_LS_SNOW = RANDOM_NUM_SD*(24.0-DUR_LS_SNOW)
          END_HOUR_LS_SNOW = INIT_HOUR_LS_SNOW+DUR_LS_SNOW

          N_TALLY = 0

          DO ISTEP=1,STEP_DAY
            HOUREVENT = (REAL(ISTEP)-0.5)*PERIOD_LEN
            IF(HOUREVENT >= INIT_HOUR_LS_SNOW .AND.               &
               HOUREVENT < END_HOUR_LS_SNOW) THEN
              N_EVENT(L,ISTEP) = 1
              N_TALLY = N_TALLY + 1
            ENDIF
          ENDDO

          DO ISTEP=1,STEP_DAY
            IF(N_EVENT(L,ISTEP) == 1) THEN       !Rains on this da
              LS_SNOW_SD(L,ISTEP) =                               &
                        (REAL(STEP_DAY)/REAL(N_TALLY))*PRECIP_L(L)
              PREC_LOC(ISTEP) = LS_SNOW_SD(L,ISTEP)
              N_EVENT_LOCAL(ISTEP) = N_EVENT(L,ISTEP)
            ENDIF
          ENDDO

! Check that no large scale snow periods exceed MAX_PRECIP_RATE,
! or if so, then redistribute.
          CALL REDIS(                                             &
            NSDMAX,STEP_DAY,MAX_PRECIP_RATE,PREC_LOC,             &
            N_EVENT_LOCAL,N_TALLY                                 &
          )

          DO ISTEP=1,STEP_DAY
            LS_SNOW_SD(L,ISTEP) = PREC_LOC(ISTEP)
          ENDDO

        ENDIF

      ENDDO        ! End of large loop over different land points.
!                    in calculation of different rainfall behaviours

    ELSE           ! Now case where no subdaily variation (TSTEP=1)

      DO L=1,POINTSM
        PSTAR_SD_LOCAL(L) = PSTAR_SD(L,1)
        T_SD_LOCAL(L) = T_SD(L,1)
      ENDDO

      CALL QSAT(QS_SD_LOCAL,T_SD_LOCAL,PSTAR_SD_LOCAL,POINTSM)

      DO L = 1,POINTSM
        SW_SD(L,1) = SW_L(L)
        T_SD(L,1) = T_L(L)
        LW_SD(L,1) = LW_L(L)
        PSTAR_SD(L,1) = PSTAR_L(L)
        WIND_SD(L,1) = WIND_L(L)
        QHUM_SD(L,1) = 0.01*RH15M_L(L)*QS_SD_LOCAL(L)

        IF(T_L(L) >= TEMP_CONV) THEN
          CONV_RAIN_SD(L,1) = PRECIP_L(L)
        ELSE IF(T_L(L) < TEMP_CONV .AND. T_L(L) >= TEMP_SNOW) THEN
          LS_RAIN_SD(L,1) = PRECIP_L(L)
        ELSE
          LS_SNOW_SD(L,1) = PRECIP_L(L)
        ENDIF

      ENDDO

    ENDIF !End of loop to chose whether sub-daily is required

    RETURN

  END SUBROUTINE DAY_CALC
#endif
