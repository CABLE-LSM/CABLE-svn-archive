#if !defined(UM_JULES)
!-----------------------------------------------------------------
! GCM ANALOGUE Model     (C.Huntingford, P.Cox, 9/98)
! Simplified 3/04 to supply just the anomalies for change in radiative
! forcing. CO2 scenarios etc are calculated outside this routine.
!-----------------------------------------------------------------
  SUBROUTINE GCM_ANLG(                                                        &
    Q,LAND_PTS,T_ANOM_AM,PRECIP_ANOM_AM,RH15M_ANOM_AM,                        &
    UWIND_ANOM_AM,VWIND_ANOM_AM,DTEMP_ANOM_AM,PSTAR_HA_ANOM_AM,               &
    SW_ANOM_AM,LW_ANOM_AM,N_OLEVS,DIR_PATT,F_OCEAN,KAPPA_O,                   &
    LAMBDA_L,LAMBDA_O,MU,DTEMP_O,LONGMIN_AM,LATMIN_AM,LONGMAX_AM,             &
    LATMAX_AM,MM                                                              &
  )

    USE io_constants, ONLY : IMOGEN_UNIT

    USE IMOGEN_MAP, ONLY : SGINDINV

    USE imogen_constants, ONLY : DRIVE_MONTH,N_IMOGEN_LAND

    IMPLICIT NONE

    INTEGER ::                                                                &
      LAND_PTS,                                                               &
                     ! IN Number of land points.
      IM,                                                                     &
                     ! IN Loop over months
      N_OLEVS,                                                                &
                     ! IN Number of ocean thermal layers.
      MM             ! In Number of months in a year

    CHARACTER*180 ::                                                          &
      DIR_PATT,                                                               &
                     ! IN Directory containing anomaly patterns.
      DRIVER_PATT    ! WORK File containing anomaly
                     !      patterns from analogue model.

    REAL ::                                                                   &
      LAT_AM(LAND_PTS),                                                       &
                      ! Latitude read from file.
      LONG_AM(LAND_PTS)
                      !Longitude read from file.

!-----------------------------------------------------------------
! Climatological forcing variables.
!-----------------------------------------------------------------
    REAL ::                                                                   &
      F_OCEAN,                                                                &
                     ! IN Fractional coverage of the ocean.
      KAPPA_O,                                                                &
                     ! IN Ocean eddy diffusivity (W/m/K).
      LAMBDA_L,                                                               &
                     ! IN Inverse climate sensitivity over
                     !    land (W/m2/K).
      LAMBDA_O,                                                               &
                     ! IN Inverse climate sensitivity over
                     !    ocean (W/m2/K).
      MU,                                                                     &
                     ! IN Ratio of land to ocean
                     !    temperature anomalies.
      DTEMP_O(N_OLEVS),                                                       &
                     ! INOUT Ocean mean temperature anomaly (K).
      DTEMP_L,                                                                &
                     ! OUT Land mean temperature anomaly (K).
      DTEMP_ANOM_AM(LAND_PTS,MM),                                             &
                     ! OUT Diurnal temperature range (K).
      LW_ANOM_AM(LAND_PTS,MM),                                                &
                     ! OUT Downward surface longwave
                     !     radiation anomaly (W/m).
      RAINFALL_ANOM_AM(LAND_PTS,MM),                                          &
                     ! OUT Rainfall rate anomaly (mm/day)
      RH15M_ANOM_AM(LAND_PTS,MM),                                             &
                     ! OUT Relative humidity at 1.5m anom
      SNOWFALL_ANOM_AM(LAND_PTS,MM),                                          &
                     ! OUT Snowfall rate anomaly (mm/day)
      SW_ANOM_AM(LAND_PTS,MM),                                                &
                     ! OUT Downward surface shortwave
                     !     radiation anomaly (W/m2).
      T_ANOM_AM(LAND_PTS,MM),                                                 &
                     ! OUT Air temperature anomaly (K).
      PRECIP_ANOM_AM(LAND_PTS,MM),                                            &
                     ! OUT Precipitation (snowfall
                     !     plus rainfall) anomaly (mm/day)
      UWIND_ANOM_AM(LAND_PTS,MM),                                             &
                     ! OUT Wind speed anomaly (m/s).
      VWIND_ANOM_AM(LAND_PTS,MM),                                             &
                     ! OUT Wind speed anomaly (m/s).
      LATMIN_AM,LATMAX_AM,                                                    &
                     ! WORK Latitudinal limits of the area
                     !      (degrees).
      LONGMIN_AM,LONGMAX_AM,                                                  &
                     ! WORK Longitudinal limits of the area
                     !      (degrees).
      PSTAR_HA_ANOM_AM(LAND_PTS,MM),                                          &
                     ! OUT Surface pressure (hPa).
      Q
                     ! WORK Increase in radiative forcing (W/m2)

!-----------------------------------------------------------------
! Anomaly patterns scaled to land mean temperature anomalies.
!-----------------------------------------------------------------
    REAL ::                                                                   &
      DDTEMP_DAY_PAT,                                                         &
                     ! WORK Diurnal temperature range (.).
      DLW_C_PAT,                                                              &
                     ! WORK Downward surface longwave
                     !      radiation (W/m2/K).
      DPSTAR_C_PAT,                                                           &
                     ! WORK Surface pressure (hPa/K).
      DRAINFALL_PAT,                                                          &
                     ! WORK Rainfall rate (mm/day/K).
      DRH15M_PAT,                                                             &
                     ! WORK Relative humidity at 1.5m (%/K).
      DSNOWFALL_PAT,                                                          &
                     ! WORK Snowfall rate (mm/day/K).
      DSW_C_PAT,                                                              &
                     ! WORK Surface shortwave radiation (W/m2/K).
      DT_C_PAT,                                                               &
                     ! WORK Air temperature (.).
      DUWIND_PAT,DVWIND_PAT
                     ! WORK Windspeed components (m/s/K).

!-----------------------------------------------------------------
! Loop counters.
!-----------------------------------------------------------------
    INTEGER ::                                                                &
      I,L    ! WORK


!-----------------------------------------------------------------
! Loop over months
!-----------------------------------------------------------------
    DO IM=1,MM

!-----------------------------------------------------------------
! Calculate new area mean temperature anomalies
!-----------------------------------------------------------------
      IF (IM == 1) THEN
        CALL DELTA_TEMP(                                                      &
          N_OLEVS,F_OCEAN,KAPPA_O,LAMBDA_L,LAMBDA_O,MU,Q,                     &
          DTEMP_L,DTEMP_O                                                     &
        )
      ENDIF

!-----------------------------------------------------------------
! Define the anomaly patterns and read the header
!-----------------------------------------------------------------
      DRIVER_PATT = TRIM(DIR_PATT) // DRIVE_MONTH(IM)

      OPEN(IMOGEN_UNIT, FILE=DRIVER_PATT,                                     &
                        STATUS='old', POSITION='rewind', ACTION='read')

      READ(IMOGEN_UNIT,*) LONGMIN_AM,LATMIN_AM,LONGMAX_AM,LATMAX_AM

!-----------------------------------------------------------------
! Read in initial climatology and then define the new climate data.
!-----------------------------------------------------------------
      DO I=1,N_IMOGEN_LAND
        IF (SGINDINV(i) > 0) THEN
          L = SGINDINV(I)

          READ(IMOGEN_UNIT,*)                                                 &
            LONG_AM(L),LAT_AM(L),DT_C_PAT,DRH15M_PAT,DUWIND_PAT,              &
            DVWIND_PAT,DLW_C_PAT,DSW_C_PAT,DDTEMP_DAY_PAT,                    &
            DRAINFALL_PAT,DSNOWFALL_PAT,DPSTAR_C_PAT

          T_ANOM_AM(L,IM) = DT_C_PAT * DTEMP_L
          RH15M_ANOM_AM(L,IM) = DRH15M_PAT * DTEMP_L
          UWIND_ANOM_AM(L,IM) = DUWIND_PAT * DTEMP_L
          VWIND_ANOM_AM(L,IM) = DVWIND_PAT * DTEMP_L
          RAINFALL_ANOM_AM(L,IM) = DRAINFALL_PAT * DTEMP_L
          SNOWFALL_ANOM_AM(L,IM) = DSNOWFALL_PAT * DTEMP_L
          DTEMP_ANOM_AM(L,IM) = DDTEMP_DAY_PAT * DTEMP_L
          LW_ANOM_AM(L,IM) = DLW_C_PAT * DTEMP_L
          SW_ANOM_AM(L,IM) = DSW_C_PAT * DTEMP_L
          PSTAR_HA_ANOM_AM(L,IM) = DPSTAR_C_PAT * DTEMP_L

          PRECIP_ANOM_AM(L,IM) = RAINFALL_ANOM_AM(L,IM) +                     &
                                 SNOWFALL_ANOM_AM(L,IM)
        ELSE
          READ(IMOGEN_UNIT,*)
        ENDIF
      ENDDO     !End of loop over land points

      CLOSE(IMOGEN_UNIT)
    ENDDO     !End of loop over months

    RETURN
    END SUBROUTINE GCM_ANLG
#endif
