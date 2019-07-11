#if !defined(UM_JULES)
!-----------------------------------------------------------------------
! Reads driving data for the stand-alone Interactive Vegetation Model
!                                                       (P.Cox 11/97).
! Adjusted by Chris Huntingford (08/04) for inclusion in IMOGEN
!-----------------------------------------------------------------------
  SUBROUTINE DRDAT(                                               &
    IYEAR,GPOINTS,T_ANOM_DAT,PRECIP_ANOM_DAT,RH15M_ANOM_DAT,      &
    UWIND_ANOM_DAT,VWIND_ANOM_DAT,DTEMP_ANOM_DAT,                 &
    PSTAR_HA_ANOM_DAT,SW_ANOM_DAT,LW_ANOM_DAT,DIR_ANOM,           &
    LONGMIN_DAT,LATMIN_DAT,LONGMAX_DAT,LATMAX_DAT,MM,DRIVE_MONTH  &
  )

    USE io_constants, ONLY : IMOGEN_UNIT

    IMPLICIT NONE

    INTEGER ::                                                    &
      IYEAR,                                                      &
                  ! WORK Year of interest.
      GPOINTS,                                                    &
                  ! IN Number of land points.
      MM,                                                         &
                  ! IN Number months in year.
      IM          ! WORK Month loop parameter

    CHARACTER*4 ::                                                &
      DRIVE_YEAR,                                                 &
                  ! IN Label for year of driving data.
      DRIVE_MONTH(12)
                  ! Labels for month of driving data.

    CHARACTER*180 ::                                              &
      DIR_ANOM,                                                   &
                  ! IN Directory containing anomalies.
      DRIVER_ANOM
                  ! IN Directory containing anomalies.

!-----------------------------------------------------------------------
! Climatological forcing variables.
!-----------------------------------------------------------------------
    REAL ::                                                       &
      LATMIN_DAT,LATMAX_DAT,                                      &
                 ! IN Latitudinal limits of the area
      LONGMIN_DAT,LONGMAX_DAT,                                    &
                 ! IN Longitudinal limits of the area
      LAT(GPOINTS),                                               &
                 ! WORK Latitude (degrees).
      LONG(GPOINTS),                                              &
                 ! WORK Longitude (degrees).
      DTEMP_ANOM_DAT(GPOINTS,MM),                                 &
                 ! OUT Diurnal temperature range (K).
      LW_ANOM_DAT(GPOINTS,MM),                                    &
                 ! OUT Downward surface longwave radiation (W/m).
      PSTAR_HA_ANOM_DAT(GPOINTS,MM),                              &
                 ! OUT Surface pressure (Pa).
      RAINFALL_ANOM_DAT(GPOINTS,MM),                              &
                 ! WORK Rainfall rate (mm/day).
      RH15M_ANOM_DAT(GPOINTS,MM),                                 &
                 ! OUT Relative humidity at 1.5m (%).
      SNOWFALL_ANOM_DAT(GPOINTS,MM),                              &
                 ! WORK Snowfall rate (mm/day).
      SW_ANOM_DAT(GPOINTS,MM),                                    &
                 ! OUT Downward surface shortwave radiation (W/m2).
      T_ANOM_DAT(GPOINTS,MM),                                     &
                 ! OUT Air temperature (K).
      PRECIP_ANOM_DAT(GPOINTS,MM),                                &
                 ! OUT Precipitation (mm/day)
      UWIND_ANOM_DAT(GPOINTS,MM),                                 &
                 ! OUT Windspeed u-components (m/s).
      VWIND_ANOM_DAT(GPOINTS,MM)
                 ! OUT Windspeed v-components (m/s).

!-----------------------------------------------------------------------
! Loop counters.
!-----------------------------------------------------------------------
    INTEGER :: L ! Loop counters

!-----------------------------------------------------------------------
! Convert year to string.
!-----------------------------------------------------------------------
    WRITE(DRIVE_YEAR,'(I4)') IYEAR

    DO IM=1,MM
      DRIVER_ANOM=TRIM(DIR_ANOM) // DRIVE_YEAR // DRIVE_MONTH(MM)

      OPEN(IMOGEN_UNIT, FILE=DRIVER_ANOM,                                     &
                        STATUS='old', POSITION='rewind', ACTION='read')

      READ(IMOGEN_UNIT,*) LONGMIN_DAT,LATMIN_DAT,LONGMAX_DAT,LATMAX_DAT
      DO L=1,GPOINTS
        READ(IMOGEN_UNIT,*) LONG(L),LAT(L),T_ANOM_DAT(L,IM),                  &
                            RH15M_ANOM_DAT(L,IM),UWIND_ANOM_DAT(L,IM),        &
                            VWIND_ANOM_DAT(L,IM),LW_ANOM_DAT(L,IM),           &
                            SW_ANOM_DAT(L,IM),DTEMP_ANOM_DAT(L,IM),           &
                            RAINFALL_ANOM_DAT(L,IM),                          &
                            SNOWFALL_ANOM_DAT(L,IM),PSTAR_HA_ANOM_DAT(L,IM)

        PRECIP_ANOM_DAT(L,IM) = RAINFALL_ANOM_DAT(L,IM)                       &
                              + SNOWFALL_ANOM_DAT(L,IM)
      ENDDO
    ENDDO
    CLOSE(IMOGEN_UNIT)

    RETURN

  END SUBROUTINE DRDAT
#endif
