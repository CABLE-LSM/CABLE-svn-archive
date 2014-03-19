#if !defined(UM_JULES)
SUBROUTINE imogen_update_clim

  USE io_constants, ONLY : IMOGEN_UNIT

  USE datetime_mod, ONLY : SECS_IN_DAY

  USE aero, ONLY : co2_mmr

  USE ancil_info, ONLY : land_pts

  USE trifctl, ONLY : cv

  USE prognostics, ONLY : cs

  USE imogen_progs, ONLY :                                        &
    CO2_PPMV,CO2_START_PPMV,DTEMP_O,SEED_RAIN

  USE imogen_run, ONLY :                                          &
    INCLUDE_CO2,C_EMISSIONS,FILE_SCEN_CO2_PPMV,ANOM,ANLG,         &
    CO2_INIT_PPMV,INCLUDE_NON_CO2,FILE_NON_CO2_VALS,WGEN,         &
    LAND_FEED

  USE imogen_time, ONLY :                                         &
    IYEAR,MM,MD,NSDMAX,STEP_DAY

  USE imogen_anlg_vals, ONLY :                                    &
    Q2CO2,NYR_NON_CO2,FILE_NON_CO2,DIR_PATT,F_OCEAN,KAPPA_O,      &
    LAMBDA_L,LAMBDA_O,MU,DIR_ANOM

  USE imogen_constants, ONLY :                                    &
    N_OLEVS,DRIVE_MONTH,MDI

  USE imogen_clim, ONLY :                                         &
!   Scalars
    LATMIN_CLIM,LATMAX_CLIM,LONGMIN_CLIM,LONGMAX_CLIM,            &
!   Arrays
    T_CLIM,RAINFALL_CLIM,SNOWFALL_CLIM,RH15M_CLIM,UWIND_CLIM,     &
    VWIND_CLIM,DTEMP_CLIM,PSTAR_HA_CLIM,SW_CLIM,LW_CLIM,          &
    F_WET_CLIM,LAT,LONG,DCTOT

  USE imogen_drive_vars, ONLY :                                   &
    T_OUT,CONV_RAIN_OUT,CONV_SNOW_OUT,LS_RAIN_OUT,LS_SNOW_OUT,    &
    QHUM_OUT,WIND_OUT,PSTAR_OUT,SW_OUT,LW_OUT

  USE logging_mod, ONLY : log_fatal

  USE jules_print_mgr, ONLY :                                       &
    jules_message,                                                 &
    jules_print
  IMPLICIT NONE

  INTEGER :: i ! Loop counter

  INTEGER ::                                                      &
    TALLY_CO2_FILE,                                               &
                   !WORK If used, checks CO2 value available in
                   !     file "FILE_SCEN_CO2_PPMV"
    YR_CO2_FILE,                                                  &
                   !WORK If used, reads in years available in
                   !     file "FILE_SCEN_CO2_PPMV"
    kode           ! Used to hold IOSTAT while reading file

  REAL ::                                                         &
    CO2_FILE_PPMV  ! If used, is prescribed CO2 value for read
                   ! in years in "FILE_SCEN_CO2_PPMV"

  REAL ::                                                         &
    Q_CO2,                                                        &
                 !WORK Radiative forcing due to CO2 (W/m2)
    Q_NON_CO2,                                                    &
                 !WORK Radiative forcing due to non-CO2 (W/m2)
    Q            !WORK Total radiative forcing, both CO2 and NON CO2

  REAL :: LATMIN,LATMAX,LONGMIN,LONGMAX
                 !WORK Max and min lat and long for files

! Dummy variables for weather generator which is not available in
! this version
  REAL ::                                                         &
    PRECIP_WG(LAND_PTS,MM,MD),                                    &
    TMIN_WG(LAND_PTS,MM,MD),                                      &
    TMAX_WG(LAND_PTS,MM,MD),                                      &
    SW_WG(LAND_PTS,MM,MD),                                        &
    RH15M_WG(LAND_PTS,MM,MD)

!-----------------------------------------------------------------
! Variables to hold calculated anomalies
!-----------------------------------------------------------------
  REAL ::                                                         &
    T_ANOM(LAND_PTS,MM),                                          &
                  ! WORK Temperature anomalies (K)
    PRECIP_ANOM(LAND_PTS,MM),                                     &
                  ! WORK Precip anomalies (mm/day)
    RH15M_ANOM(LAND_PTS,MM),                                      &
                  ! WORK Relative humidity anomalies
    UWIND_ANOM(LAND_PTS,MM),                                      &
                  ! WORK u-wind anomalies (m/s)
    VWIND_ANOM(LAND_PTS,MM),                                      &
                  ! WORK v-wind anomalies (m/s)
    DTEMP_ANOM(LAND_PTS,MM),                                      &
                  ! WORK Diurnal Temperature (K)
    PSTAR_HA_ANOM(LAND_PTS,MM),                                   &
                  ! WORK Pressure anomalies (hPa)
    SW_ANOM(LAND_PTS,MM),                                         &
                  ! WORK Shortwave radiation anomalie
    LW_ANOM(LAND_PTS,MM)
                  ! WORK Longwave radiation anomalies

!------------------------------------------------------------------------
! Initialisation
  Q_CO2     = 0.0
  Q_NON_CO2 = 0.0

  T_ANOM(:,:)        = 0.0
  PRECIP_ANOM(:,:)   = 0.0
  RH15M_ANOM(:,:)    = 0.0
  UWIND_ANOM(:,:)    = 0.0
  VWIND_ANOM(:,:)    = 0.0
  DTEMP_ANOM(:,:)    = 0.0
  PSTAR_HA_ANOM(:,:) = 0.0
  SW_ANOM(:,:)       = 0.0
  LW_ANOM(:,:)       = 0.0

  WRITE(jules_message,*) 'Updating IMOGEN climate'
  CALL jules_print('imogen_update_clim',jules_message)


! Capture CO2 concentration at beginning of year
  IF(INCLUDE_CO2) THEN
    CO2_START_PPMV = CO2_PPMV
  ENDIF


! Hydrology 20th Century simulations (note also check for this run in
! subroutine IMOGEN_CONFIRMED_RUN which includes more stringent checks)
! OR analogue model simulations with CO2 prescribed.
  IF((.NOT.C_EMISSIONS) .AND. INCLUDE_CO2) THEN
! This works by reading in a file of CO2 concentrations, and checks that
! year is represented
    OPEN(IMOGEN_UNIT, FILE=FILE_SCEN_CO2_PPMV,                                &
                        STATUS='old', POSITION='rewind', ACTION='read')

    TALLY_CO2_FILE = 0
    kode = 0

    DO WHILE (.TRUE.)
      READ(IMOGEN_UNIT,FMT=*,IOSTAT=kode) YR_CO2_FILE,CO2_FILE_PPMV
! Check for end of file (or just any error) and exit loop if found
! iostat > 0 means illegal data
! iostat < 0 means end of record or end of file
      IF(kode /= 0) EXIT

      IF(YR_CO2_FILE == IYEAR) THEN
        CO2_PPMV = CO2_FILE_PPMV
        TALLY_CO2_FILE = TALLY_CO2_FILE + 1
! We have found the correct year so exit loop
        EXIT
      ENDIF
    ENDDO

    CLOSE(IMOGEN_UNIT)

! Check that value has been found.
    IF(TALLY_CO2_FILE /= 1)                                                   &
      CALL log_fatal("imogen_update_clim", 'CO2 value not found in file')
  ENDIF

! Now calculate the added monthly anomalies, either from analogue model
! or prescribed directly.
  IF(ANOM) THEN
    IF(ANLG) THEN
! This call is to the GCM analogue model. It is prescribed CO2
! concentration, and calculates non-CO2, and returns total change in
! radiative forcing, Q. Recall that the AM has a "memory" through
! DTEMP_0 - ie the ocean temperatures. Note that in this version of the
! code, the AM is updated yearly.

! Calculate the CO2 forcing
      IF(INCLUDE_CO2) THEN
        CALL RADF_CO2(CO2_PPMV,CO2_INIT_PPMV,Q2CO2,Q_CO2)
      ENDIF

! Calculate the non CO2 forcing
      IF(INCLUDE_NON_CO2) THEN
        CALL RADF_NON_CO2(                                        &
          IYEAR,Q_NON_CO2,NYR_NON_CO2,FILE_NON_CO2,               &
          FILE_NON_CO2_VALS                                       &
        )
      ENDIF

! Calculate the total forcing
      Q = Q_CO2 + Q_NON_CO2

! Call the GCM analogue model that responds to this forcing
      IF(INCLUDE_CO2 .OR. INCLUDE_NON_CO2) THEN
        CALL GCM_ANLG(                                            &
          Q,LAND_PTS,T_ANOM,PRECIP_ANOM,RH15M_ANOM,UWIND_ANOM,    &
          VWIND_ANOM,DTEMP_ANOM,PSTAR_HA_ANOM,SW_ANOM,LW_ANOM,    &
          N_OLEVS,DIR_PATT,F_OCEAN,KAPPA_O,LAMBDA_L,LAMBDA_O,     &
          MU,DTEMP_O,LONGMIN,LATMIN,LONGMAX,LATMAX,MM             &
        )

! Check driving files are compatible.
        IF((ABS(LONGMIN_CLIM - LONGMIN) >= 1.0E-6) .OR.           &
           (ABS(LATMIN_CLIM - LATMIN) >= 1.0E-6) .OR.             &
           (ABS(LONGMAX_CLIM - LONGMAX) >= 1.0E-6) .OR.           &
           (ABS(LATMAX_CLIM - LATMAX) >= 1.0E-6)) THEN
          CALL log_fatal("imogen_update_clim",                    &
                         'Driving files are incompatible')
        ENDIF
      ENDIF

! Option where anomalies are prescribed directly not using the analogue
    ELSE
      CALL DRDAT(                                                 &
        IYEAR,LAND_PTS,T_ANOM,PRECIP_ANOM,RH15M_ANOM,UWIND_ANOM,  &
        VWIND_ANOM,DTEMP_ANOM,PSTAR_HA_ANOM,SW_ANOM,LW_ANOM,      &
        DIR_ANOM,LONGMIN,LATMIN,LONGMAX,LATMAX,MM,DRIVE_MONTH     &
      )

! Check driving files are compatible.
      IF((ABS(LONGMIN_CLIM - LONGMIN) >= 1.0E-6) .OR.             &
         (ABS(LATMIN_CLIM - LATMIN) >= 1.0E-6) .OR.               &
         (ABS(LONGMAX_CLIM - LONGMAX) >= 1.0E-6) .OR.             &
         (ABS(LATMAX_CLIM - LATMAX) >= 1.0E-6)) THEN
        CALL log_fatal("imogen_update_clim",                      &
                       'Driving files are incompatible')
      ENDIF
    ENDIF
  ELSE
! Set anomalies to zero.
    T_ANOM(:,:) = 0.0
    SW_ANOM(:,:) = 0.0
    LW_ANOM(:,:) = 0.0
    PSTAR_HA_ANOM(:,:) = 0.0
    RH15M_ANOM(:,:) = 0.0
    PRECIP_ANOM(:,:) = 0.0
    UWIND_ANOM(:,:) = 0.0
    VWIND_ANOM(:,:) = 0.0
    DTEMP_ANOM(:,:) = 0.0
  ENDIF         !End of where anomalies are calculated

! Now incorporate anomalies with climate data.
! At this point, have climatology, WG if called and anomalies of either
! "_AM" or "_DAT". Now calculate the daily values of the driving data.
! This is calculated using subroutine CLIM_CALC
  CALL CLIM_CALC(                                                 &
    ANOM,ANLG,LAND_PTS,WGEN,MM,MD,T_CLIM,SW_CLIM,LW_CLIM,         &
    PSTAR_HA_CLIM,RH15M_CLIM,RAINFALL_CLIM,SNOWFALL_CLIM,         &
    UWIND_CLIM,VWIND_CLIM,DTEMP_CLIM,F_WET_CLIM,TMIN_WG,          &
    TMAX_WG,SW_WG,RH15M_WG,PRECIP_WG,T_ANOM,SW_ANOM,LW_ANOM,      &
    PSTAR_HA_ANOM,RH15M_ANOM,PRECIP_ANOM,UWIND_ANOM,              &
    VWIND_ANOM,DTEMP_ANOM,SW_OUT,T_OUT,LW_OUT,CONV_RAIN_OUT,      &
    CONV_SNOW_OUT,LS_RAIN_OUT,LS_SNOW_OUT,PSTAR_OUT,WIND_OUT,     &
    QHUM_OUT,NSDMAX,STEP_DAY,SEED_RAIN,SECS_IN_DAY,LAT,LONG,MDI   &
  )

! Compute current land carbon.
  IF(LAND_FEED) THEN
    DO I=1,LAND_PTS
      DCTOT(I) = CV(I) + CS(I,1)
    ENDDO
  ENDIF

  co2_mmr = CO2_PPMV * 44.0 / 28.97 * 1.0e-6


  RETURN

END SUBROUTINE imogen_update_clim
#endif
