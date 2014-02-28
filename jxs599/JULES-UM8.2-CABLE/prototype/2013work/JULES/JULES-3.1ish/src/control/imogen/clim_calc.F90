#if !defined(UM_JULES)
  SUBROUTINE CLIM_CALC(                                           &
    ANOM,ANLG,GPOINTS,WGEN,MM,MD,T_CLIM,SW_CLIM,LW_CLIM,          &
    PSTAR_HA_CLIM,RH15M_CLIM,RAINFALL_CLIM,SNOWFALL_CLIM,         &
    UWIND_CLIM,VWIND_CLIM,DTEMP_CLIM,F_WET_CLIM,TMIN_WG,TMAX_WG,  &
    SW_WG,RH15M_WG,PRECIP_WG,T_ANOM,SW_ANOM,LW_ANOM,              &
    PSTAR_HA_ANOM,RH15M_ANOM,PRECIP_ANOM,UWIND_ANOM,VWIND_ANOM,   &
    DTEMP_ANOM,SW_OUT,T_OUT,LW_OUT,CONV_RAIN_OUT,CONV_SNOW_OUT,   &
    LS_RAIN_OUT,LS_SNOW_OUT,PSTAR_OUT,WIND_OUT,QHUM_OUT,NSDMAX,   &
    STEP_DAY,SEED_RAIN,SEC_DAY,LAT,LONG,MDI                       &
  )

    IMPLICIT NONE                                                     
                                                                        
    INTEGER ::                                                    &
      MM,                                                         &
               !IN Number of months
      MD,                                                         &
               !IN Number of days in (GCM) month
      STEP_DAY,                                                   &
               !IN Number of daily timesteps of IMPACTS_MODEL
      ISTEP,                                                      &
               !Looping parameter over suib-daily periods
      SEC_DAY,                                                    &
               !IN Number of seconds in each day
      NSDMAX,                                                     &
               !IN Maximum number of possible subdaily increments
      L,J,K,                                                      &
               !Loop parameters
      GPOINTS,                                                    &
               !Number of points, not including Antartica
      SEED_RAIN(4)
               !WORK Seeding number for subdaily rainfall.
                                                                        
    LOGICAL ::                                                    &
      ANLG,                                                       &
               !IN If true, then use the GCM analogue model
      ANOM,                                                       &
               !IN If true, then use the GCM analogue model
      WGEN     !IN Is the weather generator switched on.

    REAL ::                                                       &
      T_ANOM(GPOINTS,MM),                                         &
               ! Temperature anomalies fro
      PRECIP_ANOM(GPOINTS,MM),                                    &
               ! Precip anomalies from AM 
      RH15M_ANOM(GPOINTS,MM),                                     &
               ! Relative humidity anomali
      UWIND_ANOM(GPOINTS,MM),                                     &
               ! u-wind anomalies from AM 
      VWIND_ANOM(GPOINTS,MM),                                     &
               ! v-wind anomalies from AM 
      DTEMP_ANOM(GPOINTS,MM),                                     &
               ! Diurnal Temperature (K)  
      PSTAR_HA_ANOM(GPOINTS,MM),                                  &
               ! Pressure anomalies from A
      SW_ANOM(GPOINTS,MM),                                        &
               ! Shortwave radiation anoma
      LW_ANOM(GPOINTS,MM),                                        &
               ! Longwave radiation anomal
      LAT(GPOINTS),                                               &
               ! Latitudinal position of l
      LONG(GPOINTS)
               ! Longitudinal position of 

! Driving "control" climatology
    REAL ::                                                       &
      T_CLIM(GPOINTS,MM),                                         &
               !IN Control climate temperature
      RAINFALL_CLIM(GPOINTS,MM),                                  &
               !IN Control climate rainfall
      SNOWFALL_CLIM(GPOINTS,MM),                                  &
               !IN Control climate snowfall
      RH15M_CLIM(GPOINTS,MM),                                     &
               !IN Control climate relative humidity
      UWIND_CLIM(GPOINTS,MM),                                     &
               !IN Control climate u-wind
      VWIND_CLIM(GPOINTS,MM),                                     &
               !IN Control climate v-wind
      DTEMP_CLIM(GPOINTS,MM),                                     &
               !IN Control climate diurnal Tem
      PSTAR_HA_CLIM(GPOINTS,MM),                                  &
               !IN Control climate pressure
      SW_CLIM(GPOINTS,MM),                                        &
               !IN Control climate shortwave radiation
      LW_CLIM(GPOINTS,MM),                                        &
               !IN Control climate longwave radiation
      F_WET_CLIM(GPOINTS,MM)
               !IN Control climate fraction we

! Output from the weather generator when called                         
    REAL ::                                                       &
      PRECIP_WG(GPOINTS,MM,MD),                                   &
               ! Daily precipitation
      TMIN_WG(GPOINTS,MM,MD),                                     &
               ! Daily minimum temperature
      TMAX_WG(GPOINTS,MM,MD),                                     &
               ! Daily maximum temperature
      SW_WG(GPOINTS,MM,MD),                                       &
               ! Daily shortwave radiation
      RH15M_WG(GPOINTS,MM,MD)
               ! Daily relative humidity

! Create "local" (in time) values of the arrays below                   
    REAL ::                                                       &
      T_OUT_LOCAL(GPOINTS,NSDMAX),                                &
               !WORK temperature (K)
      CONV_RAIN_OUT_LOCAL(GPOINTS,NSDMAX),                        &
               !WORK temperature (mm/day)
      LS_RAIN_OUT_LOCAL(GPOINTS,NSDMAX),                          &
               !WORK temperature (mm/day)
      LS_SNOW_OUT_LOCAL(GPOINTS,NSDMAX),                          &
               !WORK temperature (mm/day)
      QHUM_OUT_LOCAL(GPOINTS,NSDMAX),                             &
               !WORK humidity (kg/kg)
      WIND_OUT_LOCAL(GPOINTS,NSDMAX),                             &
               !WORK wind  (m/s)
      PSTAR_OUT_LOCAL(GPOINTS,NSDMAX),                            &
               !WORK pressure (Pa)
      SW_OUT_LOCAL(GPOINTS,NSDMAX),                               &
               !WORK shortwave radiation
      LW_OUT_LOCAL(GPOINTS,NSDMAX)
               !WORK longwave radiation

! Create fine temperal resolution year of climatology (to be used by imp
! studies or DGVMs).                                                    
    REAL ::                                                       &
      T_OUT(GPOINTS,MM,MD,NSDMAX),                                &
               !OUT Calculated temperature 
      CONV_RAIN_OUT(GPOINTS,MM,MD,NSDMAX),                        &
               !OUT Calculated convective rainfall (mm/day)
      CONV_SNOW_OUT(GPOINTS,MM,MD,NSDMAX),                        &
               !OUT Calculated convective rainfall (mm/day)
      LS_RAIN_OUT(GPOINTS,MM,MD,NSDMAX),                          &
               !OUT Calculated large scale rainfall (mm/day)
      LS_SNOW_OUT(GPOINTS,MM,MD,NSDMAX),                          &
               !OUT Calculated large scale snowfall (mm/day)
      QHUM_OUT(GPOINTS,MM,MD,NSDMAX),                             &
               !OUT Calculated humidity
      WIND_OUT(GPOINTS,MM,MD,NSDMAX),                             &
               !OUT Calculated wind  (m/s)
      PSTAR_OUT(GPOINTS,MM,MD,NSDMAX),                            &
               !OUT Calculated pressure
      SW_OUT(GPOINTS,MM,MD,NSDMAX),                               &
               !OUT Calculated shortwave radiation
      LW_OUT(GPOINTS,MM,MD,NSDMAX)
               !OUT Calculated longwave radiation

! Variables for daily climatology                                       
    REAL ::                                                       &
      T_DAILY(GPOINTS,MM,MD),                                     &
               !WORK Calculated temperature
      PRECIP_DAILY(GPOINTS,MM,MD),                                &
               !WORK Calculated temperature
      UWIND_DAILY(GPOINTS,MM,MD),                                 &
               !WORK Calculated "u"-wind
      VWIND_DAILY(GPOINTS,MM,MD),                                 &
               !WORK Calculated "v"-wind
      WIND_DAILY(GPOINTS,MM,MD),                                  &
               !WORK Calculated wind  (m/s)
      DTEMP_DAILY(GPOINTS,MM,MD),                                 &
               !WORK Calculated diurnal Temperature
      PSTAR_DAILY(GPOINTS,MM,MD),                                 &
               !WORK Calculated pressure (Pa) 
      SW_DAILY(GPOINTS,MM,MD),                                    &
               !WORK Calculated shortwave radiation
      LW_DAILY(GPOINTS,MM,MD)
               !WORK Calculated longwave radiation

! And "local" (in time) values of the DAILY variables.                  
    REAL ::                                                       &
      T_DAILY_LOCAL(GPOINTS),                                     &
               !WORK Calculated temperature (K)
      PRECIP_DAILY_LOCAL(GPOINTS),                                &
               !WORK Calculated precip
      RH15M_DAILY_LOCAL(GPOINTS),                                 &
               !WORK Calculated humidity (%)
      WIND_DAILY_LOCAL(GPOINTS),                                  &
               !WORK Calculated wind  (m/s)
      DTEMP_DAILY_LOCAL(GPOINTS),                                 &
               !WORK Calculated diurnal Temperature
      PSTAR_DAILY_LOCAL(GPOINTS),                                 &
               !WORK Calculated pressure (Pa) 
      SW_DAILY_LOCAL(GPOINTS),                                    &
               !WORK Calculated shortwave radiation
      LW_DAILY_LOCAL(GPOINTS)
               !WORK Calculated longwave radiation

    REAL ::                                                       &
      RH15M_DAILY(GPOINTS,MM,MD),                                 &
               !WORK Calculated relative humidity
      MDI      !WORK Missing data indicator

    REAL ::                                                       &
      TLOCAL,                                                     &
               !WORK Local value of T_DAILY in
      PLOCAL   !WORK Local value of PSTAR_DAIL
                                                                        
!---------------------------------------------------------------------
! Variables required to split rainfall up so that it rains roughly
! the correct no. of days/month when weather generator is switched off.
!    INTEGER ::                                                    &
!      NO_RAINDAY,                                                 &
               ! WORK No. of rainy days in the month.
!      INT_RAINDAY
               ! WORK Rain day interval.
!      IC_RAINDAY
               ! WORK No. of rain days counter.
!    REAL ::                                                       &
!      TOT_RAIN ! WORK Total rain counter.
!---------------------------------------------------------------------


! Calculate monthly means and add anomalies.
! Anomalies will be zero, if set ANOM=.F.
    DO J=1,MM                 !Loop over the months                  
      DO L=1,GPOINTS         !Loop over land points
        DO K=1,MD           !Loop over the days                    

          IF(WGEN) THEN                                            
!     Temperature (K)                                                   
            T_DAILY(L,J,K) = (0.5 *                               &
                              (TMIN_WG(L,J,K) + TMAX_WG(L,J,K))   &
                             ) + T_ANOM(L,J)

!     Shortwave radiation (W/m2)                                        
            SW_DAILY(L,J,K) = SW_WG(L,J,K) + SW_ANOM(L,J)

!     Relative humidity (kg/kg)
            RH15M_DAILY(L,J,K) = RH15M_WG(L,J,K) + RH15M_ANOM(L,J)
            DTEMP_DAILY(L,J,K) = (TMAX_WG(L,J,K)-TMIN_WG(L,J,K))  &
                               + DTEMP_ANOM(L,J)

!     Precip : At present, the Hadley Centre GCM AM has only PRECIP in t
!     column (ie include snowfall).                                     
            PRECIP_DAILY(L,J,K) = PRECIP_WG(L,J,K)+PRECIP_ANOM(L,J)
          ELSE
            T_DAILY(L,J,K)     = T_CLIM(L,J) + T_ANOM(L,J)
            SW_DAILY(L,J,K)    = SW_CLIM(L,J) + SW_ANOM(L,J)
            RH15M_DAILY(L,J,K) = RH15M_CLIM(L,J) + RH15M_ANOM(L,J)
            DTEMP_DAILY(L,J,K) = DTEMP_CLIM(L,J) + DTEMP_ANOM(L,J)

            PRECIP_DAILY(L,J,K) = RAINFALL_CLIM(L,J)              &
                                + SNOWFALL_CLIM(L,J)              &
                                + PRECIP_ANOM(L,J)

!CH-EDITTED OUT LINES BELOW BUT WILL BE RE_IMPLEMENTED WITH DOCUMENTATIO
! To correct the no. of rain days per month if WGEN is off:
!                 IF(K.EQ.1)THEN
!                    NO_RAINDAY=NINT(F_WET_CLIM(L,J)*MD)
!                    IF(NO_RAINDAY.LT.1)NO_RAINDAY=1
!                    IF(NO_RAINDAY.GT.MD)NO_RAINDAY=MD
!                    INT_RAINDAY=INT(FLOAT(MD)/FLOAT(NO_RAINDAY))
!                    IF(INT_RAINDAY.GT.MD)INT_RAINDAY=MD
!                    IC_RAINDAY=0
!                    TOT_RAIN=0.0
!                 ENDIF
!
!                 IF(MOD(K,INT_RAINDAY).EQ.0
!    &                 .AND.IC_RAINDAY.LT.NO_RAINDAY)THEN
!                    PRECIP_DAILY(L,J,K) = (RAINFALL_CLIM(L,J)+         
!    &                 SNOWFALL_CLIM(L,J)+PRECIP_ANOM(L,J))
!    &                    *FLOAT(MD)/FLOAT(NO_RAINDAY)
!                    IC_RAINDAY=IC_RAINDAY+1
!                 ELSE
!                    PRECIP_DAILY(L,J,K) = 0.0
!                 ENDIF
!                 TOT_RAIN=TOT_RAIN+PRECIP_DAILY(L,J,K)
!CH-END OF EDITTING OUT
          ENDIF

!     Make sure precip anomalies do not produce negative rainfall   
          PRECIP_DAILY(L,J,K) = MAX(PRECIP_DAILY(L,J,K), 0.0)

!     Pressure (Pa) - Include unit conversion from HPa to Pa
          PSTAR_DAILY(L,J,K) =                                    &
                    100.0*(PSTAR_HA_CLIM(L,J)+PSTAR_HA_ANOM(L,J))

!     Check on humidity bounds
          RH15M_DAILY(L,J,K) = MIN(RH15M_DAILY(L,J,K), 100.0)
          RH15M_DAILY(L,J,K) = MAX(RH15M_DAILY(L,J,K), 0.0)
!     Now convert humidity units into required (g/kg)
          TLOCAL = T_DAILY(L,J,K)
          PLOCAL = PSTAR_DAILY(L,J,K)

!     Check to make sure anomalies do not produce negative values
          DTEMP_DAILY(L,J,K) = MAX(DTEMP_DAILY(L,J,K), 0.0)

!     Longwave radiation (W/m2)                                         
          LW_DAILY(L,J,K) = LW_CLIM(L,J) + LW_ANOM(L,J)

!     Wind
          UWIND_DAILY(L,J,K) = UWIND_CLIM(L,J) + UWIND_ANOM(L,J)
          VWIND_DAILY(L,J,K) = VWIND_CLIM(L,J) + VWIND_ANOM(L,J)
          WIND_DAILY(L,J,K)  = SQRT(                              &
            (UWIND_DAILY(L,J,K)**2) + (VWIND_DAILY(L,J,K)**2)     &
          )

!     Check to make sure anomalies do not produce zero windspeed (ie  
!     below measurement level).                                         
          WIND_DAILY(L,J,K) = MAX(WIND_DAILY(L,J,K), 0.01)
        ENDDO
      ENDDO
    ENDDO

!     Disaggregate down to sub-daily (ie TSTEP values)                  
!     This calls subroutine DAY_CALC which for each day converts values 
!     "_daily" to "_out".                                               

!     Variables going in (with units) are SW(W/m2), Precip (mm/day), Tem
!     (K), DTemp (K), LW (W/m2), PSTAR(Pa), Wind (m/2) and QHUM (kg/kg) 
      
!     Variables coming out are subdaily estimates of above variables, ex
!     for DTEMP (which no longer has meaning), and temperature dependent
!     the splitting of precipitation back into LS_CONV, LS_SNOW and     
!     Convective rainfall (there is no convective snow, and so below, th
!     set of have a zero value).                                        
                                                                        
    DO J=1,MM                !Loop over the months                    
      DO K=1,MD              !Loop over the days                      
        DO L = 1,GPOINTS
          SW_DAILY_LOCAL(L)     = SW_DAILY(L,J,K)
          PRECIP_DAILY_LOCAL(L) = PRECIP_DAILY(L,J,K)
          T_DAILY_LOCAL(L)      = T_DAILY(L,J,K)
          DTEMP_DAILY_LOCAL(L)  = DTEMP_DAILY(L,J,K)
          LW_DAILY_LOCAL(L)     = LW_DAILY(L,J,K)
          PSTAR_DAILY_LOCAL(L)  = PSTAR_DAILY(L,J,K)
          WIND_DAILY_LOCAL(L)   = WIND_DAILY(L,J,K)
          RH15M_DAILY_LOCAL(L)  = RH15M_DAILY(L,J,K)
        ENDDO                                                         
                                                                        
        CALL DAY_CALC(                                            &
          GPOINTS,STEP_DAY,MD,SW_DAILY_LOCAL,PRECIP_DAILY_LOCAL,  &
          T_DAILY_LOCAL,DTEMP_DAILY_LOCAL,LW_DAILY_LOCAL,         &
          PSTAR_DAILY_LOCAL,WIND_DAILY_LOCAL,RH15M_DAILY_LOCAL,   &
          SW_OUT_LOCAL,T_OUT_LOCAL,LW_OUT_LOCAL,                  &
          CONV_RAIN_OUT_LOCAL,LS_RAIN_OUT_LOCAL,                  &
          LS_SNOW_OUT_LOCAL,PSTAR_OUT_LOCAL,WIND_OUT_LOCAL,       &
          QHUM_OUT_LOCAL,J,K,LAT,LONG,SEC_DAY,NSDMAX,SEED_RAIN    &
        )                       

! Finalise value and set unused output values as MDI as a precaution.   
        DO L = 1,GPOINTS                                              
          DO ISTEP = 1,STEP_DAY
            SW_OUT(L,J,K,ISTEP)        = SW_OUT_LOCAL(L,ISTEP)
            T_OUT(L,J,K,ISTEP)         = T_OUT_LOCAL(L,ISTEP)
            LW_OUT(L,J,K,ISTEP)        = LW_OUT_LOCAL(L,ISTEP)
            CONV_RAIN_OUT(L,J,K,ISTEP) =                          &
                                     CONV_RAIN_OUT_LOCAL(L,ISTEP)
            CONV_SNOW_OUT(L,J,K,ISTEP) = 0.0
            LS_RAIN_OUT(L,J,K,ISTEP)   = LS_RAIN_OUT_LOCAL(L,ISTEP)
            LS_SNOW_OUT(L,J,K,ISTEP)   = LS_SNOW_OUT_LOCAL(L,ISTEP)
            PSTAR_OUT(L,J,K,ISTEP)     = PSTAR_OUT_LOCAL(L,ISTEP)
            WIND_OUT(L,J,K,ISTEP)      = WIND_OUT_LOCAL(L,ISTEP)
            QHUM_OUT(L,J,K,ISTEP)      = QHUM_OUT_LOCAL(L,ISTEP)
          ENDDO                                                       
                                                                        
          DO ISTEP = STEP_DAY+1,NSDMAX                                
            SW_OUT(L,J,K,ISTEP)        = MDI
            T_OUT(L,J,K,ISTEP)         = MDI
            LW_OUT(L,J,K,ISTEP)        = MDI
            CONV_RAIN_OUT(L,J,K,ISTEP) = MDI
            CONV_SNOW_OUT(L,J,K,ISTEP) = MDI
            LS_RAIN_OUT(L,J,K,ISTEP)   = MDI
            LS_SNOW_OUT(L,J,K,ISTEP)   = MDI
            PSTAR_OUT(L,J,K,ISTEP)     = MDI
            WIND_OUT(L,J,K,ISTEP)      = MDI
            QHUM_OUT(L,J,K,ISTEP)      = MDI
          ENDDO
        ENDDO
      ENDDO                  !End of loop over days                   
    ENDDO                    !End of loop over months                 
                                                                        
    RETURN
  
  END SUBROUTINE CLIM_CALC
#endif
