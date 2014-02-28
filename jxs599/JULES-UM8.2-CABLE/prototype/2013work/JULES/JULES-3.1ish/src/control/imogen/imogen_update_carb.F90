#if !defined(UM_JULES)
SUBROUTINE imogen_update_carb

  USE ancil_info, ONLY : land_pts
  
  USE trifctl, ONLY : cv
  
  USE prognostics, ONLY : cs

  USE imogen_run, ONLY :                                          &
    INCLUDE_CO2,C_EMISSIONS,ANOM,ANLG,CO2_INIT_PPMV,LAND_FEED,    &
    NYR_EMISS,OCEAN_FEED
    
  USE imogen_clim, ONLY :                                         &
    DCTOT,LAT
    
  USE imogen_constants, ONLY :                                    &
    CONV,OCEAN_AREA,NFARRAY
    
  USE imogen_time, ONLY :                                         &
    IYEAR,YEAR1,IYEND
    
  USE imogen_progs, ONLY :                                        &
    CO2_PPMV,CO2_START_PPMV,DTEMP_O,CO2_CHANGE_PPMV,FA_OCEAN
    
  USE imogen_anlg_vals, ONLY :                                    &
    T_OCEAN_INIT
    
  USE imogen_io_vars, ONLY :                                      &
    YR_EMISS,C_EMISS

  IMPLICIT NONE
  
  INTEGER ::                                                      &
    EMISS_TALLY ! Checks that datafile of emissions includes
                ! year of interest
                
  INTEGER :: i,n ! loop counters

  REAL ::                                                         &
    C_EMISS_LOCAL  ! Local value of C_EMISS
    
  REAL ::                                                         &
    D_LAND_ATMOS,                                                 &
                ! Change in atmospheric CO2 concentration due to
                ! feedbacks (ppm/year).
    D_OCEAN_ATMOS
                ! Change in atmospheric CO2 concentration due to
                ! feedbacks (ppm/year)
      
!-----------------------------------------------------------------


!-----------------------------------------------------------------
! Calculate land atmosphere carbon exchange.
!-----------------------------------------------------------------
  IF(LAND_FEED) THEN 
    DO I=1,LAND_PTS
      DCTOT(I) = CV(I) + CS(I,1) - DCTOT(I)
    ENDDO

    CALL DIFFCARB_LAND(LAND_PTS,D_LAND_ATMOS,LAT,DCTOT,CONV)
  ENDIF

!-----------------------------------------------------------------
! Now do the carbon cycling update.
!-----------------------------------------------------------------
  IF(INCLUDE_CO2 .AND. C_EMISSIONS .AND. ANOM .AND. ANLG) THEN

! Include anthropogenic carbon emissions.
    EMISS_TALLY=0
    DO N = 1,NYR_EMISS
      IF(YR_EMISS(N) == IYEAR) THEN
        C_EMISS_LOCAL = C_EMISS(N)
        CO2_PPMV = CO2_PPMV + CONV * C_EMISS_LOCAL 
        EMISS_TALLY = EMISS_TALLY + 1
! We have found the right year so we can exit the loop
        EXIT
      ENDIF
    ENDDO

    IF(EMISS_TALLY /= 1) THEN
      WRITE(*,*) 'IMOGEN: Emission dataset does not match run'
      STOP
    ENDIF

    IF(LAND_FEED) THEN 
! Update with land feedbacks if required
      CO2_PPMV = CO2_PPMV + D_LAND_ATMOS
    ENDIF

    IF(OCEAN_FEED) THEN 
! Update with ocean feedbacks if required
      CALL OCEAN_CO2(                                             &
        IYEAR-YEAR1+1,1,CO2_PPMV,CO2_INIT_PPMV,DTEMP_O(1),        &
        FA_OCEAN,OCEAN_AREA,CO2_CHANGE_PPMV,IYEND-YEAR1,          &
        T_OCEAN_INIT,NFARRAY,D_OCEAN_ATMOS                        &
      )
      CO2_PPMV = CO2_PPMV + D_OCEAN_ATMOS    
    ENDIF
  ENDIF

  IF(INCLUDE_CO2) CO2_CHANGE_PPMV = CO2_PPMV - CO2_START_PPMV


  RETURN

END SUBROUTINE imogen_update_carb
#endif
