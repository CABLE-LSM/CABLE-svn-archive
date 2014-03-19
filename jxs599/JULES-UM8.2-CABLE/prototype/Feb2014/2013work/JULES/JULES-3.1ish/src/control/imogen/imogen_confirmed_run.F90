#if !defined(UM_JULES)
! Unlike IMOGEN_CHECK, this is a far more rigorous check, only allowing 
! have been checked and believed to function as expected. 
! In time as all the configurations are checked, then this routine will 
! identical checks as those in "IMOGEN_CHECK"
! C. Huntingford (5th October 2004)

  SUBROUTINE IMOGEN_CONFIRMED_RUN(                                &
    C_EMISSIONS,INCLUDE_CO2,INCLUDE_NON_CO2,LAND_FEED,OCEAN_FEED, &
    L_IMPACTS,WGEN,ANLG,ANOM                                      &
  )

    IMPLICIT NONE

    LOGICAL ::                                                    &
      C_EMISSIONS,                                                &
                  !IN If true, means CO2 concentration is calcula
      INCLUDE_CO2,                                                &
                  !IN Are adjustments to CO2 values allowed?
      INCLUDE_NON_CO2,                                            &
                  !IN Are adjustments to non-CO2 values allowed?
      LAND_FEED,                                                  &
                  !IN Are land feedbacks allowed on atmospheric C
      OCEAN_FEED,                                                 &
                  !IN Are ocean feedbacks allowed on atmospheric 
      CHECK_FLAG_EXTRA,                                           &
                  !WORK Confirms that configuration is OK  
      L_IMPACTS,                                                  &
                  !IN If true an impacts model used
      WGEN,                                                       &
                  !IN Weather generator
      ANLG,                                                       &
                  !IN True if the analogue model is used
      ANOM        !IN True if the anomalies are used 



    CHECK_FLAG_EXTRA=.FALSE.

! Transient run with carbon cycle 
!      IF(C_EMISSIONS.AND.INCLUDE_CO2.AND.INCLUDE_NON_CO2.AND. 
    IF(C_EMISSIONS .AND. INCLUDE_CO2 .AND. LAND_FEED   .AND.      &
       OCEAN_FEED  .AND. L_IMPACTS   .AND. (.NOT.WGEN) .AND.      &
       ANLG        .AND. ANOM) THEN
      WRITE(*,*) ' '
      WRITE(*,*) '***************************************************'
      WRITE(*,*) 'This is a full analogue model simulation' 
      WRITE(*,*) 'Carbon cycle driven by emissions' 
!        WRITE(*,*) 'Non-co2 forcings are prescribed'
      WRITE(*,*) 'Land and ocean feedbacks on the carbon cycle included'
      WRITE(*,*) 'Terrestrial carbon content calculated with user DGVM' 
      WRITE(*,*) 'Weather generator is switched off'
      WRITE(*,*) '***************************************************'
      WRITE(*,*) ' '
      CHECK_FLAG_EXTRA=.TRUE.
    ENDIF

! Spin-up to the above
!      IF(C_EMISSIONS.AND.INCLUDE_CO2.AND.INCLUDE_NON_CO2.AND. 
    IF(C_EMISSIONS .AND. INCLUDE_CO2 .AND. LAND_FEED   .AND.      &
       OCEAN_FEED  .AND. L_IMPACTS   .AND. (.NOT.WGEN) .AND.      &
       ANLG        .AND. (.NOT.ANOM)) THEN
      WRITE(*,*) ' '
      WRITE(*,*) '***************************************************'
      WRITE(*,*) 'SPIN UP for full analogue model simulation below:' 
      WRITE(*,*) 'Carbon cycle driven by emissions' 
!        WRITE(*,*) 'Non-co2 forcings are prescribed'
      WRITE(*,*) 'Land and ocean feedbacks on the carbon cycle included'
      WRITE(*,*) 'Terrestrial carbon content calculated with user DGVM' 
      WRITE(*,*) 'Weather generator is switched off'
      WRITE(*,*) '***************************************************'
      WRITE(*,*) ' '
      CHECK_FLAG_EXTRA=.TRUE.
    ENDIF

! "Hydrology 20th century simulations" 
    IF((.NOT.C_EMISSIONS)     .AND. INCLUDE_CO2      .AND.        &
       (.NOT.INCLUDE_NON_CO2) .AND. (.NOT.LAND_FEED) .AND.        &
       (.NOT.OCEAN_FEED)      .AND. L_IMPACTS        .AND.        &
       (.NOT.WGEN)            .AND. (.NOT.ANLG)      .AND.        &
       (.NOT.ANOM)) THEN
      WRITE(*,*) ' '
      WRITE(*,*) '***************************************************'
      WRITE(*,*) 'Hydrology of 20th Century simulation' 
      WRITE(*,*) 'Allows user to provide a file of CO2 concs.' 
      WRITE(*,*) 'Terrestrial response calculated with user SVAT' 
      WRITE(*,*) 'Weather generator is switched off'
      WRITE(*,*) '***************************************************'
      WRITE(*,*) ' '
      CHECK_FLAG_EXTRA=.TRUE.
    ENDIF

    IF((.NOT.C_EMISSIONS)     .AND. (.NOT.INCLUDE_CO2) .AND.      &
       (.NOT.INCLUDE_NON_CO2) .AND. (.NOT.LAND_FEED)   .AND.      &
       (.NOT.OCEAN_FEED)      .AND. L_IMPACTS          .AND.      &
       (.NOT.WGEN)            .AND. (.NOT.ANLG)        .AND.      &
       (.NOT.ANOM)) THEN
      WRITE(*,*) ' '
      WRITE(*,*) '***************************************************'
      WRITE(*,*) 'SPIN TO or CO2 FIXED for:'
      WRITE(*,*) 'Hydrology of 20th Century simulation' 
      WRITE(*,*) 'Allows user to provide a file of CO2 concs.' 
      WRITE(*,*) 'Terrestrial response calculated with user SVAT' 
      WRITE(*,*) 'Weather generator is switched off'
      WRITE(*,*) '***************************************************'
      WRITE(*,*) ' '
      CHECK_FLAG_EXTRA=.TRUE.
    ENDIF

! Analogue model simulations with CO2 prescribed
    IF((.NOT.C_EMISSIONS) .AND. (INCLUDE_CO2)    .AND.            &
       (INCLUDE_NON_CO2)  .AND. (.NOT.LAND_FEED) .AND.            &
       (.NOT.OCEAN_FEED)  .AND. L_IMPACTS        .AND.            &
       (.NOT.WGEN)        .AND. ANLG             .AND.            &
       ANOM) THEN
      WRITE(*,*) ' '
      WRITE(*,*) '***************************************************'
      WRITE(*,*) 'Transient analogue model simulation'
      WRITE(*,*) 'Allows user to provide a file of CO2 concs.' 
      WRITE(*,*) 'Terrestrial response calculated with user SVAT' 
      WRITE(*,*) 'Weather generator is switched off'
      WRITE(*,*) '***************************************************'
      WRITE(*,*) ' '
      CHECK_FLAG_EXTRA=.TRUE.
    ENDIF

    IF((.NOT.C_EMISSIONS)     .AND. (.NOT.INCLUDE_CO2) .AND.      &
       (.NOT.INCLUDE_NON_CO2) .AND. (.NOT.LAND_FEED)   .AND.      &
       (.NOT.OCEAN_FEED)      .AND. L_IMPACTS          .AND.      &
       (.NOT.WGEN)            .AND. ANLG               .AND.      &
       (.NOT.ANOM)) THEN
      WRITE(*,*) ' '
      WRITE(*,*) '***************************************************'
      WRITE(*,*) 'SPIN-UP solution to:'
      WRITE(*,*) 'Transient analogue model simulation'
      WRITE(*,*) 'Allows user to provide a file of CO2 concs.' 
      WRITE(*,*) 'Terrestrial response calculated with user SVAT' 
      WRITE(*,*) 'Weather generator is switched off'
      WRITE(*,*) '***************************************************'
      WRITE(*,*) ' '
      CHECK_FLAG_EXTRA=.TRUE.
    ENDIF

    IF (.NOT.ANOM) THEN
      CHECK_FLAG_EXTRA=.TRUE.
    ENDIF

    IF(.NOT.CHECK_FLAG_EXTRA) THEN
      WRITE(*,*) 'This configuration of the AM will be available' &
              // ' at some point, but is not yet checked' 
      STOP
    ENDIF

    RETURN

  END SUBROUTINE IMOGEN_CONFIRMED_RUN
#endif
