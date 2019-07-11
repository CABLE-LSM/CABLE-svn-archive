#if !defined(UM_JULES)
  SUBROUTINE RESPONSE(NCALLYR,YEAR_RUN,RS)

    IMPLICIT NONE

    INTEGER ::                                                    &
      I,                                                          &
              !WORK looping parameter
      YEAR_RUN,                                                   &
              !IN Number of years in the simulation
      NCALLYR
              !IN Number of calls per year to the ocean routine
              !   (with CO2 conc. fixed)

    REAL ::                                                       &
      TIME_RS,                                                    &
              !IN Time delay required by the response function (yr)
      RS(NCALLYR*YEAR_RUN)
              !OUT Value of the response function (.)

 
    DO I = 1,NCALLYR*YEAR_RUN
      TIME_RS = FLOAT(I) / FLOAT(NCALLYR)
      IF(TIME_RS >= 1.0) THEN
        RS(I) = 0.014819 + 0.70367*EXP(-TIME_RS/0.70177)          &
                         + 0.24966*EXP(-TIME_RS/2.3488)           &
                         + 0.066485*EXP(-TIME_RS/15.281)          &
                         + 0.038344*EXP(-TIME_RS/65.359)          &
                         + 0.019439*EXP(-TIME_RS/347.55)
      ELSE
        RS(I) = 1.0 - 2.2617*TIME_RS + 14.002*(TIME_RS**2)        &
                    - 48.770*(TIME_RS**3) + 82.986*(TIME_RS**4)   &
                    - 67.527*(TIME_RS**5) + 21.037*(TIME_RS**6)
      ENDIF 
    ENDDO

    RETURN

  END SUBROUTINE RESPONSE
#endif
