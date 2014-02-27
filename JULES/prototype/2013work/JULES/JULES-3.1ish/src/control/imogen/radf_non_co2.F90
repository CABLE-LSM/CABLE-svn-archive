#if !defined(UM_JULES)
!**********************************************************************
! Calculates the radiative forcing due to non-CO2 GHGs. This version
! is consistent with the HadCM3 GHG run (AAXZE).
!
! Written by Peter Cox (Sept 1998)
! Adjusted by Chris Huntingford (Dec 1999) to include 
! other scenarios copied from file.
!**********************************************************************
  SUBROUTINE RADF_NON_CO2(                                        &
    YEAR,Q,NYR_NON_CO2,FILE_NON_CO2,FILE_NON_CO2_VALS             &
  )

    USE io_constants, ONLY : IMOGEN_UNIT

    IMPLICIT NONE

    INTEGER ::                                                    &
      YEAR,                                                       &
               ! IN Julian year.
      I,                                                          &
               ! WORK Loop counter.
      NYR_NON_CO2
               ! IN Number of years for which NON_CO2
               !    forcing is prescribed.

    CHARACTER(len=180) ::                                         &
      FILE_NON_CO2_VALS   ! IN File of non-co2 radiative forcings

    LOGICAL ::                                                    &
      FILE_NON_CO2        ! .T. if non_co2 forcings are read in from
                          ! a data file

    REAL :: Q             ! OUT non-CO2 radiative forcing (W/m2).

!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------

    INTEGER ::                                                    &
      YEARS(300)          ! Years for which radiative is
                          ! prescribed.

    REAL ::                                                       &
      Q_NON_CO2(300),                                             &
                          ! Specified radiative forcings (W/m2).
      GROWTH_RATE         ! Growth rate in the forcing after the
                          ! last prescribed value (%/pa).
    PARAMETER (GROWTH_RATE=0.0)

    DATA YEARS   / 1859,  1875,  1890,  1900,  1917,              &
                   1935,  1950,  1960,  1970,  1980,              &
                   1990,  2005,  2020,  2030,  2040,              &
                   2050,  2060,  2070,  2080,  2090,              &  
                   2100,  279*2100 /

    DATA Q_NON_CO2/ 0.0344,0.0557,0.0754,0.0912,0.1176,           &
                    0.1483,0.1831,0.2387,0.3480,0.4987,           &
                    0.6627,0.8430,0.9225,0.9763,1.0575,           &
                    1.1486,1.2316,1.3025,1.3604,1.4102,           &
                    1.4602,279*1.4602 /
                    
                    


    IF(.NOT.FILE_NON_CO2 .AND. NYR_NON_CO2 /= 21) THEN
      PRINT *,'Reset value of NYR_NON_CO2'
      STOP
    ENDIF

    IF(NYR_NON_CO2 > 300) THEN
      PRINT *,'NYR_NON_CO2 too large'
      STOP
    ENDIF


!-----------------------------------------------------------------------
! File of non-co2 forcings read in if required
!-----------------------------------------------------------------------
    IF(FILE_NON_CO2) THEN
      OPEN(IMOGEN_UNIT, FILE=FILE_NON_CO2_VALS,                               &
                        STATUS='old', POSITION='rewind', ACTION='read')
      DO I = 1,NYR_NON_CO2
        READ(IMOGEN_UNIT,*) YEARS(I),Q_NON_CO2(I)
      ENDDO
      CLOSE(IMOGEN_UNIT)
    ENDIF

!-----------------------------------------------------------------------
! Now calculate the non_co2 forcing 
!-----------------------------------------------------------------------
    IF (YEAR < YEARS(1)) THEN
      Q = 0.0
    ELSEIF (YEAR > YEARS(NYR_NON_CO2)) THEN
      Q = Q_NON_CO2(NYR_NON_CO2)                                  &
        * ((1.0 + 0.01*GROWTH_RATE)**(YEAR - YEARS(NYR_NON_CO2)))
    ELSE
      DO I=1,NYR_NON_CO2 - 1
        IF((YEAR >= YEARS(I)) .AND. (YEAR <= YEARS(I+1))) THEN
          Q = Q_NON_CO2(I) + (YEAR-YEARS(I))                      & 
            * (Q_NON_CO2(I+1)-Q_NON_CO2(I))/ (YEARS(I+1)-YEARS(I))
        ENDIF
      ENDDO
    ENDIF

    RETURN

  END SUBROUTINE RADF_NON_CO2
#endif
