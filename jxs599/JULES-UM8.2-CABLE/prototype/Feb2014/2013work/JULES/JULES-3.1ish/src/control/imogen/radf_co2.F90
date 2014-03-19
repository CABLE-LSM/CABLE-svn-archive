#if !defined(UM_JULES)
!**********************************************************************
! Calculates the radiative forcing due to a given CO2 concentration 
!
! Written by Peter Cox (August 1998)
!**********************************************************************
  SUBROUTINE RADF_CO2(CO2,CO2REF,Q2CO2,Q_CO2)

    IMPLICIT NONE


    REAL ::                                                       &
      CO2,                                                        &
                ! IN CO2 concentration (ppmv).
      CO2REF,                                                     &
                ! IN Reference CO2 concentration (ppmv).
      Q2CO2,                                                      &
                ! IN Radiative forcing due to doubling CO2 (W/m2).
      Q_CO2     ! OUT Radiative forcing due to CO2 (W/m2).

    Q_CO2=(Q2CO2/(LOG(2.0)))*LOG(CO2/CO2REF)

    RETURN
  END SUBROUTINE RADF_CO2
#endif
