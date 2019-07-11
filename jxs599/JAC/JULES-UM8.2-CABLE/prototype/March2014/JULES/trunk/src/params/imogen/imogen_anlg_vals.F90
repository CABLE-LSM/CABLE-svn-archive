#if !defined(UM_JULES)
  MODULE imogen_anlg_vals
  
    IMPLICIT NONE
    
    REAL ::                                                       &
      Q2CO2=3.74,                                                 &
                  ! Radiative forcing due to doubling CO2 (W/m2
      F_OCEAN=0.711,                                              &
                  ! Fractional coverage of the ocean
      KAPPA_O=383.8,                                              &
                  ! Ocean eddy diffusivity (W/m/K) 
      LAMBDA_L=0.52,                                              &
                  ! Inverse of climate sensitivity
                  !   over land (W/m2/K)
      LAMBDA_O=1.75,                                              &
                  ! Inverse of climate sensitivity
                  !   over ocean (W/m2/K)
      MU=1.87,                                                    &
                  ! Ratio of land to ocean temperature 
                  !   anomalies
      T_OCEAN_INIT=289.28
                  ! Initial ocean temperature (K)

    INTEGER ::                                                    &
      NYR_NON_CO2 = 21
                  ! Number of years for which NON_CO2 forcing
                  ! is prescribed.
                   
    CHARACTER(len=180) ::                                         &
      DIR_PATT = '',                                              &
                  ! Directory containing the patterns
      DIR_CLIM = '',                                              &
                  ! Directory containing initialising climatology.
      DIR_ANOM = ''
                  ! Directory containing prescribed anomalies
                  
    LOGICAL ::                                                    &
      FILE_NON_CO2=.FALSE.
                  ! If true, then non-CO2 radiative forcings ar
                  ! contained within a file.


    NAMELIST /IMOGEN_ANLG_VALS_LIST/ Q2CO2,F_OCEAN,KAPPA_O,       &
                                     LAMBDA_L,LAMBDA_O,MU,        &
                                     T_OCEAN_INIT,DIR_PATT,       &
                                     DIR_CLIM,NYR_NON_CO2,        &
                                     FILE_NON_CO2,DIR_ANOM
  
  END MODULE imogen_anlg_vals
#endif
