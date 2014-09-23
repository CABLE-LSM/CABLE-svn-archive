#if !defined(UM_JULES)
MODULE imogen_progs

  USE imogen_constants, ONLY : n_olevs,nfarray

  IMPLICIT NONE
  
  REAL ::                                                         &
    CO2_PPMV,                                                     &
              ! Atmospheric CO2 concentration (ppmv)
    CO2_START_PPMV,                                               &
              ! Atmospheric CO2 concentration at start of year (ppmv)
    CO2_CHANGE_PPMV
              ! Change in CO2 between restarts (ppmv)
                  
  REAL ::                                                         &
    DTEMP_O(N_OLEVS),                                             &
              ! Ocean temperature anomalies (K)
    FA_OCEAN(NFARRAY)
              ! CO2 fluxes from the atmosphere to the
              ! ocean (ie positive upwards) (ppm/m2/yr)
              
  INTEGER ::                                                      &
    SEED_WG(4),                                                   &
              ! Seeding number for weather generator
              ! random number
    SEED_RAIN(4)
              ! Seeding number for subdaily rainfall.
  

END MODULE imogen_progs
#endif
