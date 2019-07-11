#if !defined(UM_JULES)
MODULE imogen_constants

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: N_OLEVS = 254
              ! Number of ocean levels in thermal calculation
  INTEGER, PARAMETER :: NFARRAY = 10000
              ! Array size for FA_OCEAN
  INTEGER, PARAMETER :: N_IMOGEN_LAND = 1631
              ! Number of land points in the IMOGEN grid

  REAL, PARAMETER :: MDI=999.9
              ! Missing data indicator
  REAL, PARAMETER :: OCEAN_AREA=3.627E14
              ! Ocean area (m2) 
  REAL, PARAMETER :: CONV=0.471
              ! Converts global emission of C (Gt) into
              !     change in atmospheric CO2 (ppm)
  
  CHARACTER(len=4) ::                                             &
    DRIVE_MONTH(12)   ! Month labels for reading in files
  
  DATA DRIVE_MONTH /'/jan','/feb','/mar','/apr','/may','/jun',    &
                    '/jul','/aug','/sep','/oct','/nov','/dec' /

END MODULE imogen_constants
#endif
