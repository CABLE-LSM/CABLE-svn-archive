   
SUBROUTINE cable_hyd_driver( SNOW_TILE, LYING_SNOW, SURF_ROFF, SUB_SURF_ROFF,  &
                             TOT_TFALL )

   USE cable_data_module,   ONLY : PHYS, OTHER
   USE cable_common_module!, only : cable_runtime, cable_user
   USE cable_um_tech_mod, only : um1, ssnow, canopy, veg
   IMPLICIT NONE

   REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS,um1%NTILES) ::                    &
      SNOW_TILE   ! IN Lying snow on tiles (kg/m2)        

   REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS) ::                               &
      LYING_SNOW,    & ! OUT Gridbox snowmass (kg/m2)        
      SUB_SURF_ROFF, & !
      SURF_ROFF,     & !
      TOT_TFALL        !

   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                                 &
      SURF_CAB_ROFF,    &
      TOT_TFALL_TILE                

   REAL :: miss =0. 
   REAL, POINTER :: TFRZ
      
      TFRZ => PHYS%TFRZ
   
      IF( cable_user%RUN_DIAG_LEVEL == 'BASIC' )                               &
         CALL cable_stat('cable_hyd_driver')
      
      cable_runtime%um_hydrology = .true.

      SNOW_TILE= UNPACK(ssnow%snowd, um1%L_TILE_PTS, miss)
      LYING_SNOW = SUM(um1%TILE_FRAC * SNOW_TILE,2) !gridbox snow mass

      SURF_CAB_ROFF  = UNPACK(ssnow%rnof1, um1%L_TILE_PTS, miss)
      SURF_ROFF      = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)
      
      ! Don't include sub-soil drainage for lakes
      WHERE( veg%iveg == 16 ) ssnow%rnof2 = 0.0
  
      SURF_CAB_ROFF  = UNPACK(ssnow%rnof2, um1%L_TILE_PTS, miss)
      SUB_SURF_ROFF  = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)

      TOT_TFALL_TILE = UNPACK(canopy%through, um1%L_TILE_PTS, miss)
      TOT_TFALL      = SUM(um1%TILE_FRAC * TOT_TFALL_TILE,2)
      
      cable_runtime%um_hydrology = .FALSE.

END SUBROUTINE cable_hyd_driver
      




