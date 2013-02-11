   
   subroutine cable_hyd_driver( &
        SNOW_TILE, LYING_SNOW, SURF_ROFF, SUB_SURF_ROFF, TOT_TFALL)

      use cable_common_module!, only : cable_runtime, cable_user
      use cable_um_tech_mod, only : um1, ssoil, canopy, veg
      use cable_diag_module, only : cable_stat
      implicit none

      REAL, intent(out), dimension(um1%LAND_PTS,um1%NTILES) :: &  
            SNOW_TILE   ! IN Lying snow on tiles (kg/m2)        

      REAL, intent(out), dimension(um1%LAND_PTS) :: &  
            LYING_SNOW       &  ! OUT Gridbox snowmass (kg/m2)        
            ,SUB_SURF_ROFF          &
            ,SURF_ROFF  &
            ,TOT_TFALL

      REAL, dimension(um1%LAND_PTS,um1%NTILES) :: &  
            SURF_CAB_ROFF               &
            ,TOT_TFALL_TILE                

      real :: miss =0. 
      
         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('cable_hyd_driver')
         
         cable_runtime%um_hydrology = .true.

         SNOW_TILE= unpack(ssoil%snowd, um1%L_TILE_PTS, miss)
         LYING_SNOW = SUM(um1%TILE_FRAC * SNOW_TILE,2) !gridbox snow mass

         SURF_CAB_ROFF  = unpack(ssoil%rnof1, um1%L_TILE_PTS, miss)
         SURF_ROFF      = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)
         ! Don't include sub-soil drainage for lakes
         where ( veg%iveg == 16 ) &
           ssoil%rnof2 = 0.0
     
         SURF_CAB_ROFF  = unpack(ssoil%rnof2, um1%L_TILE_PTS, miss)
         SUB_SURF_ROFF  = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)

         TOT_TFALL_TILE = unpack(canopy%through, um1%L_TILE_PTS, miss)
         TOT_TFALL      = SUM(um1%TILE_FRAC * TOT_TFALL_TILE,2)
         
         cable_runtime%um_hydrology = .false.

      return
   end subroutine cable_hyd_driver
      




