
         !___UM vars returned from land-surface (CABLE)
         !----------------------------------------------------------------------------   
      
         !___return fluxes
         REAL, intent(out), dimension(um1%land_pts) ::   &
            FTL_CAB, &
            LE_CAB
         REAL, intent(out), dimension(um1%land_pts,um1%ntiles) :: &
            FTL_TILE_CAB, &
            FTL_TILE,   &  ! Surface FTL for land tiles     
            FQW_TILE,   &  ! Surface FQW for land tiles     
            LE_TILE_CAB
      
         !___return temp and roughness
         real, intent(out), dimension(um1%land_pts,um1%ntiles) :: &
            TSTAR_TILE_CAB, TSTAR_TILE,  Z0H_TILE, Z0M_TILE
         real, intent(out), dimension(um1%land_pts) ::                  &
            TSTAR_CAB
      
         !___return friction velocities/drags/ etc
         real, intent(out), dimension(um1%land_pts,um1%ntiles) :: &
            CD_TILE,    &     ! Drag coefficient
            CH_TILE,    &     ! Transfer coefficient for heat & moisture
            U_S_STD_TILE      ! Surface friction velocity
         real, intent(out), dimension(um1%row_length,um1%rows)  :: &
            U_S               ! Surface friction velocity (m/s)
         real, intent(out), dimension(um1%land_pts) ::                  &
            CH_CAB,  &  ! Turbulent surface exchange
            CD_CAB,  &  ! Turbulent surface exchange
            U_S_CAB     ! Surface friction velocity (m/s)
      
         !___return miscelaneous 
         real, intent(out), dimension(um1%land_pts,um1%ntiles) :: &
            RADNET_TILE,   &  ! Surface net radiation
            RESFS,   &        ! Combined soil, stomatal & aerodynamic resistance
                              ! factor for fraction (1-FRACA) of snow-free land tiles
            RESFT,   &        ! Total resistance factor.
                              ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts, 1 for snow.    
            FRACA,   &        ! Fraction of surface moisture
            RECIP_L_MO_TILE,  & ! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
            EPOT_TILE
         
         logical,dimension(um1%land_pts,um1%ntiles) :: l_tile_pts
         !----------------------------------------------------------------------------   
      
         !___UM vars used but NOT returned 
         real, intent(in), dimension(um1%land_pts) ::   &
            FLAND(um1%land_pts)              ! IN Land fraction on land tiles.
      
