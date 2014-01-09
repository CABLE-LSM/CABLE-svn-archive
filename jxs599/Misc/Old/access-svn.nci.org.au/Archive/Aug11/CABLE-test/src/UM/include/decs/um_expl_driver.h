      !___IN: UM dimensions, array indexes, flags
      integer, intent(in) :: row_length, rows, land_pts, ntiles, npft, sm_levels
      integer, intent(in), dimension(land_pts) :: land_index 
      integer, intent(in), dimension(ntiles) :: tile_pts 
      integer, intent(in), dimension(land_pts, ntiles) :: tile_index ,isnow_flg3l 
      !--- T if L_TILE_PTS, F elsewhere.
      logical,dimension(row_length,rows) :: land_mask   

      !___UM parameters 
      real, intent(in) :: timestep
      real, intent(in) :: rho_water 
      real, intent(in), dimension(sm_levels) :: dzsoil

      !___UM soil/snow/radiation/met vars
      real, intent(in), dimension(land_pts) :: bexp, hcon, satcon, sathh, smvcst,&
            smvcwt, smvccl, albsoil, fland 
      real, intent(inout), dimension(row_length,rows) :: sw_down,  cos_zenith_angle
      real, intent(in), dimension(row_length,rows) :: latitude, longitude,       &
            lw_down, ls_rain, ls_snow, tl_1, qw_1,    &
            vshr_land, pstar, z1_tq, z1_uv
      real, intent(inout), dimension(land_pts, ntiles) :: snow_tile
      real, intent(in), dimension(land_pts, ntiles) :: tile_frac, &   
            snow_rho1l, snage_tile
      real, intent(in), dimension(row_length, rows, 4) :: surf_down_sw 
      real, intent(in), dimension(land_pts, npft) :: canht_ft, lai_ft 
      real, intent(in),dimension(land_pts, ntiles) :: canopy_tile
      real, intent(inout), dimension(land_pts, ntiles,3) :: snow_cond
      real, intent(in), dimension(land_pts, ntiles,3) :: &
            snow_rho3l, snow_depth3l, snow_mass3l, snow_tmp3l
      real, intent(in), dimension(land_pts, sm_levels) :: sthu 
      real, intent(in), dimension(land_pts, ntiles, sm_levels) :: sthu_tile,     &
            sthf_tile,smcl_tile, tsoil_tile
      real, intent(in) :: co2_mmr

      !___true IF vegetation (tile) fraction is greater than 0
      logical, intent(inout),dimension(land_pts, ntiles) :: L_tile_pts
     
      real :: sin_theta_latitude(row_length,rows) 

        
         !___return fluxes
         REAL, intent(out), dimension(land_pts) ::   &
            FTL_CAB, &
            LE_CAB
         REAL, intent(out), dimension(land_pts,ntiles) :: &
            FTL_TILE_CAB, &
            FTL_TILE,   &  ! Surface FTL for land tiles     
            FQW_TILE,   &  ! Surface FQW for land tiles     
            LE_TILE_CAB
      
         !___return temp and roughness
         real, intent(out), dimension(land_pts,ntiles) :: &
            TSTAR_TILE_CAB, TSTAR_TILE,  Z0H_TILE, Z0M_TILE
         real, intent(out), dimension(land_pts) ::                  &
            TSTAR_CAB
      
         !___return friction velocities/drags/ etc
         real, intent(out), dimension(land_pts,ntiles) :: &
            CD_TILE,    &     ! Drag coefficient
            CH_TILE,    &     ! Transfer coefficient for heat & moisture
            U_S_STD_TILE      ! Surface friction velocity
         real, intent(out), dimension(row_length,rows)  :: &
            U_S               ! Surface friction velocity (m/s)
         real, intent(out), dimension(land_pts) ::                  &
            CH_CAB,  &  ! Turbulent surface exchange
            CD_CAB,  &  ! Turbulent surface exchange
            U_S_CAB     ! Surface friction velocity (m/s)
      
         !___return miscelaneous 
         real, intent(out), dimension(land_pts,ntiles) :: &
            RADNET_TILE,   &  ! Surface net radiation
            RESFS,   &        ! Combined soil, stomatal & aerodynamic resistance
                              ! factor for fraction (1-FRACA) of snow-free land tiles
            RESFT,   &        ! Total resistance factor.
                              ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts, 1 for snow.    
            FRACA,   &        ! Fraction of surface moisture
            RECIP_L_MO_TILE,  & ! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
            EPOT_TILE
         !----------------------------------------------------------------------------   
      
   
     
  
