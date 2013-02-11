!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++ CABLE HEADER 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
!=== CABLE is called in the UM from 4 locations. this is the first call at  
!=== the beginning of the timestep from sf_exch(). others are:
!===     cable_implicit_driver() from sf_impl(). 
!===     cable_rad_driver from glue_ctl_rad()
!===     cable_hyd_driver from hydrol()
!=== the main purpose of this call is to calculate the surface-atmosphere 
!=== exchange co-effs so that the UM can convect etc. at mid-timestep calls
!=== cable_implicit() where it uses updated forcings. 
!===============================================================================


subroutine cable_explicit_driver(  row_length, rows, land_pts, ntiles,npft,  & 
           sm_levels, timestep, latitude, longitude, land_index, tile_frac,  &
           tile_pts, tile_index, bexp, hcon, satcon, sathh, smvcst, smvcwt,  &
           smvccl, albsoil, snow_tile, snow_rho1l, snage_tile, isnow_flg3l,  &
           snow_rho3l, snow_cond, snow_depth3l, snow_tmp3l, snow_mass3l,     &
           sw_down, lw_down, cos_zenith_angle, surf_down_sw, ls_rain,        &
           ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, rho_water,   &
           L_tile_pts, canopy_tile, Fland, CO2_MMR, sthu_tile, smcl_tile,    &
           sthf_tile, sthu, tsoil_tile, canht_ft, lai_ft, sin_theta_latitude,&
           dzsoil, LAND_MASK, FTL_TILE_CAB, FTL_CAB, FTL_TILE, FQW_TILE,     &
           LE_TILE_CAB, LE_CAB, TSTAR_TILE, TSTAR_TILE_CAB, TSTAR_CAB, U_S,  &
           U_S_STD_TILE, U_S_CAB, CH_CAB, CD_CAB, CD_TILE, CH_TILE,          &
           RADNET_TILE, FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,             &
           RECIP_L_MO_TILE, EPOT_TILE, endstep, timestep_number, mype )    
   
   !--- reads runtime and user switches and reports
   use cable_um_tech_mod, only : cable_um_runtime_vars, air, bgc, canopy,     &
                              met, bal, rad, rough, soil, ssoil, sum_flux, veg 
   
   !--- vars common to CABLE declared 
   use cable_common_module, only : cable_runtime, cable_user, ktau_gl, knode_gl, &
      kwidth_gl, kend_gl
   
   !--- CABLE diagnostic tools 
   use cable_diag_module, only : cable_stat
   
   !--- subr to (manage)interface UM data to CABLE
   use cable_um_init_mod, only : interface_UM_data
   
   !--- subr to call CABLE model
   use cbm_module, only : cbm

   use define_dimensions, only : mp,ms

   implicit none
 
 
 
  
   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM sf_exch() --------------------------------------------
   !-------------------------------------------------------------------------- 
   
   !___IN: UM dimensions, array indexes, flags
   integer, intent(in) :: & 
      row_length, rows, & ! UM grid resolution
      land_pts,         & ! # of land points being processed
      ntiles,           & ! # of tiles 
      npft,             & ! # of plant functional types
      sm_levels           ! # of soil layers 

   ! index of land points being processed
   integer, intent(in), dimension(land_pts) :: land_index 

   ! # of land points on each tile
   integer, intent(in), dimension(ntiles) :: tile_pts 
   
   integer, intent(in), dimension(land_pts, ntiles) :: & 
      tile_index ,& ! index of tile points being processed
      isnow_flg3l   ! 3 layer snow flag

   !--- TRUE if land, F elsewhere.
   !jhan:rm land_mask
   logical,dimension(row_length,rows) :: land_mask   

   !___UM parameters: water density, soil layer thicknesses 
   real, intent(in) :: rho_water 
   real, intent(in), dimension(sm_levels) :: dzsoil

   !___UM soil/snow/radiation/met vars
   real, intent(in), dimension(land_pts) :: & 
      bexp,    & ! => parameter b in Campbell equation 
      hcon,    & ! Soil thermal conductivity (W/m/K).
      satcon,  & ! hydraulic conductivity @ saturation [mm/s]
      sathh,   &
      smvcst,  &
      smvcwt,  &
      smvccl,  &
      albsoil, &
      fland 
   
   real, intent(inout), dimension(row_length,rows) :: &
      sw_down,          & 
      cos_zenith_angle
   
   real, intent(in), dimension(row_length,rows) :: latitude, longitude,       &
      lw_down,    &
      ls_rain,    &
      ls_snow,    &
      tl_1,       &
      qw_1,       &  
      vshr_land,  &
      pstar,      &
      z1_tq,      &
      z1_uv

   real, intent(inout), dimension(land_pts, ntiles) :: snow_tile

   real, intent(in), dimension(land_pts, ntiles) :: tile_frac, &    
      snow_rho1l, &
      snage_tile
   
   real, intent(in), dimension(row_length, rows, 4) :: surf_down_sw 
   
   real, intent(in), dimension(land_pts, npft) :: canht_ft, lai_ft 
   
   real, intent(in),dimension(land_pts, ntiles) :: canopy_tile
   
   real, intent(inout), dimension(land_pts, ntiles,3) :: snow_cond
   
   real, intent(in), dimension(land_pts, ntiles,3) :: &
      snow_rho3l,    &
      snow_depth3l,  &
      snow_mass3l,   &
      snow_tmp3l
   
   real, intent(in), dimension(land_pts, sm_levels) :: sthu 
   
   real, intent(in), dimension(land_pts, ntiles, sm_levels) :: & 
      sthu_tile, &
      sthf_tile, &
      smcl_tile, &
      tsoil_tile
   
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
      TSTAR_TILE_CAB,   &
      TSTAR_TILE,       & 
      Z0H_TILE,         &
      Z0M_TILE

   real, intent(out), dimension(land_pts) ::  TSTAR_CAB

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

   ! end step of experiment, this step, step width, processor num
   integer, intent(in) :: endstep, timestep_number, mype
   real, intent(in) ::  timestep     
   
   integer:: itimestep
    
   !___return miscelaneous 
   real, intent(out), dimension(land_pts,ntiles) :: &
      RADNET_TILE,   & ! Surface net radiation
      RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
                       ! factor for fraction (1-FRACA) of snow-free land tiles
      RESFT,         & ! Total resistance factor.
                       ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts, 1 for snow.    
      FRACA,         & ! Fraction of surface moisture
      RECIP_L_MO_TILE,  & ! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
      EPOT_TILE
     
   !-------------------------------------------------------------------------- 
   !--- end INPUT ARGS FROM sf_exch() ----------------------------------------
   !-------------------------------------------------------------------------- 
   


   
   !___ declare local vars 
   
   !___ location of namelist file defining runtime vars
   character(len=200), parameter ::   & 
         runtime_vars_file = '/home/599/ewk599/CABLE-UM/cable.nml'

   !___ 1st call in RUN (!=ktau_gl -see below) 
   logical, save :: first_cable_call = .true.
 



      if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         call cable_stat('cable_explicit_driver')

      !--- initialize cable_runtime% switches 
      if(first_cable_call) & 
         cable_runtime%um = .true.
      
      !--- basic info from global model passed to cable_common_module 
      !--- vars so don't need to be passed around, just USE _module
      ktau_gl = timestep_number     !timestep of EXPERIMENT not necesarily 
                                    !the same as timestep of particular RUN
      knode_gl = mype               !which processor am i on?
      itimestep = int(timestep)    !realize for 'call cbm' pass
      kwidth_gl = itimestep          !width of timestep (secs)
      kend_gl = endstep             !timestep of EXPERIMENT not necesarily 

      !--- internal FLAGS def. specific call of CABLE from UM
      !--- from cable_common_module
      cable_runtime%um_explicit = .true.

      !--- user FLAGS, variables etc def. in cable.nml is read on 
      !--- first time step of each run. these variables are read at 
      !--- runtime and for the most part do not require a model rebuild.
      if(first_cable_call) then
         call cable_um_runtime_vars(runtime_vars_file) 
         first_cable_call = .false.
      endif      




      !---------------------------------------------------------------------!
      !--- initialize CABLE using UM forcings etc. these args are passed ---!
      !--- down from UM.                                                 ---! 
      !---------------------------------------------------------------------!
      call interface_UM_data( row_length, rows, land_pts, ntiles,    & 
               npft, sm_levels, itimestep, latitude, longitude, land_index,       &
               tile_frac, tile_pts, tile_index, bexp, hcon, satcon, sathh,       &
               smvcst, smvcwt, smvccl, albsoil, snow_tile, snow_rho1l,           &
               snage_tile, isnow_flg3l, snow_rho3l, snow_cond, snow_depth3l,     &
               snow_tmp3l, snow_mass3l, sw_down, lw_down, cos_zenith_angle,      &
               surf_down_sw, ls_rain, ls_snow, tl_1, qw_1, vshr_land, pstar,     &
               z1_tq, z1_uv, rho_water, L_tile_pts, canopy_tile, Fland,          &
               CO2_MMR, sthu_tile, smcl_tile, sthf_tile, sthu, tsoil_tile,       &
               canht_ft, lai_ft, sin_theta_latitude, dzsoil )                        

      canopy%oldcansto=canopy%cansto


      !---------------------------------------------------------------------!
      !--- real(timestep) width, CABLE types passed to CABLE "engine" as ---!  
      !--- req'd by Mk3L  --------------------------------------------------!
      !---------------------------------------------------------------------!
      CALL cbm(timestep, air, bgc, canopy, met, bal, &
                  rad, rough, soil, ssoil, sum_flux, veg )




      !---------------------------------------------------------------------!
      !--- pass land-surface quantities calc'd by CABLE in explicit call ---!
      !--- back to UM.                                                   ---!
      !---------------------------------------------------------------------!
      call cable_expl_unpack( FTL_TILE_CAB, FTL_CAB, FTL_TILE, FQW_TILE,         &
               LE_TILE_CAB, LE_CAB, TSTAR_TILE, TSTAR_TILE_CAB, TSTAR_CAB,       &
               U_S, U_S_STD_TILE, U_S_CAB, CH_CAB, CD_CAB, CD_TILE, CH_TILE,     &
               FLAND, RADNET_TILE, FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,      &
               RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts, ssoil%snowd, ssoil%cls,   &
               air%rlam, air%rho, canopy%fe, canopy%fh, canopy%us, canopy%cdtq,  &
               canopy%fwet, canopy%wetfac_cs, canopy%rnet, canopy%zetar,         &
               canopy%epot, met%ua, rad%trad, rad%transd, rough%z0m,             &
               rough%zref_tq )

      cable_runtime%um_explicit = .false.




   return
end subroutine cable_explicit_driver




!---------------------------------------------------------------------!
!--- pass land-surface quantities calc'd by CABLE in explicit call ---!
!--- back to UM.                                                   ---!
!---------------------------------------------------------------------!
subroutine cable_expl_unpack( FTL_TILE_CAB, FTL_CAB, FTL_TILE, FQW_TILE,         &
               LE_TILE_CAB, LE_CAB, TSTAR_TILE, TSTAR_TILE_CAB, TSTAR_CAB,       &
               U_S, U_S_STD_TILE, U_S_CAB, CH_CAB, CD_CAB, CD_TILE, CH_TILE,     &
               FLAND, RADNET_TILE, FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,      &
               RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts, ssoil_snowd, ssoil_cls,   &
               air_rlam, air_rho, canopy_fe, canopy_fh, canopy_us, canopy_cdtq,  & 
               canopy_fwet, canopy_wetfac_cs, canopy_rnet, canopy_zetar,         &
               canopy_epot, met_ua, rad_trad, rad_transd, rough_z0m,             &
               rough_zref_tq )

   use define_dimensions, only : mp 
   use physical_constants, only : niter, capp
   use cable_um_tech_mod, only : um1
   use cable_common_module, only : cable_runtime, cable_user, &
                                   ktau_gl, knode_gl 
   use cable_diag_module, only : cable_stat, cable_diag
   implicit none         




   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
   !-------------------------------------------------------------------------- 



   
   !___ UM variables to recieve unpacked CABLE vars

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

   !___UM vars used but NOT returned 
   real, intent(in), dimension(um1%land_pts) ::   &
      FLAND(um1%land_pts)              ! IN Land fraction on land tiles.




   !___ decs of intent(in) CABLE variables to be unpacked

   ! snow depth (liquid water), factor for latent heat
   real, intent(in), dimension(mp) :: ssoil_snowd, ssoil_cls
   ! surface wind speed (m/s)
   real, intent(in), dimension(mp) :: met_ua 
   ! latent heat for water (j/kg), dry air density (kg m-3)
   real, intent(in), dimension(mp) :: air_rlam, air_rho 
   ! frac SW diffuse transmitted thru canopy, rad. temp. (soil and veg)
   real, intent(in), dimension(mp) :: rad_trad,rad_transd 
   ! total latent heat (W/m2), total sensible heat (W/m2)
   real, intent(in), dimension(mp) :: canopy_fe, canopy_fh  
   ! fraction of canopy wet
   real, intent(in), dimension(mp) :: canopy_fwet, canopy_wetfac_cs
   ! friction velocity, drag coefficient for momentum
   real, intent(in), dimension(mp) :: canopy_us, canopy_cdtq
   ! net rad. absorbed by surface (W/m2), total potential evaporation 
   real, intent(in), dimension(mp) :: canopy_rnet, canopy_epot        
   ! stability correction
   real, intent(in), dimension(mp,niter) :: canopy_zetar
   ! roughness length, Reference height for met forcing
   real, intent(in), dimension(mp) :: rough_z0m, rough_zref_tq 
 
   !-------------------------------------------------------------------------- 
   !--- end INPUT ARGS FROM cable_explicit_driver() ------------------------------
   !-------------------------------------------------------------------------- 


   
        
   !___vars in local calc. of latent heat fluxes
   REAL, dimension(um1%land_pts,um1%ntiles) ::                  &
      FQW_TILE_CAB,  &
      LE_TILE

   !___vars in local calc of Surface friction velocities
   REAL, dimension(um1%land_pts,um1%ntiles) ::                  &
      CD_CAB_TILE,   &  
      CH_CAB_TILE,   &  ! (bulk transfer) coeff. for momentum
      U_S_TILE
   REAL, dimension(mp)  :: &
      CDCAB,CHCAB

   !___local miscelaneous
   real, dimension(mp)  :: &
   THETAST,fraca_cab,rfsfs_cab, RECIPLMOTILE, fe_dlh
   integer :: i,j,k,N,L
   real :: miss = 0.0
   logical, save :: first_cable_call = .true.




      if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         call cable_stat('cable_expl_unpack')



      
      !___return fluxes
      FTL_TILE_CAB = unpack(canopy_fh,  um1%l_tile_pts, miss)
      FTL_CAB = SUM(um1%TILE_FRAC * FTL_TILE_CAB,2)
      FQW_TILE_CAB = unpack(canopy_fe,  um1%l_tile_pts, miss)
      LE_TILE_CAB = unpack(canopy_fe,  um1%l_tile_pts, miss)
      LE_CAB = SUM(um1%TILE_FRAC * LE_TILE_CAB,2)
      fe_dlh = canopy_fe/(air_rlam*ssoil_cls)
      FTL_TILE = unpack(canopy_fh,  um1%l_tile_pts, miss)
      FTL_TILE = FTL_TILE / capp
      FQW_TILE = unpack(fe_dlh, um1%l_tile_pts, miss)

      !jhan: for testing purposes only
      !call cable_diag( 1, 'canopy_flux', um1%land_pts, 48, ktau_gl,  & 
      !      knode_gl, 'cnp_flux', ftl_cab + le_cab )
      
      
      !___return temp and roughness
      TSTAR_TILE_CAB = unpack(rad_trad, um1%l_tile_pts, miss)
      TSTAR_CAB = SUM(um1%TILE_FRAC * TSTAR_TILE_CAB,2)
      TSTAR_TILE = unpack(rad_trad,  um1%l_tile_pts, miss)
      Z0M_TILE = unpack(rough_z0m,  um1%l_tile_pts, miss)
      Z0H_TILE = Z0M_TILE
      
      !___return friction velocities/drags/ etc
      U_S_TILE  =  unpack(canopy_us, um1%l_tile_pts, miss)
      U_S_CAB  = SUM(um1%TILE_FRAC *  U_S_TILE,2)
      CDCAB = canopy_us**2/met_ua**2   ! met%ua is always above umin = 0.1m/s
      ! for Cable CD*
      CD_CAB_TILE =  unpack(CDCAB,um1%l_tile_pts, miss)
      CD_CAB= SUM(um1%TILE_FRAC * CD_CAB_TILE,2)
      ! for Cable CH*
      CH_CAB_TILE =  unpack(canopy_cdtq,um1%l_tile_pts, miss)
      CH_CAB= SUM(um1%TILE_FRAC * CH_CAB_TILE,2)

      U_S_STD_TILE=U_S_TILE
      CD_TILE = CD_CAB_TILE
      CH_TILE = CH_CAB_TILE

      U_S = 0.
      DO N=1,um1%ntiles
      DO K=1,um1%TILE_PTS(N)
        L = um1%TILE_INDEX(K,N)
        J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
        I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
        U_S(I,J) = U_S(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*U_S_TILE(L,N)
      ENDDO
      ENDDO




      !___return miscelaneous 
      fraca_cab = canopy_fwet * (1.-rad_transd)
      where ( ssoil_snowd > 1.0) fraca_cab = 1.0
      rfsfs_cab = min(1.,max(0.01,canopy_wetfac_cs - fraca_cab)/max(0.01,1. - fraca_cab) )
      FRACA = unpack(fraca_cab, um1%l_tile_pts, miss)
      RESFT =  unpack( canopy_wetfac_cs,um1%l_tile_pts, miss)
      RESFS = unpack( rfsfs_cab , um1%l_tile_pts, miss)

      RADNET_TILE = unpack( canopy_rnet , um1%l_tile_pts, miss)
      THETAST = abs(canopy_fh)/(air_rho*capp*canopy_us)
      RECIPLMOTILE =  canopy_zetar(:,niter) / rough_zref_tq
      RECIP_L_MO_TILE = unpack(RECIPLMOTILE,um1%l_tile_pts, miss)
      EPOT_TILE = unpack(canopy_epot,um1%l_tile_pts, miss)
      



      if(first_cable_call) then 
         l_tile_pts = um1%l_tile_pts
         first_cable_call = .false.
      endif




   return
end subroutine cable_expl_unpack
    
!============================================================================
!============================================================================
!============================================================================

