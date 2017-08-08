!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose:
!
! Called from: JULES: surf_couple_
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE-JULES coupling in UM 10.5
!
!
! ==============================================================================

module cable_explicit_main_mod
  
contains

SUBROUTINE cable_explicit_main(                                                &
            mype, timestep, timestep_number, cycleno, numcycles,               &
            ! grid, model, dimensions. PFT frac per landpoint    
            row_length, rows, land_pts, ntiles,                                &
            npft, sm_levels,                                                   &
            latitude, longitude,                                               &
            land_index, tile_frac, tile_pts, tile_index,                       &
            ! Soil parameters **jhan: could be undef.@some point  issue here
            bexp, hcon, satcon, sathh,                                         &
            smvcst, smvcwt, smvccl, albsoil,                                   & 
            ! packs/ unpacked ssnow% snowd 
            snow_tile, lw_down, cosine_zenith_angle,                           &
            surf_down_sw, ls_rain, ls_snow,                                    &          
            tl_1, qw_1, vshr_land, pstar, z1_tq, z1_uv, canopy_tile,           &
            Fland, CO2_MMR, sthu, canht_ft, lai_ft ,                           &
            sin_theta_latitude, dzsoil, FTL_TILE, FQW_TILE,                    &
            TSTAR_TILE, U_S, U_S_STD_TILE, CD_TILE, CH_TILE,                   &
            RADNET_TILE, FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,              &
            RECIP_L_MO_TILE, EPOT_TILE )
 
  USE cable_common_module, ONLY : knode_gl,        & ! processor number
                                  ktau_gl,         & ! number
                                  kend_gl,         & ! number
                                  kwidth_gl          ! width in S 
  
  USE atm_fields_real_mod, ONLY : soil_temp_cable, soil_moist_cable,           &
                                  soil_froz_frac_cable, snow_dpth_cable,       & 
                                  snow_mass_cable, snow_temp_cable,            &
                                  snow_rho_cable, snow_avg_rho_cable,          &   
                                  snow_age_cable, snow_flg_cable

  USE cable_explicit_driv_mod, ONLY : cable_explicit_driver

  USE cable_expl_unpack_mod, ONLY : cable_expl_unpack
  
  USE cable_decs_mod, only : L_tile_pts, rho_water

  USE cable_um_tech_mod, ONLY : air, bgc, canopy, met, bal, rad, rough, soil,  &
                                ssnow, sum_flux, veg
  
  implicit none
 
  !--- IN ARGS FROM sf_exch_cable, passed from surf_couple_explicit() down ----
   INTEGER ::                                                      & 
     mype, timestep_number, cycleno, numcycles
   real :: timestep
   INTEGER ::                                                      & 
      row_length, rows, & ! UM grid resolution
      land_pts,         & ! # of land points being processed
      ntiles,           & ! # of tiles 
      npft,             & ! # of plant functional types
      sm_levels          ! # of soil layers 

   REAL,  DIMENSION(row_length,rows) ::                             &
      latitude,   &
      longitude
   
   ! index of land points being processed
   INTEGER, DIMENSION(land_pts) :: land_index 

   REAL, DIMENSION(land_pts, ntiles) ::                            &
      tile_frac
   
   ! # of land points on each tile
   INTEGER,  DIMENSION(ntiles) :: tile_pts 

   INTEGER,  DIMENSION(land_pts, ntiles) ::                         & 
      tile_index   ! index of tile points being processed

   !___UM soil/snow/radiation/met vars
   REAL,  DIMENSION(land_pts) :: & 
      bexp,    & ! => parameter b in Campbell equation 
      hcon,    & ! Soil thermal conductivity (W/m/K).
      satcon,  & ! hydraulic conductivity @ saturation [mm/s]
      sathh,   &
      smvcst,  &
      smvcwt,  &
      smvccl,  &
      albsoil

   REAL,  DIMENSION(land_pts, ntiles) ::                         &
      snow_tile
   
   REAL,  DIMENSION(row_length,rows) ::                             &
      lw_down,    &
      cosine_zenith_angle
 
   REAL, DIMENSION(row_length, rows, 4) ::                         &
      surf_down_sw 
   
   REAL,  DIMENSION(row_length,rows) ::                             &
      ls_rain,    &
      ls_snow

   REAL,  DIMENSION(row_length,rows) ::                             &
      tl_1,       &
      qw_1,       &  
      vshr_land,  &
      pstar,      &
      z1_tq,      &
      z1_uv

   REAL, DIMENSION(land_pts, ntiles) ::                             &
      canopy_tile
   
   REAL,  DIMENSION(land_pts) :: & 
      fland
   
   REAL :: co2_mmr

   REAL, DIMENSION(land_pts, sm_levels) ::                         &
      sthu 
   
   REAL, DIMENSION(land_pts, npft) ::                              &
      canht_ft, lai_ft 
   
   REAL :: sin_theta_latitude(row_length,rows) 
     
   REAL,  DIMENSION(sm_levels) :: dzsoil

   REAL, DIMENSION(land_pts,ntiles) :: &
      FTL_TILE,   &  ! Surface FTL for land tiles     
      FQW_TILE       ! Surface FQW for land tiles     
   
   !___return temp and roughness, friction velocities/drags/ etc

   REAL, DIMENSION(land_pts,ntiles) :: &
      TSTAR_TILE
   
   REAL, DIMENSION(row_length,rows)  :: &
      U_S               ! Surface friction velocity (m/s)

   REAL, DIMENSION(land_pts,ntiles) :: &
      U_S_STD_TILE,     & ! Surface friction velocity
      CD_TILE,    &     ! Drag coefficient
      CH_TILE           ! Transfer coefficient for heat & moisture

   !___return miscelaneous 
   REAL,  DIMENSION(land_pts,ntiles) :: &
      RADNET_TILE,    & ! Surface net radiation
      FRACA,          & ! Fraction of surface moisture
      RESFS,          & ! Combined soil, stomatal & aerodynamic resistance
                        ! factor for fraction (1-FRACA) of snow-free land tiles
      RESFT,          & ! Total resistance factor.
                        ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
                        ! 1 for snow.    
      Z0H_TILE,       &
      Z0M_TILE,       &
      RECIP_L_MO_TILE,& ! Reciprocal of the Monin-Obukhov length for tiles (m^-1)
      EPOT_TILE

!real :: SOIL_ORDER_casa(:)
!real :: LAI_casa(:,:)
!real :: PHENPHASE_casa(:,:)
!real :: C_pool_casa(:,:,:)
!real :: N_pool_casa(:,:,:)
!real :: N_dep_casa(:)
!real :: N_fix_casa(:)
!real :: P_pool_casa(:,:,:)
!real :: P_dust_casa(:)
!real :: P_weath_casa(:)

!REAL,  DIMENSION(land_pts, ntiles,3) ::                       &
!   snow_cond
!
  
! rml 2/7/13 Extra atmospheric co2 variables
!LOGICAL, INTENT(IN) :: L_CO2_INTERACTIVE
!INTEGER, INTENT(IN) ::                              &
!  CO2_DIM_LEN                                      &
!  ,CO2_DIM_ROW
!REAL, INTENT(IN) :: CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)  ! co2 mass mixing ratio

  !--- End IN ARGS ----------------------------------------------------------

  !___true IF vegetation (tile) fraction is greater than 0
  integer,  DIMENSION(land_pts, ntiles) :: isnow_flg_cable
  
!   !r825 adds CASA vars here
!   REAL, DIMENSION(land_pts,ntiles,10) :: &
!      CPOOL_TILE,    & ! Carbon Pools
!      NPOOL_TILE       ! Nitrogen Pools
!
!   REAL, DIMENSION(land_pts,ntiles,12) :: &
!      PPOOL_TILE       ! Phosphorus Pools
!
!   REAL, DIMENSION(land_pts) :: &
!      SOIL_ORDER,    & ! Soil Order (1 to 12)
!      NIDEP,         & ! Nitrogen Deposition
!      NIFIX,         & ! Nitrogen Fixation
!      PWEA,          & ! Phosphorus from Weathering
!      PDUST            ! Phosphorus from Dust
!
!   REAL, DIMENSION(land_pts,ntiles) :: &
!      GLAI, &          ! Leaf Area Index for Prognostics LAI
!      PHENPHASE        ! Phenology Phase for Casa-CNP
!                                  
!   REAL, DIMENSION(land_pts,ntiles) :: &
!      NPP_FT_ACC,     &
!      RESP_W_FT_ACC

  !--- declare local vars ------------------------------------------------------ 

  character(len=*), parameter :: subr_name = "cable_explicit_main"
  logical, save :: first_call = .true.
  integer :: endstep = 0 !dummy 
  real :: radians_degrees
  REAL,  DIMENSION(row_length,rows) ::                             &
      latitude_deg,   &
      longitude_deg
   
  !--- End header -------------------------------------------------------------

19 format(  "CABLE_LSM: ", A20, " @", I8.1, " on ",I3.1  )
  write (6, 19)  subr_name, timestep_number, mype 

  !----------------------------------------------------------------------------
  !--- Organize report writing for CABLE.                         -------------
  !--- Progress log and IN args @ timestep X,Y,Z                  -------------
  !----------------------------------------------------------------------------
  !--- Convert lat/long to degrees
  radians_degrees = 180.0 / ( 4.0*atan(1.0) ) ! * 180 / PI
  latitude_deg  = latitude * radians_degrees
  longitude_deg  = longitude * radians_degrees
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  !initialize processor number, timestep width adn number here as we know this 
  !is first CABLE call in model. However make this more generic
  if( first_call ) then
    knode_gl  = mype 
    kwidth_gl = int(timestep)
  endif
  ktau_gl   = timestep_number
  kend_gl   = endstep   !dummy initialization
   
  if( .NOT. allocated(L_tile_pts) ) allocate( L_tile_pts(land_pts, ntiles) ) 

  !----------------------------------------------------------------------------
  !--- CALL _driver to run specific and necessary components of CABLE with IN -
  !--- args PACKED to force CABLE
  !----------------------------------------------------------------------------
  isnow_flg_cable = int(snow_flg_cable)

  call cable_explicit_driver( row_length, rows, land_pts, ntiles,npft,         &
                              sm_levels, timestep, latitude_deg, longitude_deg,&
                              land_index, tile_frac,  tile_pts, tile_index,    &
                              bexp, hcon, satcon, sathh, smvcst,               &
                              smvcwt,  smvccl, albsoil, snow_tile,             &
                              snow_avg_rho_cable, snow_age_cable,              &
                              isnow_flg_cable, snow_rho_cable, snow_dpth_cable,&
                              snow_temp_cable, snow_mass_cable,                &
                              lw_down, cosine_zenith_angle, surf_down_sw,      &
                              ls_rain, ls_snow, tl_1, qw_1, vshr_land, pstar,  &
                              z1_tq, z1_uv,  canopy_tile, Fland, CO2_MMR,      &
                              soil_moist_cable, soil_froz_frac_cable, sthu,    &
                              soil_temp_cable, canht_ft, lai_ft,               &
                              sin_theta_latitude, dzsoil, FTL_TILE, FQW_TILE,  &
                              TSTAR_TILE, U_S, U_S_STD_TILE, CD_TILE, CH_TILE, &
                              RADNET_TILE, FRACA, RESFS, RESFT,                &
                              Z0H_TILE, Z0M_TILE,                              &
                              RECIP_L_MO_TILE, EPOT_TILE,                      &
!                  !CPOOL_TILE, NPOOL_TILE, PPOOL_TILE,          &
!                  !SOIL_ORDER, NIDEP, NIFIX, PWEA, PDUST,       &
!                  !GLAI, PHENPHASE, NPP_FT_ACC, RESP_W_FT_ACC,  &
                              endstep, timestep_number, mype )    

  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !--- CALL _unpack to unpack variables from CABLE back to UM format to return
  !----------------------------------------------------------------------------
!*!    call cable_expl_unpack( latitude, longitude, FTL_TILE, FQW_TILE, TSTAR_TILE, &
!*!                           U_S, U_S_STD_TILE, &
!*!                           CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
!*!                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
!*!                           RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts,             &
!*!                           ssnow%snowd, ssnow%cls, air%rlam, air%rho,          &
!*!                           canopy%fe, canopy%fh, canopy%us, canopy%cdtq,       &
!*!                           canopy%fwet, canopy%wetfac_cs, canopy%rnet,         &
!*!                           canopy%zetar, canopy%epot, met%ua, rad%trad,        &
!*!                           rad%transd, rough%z0m, rough%zref_tq )
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  !----------------------------------------------------------------------------
  !--- Organize report writing for CABLE.                         -------------
  !--- OUT args @ timestep X,Y,Z                                  -------------
  !----------------------------------------------------------------------------

  !jhan: call checks as required by namelis      
  
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  
  first_call = .false.        

return

End subroutine cable_explicit_main
  
End module cable_explicit_main_mod











































