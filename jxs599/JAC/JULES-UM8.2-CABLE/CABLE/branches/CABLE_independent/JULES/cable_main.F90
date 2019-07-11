!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
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
! Purpose: Passes UM variables to CABLE, calls cbm, passes CABLE variables 
!          back to UM. 'Explicit' is the first of two routines that call cbm at 
!          different parts of the UM timestep.
!
! Called from: UM code sf_exch
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
!
! ==============================================================================

MODULE cable_main_mod

CONTAINS

  SUBROUTINE cable_main(                                           &
           ! vars native to JULES/UM 
           !
           ! # grid cells 
           row_length, rows, &
           ! # AND index of land points
           land_pts, land_index, &
           ! model levels (modified for CABLE)
           ntiles, npft, sm_levels,   &
           ! time info 
           timestep, timestep_number, &
           ! time info in UM (passing dummy in JULES) 
           endstep, &
           ! # of processor in UM (passing dummy in JULES) 
           mype, &
           ! grid cell data
           latitude, longitude, &
           ! fraction of land on each grid cell 
           Fland, & 
           ! fraction of each tile (modified for CABLE)
           tile_frac,  &
           ! # AND index of tile points following tile_frac
           tile_pts, tile_index,&
           ! canopy height, LAI, soil levels 
           ! as prescribed in JULES 
           canht_ft,       &
           lai_ft, dzsoil,       &
           ! soil vars used in initialization of CABLE vars
           ! INTENT(IN)  used on first CALL only
           !jhan: ultimately read in as tiled vars
           bexp, hcon, satcon, sathh, smvcst,           &
           smvcwt, smvccl, albsoil, &
           ! JULES non-tiled unfrozen soil. returned from CABLE 
           sthu, &
           ! used in initialization of CABLE var every step
           ! returned to JULES from cable_hydrolog AND cable_implicit
           snow_tile,          &
           ! canopy water storage - reieved as canopy_tile in CABLE
           !unpacked from CABLE - but only for dumping?
           canopy_tile,   &
           ! CO2 mass mixing ratio INTENT(IN) 
           CO2_MMR, &
           !
           !
           ! JULES forcing  
           sw_down, &
           lw_down,   &
           ls_rain, ls_snow, &
           cos_zenith_angle, & ! constructed for CABLE
           tl_1, qw_1, vshr_land, pstar, z1_tq,&
           z1_uv, &
           !
           ! End - vars native to JULES/UM 
           !
           !from IMPLICIT CALL
           dtl_1, dqw_1, &
           TSOIL, SMCL, STHF,  &
           FTL_1,FQW_1,  &
           SURF_HT_FLUX_LAND, &
           ECAN_TILE,ESOIL_TILE,EI_TILE,&
           GS, T1P5M_TILE, Q1P5M_TILE, &
           CANOPY_GB, MELT_TILE, DIM_CS1,DIM_CS2, NPP, NPP_FT, &
           GPP, GPP_FT, RESP_S, &
           RESP_S_TOT, RESP_P, RESP_P_FT,  &
           G_LEAF, &
           !
           ! new CABLE vars
           !
           ! CABLE_vars initialized from ancillaries
           !
           ! snow vars - passed as read 
           isnow_flg3l, snow_rho1l, snage_tile,          &
           ! snow vars - read as separate var per layer and pre-packed 
           snow_rho3l, snow_depth3l, snow_tmp3l, snow_mass3l, &
           ! soil vars - read as separate var per layer and pre-packed
           sthu_tile, smcl_tile, sthf_tile, tsoil_tile, &
           ! snow_cond requires no init from file?
           snow_cond, &
           FTL_TILE,  &
           FQW_TILE, TSTAR_TILE,   &
           U_S, U_S_STD_TILE,&
           RADNET_TILE, FRACA, rESFS, RESFT, Z0H_TILE,  &
           Z0M_TILE, EPOT_TILE & 
                                  )    

   
   !--- reads runtime and user switches and reports
   USE cable_um_tech_mod, ONLY : cable_um_runtime_vars
   USE cable_um_tech_mod, ONLY : air, bgc, canopy,    &
                                 met, bal, rad, rough, soil, ssnow, sum_flux, veg 
   
   !--- vars common to CABLE declared 
   USE cable_common_module, ONLY : cable_runtime, cable_user, ktau_gl,          &
                                   knode_gl, kwidth_gl, kend_gl
   
   !--- subr to (manage)interface UM data to CABLE
   USE cable_um_init_mod!, ONLY : interface_UM_data
   
   !--- subr to call CABLE model
   USE cable_cbm_module, ONLY : cbm

!!   USE cable_def_types_mod, ONLY : mp, ms, ssnow, rough, canopy, air, rad,     &
!!                                   met

   IMPLICIT NONE
 
 
  
   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM sf_exch() --------------------------------------------
   !-------------------------------------------------------------------------- 
   !___IN: UM dimensions, array indexes, flags
   INTEGER, INTENT(IN) ::                                                      & 
      row_length, rows, & ! UM grid resolution
      land_pts,         & ! # of land points being processed
      ntiles,           & ! # of tiles 
      npft,             & ! # of plant functional types
      sm_levels           ! # of soil layers 

   ! end step of experiment, this step, step width, processor num
   INTEGER, INTENT(IN) :: timestep, endstep, timestep_number, mype

   ! index of land points being processed
   INTEGER, INTENT(IN), DIMENSION(land_pts) :: land_index 

   REAL, INTENT(IN), DIMENSION(row_length,rows) ::                             &
      latitude,   &
      longitude

   REAL, INTENT(IN), DIMENSION(land_pts, ntiles) ::                            &
      tile_frac

   ! # of land points on each tile
   INTEGER, INTENT(IN), DIMENSION(ntiles) :: tile_pts 
   
   INTEGER, INTENT(IN), DIMENSION(land_pts, ntiles) ::                         & 
      tile_index   ! index of tile points being processed

   !___UM parameters: soil layer thicknesses 
   REAL, INTENT(IN), DIMENSION(sm_levels) :: dzsoil

   !___UM soil/snow/radiation/met vars
   REAL, INTENT(IN), DIMENSION(land_pts) :: & 
      bexp,    & ! => parameter b in Campbell equation 
      hcon,    & ! Soil thermal conductivity (W/m/K).
      satcon,  & ! hydraulic conductivity @ saturation [mm/s]
      sathh,   &
      smvcst,  &
      smvcwt,  &
      smvccl,  &
      albsoil

   REAL, INTENT(IN), DIMENSION(land_pts) :: & 
      fland 
   
   REAL, INTENT(INOUT), DIMENSION(row_length,rows) :: &
      sw_down,          & 
      cos_zenith_angle
   
   REAL, INTENT(IN), DIMENSION(row_length,rows) ::                             &
      lw_down,    &
      ls_rain,    &
      ls_snow,    &
      tl_1,       &
      qw_1,       &  
      vshr_land,  &
      pstar,      &
      z1_tq,      &
      z1_uv

   REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles) ::                         &
      snow_tile

   REAL, INTENT(IN), DIMENSION(land_pts, npft) ::                              &
      canht_ft, lai_ft 
   
   REAL, INTENT(IN),DIMENSION(land_pts, ntiles) ::                             &
      canopy_tile
 
    INTEGER, INTENT(IN), DIMENSION(land_pts, ntiles) ::                        & 
      isnow_flg3l   ! 3 layer snow flag

   REAL, INTENT(IN), DIMENSION(land_pts, ntiles) ::                            &
      snow_rho1l, &
      snage_tile
   
   REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles,3) ::                       &
      snow_cond
   
   REAL, INTENT(IN), DIMENSION(land_pts, ntiles,3) ::                          &
      snow_rho3l,    &
      snow_depth3l,  &
      snow_mass3l,   &
      snow_tmp3l
   
   REAL, INTENT(IN), DIMENSION(land_pts, sm_levels) ::                         &
      sthu 
   
   REAL, INTENT(IN), DIMENSION(land_pts, ntiles, sm_levels) :: & 
      sthu_tile, &
      sthf_tile, &
      smcl_tile, &
      tsoil_tile
   
   REAL, INTENT(IN) :: co2_mmr

   !___return fluxes

   REAL, DIMENSION(land_pts,ntiles) :: &
      FTL_TILE,   &  ! Surface FTL for land tiles     
      FQW_TILE       ! Surface FQW for land tiles     

   !___return temp and roughness
   REAL, DIMENSION(land_pts,ntiles) :: &
      TSTAR_TILE,       & 
      Z0H_TILE,         &
      Z0M_TILE

   !___return friction velocities/drags/ etc
   REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
      U_S_STD_TILE      ! Surface friction velocity

   REAL, INTENT(OUT), DIMENSION(row_length,rows)  :: &
      U_S               ! Surface friction velocity (m/s)
   
   !___return miscelaneous 
   REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
      RADNET_TILE,   & ! Surface net radiation
      RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
                       ! factor for fraction (1-FRACA) of snow-free land tiles
      RESFT,         & ! Total resistance factor.
                       ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
                       ! 1 for snow.    
      FRACA,         & ! Fraction of surface moisture
      EPOT_TILE
     
   !-------------------------------------------------------------------------- 
   !--- end INPUT ARGS FROM sf_exch() ----------------------------------------
   !-------------------------------------------------------------------------- 
   

   !From IMPLICIT call
   REAL, DIMENSION(ROW_LENGTH,ROWS) ::                                 &
      DTL_1,    & ! IN Level 1 increment to T field 
      DQW_1       ! IN Level 1 increment to q field 

   INTEGER ::                                                                  &
      DIM_CS1, DIM_CS2 

   REAL, DIMENSION(land_pts) ::                                            &
      GS        ! OUT "Stomatal" conductance to
   
   REAL, DIMENSION(ROW_LENGTH,ROWS) ::                                 &
      !--- Net downward heat flux at surface over land.
      !--- fraction of gridbox (W/m2).
      SURF_HT_FLUX_LAND,                                                       &   
      !--- Moisture flux between layers. (kg/m^2/sec).
      !--- FQW(,1) is total water flux from surface, 'E'.
      FQW_1,                                                                   &  
      !--- FTL(,K) =net turbulent sensible heat flux into layer K
      !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
      FTL_1         

   REAL, DIMENSION(land_pts,ntiles) ::                              &
      !___(tiled) stomatatal conductance
     GS_TILE,                                           &
     !___ INOUT Surface net radiation on tiles (W/m2)
     EI_TILE,     & ! OUT EI for land tiles.
     ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
     ESOIL_TILE     ! evapotranspiration from soil moisture store (kg/m2/s) 

   !___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
   !___ runoff ??
   REAL, dimension(land_pts,sm_levels) ::                           &
      SMCL,       & ! 
      STHF,       & !
      TSOIL         !

   REAL, DIMENSION(land_pts,ntiles) ::                              &
      MELT_TILE,   & 
      FRS_TILE,   & ! Local
      NEE_TILE,   & ! Local
      NPP_TILE,   & ! Local
      GPP_TILE,   & ! Local
      GLEAF_TILE, & ! Local, kdcorbin, 10/10
      FRP_TILE,   &
      NPP_FT,     &
      NPP_FT_old, &
      GPP_FT,     &
      GPP_FT_old       

   REAL, DIMENSION(land_pts) ::                                         &
      CANOPY_GB,   & !
      RESP_P,      & !
      NPP,         & !
      GPP            !
      
   REAL, DIMENSION( land_pts,ntiles ) ::                               &
      T1P5M_TILE,    &
      Q1P5M_TILE,    &
      RESP_S_TILE,   & 
      RESP_P_FT,     &
      RESP_P_FT_old, &
      G_LEAF

   REAL ::                                                                     &
      RESP_S(LAND_PTS,DIM_CS1),     &
      RESP_S_old(LAND_PTS,DIM_CS1), &
      RESP_S_TOT(DIM_CS2)    
     

   
   !___ declare local vars 
   
   REAL, DIMENSION(ROW_LENGTH,ROWS) ::                                 &
      CONV_RAIN, & ! IN Convective rain
      CONV_SNOW   ! IN Convective snow
   
   REAL, DIMENSION(land_pts,ntiles) :: &
      CD_TILE,    &     ! Drag coefficient
      CH_TILE,    &     ! Transfer coefficient for heat & moisture
      RECIP_L_MO_TILE   ! Reciprocal of the Monin-Obukhov length for tiles (m^-1)

   !___ location of namelist file defining runtime vars
   CHARACTER(LEN=200), PARAMETER ::                                            & 
      runtime_vars_file = 'cable.nml'


   !___ 1st call in RUN (!=ktau_gl -see below) 
   LOGICAL, SAVE :: first_call = .TRUE.
   

!jhan: use after defined cable%
   !--- basic info from global model passed to cable_common_module 
   !--- vars so don't need to be passed around, just USE _module
   ktau_gl = timestep_number     !timestep of EXPERIMENT not necesarily 
                                 !the same as timestep of particular RUN
   knode_gl = mype               !which processor am i on?
   kwidth_gl = timestep          !width of timestep (secs)
   kend_gl = endstep             !timestep of EXPERIMENT not necesarily 

   !--- user FLAGS, variables etc def. in cable.nml is read on 
   !--- first time step of each run. these variables are read at 
   !--- runtime and for the most part do not require a model rebuild.
   IF(first_call) THEN
      CALL cable_um_runtime_vars(runtime_vars_file) 
      first_call = .FALSE.
   ENDIF      

   !---------------------------------------------------------------------!
   !--- initialize CABLE using UM forcings etc. these args are passed ---!
   !--- down from ~UM. explicit_call                                  ---! 
   !---------------------------------------------------------------------!
   CALL explicit_call_initialization(                                          & 
            ! cable% type yet to be officailly implemented
            row_length, & ! -> cable%row_length 
            rows,       & ! -> cable%rows
            land_pts,   & ! -> cable%land_pts
            ntiles,     & ! -> cable%ntiles
            npft,       & ! -> cable%npft 
            sm_levels,  & ! -> cable%ms
            timestep,   & ! -> cable%timestep_width
            latitude,   & ! -> cable%latitude
            longitude,  & ! -> cable%longitude
            land_index, & ! -- necessary for packing 
            tile_frac,  & ! -> cable%tile_frac
            tile_pts,   & ! -> cable%
            tile_index, & ! -- necessary for packing
            dzsoil,     & ! -> soil%zse                        
            ! soil properties from UM/JULES
            bexp,       & ! -> soil%bch
            hcon,       & ! ~> soil%cnsd
            satcon,     & ! ~> soil%hyds
            sathh,      & ! -> soil%sucs
            smvcst,     & ! -> soil%ssat
            smvcwt,     & ! -> soil%swilt
            smvccl,     & ! -> soil%sfc
            albsoil,    & ! -> soil%albsoil
            ! canopy properties from UM/JULES
            canht_ft,   & ! ~> veg%hc
            lai_ft,     & ! ~> veg%lai
            ! forcing from JULES
            sw_down,    & ! ~> met%fsd
            lw_down,    & ! -> met%fld 
            ls_rain,    & ! ~> met%precip
            ls_snow,    & ! ~> met%precip_sn
            tl_1,       & ! -> met%tk
            qw_1,       & ! -> met%qv
            vshr_land,  & ! -> met%ua
            pstar,      & ! ~> met%pmb
            z1_tq,      & ! -> rough%za_tq
            z1_uv,      & ! -> rough%za_uv
            canopy_tile,& ! -> canopy%cansto
            Fland,      & ! -> ssnow%fland
            CO2_MMR,    & ! ~> met%ca
         
            !jhan:adapted from JULES var cosz, done elsewhere, move to here 
            ! and make switchable
            cos_zenith_angle,    & ! ->met%coszen
         
            ! snow properties from UM/JULES
            snow_tile,     & ! -> ssnow%snowd
         
            ! snow properties from CABLE vars 
            snage_tile,    & ! -> ssnow%snage
            snow_rho1l,    & ! -> ssnow%ssdnn
            isnow_flg3l,   & ! -> ssnow%isflag
            snow_rho3l,    & ! -> ssnow%ssdn
            snow_depth3l,  & ! -> ssnow%sdepth
            snow_tmp3l,    & ! -> ssnow%tggsn
            snow_mass3l,   & ! -> ssnow%smass
            snow_cond,     & ! -> ssnow%sconds
         
         
            !soil properties from CABLE vars 
            smcl_tile,     & ! ~> soil%wb
            sthf_tile,     & ! ~> soil%wbice
            tsoil_tile     & ! -> ssnow%tgg
   )                         

   CALL implicit_call_initialization(                                          & 
            row_length, & ! -> cable%row_length 
            rows,       & ! -> cable%rows
            ls_rain,    & ! ~>  met%precip
            ls_snow,    & ! ~>  met%precip_sn
            conv_rain,  & ! ~~> met%precip
            conv_snow,  & ! ~~> met%precip_sn
            dtl_1,      & ! ~~> met%precip_sn
            dqw_1       & ! ~~> met%precip_sn
   )
    


   !---------------------------------------------------------------------!
   !--- real(timestep) width, CABLE types passed to CABLE "engine" as ---!  
   !--- req'd by Mk3L  --------------------------------------------------!
   !---------------------------------------------------------------------!
   CALL cbm( REAL(timestep), air, bgc, canopy, met, bal,                       &
             rad, rough, soil, ssnow, sum_flux, veg )


   CALL cable_expl_unpack( FTL_TILE, FQW_TILE,       &
                           TSTAR_TILE, &
                           U_S, U_S_STD_TILE, &
                           CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           RECIP_L_MO_TILE, EPOT_TILE, &
                           ssnow%snowd, ssnow%cls, air%rlam, air%rho,          &
                           canopy%fe, canopy%fh, canopy%us, canopy%cdtq,       &
                           canopy%fwet, canopy%wetfac_cs, canopy%rnet,         &
                           canopy%zetar, canopy%epot, met%ua, rad%trad,        &
                           rad%transd, rough%z0m, rough%zref_tq )



   CALL implicit_unpack( TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,                &
                            SMVCST, STHF, STHF_TILE, STHU, STHU_TILE,          &
                            snow_tile, SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L,  &
                            SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,    &
                            FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
                            SURF_HT_FLUX_LAND, ECAN_TILE,        &
                            ESOIL_TILE, EI_TILE, RADNET_TILE, &
                            SNAGE_TILE, CANOPY_TILE, GS, T1P5M_TILE,           &
                            Q1P5M_TILE, CANOPY_GB, FLAND, MELT_TILE, DIM_CS1,  &
                            DIM_CS2, NPP, NPP_FT, GPP, GPP_FT, RESP_S,         &
                            RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT, G_LEAF )
 


END SUBROUTINE cable_main



!---------------------------------------------------------------------!
!--- pass land-surface quantities calc'd by CABLE in explicit call ---!
!--- back to UM.                                                   ---!
!---------------------------------------------------------------------!
SUBROUTINE cable_expl_unpack( FTL_TILE, FQW_TILE,       &
                           TSTAR_TILE, &
                           U_S, U_S_STD_TILE, &
                           CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           RECIP_L_MO_TILE, EPOT_TILE, &
                           ssnow_snowd, ssnow_cls, air_rlam, air_rho,          &
                           canopy_fe, canopy_fh, canopy_us, canopy_cdtq,       &
                           canopy_fwet, canopy_wetfac_cs, canopy_rnet,         &
                           canopy_zetar, canopy_epot, met_ua, rad_trad,        &
                           rad_transd, rough_z0m, rough_zref_tq )

   USE cable_def_types_mod, ONLY : mp, NITER 
   USE cable_data_module,   ONLY : PHYS
   USE cable_um_tech_mod,   ONLY : um1
   USE cable_common_module, ONLY : cable_runtime, cable_user, &
                                   ktau_gl, knode_gl 
   IMPLICIT NONE         


   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
   !-------------------------------------------------------------------------- 


   !___ UM variables to recieve unpacked CABLE vars

   !___return fluxes
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      FQW_TILE       ! Surface FQW for land tiles     

   !___return temp and roughness
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      TSTAR_TILE,  Z0H_TILE, Z0M_TILE

   !___return friction velocities/drags/ etc
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      CD_TILE,    &     ! Drag coefficient
      CH_TILE,    &     ! Transfer coefficient for heat & moisture
      U_S_STD_TILE      ! Surface friction velocity
   REAL, INTENT(OUT), DIMENSION(um1%row_length,um1%rows)  :: &
      U_S               ! Surface friction velocity (m/s)

   !___return miscelaneous 
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      RADNET_TILE,   & ! Surface net radiation
      RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
                       ! factor for fraction (1-FRACA) of snow-free land tiles
      RESFT,         & ! Total resistance factor.
                       ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,
                       ! 1 for snow.    
      FRACA,         & ! Fraction of surface moisture
      RECIP_L_MO_TILE,&! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
      EPOT_TILE
   
   LOGICAL,DIMENSION(um1%land_pts,um1%ntiles) :: l_tile_pts

   !___UM vars used but NOT returned 
   REAL, INTENT(IN), DIMENSION(um1%land_pts) ::   &
      FLAND(um1%land_pts)              ! IN Land fraction on land tiles.




   !___ decs of intent(in) CABLE variables to be unpacked

   ! snow depth (liquid water), factor for latent heat
   REAL, INTENT(IN), DIMENSION(mp) :: ssnow_snowd, ssnow_cls
   
   ! surface wind speed (m/s)
   REAL, INTENT(IN), DIMENSION(mp) :: met_ua 
   
   ! latent heat for water (j/kg), dry air density (kg m-3)
   REAL, INTENT(IN), DIMENSION(mp) :: air_rlam, air_rho 
   
   ! frac SW diffuse transmitted thru canopy, rad. temp. (soil and veg)
   REAL, INTENT(IN), DIMENSION(mp) :: rad_trad,rad_transd 
   
   ! total latent heat (W/m2), total sensible heat (W/m2)
   REAL, INTENT(IN), DIMENSION(mp) :: canopy_fe, canopy_fh  
   
   ! fraction of canopy wet
   REAL, INTENT(IN), DIMENSION(mp) :: canopy_fwet, canopy_wetfac_cs
   
   ! friction velocity, drag coefficient for momentum
   REAL, INTENT(IN), DIMENSION(mp) :: canopy_us, canopy_cdtq
   
   ! net rad. absorbed by surface (W/m2), total potential evaporation 
   REAL, INTENT(IN), DIMENSION(mp) :: canopy_rnet, canopy_epot        
   
   ! stability correction
   REAL, INTENT(IN), DIMENSION(mp,niter) :: canopy_zetar
   
   ! roughness length, Reference height for met forcing
   REAL, INTENT(IN), DIMENSION(mp) :: rough_z0m, rough_zref_tq 
 
   !-------------------------------------------------------------------------- 
   !--- end INPUT ARGS FROM cable_explicit_driver() ------------------------------
   !-------------------------------------------------------------------------- 


   
        
   !___vars in local calc. of latent heat fluxes
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
      FQW_TILE_CAB,  &
      LE_TILE

   REAL, DIMENSION(um1%land_pts) ::   &
      FTL_CAB, &
      LE_CAB

   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: &
      FTL_TILE_CAB, &
      FTL_TILE,   &  ! Surface FTL for land tiles     
      LE_TILE_CAB
   
   !___vars in local calc of Surface friction velocities
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
      CD_CAB_TILE,   &  
      CH_CAB_TILE,   &  ! (bulk transfer) coeff. for momentum
      U_S_TILE

   REAL, DIMENSION(um1%land_pts) ::                  &
      CH_CAB,  &  ! Turbulent surface exchange
      CD_CAB,  &  ! Turbulent surface exchange
      U_S_CAB     ! Surface friction velocity (m/s)

   REAL, DIMENSION(mp)  :: &
      CDCAB,CHCAB

   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: &
      TSTAR_TILE_CAB
   
   REAL, DIMENSION(um1%land_pts) ::                  &
      TSTAR_CAB
   

   !___local miscelaneous
   REAL, DIMENSION(mp)  :: &
   THETAST,fraca_cab,rfsfs_cab, RECIPLMOTILE, fe_dlh
   INTEGER :: i,j,k,N,L
   REAL :: miss = 0.0
   LOGICAL, SAVE :: first_call = .true.
   REAL, POINTER :: CAPP 
   
      CAPP => PHYS%CAPP
      
      !___return fluxes
      FTL_TILE_CAB = UNPACK(canopy_fh,  um1%l_tile_pts, miss)
      FTL_CAB = SUM(um1%TILE_FRAC * FTL_TILE_CAB,2)
      FQW_TILE_CAB = UNPACK(canopy_fe,  um1%l_tile_pts, miss)
      LE_TILE_CAB = UNPACK(canopy_fe,  um1%l_tile_pts, miss)
      LE_CAB = SUM(um1%TILE_FRAC * LE_TILE_CAB,2)
      fe_dlh = canopy_fe/(air_rlam*ssnow_cls)
      FTL_TILE = UNPACK(canopy_fh,  um1%l_tile_pts, miss)
      FTL_TILE = FTL_TILE / capp
      FQW_TILE = UNPACK(fe_dlh, um1%l_tile_pts, miss)

      !___return temp and roughness
      TSTAR_TILE_CAB = UNPACK(rad_trad, um1%l_tile_pts, miss)
      TSTAR_CAB = SUM(um1%TILE_FRAC * TSTAR_TILE_CAB,2)
      TSTAR_TILE = UNPACK(rad_trad,  um1%l_tile_pts, miss)
      Z0M_TILE = UNPACK(rough_z0m,  um1%l_tile_pts, miss)
      Z0H_TILE = Z0M_TILE
      
      !___return friction velocities/drags/ etc
      U_S_TILE  =  UNPACK(canopy_us, um1%l_tile_pts, miss)
      U_S_CAB  = SUM(um1%TILE_FRAC *  U_S_TILE,2)
      CDCAB = canopy_us**2/met_ua**2   ! met%ua is always above umin = 0.1m/s
      ! for Cable CD*
      CD_CAB_TILE =  UNPACK(CDCAB,um1%l_tile_pts, miss)
      CD_CAB= SUM(um1%TILE_FRAC * CD_CAB_TILE,2)
      ! for Cable CH*
      CH_CAB_TILE =  UNPACK(canopy_cdtq,um1%l_tile_pts, miss)
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
      WHERE( ssnow_snowd > 1.0 ) fraca_cab = 1.0
      rfsfs_cab = MIN( 1., MAX( 0.01, canopy_wetfac_cs - fraca_cab ) /         &
                  MAX( 0.01,1. - fraca_cab ) )
      FRACA = UNPACK( fraca_cab, um1%l_tile_pts, miss )
      RESFT = UNPACK( canopy_wetfac_cs,um1%l_tile_pts, miss )
      RESFS = UNPACK( rfsfs_cab , um1%l_tile_pts, miss )

      RADNET_TILE = UNPACK( canopy_rnet , um1%l_tile_pts, miss )
      THETAST = ABS( canopy_fh ) / ( air_rho * capp*canopy_us )
      RECIPLMOTILE =  canopy_zetar(:,niter) / rough_zref_tq
      RECIP_L_MO_TILE = UNPACK( RECIPLMOTILE, um1%l_tile_pts, miss )
      EPOT_TILE = UNPACK( canopy_epot, um1%l_tile_pts, miss )
      

      IF(first_call) THEN 
         l_tile_pts = um1%l_tile_pts
         first_call = .FALSE.
      ENDIF

   
END SUBROUTINE cable_expl_unpack


!========================================================================= 
!========================================================================= 
!========================================================================= 
        
SUBROUTINE implicit_unpack( TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,                &
                            SMVCST, STHF, STHF_TILE, STHU, STHU_TILE,          &
                            snow_tile, SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L,  &
                            SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,    &
                            FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
                            SURF_HT_FLUX_LAND, ECAN_TILE,        &
                            ESOIL_TILE, EI_TILE, RADNET_TILE, &
                            SNAGE_TILE, CANOPY_TILE, GS, T1P5M_TILE,           &
                            Q1P5M_TILE, CANOPY_GB, FLAND, MELT_TILE, DIM_CS1,  &
                            DIM_CS2, NPP, NPP_FT, GPP, GPP_FT, RESP_S,         &
                            RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT, G_LEAF )
 
   USE cable_def_types_mod, ONLY : mp
   USE cable_data_module,   ONLY : PHYS
   USE cable_um_tech_mod,   ONLY : um1 ,canopy, rad, soil, ssnow, air
   USE cable_common_module, ONLY : cable_runtime, cable_user
   IMPLICIT NONE
 
   !jhan:these need to be cleaned out to what is actualllly passed
   INTEGER :: DIM_CS1 ,DIM_CS2 

   REAL, DIMENSION(um1%land_pts) ::                                            &
      GS,         &  ! OUT "Stomatal" conductance to
      SMVCST,     &  ! IN Volumetric saturation point
      FLAND          ! IN Land fraction on land tiles
   
   real, dimension(um1%ROW_LENGTH,um1%ROWS) ::                                 &
      !--- Net downward heat flux at surface over land.
      !--- fraction of gridbox (W/m2).
      SURF_HT_FLUX_LAND,           &
      !--- Moisture flux between layers. (kg/m^2/sec).
      !--- FQW(,1) is total water flux from surface, 'E'.
      FQW_1,       &  
      !--- FTL(,K) =net turbulent sensible heat flux into layer K
      !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
      FTL_1         

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
      !___Surface FTL, FQL for land tiles
      FTL_TILE, FQW_TILE, FQW_TILE_CAB,   &  
      !___(tiled) latent heat flux, melting, stomatatal conductance
      LE_TILE, MELT_TILE, GS_TILE,     &  
      RADNET_TILE, & ! INOUT Surface net radiation on tiles (W/m2)
      EI_TILE,     & ! OUT EI for land tiles.
      ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
      ESOIL_TILE     ! evapotranspiration from soil moisture store (kg/m2/s) 

   !___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
   !___ runoff ??
   REAL, DIMENSION(um1%land_pts,um1%sm_levels) ::                              &
      SMCL,       & !
      STHF,       &
      STHU,       &
      SMCL_CAB,   &
      TSOIL_CAB,  &
      TSOIL,      &
      SURF_CAB_ROFF !      

   !___(tiled) soil prognostics: as above 
   REAL, DIMENSION(um1%land_pts,um1%ntiles,um1%sm_levels) ::                   &
      SMCL_TILE,  & 
      STHU_TILE,  &
      TSOIL_TILE, &
      STHF_TILE  

   !___flag for 3 layer snow pack
   INTEGER :: ISNOW_FLG3L(um1%LAND_PTS,um1%NTILES)
   
   !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
   REAL, DIMENSION(um1%land_pts,um1%ntiles,3) ::                               &
      SNOW_DEPTH3L,  &
      SNOW_MASS3L,   &
      SNOW_RHO3L,    &
      SNOW_TMP3L,    &
      SNOW_COND 

   REAL, dimension(um1%land_pts,um1%ntiles) ::                                 &
      FRS_TILE,   & ! Local
      NEE_TILE,   & ! Local
      NPP_TILE,   & ! Local
      GPP_TILE,   & ! Local
      GLEAF_TILE, & ! Local, kdcorbin, 10/10
      FRP_TILE,   &
      NPP_FT,     &
      NPP_FT_old, &
      GPP_FT,     &
      GPP_FT_old       

   REAL, dimension(um1%land_pts) ::                                            &
      SNOW_GRD,   &  
      CANOPY_GB,  &
      FTL_CAB,    &
      LE_CAB,     &
      TSTAR_CAB,  &
      SURF_HTF_CAB,&
      RESP_P,     & 
      NPP,        & 
      GPP
   
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
      SNOW_TILE,     & !
      SNOW_RHO1L,    & ! Mean snow density
      SNAGE_TILE,    & !
      CANOPY_TILE,   & !
      FTL_TILE_CAB,  & !
      LE_TILE_CAB,   & !
      T1P5M_TILE,    &
      Q1P5M_TILE,    &
      TSTAR_TILE_CAB,&
      TSTAR_TILE,    &
      SURF_HTF_T_CAB,& 
      RESP_S_TILE,   & 
      RESP_P_FT,     &
      RESP_P_FT_old, &
      G_LEAF

   REAL ::                                                                     &
      RESP_S(um1%LAND_PTS,DIM_CS1),    & !
      RESP_S_old(um1%LAND_PTS,DIM_CS1),& !
      RESP_S_TOT(DIM_CS2)                !
  
   REAL, DIMENSION(mp) ::                                                                     &
      fe_dlh,    & !
      fes_dlh,   & !
      fev_dlh      !

   !--- Local vars
   INTEGER :: i,j,l,k,n

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
         !--- Local buffer surface FTL, FQL @ prev dt
         FTL_TILE_old, FQW_TILE_old

   INTEGER:: i_miss = 0
   REAL :: miss = 0.0
   
   REAL, POINTER :: TFRZ
   
   LOGICAL, SAVE :: first_call = .TRUE. 
   
      TFRZ => PHYS%TFRZ
  
      !--- set UM vars to zero
      TSOIL_CAB = 0.; SMCL_CAB = 0.; TSOIL_TILE = 0.; 
      SMCL_TILE = 0.; STHF_TILE = 0.; STHU_TILE = 0.

      DO j = 1,um1%SM_LEVELS
         TSOIL_TILE(:,:,j)= UNPACK(ssnow%tgg(:,j), um1%L_TILE_PTS, miss)
         TSOIL_CAB(:,j) = SUM(um1%TILE_FRAC * TSOIL_TILE(:,:,j),2)
         SMCL_TILE(:,:,j)= UNPACK(REAL(ssnow%wb(:,j)), um1%L_TILE_PTS, miss)
         SMCL_TILE(:,:,j)=SMCL_TILE(:,:,j)*soil%zse(j)*PHYS%RHOW
         SMCL_CAB(:,j) = SUM(um1%TILE_FRAC * SMCL_TILE(:,:,j),2)
         STHF_TILE(:,:,j)= UNPACK(REAL(ssnow%wbice(:,j)), um1%L_TILE_PTS, miss)
         SMCL(:,j) = SUM(um1%TILE_FRAC * SMCL_TILE(:,:,j),2)
         TSOIL(:,j) = SUM(um1%TILE_FRAC * TSOIL_TILE(:,:,j),2)
         
         DO N=1,um1%NTILES
            DO K=1,um1%TILE_PTS(N)
               I = um1%TILE_INDEX(K,N)
               IF ( SMVCST(I) > 0. ) THEN ! Exclude permanent ice - mrd
                  STHF_TILE(I,N,J)= STHF_TILE(I,N,J)/SMVCST(I)
                  STHU_TILE(I,N,J)= MAX( 0., SMCL_TILE(I,N,J) -                &
                                    STHF_TILE(I,N,J) * SMVCST(I) * soil%zse(J) &
                                    * PHYS%RHOW ) / ( soil%zse(J) *        &
                                    PHYS%RHOW * SMVCST(I) )
               ENDIF
            ENDDO
         ENDDO

         STHF(:,J) = SUM(um1%TILE_FRAC * STHF_TILE(:,:,J),2)
         STHU(:,J) = SUM(um1%TILE_FRAC * STHU_TILE(:,:,J),2)
      ENDDO

!IF(first_call) THEN 
!   print *, " "
!   print *, "jhan:impl::tgg", ssnow%tgg
!   print *, " "
!   print *, "jhan:impl::tsoil_tile ", tsoil_tile
!   print *, " "
!   print *, "jhan:impl::tsoil", tsoil
!ENDIF

!IF(first_call) THEN 
!   print *, 'tsoil', tsoil, tsoil_tile
!   first_call=.FALSE.
!ENDIF

      !--- unpack snow vars 
      SNOW_RHO1L  = UNPACK(ssnow%ssdnn, um1%L_TILE_PTS, miss)
      ISNOW_FLG3L = UNPACK(ssnow%isflag, um1%L_TILE_PTS, i_miss)
      MELT_TILE   = UNPACK(ssnow%smelt, um1%L_TILE_PTS, miss)
      SNOW_TILE= UNPACK(ssnow%snowd, um1%L_TILE_PTS, miss)
      SNOW_GRD=  SUM(um1%TILE_FRAC * SNOW_TILE,2)  ! gridbox snow mass & snow below canopy 

      !--- unpack layered snow vars 
      do k = 1,3
        SNOW_TMP3L(:,:,k) = UNPACK(ssnow%tggsn(:,k), um1%L_TILE_PTS, miss)
        SNOW_MASS3L(:,:,k)= UNPACK(ssnow%smass(:,k), um1%L_TILE_PTS, miss)
        SNOW_RHO3L(:,:,k) = UNPACK(ssnow%ssdn(:,k), um1%L_TILE_PTS, miss)
        SNOW_COND(:,:,k)  = UNPACK(ssnow%sconds(:,k),um1%L_TILE_PTS,miss)
        SNOW_DEPTH3L(:,:,k)  = UNPACK(ssnow%sdepth(:,k),um1%L_TILE_PTS,miss)
      enddo

      !---???
      GS_TILE = UNPACK(canopy%gswx_T,um1%L_TILE_PTS,miss)
      GS =  SUM(um1%TILE_FRAC * GS_TILE,2)

      !___return fluxes
      FTL_TILE_CAB = UNPACK(canopy%fh,  um1%l_tile_pts, miss)
      FTL_CAB = SUM(um1%TILE_FRAC * FTL_TILE_CAB,2)
      FQW_TILE_CAB = UNPACK(canopy%fe,  um1%l_tile_pts, miss)
      LE_TILE_CAB = FQW_TILE_CAB 
      LE_CAB = SUM(um1%TILE_FRAC * LE_TILE_CAB,2)
      fe_dlh = canopy%fe/(air%rlam*ssnow%cls)
      where ( fe_dlh .ge. 0.0 ) fe_dlh = MAX ( 1.e-6, fe_dlh )
      where ( fe_dlh .lt. 0.0 ) fe_dlh = MIN ( -1.e-6, fe_dlh )
      fes_dlh = canopy%fes/(air%rlam*ssnow%cls)
      fev_dlh = canopy%fev/air%rlam

      !---preserve fluxes from the previous time step for the coastal grids
      FTL_TILE_old = FTL_TILE
      FQW_TILE_old = FQW_TILE
      
      !---update fluxes 
      FTL_TILE = FTL_TILE_CAB 
      FQW_TILE = UNPACK(fe_dlh, um1%l_tile_pts, miss)

      !___return temp and roughness
      TSTAR_TILE_CAB = UNPACK(rad%trad, um1%l_tile_pts, miss)
      TSTAR_CAB = SUM(um1%TILE_FRAC * TSTAR_TILE_CAB,2)
      TSTAR_TILE = TSTAR_TILE_CAB 

      !___return miscelaneous 
      RADNET_TILE = unpack( canopy%rnet , um1%l_tile_pts, miss)

      SURF_HTF_T_CAB = UNPACK(canopy%ga,um1%L_TILE_PTS,miss)
      SURF_HTF_CAB = SUM(um1%TILE_FRAC * SURF_HTF_T_CAB,2)

     ESOIL_TILE = UNPACK(fes_dlh, um1%L_TILE_PTS, miss)
     ECAN_TILE = UNPACK(fev_dlh,  um1%L_TILE_PTS, miss)
     EI_TILE = 0.
     SNAGE_TILE = UNPACK(ssnow%snage, um1%L_TILE_PTS, miss) 

     !unpack screen level (1.5m) variables
     !Convert back to K 
     t1p5m_tile     = UNPACK(canopy%tscrn+tfrz, um1%L_TILE_PTS, miss)
     q1p5m_tile     = UNPACK(canopy%qscrn, um1%L_TILE_PTS, miss)
     CANOPY_TILE    = UNPACK(canopy%cansto, um1%L_TILE_PTS, miss)
     CANOPY_GB      = SUM(um1%TILE_FRAC * CANOPY_TILE,2)

     ! Lestevens - Passing CO2 from CABLE to bl_trmix_dd.F90
     FRS_TILE       = UNPACK(canopy%frs, um1%L_TILE_PTS, miss)
     NEE_TILE       = UNPACK(canopy%fnee, um1%L_TILE_PTS, miss)
     NPP_TILE       = UNPACK(canopy%fnpp, um1%L_TILE_PTS, miss)
     GLEAF_TILE     = UNPACK(canopy%frday,um1%L_TILE_PTS, miss)

      IF( cable_user%leaf_respiration == 'on' .OR.                             &
           cable_user%leaf_respiration == 'ON') THEN
         GPP_TILE = UNPACK(canopy%fnpp+canopy%frp, um1%L_TILE_PTS, miss)
      ELSE 
         GPP_TILE = UNPACK(canopy%fnpp+canopy%frp+canopy%frday,  &
                            um1%L_TILE_PTS, miss)
      ENDIF

     FRP_TILE       = UNPACK(canopy%frp, um1%L_TILE_PTS, miss)
     NPP_FT_old     = NPP_FT
     GPP_FT_old     = GPP_FT
     RESP_P_FT_old  = RESP_P_FT
     RESP_S_old     = RESP_S

     !initialse full land grids and retain coastal grid fluxes
      DO N=1,um1%NTILES
         DO K=1,um1%TILE_PTS(N)
           L = um1%TILE_INDEX(K,N)
           J=(um1%LAND_INDEX(L)-1)/um1%ROW_LENGTH + 1
           I = um1%LAND_INDEX(L) - (J-1)*um1%ROW_LENGTH
           IF( FLAND(L) == 1.0) THEN 
             FTL_1(I,J) =  0.0
             FQW_1(I,J) =  0.0
           ELSE
             !retain sea/ice contribution and remove land contribution
             FTL_1(I,J) = FTL_1(I,J) - FLAND(L) * um1%TILE_FRAC(L,N) *         &
                          FTL_TILE_old(L,N)
             FQW_1(I,J) = FQW_1(I,J) - FLAND(L) * um1%TILE_FRAC(L,N) *         &
                          FQW_TILE_old(L,N)
           ENDIF
           SURF_HT_FLUX_LAND(I,J) = 0.
         ENDDO
     ENDDO

      DO N=1,um1%NTILES
         DO K=1,um1%TILE_PTS(N)
           L = um1%TILE_INDEX(K,N)
           J=(um1%LAND_INDEX(L)-1)/um1%ROW_LENGTH + 1
           I = um1%LAND_INDEX(L) - (J-1)*um1%ROW_LENGTH
           FTL_1(I,J) = FTL_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FTL_TILE(L,N)
           FQW_1(I,J) = FQW_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FQW_TILE(L,N)
           SURF_HT_FLUX_LAND(I,J) = SURF_HT_FLUX_LAND(I,J) +                   &
                                    FLAND(L)*um1%TILE_FRAC(L,N) *              &
                                    SURF_HTF_T_CAB(L,N)
         ENDDO
      ENDDO

      DO N=1,um1%NTILES
         DO K=1,um1%TILE_PTS(N)
            L = um1%TILE_INDEX(K,N)
            IF( FLAND(L) == 1.0) THEN
               NPP(L)=0.; NPP_FT(L,N)=0.; GPP(L)=0.; GPP_FT(L,N)=0.
               RESP_P(L)=0.; RESP_P_FT(L,N)=0.; RESP_S(L,:)=0.; G_LEAF(L,N)=0.   
            ELSE
               ! For coastal points: currently no contribution
               NPP(L)=NPP(L)-FLAND(L)*um1%TILE_FRAC(L,N)*NPP_FT_old(L,N)
               GPP(L)=GPP(L)-FLAND(L)*um1%TILE_FRAC(L,N)*GPP_FT_old(L,N)
               RESP_P(L)=RESP_P(L)-FLAND(L)*um1%TILE_FRAC(L,N)*RESP_P_FT_old(L,N)
               !--- loop for soil respiration
               DO I=1,DIM_CS1
                  RESP_S(L,I)=RESP_S(L,I)-FLAND(L)*RESP_S_old(L,I)
               ENDDO
               RESP_S_TOT(L)=sum(RESP_S(L,:))
            ENDIF
         ENDDO
      ENDDO

     RESP_S_TILE=FRS_TILE*1.e-3

      DO N=1,um1%NTILES 
         DO K=1,um1%TILE_PTS(N)
            L = um1%TILE_INDEX(K,N)
            !add leaf respiration to output
            G_LEAF(L,N)=GLEAF_TILE(L,N)*1.e-3
            NPP_FT(L,N)=NPP_TILE(L,N)*1.e-3
            NPP(L)=NPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*NPP_FT(L,N)
            GPP_FT(L,N)=GPP_TILE(L,N)*1.e-3
            GPP(L)=GPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*GPP_FT(L,N)

            !loop for soil resp. - all UM levels = single CABLE output 
            DO I=1,DIM_CS1
               RESP_S(L,I) = RESP_S(L,I) + &
                             FLAND(L)*um1%TILE_FRAC(L,N)*FRS_TILE(L,N)*1.e-3
            ENDDO

            RESP_S_TOT(L)=sum(RESP_S(L,:))
            RESP_P_FT(L,N)=FRP_TILE(L,N)*1.e-3
            RESP_P(L)=RESP_P(L)+FLAND(L)*um1%TILE_FRAC(L,N)*RESP_P_FT(L,N)
         ENDDO
      ENDDO

END SUBROUTINE Implicit_unpack









 
!
!SUBROUTINE cable_rad_driver(                                                   &
!                             ! IN atmospheric forcing
!                             surf_down_sw, cos_zenith_angle,                   &
!                             ! IN soil/vegetation/land surface data :
!                             NOW_TILE, SNOW_TMP3L, SNOW_RHO1L, TSOIL_TILE,    &
!                             ISNOW_FLG3L, ALBSOIL,                             &
!                             ! OUT
!                             LAND_ALBEDO_CABLE, ALB_TILE, LAND_ALB_CABLE ) 
!
!   USE cable_def_types_mod, ONLY : mp
!   USE cable_albedo_module, ONLY : surface_albedo
!   USE cable_um_tech_mod,   ONLY : kblum_rad, um1, soil, ssnow, rad, veg,      &
!                                   met, canopy
!   USE cable_um_init_subrs_mod, ONLY : update_kblum_radiation,  um2cable_met_rad,  &
!                                   um2cable_lp 
!   USE cable_common_module, ONLY : cable_runtime, cable_user
!   
!   IMPLICIT NONE                     
!
!   INTEGER, DIMENSION(um1%LAND_PTS,um1%NTILES) :: isnow_flg3l    
!   
!   REAL :: ALBSOIL(um1%LAND_PTS)          !          &! IN soil albedo 
!   
!   REAL, DIMENSION(um1%row_length,um1%rows) ::                                 &
!      LAND_ALB_CABLE,      & ! Land albedo calculated by Cable
!      SW_DOWN,             & ! Surface downward SW radiation (W/m2).
!      cos_zenith_angle
!
!   REAL, DIMENSION(um1%row_length,um1%rows,4) ::                               &
!      LAND_ALBEDO_CABLE, & ! Land albedo calculated by Cable NIR/VIS/Beam/Diffuse
!      surf_down_sw         ! IN Surface downward SW radiation
!
!   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                                 &
!      LAND_ALB_CABLE_TILE,  & ! Land albedo calculated by Cable
!      SNOW_TILE,            & ! IN Lying snow on tiles (kg/m2) 
!      SNOW_RHO1L             ! snow cover in the ice tile.
!   
!   REAL, DIMENSION( um1%LAND_PTS, um1%NTILES, 4 ) ::                           &
!      ALB_TILE    ! Land albedo calculated by Cable
!      
!   REAL, DIMENSION( um1%LAND_PTS, um1%NTILES, um1%SM_LEVELS ) ::               &
!      TSOIL_TILE  ! Mean snow density  (or 1 layer)
!
!   REAL, DIMENSION( um1%LAND_PTS, um1%NTILES, 3 ) ::                           &
!      SNOW_TMP3L  ! Snow temperature (3 layer)
!                                                                   
!   INTEGER :: i,J,N,K,L
!   REAL :: miss = 0.0
!   LOGICAL :: skip =.TRUE. 
!   
!   REAL :: rad_vis(mp), rad_nir(mp), met_fsd_tot_rel(mp), rad_albedo_tot(mp) 
!
!      !jhan:check that these are reset after call done
!      cable_runtime%um_radiation= .TRUE.
!      
!      !     **** surf_down_sw is from the previous time step  ****
!      !--- re-set UM rad. forcings to suit CABLE. also called in explicit call to 
!      !--- CABLE from subr cable_um_expl_update() 
!      CALL update_kblum_radiation( sw_down, cos_zenith_angle, surf_down_sw )
!   
!      !--- set met. and rad. forcings to CABLE. also called in explicit call to 
!      !--- CABLE from subr update_explicit_vars() 
!      !--- subr.  um2cable_met_rad_alb() USES CABLE types met%, rad%, soil%
!      !--- and kblum%rad. calculated in  update_kblum_radiation() above 
!      CALL um2cable_met_rad( cos_zenith_angle)
!      
!      !--- CABLE soil albedo forcing
!      CALL um2cable_lp( albsoil, albsoil, soil%albsoil(:,1), soil%isoilm, skip )
!      !-------------------------------------------------------------------
!  
!      !At present only single value is used for each land point 
!      soil%albsoil(:,2) = REAL(0) 
!      soil%albsoil(:,3) = REAL(0) 
!
!      ! soil%albsoil should be set to geograpically explicit data for 
!      ! snow free soil albedo in VIS and NIR  
!      ! get the latest surface condtions
!      ssnow%snowd  =     PACK( SNOW_TILE, um1%L_TILE_PTS )
!      ssnow%ssdnn  =     PACK( SNOW_RHO1L, um1%L_TILE_PTS )
!      ssnow%isflag =     PACK( ISNOW_FLG3L, um1%L_TILE_PTS )
!      ssnow%tggsn(:,1) = PACK( SNOW_TMP3L(:,:,1), um1%L_TILE_PTS )
!      ssnow%tgg(:,1) =   PACK( TSOIL_TILE(:,:,1), um1%L_TILE_PTS )
!
!
!      CALL surface_albedo(ssnow, veg, met, rad, soil, canopy)
!
!      ! only for land points, at present do not have a method for treating 
!      ! mixed land/sea or land/seaice points as yet.
!      ALB_TILE(:,:,1) = UNPACK(rad%reffbm(:,1),um1%L_TILE_PTS, miss)
!      ALB_TILE(:,:,2) = UNPACK(rad%reffdf(:,1),um1%L_TILE_PTS, miss)
!      ALB_TILE(:,:,3) = UNPACK(rad%reffbm(:,2),um1%L_TILE_PTS, miss)
!      ALB_TILE(:,:,4) = UNPACK(rad%reffdf(:,2),um1%L_TILE_PTS, miss)
!
!      rad_vis = ( (1.0-rad%fbeam(:,1) ) * rad%reffdf(:,1) + rad%fbeam(:,1) *   &
!                rad%reffbm(:,1) ) 
!      rad_nir = ( (1.0-rad%fbeam(:,2) ) * rad%reffdf(:,2) + rad%fbeam(:,2) *   &
!                rad%reffbm(:,2)) 
!
!      met_fsd_tot_rel = met%fsd(:,1) / MAX( 0.1, ( met%fsd(:,1)+met%fsd(:,2) ) )
!      rad_albedo_tot = met_fsd_tot_rel  * rad_vis                              &
!                       + ( 1.- met_fsd_tot_rel ) * rad_nir
!
!      LAND_ALBEDO_CABLE =0.
!      LAND_ALB_CABLE =0.
!      LAND_ALB_CABLE_TILE = UNPACK( rad_albedo_tot, um1%L_TILE_PTS, miss )
!
!      DO N=1,um1%NTILES
!         DO K=1,um1%TILE_PTS(N)
!            L = um1%TILE_INDEX(K,N)
!            J=(um1%LAND_INDEX(L)-1)/um1%ROW_LENGTH + 1
!            I = um1%LAND_INDEX(L) - (J-1)*um1%ROW_LENGTH
!            
!            ! direct beam visible
!            LAND_ALBEDO_CABLE(I,J,1) = LAND_ALBEDO_CABLE(I,J,1) +              &
!                                       um1%TILE_FRAC(L,N)*ALB_TILE(L,N,1)
!
!            ! diffuse beam visible
!            LAND_ALBEDO_CABLE(I,J,2) = LAND_ALBEDO_CABLE(I,J,2) +              &
!                                       um1%TILE_FRAC(L,N)*ALB_TILE(L,N,2)
!
!            ! direct beam nearinfrared 
!            LAND_ALBEDO_CABLE(I,J,3) = LAND_ALBEDO_CABLE(I,J,3) +              &
!                                       um1%TILE_FRAC(L,N)*ALB_TILE(L,N,3)
!
!            ! diffuse beam nearinfrared
!            LAND_ALBEDO_CABLE(I,J,4) = LAND_ALBEDO_CABLE(I,J,4) +              &
!                                       um1%TILE_FRAC(L,N)*ALB_TILE(L,N,4)
!            LAND_ALB_CABLE(I,J) = LAND_ALB_CABLE(I,J) +                        &
!                                um1%TILE_FRAC(L,N)*LAND_ALB_CABLE_TILE(L,N)
!         ENDDO
!      ENDDO
!
!      cable_runtime%um_radiation= .FALSE.
!
!END SUBROUTINE cable_rad_driver
!
!  
!
!
!
!SUBROUTINE cable_hyd_driver( SNOW_TILE, LYING_SNOW, SURF_ROFF, SUB_SURF_ROFF,  &
!                             TOT_TFALL )
!
!   USE cable_data_module,   ONLY : PHYS, OTHER
!   USE cable_common_module!, only : cable_runtime, cable_user
!   USE cable_um_tech_mod, only : um1, ssnow, canopy, veg
!   IMPLICIT NONE
!
!   REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS,um1%NTILES) ::                    &
!      SNOW_TILE   ! IN Lying snow on tiles (kg/m2)        
!
!   REAL, INTENT(OUT), DIMENSION(um1%LAND_PTS) ::                               &
!      LYING_SNOW,    & ! OUT Gridbox snowmass (kg/m2)        
!      SUB_SURF_ROFF, & !
!      SURF_ROFF,     & !
!      TOT_TFALL        !
!
!   REAL, DIMENSION(um1%LAND_PTS,um1%NTILES) ::                                 &
!      SURF_CAB_ROFF,    &
!      TOT_TFALL_TILE                
!
!   REAL :: miss =0. 
!   REAL, POINTER :: TFRZ
!      
!      TFRZ => PHYS%TFRZ
!   
!      SNOW_TILE= UNPACK(ssnow%snowd, um1%L_TILE_PTS, miss) 
!      LYING_SNOW = SUM(um1%TILE_FRAC * SNOW_TILE,2) !gridbox snow mass
!
!      SURF_CAB_ROFF  = UNPACK(ssnow%rnof1, um1%L_TILE_PTS, miss)
!      SURF_ROFF      = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)
!      
!      ! Don't include sub-soil drainage for lakes
!      ! NB: Hard-wired type to be removed in future version
!      WHERE( veg%iveg == 16 ) ssnow%rnof2 = 0.0
!  
!      SURF_CAB_ROFF  = UNPACK(ssnow%rnof2, um1%L_TILE_PTS, miss)
!      SUB_SURF_ROFF  = SUM(um1%TILE_FRAC * SURF_CAB_ROFF,2)
!
!      TOT_TFALL_TILE = UNPACK(canopy%through, um1%L_TILE_PTS, miss)
!      TOT_TFALL      = SUM(um1%TILE_FRAC * TOT_TFALL_TILE,2)
!      
!END SUBROUTINE cable_hyd_driver
!      



END MODULE cable_main_mod
