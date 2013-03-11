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
           ! constructed for CABLE
           cos_zenith_angle, & 
           ! these exist in JULES/UM 
           tl_1, qw_1, vshr_land, pstar, z1_tq,&
           z1_uv, &
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
           !jhan: make sure NOT used and delete
           ! CABLE uses patitioned surf_down_sw in UM, 
           ! however is not available in JULES
           surf_down_sw, &
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

   
!!   !--- reads runtime and user switches and reports
!!   USE cable_um_tech_mod, ONLY : cable_um_runtime_vars, air, bgc, canopy,      &
!!                                 met, bal, rad, rough, soil, ssnow, sum_flux, veg 
!!   
!!   !--- vars common to CABLE declared 
!!   USE cable_common_module, ONLY : cable_runtime, cable_user, ktau_gl,          &
!!                                   knode_gl, kwidth_gl, kend_gl
!!   
!!   !--- subr to (manage)interface UM data to CABLE
!!   USE cable_um_init_mod, ONLY : interface_UM_data
!!   
!!   !--- subr to call CABLE model
!!   USE cable_cbm_module, ONLY : cbm
!!
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
      tile_index,& ! index of tile points being processed
      isnow_flg3l   ! 3 layer snow flag

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

   REAL, INTENT(IN), DIMENSION(land_pts, ntiles) ::                            &
      snow_rho1l, &
      snage_tile
   
   REAL, INTENT(IN), DIMENSION(row_length, rows, 4) ::                         &
      surf_down_sw 
   
   REAL, INTENT(IN), DIMENSION(land_pts, npft) ::                              &
      canht_ft, lai_ft 
   
   REAL, INTENT(IN),DIMENSION(land_pts, ntiles) ::                             &
      canopy_tile
   
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
   LOGICAL, SAVE :: first_cable_call = .TRUE.
   
      CONV_RAIN = 0.0 
      CONV_SNOW = 0.0
   !jhan: because intent(out) need init where dependencies not veiwed in testing
      FTL_TILE=0.0
      FQW_TILE=0.
      TSTAR_TILE=0.
      U_S=0.
      U_S_STD_TILE=0.
      CH_TILE=0.
      CD_TILE=0.
      RADNET_TILE=0.
      FRACA=0.
      RESFS=0.
      RESFT=0.
      Z0H_TILE=0.
      Z0M_TILE=0.
      RECIP_L_MO_TILE=0.
      EPOT_TILE=0.
 

   !jhan: split the interface_data


   print *, " "
   print *, "JHAN in CABLE ", timestep_number
   print *, " "

!jhan: use after defined cable%
!   !--- basic info from global model passed to cable_common_module 
!   !--- vars so don't need to be passed around, just USE _module
!   ktau_gl = timestep_number     !timestep of EXPERIMENT not necesarily 
!                                 !the same as timestep of particular RUN
!   knode_gl = mype               !which processor am i on?
!   kwidth_gl = timestep          !width of timestep (secs)
!   kend_gl = endstep             !timestep of EXPERIMENT not necesarily 
!
!   !--- user FLAGS, variables etc def. in cable.nml is read on 
!   !--- first time step of each run. these variables are read at 
!   !--- runtime and for the most part do not require a model rebuild.
!   IF(first_cable_call) THEN
!      CALL cable_um_runtime_vars(runtime_vars_file) 
!      first_cable_call = .FALSE.
!   ENDIF      
!
!
!   !---------------------------------------------------------------------!
!   !--- initialize CABLE using UM forcings etc. these args are passed ---!
!   !--- down from ~UM. explicit_call                                  ---! 
!   !---------------------------------------------------------------------!
!   CALL explicit_call_initialization(                                          & 
!            ! cable% type yet to be officailly implemented
!            row_length, & ! -> cable%row_length 
!            rows,       & ! -> cable%rows
!            land_pts,   & ! -> cable%land_pts
!            ntiles,     & ! -> cable%ntiles
!            npft,       & ! -> cable%npft 
!            sm_levels,  & ! -> cable%ms
!            timestep,   & ! -> cable%timestep_width
!            latitude,   & ! -> cable%latitude
!            longitude,  & ! -> cable%longitude
!            land_index, & ! -- necessary for packing 
!            tile_frac,  & ! -> cable%tile_frac
!            tile_pts,   & ! -> cable%
!            tile_index, & ! -- necessary for packing
!            dzsoil,     & ! -> soil%zse                        
!            ! soil properties from UM/JULES
!            bexp,       & ! -> soil%bch
!            hcon,       & ! ~> soil%cnsd
!            satcon,     & ! ~> soil%hyds
!            sathh,      & ! -> soil%sucs
!            smvcst,     & ! -> soil%ssat
!            smvcwt,     & ! -> soil%swilt
!            smvccl,     & ! -> soil%sfc
!            albsoil,    & ! -> soil%albsoil
!            ! canopy properties from UM/JULES
!            canht_ft,   & ! ~> veg%hc
!            lai_ft,     & ! ~> veg%lai
!            ! forcing from JULES
!            sw_down,    & ! ~> met%fsd
!            lw_down,    & ! -> met%fld 
!            ls_rain,    & ! ~> met%precip
!            ls_snow,    & ! ~> met%precip_sn
!            tl_1,       & ! -> met%tk
!            qw_1,       & ! -> met%qv
!            vshr_land,  & ! -> met%ua
!            pstar,      & ! ~> met%pmb
!            z1_tq,      & ! -> rough%za_tq
!            z1_uv,      & ! -> rough%za_uv
!            canopy_tile,& ! -> canopy%cansto
!            Fland,      & ! -> ssnow%fland
!            CO2_MMR,    & ! ~> met%ca
!         
!            !jhan:adapted from JULES var cosz, done elsewhere, move to here 
!            ! and make switchable
!            cos_zenith_angle,    & ! ->met%coszen
!         
!            ! snow properties from UM/JULES
!            snow_tile,     & ! -> ssnow%snowd
!         
!            ! snow properties from CABLE vars 
!            snage_tile,    & ! -> ssnow%snage
!            snow_rho1l,    & ! -> ssnow%ssdnn
!            isnow_flg3l,   & ! -> ssnow%isflag
!            snow_rho3l,    & ! -> ssnow%ssdn
!            snow_depth3l,  & ! -> ssnow%sdepth
!            snow_tmp3l,    & ! -> ssnow%tggsn
!            snow_mass3l,   & ! -> ssnow%smass
!            snow_cond,     & ! -> ssnow%sconds
!         
!         
!            !soil properties from CABLE vars 
!            smcl_tile,     & ! ~> soil%wb
!            sthf_tile,     & ! ~> soil%wbice
!            tsoil_tile     & ! -> ssnow%tgg
!   )                         

!   CALL implicit_call_initialization(                                          & 
!            ls_rain,    & ! ~>  met%precip
!            ls_snow,    & ! ~>  met%precip_sn
!            conv_rain,  & ! ~~> met%precip
!            conv_snow,  & ! ~~> met%precip_sn
!            dtl_1,      & ! ~~> met%precip_sn
!            dqw_1       & ! ~~> met%precip_sn
!   )
!    



!
!   !---------------------------------------------------------------------!
!   !--- real(timestep) width, CABLE types passed to CABLE "engine" as ---!  
!   !--- req'd by Mk3L  --------------------------------------------------!
!   !---------------------------------------------------------------------!
!   CALL cbm( REAL(timestep), air, bgc, canopy, met, bal,                       &
!             rad, rough, soil, ssnow, sum_flux, veg )
!
!
!
!
!   !---------------------------------------------------------------------!
!   !--- pass land-surface quantities calc'd by CABLE in explicit call ---!
!   !--- back to UM.                                                   ---!
!   !---------------------------------------------------------------------!
!   call cable_expl_unpack( FTL_TILE_cable, FTL_cable, FTL_TILE, FQW_TILE,          &
!                           LE_TILE_cable, LE_cable, TSTAR_TILE, TSTAR_TILE_cable,    &
!                           TSTAR_cable, U_S, U_S_STD_TILE, U_S_cable, CH_cable,      &
!                           CD_cable, CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
!                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
!                           RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts,             &
!                           ssnow%snowd, ssnow%cls, air%rlam, air%rho,          &
!                           canopy%fe, canopy%fh, canopy%us, canopy%cdtq,       &
!                           canopy%fwet, canopy%wetfac_cs, canopy%rnet,         &
!                           canopy%zetar, canopy%epot, met%ua, rad%trad,        &
!                           rad%transd, rough%z0m, rough%zref_tq )

!      CALL implicit_unpack( TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,                &
!                            SMVCST, STHF, STHF_TILE, STHU, STHU_TILE,          &
!                            snow_tile, SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L,  &
!                            SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,    &
!                            FTL_TILE_CAB, FTL_CAB, LE_TILE_CAB, LE_CAB,        &
!                            FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
!                            TSTAR_TILE_CAB, TSTAR_CAB, SMCL_CAB, TSOIL_CAB,    &
!                            SURF_HTF_CAB, SURF_HT_FLUX_LAND, ECAN_TILE,        &
!                            ESOIL_TILE, EI_TILE, RADNET_TILE, TOT_ALB,         &
!                            SNAGE_TILE, CANOPY_TILE, GS, T1P5M_TILE,           &
!                            Q1P5M_TILE, CANOPY_GB, FLAND, MELT_TILE, DIM_CS1,  &
!                            DIM_CS2, NPP, NPP_FT, GPP, GPP_FT, RESP_S,         &
!                            RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT, G_LEAF )
! 

END SUBROUTINE cable_main


!SUBROUTINE implicit_call_initialization(                                    & 
!            ls_rain,    & ! ~>  met%precip
!            ls_snow,    & ! ~>  met%precip_sn
!            conv_rain,  & ! ~~> met%precip
!            conv_snow,  & ! ~~> met%precip_sn
!            dtl_1,      & ! ~~> met%precip_sn
!            dqw_1       & ! ~~> met%precip_sn
!   )
!    
!   REAL, INTENT(IN), DIMENSION(row_length,rows) ::                             &
!      ls_rain,    &
!      ls_snow     
!
!   REAL, DIMENSION(ROW_LENGTH,ROWS) ::                                 &
!      CONV_RAIN, & ! IN Convective rain
!      CONV_SNOW   ! IN Convective snow
!   
!   REAL, DIMENSION(ROW_LENGTH,ROWS) ::                                 &
!      DTL_1,    & ! IN Level 1 increment to T field 
!      DQW_1       ! IN Level 1 increment to q field 
!
!   !jhan:mp must be foundbefore we can do this
!   REAL, DIMENSION(mp) ::                                                      & 
!      dtlc, & 
!      dqwc
!
!      dtlc = 0. ; dqwc = 0.
!
!      !--- All these subrs do is pack a CABLE var with a UM var.
!      !-------------------------------------------------------------------
!      !--- UM met forcing vars needed by CABLE which have UM dimensions
!      !---(rowlength,rows)[_rr], which is no good to CABLE. These have to be 
!      !--- re-packed in a single vector of active tiles. Hence we use 
!      !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
!      !--- if the land point is/has an active tile
!      !--- generic format:
!      !--- um2cable_rr( UM var, default value for snow tile, CABLE var, mask )
!      !--- where mask tells um2cable_rr whether or not to use default value 
!      !--- for snow tile 
!      !-------------------------------------------------------------------
!      CALL um2cable_rr( (LS_RAIN+CONV_RAIN)*um1%TIMESTEP, met%precip)
!      CALL um2cable_rr( (LS_SNOW+CONV_SNOW)*um1%TIMESTEP, met%precip_sn)
!      CALL um2cable_rr( dtl_1, dtlc)
!      CALL um2cable_rr( dqw_1, dqwc)
!      
!      !--- conv_rain(snow)_prevstep are added to precip. in explicit call
!      CALL um2cable_rr( (CONN_RAIN)*um1%TIMESTEP, conv_rain_prevstep)
!      CALL um2cable_rr( (CONV_snow)*um1%TIMESTEP, conv_snow_prevstep)
!      
!      met%precip   =  met%precip + met%precip_sn
!      met%tk = met%tk + dtlc
!      met%qv = met%qv + dqwc
!      met%tvair = met%tk
!      met%tvrad = met%tk
!   END SUBROUTINE implicit_call_initialization
 


!
!
!
!
!!---------------------------------------------------------------------!
!!--- pass land-surface quantities calc'd by CABLE in explicit call ---!
!!--- back to UM.                                                   ---!
!!---------------------------------------------------------------------!
!SUBROUTINE cable_expl_unpack( FTL_TILE_cable, FTL_cable, FTL_TILE, FQW_TILE,       &
!                           LE_TILE_cable, LE_cable, TSTAR_TILE, TSTAR_TILE_cable,    &
!                           TSTAR_cable, U_S, U_S_STD_TILE, U_S_cable, CH_cable,      &
!                           CD_cable, CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
!                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
!                           RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts,             &
!                           ssnow_snowd, ssnow_cls, air_rlam, air_rho,          &
!                           canopy_fe, canopy_fh, canopy_us, canopy_cdtq,       &
!                           canopy_fwet, canopy_wetfac_cs, canopy_rnet,         &
!                           canopy_zetar, canopy_epot, met_ua, rad_trad,        &
!                           rad_transd, rough_z0m, rough_zref_tq )
!
!   USE cable_def_types_mod, ONLY : mp, NITER 
!   USE cable_data_module,   ONLY : PHYS
!   USE cable_um_tech_mod,   ONLY : um1
!   USE cable_common_module, ONLY : cable_runtime, cable_user, &
!                                   ktau_gl, knode_gl 
!   IMPLICIT NONE         
!
!
!   !-------------------------------------------------------------------------- 
!   !--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
!   !-------------------------------------------------------------------------- 
!
!
!   !___ UM variables to recieve unpacked CABLE vars
!
!   !___return fluxes
!   REAL, INTENT(OUT), DIMENSION(land_pts) ::   &
!      FTL_cable, &
!      LE_cable
!   REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
!      FTL_TILE_cable, &
!      FTL_TILE,   &  ! Surface FTL for land tiles     
!      FQW_TILE,   &  ! Surface FQW for land tiles     
!      LE_TILE_cable
!
!   !___return temp and roughness
!   REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
!      TSTAR_TILE_cable, TSTAR_TILE,  Z0H_TILE, Z0M_TILE
!   REAL, INTENT(OUT), DIMENSION(land_pts) ::                  &
!      TSTAR_cable
!
!   !___return friction velocities/drags/ etc
!   REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
!      U_S_STD_TILE      ! Surface friction velocity
!   REAL, INTENT(OUT), DIMENSION(row_length,rows)  :: &
!      U_S               ! Surface friction velocity (m/s)
!   REAL, INTENT(OUT), DIMENSION(land_pts) ::                  &
!      CH_cable,  &  ! Turbulent surface exchange
!      CD_cable,  &  ! Turbulent surface exchange
!      U_S_cable     ! Surface friction velocity (m/s)
!
!   !___return miscelaneous 
!   REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
!      RADNET_TILE,   & ! Surface net radiation
!      RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
!                       ! factor for fraction (1-FRACA) of snow-free land tiles
!      RESFT,         & ! Total resistance factor.
!                       ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,
!                       ! 1 for snow.    
!      FRACA,         & ! Fraction of surface moisture
!      EPOT_TILE
!   
!   LOGICAL,DIMENSION(land_pts,ntiles) :: l_tile_pts
!
!   !___UM vars used but NOT returned 
!   REAL, INTENT(IN), DIMENSION(land_pts) ::   &
!      FLAND(land_pts)              ! IN Land fraction on land tiles.
!
!
!
!
!   !___ decs of intent(in) CABLE variables to be unpacked
!
!   ! snow depth (liquid water), factor for latent heat
!   REAL, INTENT(IN), DIMENSION(mp) :: ssnow_snowd, ssnow_cls
!   
!   ! surface wind speed (m/s)
!   REAL, INTENT(IN), DIMENSION(mp) :: met_ua 
!   
!   ! latent heat for water (j/kg), dry air density (kg m-3)
!   REAL, INTENT(IN), DIMENSION(mp) :: air_rlam, air_rho 
!   
!   ! frac SW diffuse transmitted thru canopy, rad. temp. (soil and veg)
!   REAL, INTENT(IN), DIMENSION(mp) :: rad_trad,rad_transd 
!   
!   ! total latent heat (W/m2), total sensible heat (W/m2)
!   REAL, INTENT(IN), DIMENSION(mp) :: canopy_fe, canopy_fh  
!   
!   ! fraction of canopy wet
!   REAL, INTENT(IN), DIMENSION(mp) :: canopy_fwet, canopy_wetfac_cs
!   
!   ! friction velocity, drag coefficient for momentum
!   REAL, INTENT(IN), DIMENSION(mp) :: canopy_us, canopy_cdtq
!   
!   ! net rad. absorbed by surface (W/m2), total potential evaporation 
!   REAL, INTENT(IN), DIMENSION(mp) :: canopy_rnet, canopy_epot        
!   
!   ! stability correction
!   REAL, INTENT(IN), DIMENSION(mp,niter) :: canopy_zetar
!   
!   ! roughness length, Reference height for met forcing
!   REAL, INTENT(IN), DIMENSION(mp) :: rough_z0m, rough_zref_tq 
! 
!   !___return friction velocities/drags/ etc
!   REAL, DIMENSION(land_pts,ntiles) :: &
!      CD_TILE,    &     ! Drag coefficient
!      CH_TILE,    &     ! Transfer coefficient for heat & moisture
!      RECIP_L_MO_TILE ! Reciprocal of the Monin-Obukhov length for tiles (m^-1)
!   
!   !-------------------------------------------------------------------------- 
!   !--- end INPUT ARGS FROM cable_explicit_driver() --------------------------
!   !-------------------------------------------------------------------------- 
!
!
!   
!        
!   !___vars in local calc. of latent heat fluxes
!   REAL, DIMENSION(land_pts,ntiles) ::                  &
!      FQW_TILE_cable,  &
!      LE_TILE
!
!   !___vars in local calc of Surface friction velocities
!   REAL, DIMENSION(land_pts,ntiles) ::                  &
!      CD_cable_TILE,   &  
!      CH_cable_TILE,   &  ! (bulk transfer) coeff. for momentum
!      U_S_TILE
!   REAL, DIMENSION(mp)  :: &
!      CDCAB,CHCAB
!
!   !___local miscelaneous
!   REAL, DIMENSION(mp)  :: &
!   THETAST,fraca_cab,rfsfs_cab, RECIPLMOTILE, fe_dlh
!   INTEGER :: i,j,k,N,L
!   REAL :: miss = 0.0
!   LOGICAL, SAVE :: first_cable_call = .true.
!   REAL, POINTER :: CAPP 
!   
!      CAPP => PHYS%CAPP
!
!   !jhan: because intent(out) need init where dependencies not veiwed in testing
!      FTL_cable=0.0
!      FTL_TILE_cable=0.0
!      FTL_TILE=0.0
!      LE_cable=0.
!      LE_TILE_cable=0.
!      FQW_TILE=0.
!      TSTAR_TILE=0.
!      TSTAR_TILE_cable=0.
!      TSTAR_cable=0.
!      U_S=0.
!      U_S_STD_TILE=0.
!      U_S_cable=0.
!      CH_cable=0.
!      CD_cable=0.
!      CH_TILE=0.
!      CD_TILE=0.
!      RADNET_TILE=0.
!      FRACA=0.
!      RESFS=0.
!      RESFT=0.
!      Z0H_TILE=0.
!      Z0M_TILE=0.
!      RECIP_L_MO_TILE=0.
!      EPOT_TILE=0.
! 
!
!      
!      !___return fluxes
!      FTL_TILE_cable = UNPACK(canopy_fh,  l_tile_pts, miss)
!      FTL_cable = SUM(TILE_FRAC * FTL_TILE_cable,2)
!      FQW_TILE_cable = UNPACK(canopy_fe,  l_tile_pts, miss)
!      LE_TILE_cable = UNPACK(canopy_fe,  l_tile_pts, miss)
!      LE_cable = SUM(TILE_FRAC * LE_TILE_cable,2)
!      fe_dlh = canopy_fe/(air_rlam*ssnow_cls)
!      FTL_TILE = UNPACK(canopy_fh,  l_tile_pts, miss)
!      FTL_TILE = FTL_TILE / capp
!      FQW_TILE = UNPACK(fe_dlh, l_tile_pts, miss)
!
!      !___return temp and roughness
!      TSTAR_TILE_cable = UNPACK(rad_trad, l_tile_pts, miss)
!      TSTAR_cable = SUM(TILE_FRAC * TSTAR_TILE_cable,2)
!      TSTAR_TILE = UNPACK(rad_trad,  l_tile_pts, miss)
!!      Z0M_TILE = UNPACK(rough_z0m,  l_tile_pts, miss)
!!      Z0H_TILE = Z0M_TILE
!      
!      !___return friction velocities/drags/ etc
!!      U_S_TILE  =  UNPACK(canopy_us, l_tile_pts, miss)
!      U_S_cable  = SUM(TILE_FRAC *  U_S_TILE,2)
!      CDCAB = canopy_us**2/met_ua**2   ! met%ua is always above umin = 0.1m/s
!      ! for Cable CD*
!      CD_cable_TILE =  UNPACK(CDCAB,l_tile_pts, miss)
!      CD_cable= SUM(TILE_FRAC * CD_cable_TILE,2)
!      ! for Cable CH*
!      CH_cable_TILE =  UNPACK(canopy_cdtq,l_tile_pts, miss)
!      CH_cable= SUM(TILE_FRAC * CH_cable_TILE,2)
!
!!      U_S_STD_TILE=U_S_TILE
!!      CD_TILE = CD_cable_TILE
!!      CH_TILE = CH_cable_TILE
!!
!!      U_S = 0.
!!      DO N=1,ntiles
!!         DO K=1,TILE_PTS(N)
!!            L = TILE_INDEX(K,N)
!!            J=(LAND_INDEX(L)-1)/row_length + 1
!!            I = LAND_INDEX(L) - (J-1)*row_length
!!            U_S(I,J) = U_S(I,J)+FLAND(L)*TILE_FRAC(L,N)*U_S_TILE(L,N)
!!         ENDDO
!!      ENDDO
!!
!!
!!
!!
!!      !___return miscelaneous 
!!      fraca_cab = canopy_fwet * (1.-rad_transd)
!!      WHERE( ssnow_snowd > 1.0 ) fraca_cab = 1.0
!!      rfsfs_cab = MIN( 1., MAX( 0.01, canopy_wetfac_cs - fraca_cab ) /         &
!!                  MAX( 0.01,1. - fraca_cab ) )
!!      FRACA = UNPACK( fraca_cab, l_tile_pts, miss )
!!      RESFT = UNPACK( canopy_wetfac_cs,l_tile_pts, miss )
!!      RESFS = UNPACK( rfsfs_cab , l_tile_pts, miss )
!!
!!      RADNET_TILE = UNPACK( canopy_rnet , l_tile_pts, miss )
!!      THETAST = ABS( canopy_fh ) / ( air_rho * capp*canopy_us )
!!      RECIPLMOTILE =  canopy_zetar(:,niter) / rough_zref_tq
!!      RECIP_L_MO_TILE = UNPACK( RECIPLMOTILE, l_tile_pts, miss )
!!      EPOT_TILE = UNPACK( canopy_epot, l_tile_pts, miss )
!      
!
!      IF(first_cable_call) THEN 
!         l_tile_pts = l_tile_pts
!         first_cable_call = .FALSE.
!      ENDIF
!
!   
!END SUBROUTINE cable_expl_unpack
 

END MODULE cable_main_mod
