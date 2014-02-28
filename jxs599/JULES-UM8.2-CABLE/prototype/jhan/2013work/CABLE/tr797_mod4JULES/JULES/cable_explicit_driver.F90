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

  SUBROUTINE cable_explicit(                                           &
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
   REAL, INTENT(IN) :: timestep
   INTEGER, INTENT(IN) :: endstep, timestep_number, mype

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

!jha
   !___JULES soil/snow/radiation/met vars
   !___JULES includes this as (:,:) 
   REAL, INTENT(IN), DIMENSION(land_pts, sm_levels) :: & 
      bexp     ! => parameter b in Campbell equation
      
       
   !___UM soil/snow/radiation/met vars
   REAL, INTENT(IN), DIMENSION(land_pts) :: & 
      !bexp,    & ! => parameter b in Campbell equation 
      hcon,    & ! Soil thermal conductivity (W/m/K).
      satcon,  & ! hydraulic conductivity @ saturation [mm/s]
      sathh,   &
      smvcst,  &
      smvcwt,  &
      smvccl,  &
      albsoil

   REAL, INTENT(IN), DIMENSION(land_pts) :: & 
      fland 
   
   !REAL, INTENT(INOUT), DIMENSION(row_length,rows) :: &
   REAL, DIMENSION(row_length,rows) :: &
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

   !REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles) ::                         &
   REAL, DIMENSION(land_pts, ntiles) ::                         &
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
   
   !REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles,3) ::                       &
   REAL, DIMENSION(land_pts, ntiles,3) ::                       &
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
   !REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
   REAL, DIMENSION(land_pts,ntiles) :: &
      U_S_STD_TILE      ! Surface friction velocity

   !REAL, INTENT(OUT), DIMENSION(row_length,rows)  :: &
   REAL, DIMENSION(row_length,rows)  :: &
      U_S               ! Surface friction velocity (m/s)
   
   !___return miscelaneous 
   !REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
   REAL, DIMENSION(land_pts,ntiles) :: &
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
   kwidth_gl = INT(timestep)     !width of timestep (secs)
   kend_gl = endstep             !timestep of EXPERIMENT not necesarily 

   !--- user FLAGS, variables etc def. in cable.nml is read on 
   !--- first time step of each run. these variables are read at 
   !--- runtime and for the most part do not require a model rebuild.
   
   cable_runtime%um = .TRUE. 
   cable_runtime%um_explicit = .TRUE. 
   
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


   canopy%oldcansto=canopy%cansto

   !---------------------------------------------------------------------!
   !--- real(timestep) width, CABLE types passed to CABLE "engine" as ---!  
   !--- req'd by Mk3L  --------------------------------------------------!
   !---------------------------------------------------------------------!
   CALL cbm( timestep, air, bgc, canopy, met, bal,                             &
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

print *, "jhan: end cable_explicit"

END SUBROUTINE cable_explicit



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
   !REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: &
      FQW_TILE       ! Surface FQW for land tiles     

   !___return temp and roughness
   !REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: &
      TSTAR_TILE,  Z0H_TILE, Z0M_TILE

   !___return friction velocities/drags/ etc
   !REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: &
      CD_TILE,    &     ! Drag coefficient
      CH_TILE,    &     ! Transfer coefficient for heat & moisture
      U_S_STD_TILE      ! Surface friction velocity
   !REAL, INTENT(OUT), DIMENSION(um1%row_length,um1%rows)  :: &
   REAL, DIMENSION(um1%row_length,um1%rows)  :: &
      U_S               ! Surface friction velocity (m/s)

   !___return miscelaneous 
   !REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: &
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

!========================================================================= 
!========================================================================= 
!========================================================================= 

