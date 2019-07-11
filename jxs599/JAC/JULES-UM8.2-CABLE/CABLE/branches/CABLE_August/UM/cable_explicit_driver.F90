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


   SUBROUTINE cable_explicit_driver(                                        &
               ! vars native to JULES/UM 
               ! ........................
               ! # grid cells 
               row_length, rows, &
               ! # AND index of land points
               land_pts, land_index, &
               ! model levels (modified for CABLE)
               ntiles, npft, sm_levels,   &
               ! 
               dim_cs1, dim_cs2,   &
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
               lai_ft,         &
               dzsoil,       &
               surf_down_sw,       &
               ! soil vars used in initialization of CABLE vars
               ! INTENT(IN)  used on first CALL only
               !jhan: ultimately read in as tiled vars
                bexp, & 
                !currently being passed 
                hcon, &
                satcon, sathh, &
                smvcst,           &
                smvcwt, smvccl, &
                albsoil, & ! SOIL_ALB (atm_step)
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
                ! JULES forcing  
                !sw_down_cable, &
                !lw_down_cable,   &
                lw_down, &
                ! large scale rain, snow
                ls_rain, ls_snow, &
                cos_zenith_angle, & 
                tl_1, qw_1, vshr_land, pstar, z1_tq,&
                z1_uv, &
                !
                ! End - vars native to JULES/UM 
                !
                ! new CABLE vars
                !
                ! CABLE_vars initialized from ancillaries
                !
                ! snow vars - passed as read 
                isnow_flg3l, snow_rho1l, snow_age,          &
                ! snow vars - (JULES-read as separate var per layer and pre-packed) 
                snow_rho3l, snow_depth3l, snow_tmp3l, snow_mass3l, &
                ! soil vars - (JULES-read as separate var per layer and pre-packed)
                !sthf_tile is set_atm FIELD, but sthu_tile is simply dec in atm_step
                sthu_tile, smcl_tile, sthf_tile, tsoil_tile, &
                !
                ! snow_cond requires no init from file
                snow_cond, &
                ! 
                FTL_TILE,  &
                FQW_TILE, TSTAR_TILE,   &
                U_S, &
                U_S_STD_TILE,&
                RADNET_TILE, FRACA, rESFS, RESFT, Z0H_TILE,  &
                Z0M_TILE, EPOT_TILE, cd_tile & 
                )

   !--- reads runtime and user switches and reports
   USE cable_um_tech_mod, ONLY : cable_um_runtime_vars, air, bgc, canopy,      &
                                 met, bal, rad, rough, soil, ssnow, sum_flux, veg 
   
   !--- vars common to CABLE declared 
   USE cable_common_module, ONLY : cable_runtime, cable_user, ktau_gl,         &
                                   knode_gl, kwidth_gl, kend_gl,               &
                                   report_version_no
   
   !--- subr to (manage)interface UM data to CABLE
   USE cable_um_init_mod, ONLY : interface_UM_data
   
   !--- subr to call CABLE model
   USE cable_cbm_module, ONLY : cbm

   USE cable_def_types_mod, ONLY : mp, ms, ssnow, rough, canopy, air, rad,     &
                                   met

   USE cable_data_module,   ONLY : PHYS
   
   !--- include subr called to , write data for testing purposes 
   USE cable_diag_module

   IMPLICIT NONE
 
 
  
   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM sf_exch() --------------------------------------------
   !-------------------------------------------------------------------------- 
   
   !___IN: UM dimensions, array indexes, flags
   !INTEGER, INTENT(IN) ::                                                      & 
   INTEGER ::                                                      & 
      row_length, rows, & ! UM grid resolution
      land_pts,         & ! # of land points being processed
      ntiles,           & ! # of tiles 
      npft,             & ! # of plant functional types
      sm_levels,        & ! # of soil layers 
      dim_cs1, dim_cs2

   ! index of land points being processed
   !INTEGER, INTENT(IN), DIMENSION(land_pts) :: land_index 
   INTEGER, DIMENSION(land_pts) :: land_index 

   ! # of land points on each tile
   !INTEGER, INTENT(IN), DIMENSION(ntiles) :: tile_pts 
   INTEGER, DIMENSION(ntiles) :: tile_pts 
   
   !INTEGER, INTENT(IN), DIMENSION(land_pts, ntiles) ::                         & 
   INTEGER, DIMENSION(land_pts, ntiles) ::                         & 
      tile_index ,& ! index of tile points being processed
      isnow_flg3l   ! 3 layer snow flag

   !--- TRUE if land, F elsewhere.
   !jhan:rm land_mask
   LOGICAL,DIMENSION(row_length,rows) :: land_mask   

   !___UM parameters: water density, soil layer thicknesses 
   !REAL, INTENT(IN), DIMENSION(sm_levels) :: dzsoil
   REAL, DIMENSION(sm_levels) :: dzsoil

   !___UM soil/snow/radiation/met vars
   !REAL, INTENT(IN), DIMENSION(land_pts) :: & 
   REAL, DIMENSION(land_pts) :: & 
      bexp,    & ! => parameter b in Campbell equation 
      hcon,    & ! Soil thermal conductivity (W/m/K).
      satcon,  & ! hydraulic conductivity @ saturation [mm/s]
      sathh,   &
      smvcst,  &
      smvcwt,  &
      smvccl,  &
      albsoil, &
      fland 
   
   !REAL, INTENT(INOUT), DIMENSION(row_length,rows) :: &
   REAL, DIMENSION(row_length,rows) :: &
      cos_zenith_angle
   
   !REAL, INTENT(IN), DIMENSION(row_length,rows) ::                             &
   REAL, DIMENSION(row_length,rows) ::                             &
      latitude,   &
      longitude,  &
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

   !REAL, INTENT(IN), DIMENSION(land_pts, ntiles) ::                            &
   REAL, DIMENSION(land_pts, ntiles) ::                            &
      tile_frac,  &    
      snow_rho1l, &
      snow_age
   
   !REAL, INTENT(IN), DIMENSION(row_length, rows, 4) ::                        &
   REAL,  DIMENSION(row_length, rows, 4) ::                        &
      surf_down_sw 
   
   !REAL, INTENT(IN), DIMENSION(land_pts, npft) ::                              &
   REAL, DIMENSION(land_pts, npft) ::                              &
      canht_ft, lai_ft 
   
   !REAL, INTENT(IN),DIMENSION(land_pts, ntiles) ::                             &
   REAL, DIMENSION(land_pts, ntiles) ::                             &
      canopy_tile
   
   !REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles,3) ::                       &
   REAL, DIMENSION(land_pts, ntiles,3) ::                       &
      snow_cond
   
   !REAL, INTENT(IN), DIMENSION(land_pts, ntiles,3) ::                          &
   REAL, DIMENSION(land_pts, ntiles,3) ::                          &
      snow_rho3l,    &
      snow_depth3l,  &
      snow_mass3l,   &
      snow_tmp3l
   
   !REAL, INTENT(IN), DIMENSION(land_pts, sm_levels) ::                         &
   REAL, DIMENSION(land_pts, sm_levels) ::                         &
      sthu 
   
   !REAL, INTENT(IN), DIMENSION(land_pts, ntiles, sm_levels) :: & 
   REAL, DIMENSION(land_pts, ntiles, sm_levels) :: & 
      sthu_tile, &
      sthf_tile, &
      smcl_tile, &
      tsoil_tile
   
   !REAL, INTENT(IN) :: co2_mmr
   REAL :: co2_mmr

   REAL :: sin_theta_latitude(row_length,rows) 
     
   !___return fluxes
   !REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
   REAL,  DIMENSION(land_pts,ntiles) :: &
      FTL_TILE,   &  ! Surface FTL for land tiles     
      FQW_TILE       ! Surface FQW for land tiles     

   !___return temp and roughness
   !REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
   REAL, DIMENSION(land_pts,ntiles) :: &
      TSTAR_TILE,       & 
      Z0H_TILE,         &
      Z0M_TILE


   !___return friction velocities/drags/ etc
   !REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
   REAL,  DIMENSION(land_pts,ntiles) :: &
      CD_TILE,    &     ! Drag coefficient
!      CH_TILE,    &     ! Transfer coefficient for heat & moisture
      U_S_STD_TILE      ! Surface friction velocity

   !REAL, INTENT(OUT), DIMENSION(row_length,rows)  :: &
   REAL,  DIMENSION(row_length,rows)  :: &
      U_S               ! Surface friction velocity (m/s)
   
   ! end step of experiment, this step, step width, processor num
   !INTEGER, INTENT(IN) :: endstep, timestep_number, mype
   INTEGER :: endstep, timestep_number, mype
   !REAL, INTENT(IN) ::  timestep     
   REAL ::  timestep     
   
   INTEGER:: itimestep
    
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
   
   !___ declare local vars 
   
   !___true IF vegetation (tile) fraction is greater than 0
   LOGICAL, DIMENSION(land_pts, ntiles) :: L_tile_pts
  
   !___ location of namelist file defining runtime vars
   CHARACTER(LEN=200), PARAMETER ::                                            & 
      runtime_vars_file = 'cable.nml'

   !___ 1st call in RUN (!=ktau_gl -see below) 
   LOGICAL, SAVE :: first_cable_call = .TRUE.
 
   character(len=30), dimension(:), allocatable :: expl1CheckNames 
   real, dimension(:,:), allocatable :: expl1CheckFields

   character(len=30), dimension(:), allocatable :: expl2CheckNames 
   real, dimension(:,:,:), allocatable :: expl2CheckFields

   character(len=30), dimension(:), allocatable :: expl3CheckNames 
   real, dimension(:,:,:), allocatable :: expl3CheckFields


!      print *, "jhan:cable_expl row_length ",row_length 
!      print *, "jhan:cable_expl rows ",      rows
!      print *, "jhan:cable_expl land_pts  ", land_pts 
!      print *, "jhan:cable_expl land_index ",shape(land_index)
!      print *, "jhan:cable_expl ntiles  ",   ntiles 
!      print *, "jhan:cable_expl npft  ",     npft 
!      print *, "jhan:cable_expl sm_levels ", sm_levels
!      print *, "jhan:cable_expl dim_cs1  ",  dim_cs1 
!      print *, "jhan:cable_expl dim_cs2 ",   dim_cs2
!      print *, "jhan:cable_expl timestep  ", timestep 
!      print *, "jhan:cable_expl timestep_number ", timestep_number
!      print *, "jhan:cable_expl endstep ",   endstep
!      print *, "jhan:cable_expl mype ",      mype
!
!      print *, "jhan:cable_ex:lw_down ",          lw_down
!
!   CALL cable_farray( land_pts, expl1CheckNames, expl1CheckFields,   &
!      "bexp",   bexp,    & 
!      "hcon",   hcon,    &
!      "satcon", satcon,  &
!      "sathh",  sathh,   &
!      "smvcst", smvcst,  &
!      "smvcwt", smvcwt,  &
!      "smvccl", smvccl,  &
!      "albsoil",albsoil, &
!      "fland",   fland    ) 
!   
!   CALL cable_NaN(expl1CheckNames,expl1CheckFields, knode_gl )
!
!   CALL cable_farray( row_length,rows, expl2CheckNames, expl2CheckFields,   &
!      "cos_zenith_angle", cos_zenith_angle, &
!      "latitude", latitude,   &
!      "longitude", longitude,  &
!      "lw_down",          lw_down,    &
!      "ls_rain",         ls_rain,    &
!      "ls_snow",      ls_snow,    &
!      "tl_1",         tl_1,       &
!      "qw_1",        qw_1,       &  
!      "vshr_land",        vshr_land,  &
!      "pstar",         pstar,      &
!      "z1_tq",        z1_tq,      &
!      "z1_uv", z1_uv )
!
!   CALL cable_NaN(expl2CheckNames,expl2CheckFields, knode_gl )
!   CALL cable_extremes2(expl2CheckNames,expl2CheckFields, knode_gl )
!
!   CALL cable_farray( land_pts, ntiles, expl3CheckNames, expl3CheckFields,   &
!      "snow_tile",snow_tile,  &
!      "tile_frac",tile_frac,  &    
!      "snow_rho1l",snow_rho1l, &
!      "snow_age",  snow_age )
!   
!   CALL cable_NaN(expl3CheckNames,expl3CheckFields, knode_gl )
!   CALL cable_extremes2(expl3CheckNames,expl3CheckFields, knode_gl )

!   !REAL, INTENT(IN), DIMENSION(row_length, rows, 4) ::                        &
!   REAL,  DIMENSION(row_length, rows, 4) ::                        &
!      surf_down_sw 
!   
!   !REAL, INTENT(IN), DIMENSION(land_pts, npft) ::                              &
!   REAL, DIMENSION(land_pts, npft) ::                              &
!      canht_ft, lai_ft 
!   
!   !REAL, INTENT(IN),DIMENSION(land_pts, ntiles) ::                             &
!   REAL, DIMENSION(land_pts, ntiles) ::                             &
!      canopy_tile
!   
!   !REAL, INTENT(INOUT), DIMENSION(land_pts, ntiles,3) ::                       &
!   REAL, DIMENSION(land_pts, ntiles,3) ::                       &
!      snow_cond
!   
!   !REAL, INTENT(IN), DIMENSION(land_pts, ntiles,3) ::                          &
!   REAL, DIMENSION(land_pts, ntiles,3) ::                          &
!      snow_rho3l,    &
!      snow_depth3l,  &
!      snow_mass3l,   &
!      snow_tmp3l
!   
!   !REAL, INTENT(IN), DIMENSION(land_pts, sm_levels) ::                         &
!   REAL, DIMENSION(land_pts, sm_levels) ::                         &
!      sthu 
!   
!   !REAL, INTENT(IN), DIMENSION(land_pts, ntiles, sm_levels) :: & 
!   REAL, DIMENSION(land_pts, ntiles, sm_levels) :: & 
!      sthu_tile, &
!      sthf_tile, &
!      smcl_tile, &
!      tsoil_tile
!   
!   !REAL, INTENT(IN) :: co2_mmr
!   REAL :: co2_mmr
!
!   REAL :: sin_theta_latitude(row_length,rows) 
!     
!   !___return fluxes
!   !REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
!   REAL,  DIMENSION(land_pts,ntiles) :: &
!      FTL_TILE,   &  ! Surface FTL for land tiles     
!      FQW_TILE       ! Surface FQW for land tiles     
!
!   !___return temp and roughness
!   !REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
!   REAL, DIMENSION(land_pts,ntiles) :: &
!      TSTAR_TILE,       & 
!      Z0H_TILE,         &
!      Z0M_TILE
!
!
!   !___return friction velocities/drags/ etc
!   !REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
!   REAL,  DIMENSION(land_pts,ntiles) :: &
!      CD_TILE,    &     ! Drag coefficient
!!      CH_TILE,    &     ! Transfer coefficient for heat & moisture
!      U_S_STD_TILE      ! Surface friction velocity
!
!   !REAL, INTENT(OUT), DIMENSION(row_length,rows)  :: &
!   REAL,  DIMENSION(row_length,rows)  :: &
!      U_S               ! Surface friction velocity (m/s)
!   
!   ! end step of experiment, this step, step width, processor num
!   !INTEGER, INTENT(IN) :: endstep, timestep_number, mype
!   INTEGER :: endstep, timestep_number, mype
!   !REAL, INTENT(IN) ::  timestep     
!   REAL ::  timestep     
!   
!   INTEGER:: itimestep
!    
!   !___return miscelaneous 
!   !REAL, INTENT(OUT), DIMENSION(land_pts,ntiles) :: &
!   REAL, DIMENSION(land_pts,ntiles) :: &
!      RADNET_TILE,   & ! Surface net radiation
!      RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
!                       ! factor for fraction (1-FRACA) of snow-free land tiles
!      RESFT,         & ! Total resistance factor.
!                       ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
!                       ! 1 for snow.    
!      FRACA,         & ! Fraction of surface moisture
!      EPOT_TILE
! 


!      print *, "jhan:cable_expl shlw_down", shape( lw_down  )
!      print *, "jhan:cable_expl lw_down",  lw_down 

!      print *, "jhan:cable_expl latitude  ", shape( latitude  )
!      print *, "jhan:cable_expl longitude ",  shape(longitude)
!      print *, "jhan:cable_expl Fland ",  shape(Fland)
!      print *, "jhan:cable_expl tile_frac ",  shape( tile_frac)
!      print *, "jhan:cable_expl tile_pts  ",  shape( tile_pts )
!      print *, "jhan:cable_expl tile_index ",  shape(tile_index)
!      print *, "jhan:cable_expl canht_ft ",  shape(  canht_ft)
!      print *, "jhan:cable_expl lai_ft ",  shape( lai_ft)
!      print *, "jhan:cable_expl dzsoil ",  shape( dzsoil)
!      print *, "jhan:cable_expl  bexp ",  shape( bexp)
!      print *, "jhan:cable_expl  hcon  ",  shape(hcon )
!      print *, "jhan:cable_expl satcon ", shape( satcon ) 
!      print *, "jhan:cable_expl sathh ", shape(  sathh)
!      print *, "jhan:cable_expl smvcst ", shape( smvcst)
!      print *, "jhan:cable_expl smvcwt  ", shape( smvcwt) 
!      print *, "jhan:cable_expl smvccl ", shape(  smvccl)
!      print *, "jhan:cable_expl albsoil ", shape(albsoil)
!      print *, "jhan:cable_expl sthu ", shape( sthu )
!      print *, "jhan:cable_expl snow_tile ", shape( snow_tile )
!      print *, "jhan:cable_expl canopy_tile ", shape( canopy_tile )
!      print *, "jhan:cable_expl CO2_MMR ", shape( CO2_MMR )
!      print *, "jhan:cable_expl lw_down ", shape( lw_down )
!      print *, "jhan:cable_expl ls_rain  ", shape(ls_rain  )
!      print *, "jhan:cable_expl ls_snow ", shape( ls_snow )
!      print *, "jhan:cable_expl cos_zenith_angle ", shape(   cos_zenith_angle )
!      print *, "jhan:cable_expl tl_1  ", shape( tl_1  )
!      print *, "jhan:cable_expl qw_1  ", shape(  qw_1  )
!      print *, "jhan:cable_expl vshr_land  ", shape( vshr_land  )
!      print *, "jhan:cable_expl pstar  ", shape(  pstar  )
!      print *, "jhan:cable_expl z1_tq ", shape( z1_tq )
!      print *, "jhan:cable_expl z1_uv ", shape( z1_uv )
!      print *, "jhan:cable_expl isnow_flg3l  ", shape( isnow_flg3l  )
!      print *, "jhan:cable_expl snow_rho1l  ", shape( snow_rho1l  )
!      print *, "jhan:cable_expl snow_age ", shape( snow_age )
!      print *, "jhan:cable_expl snow_rho3l  ", shape( snow_rho3l  )
!      print *, "jhan:cable_expl snow_depth3l  ", shape( snow_depth3l  )
!      print *, "jhan:cable_expl snow_tmp3l  ", shape( snow_tmp3l  )
!      print *, "jhan:cable_expl snow_mass3l ", shape( snow_mass3l )
!      print *, "jhan:cable_expl sthu_tile  ", shape( sthu_tile  )
!      print *, "jhan:cable_expl smcl_tile  ", shape( smcl_tile  )
!      print *, "jhan:cable_expl sthf_tile  ", shape( sthf_tile  )
!      print *, "jhan:cable_expl tsoil_tile ", shape( tsoil_tile )
!      print *, "jhan:cable_expl snow_cond ", shape( snow_cond )
!      print *, "jhan:cable_expl FTL_TILE ", shape(  FTL_TILE )
!      print *, "jhan:cable_expl FQW_TILE  ", shape(FQW_TILE  )
!      print *, "jhan:cable_expl TSTAR_TILE ", shape( TSTAR_TILE )
!      print *, "jhan:cable_expl U_S ", shape( U_S )
!      print *, "jhan:cable_expl U_S_STD_TILE ", shape( U_S_STD_TILE )
!      print *, "jhan:cable_expl RADNET_TILE  ", shape( RADNET_TILE  )
!      print *, "jhan:cable_expl FRACA  ", shape( FRACA  )
!      print *, "jhan:cable_expl rESFS  ", shape( rESFS  )
!      print *, "jhan:cable_expl RESFT  ", shape( RESFT  )
!      print *, "jhan:cable_expl Z0H_TILE ", shape( Z0H_TILE )
!      print *, "jhan:cable_expl Z0M_TILE ", shape( Z0M_TILE )
!      print *, "jhan:cable_expl EPOT_TILE ", shape(EPOT_TILE )
!      print *, "jhan:cable_expl cd_tile  ", shape( cd_tile  )
   
 

   !--- initialize cable_runtime% switches 
   IF(first_cable_call) THEN
      cable_runtime%um = .TRUE.
      write(6,*) ""
      write(6,*) "CABLE_log"
      CALL report_version_no(6) ! wriite revision number to stdout(6)
   ENDIF
      
      print *, "jhan:cable_expl AA"

   !--- basic info from global model passed to cable_common_module 
   !--- vars so don't need to be passed around, just USE _module
   ktau_gl = timestep_number     !timestep of EXPERIMENT not necesarily 
                                 !the same as timestep of particular RUN
   knode_gl = mype               !which processor am i on?
   itimestep = INT(timestep)    !realize for 'call cbm' pass
   kwidth_gl = itimestep          !width of timestep (secs)
   kend_gl = endstep             !timestep of EXPERIMENT not necesarily 

   !--- internal FLAGS def. specific call of CABLE from UM
   !--- from cable_common_module
   cable_runtime%um_explicit = .TRUE.

   !--- user FLAGS, variables etc def. in cable.nml is read on 
   !--- first time step of each run. these variables are read at 
   !--- runtime and for the most part do not require a model rebuild.
   IF(first_cable_call) THEN
      CALL cable_um_runtime_vars(runtime_vars_file) 
      first_cable_call = .FALSE.
   ENDIF      

      call cable_diag( 1, "lw_down", rows, 1, ktau_gl, knode_gl,            &
                          "lw_down", lw_down(1,:) )
                

      print *, "jhan:cable_expl 1"
! shapes are wrong 
! CALL cable_farray( mp,                      &  
!               'tile_frac',  tile_frac,         & 
!               'bexp',           bexp, &
!               'hcon',           hcon, &
!               'satcon', satcon,& 
!               'sathh', sathh, &
!               'smvcst', smvcst, &
!               'smvcwt', smvcwt,          &
!               'smvccl', smvccl, &
!               'albsoil', albsoil, &
!               'snow_tile', snow_tile,& 
!               'snow_rho1l',snow_rho1l,             &
!               'snow_age', snow_age, &
!               'isnow_flg3l',isnow_flg3l, &
!               'snow_rho3l', snow_rho3l,& 
!               'snow_cond', snow_cond,     &
!               'snow_depth3l', snow_depth3l,& 
!               'snow_tmp3l', snow_tmp3l, &
!               'snow_mass3l', snow_mass3l, &
!               'lw_down', lw_down, &
!               'cos_zenith_angle',  cos_zenith_angle, &
!               'surf_down_sw', surf_down_sw, &
!               'ls_rain', ls_rain,   &
!               'ls_snow', ls_snow, &
!               'tl_1',tl_1, &
!               'qw_1',qw_1, &
!               'vshr_land', vshr_land, &
!               'pstar', pstar, &
!               'z1_tq',z1_tq,       &
!               'z1_uv',z1_uv, &
!               'canopy_tile', canopy_tile, &
!               'Fland',  Fland,   &
!               'sthu_tile',sthu_tile,& 
!               'smcl_tile',smcl_tile, &
!               'sthf_tile',sthf_tile,           &
!               'sthu', sthu, &
!               'tsoil_tile', tsoil_tile, &
!               'canht_ft', canht_ft, &
!               'lai_ft',   lai_ft,                 &
!               'sin_theta_latitude', sin_theta_latitude, &
!               'dzsoil', dzsoil &
!               )                         
!
  !---------------------------------------------------------------------!
  !--- initialize CABLE using UM forcings etc. these args are passed ---!
  !--- down from UM.                                                 ---! 
  !---------------------------------------------------------------------!
   CALL interface_UM_data( row_length, rows, land_pts, ntiles, npft,           & 
                           sm_levels, itimestep, latitude, longitude,          &
                           land_index, tile_frac, tile_pts, tile_index,        &
                           bexp, hcon, satcon, sathh, smvcst, smvcwt,          &
                           smvccl, albsoil, snow_tile, snow_rho1l,             &
                           snow_age, isnow_flg3l, snow_rho3l, snow_cond,     &
                           snow_depth3l, snow_tmp3l, snow_mass3l, &
                           lw_down, cos_zenith_angle, surf_down_sw, ls_rain,   &
                           ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq,       &
                           z1_uv, PHYS%rhow, L_tile_pts, canopy_tile, Fland,   &
                           CO2_MMR, sthu_tile, smcl_tile, sthf_tile,           &
                           sthu, tsoil_tile, canht_ft, lai_ft,                 &
                           sin_theta_latitude, dzsoil )                         

   canopy%oldcansto=canopy%cansto


      print *, "jhan:cable_expl 2"
   !---------------------------------------------------------------------!
   !--- real(timestep) width, CABLE types passed to CABLE "engine" as ---!  
   !--- req'd by Mk3L  --------------------------------------------------!
   !---------------------------------------------------------------------!
   CALL cbm( timestep, air, bgc, canopy, met, bal,                             &
             rad, rough, soil, ssnow, sum_flux, veg )

      print *, "jhan:cable_expl 3 pre-unpack"



   !---------------------------------------------------------------------!
   !--- pass land-surface quantities calc'd by CABLE in explicit call ---!
   !--- back to UM.                                                   ---!
   !---------------------------------------------------------------------!
   call cable_expl_unpack( FTL_TILE, FQW_TILE,          &
                           TSTAR_TILE,     &
                           U_S, U_S_STD_TILE,       &
                           CD_TILE, &
                           !CH_TILE, &
                           FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           , EPOT_TILE, l_tile_pts,             &
                           ssnow%snowd, ssnow%cls, air%rlam, air%rho,          &
                           canopy%fe, canopy%fh, canopy%us, canopy%cdtq,       &
                           canopy%fwet, canopy%wetfac_cs, canopy%rnet,         &
                           canopy%zetar, canopy%epot, met%ua, rad%trad,        &
                           rad%transd, rough%z0m, rough%zref_tq )


      print *, "jhan:cable_expl 4"
   ! dump bitwise reproducible testing data
   IF( cable_user%RUN_DIAG_LEVEL == 'zero')                                    &
      call cable_diag( 1, "FLUXES", mp, kend_gl, ktau_gl, knode_gl,            &
                          "FLUXES", canopy%fe + canopy%fh )
                

   cable_runtime%um_explicit = .FALSE.

      print *, "jhan:cable_expl 5"

END SUBROUTINE cable_explicit_driver




!---------------------------------------------------------------------!
!--- pass land-surface quantities calc'd by CABLE in explicit call ---!
!--- back to UM.                                                   ---!
!---------------------------------------------------------------------!
SUBROUTINE cable_expl_unpack( FTL_TILE, FQW_TILE,       &
                           TSTAR_TILE,     &
                           U_S, U_S_STD_TILE,       &
                           CD_TILE, FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           EPOT_TILE, l_tile_pts,             &
                           ssnow_snowd, ssnow_cls, air_rlam, air_rho,          &
                           canopy_fe, canopy_fh, canopy_us, canopy_cdtq,       &
                           canopy_fwet, canopy_wetfac_cs, canopy_rnet,         &
                           canopy_zetar, canopy_epot, met_ua, rad_trad,        &
                           rad_transd, rough_z0m, rough_zref_tq )

   USE cable_diag_module
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
   REAL, DIMENSION(um1%land_pts) ::   &
      FTL_CAB, &
      LE_CAB
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      FTL_TILE,   &  ! Surface FTL for land tiles     
      FQW_TILE       ! Surface FQW for land tiles     
   
   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: &
      FTL_TILE_CAB, &
      LE_TILE_CAB

   !___return temp and roughness
   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: &
      TSTAR_TILE_CAB, TSTAR_TILE,  Z0H_TILE, Z0M_TILE
   REAL, DIMENSION(um1%land_pts) ::                  &
      TSTAR_CAB

   !___return friction velocities/drags/ etc
   !REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: &
      CD_TILE,    &     ! Drag coefficient
      CH_TILE,    &     ! Transfer coefficient for heat & moisture
      U_S_STD_TILE      ! Surface friction velocity
   REAL, INTENT(OUT), DIMENSION(um1%row_length,um1%rows)  :: &
      U_S               ! Surface friction velocity (m/s)
   !REAL, INTENT(OUT), DIMENSION(um1%land_pts) ::                  &
   REAL, DIMENSION(um1%land_pts) ::                  &
      CH_CAB,  &  ! Turbulent surface exchange
      CD_CAB,  &  ! Turbulent surface exchange
      U_S_CAB     ! Surface friction velocity (m/s)

   !___return miscelaneous 
   REAL, DIMENSION(um1%land_pts,um1%ntiles) :: &
      RECIP_L_MO_TILE  ! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      RADNET_TILE,   & ! Surface net radiation
      RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
                       ! factor for fraction (1-FRACA) of snow-free land tiles
      RESFT,         & ! Total resistance factor.
                       ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,
                       ! 1 for snow.    
      FRACA,         & ! Fraction of surface moisture
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

   !___vars in local calc of Surface friction velocities
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
      CD_CAB_TILE,   &  
      CH_CAB_TILE,   &  ! (bulk transfer) coeff. for momentum
      U_S_TILE
   REAL, DIMENSION(mp)  :: &
      CDCAB,CHCAB

   !___local miscelaneous
   REAL, DIMENSION(mp)  :: &
   THETAST,fraca_cab,rfsfs_cab, RECIPLMOTILE, fe_dlh
   INTEGER :: i,j,k,N,L
   REAL :: miss = 0.0
   LOGICAL, SAVE :: first_cable_call = .true.
   REAL, POINTER :: CAPP 
   
   character(len=30), dimension(:), allocatable :: explUnpackCheckNames 
   real, dimension(:,:), allocatable :: explUnpackCheckFields

      CAPP => PHYS%CAPP

     print*, "jhan:cable_expl_unpack 1" 

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

      !CALL cable_farray( mp, explUnpackCheckNames, explUnpackCheckFields,      &
      !                   "canopy%fh", canopy_fh,        &                        
      !                   "canopy%us", canopy_us,        &                        
      !                   "met%ua", met_ua,        &                        
      !                   "cdcab", cdcab,        &                        
      !                   "canopy%fe", canopy_fe ) !,        &                        
      
      !CALL cable_NaN( explUnpackCheckNames, explUnpackCheckFields, knode_gl )
                         
      
     print*, "jhan:cable_expl_unpack 2" 

      U_S_STD_TILE=U_S_TILE
     print*, "jhan:cable_expl_unpack 3" 

      CD_TILE = CD_CAB_TILE
     print*, "jhan:cable_expl_unpack 4" 

      CH_TILE = CH_CAB_TILE
     print*, "jhan:cable_expl_unpack 5" 


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
      

      IF(first_cable_call) THEN 
         l_tile_pts = um1%l_tile_pts
         first_cable_call = .FALSE.
      ENDIF

   
END SUBROUTINE cable_expl_unpack
    
!============================================================================
!============================================================================
!============================================================================

