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

SUBROUTINE cable_explicit_driver( &
   row_length, rows, land_pts, ntiles,npft,     &
                                  sm_levels, timestep, latitude, longitude,    &
                                  land_index, tile_frac,  tile_pts, tile_index,&
                                  bexp, hcon, satcon, sathh, smvcst,           &
                                  smvcwt,  smvccl, albsoil, snow_tile,         &
                                  snow_rho1l, snage_tile, isnow_flg3l,         &
                                  snow_rho3l, snow_cond, snow_depth3l,         &
                                  snow_tmp3l, snow_mass3l, sw_down, lw_down,   &
                                  cos_zenith_angle, surf_down_sw, ls_rain,     &
                                  ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq,&
                                  z1_uv, rho_water, L_tile_pts, canopy_tile,   &
                                  Fland, CO2_MMR, sthu_tile, smcl_tile,        &
                                  sthf_tile, sthu, tsoil_tile, canht_ft,       &
                                  lai_ft, sin_theta_latitude, dzsoil,          &
                                  LAND_MASK, FTL_TILE,  &
                                  FQW_TILE, TSTAR_TILE,   &
                                  U_S, U_S_STD_TILE,&
                                  CD_TILE, CH_TILE,   &
                                  RADNET_TILE, FRACA, rESFS, RESFT, Z0H_TILE,  &
                                  Z0M_TILE, RECIP_L_MO_TILE, EPOT_TILE,        &
                                  endstep, timestep_number, mype & 
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

   !--- include subr called to write data for testing purposes 
   USE cable_diag_module

   IMPLICIT NONE
 
 
  
   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM sf_exch() --------------------------------------------
   !-------------------------------------------------------------------------- 
   
   !___IN: UM dimensions, array indexes, flags
   INTEGER ::                                                      & 
      row_length, rows, & ! UM grid resolution
      land_pts,         & ! # of land points being processed
      ntiles,           & ! # of tiles 
      npft,             & ! # of plant functional types
      sm_levels           ! # of soil layers 

   ! index of land points being processed
   INTEGER, DIMENSION(land_pts) :: land_index 

   ! # of land points on each tile
   INTEGER,  DIMENSION(ntiles) :: tile_pts 
   
   INTEGER,  DIMENSION(land_pts, ntiles) ::                         & 
      tile_index ,& ! index of tile points being processed
      isnow_flg3l   ! 3 layer snow flag

   !--- TRUE if land, F elsewhere.
   !jhan:rm land_mask
   LOGICAL,DIMENSION(row_length,rows) :: land_mask   

   !___UM parameters: water density, soil layer thicknesses 
   REAL :: rho_water 
   REAL,  DIMENSION(sm_levels) :: dzsoil

   !___UM soil/snow/radiation/met vars
   REAL,  DIMENSION(land_pts) :: & 
      bexp,    & ! => parameter b in Campbell equation 
      hcon,    & ! Soil thermal conductivity (W/m/K).
      satcon,  & ! hydraulic conductivity @ saturation [mm/s]
      sathh,   &
      smvcst,  &
      smvcwt,  &
      smvccl,  &
      albsoil, &
      fland 
   
   REAL,  DIMENSION(row_length,rows) :: &
      sw_down,          & 
      cos_zenith_angle
   
   REAL,  DIMENSION(row_length,rows) ::                             &
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

   REAL,  DIMENSION(land_pts, ntiles) ::                         &
      snow_tile

   REAL, DIMENSION(land_pts, ntiles) ::                            &
      tile_frac,  &    
      snow_rho1l, &
      snage_tile
   
   REAL, DIMENSION(row_length, rows, 4) ::                         &
      surf_down_sw 
   
   REAL, DIMENSION(land_pts, npft) ::                              &
      canht_ft, lai_ft 
   
   REAL, DIMENSION(land_pts, ntiles) ::                             &
      canopy_tile
   
   REAL,  DIMENSION(land_pts, ntiles,3) ::                       &
      snow_cond
   
   REAL, DIMENSION(land_pts, ntiles,3) ::                          &
      snow_rho3l,    &
      snow_depth3l,  &
      snow_mass3l,   &
      snow_tmp3l
   
   REAL, DIMENSION(land_pts, sm_levels) ::                         &
      sthu 
   
   REAL, DIMENSION(land_pts, ntiles, sm_levels) :: & 
      sthu_tile, &
      sthf_tile, &
      smcl_tile, &
      tsoil_tile
   
   REAL :: co2_mmr

   !___true IF vegetation (tile) fraction is greater than 0
   LOGICAL, DIMENSION(land_pts, ntiles) :: L_tile_pts
  
   REAL :: sin_theta_latitude(row_length,rows) 
     
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
   REAL, DIMENSION(land_pts,ntiles) :: &
      CD_TILE,    &     ! Drag coefficient
      CH_TILE,    &     ! Transfer coefficient for heat & moisture
      U_S_STD_TILE      ! Surface friction velocity

   REAL, DIMENSION(row_length,rows)  :: &
      U_S               ! Surface friction velocity (m/s)

   ! end step of experiment, this step, step width, processor num
   INTEGER :: endstep, timestep_number, mype
   REAL ::  timestep     
   
   INTEGER:: itimestep
    
   !___return miscelaneous 
   REAL,  DIMENSION(land_pts,ntiles) :: &
      RADNET_TILE,   & ! Surface net radiation
      RESFS,         & ! Combined soil, stomatal & aerodynamic resistance
                       ! factor for fraction (1-FRACA) of snow-free land tiles
      RESFT,         & ! Total resistance factor.
                       ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,        
                       ! 1 for snow.    
      FRACA,         & ! Fraction of surface moisture
      RECIP_L_MO_TILE,& ! Reciprocal of the Monin-Obukhov length for tiles (m^-1)
      EPOT_TILE
     
   !-------------------------------------------------------------------------- 
   !--- end INPUT ARGS FROM sf_exch() ----------------------------------------
   !-------------------------------------------------------------------------- 
   
   

   
   !___ declare local vars 
   
   !___ location of namelist file defining runtime vars
   CHARACTER(LEN=200), PARAMETER ::                                            & 
      runtime_vars_file = 'cable.nml'


   !___ 1st call in RUN (!=ktau_gl -see below) 
   LOGICAL, SAVE :: first_cable_call = .TRUE.
 
   INTEGER, SAVE ::  iDiag0=0,iDiag1=0, iDiag2=0
!umoutput print_args.1
!print *, "row_length ", row_length 
!print *, "rows " , rows 
!print *, "land_pts " , land_pts 
!print *, "ntiles,npft " , ntiles,npft
!print *, "sm_levels " , sm_levels

!print *, "jhan:_explicit: args read"
!umoutput print_args.2
!print *, ""
!print *,"timestep ", timestep
!print *, ""
!print *,"latitude ",  latitude 
!print *, ""
!print *,"longitude ", longitude
!print *, ""
!print *,"land_index ", land_index
!print *, ""
!print *,"tile_frac ",tile_frac
!print *, ""
!print *,"tile_pts ", tile_pts
!print *, ""
!print *,"tile_index ", tile_index
!print *, ""
!print *,"bexp ", bexp
!print *, ""
!print *,"hcon ", hcon
!print *, ""
!print *,"satcon ",satcon
!print *, ""
!print *,"sathh ", sathh
!print *, ""
!print *,"smvcst ", smvcst
!stop
!print *, ""
!print *,"smvcwt ", smvcwt
!print *, ""
!print *,"smvccl ", smvccl
!print *, ""
!print *,"albsoil ", albsoil
!print *, ""
!print *,"snow_tile ", snow_tile
!print *, ""
!print *,"snow_rho1l ", snow_rho1l
!print *, ""
!print *,"snage_tile ", snage_tile
!print *, ""
!print *,"isnow_flg3l ", isnow_flg3l
!print *, ""
!print *,"snow_rho3l ", snow_rho3l
!print *, ""
!print *,"snow_cond ", snow_cond
!print *, ""
!print *,"snow_depth3l ", snow_depth3l
!print *, ""
!print *,"snow_tmp3l ", snow_tmp3l
!print *, ""
!print *,"snow_mass3l ", snow_mass3l
!print *, ""
!print *,"sw_down ", sw_down
!print *, ""
!print *,"lw_down ", lw_down
!print *, ""
!print *,"cos_zenith_angle ", cos_zenith_angle
!print *, ""
!print *,"surf_down_sw ", surf_down_sw
!print *, ""
!print *,"ls_rain ", ls_rain
!print *, ""
!print *,"ls_snow ", ls_snow
!print *, ""
!print *,"tl_1 ", tl_1
!print *, ""
!print *,"qw_1 ", qw_1
!print *, ""
!print *,"vshr_land ", vshr_land
!print *, ""
!print *,"pstar ", pstar
!print *, ""
!print *,"z1_tq ", z1_tq
!print *, ""
!print *,"z1_uv ", z1_uv
!print *, ""
!print *,"rho_water ", rho_water
!print *, ""
!print *,"L_tile_pts ", L_tile_pts
!print *, ""
!print *,"canopy_tile ", canopy_tile
!print *, ""
!print *,"Fland ", Fland
!print *, ""
!print *,"CO2_MMR ", CO2_MMR
!print *, ""
!print *,"sthu_tile ", sthu_tile
!print *, ""
!print *,"smcl_tile ", smcl_tile
!print *, ""
!print *,"sthf_tile ", sthf_tile
!print *, ""
!print *,"sthu ", sthu
!print *, ""
!print *,"tsoil_tile ", tsoil_tile
!print *, ""
!print *,"canht_ft ", canht_ft
!print *, ""
!print *,"lai_ft ", lai_ft
!print *, ""
!print *,"sin_theta_latitude ", sin_theta_latitude
!print *, ""
!print *,"dzsoil ", dzsoil
!print *, ""
!print *,"LAND_MASK ", LAND_MASK
!print *, ""
!print *,"FTL_TILE ", FTL_TILE
!print *, ""
!print *,"FQW_TILE ", FQW_TILE
!print *, ""
!print *,"TSTAR_TILE ", TSTAR_TILE
!print *, ""
!print *,"U_S ", U_S
!print *, ""
!print *,"U_S_STD_TILE ", U_S_STD_TILE
!print *, ""
!print *,"CD_TILE ", CD_TILE
!print *, ""
!print *,"CH_TILE ", CH_TILE
!print *, ""
!print *,"RADNET_TILE ", RADNET_TILE
!print *, ""
!print *,"FRACA ", FRACA
!print *, ""
!print *,"rESFS ", rESFS
!print *, ""
!print *,"RESFT ", RESFT
!print *, ""
!print *,"Z0H_TILE ", Z0H_TILE
!print *, ""
!print *,"Z0M_TILE ", Z0M_TILE
!print *, ""
!print *,"RECIP_L_MO_TILE ", RECIP_L_MO_TILE
!print *, ""
!print *,"EPOT_TILE ", EPOT_TILE
!print *, ""
!print *,"endstep ", endstep
!print *, ""
!print *,"timestep_number ", timestep_number
!print *, ""
!print *,"mype ", mype

   !--- initialize cable_runtime% switches 
   IF(first_cable_call) THEN
      cable_runtime%um = .TRUE.
      write(6,*) ""
      write(6,*) "CABLE_log"
      CALL report_version_no(6) ! wriite revision number to stdout(6)
   ENDIF
      
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

         write (6,*) 'CABLE about to write file'
   IF(first_cable_call) & 
open(unit=7777771,file="/home/599/jxs599/cable.txt",status="unknown", &
  action="write",  form="formatted", &
  position='append' )
         write (7777771,*) ktau_gl 
   !--- user FLAGS, variables etc def. in cable.nml is read on 
   !--- first time step of each run. these variables are read at 
   !--- runtime and for the most part do not require a model rebuild.
   IF(first_cable_call) THEN
      CALL cable_um_runtime_vars(runtime_vars_file) 
      first_cable_call = .FALSE.
   ENDIF      

      

!print *, "jhan:_explicit: pre interface", mype, shape(tile_frac)

   !---------------------------------------------------------------------!
   !--- initialize CABLE using UM forcings etc. these args are passed ---!
   !--- down from UM.                                                 ---! 
   !---------------------------------------------------------------------!
   CALL interface_UM_data( row_length, rows, land_pts, ntiles, npft,           & 
                           sm_levels, itimestep, latitude, longitude,          &
                           land_index, tile_frac, tile_pts, tile_index,        &
                           bexp, hcon, satcon, sathh, smvcst, smvcwt,          &
                           smvccl, albsoil, snow_tile, snow_rho1l,             &
                           snage_tile, isnow_flg3l, snow_rho3l, snow_cond,     &
                           snow_depth3l, snow_tmp3l, snow_mass3l, sw_down,     &
                           lw_down, cos_zenith_angle, surf_down_sw, ls_rain,   &
                           ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq,       &
                           z1_uv, rho_water, L_tile_pts, canopy_tile, Fland,   &
                           CO2_MMR, sthu_tile, smcl_tile, sthf_tile,           &
                           sthu, tsoil_tile, canht_ft, lai_ft,                 &
                           sin_theta_latitude, dzsoil )                         

   canopy%oldcansto=canopy%cansto

!print *, "jhan:_explicit: pre cbm", mype

   !---------------------------------------------------------------------!
   !--- real(timestep) width, CABLE types passed to CABLE "engine" as ---!  
   !--- req'd by Mk3L  --------------------------------------------------!
   !---------------------------------------------------------------------!
   CALL cbm( timestep, air, bgc, canopy, met, bal,                             &
             rad, rough, soil, ssnow, sum_flux, veg )


!print *, "jhan:_explicit: pre unpack"


   !---------------------------------------------------------------------!
   !--- pass land-surface quantities calc'd by CABLE in explicit call ---!
   !--- back to UM.                                                   ---!
   !---------------------------------------------------------------------!
   call cable_expl_unpack( FTL_TILE, FQW_TILE,          &
                           TSTAR_TILE, &
                           U_S, U_S_STD_TILE, &
                           CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts,             &
                           ssnow%snowd, ssnow%cls, air%rlam, air%rho,          &
                           canopy%fe, canopy%fh, canopy%us, canopy%cdtq,       &
                           canopy%fwet, canopy%wetfac_cs, canopy%rnet,         &
                           canopy%zetar, canopy%epot, met%ua, rad%trad,        &
                           rad%transd, rough%z0m, rough%zref_tq )


   ! dump bitwise reproducible testing data
   IF( cable_user%RUN_DIAG_LEVEL == 'zero')                                    &
      call cable_diag( iDiag0, "FLUXES", mp, kend_gl, ktau_gl, knode_gl,            &
                          "FLUXES", canopy%fe + canopy%fh )
                

   cable_runtime%um_explicit = .FALSE.

!print *, "jhan:_explicit: end _explicit"

END SUBROUTINE cable_explicit_driver




!---------------------------------------------------------------------!
!--- pass land-surface quantities calc'd by CABLE in explicit call ---!
!--- back to UM.                                                   ---!
!---------------------------------------------------------------------!
SUBROUTINE cable_expl_unpack( FTL_TILE, FQW_TILE,       &
                           TSTAR_TILE, &
                           U_S, U_S_STD_TILE, &
                           CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts,             &
                           ssnow_snowd, ssnow_cls, air_rlam, air_rho,          &
                           canopy_fe, canopy_fh, canopy_us, canopy_cdtq,       &
                           canopy_fwet, canopy_wetfac_cs, canopy_rnet,         &
                           canopy_zetar, canopy_epot, met_ua, rad_trad,        &
                           rad_transd, rough_z0m, rough_zref_tq )

   USE cable_def_types_mod, ONLY : mp, NITER 
   USE cable_data_module,   ONLY : PHYS
   USE cable_um_tech_mod,   ONLY : um1
   USE cable_common_module, ONLY : cable_runtime, cable_user, &
                                   ktau_gl, knode_gl, kend_gl 
   USE cable_diag_module
   IMPLICIT NONE         


   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
   !-------------------------------------------------------------------------- 


   !___ UM variables to recieve unpacked CABLE vars

   !___return fluxes
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      FTL_TILE,   &  ! Surface FTL for land tiles     
      FQW_TILE      ! Surface FQW for land tiles     

   !___return temp and roughness
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      Z0H_TILE,         &
      Z0M_TILE

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


   
   REAL, INTENT(OUT), DIMENSION(um1%land_pts,um1%ntiles) :: &
      TSTAR_TILE
        
   !___vars in local calc. of latent heat fluxes
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
      LE_TILE

   !___vars in local calc of Surface friction velocities
   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                  &
      U_S_TILE
   REAL, DIMENSION(mp)  :: &
      CDCAB

   !___local miscelaneous
   REAL, DIMENSION(mp)  :: &
   THETAST,fraca_cab,rfsfs_cab, RECIPLMOTILE, fe_dlh
   INTEGER :: i,j,k,N,L
   REAL :: miss = 0.0
   LOGICAL, SAVE :: first_cable_call = .true.
   REAL, POINTER :: CAPP 
   
   INTEGER, SAVE ::  iDiag0=0,iDiag1=0, iDiag2=0
   
      CAPP => PHYS%CAPP
      
      !___return fluxes
      fe_dlh = canopy_fe/(air_rlam*ssnow_cls)
      FTL_TILE = UNPACK(canopy_fh,  um1%l_tile_pts, miss)
      FTL_TILE = FTL_TILE / capp
      FQW_TILE = UNPACK(fe_dlh, um1%l_tile_pts, miss)

      !___return temp and roughness
      TSTAR_TILE = UNPACK(rad_trad,  um1%l_tile_pts, miss)
      Z0M_TILE = UNPACK(rough_z0m,  um1%l_tile_pts, miss)
      Z0H_TILE = Z0M_TILE
      
      !___return friction velocities/drags/ etc
      U_S_TILE  =  UNPACK(canopy_us, um1%l_tile_pts, miss)
      CDCAB = canopy_us**2/met_ua**2   ! met%ua is always above umin = 0.1m/s
      ! for Cable CD*
      CD_TILE =  UNPACK(CDCAB,um1%l_tile_pts, miss)
      ! for Cable CH*
      CH_TILE =  UNPACK(canopy_cdtq,um1%l_tile_pts, miss)

      U_S_STD_TILE=U_S_TILE

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

!cable_diag all of these and plot -actually do the cable (mp) version
       
      call cable_diag( iDiag0, "FTL_tile", mp, kend_gl, ktau_gl, knode_gl,            &
                          "FTL_tile", canopy_fh/capp )
                
!print *, "FTL_TILE ", FTL_TILE

!print *, "FTL_TILE ", FTL_TILE
!print *, "FQW_TILE ", FQW_TILE
!print *, "TSTAR_TILE ", TSTAR_TILE
!print *, "Z0M_TILE ", Z0M_TILE
!print *, "U_S_TILE ", U_S_TILE
!print *, "CD_TILE  ", CD_TILE 
!print *, "CH_TILE ", CH_TILE
!print *, "FRACA ", FRACA
!print *, "RESFT ", RESFT
!print *, "RESFS ", RESFS
!print *, "RADNET_TILE ", RADNET_TILE
!print *, "RECIP_L_MO_TILE ", RECIP_L_MO_TILE
!print *, "EPOT_TILE ", EPOT_TILE
!stop

      IF(first_cable_call) THEN 
         l_tile_pts = um1%l_tile_pts
         first_cable_call = .FALSE.
      ENDIF

   
END SUBROUTINE cable_expl_unpack
    
!============================================================================
!============================================================================
!============================================================================

