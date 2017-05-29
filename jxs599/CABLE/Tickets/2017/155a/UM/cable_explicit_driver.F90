!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
! 
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located 
! in each directory containing CABLE code.
!
! ==============================================================================
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
  
 
module cable_explicit_driv_mod
  
contains

SUBROUTINE cable_explicit_driver( row_length, rows, land_pts, ntiles,npft,     &
                                  sm_levels, timestep, latitude, longitude,    &
                                  land_index, tile_frac,  tile_pts, tile_index,&
                                  bexp, hcon, satcon, sathh, smvcst,           &
                                  smvcwt,  smvccl, albsoil,                    &
                                  snow_tile,                                   &
                                  snow_rho1l,                                  &
                                  snow_age,                                    &
                                  snow_flg3l,                                  &
                                  snow_rho3l,                                  &
                                  !snow_cond, 
                                  snow_depth3l,                                &
                                  snow_tmp3l,                                  &
                                  snow_mass3l,                                 & 
                                  !sw_down, 
                                  lw_down,                                     &
                                  cos_zenith_angle,                            &
                                  surf_down_sw,                                &
                                  ls_rain,                                     &
                                  ls_snow,                                     &
                                  tl_1, qw_1, &
                                  vshr_land, pstar, z1_tq, z1_uv,  &
                                  canopy_tile,                                 &
                                  Fland,                                       &
                                  CO2_MMR,                                     & 
                                  ! r935 rml 2/7/13 pass 3d co2 through to cable if required
                                  !CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,   &
                                  !
                                  !sthu_tile, 
                                  smcl_tile,                                   &
                                  sthf_tile,                                   &
                                  sthu,                                        &
                                  tsoil_tile,                                  &
                                  canht_ft,                                    &
                                  lai_ft, sin_theta_latitude, dzsoil,          &
                                  FTL_TILE, FQW_TILE,                          &
                                  TSTAR_TILE,                                  &
                                  U_S, U_S_STD_TILE,                           &
                                  CD_TILE, CH_TILE,                            &
                                  RADNET_TILE, FRACA,                          &
                                  RESFS, RESFT,                                &
                                  Z0H_TILE,Z0M_TILE,                           &
                                  RECIP_L_MO_TILE, EPOT_TILE,                  &
                                  ! r825 adds CASA vars here
                                  !CPOOL_TILE, NPOOL_TILE, PPOOL_TILE,          &
                                  !SOIL_ORDER, NIDEP, NIFIX, PWEA, PDUST,       &
                                  !GLAI, PHENPHASE, NPP_FT_ACC, RESP_W_FT_ACC,  &
                                  endstep, timestep_number, mype )    
   
   !--- reads runtime and user switches and reports
   USE cable_um_tech_mod, ONLY : cable_um_runtime_vars, air, bgc, canopy,      &
                                 met, bal, rad, rough, soil, ssnow, sum_flux,  &
                                 veg, basic_diag 
   
   !--- vars common to CABLE declared 
   USE cable_common_module!, ONLY : cable_runtime, cable_user, ktau_gl,         &
                          !         knode_gl, kwidth_gl, kend_gl,               &
                          !         report_version_no,                          & 
                          !         l_vcmaxFeedbk, l_laiFeedbk
   
   !--- subr to (manage)interface UM data to CABLE
   USE cable_um_init_mod, ONLY : interface_UM_data
   
   !--- subr to call CABLE model
   USE cable_cbm_module, ONLY : cbm

   USE cable_def_types_mod, ONLY : mp, ms

   USE cable_expl_unpack_mod, ONLY : cable_expl_unpack
  
  USE cable_decs_mod, only : L_tile_pts, rho_water

   USE casavariable
   USE casa_types_mod
  !fprintf
  !USE cable_common_module, ONLY :    fprintf_dir_root, fprintf_dir
  USE cable_diag_module
  USE cable_climate_mod
  !made this a module so can build in the UM re:dependencies etc
  !did the same dor sli_main. promote everything to modules
  USE casa_cable
   
   IMPLICIT NONE
 
  character(len=*), parameter :: subr_name = "cable_explicit_driver"

   !jhan: this can be moved and USEd ?
   TYPE (climate_type)	:: climate     ! climate variables
   !-------------------------------------------------------------------------- 
   !--- INPUT ARGS FROM () --------------------------------------------
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
      snow_flg3l   ! 3 layer snow flag

   !___UM parameters: water density, soil layer thicknesses 
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
      snow_age
   
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
! rml 2/7/13 Extra atmospheric co2 variables
   !LOGICAL, INTENT(IN) :: L_CO2_INTERACTIVE
   !INTEGER, INTENT(IN) ::                              &
   !  CO2_DIM_LEN                                      &
   !  ,CO2_DIM_ROW
   !REAL, INTENT(IN) :: CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)  ! co2 mass mixing ratio
  
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
   !INTEGER:: itimestep
    
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

   !r825 adds CASA vars here
   REAL, DIMENSION(land_pts,ntiles,10) :: &
      CPOOL_TILE,    & ! Carbon Pools
      NPOOL_TILE       ! Nitrogen Pools

   REAL, DIMENSION(land_pts,ntiles,12) :: &
      PPOOL_TILE       ! Phosphorus Pools

   REAL, DIMENSION(land_pts) :: &
      SOIL_ORDER,    & ! Soil Order (1 to 12)
      NIDEP,         & ! Nitrogen Deposition
      NIFIX,         & ! Nitrogen Fixation
      PWEA,          & ! Phosphorus from Weathering
      PDUST            ! Phosphorus from Dust

   REAL, DIMENSION(land_pts,ntiles) :: &
      GLAI, &          ! Leaf Area Index for Prognostics LAI
      PHENPHASE        ! Phenology Phase for Casa-CNP
                                  
   REAL, DIMENSION(land_pts,ntiles) :: &
      NPP_FT_ACC,     &
      RESP_W_FT_ACC
     
   !-------------------------------------------------------------------------- 
   !--- end INPUT ARGS FROM sf_exch() ----------------------------------------
   !-------------------------------------------------------------------------- 
   
   !___ declare local vars 
   
   !___ location of namelist file defining runtime vars
   CHARACTER(LEN=200), PARAMETER ::                                            & 
      runtime_vars_file = 'cable.nml'


   !___ 1st call in RUN (!=ktau_gl -see below) 
   LOGICAL, SAVE :: first_cable_call = .TRUE.
 
   !___ unique unit/file identifiers for cable_diag: arbitrarily 5 here 
   INTEGER, SAVE :: iDiagZero=0, iDiag1=0, iDiag2=0, iDiag3=0, iDiag4=0

   ! Vars for standard for quasi-bitwise reproducability b/n runs
   ! Check triggered by cable_user%consistency_check = .TRUE. in cable.nml
   CHARACTER(len=30), PARAMETER ::                                             &
      Ftrunk_sumbal  = ".trunk_sumbal",                                        &
      Fnew_sumbal    = "new_sumbal"

   DOUBLE PRECISION, save ::                                                         &
      trunk_sumbal = 0.0, & !
      new_sumbal = 0.0

   INTEGER :: ioerror=0
integer ::  i,j

  !fprintf{_____________________________________________________________________  
  !USE cable_common_module, ONLY :    fprintf_dir_root, fprintf_dir
  !USE cable_diag_module
  character(len=25) :: vname
  character(len=70) :: dir
 
  INTEGER, SAVE ::                                                            &
   cDiag00=0, cDiag0=0, cDiag1=0, cDiag2=0, cDiag3=0, cDiag4=0,                &
   cDiag5=0, cDiag6=0, cDiag7=0, cDiag8=0, cDiag9=0, cDiag10=0,    cDiag11=0,  &  
   cDiag12=0, cDiag13=0, cDiag14=0, cDiag15=0, cDiag16=0, cDiag17=0, cDiag18=0,& 
   cDiag19=0,cDiag20=0, cDiag21=0, cDiag22=0, cDiag23=0, cDiag24=0, cDiag25=0, &
   cDiag26=0, cDiag27=0, cDiag28=0, cDiag29=0, cDiag30=0, cDiag31=0, cDiag32=0,&  
   cDiag33=0, cDiag34=0, cDiag35=0, cDiag36=0, cDiag37=0, cDiag38=0, cDiag39=0,& 
   cDiag40=0, cDiag41=0, cDiag42=0, cDiag43=0, cDiag44=0, cDiag45=0, cDiag46=0,& 
   cDiag47=0, cDiag48=0, cDiag49=0, cDiag50=0, cDiag51=0, cDiag52=0, cDiag53=0,& 
   cDiag54=0, cDiag55=0, cDiag56=0, cDiag57=0, cDiag58=0, cDiag59=0, cDiag60=0,& 
   cDiag61=0, cDiag62=0, cDiag63=0, cDiag64=0, cDiag65=0, cDiag66=0, cDiag67=0,& 
   cDiag68=0, cDiag69=0

  logical, parameter :: L_fprint_HW = .false.
  logical :: L_fprint
  
  L_fprint = .false. ! default
  
  !if( L_fprint_HW ) then
  !  if ( ktau_gl==1 .OR. ktau_gl==54 .OR. ktau_gl==154 .OR. &
  !       ktau_gl==154 .OR. ktau_gl==154 .OR. mod(ktau_gl,10)==0. ) then
  !    L_fprint = .true.
  !  endif  
  !endif  

  fprintf_dir=trim(fprintf_dir_root)//trim("expl_unpack")//"/"
  
  !vname='' 
  !call cable_fprintf( cDiagX, vname, %, mp, L_fprint )
  !fprintf_____________________________________________________________________} 
 
!hard-wired for now as not reasing namelists
cable_user%cable_runtime_coupled=.FALSE.
cable_user%diag_soil_resp='ON '
cable_user%fwsoil_switch='standard'
cable_user%l_new_reduce_soilevp=.FALSE.
cable_user%l_new_roughness_soil=.FALSE.
cable_user%l_new_runoff_speed=.FALSE.
cable_user%leaf_respiration='OFF'
cable_user%run_diag_level='BASIC'
cable_user%ssnow_potev=''
!casafile%cnpbiome='~/CABLE-AUX/core/biogeochem/pftlookup_csiro_v16_17tiles.csv'
!casafile%phen='~/CABLE-AUX/core/biogeochem/modis_phenology_csiro.txt'
!filename%soil='/home/599/lxs599/CABLE-AUX/core/biogeophys/def_soil_params.txt'
!filename%veg='/home/599/lxs599/CABLE-AUX/core/biogeophys/def_veg_params.txt'
icycle=0
redistrb=.FALSE.
satuparam=0.8
wiltparam=0.5

   !IF(cable_user%run_diag_level == "BASIC")                                    &     
   !   CALL basic_diag(subr_name, "Called.") 

   !--- initialize cable_runtime% switches 
      cable_runtime%um = .TRUE.
   cable_runtime%um_explicit = .TRUE.
   
   !--- UM7.3 latitude is not passed correctly. hack 
   !IF(first_cable_call) latitude = sin_theta_latitude

   !--- user FLAGS, variables etc def. in cable.nml is read on 
   !--- first time step of each run. these variables are read at 
   !--- runtime and for the most part do not require a model rebuild.
!Hardwire these aswell
   IF(first_cable_call) THEN
      !*!CALL cable_um_runtime_vars(runtime_vars_file) 
      first_cable_call = .FALSE.
   ENDIF      

   !---------------------------------------------------------------------!
   !--- initialize CABLE using UM forcings etc. these args are passed ---!
   !--- down from UM.                                                 ---! 
   !---------------------------------------------------------------------!

   CALL interface_UM_data( row_length, rows, land_pts, ntiles, npft,           & 
                           sm_levels, ktau_gl, &
                           latitude, &
                           longitude,         &
                           land_index, tile_frac, tile_pts, tile_index,        &
                           bexp, hcon, satcon, sathh, smvcst, smvcwt,          &
                           smvccl, albsoil, snow_tile, snow_rho1l,             &
                           snow_age, snow_flg3l, snow_rho3l, snow_cond,     &
                           snow_depth3l, snow_tmp3l, snow_mass3l, sw_down,     &
                           lw_down, cos_zenith_angle, surf_down_sw, ls_rain,   &
                           ls_snow, tl_1, qw_1, vshr_land, pstar, z1_tq,       &
                           z1_uv, rho_water, L_tile_pts, canopy_tile, Fland,   &
                           CO2_MMR, &
! r935 rml 2/7/13 pass 3d co2 through to cable if required
                   !CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,   &
                           sthu_tile, smcl_tile, sthf_tile,                    &
                           sthu, tsoil_tile, canht_ft, lai_ft,                 &
                           sin_theta_latitude, dzsoil )!,                         &
                           ! r825	
                           !CPOOL_TILE, NPOOL_TILE, PPOOL_TILE, SOIL_ORDER,     &
                           !NIDEP, NIFIX, PWEA, PDUST, GLAI, PHENPHASE,         &
                           !NPP_FT_ACC,RESP_W_FT_ACC )

  fprintf_dir=trim(fprintf_dir_root)//trim("expl_driver")//"/"

   !---------------------------------------------------------------------!
   !--- Feedback prognostic vcmax and daily LAI from casaCNP to CABLE ---!
   !---------------------------------------------------------------------!
   !IF(l_vcmaxFeedbk) call casa_feedback(ktau_gl,veg,casabiome,casapool,casamet)
   !IF(l_laiFeedbk) veg%vlai(:) = casamet%glai(:)

   canopy%oldcansto=canopy%cansto

   IF(cable_user%run_diag_level == "BASIC")                                    &     
      CALL basic_diag(subr_name, "call cbm.") 
   !---------------------------------------------------------------------!
   !--- real(timestep) width, CABLE types passed to CABLE "engine" as ---!  
   !--- req'd by Mk3L  --------------------------------------------------!
   !---------------------------------------------------------------------!
   CALL cbm( ktau_gl,timestep, air, bgc, canopy, met, bal,                             &
             rad, rough, soil, ssnow, sum_flux, veg, climate )

   !---------------------------------------------------------------------!
   ! Check this run against standard for quasi-bitwise reproducability   !  
   ! Check triggered by cable_user%consistency_check=.TRUE. in cable.nml !
   !---------------------------------------------------------------------!
   IF(cable_user%consistency_check) THEN 
         
      IF( knode_gl==1 ) &
         new_sumbal = new_sumbal + ( SUM(canopy%fe) + SUM(canopy%fh)           &
                    + SUM(ssnow%wb(:,1)) + SUM(ssnow%tgg(:,1)) )
     
      IF( knode_gl==1 .and. ktau_gl==kend_gl ) then 
         
         IF( abs(new_sumbal-trunk_sumbal) < 1.e-7 ) THEN
   
            print *, ""
            print *, &
            "Internal check shows this version reproduces the trunk sumbal"
         
         ELSE
   
            print *, ""
            print *, &
            "Internal check shows in this version new_sumbal != trunk sumbal"
            print *, "The difference is: ", new_sumbal - trunk_sumbal
            print *, &
            "Writing new_sumbal to the file:", TRIM(Fnew_sumbal)
                  
            OPEN( 12, FILE = Fnew_sumbal )
               WRITE( 12, '(F20.7)' ) new_sumbal  ! written by previous trunk version
            CLOSE(12)
         
         ENDIF   
      
      ENDIF   
   
   ENDIF


   !---------------------------------------------------------------------!
   !--- pass land-surface quantities calc'd by CABLE in explicit call ---!
   !--- back to UM.                                                   ---!
   !---------------------------------------------------------------------!

   call cable_expl_unpack( latitude, longitude, FTL_TILE, FQW_TILE, TSTAR_TILE, &
                           U_S, U_S_STD_TILE, &
                           CD_TILE, CH_TILE, FLAND, RADNET_TILE,       &
                           FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,            &
                           RECIP_L_MO_TILE, EPOT_TILE, l_tile_pts,             &
                           ssnow%snowd, ssnow%cls, air%rlam, air%rho,          &
                           canopy%fe, canopy%fh, canopy%us, canopy%cdtq,       &
                           canopy%fwet, canopy%wetfac_cs, canopy%rnet,         &
                           canopy%zetar, canopy%epot, met%ua, rad%trad,        &
                           rad%transd, rough%z0m, rough%zref_tq )

   cable_runtime%um_explicit = .FALSE.

END SUBROUTINE cable_explicit_driver

End module cable_explicit_driv_mod
