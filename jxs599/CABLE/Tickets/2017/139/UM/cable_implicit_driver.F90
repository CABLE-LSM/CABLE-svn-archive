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
! Purpose: Updates CABLE variables (as altered by first pass through boundary 
!          layer and convection scheme), calls cbm, passes CABLE variables back 
!          to UM. 'Implicit' is the second call to cbm in each UM timestep.
!
! Called from: UM code sf_impl
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
! ==============================================================================

module cable_implicit_driv_mod
  
contains

subroutine cable_implicit_driver( i_day_number, &
cycleno,                                                         & 
row_length,rows, land_pts, ntiles, npft, sm_levels,              &
dim_cs1, dim_cs2,                                                & 
Fland,                                                           &
LS_RAIN, CON_RAIN, LS_SNOW, CONV_SNOW,       &
                                  DTL_1,DQW_1, TSOIL, TSOIL_TILE, SMCL,        &
                                  SMCL_TILE, timestep, SMVCST,STHF, STHF_TILE, &
                                  STHU,&
!STHU_TILE, &
snow_tile, SNOW_RHO1L,       &
                                  ISNOW_FLG3L, SNOW_DEPTH3L, SNOW_MASS3L,      &
                                  SNOW_RHO3L, SNOW_TMP3L, &
!SNOW_COND,           &
                                  FTL_1, FTL_TILE, FQW_1, FQW_TILE,    &
                                  TSTAR_TILE, &
                                  SURF_HT_FLUX_LAND, ECAN_TILE, ESOIL_TILE,    &
                                  EI_TILE, RADNET_TILE, &
!TOT_ALB, 
                                  SNOW_AGE,   &
                                  CANOPY_TILE, GS,GS_TILE,T1P5M_TILE, Q1P5M_TILE,&
                                  CANOPY_GB, MELT_TILE, &
                                  NPP, NPP_FT, GPP, GPP_FT, RESP_S,   &
                                  RESP_S_TOT, &
                                  RESP_S_TILE, &
                                  RESP_P, RESP_P_FT,  &
                     ! r825 added casa vars after G_LEAF, but we need
                    ! need. vars here satisfy _hydrol CALL that is now from _impl
                    !idoy added r1164+
                                  G_LEAF, & 
!LYING_SNOW, SURF_ROFF, SUB_SURF_ROFF,  &
!TOT_TFALL, 
TL_1, QW_1, SURF_HTF_TILE, & 
                                  !TOT_TFALL )
                                  !G_LEAF, TRANSP_TILE, 
                                  CPOOL_TILE, NPOOL_TILE, PPOOL_TILE,  &
                                  GLAI, PHENPHASE, & 
                                  NPP_FT_ACC, RESP_W_FT_ACC )

  USE cable_implicit_unpack_mod, ONLY : implicit_unpack
   USE cable_def_types_mod, ONLY : mp, msn
   USE cable_data_module,   ONLY : PHYS
   !USE cable_data_mod,      ONLY : cable
   USE cable_um_tech_mod,   ONLY : um1, conv_rain_prevstep, conv_snow_prevstep,&
                                  air, bgc, canopy, met, bal, rad, rough, soil,&
                                  ssnow, sum_flux, veg, basic_diag
   USE cable_common_module, ONLY : cable_runtime, cable_user, l_casacnp,       &
                                   l_vcmaxFeedbk, knode_gl, ktau_gl, kend_gl
  
   USE cable_um_init_subrs_mod, ONLY : um2cable_rr
   USE cable_cbm_module,    ONLY : cbm
   USE casavariable
   USE phenvariable
   USE casa_types_mod
   USE casa_cable, only : bgcdriver, sumcflux
   USE casa_um_inout_mod
  USE cable_common_module, ONLY :    fprintf_dir_root, fprintf_dir
  USE cable_diag_module
   USE cable_climate_mod
   use POP_TYPES, only : pop_type

   IMPLICIT NONE

   TYPE (climate_type)  :: climate     ! climate variables
  !necessary as arg checking is enforce in modular structure that now present 
  ! - HOWEVER *NB*  this POP is not initialized anywhere
  TYPE(POP_TYPE) :: POP 
  
  character(len=*), parameter :: subr_name = "cable_implicit_driver"
        
  integer :: cycleno
  integer :: row_length,rows, land_pts, ntiles, npft, sm_levels
  integer :: dim_cs1, dim_cs2 
  
  REAL,  DIMENSION(land_pts) :: & 
      fland       ! IN Land fraction on land tiles
   
   REAL, DIMENSION(ROW_LENGTH,ROWS) ::                                 &
      LS_RAIN,  & ! IN Large scale rain
      LS_SNOW,  & ! IN Large scale snow
      CON_RAIN, & ! IN Convective rain
      CONV_SNOW,& ! IN Convective snow
      TL_1,     & !
      QW_1,     & !
      DTL_1,    & ! IN Level 1 increment to T field 
      DQW_1       ! IN Level 1 increment to q field 

   REAL :: timestep

   REAL, DIMENSION(land_pts) ::                                            &
      GS,      &  ! OUT "Stomatal" conductance to
      SMVCST
   
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
      SURF_HTF_TILE, &
      !___Surface FTL, FQL for land tiles
      FTL_TILE, FQW_TILE,                                                   &
      
      !___(tiled) latent heat flux, melting, stomatatal conductance
     LE_TILE, MELT_TILE, GS_TILE,                                           &
     
     !___ INOUT Surface net radiation on tiles (W/m2)
     RADNET_TILE, &
     TOT_ALB,     & ! total albedo
     EI_TILE,     & ! OUT EI for land tiles.
     ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
     ESOIL_TILE     ! evapotranspiration from soil moisture store (kg/m2/s) 

   !___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
   !___ runoff ??
   REAL, dimension(land_pts,sm_levels) ::                           &
      SMCL,       & ! 
      STHF,       & !
      STHU,       & !
      TSOIL         !

   !___(tiled) soil prognostics: as above 
   REAL, dimension(land_pts,ntiles,sm_levels) ::                &
      SMCL_TILE, & !
      STHU_TILE, & !
      TSOIL_TILE,& !
      STHF_TILE    !

   !___flag for 3 layer snow pack
   INTEGER :: ISNOW_FLG3L(LAND_PTS,NTILES)
   
   !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
   REAL, dimension(land_pts,ntiles,3) :: &
      SNOW_DEPTH3L,  & ! 
      SNOW_MASS3L,   & !
      SNOW_RHO3L,    & !
      SNOW_TMP3L,    & !
      SNOW_COND        !
  !NOT declared as Ntiles in UM yet: reverted
   !REAL, DIMENSION(land_pts,ntiles) ::                              &
   REAL, DIMENSION(land_pts,ntiles) ::                              &
      RESP_P_FT,     &
      G_LEAF,        &
!     FRS_TILE,   & ! Local
!     NEE_TILE,   & ! Local
!     NPP_TILE,   & ! Local
!     GPP_TILE,   & ! Local
!     GLEAF_TILE, & ! Local, kdcorbin, 10/10
!     FRP_TILE,   &
      NPP_FT,     &
!     NPP_FT_old, &
      GPP_FT      
!     GPP_FT_old       

   REAL, DIMENSION(land_pts) ::                                         &
      SNOW_GRD,    & !
      CANOPY_GB,   & !
      RESP_P,      & !
      NPP,         & !
      GPP            !
      
   REAL, DIMENSION( land_pts,ntiles ) ::                               &
      SNOW_TILE,     &
      SNOW_RHO1L,    &  ! Mean snow density
      SNOW_AGE,    &
      CANOPY_TILE,   &
      T1P5M_TILE,    &
      Q1P5M_TILE,    &
      TSTAR_TILE,    &
      RESP_S_TILE,   & 
!     RESP_P_FT_old, &
      TRANSP_TILE

   REAL ::                                                                     &
      RESP_S(LAND_PTS,DIM_CS1),     &
!     RESP_S_old(LAND_PTS,DIM_CS1), &
      RESP_S_TOT(DIM_CS2)    

   REAL, DIMENSION(LAND_PTS,NTILES,10) ::                              &
      CPOOL_TILE, &
      NPOOL_TILE     
   REAL, DIMENSION(LAND_PTS,NTILES,12) ::                              &
      PPOOL_TILE
   REAL, DIMENSION(LAND_PTS,NTILES) ::                                 &
      GLAI, &
   !INTEGER, DIMENSION(LAND_PTS,NTILES) ::                              &
      PHENPHASE

   ! Lestevens 23apr13
   REAL, DIMENSION(LAND_PTS,NTILES) ::                                 &
      NPP_FT_ACC, &
      RESP_W_FT_ACC

   INTEGER ::     &
      ktauday,      & ! day counter for CASA-CNP
      i_day_number, & ! day of year (1:365) counter for CASA-CNP
      idoy            ! day of year (1:365) counter for CASA-CNP
   INTEGER, SAVE :: &
      kstart = 1

   REAL, DIMENSION(mp) ::                                                      & 
      dtlc, & 
      dqwc
   
   REAL, DIMENSION(LAND_PTS) ::                               &
      LYING_SNOW,    & ! OUT Gridbox snowmass (kg/m2)        
      SUB_SURF_ROFF, & !
      SURF_ROFF,     & !
      TOT_TFALL        !

   REAL, POINTER :: TFRZ

   !This is a quick fix. These can be organised through namelists
   logical :: spinup=.false., spinconv=.false.,                   &
              dump_read=.false., dump_write=.false.
   integer :: loy=365, lalloc=0
   
   !___ 1st call in RUN (!=ktau_gl -see below) 
   LOGICAL, SAVE :: first_cable_call = .TRUE.
   REAL, ALLOCATABLE:: fwork(:,:,:)

   type ProgBank
      
      real, dimension(:,:), allocatable ::                &
         TSOIL, SMCL, STHF 
         
      real, dimension(:,:), allocatable ::                &
         SNOW_DEPTH,                                        &
         SNOW_MASS, SNOW_TMP, SNOW_RHO
      
      real, dimension(:), allocatable ::                &
         SNOW_RHO1L, SNOW_AGE, SNOW_TILE,&
         OCANOPY
      
      integer, dimension(:), allocatable ::                &
         SNOW_FLG3L

      real, dimension(:), allocatable ::                &
         fes_cor,fhs_cor

      real, dimension(:), allocatable ::                &
         osnowd,owetfac,otss

   End type ProgBank

   !integer :: cpb
   integer, parameter :: cpb =2
   type (ProgBank), dimension(cpb), save :: PB
   
   integer :: i2, j2, n2, l2, m2
   integer :: ipb, i, j, k , n, l, m
   integer :: flpt,ftile,flev
   integer :: lpt,tile,lev, itile
   !decs to write text files 
   character(len=*), parameter :: hfmt1 = '("tile= ", I2, "    ", ES8.2)'
   character(len=30) :: chnode, chipb, chlev
   !character(len=25), dimension(9) :: filename
   integer, dimension(9) :: iunit 
   character(len=25), dimension(9) :: basename
  
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

  L_fprint = .false. !default

  IF(cable_user%run_diag_level == "fprint")                                    &     
  fprintf_dir=trim(fprintf_dir_root)//trim("impl_driver")//"/"
  
  !if( L_fprint_HW ) then
  !  if ( ktau_gl==1 .OR. ktau_gl==54 .OR. ktau_gl==154 .OR. &
  !       ktau_gl==154 .OR. ktau_gl==154 .OR. mod(ktau_gl,10)==0. ) then
  !    L_fprint = .true.
  !  endif  
  !endif  

  !vname='' 
  !call cable_fprintf( cDiagX, vname, %, mp, L_fprint )
  !fprintf_____________________________________________________________________} 

  if(knode_gl==0) then
    write (6, *) "CABLE_LSM:Start Subr: ", subr_name
    write (6, *) "@ knode_gl: ", knode_gl 
    write (6, *) "@ timestep: ", ktau_gl         
    write (6, *) "@ cycleno: ", cycleno         
    write (6, *) "============================="
  endif
     

   !cpb = cable% um% numcycles
   flpt  = land_pts
   ftile = ntiles
   flev  = sm_levels
   if (.NOT. allocated(PB(1) %tsoil) ) then
      do ipb = 1, cpb 
         allocate( PB(ipb)%TSOIL(mp,sm_levels) )
         allocate( PB(ipb)%SMCL(mp,sm_levels) )
         allocate( PB(ipb)%STHF(mp,sm_levels) )
         allocate( PB(ipb)%snow_depth(mp,3) )
         allocate( PB(ipb)%snow_mass(mp,3) )
         allocate( PB(ipb)%snow_tmp(mp,3) )
         allocate( PB(ipb)%snow_rho(mp,3) )
         allocate( PB(ipb)%snow_rho1l(mp) )
         allocate( PB(ipb)%snow_age(mp) )
         allocate( PB(ipb)%snow_flg3l(mp) )
         allocate( PB(ipb)%snow_tile(mp) )
         allocate( PB(ipb)%ocanopy(mp) )
      enddo
   endif   

   ipb = cycleno

   PB(ipb)%tsoil     = ssnow%tgg
   PB(ipb)%smcl      = ssnow%wb
   PB(ipb)%sthf      = ssnow%wbice
   PB(ipb)%snow_depth= ssnow%sdepth
   PB(ipb)%snow_mass = ssnow%smass
   PB(ipb)%snow_tmp  = ssnow%tggsn
   PB(ipb)%snow_rho  = ssnow%ssdn
   PB(ipb)%snow_rho1l= ssnow%ssdnn
   PB(ipb)%snow_age  = ssnow%snage
   PB(ipb)%snow_flg3l= ssnow%isflag
   PB(ipb)%snow_tile = ssnow%snowd
   PB(ipb)%ocanopy   = canopy%oldcansto

   if (ipb == cpb) then
      ssnow%tgg     = PB(1)%tsoil
      ssnow%wb      = PB(1)%smcl
      ssnow%wbice   = PB(1)%sthf
      ssnow%sdepth  = PB(1)%snow_depth
      ssnow%smass   = PB(1)%snow_mass
      ssnow%tggsn   = PB(1)%snow_tmp
      ssnow%ssdn    = PB(1)%snow_rho
      ssnow%ssdnn   = PB(1)%snow_rho1l
      ssnow%snage   = PB(1)%snow_age
      ssnow%isflag  = PB(1)%snow_flg3l
      ssnow%snowd   = PB(1)%snow_tile
      canopy%oldcansto = PB(1)%ocanopy
    endif

      IF(cable_user%run_diag_level == "BASIC") &     
         CALL basic_diag(subr_name, "Called.") 

      TFRZ => PHYS%TFRZ
   
      ! FLAGS def. specific call to CABLE from UM
      cable_runtime%um_explicit = .FALSE.
      cable_runtime%um_implicit = .TRUE.
   
      dtlc = 0. ; dqwc = 0.

      !--- All these subrs do is pack a CABLE var with a UM var.
      !-------------------------------------------------------------------
      !--- UM met forcing vars needed by CABLE which have UM dimensions
      !---(rowlength,rows)[_rr], which is no good to CABLE. These have to be 
      !--- re-packed in a single vector of active tiles. Hence we use 
      !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
      !--- if the land point is/has an active tile
      !--- generic format:
      !--- um2cable_rr( UM var, default value for snow tile, CABLE var, mask )
      !--- where mask tells um2cable_rr whether or not to use default value 
      !--- for snow tile 
      !-------------------------------------------------------------------

  vname='met_precip' 
  !-----------------------------------------------------------------------
  CALL um2cable_rr( (LS_RAIN+CON_RAIN)*TIMESTEP, met%precip)
  call cable_fprintf( cDiag0, vname, met%precip, mp, L_fprint )
  !-----------------------------------------------------------------------
  
  vname='met_precip_sn' 
  !-----------------------------------------------------------------------
  CALL um2cable_rr( (LS_SNOW+CONV_SNOW)*TIMESTEP, met%precip_sn)
  call cable_fprintf( cDiag1, vname, met%precip_sn, mp, L_fprint )
  !-----------------------------------------------------------------------
      
  if( L_fprint_HW ) then
    if ( ipb == cpb ) L_fprint = .true.
  endif  
  vname='tl_1' 
  !-----------------------------------------------------------------------
  CALL um2cable_rr( TL_1, met%tk)
  call cable_fprintf( cDiag2, vname, met%tk, mp, L_fprint )
  !-----------------------------------------------------------------------
  L_fprint = .false.
      
  vname='met_qv' 
  !-----------------------------------------------------------------------
  CALL um2cable_rr( QW_1, met%qv)
  call cable_fprintf( cDiag3, vname, met%qv, mp, L_fprint )
  !-----------------------------------------------------------------------
      
  if( L_fprint_HW ) then
    if ( ipb == cpb ) L_fprint = .true.
  endif  
  vname='dtlc' 
  !-----------------------------------------------------------------------
  CALL um2cable_rr( dtl_1, dtlc)
  call cable_fprintf( cDiag4, vname, dtlc, mp, L_fprint )
  !-----------------------------------------------------------------------
  L_fprint = .false.
      
  vname='dqwc' 
  !-----------------------------------------------------------------------
  CALL um2cable_rr( dqw_1, dqwc)
  call cable_fprintf( cDiag5, vname,dqwc, mp, L_fprint )
  !-----------------------------------------------------------------------
      
  !--- conv_rain(snow)_prevstep are added to precip. in explicit call
  vname='conv_rain_prevstep' 
  !-----------------------------------------------------------------------
  CALL um2cable_rr( (CON_RAIN)*TIMESTEP, conv_rain_prevstep)
  call cable_fprintf( cDiag6, vname, conv_rain_prevstep, mp, L_fprint )
  !-----------------------------------------------------------------------
  
  vname='conv_snow_prevstep' 
  !-----------------------------------------------------------------------
  CALL um2cable_rr( (CONV_snow)*TIMESTEP, conv_snow_prevstep)
  call cable_fprintf( cDiag7, vname, conv_snow_prevstep, mp, L_fprint )
  !-----------------------------------------------------------------------

  met%precip   =  met%precip + met%precip_sn
  met%tk = met%tk + dtlc
  
  if( L_fprint_HW ) then
    if ( ipb == cpb ) L_fprint = .true.
  endif  
  call cable_fprintf( cDiag8, vname, met%tk, mp, L_fprint )
  L_fprint = .false.

      met%qv = met%qv + dqwc
      met%tvair = met%tk
      met%tvrad = met%tk

      canopy%cansto = canopy%oldcansto

   CALL cbm( ktau_gl,timestep, air, bgc, canopy, met, bal,                             &
             rad, rough, soil, ssnow, sum_flux, veg, climate )

      ! Lestevens - temporary ?
      ktauday = int(24.0*3600.0/TIMESTEP)
      idoy=i_day_number
!      IF(idoy==0) idoy =365
  
!Call CASA-CNP
      if (l_casacnp) then
      CALL bgcdriver(ktau_gl,kstart,kend_gl,timestep,met,ssnow,canopy,veg,soil, &
                     climate,casabiome,casapool,casaflux,casamet,casabal,phen, &
                     pop, spinConv,spinup, ktauday, idoy,loy, dump_read,   &
                     dump_write, LALLOC)
      endif

      CALL sumcflux(ktau_gl,kstart,kend_gl,TIMESTEP,bgc,canopy,soil,ssnow,      &
                    sum_flux,veg,met,casaflux,l_vcmaxFeedbk)

     CALL implicit_unpack(& 
                            cycleno,                                         & 
                            row_length,rows, land_pts, ntiles, npft, sm_levels, &
                            dim_cs1, dim_cs2,                                   & 
                            TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,                &
                           SMVCST, STHF, STHF_TILE, STHU, STHU_TILE,          &
                           snow_tile, SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L,  &
                           SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,    &
                           FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
                           SURF_HT_FLUX_LAND, ECAN_TILE,        &
                           ESOIL_TILE, EI_TILE, RADNET_TILE, TOT_ALB,         &
                           SNOW_AGE, CANOPY_TILE, GS,GS_TILE, T1P5M_TILE,        &
                           Q1P5M_TILE, CANOPY_GB, FLAND, MELT_TILE,           &
                           NPP, NPP_FT, GPP, GPP_FT, RESP_S,         &
                           RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT, G_LEAF,& 
                           TRANSP_TILE, NPP_FT_ACC, RESP_W_FT_ACC,SURF_HTF_TILE )

! Call CASA-CNP
      if (l_casacnp) then
        if (knode_gl==0 .and. ktau_gl==kend_gl) then
         print *, '  '; print *, 'CASA_log:'
         print *, '  Calling CasaCNP - Poolout '
         print *, '  l_casacnp = ',l_casacnp
         print *, '  ktau_gl, kend_gl = ',ktau_gl,kend_gl
         print *, 'End CASA_log:'; print *, '  '
        endif
       CALL casa_poolout_unpk(casapool,casaflux,casamet,casabal,phen,  &
                              CPOOL_TILE,NPOOL_TILE,PPOOL_TILE, &
                              GLAI,PHENPHASE)
      endif

      cable_runtime%um_implicit = .FALSE.

  return

END SUBROUTINE cable_implicit_driver


!========================================================================= 
!========================================================================= 
!========================================================================= 
        

End module cable_implicit_driv_mod
