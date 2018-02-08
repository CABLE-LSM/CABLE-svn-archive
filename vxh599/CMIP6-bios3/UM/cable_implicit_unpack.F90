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
! Purpose: Updates CABLE variables (as altered by first pass through boundary 
!          layer and convection scheme), calls cbm, passes CABLE variables back 
!          to UM. 'Implicit' is the second call to cbm in each UM timestep.
!
! Called from: UM/JULES sf_impl
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
! ==============================================================================

module cable_implicit_unpack_mod
  
contains

! rml 30/10/15 rml: rewrite carbon fluxes - previous version seemed
! unnecessarily complicated by allowing for fluxes from sea-ice when don't
! need to account for any carbon flux not from land fraction.
SUBROUTINE implicit_unpack( & 
                            cycleno,                                                         & 
                            row_length,rows, land_pts, ntiles, npft, sm_levels,              &
                            dim_cs1, dim_cs2,                                                & 
                            TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,                &
                            SMVCST, STHF, STHF_TILE, STHU, STHU_TILE,          &
                            snow_tile, SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L,  &
                            SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,    &
                            FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
                            SURF_HT_FLUX_LAND, ECAN_TILE,        &
                            ESOIL_TILE, EI_TILE, RADNET_TILE, TOT_ALB,         &
                            SNOW_AGE, CANOPY_TILE, GS, gs_tile, T1P5M_TILE,  &
                            !args increased following merge
                            Q1P5M_TILE, CANOPY_GB, FLAND, MELT_TILE, &
                            NPP, NPP_FT, GPP, GPP_FT, RESP_S,         &
                            RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT, G_LEAF,&
                            TRANSP_TILE, NPP_FT_ACC, RESP_W_FT_ACC,SURF_HTF_TILE )
 
   USE cable_def_types_mod, ONLY : mp
   USE cable_data_module,   ONLY : PHYS
   USE cable_um_tech_mod,   ONLY : um1 ,canopy, rad, soil, ssnow, air,         &
                                   basic_diag, veg
   USE cable_common_module!, ONLY : cable_runtime, cable_user, fudge_out,       &
                          !         L_fudge, ktau_gl

  USE cable_decs_mod, ONLY : L_tile_pts, rho_water

  USE cable_um_tech_mod,   ONLY : um1
  !fprintf
  USE cable_common_module, ONLY :    fprintf_dir_root, fprintf_dir
  USE cable_diag_module
use arraydiag_m

!C!use cable_implicit_nunpack_mod
   
   IMPLICIT NONE
 
   character(len=*), parameter :: subr_name = "cable_implicit_unpack"
        
  integer :: cycleno
  integer :: row_length,rows, land_pts, ntiles, npft, sm_levels
  integer :: dim_cs1, dim_cs2 

   REAL, DIMENSION(land_pts) ::                                            &
      GS,         &  ! OUT "Stomatal" conductance to
      SMVCST,     &  ! IN Volumetric saturation point
      FLAND          ! IN Land fraction on land tiles
   
   real, dimension(ROW_LENGTH,ROWS) ::                                 &
      !--- Net downward heat flux at surface over land.
      !--- fraction of gridbox (W/m2).
      SURF_HT_FLUX_LAND,           &
      !--- Moisture flux between layers. (kg/m^2/sec).
      !--- FQW(,1) is total water flux from surface, 'E'.
      FQW_1,       &  
      !--- FTL(,K) =net turbulent sensible heat flux into layer K
      !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
      FTL_1         

   REAL, DIMENSION(land_pts,ntiles) ::                                 &
      SURF_HTF_TILE,&
      !___Surface FTL, FQL for land tiles
      FTL_TILE, FQW_TILE,                 &  
      !___(tiled) latent heat flux, melting, stomatatal conductance
      LE_TILE, MELT_TILE, GS_TILE,     &  
      RADNET_TILE, & ! INOUT Surface net radiation on tiles (W/m2)
      TOT_ALB,     & ! total albedo
      EI_TILE,     & ! OUT EI for land tiles.
      ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
      ESOIL_TILE     ! evapotranspiration from soil moisture store (kg/m2/s) 

   !___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
   !___ runoff ??
   REAL, DIMENSION(land_pts,sm_levels) ::                              &
      SMCL,       & !
      STHF,       &
      STHU,       &
      TSOIL       

   !___(tiled) soil prognostics: as above 
   REAL, DIMENSION(land_pts,ntiles,sm_levels) ::                   &
      SMCL_TILE,  & 
      STHU_TILE,  &
      TSOIL_TILE, &
      STHF_TILE  

   !___flag for 3 layer snow pack
   INTEGER :: ISNOW_FLG3L(LAND_PTS,NTILES)
   
   !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
   REAL, DIMENSION(land_pts,ntiles,3) ::                               &
      SNOW_DEPTH3L,  &
      SNOW_MASS3L,   &
      SNOW_RHO3L,    &
      SNOW_TMP3L,    &
      SNOW_COND 

!   REAL, dimension(land_pts,ntiles) ::                                 &
!     FRS_TILE,   & ! Local
!     NEE_TILE,   & ! Local
!     NPP_TILE,   & ! Local
!     GPP_TILE,   & ! Local
!     SURF_HTF_T_CAB
!     GLEAF_TILE, & ! Local, kdcorbin, 10/10
!     FRP_TILE

   REAL, dimension(land_pts,ntiles) ::                           &
   !REAL, dimension(land_pts,ntiles) ::                           &
      RESP_P_FT,     &
      G_LEAF,        &
      NPP_FT,     &
!     NPP_FT_old, &
      GPP_FT       
!     GPP_FT_old

   REAL, dimension(land_pts) ::                                            &
      RESP_P,     & 
      NPP,        & 
      GPP

   REAL, dimension(land_pts) ::                                            &
      SNOW_GRD,   &  
      CANOPY_GB

   REAL, DIMENSION(land_pts,ntiles) ::                                 &
      NPP_FT_ACC,    & ! sresp for CASA-CNP
      RESP_W_FT_ACC    ! presp for CASA-CNP
   
   REAL, DIMENSION(land_pts,ntiles) ::                                 &
      SNOW_TILE,     & !
      SNOW_RHO1L,    & ! Mean snow density
      SNOW_AGE,    & !
      CANOPY_TILE,   & !
      T1P5M_TILE,    &
      Q1P5M_TILE,    &
      TSTAR_TILE,    &
      RESP_S_TILE,   & 
!     RESP_P_FT_old, &
      TRANSP_TILE

   REAL ::                                                                     &
      RESP_S(LAND_PTS,DIM_CS1),    & !
      RESP_S_old(LAND_PTS,DIM_CS1),& !
      RESP_S_TOT(DIM_CS2)                !
  
   REAL, DIMENSION(mp) ::                                                                     &
      fe_dlh,    & !
      fes_dlh,   & !
      fev_dlh      !

   !--- Local vars
   INTEGER :: i,j,l,k,n,m

   REAL, DIMENSION(land_pts,ntiles) ::                                 &
         !--- Local buffer surface FTL, FQL @ prev dt
         FTL_TILE_old, FQW_TILE_old, &
         lpts_ntiles

   INTEGER:: i_miss = 0
   REAL :: miss = 0.0
   
   REAL, POINTER :: TFRZ

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
    fprintf_dir=trim(fprintf_dir_root)//trim("impl_unpack")//"/"
  
  !if( L_fprint_HW ) then
  !  if ( ktau_gl==1 .OR. ktau_gl==54 .OR. ktau_gl==154 .OR. &
  !       ktau_gl==154 .OR. ktau_gl==154 .OR. mod(ktau_gl,10)==0. ) then
  !    L_fprint = .true.
  !  endif  
  !endif  

  !vname='' 
  !call cable_fprintf( cDiagX, vname, %, mp, L_fprint )
  !fprintf_____________________________________________________________________} 

!  LOGICAL, SAVE :: first_call = .TRUE.
   
  IF(cable_user%run_diag_level == "BASIC")                                    &     
    CALL basic_diag(subr_name, "Called.") 

  TFRZ => PHYS%TFRZ
  
  !--- set UM vars to zero
  SMCL_TILE = 0.; STHF_TILE = 0.; STHU_TILE = 0.
  TSOIL_TILE = 0.

  if( L_fprint_HW ) L_fprint = .true.

  vname='ssnow_wb_1' 
  call cable_fprintf( cDiag0, vname, ssnow%wb(:,1), mp, L_fprint )
  
  vname='ssnow_wb_6' 
  call cable_fprintf( cDiag1, vname, ssnow%wb(:,6), mp, L_fprint )
  
  vname='ssnow_tgg_1' 
  call cable_fprintf( cDiag2, vname, ssnow%tgg(:,1), mp, L_fprint )
  
  vname='ssnow_tgg_6' 
  call cable_fprintf( cDiag3, vname, ssnow%tgg(:,6), mp, L_fprint )

  DO j = 1,SM_LEVELS
     TSOIL_TILE(:,:,j)= UNPACK(ssnow%tgg(:,j), L_TILE_PTS, miss)
     SMCL_TILE(:,:,j)= UNPACK(REAL(ssnow%wb(:,j)), L_TILE_PTS, miss)
     SMCL_TILE(:,:,j)=SMCL_TILE(:,:,j)*soil%zse(j)*RHO_WATER
     STHF_TILE(:,:,j)= UNPACK(REAL(ssnow%wbice(:,j)), L_TILE_PTS, miss)
     SMCL(:,j) = SUM(um1%TILE_FRAC * SMCL_TILE(:,:,j),2)
     TSOIL(:,j) = SUM(um1%TILE_FRAC * TSOIL_TILE(:,:,j),2)
     
     DO N=1,NTILES
        DO K=1,um1%TILE_PTS(N)
           I = um1%TILE_INDEX(K,N)
           IF ( SMVCST(I) > 0. ) THEN ! Exclude permanent ice - mrd
              STHF_TILE(I,N,J)= STHF_TILE(I,N,J)/SMVCST(I)
              STHU_TILE(I,N,J)= MAX( 0., SMCL_TILE(I,N,J) -                &
                                STHF_TILE(I,N,J) * SMVCST(I) * soil%zse(J) &
                                * RHO_WATER ) / ( soil%zse(J) *        &
                                RHO_WATER * SMVCST(I) )
           ENDIF
        ENDDO
     ENDDO

     STHF(:,J) = SUM(um1%TILE_FRAC * STHF_TILE(:,:,J),2)
     STHU(:,J) = SUM(um1%TILE_FRAC * STHU_TILE(:,:,J),2)
  ENDDO





      !--- unpack snow vars 
      SNOW_RHO1L  = UNPACK(ssnow%ssdnn, L_TILE_PTS, miss)
      ISNOW_FLG3L = UNPACK(ssnow%isflag, L_TILE_PTS, i_miss)
      MELT_TILE   = UNPACK(ssnow%smelt, L_TILE_PTS, miss)
      SNOW_TILE= UNPACK(ssnow%snowd, L_TILE_PTS, miss)
      SNOW_GRD=  SUM(um1%TILE_FRAC * SNOW_TILE,2)  ! gridbox snow mass & snow below canopy 
      !--- unpack layered snow vars 
      do k = 1,3
        SNOW_TMP3L(:,:,k) = UNPACK(ssnow%tggsn(:,k), L_TILE_PTS, miss)
        SNOW_MASS3L(:,:,k)= UNPACK(ssnow%smass(:,k), L_TILE_PTS, miss)
        SNOW_RHO3L(:,:,k) = UNPACK(ssnow%ssdn(:,k), L_TILE_PTS, miss)
        SNOW_COND(:,:,k)  = UNPACK(ssnow%sconds(:,k),L_TILE_PTS,miss)
        SNOW_DEPTH3L(:,:,k)  = UNPACK(ssnow%sdepth(:,k),L_TILE_PTS,miss)
      enddo

      
      canopy%gswx_T = canopy%gswx_T/air%cmolar
      GS_TILE = UNPACK(canopy%gswx_T,L_TILE_PTS,miss)
      GS =  SUM(um1%TILE_FRAC * GS_TILE,2)

      !---preserve fluxes from the previous time step for the coastal grids
      FTL_TILE_old = FTL_TILE
      FQW_TILE_old = FQW_TILE
      !___return fluxes
      FTL_TILE = UNPACK(canopy%fh,  l_tile_pts, miss)
!      fe_dlh = canopy%fe/(air%rlam*ssnow%cls)
      fes_dlh = canopy%fes/(air%rlam*ssnow%cls)
      fev_dlh = canopy%fev/air%rlam
      fe_dlh =  fev_dlh + fes_dlh

      !---update fluxes 
      FQW_TILE = UNPACK(fe_dlh, l_tile_pts, miss)


      !___return temp and roughness
      lpts_ntiles= tstar_tile
  vname='rad_trad' 
  call cable_fprintf( cDiag4, vname, rad%trad, mp, L_fprint )
!testing
!TSTAR_TILE = UNPACK(rad%trad, l_tile_pts, lpts_ntiles )
      TSTAR_TILE = UNPACK(rad%trad, l_tile_pts, miss)

      !___return miscelaneous 
      RADNET_TILE = unpack( canopy%rnet , l_tile_pts, miss)

     SURF_HTF_TILE = UNPACK(canopy%ga,L_TILE_PTS,miss)

     TOT_ALB=UNPACK(rad%albedo_T,L_TILE_PTS, miss) 

     ESOIL_TILE = UNPACK(fes_dlh, L_TILE_PTS, miss)

     ECAN_TILE = UNPACK(fev_dlh,  L_TILE_PTS, miss)

     EI_TILE = 0.

     SNOW_AGE = UNPACK(ssnow%snage, L_TILE_PTS, miss) 

     TRANSP_TILE = UNPACK(canopy%fevc, L_TILE_PTS, miss) 

     !unpack screen level (1.5m) variables
     !Convert back to K 
     t1p5m_tile     = UNPACK(canopy%tscrn+tfrz, L_TILE_PTS, miss)

     q1p5m_tile     = UNPACK(canopy%qscrn, L_TILE_PTS, miss)

     CANOPY_TILE    = UNPACK(canopy%cansto, L_TILE_PTS, miss)

     CANOPY_GB      = SUM(um1%TILE_FRAC * CANOPY_TILE,2)

     !initialse full land grids and retain coastal grid fluxes
      DO N=1,NTILES
         DO K=1,um1%TILE_PTS(N)
           L = um1%TILE_INDEX(K,N)
           J=(um1%LAND_INDEX(L)-1)/ROW_LENGTH + 1
           I = um1%LAND_INDEX(L) - (J-1)*ROW_LENGTH
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

      DO N=1,NTILES
         DO K=1,um1%TILE_PTS(N)
           L = um1%TILE_INDEX(K,N)
           J=(um1%LAND_INDEX(L)-1)/ROW_LENGTH + 1
           I = um1%LAND_INDEX(L) - (J-1)*ROW_LENGTH
           FTL_1(I,J) = FTL_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FTL_TILE(L,N)
           FQW_1(I,J) = FQW_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FQW_TILE(L,N)
           SURF_HT_FLUX_LAND(I,J) = SURF_HT_FLUX_LAND(I,J) +                   &
                                    FLAND(L)*um1%TILE_FRAC(L,N) *              &
                                    SURF_HTF_TILE(L,N)
         ENDDO
      ENDDO

!    rml 30/10/15 Initialise grid-cell carbon fields that are accumulated over tiles
     RESP_P = 0.
     NPP = 0.
     GPP = 0.
     RESP_S = 0.

     ! Lestevens - Passing CO2 from CABLE to bl_trmix_dd.F90
!    rml 30/10/15 unpack cable carbon variables straight into um variables
!    (previously to local)
     RESP_S_TILE    = UNPACK(canopy%frs, L_TILE_PTS, miss)

!    NEE_TILE       = UNPACK(canopy%fnee, L_TILE_PTS, miss)
!May1! reverted to UM declared as npft - therefore unpacking wont work
!May1!     NPP_FT         = UNPACK(canopy%fnpp, L_TILE_PTS, miss)
!May1!
 !May1!    G_LEAF = UNPACK(canopy%frday,L_TILE_PTS, miss)
!May1!
 !May1!     RESP_P_FT = UNPACK(canopy%frp, L_TILE_PTS, miss)

!    rml - probably want to retire this switch and always have GPP as real GPP
!    (off case)
!May1!     IF( cable_user%leaf_respiration == 'on' .OR.                             &
!May1!           cable_user%leaf_respiration == 'ON') THEN
!May1!        GPP_FT = UNPACK(canopy%fnpp+canopy%frp, L_TILE_PTS, miss)
!May1!     ELSE
!May1!        GPP_FT = UNPACK(canopy%fnpp+canopy%frp+canopy%frday,  &
!May1!                            L_TILE_PTS, miss)
!May1!     ENDIF

!    convert from CABLE units (gC/m2/s) to UM units (kgC/m2/s)
     RESP_S_TILE = 0.0
     G_LEAF      = 0.0  
     NPP_FT      = 0.0 
     GPP_FT      = 0.0 
     RESP_P_FT   = 0.0 
     RESP_S_TILE = RESP_S_TILE*1.e-3
     G_LEAF      = G_LEAF*1.e-3
     NPP_FT      = NPP_FT*1.e-3
     GPP_FT      = GPP_FT*1.e-3
     RESP_P_FT   = RESP_P_FT*1.e-3

     ! Lestevens 23apr13 - possible miss match ntiles<-->npft
!    If CASA-CNP used, plant and soil resp need to be passed into 
!    variables that are dumped to restart, because CASA-CNP only run daily
     NPP_FT_ACC = RESP_S_TILE
     RESP_W_FT_ACC = RESP_P_FT

!May1!      DO N=1,NTILES 
!May1!         DO K=1,um1%TILE_PTS(N)
!May1!            L = um1%TILE_INDEX(K,N)
!May1!            NPP(L)=NPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*NPP_FT(L,N)
!May1!            GPP(L)=GPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*GPP_FT(L,N)
!May1!            RESP_P(L)=RESP_P(L)+FLAND(L)*um1%TILE_FRAC(L,N)*RESP_P_FT(L,N)
!May1!
!May1!            !loop for soil resp. although DIM_CS1=1 (not 1 for triffid)
!May1!            DO I=1,DIM_CS1
!May1!               RESP_S(L,I) = RESP_S(L,I) + &
!May1!                             FLAND(L)*um1%TILE_FRAC(L,N)*RESP_S_TILE(L,N)
!May1!            ENDDO
!May1!            RESP_S_TOT(L)=sum(RESP_S(L,:))
!May1!         ENDDO
!May1!      ENDDO


!*!call implicit_nunpack( & 
!*!                            cycleno,                                                         & 
!*!                            row_length,rows, land_pts, ntiles, npft, sm_levels,              &
!*!                            dim_cs1, dim_cs2,                                                & 
!*!                            TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,                &
!*!                            SMVCST, STHF, STHF_TILE, STHU, STHU_TILE,          &
!*!                            snow_tile, SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L,  &
!*!                            SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND,    &
!*!                            FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
!*!                            SURF_HT_FLUX_LAND, ECAN_TILE,        &
!*!                            ESOIL_TILE, EI_TILE, RADNET_TILE, TOT_ALB,         &
!*!                            SNOW_AGE, CANOPY_TILE, GS, gs_tile, T1P5M_TILE,  &
!*!                            !args increased following merge
!*!                            Q1P5M_TILE, CANOPY_GB, FLAND, MELT_TILE, &
!*!                            NPP, NPP_FT, GPP, GPP_FT, RESP_S,         &
!*!                            RESP_S_TOT, RESP_S_TILE, RESP_P, RESP_P_FT, G_LEAF,&
!*!                            TRANSP_TILE, NPP_FT_ACC, RESP_W_FT_ACC,SURF_HTF_TILE )
!*! 

END SUBROUTINE Implicit_unpack


End module cable_implicit_unpack_mod
