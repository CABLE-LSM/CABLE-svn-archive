#include "../include/cable_directives.h"

   subroutine cable_implicit_driver(                      &
                  LS_RAIN, CON_RAIN, LS_SNOW, CONV_SNOW, DTL_1,DQW_1,&
                  TSOIL, TSOIL_TILE, SMCL, SMCL_TILE, timestep,                             &
                  SMVCST,STHF, STHF_TILE, STHU,STHU_TILE, snow_tile,                         &
                  SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L, SNOW_MASS3L, &
                  SNOW_RHO3L,SNOW_TMP3L,SNOW_COND, FTL_TILE_CAB,     &
                  FTL_CAB,LE_TILE_CAB,LE_CAB, FTL_1,FTL_TILE,FQW_1,  &
                  FQW_TILE, TSTAR_TILE, TSTAR_TILE_CAB, TSTAR_CAB, &
                  SMCL_CAB,TSOIL_CAB, SURF_HTF_CAB, SURF_HT_FLUX_LAND, &
                  ECAN_TILE,ESOIL_TILE,EI_TILE,RADNET_TILE,TOT_ALB,  &
                  SNAGE_TILE,CANOPY_TILE, GS, T1P5M_TILE, Q1P5M_TILE, &
                  CANOPY_GB,FLAND,MELT_TILE, DIM_CS1,DIM_CS2, NPP, NPP_FT, &
                  GPP, GPP_FT, RESP_S, RESP_S_TOT, RESP_S_TILE,     & 
                  RESP_P,RESP_P_FT, G_LEAF)  
 
      use define_dimensions, only : mp
      use physical_constants, only : tfrz
      use cable_variables
      use cable_um_tech_mod, only : um1
      use cable_common_module, only : cable_runtime, cable_user
      use cable_um_init_subrs, only : um2cable_rr
      use cbm_module, only : cbm
      use cable_diag_module, only : cable_stat

      implicit none
           
      real, dimension(um1%ROW_LENGTH,um1%ROWS) :: &
       LS_RAIN  & ! IN Large scale rain
      ,LS_SNOW  & ! IN Large scale snow
      ,CON_RAIN & ! IN Convective rain
      ,CONV_SNOW & ! IN Convective snow
      ,DTL_1,    &     ! IN Level 1 increment to T field 
       DQW_1  ! IN Level 1 increment to q field 

      real :: timestep

      integer :: DIM_CS1 ,DIM_CS2 

      real, dimension(um1%land_pts) :: &
            GS,         &  ! OUT "Stomatal" conductance to
            SMVCST,     &  ! IN Volumetric saturation point
            FLAND          ! IN Land fraction on land tiles
      
      real, dimension(um1%ROW_LENGTH,um1%ROWS) :: &
            !--- Net downward heat flux at surface over land.
            !--- fraction of gridbox (W/m2).
            SURF_HT_FLUX_LAND,           &
            !--- Moisture flux between layers. (kg/m^2/sec).
            !--- FQW(,1) is total water flux from surface, 'E'.
            FQW_1,       &  
            !--- FTL(,K) =net turbulent sensible heat flux into layer K
            !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
            FTL_1         

      REAL, dimension(um1%land_pts,um1%ntiles) :: &
            !___Surface FTL, FQL for land tiles
            FTL_TILE, FQW_TILE, FQW_TILE_CAB,   &  
            !___(tiled) latent heat flux, melting, stomatatal conductance
            LE_TILE, MELT_TILE, GS_TILE,     &  
            !___ INOUT Surface net radiation on tiles (W/m2)
            RADNET_TILE,   &
            !___total albedo
            TOT_ALB,    &
            EI_TILE,    & ! OUT EI for land tiles.
            ECAN_TILE   & ! OUT ECAN for snow-free land tiles
            !___evapotranspiration from soil moisture store (kg/m2/s) 
            ,ESOIL_TILE    

      !___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
      !___ runoff ??
      REAL, dimension(um1%land_pts,um1%sm_levels) :: &
            SMCL, STHF, STHU, SMCL_CAB, TSOIL_CAB, TSOIL, &
            SURF_CAB_ROFF       

      !___(tiled) soil prognostics: as above 
      REAL, dimension(um1%land_pts,um1%ntiles,um1%sm_levels) :: &
            SMCL_TILE, STHU_TILE, TSOIL_TILE, STHF_TILE  

      !___flag for 3 layer snow pack
      INTEGER :: ISNOW_FLG3L(um1%LAND_PTS,um1%NTILES)
      
      !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
      REAL, dimension(um1%land_pts,um1%ntiles,3) :: &
            SNOW_DEPTH3L ,SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND 

      REAL, dimension(um1%land_pts,um1%ntiles) :: &
      FRS_TILE         & ! Local
      ,NEE_TILE         & ! Local
      ,NPP_TILE         & ! Local
      ,GPP_TILE         & ! Local
      ,GLEAF_TILE      & ! Local, kdcorbin, 10/10
      ,FRP_TILE       &
      ,NPP_FT           &
      ,NPP_FT_old       &
      ,GPP_FT           &
      ,GPP_FT_old       

      REAL, dimension(um1%land_pts) :: &
     SNOW_GRD,                &
      CANOPY_GB        & 
      ,FTL_CAB                                             &
      ,LE_CAB                                              &
      ,TSTAR_CAB                                           &
      ,SURF_HTF_CAB                                       &
       ,RESP_P                  &
      ,NPP                   &
      ,GPP
      
      REAL, dimension(um1%land_pts,um1%ntiles) :: &
     SNOW_TILE,       &
     SNOW_RHO1L       &  ! Mean snow density
      ,SNAGE_TILE       &
      ,CANOPY_TILE      &
       ,FTL_TILE_CAB                                 &
      ,LE_TILE_CAB                                  &
      ,T1P5M_TILE                                  &
      ,Q1P5M_TILE                                &
      ,TSTAR_TILE_CAB                               &
      ,TSTAR_TILE                                   &
      ,SURF_HTF_T_CAB                            & 
      ,RESP_S_TILE      & 
      ,RESP_P_FT          &
      ,RESP_P_FT_old  &
      ,G_LEAF

      real :: & 
      RESP_S(um1%LAND_PTS,DIM_CS1)          &
      ,RESP_S_old(um1%LAND_PTS,DIM_CS1)      &
      ,RESP_S_TOT(DIM_CS2)    
     
      real, dimension(mp) ::  dtlc, dqwc

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('cable_implicit_driver')
         
         !___FLAGS def. specific call of CABLE from UM
         cable_runtime%um_explicit = .false.
         cable_runtime%um_hydrology = .true.
      
         dtlc = 0. ; dqwc = 0.

         !jhan: in principle this shouln't need to be done and so far in testing doesnt
         !however, SHOULd only need to be done once in explicit, but was forgetting pointers 
         !if not re-assigned. check & also in other CALLs to CABLE pass rows, etc
         !call assign_um_fundamentals_to_um1( row_length, rows, land_pts, ntiles,  &
         !                              npft, sm_levels, timestep, latitude, longitude, &
         !                              land_index, tile_frac, tile_pts, tile_index,    &
         !                              l_tile_pts, rho_water  )
         
         !--- All these subrs do is pack a CABLE var with a UM var.
         !-------------------------------------------------------------------
         !--- UM met forcing vars needed by CABLE which have UM dimensions
         !---(land_pts,ntiles)[_lp], which is no good to CABLE. These have to be 
         !--- re-packed in a single vector of active tiles. Hence we use 
         !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
         !--- if the land point is/has an active tile
         !--- generic format:
         !--- um2cable_lp( UM var, default value for snow tile, CABLE var, mask )
         !--- where mask tells um2cable_lp whether or not to use default value 
         !--- for snow tile 
         !-------------------------------------------------------------------
         call um2cable_rr( (LS_RAIN+CON_RAIN)*um1%TIMESTEP, met%precip)
         call um2cable_rr( (LS_SNOW+CONV_SNOW)*um1%TIMESTEP, met%precip_sn)
         call um2cable_rr( dtl_1, dtlc)
         call um2cable_rr( dqw_1, dqwc)
   
         met%precip   =  met%precip + met%precip_sn
         met%tk = met%tk + dtlc
         met%tc = met%tk - tfrz
         met%qv = met%qv + dqwc
         met%tvair = met%tk
         met%tvrad = met%tk
    
         canopy%cansto = canopy%oldcansto

         CALL cbm(TIMESTEP, air, bgc, canopy, met, bal,  &
              rad, rough, soil, ssoil, sum_flux, veg)
  
        
         call implicit_unpack( TSOIL, TSOIL_TILE, SMCL, SMCL_TILE,                             &
                  SMVCST,STHF, STHF_TILE, STHU,STHU_TILE, snow_tile,                         &
                  SNOW_RHO1L ,ISNOW_FLG3L, SNOW_DEPTH3L, SNOW_MASS3L, &
                  SNOW_RHO3L,SNOW_TMP3L,SNOW_COND, FTL_TILE_CAB,     &
                  FTL_CAB,LE_TILE_CAB,LE_CAB, FTL_1,FTL_TILE,FQW_1,  &
                  FQW_TILE, TSTAR_TILE, TSTAR_TILE_CAB, TSTAR_CAB, &
                  SMCL_CAB,TSOIL_CAB, SURF_HTF_CAB, SURF_HT_FLUX_LAND, &
                  ECAN_TILE,ESOIL_TILE,EI_TILE,RADNET_TILE,TOT_ALB,  &
                  SNAGE_TILE,CANOPY_TILE, GS, T1P5M_TILE, Q1P5M_TILE, &
                  CANOPY_GB,FLAND,MELT_TILE, DIM_CS1,DIM_CS2, NPP, NPP_FT, &
                  GPP, GPP_FT, RESP_S, RESP_S_TOT, RESP_S_TILE,     & 
                  RESP_P,RESP_P_FT, G_LEAF)  
  
      return
   end subroutine cable_implicit_driver


!========================================================================= 
!========================================================================= 
!========================================================================= 

   subroutine implicit_unpack(       TSOIL,TSOIL_TILE,SMCL,SMCL_TILE,                             &
       SMVCST,STHF,STHF_TILE,STHU,STHU_TILE,                        &  
       SNOW_TILE,SNOW_RHO1L,ISNOW_FLG3L,                            &
       SNOW_DEPTH3L,SNOW_MASS3L,SNOW_RHO3L,SNOW_TMP3L,SNOW_COND,    &
       FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,                     &
       FTL_1,FTL_TILE,FQW_1,FQW_TILE,                               &
       TSTAR_TILE,                                                  &
       TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,                 &
       SURF_HTF_CAB,                                      &
       SURF_HT_FLUX_LAND,ECAN_TILE,ESOIL_TILE,EI_TILE,RADNET_TILE   &
       ,TOT_ALB &
       ,SNAGE_TILE,CANOPY_TILE,   &
       GS,T1P5M_TILE,Q1P5M_TILE,CANOPY_GB,FLAND,MELT_TILE,          & 
       DIM_CS1,DIM_CS2,                                             &
       NPP,NPP_FT,GPP,GPP_FT,RESP_S,RESP_S_TOT,RESP_S_TILE, & 
       RESP_P,RESP_P_FT,    &
       G_LEAF)   

      use define_dimensions, only : mp
      use physical_constants, only : tfrz
      use cable_variables
      use cable_um_tech_mod, only : um1
      use cable_common_module, only : cable_runtime, cable_user
      use cable_diag_module, only : cable_stat
      implicit none
    
!       real, dimension(um1%ROW_LENGTH,um1%ROWS) :: &
!       LS_RAIN  & ! IN Large scale !rain
!      ,LS_SNOW  & ! IN Large scale snow
!      ,CON_RAIN & ! IN Convective rain
!      ,CONV_SNOW & ! IN Convective snow
!      ,DTL_1,    &     ! IN Level 1 increment to T field 
!       DQW_1  ! IN Level 1 increment to q field 
!
!      real :: timestep
!jhan:these need to be cleaned out to what is actualllly passed
      integer :: DIM_CS1 ,DIM_CS2 

      real, dimension(um1%land_pts) :: &
            GS,         &  ! OUT "Stomatal" conductance to
            SMVCST,     &  ! IN Volumetric saturation point
            FLAND          ! IN Land fraction on land tiles
      
      real, dimension(um1%ROW_LENGTH,um1%ROWS) :: &
            !--- Net downward heat flux at surface over land.
            !--- fraction of gridbox (W/m2).
            SURF_HT_FLUX_LAND,           &
            !--- Moisture flux between layers. (kg/m^2/sec).
            !--- FQW(,1) is total water flux from surface, 'E'.
            FQW_1,       &  
            !--- FTL(,K) =net turbulent sensible heat flux into layer K
            !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
            FTL_1         

      REAL, dimension(um1%land_pts,um1%ntiles) :: &
            !___Surface FTL, FQL for land tiles
            FTL_TILE, FQW_TILE, FQW_TILE_CAB,   &  
            !___(tiled) latent heat flux, melting, stomatatal conductance
            LE_TILE, MELT_TILE, GS_TILE,     &  
            !___ INOUT Surface net radiation on tiles (W/m2)
            RADNET_TILE,   &
            !___total albedo
            TOT_ALB,    &
            EI_TILE,    & ! OUT EI for land tiles.
            ECAN_TILE   & ! OUT ECAN for snow-free land tiles
            !___evapotranspiration from soil moisture store (kg/m2/s) 
            ,ESOIL_TILE    

      !___ soil prognostics: moisture, frozen, unfrozen content, soil temp.
      !___ runoff ??
      REAL, dimension(um1%land_pts,um1%sm_levels) :: &
            SMCL, STHF, STHU, SMCL_CAB, TSOIL_CAB, TSOIL, &
            SURF_CAB_ROFF       

      !___(tiled) soil prognostics: as above 
      REAL, dimension(um1%land_pts,um1%ntiles,um1%sm_levels) :: &
            SMCL_TILE, STHU_TILE, TSOIL_TILE, STHF_TILE  

      !___flag for 3 layer snow pack
      INTEGER :: ISNOW_FLG3L(um1%LAND_PTS,um1%NTILES)
      
      !___(tiled, 3 layer) Snow depth (m), mass, density, temp., conductivity
      REAL, dimension(um1%land_pts,um1%ntiles,3) :: &
            SNOW_DEPTH3L ,SNOW_MASS3L, SNOW_RHO3L, SNOW_TMP3L, SNOW_COND 

      REAL, dimension(um1%land_pts,um1%ntiles) :: &
      FRS_TILE         & ! Local
      ,NEE_TILE         & ! Local
      ,NPP_TILE         & ! Local
      ,GPP_TILE         & ! Local
      ,GLEAF_TILE      & ! Local, kdcorbin, 10/10
      ,FRP_TILE       &
      ,NPP_FT           &
      ,NPP_FT_old       &
      ,GPP_FT           &
      ,GPP_FT_old       

      REAL, dimension(um1%land_pts) :: &
     SNOW_GRD,                &
      CANOPY_GB        & 
      ,FTL_CAB                                             &
      ,LE_CAB                                              &
      ,TSTAR_CAB                                           &
      ,SURF_HTF_CAB                                       &
       ,RESP_P                  &
      ,NPP                   &
      ,GPP
      
      REAL, dimension(um1%land_pts,um1%ntiles) :: &
     SNOW_TILE,       &
     SNOW_RHO1L       &  ! Mean snow density
      ,SNAGE_TILE       &
      ,CANOPY_TILE      &
       ,FTL_TILE_CAB                                 &
      ,LE_TILE_CAB                                  &
      ,T1P5M_TILE                                  &
      ,Q1P5M_TILE                                &
      ,TSTAR_TILE_CAB                               &
      ,TSTAR_TILE                                   &
      ,SURF_HTF_T_CAB                            & 
      ,RESP_S_TILE      & 
      ,RESP_P_FT          &
      ,RESP_P_FT_old  &
      ,G_LEAF

      real :: & 
      RESP_S(um1%LAND_PTS,DIM_CS1)          &
      ,RESP_S_old(um1%LAND_PTS,DIM_CS1)      &
      ,RESP_S_TOT(DIM_CS2)    
     
!      real, dimension(mp) ::  dtlc, dqwc
  
       REAL fe_dlh(mp),fes_dlh(mp),fev_dlh(mp)

      !--- Local vars
      INTEGER :: i,j,l,k,n

      REAL, dimension(um1%land_pts,um1%ntiles) :: &
            !--- Local buffer surface FTL, FQL @ prev dt
            FTL_TILE_old, FQW_TILE_old

      integer:: i_miss = 0
      real :: miss = 0.0
      
         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('implicit_unpack')
         
         !--- set UM vars to zero
         TSOIL_CAB = 0.; SMCL_CAB = 0.; TSOIL_TILE = 0.; 
         SMCL_TILE = 0.; STHF_TILE = 0.; STHU_TILE = 0.
   
         do j = 1,um1%SM_LEVELS
            TSOIL_TILE(:,:,j)= unpack(ssoil%tgg(:,j), um1%L_TILE_PTS, miss)
            TSOIL_CAB(:,j) = SUM(um1%TILE_FRAC * TSOIL_TILE(:,:,j),2)
            SMCL_TILE(:,:,j)= unpack(real(ssoil%wb(:,j)), um1%L_TILE_PTS, miss)
            SMCL_TILE(:,:,j)=SMCL_TILE(:,:,j)*soil%zse(j)*um1%RHO_WATER
            SMCL_CAB(:,j) = SUM(um1%TILE_FRAC * SMCL_TILE(:,:,j),2)
            STHF_TILE(:,:,j)= unpack(real(ssoil%wbice(:,j)), um1%L_TILE_PTS, miss)
            SMCL(:,j) = SUM(um1%TILE_FRAC * SMCL_TILE(:,:,j),2)
            TSOIL(:,j) = SUM(um1%TILE_FRAC * TSOIL_TILE(:,:,j),2)
            DO N=1,um1%NTILES
               DO K=1,um1%TILE_PTS(N)
                  I = um1%TILE_INDEX(K,N)
                  if ( SMVCST(I) > 0. ) then ! Exclude permanent ice - mrd
                     STHF_TILE(I,N,J)= STHF_TILE(I,N,J)/SMVCST(I)
                     STHU_TILE(I,N,J)= max( 0. ,SMCL_TILE(I,N,J)- STHF_TILE(I,N,J) * &
                           SMVCST(I)*soil%zse(J)*um1%RHO_WATER )/ &
                           (soil%zse(J)*um1%RHO_WATER*SMVCST(I))
                  ENDIF
               ENDDO
            ENDDO
            STHF(:,J) = SUM(um1%TILE_FRAC * STHF_TILE(:,:,J),2)
            STHU(:,J) = SUM(um1%TILE_FRAC * STHU_TILE(:,:,J),2)
         enddo

         !--- unpack snow vars 
         SNOW_RHO1L  = unpack(ssoil%ssdnn, um1%L_TILE_PTS, miss)
         ISNOW_FLG3L = unpack(ssoil%isflag, um1%L_TILE_PTS, i_miss)
         MELT_TILE   = unpack(ssoil%smelt, um1%L_TILE_PTS, miss)
         SNOW_TILE= unpack(ssoil%snowd, um1%L_TILE_PTS, miss)
         SNOW_GRD=  SUM(um1%TILE_FRAC * SNOW_TILE,2)  ! gridbox snow mass & snow below canopy 

         !--- unpack layered snow vars 
         do k = 1,3
           SNOW_TMP3L(:,:,k) = unpack(ssoil%tggsn(:,k), um1%L_TILE_PTS, miss)
           SNOW_MASS3L(:,:,k)= unpack(ssoil%smass(:,k), um1%L_TILE_PTS, miss)
           SNOW_RHO3L(:,:,k) = unpack(ssoil%ssdn(:,k), um1%L_TILE_PTS, miss)
           SNOW_COND(:,:,k)  = unpack(ssoil%sconds(:,k),um1%L_TILE_PTS,miss)
           SNOW_DEPTH3L(:,:,k)  = unpack(ssoil%sdepth(:,k),um1%L_TILE_PTS,miss)
         enddo
   
         !---???
         GS_TILE = unpack(canopy%gswx_T,um1%L_TILE_PTS,miss)
         GS =  SUM(um1%TILE_FRAC * GS_TILE,2)

         !___return fluxes
         FTL_TILE_CAB = unpack(canopy%fh,  um1%l_tile_pts, miss)
         FTL_CAB = SUM(um1%TILE_FRAC * FTL_TILE_CAB,2)
         FQW_TILE_CAB = unpack(canopy%fe,  um1%l_tile_pts, miss)
         LE_TILE_CAB = FQW_TILE_CAB 
         LE_CAB = SUM(um1%TILE_FRAC * LE_TILE_CAB,2)
         fe_dlh = canopy%fe/(air%rlam*ssoil%cls)
         fes_dlh = canopy%fes/(air%rlam*ssoil%cls)
         fev_dlh = canopy%fev/air%rlam

         !---preserve fluxes from the previous time step for the coastal grids
         FTL_TILE_old = FTL_TILE
         FQW_TILE_old = FQW_TILE
         
         !---update fluxes 
         FTL_TILE = FTL_TILE_CAB 
         FQW_TILE = unpack(fe_dlh, um1%l_tile_pts, miss)


         !___return temp and roughness
         TSTAR_TILE_CAB = unpack(rad%trad, um1%l_tile_pts, miss)
         TSTAR_CAB = SUM(um1%TILE_FRAC * TSTAR_TILE_CAB,2)
         TSTAR_TILE = TSTAR_TILE_CAB 
   
         !___return miscelaneous 
         RADNET_TILE = unpack( canopy%rnet , um1%l_tile_pts, miss)

         SURF_HTF_T_CAB = unpack(canopy%ga,um1%L_TILE_PTS,miss)
         SURF_HTF_CAB = SUM(um1%TILE_FRAC * SURF_HTF_T_CAB,2)

        TOT_ALB=unpack(rad%albedo_T,um1%L_TILE_PTS, miss) 
        ESOIL_TILE = unpack(fes_dlh, um1%L_TILE_PTS, miss)
        ECAN_TILE = unpack(fev_dlh,  um1%L_TILE_PTS, miss)
        EI_TILE = 0.
        SNAGE_TILE = unpack(ssoil%snage, um1%L_TILE_PTS, miss) 

        !unpack screen level (1.5m) variables
        !Convert back to K 
        t1p5m_tile     = unpack(canopy%tscrn+tfrz, um1%L_TILE_PTS, miss)
        q1p5m_tile     = unpack(canopy%qscrn, um1%L_TILE_PTS, miss)
        CANOPY_TILE    = unpack(canopy%cansto, um1%L_TILE_PTS, miss)
        CANOPY_GB      = SUM(um1%TILE_FRAC * CANOPY_TILE,2)

        ! Lestevens - Passing CO2 from CABLE to bl_trmix_dd.F90
        FRS_TILE       = unpack(canopy%frs, um1%L_TILE_PTS, miss)
        NEE_TILE       = unpack(canopy%fnee, um1%L_TILE_PTS, miss)
        NPP_TILE       = unpack(canopy%fnpp, um1%L_TILE_PTS, miss)
        GLEAF_TILE     = unpack(canopy%frday,um1%L_TILE_PTS, miss)

         if (cable_user%leaf_respiration == 'on') then
            GPP_TILE       = unpack(canopy%fnpp+canopy%frp, um1%L_TILE_PTS, miss)
         else 
            GPP_TILE       = unpack(canopy%fnpp+canopy%frp+canopy%frday,  &
                               um1%L_TILE_PTS, miss)
         endif

        FRP_TILE       = unpack(canopy%frp, um1%L_TILE_PTS, miss)
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
            FTL_1(I,J) = FTL_1(I,J) - FLAND(L)*um1%TILE_FRAC(L,N)*FTL_TILE_old(L,N)
            FQW_1(I,J) = FQW_1(I,J) - FLAND(L)*um1%TILE_FRAC(L,N)*FQW_TILE_old(L,N)
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
          SURF_HT_FLUX_LAND(I,J) = SURF_HT_FLUX_LAND(I,J) +                        &
                                   FLAND(L)*um1%TILE_FRAC(L,N)*SURF_HTF_T_CAB(L,N)
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

    return

end subroutine implicit_unpack


