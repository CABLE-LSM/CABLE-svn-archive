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
! Called from: cable_explicit_driver
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Developed for CABLE v1.8
!
!
! ==============================================================================

!---------------------------------------------------------------------!
!--- pass land-surface quantities calc'd by CABLE in explicit call ---!
!--- back to UM.                                                   ---!
!---------------------------------------------------------------------!
SUBROUTINE cable_expl_unpack( DIM_CS1, DIM_CS2,                               &
                              FTL_TILE, FQW_TILE,                             &
                              TSTAR_TILE,                                     &
                              U_S, U_S_STD_TILE,                              &
                              CD_TILE, CH_TILE, FLAND, RADNET_TILE,           &
                              FRACA, rESFS, RESFT, Z0H_TILE, Z0M_TILE,        &
                              RECIP_L_MO_TILE, EPOT_TILE,                     &
                              GS, NPP, NPP_FT, GPP, GPP_FT, RESP_S,           &
                              RESP_S_TOT, RESP_P, RESP_P_FT, G_LEAF,          &
                              l_tile_pts )

  USE cable_def_types_mod, ONLY : mp, NITER
  USE cable_data_module,   ONLY : PHYS
  USE cable_um_tech_mod,   ONLY : um1, ssnow, air, canopy, met, rad, rough
  USE cable_common_module, ONLY : cable_runtime, cable_user, &
                                  ktau_gl, knode_gl, kend_gl
  IMPLICIT NONE


  !--------------------------------------------------------------------------
  !--- INPUT ARGS FROM cable_explicit_driver() ------------------------------
  !--------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: DIM_CS1, DIM_CS2

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
     RADNET_TILE,     & ! Surface net radiation
     RESFS,           & ! Combined soil, stomatal & aerodynamic resistance
                        ! factor for fraction (1-FRACA) of snow-free land tiles
     RESFT,           & ! Total resistance factor.
                        ! FRACA+(1-FRACA)*RESFS for snow-free l_tile_pts,
                        ! 1 for snow.
     FRACA,           & ! Fraction of surface moisture
     RECIP_L_MO_TILE, &! Reciprocal of the Monin-Obukhov length for tiles (m^-1).
     EPOT_TILE,       &
     NPP_FT,          &
     GPP_FT,          &
     RESP_P_FT,       &
     G_LEAF,          &
     TSTAR_TILE

  REAL, INTENT(OUT), DIMENSION(um1%land_pts) ::                              &
     GS,         &
     NPP,        &
     GPP,        &
     RESP_P

  REAL, INTENT(OUT) ::                                                       &
    RESP_S(um1%LAND_PTS,DIM_CS1),    & !
    RESP_S_TOT(DIM_CS2)                !

  LOGICAL,DIMENSION(um1%land_pts,um1%ntiles) :: l_tile_pts

  !___UM vars used but NOT returned
  REAL, INTENT(IN), DIMENSION(um1%land_pts) ::   &
     FLAND(um1%land_pts)              ! IN Land fraction on land tiles.


!-----------------------------------------------------------------------------
! Local variables
  REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
     LE_TILE, GS_TILE, NPP_TILE, NPP_FT_old, GPP_TILE, GPP_FT_old, FRS_TILE,  &
     RESP_P_FT_old, FRP_TILE, GLEAF_TILE


  REAL ::                                                                     &
    RESP_S_old(um1%LAND_PTS,DIM_CS1)

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

!-----------------------------------------------------------------------------

  CAPP => PHYS%CAPP

  !___return fluxes
  fe_dlh = canopy%fe/(air%rlam*ssnow%cls)
  FTL_TILE = UNPACK(canopy%fh,  um1%l_tile_pts, miss)
  FTL_TILE = FTL_TILE / capp
  FQW_TILE = UNPACK(fe_dlh, um1%l_tile_pts, miss)

  !___return temp and roughness
  TSTAR_TILE = UNPACK(rad%trad,  um1%l_tile_pts, miss)
  Z0M_TILE = UNPACK(rough%z0m,  um1%l_tile_pts, miss)
  Z0H_TILE = Z0M_TILE

  !___return friction velocities/drags/ etc
  U_S_TILE  =  UNPACK(canopy%us, um1%l_tile_pts, miss)
  CDCAB = canopy%us**2/met%ua**2   ! met%ua is always above umin = 0.1m/s
  ! for Cable CD*
  CD_TILE =  UNPACK(CDCAB,um1%l_tile_pts, miss)
  ! for Cable CH*
  CH_TILE =  UNPACK(canopy%cdtq,um1%l_tile_pts, miss)

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
  fraca_cab = canopy%fwet * (1.-rad%transd)
  WHERE( ssnow%snowd > 1.0 ) fraca_cab = 1.0
  rfsfs_cab = MIN( 1., MAX( 0.01, canopy%wetfac_cs - fraca_cab ) /         &
              MAX( 0.01,1. - fraca_cab ) )
  FRACA = UNPACK( fraca_cab, um1%l_tile_pts, miss )
  RESFT = UNPACK( canopy%wetfac_cs,um1%l_tile_pts, miss )
  RESFS = UNPACK( rfsfs_cab , um1%l_tile_pts, miss )

  RADNET_TILE = UNPACK( canopy%rnet , um1%l_tile_pts, miss )
  THETAST = ABS( canopy%fh ) / ( air%rho * capp*canopy%us )
  RECIPLMOTILE =  canopy%zetar(:,niter) / rough%zref_tq
  RECIP_L_MO_TILE = UNPACK( RECIPLMOTILE, um1%l_tile_pts, miss )
  EPOT_TILE = UNPACK( canopy%epot, um1%l_tile_pts, miss )

    !---???
  GS_TILE = UNPACK(canopy%gswx_T,um1%L_TILE_PTS,miss)
  GS = SUM(um1%TILE_FRAC * GS_TILE,2)

  NPP_TILE   = UNPACK(canopy%fnpp, um1%L_TILE_PTS, miss)
  FRS_TILE   = UNPACK(canopy%frs, um1%L_TILE_PTS, miss)
  FRP_TILE   = UNPACK(canopy%frp, um1%L_TILE_PTS, miss)
  GLEAF_TILE = UNPACK(canopy%frday,um1%L_TILE_PTS, miss)

  IF( cable_user%leaf_respiration == 'on' .OR.                                &
      cable_user%leaf_respiration == 'ON') THEN
    GPP_TILE = UNPACK(canopy%fnpp+canopy%frp, um1%L_TILE_PTS, miss)
  ELSE
    GPP_TILE = UNPACK(canopy%fnpp+canopy%frp+canopy%frday,                    &
                                                         um1%L_TILE_PTS, miss)
  END IF

  NPP_FT_old    = NPP_FT
  GPP_FT_old    = GPP_FT
  RESP_S_old    = RESP_S
  RESP_P_FT_old = RESP_P_FT

  DO N=1,um1%NTILES
    DO K=1,um1%TILE_PTS(N)
      L = um1%TILE_INDEX(K,N)
      IF( FLAND(L) == 1.0) THEN
        NPP(L)         = 0.0
        NPP_FT(L,N)    = 0.0
        GPP(L)         = 0.0
        GPP_FT(L,N)    = 0.0
        RESP_P(L)      = 0.0
        RESP_P_FT(L,N) = 0.0
        RESP_S(L,:)    = 0.0
        G_LEAF(L,N)    = 0.0
      ELSE
        ! For coastal points: currently no contribution
        NPP(L) = NPP(L)-FLAND(L)*um1%TILE_FRAC(L,N)*NPP_FT_old(L,N)
        GPP(L) = GPP(L)-FLAND(L)*um1%TILE_FRAC(L,N)*GPP_FT_old(L,N)
        RESP_P(L) = RESP_P(L)-FLAND(L)*um1%TILE_FRAC(L,N)*RESP_P_FT_old(L,N)
        !--- loop for soil respiration
        DO I=1,DIM_CS1
          RESP_S(L,I) = RESP_S(L,I)-FLAND(L)*RESP_S_old(L,I)
        END DO
        RESP_S_TOT(L) = SUM(RESP_S(L,:))
      END IF
    END DO
  END DO

  DO N=1,um1%NTILES
    DO K=1,um1%TILE_PTS(N)
      L = um1%TILE_INDEX(K,N)
      !add leaf respiration to output
      G_LEAF(L,N) = GLEAF_TILE(L,N)*1.e-3
      NPP_FT(L,N) = NPP_TILE(L,N)*1.e-3
      NPP(L) = NPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*NPP_FT(L,N)
      GPP_FT(L,N) = GPP_TILE(L,N)*1.e-3
      GPP(L) = GPP(L)+FLAND(L)*um1%TILE_FRAC(L,N)*GPP_FT(L,N)

      !loop for soil resp. - all UM levels = single CABLE output
      DO I=1,DIM_CS1
        RESP_S(L,I) = RESP_S(L,I) +                                           &
                      FLAND(L)*um1%TILE_FRAC(L,N)*FRS_TILE(L,N)*1.e-3
      END DO

      RESP_S_TOT(L) = SUM(RESP_S(L,:))
      RESP_P_FT(L,N) = FRP_TILE(L,N)*1.e-3
      RESP_P(L) = RESP_P(L)+FLAND(L)*um1%TILE_FRAC(L,N)*RESP_P_FT(L,N)
    END DO
  END DO

  IF(first_cable_call) THEN
     l_tile_pts = um1%l_tile_pts
     first_cable_call = .FALSE.
  ENDIF

   
END SUBROUTINE cable_expl_unpack
    
!============================================================================
!============================================================================
!============================================================================

