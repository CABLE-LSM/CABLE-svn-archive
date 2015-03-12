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

subroutine cable_implicit_driver( LS_RAIN, CON_RAIN, LS_SNOW, CONV_SNOW,       &
                                  DTL_1, DQW_1, timestep,                      &
                                  FTL_1, FTL_TILE, FQW_1, FQW_TILE,            &
                                  TSTAR_TILE,                                  &
                                  SURF_HT_FLUX_LAND, ECAN_TILE, ESOIL_TILE,    &
                                  EI_TILE, RADNET_TILE,                        &
                                  T1P5M_TILE, Q1P5M_TILE,                      &
                                  FLAND, MELT_TILE )

   USE cable_def_types_mod, ONLY : mp
   USE cable_data_module,   ONLY : PHYS
   USE cable_data_mod,   ONLY : cable 
   USE cable_um_tech_mod,   ONLY : um1, conv_rain_prevstep, conv_snow_prevstep, &
                                  air, bgc, canopy, met, bal, rad, rough, soil,&
                                  ssnow, sum_flux, veg
   USE cable_common_module, ONLY : cable_runtime, cable_user, ktau_gl
   USE cable_um_init_subrs_mod, ONLY : um2cable_rr
   USE cable_cbm_module,    ONLY : cbm

   IMPLICIT NONE
        
   REAL, DIMENSION(um1%ROW_LENGTH,um1%ROWS) ::                                &
      LS_RAIN,  & ! IN Large scale rain
      LS_SNOW,  & ! IN Large scale snow
      CON_RAIN, & ! IN Convective rain
      CONV_SNOW,& ! IN Convective snow
      DTL_1,    & ! IN Level 1 increment to T field 
      DQW_1       ! IN Level 1 increment to q field 

   INTEGER :: timestep

   REAL, DIMENSION(um1%land_pts) ::                                           &
      FLAND       ! IN Land fraction on land tiles
   
   REAL, DIMENSION(um1%ROW_LENGTH,um1%ROWS) ::                                &
      !--- Net downward heat flux at surface over land.
      !--- fraction of gridbox (W/m2).
      SURF_HT_FLUX_LAND,                                                      &
      !--- Moisture flux between layers. (kg/m^2/sec).
      !--- FQW(,1) is total water flux from surface, 'E'.
      FQW_1,                                                                  &
      !--- FTL(,K) =net turbulent sensible heat flux into layer K
      !--- from below; so FTL(,1) = surface sensible heat, H.(W/m2)
      FTL_1         

   REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                &
      !___Surface FTL, FQL for land tiles
      FTL_TILE, FQW_TILE,                                                     &
      
      !___(tiled) latent heat flux, melting, stomatatal conductance
     LE_TILE, MELT_TILE,                                                      &
     
     !___ INOUT Surface net radiation on tiles (W/m2)
     RADNET_TILE, &
     EI_TILE,     & ! OUT EI for land tiles.
     ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
     ESOIL_TILE     ! evapotranspiration from soil moisture store (kg/m2/s) 

   REAL, DIMENSION( um1%land_pts,um1%ntiles ) ::                              &
      T1P5M_TILE,    &
      Q1P5M_TILE,    &
      TSTAR_TILE,    &
     
   REAL, DIMENSION(mp) ::                                                     &
      dtlc, & 
      dqwc

   REAL, POINTER :: TFRZ

   !___ 1st call in RUN (!=ktau_gl -see below) 
   LOGICAL, SAVE :: first_cable_call = .TRUE.

   integer :: itest =2
   integer :: iunit=77772
   character(len=44) :: testfname = "/home/599/jxs599/cable_implicit.txt"
    
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
      CALL um2cable_rr( (LS_RAIN+CON_RAIN)*um1%TIMESTEP, met%precip)
      CALL um2cable_rr( (LS_SNOW+CONV_SNOW)*um1%TIMESTEP, met%precip_sn)
      CALL um2cable_rr( dtl_1, dtlc)
      CALL um2cable_rr( dqw_1, dqwc)
      
      !--- conv_rain(snow)_prevstep are added to precip. in explicit call
      CALL um2cable_rr( (CON_RAIN)*um1%TIMESTEP, conv_rain_prevstep)
      CALL um2cable_rr( (CONV_snow)*um1%TIMESTEP, conv_snow_prevstep)
      
      met%precip   =  met%precip + met%precip_sn
      met%tk = met%tk + dtlc
      met%qv = met%qv + dqwc
      met%tvair = met%tk
      met%tvrad = met%tk
 
      canopy%cansto = canopy%oldcansto

      CALL cbm(real(TIMESTEP), air, bgc, canopy, met, bal,  &
           rad, rough, soil, ssnow, sum_flux, veg)
  
        
      CALL implicit_unpack( FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
                            SURF_HT_FLUX_LAND, ECAN_TILE,                      &
                            ESOIL_TILE, EI_TILE, RADNET_TILE,                  &
                            T1P5M_TILE,                                        &
                            Q1P5M_TILE, FLAND, MELT_TILE )
       
      cable_runtime%um_implicit = .FALSE.
  
END SUBROUTINE cable_implicit_driver


!========================================================================= 
!========================================================================= 
!========================================================================= 
        
SUBROUTINE implicit_unpack( FTL_1, FTL_TILE, FQW_1,  FQW_TILE, TSTAR_TILE,     &
                            SURF_HT_FLUX_LAND, ECAN_TILE,                      &
                            ESOIL_TILE, EI_TILE, RADNET_TILE,                  &
                            T1P5M_TILE,                                    &
                            Q1P5M_TILE, FLAND, MELT_TILE )
 
  USE cable_def_types_mod, ONLY : mp
  USE cable_data_module,   ONLY : PHYS
  USE cable_um_tech_mod,   ONLY : um1 ,canopy, rad, soil, ssnow, air
  USE cable_common_module, ONLY : cable_runtime, cable_user

  IMPLICIT NONE
 
  !jhan:these need to be cleaned out to what is actualllly passed
  REAL, DIMENSION(um1%land_pts) ::                                            &
    FLAND          ! IN Land fraction on land tiles
   
  REAL, DIMENSION(um1%ROW_LENGTH,um1%ROWS) ::                                 &
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
    FTL_TILE, FQW_TILE,                 &
    !___(tiled) latent heat flux, melting, stomatatal conductance
    LE_TILE, MELT_TILE,     &
    RADNET_TILE, & ! INOUT Surface net radiation on tiles (W/m2)
    EI_TILE,     & ! OUT EI for land tiles.
    ECAN_TILE,   & ! OUT ECAN for snow-free land tiles
    ESOIL_TILE     ! evapotranspiration from soil moisture store (kg/m2/s)

  REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
    SURF_HTF_T_CAB, &

  REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
    T1P5M_TILE,    &
    Q1P5M_TILE,    &
    TSTAR_TILE,    &


!--- Local vars
  REAL, DIMENSION(mp) ::                                                      &
      fe_dlh,    & !
      fes_dlh,   & !
      fev_dlh      !

  INTEGER :: i,j,l,k,n

  REAL, DIMENSION(um1%land_pts,um1%ntiles) ::                                 &
         !--- Local buffer surface FTL, FQL @ prev dt
         FTL_TILE_old, FQW_TILE_old

  INTEGER:: i_miss = 0
  REAL :: miss = 0.0
   
  REAL, POINTER :: TFRZ
   
!-----------------------------------------------------------------------------
   
  TFRZ => PHYS%TFRZ

  !--- unpack snow vars
  MELT_TILE   = UNPACK(ssnow%smelt, um1%L_TILE_PTS, miss)

  !---preserve fluxes from the previous time step for the coastal grids
  FTL_TILE_old = FTL_TILE
  FQW_TILE_old = FQW_TILE
  !___return fluxes
  FTL_TILE = UNPACK(canopy%fh,  um1%l_tile_pts, miss)
  fe_dlh = canopy%fe/(air%rlam*ssnow%cls)
  WHERE ( fe_dlh .ge. 0.0 ) fe_dlh = MAX ( 1.e-6, fe_dlh )
  WHERE ( fe_dlh .lt. 0.0 ) fe_dlh = MIN ( -1.e-6, fe_dlh )
  fes_dlh = canopy%fes/(air%rlam*ssnow%cls)
  fev_dlh = canopy%fev/air%rlam

  !---update fluxes
  FQW_TILE = UNPACK(fe_dlh, um1%l_tile_pts, miss)

  !___return temp and roughness
  TSTAR_TILE = UNPACK(rad%trad, um1%l_tile_pts, miss)

  !___return miscelaneous
  RADNET_TILE = unpack( canopy%rnet , um1%l_tile_pts, miss)
  ESOIL_TILE = UNPACK(fes_dlh, um1%L_TILE_PTS, miss)
  ECAN_TILE = UNPACK(fev_dlh,  um1%L_TILE_PTS, miss)
  EI_TILE = 0.

  !unpack screen level (1.5m) variables
  !Convert back to K
  t1p5m_tile     = UNPACK(canopy%tscrn+tfrz, um1%L_TILE_PTS, miss)
  q1p5m_tile     = UNPACK(canopy%qscrn, um1%L_TILE_PTS, miss)

  !initialse full land grids and retain coastal grid fluxes
  DO N=1,um1%NTILES
    DO K=1,um1%TILE_PTS(N)
      L = um1%TILE_INDEX(K,N)
      J = (um1%LAND_INDEX(L)-1)/um1%ROW_LENGTH + 1
      I = um1%LAND_INDEX(L) - (J-1)*um1%ROW_LENGTH
      IF( FLAND(L) == 1.0) THEN
        FTL_1(I,J) =  0.0
        FQW_1(I,J) =  0.0
      ELSE
        !retain sea/ice contribution and remove land contribution
        FTL_1(I,J) = FTL_1(I,J) - FLAND(L) * um1%TILE_FRAC(L,N) *             &
                     FTL_TILE_old(L,N)
        FQW_1(I,J) = FQW_1(I,J) - FLAND(L) * um1%TILE_FRAC(L,N) *             &
                     FQW_TILE_old(L,N)
      END IF
      SURF_HT_FLUX_LAND(I,J) = 0.0
    END DO
  END DO

  DO N=1,um1%NTILES
    DO K=1,um1%TILE_PTS(N)
      L = um1%TILE_INDEX(K,N)
      J = (um1%LAND_INDEX(L)-1)/um1%ROW_LENGTH + 1
      I = um1%LAND_INDEX(L) - (J-1)*um1%ROW_LENGTH
      FTL_1(I,J) = FTL_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FTL_TILE(L,N)
      FQW_1(I,J) = FQW_1(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*FQW_TILE(L,N)
      SURF_HT_FLUX_LAND(I,J) = SURF_HT_FLUX_LAND(I,J) +                       &
                               FLAND(L)*um1%TILE_FRAC(L,N) *                  &
                               SURF_HTF_T_CAB(L,N)
    END DO
  END DO

END SUBROUTINE implicit_unpack
