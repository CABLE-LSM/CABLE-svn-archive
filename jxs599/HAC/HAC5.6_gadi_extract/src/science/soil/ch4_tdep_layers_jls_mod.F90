! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE CH4_TDEP_LAYERS----------------------------------------

! Description:
!     Calculates methane emissions from wetland area.

MODULE ch4_tdep_layers_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CH4_TDEP_LAYERS_MOD'

CONTAINS

SUBROUTINE ch4_tdep_layers(npnts, soil_pts, timestep, dim_cs1, soil_index,    &
                       sm_levels, dim_cslayer, tsoil, cs, resp_s, resp_s_tot, &
                       npp, f_wetl, kaps_weight, dzsoil, ztot, t0_ch4,        &
                       ch4_substrate, const_tdep_cs, const_tdep_npp,          &
                       const_tdep_resps, fch4_wetl_cs, fch4_wetl_npp,         &
                       fch4_wetl_resps) 

USE water_constants_mod, ONLY: tm

USE jules_soil_biogeochem_mod, ONLY:                                          &
  ch4_substrate_npp,  ch4_substrate_soil, ch4_substrate_soil_resp,            &
  l_ch4_interactive, tau_resp, l_layeredC,                                    &
  const_ch4_cs, const_ch4_npp, const_ch4_resps

USE parkind1, ONLY: jprb, jpim
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN):
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  npnts,                                                                      &
    ! Number of gridpoints.
  soil_pts,                                                                   &
    ! Number of soil points.
  dim_cs1,                                                                    &
    ! Number of soil carbon pools
  soil_index(npnts),                                                          &
    ! Array of soil points.
  sm_levels,                                                                  &
    ! Soil layers
  dim_cslayer,                                                                &
    ! Size of soil carbon vertical dimension
  ch4_substrate
    ! Choice of substrate for CH4 production.

REAL, INTENT(IN) ::                                                           &
  tsoil(npnts,sm_levels),                                                     &
    ! Layered soil temperature (K)
  resp_s(npnts,dim_cslayer,dim_cs1),                                          &
    ! Soil respiration in pools (kg C/m2/s).
  resp_s_tot(npnts),                                                          &
    ! Soil respiration total (kg C/m2/s).
  npp(npnts),                                                                 &
    ! Gridbox mean net primary productivity (kg C/m2/s).
  f_wetl(npnts),                                                              &
    ! Wetland fraction
  timestep,                                                                   &
    ! Model timestep (s)
  t0_ch4,                                                                     &
    !  T0 value (zero celsius in kelvin)
  const_tdep_cs,                                                              &
    ! T and Q10(0) dependent function
  const_tdep_npp,                                                             &
    ! T and Q10(0) dependent function
  const_tdep_resps,                                                           &
    ! T and Q10(0) dependent function
  ztot(sm_levels),                                                            &
    ! Depth at center of each soil layer (m)
  dzsoil(sm_levels),                                                          &
    ! Thicknesses of soil layers (m)
  kaps_weight(dim_cs1)
    ! Variable for weighting by kappas

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT):
!-----------------------------------------------------------------------------
REAL, INTENT(INOUT) ::                                                        &
  fch4_wetl_cs(npnts),                                                        &
    ! Scaled methane flux (soil carbon substrate) (kg C/m2/s)
  fch4_wetl_npp(npnts),                                                       &
    ! Scaled methane flux (npp substrate) (kg C/m2/s)
  fch4_wetl_resps(npnts),                                                     &
    ! Scaled methane flux (soil respiration substrate) (kg C/m2/s)
  cs(npnts,dim_cslayer,dim_cs1)
    ! Soil carbon
    ! For RothC (dim_cs1=4), the pools are DPM, RPM, biomass and humus
    ! (kg C/m2).

!-----------------------------------------------------------------------------
! Local scalars:
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i, j, k, n

REAL ::                                                                       &
  cor_acc_cs,                                                                 &
    ! Accumulated correction based on not exhausting the soil carbon.
    ! With soil carbon as substrate.
  cor_acc_npp,                                                                &
    ! Accumulated correction based on not exhausting the soil carbon.
    ! With NPP as substrate.
  cor_acc_resp,                                                               &
    ! Accumulated correction based on not exhausting the soil carbon.
    ! With soil respiration as substrate.
  cs_tmp,                                                                     &
    ! Working variable to assign appropriate soil carbon quantity depending
    ! on l_layeredC TRUE/FALSE.
  q10t_ch4_cs,                                                                &
    ! Q10 value at T
  q10t_ch4_npp,                                                               &
    ! Q10 value at T
  q10t_ch4_resps,                                                             &
    ! Q10 value at Tcs
  resp_tmp   
    ! Working variable to assign appropriate soil respiration depending on
    ! l_layeredC TRUE/FALSE

!-----------------------------------------------------------------------------
! Local arrays:
!-----------------------------------------------------------------------------
REAL ::                                                                       &
  ch4_tmp_lys(sm_levels),                                                     &
    ! Working variable in calculation of interactive CH4.
  fch4_wetl_cs_lys(sm_levels),                                                &
    ! Working variable in calculation of CH4 with soil carbon as substrate
  fch4_wetl_resps_lys(sm_levels),                                             &
    ! Working variable in calculation of CH4 with respiration as substrate
  fch4_wetl_npp_lys(sm_levels)
    ! Working variable in calculation of CH4 with NPP as substrate

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CH4_TDEP_LAYERS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Calculate an effective soil carbon for wetland methane emission.
!-----------------------------------------------------------------------------

DO j = 1,soil_pts
  i = soil_index(j)
  ! Set accumulated variables to zero.
  cor_acc_cs         = 0.0
  cor_acc_npp        = 0.0
  cor_acc_resp       = 0.0

  IF ( ANY(tsoil(i,:) > (tm-0.1)) .AND. f_wetl(i) > 0.0 ) THEN
    !-------------------------------------------------------------------------
    ! CH4 flux using NPP as substrate
    fch4_wetl_npp_lys(:) = 0.0
    IF ( npp(i) > 0.0 ) THEN !This is the same for layered and non-layered
      DO n = 1,sm_levels
        IF ( tsoil(i,n) > (tm-0.1) ) THEN
          q10t_ch4_npp         = EXP( const_tdep_npp / tsoil(i,n) )
          fch4_wetl_npp_lys(n) = const_ch4_npp * npp(i) * f_wetl(i) *         &
                                 q10t_ch4_npp**( 0.1 *                        &
                                 ( tsoil(i,n) - t0_ch4 ) ) *                  &
                                 EXP( -tau_resp * ztot(n) ) *                 &
                                 dzsoil(n) * tau_resp
        END IF
      END DO
    END IF !npp > 0.0

    DO k = 1,dim_cs1
      fch4_wetl_cs_lys(:)    = 0.0
      fch4_wetl_resps_lys(:) = 0.0
      DO n = 1,sm_levels
        IF ( tsoil(i,n) > (tm-0.1) ) THEN
          !-------------------------------------------------------------------
          ! CH4 flux using soil carbon as substrate.
          IF ( l_layeredC ) THEN
            cs_tmp = cs(i,n,k)
          ELSE
            cs_tmp = cs(i,1,k)
          END IF
          q10t_ch4_cs         = EXP( const_tdep_cs / tsoil(i,n) )
          fch4_wetl_cs_lys(n) = const_ch4_cs * cs_tmp * f_wetl(i) *           &
                                q10t_ch4_cs**( 0.1 *                          &
                                ( tsoil(i,n) - t0_ch4 ) ) *                   &
                                EXP( -tau_resp * ztot(n) ) * kaps_weight(k)

          ! Weight so that the non-layered case is split between layers.
          IF ( .NOT. l_layeredC ) THEN
            fch4_wetl_cs_lys(n) = fch4_wetl_cs_lys(n) * dzsoil(n) * tau_resp
          END IF

          !-------------------------------------------------------------------
          ! CH4 flux using soil respiration as substrate.
          IF ( resp_s_tot(i) > 0.0 ) THEN
            IF ( l_layeredC ) THEN
              resp_tmp = resp_s(i,n,k)
            ELSE
              resp_tmp = resp_s(i,1,k)
            END IF
            q10t_ch4_resps         = EXP( const_tdep_resps / tsoil(i,n) )
            fch4_wetl_resps_lys(n) = const_ch4_resps * f_wetl(i) *            &
                                     resp_tmp * q10t_ch4_resps                &
                                     **( 0.1 * ( tsoil(i,n) - t0_ch4 ) )

            ! Weight so that the non-layered case is split between layers
            ! and descreases with depth
            IF ( .NOT. l_layeredC ) THEN
              fch4_wetl_resps_lys(n) = fch4_wetl_resps_lys(n) * tau_resp *    &
                                       EXP( -tau_resp * ztot(n) ) * dzsoil(n)
            END IF
          END IF !resp_s_tot > 0
        END IF !tsoil
      END DO !sm_levels

      !-----------------------------------------------------------------------
      ! Interactive methane (coupled into carbon cycle)
      !-----------------------------------------------------------------------
      IF ( l_ch4_interactive ) THEN 

        IF ( ch4_substrate == ch4_substrate_soil ) THEN
          ! Substrate 1: soil carbon.
          ch4_tmp_lys(:)  = fch4_wetl_cs_lys(:)
        ELSE IF ( ch4_substrate == ch4_substrate_npp ) THEN
          ! Substrate 2: NPP
          ch4_tmp_lys(:)  = fch4_wetl_npp_lys(:) * kaps_weight(k)
        ELSE IF ( ch4_substrate == ch4_substrate_soil_resp ) THEN
          ! Substrate 3: respiration
          ch4_tmp_lys(:)  = fch4_wetl_resps_lys(:)
        END IF

        !Prevent it from exhausting any of the soil carbon pools
        IF ( l_layeredC ) THEN
          DO n = 1,dim_cslayer
            ch4_tmp_lys(n) = MIN( ch4_tmp_lys(n), cs(i,n,k) / timestep )
          END DO
        ELSE
          IF ( SUM(ch4_tmp_lys) > cs(i,1,k) / timestep ) THEN
            ch4_tmp_lys(:) = ch4_tmp_lys(:) * cs(i,1,k) /                     &
                             ( timestep * SUM(ch4_tmp_lys) )
          END IF
        END IF !l_layeredC

        ! ACCUMULATE ANY CORRECTIONS SO WE CAN TAKE THEM OFF THE CH4
        ! EMISSIONS AT THE END
        IF ( ch4_substrate == ch4_substrate_soil ) THEN
          cor_acc_cs   = cor_acc_cs   + SUM( fch4_wetl_cs_lys - ch4_tmp_lys )
        ELSE IF ( ch4_substrate == ch4_substrate_npp ) THEN
          cor_acc_npp  = cor_acc_npp  + SUM( fch4_wetl_npp_lys *              &
                                             kaps_weight(k) - ch4_tmp_lys )
        ELSE IF ( ch4_substrate == ch4_substrate_soil_resp ) THEN
          cor_acc_resp = cor_acc_resp +                                       &
                         SUM( fch4_wetl_resps_lys - ch4_tmp_lys )
        END IF

        ! Subtract CH4 emissions from soil carbon.
        IF ( l_layeredC ) THEN
          DO n = 1,dim_cslayer
            cs(i,n,k) = MAX(cs(i,n,k) - ch4_tmp_lys(n) * timestep, 0.0)
          END DO !dim_cslayer
        ELSE
          cs(i,1,k) = MAX(cs(i,1,k) - SUM(ch4_tmp_lys) * timestep, 0.0)
        END IF !l_layeredC
      END IF !l_ch4_interactive

      !-----------------------------------------------------------------------
      ! Update wetland methane flux.
      !-----------------------------------------------------------------------
      fch4_wetl_cs(i)    = fch4_wetl_cs(i) + SUM( fch4_wetl_cs_lys )
      fch4_wetl_resps(i) = fch4_wetl_resps(i) + SUM( fch4_wetl_resps_lys )
    END DO !dim_cs1

    fch4_wetl_cs(i)    = fch4_wetl_cs(i)          - cor_acc_cs
    fch4_wetl_npp(i)   = SUM( fch4_wetl_npp_lys ) - cor_acc_npp
    fch4_wetl_resps(i) = fch4_wetl_resps(i)       - cor_acc_resp

  END IF !tsoil and fwetl
END DO !land_points

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ch4_tdep_layers
END MODULE ch4_tdep_layers_mod
