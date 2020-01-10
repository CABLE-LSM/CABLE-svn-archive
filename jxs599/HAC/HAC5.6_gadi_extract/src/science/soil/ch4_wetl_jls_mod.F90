! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE CH4_WETL-----------------------------------------------

! Description:
!     Calculates methane emissions from wetland area.

MODULE ch4_wetl_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CH4_WETL_MOD'

CONTAINS
SUBROUTINE ch4_wetl(npnts, soil_pts, dim_cs1, soil_index, sm_levels,          &
                    tsoil_d, cs_ch4, tsoil, cs, resp_s, npp, f_wetl,          &
                    fch4_wetl, fch4_wetl_cs, fch4_wetl_npp, fch4_wetl_resps,  &
                    timestep, l_ch4_tlayered )

!Use in relevant subroutines
USE ch4_tdep_mod, ONLY: ch4_tdep
USE ch4_tdep_layers_mod, ONLY: ch4_tdep_layers

USE jules_soil_biogeochem_mod, ONLY:                                          &
  kaps_roth, soil_model_rothc, soil_bgc_model, ch4_substrate,                 &
  ch4_substrate_npp,  ch4_substrate_soil, ch4_substrate_soil_resp,            &
  t0_ch4, q10_ch4_cs, q10_ch4_npp, q10_ch4_resps

USE ancil_info, ONLY: dim_cslayer

USE jules_soil_mod, ONLY: dzsoil

USE parkind1, ONLY: jprb, jpim
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN):
!-----------------------------------------------------------------------------
LOGICAL, INTENT(IN) ::                                                        &
  l_ch4_tlayered

INTEGER, INTENT(IN) ::                                                        &
  npnts,                                                                      &
    ! Number of gridpoints.
  soil_pts,                                                                   &
    ! Number of soil points.
  dim_cs1,                                                                    &
    ! Number of soil carbon pools
  soil_index(npnts),                                                          &
    ! Array of soil points.
  sm_levels
    ! Soil layers

REAL, INTENT(IN) ::                                                           &
  tsoil(npnts,sm_levels),                                                     &
    ! Layered soil temperature (K).
  tsoil_d(npnts),                                                             &
    ! Diagnosed soil temp to 1 metre (K).
  cs_ch4(npnts),                                                              &
    ! Soil carbon used in CH4 wetlands if single-pool model is used (kg C/m2).
  resp_s(npnts,dim_cslayer,dim_cs1),                                          &
    ! Soil respiration in pools (kg C/m2/s).
  npp(npnts),                                                                 &
    ! Gridbox mean net primary productivity (kg C/m2/s).
  f_wetl(npnts),                                                              &
    ! Wetland fraction.
  timestep
    ! Model timestep (s).

!-----------------------------------------------------------------------------
! Arguments with INTENT(INOUT):
!-----------------------------------------------------------------------------
REAL, INTENT(INOUT) ::                                                        &
  fch4_wetl(npnts),                                                           &
    ! Scaled methane flux (as used in atmos chem)
    ! (10^-9 kg C/m2/s).
  fch4_wetl_cs(npnts),                                                        &
    ! Scaled methane flux (soil carbon substrate) (kg C/m2/s).
  fch4_wetl_npp(npnts),                                                       &
    ! Scaled methane flux (npp substrate) (kg C/m2/s).
  fch4_wetl_resps(npnts),                                                     &
    ! Scaled methane flux (soil respiration substrate) (kg C/m2/s).
  cs(npnts,dim_cslayer,dim_cs1)
    ! Soil carbon
    ! For RothC (dim_cs1=4), the pools are DPM, RPM, biomass and humus
    ! (kg C/m2).

!-----------------------------------------------------------------------------
! Local scalar variables:
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i, j, k, n

REAL ::                                                                       &
  const_tdep_cs,                                                              &
    ! T and Q10(0) dependent function.
  const_tdep_npp,                                                             &
    ! T and Q10(0) dependent function.
  const_tdep_resps,                                                           &
    ! T and Q10(0) dependent function.
  sumkaps
    ! Sum of kaps_roth values.

!-----------------------------------------------------------------------------
! Local array variables:
!-----------------------------------------------------------------------------
REAL ::                                                                       &
  cs_eff(npnts),                                                              &
    ! Effective soil carbon (kg C/m2).
  resp_s_tot(npnts),                                                          &
    ! Soil respiration total (kg C/m2/s).
  ztot(sm_levels),                                                            &
    ! Depth at center of each soil layer (m).
  kaps_weight(dim_cs1)                                                         
    ! Working variable for weighting by kappas.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CH4_WETL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Calculate an effective soil carbon for wetland methane emission.
!-----------------------------------------------------------------------------

IF (l_ch4_tlayered) THEN
  ! Calculate soil layer depths from layer thicknesses.
  ztot(1) = dzsoil(1) * 0.5
  DO n = 2,sm_levels
    ztot(n) = ztot(n-1) + 0.5 * ( dzsoil(n-1) + dzsoil(n) )
  END DO
END IF

IF ( dim_cs1 > 1 ) THEN
  sumkaps = 0.0
  DO k = 1,dim_cs1
    sumkaps = sumkaps + kaps_roth(k)
  END DO
  IF ( .NOT. l_ch4_tlayered ) THEN
    ! Weight each pool by specific respiration rate.
    DO j = 1,soil_pts
      i = soil_index(j)
      cs_eff(i) = 0.0
      DO k = 1,dim_cs1
        DO n = 1,dim_cslayer
          cs_eff(i) = cs_eff(i) + cs(i,n,k) * kaps_roth(k)
        END DO
      END DO
      cs_eff(i) = cs_eff(i) / sumkaps
    END DO
  ELSE !l_ch4_tlayered=T
    ! Calculate weighting according to kappas.
    DO k = 1,dim_cs1
      kaps_weight(k) = kaps_roth(k) / sumkaps
    END DO
  END IF !l_ch4_tlayered
    
ELSE !dim_cs1==1

  IF ( .NOT. l_ch4_tlayered ) THEN
    ! Use the single soil carbon pool.
    DO j = 1,soil_pts
      i = soil_index(j)
      cs_eff(i) = cs_ch4(i)
    END DO
  ELSE
    ! l_ch4_tlayered=T
    kaps_weight(1) = 1.0
  END IF !l_ch4_tlayered

END IF  !dim_cs1

! Calculate total soil respiration
DO j = 1,soil_pts
  i = soil_index(j)
  resp_s_tot(i) = 0.0
  DO k = 1,dim_cs1
    DO n = 1,dim_cslayer
      IF ( resp_s(i,n,k) > 0.0 )                                              &
        resp_s_tot(i) = resp_s_tot(i) + resp_s(i,n,k)
    END DO
  END DO
END DO

!-----------------------------------------------------------------------------
! Calculate scaled wetland methane emission.
!-----------------------------------------------------------------------------
const_tdep_cs = t0_ch4 * LOG(q10_ch4_cs)
const_tdep_npp = t0_ch4 * LOG(q10_ch4_npp)
const_tdep_resps = t0_ch4 * LOG(q10_ch4_resps)

IF ( l_ch4_tlayered ) THEN
  CALL ch4_tdep_layers(npnts, soil_pts, timestep, dim_cs1, soil_index,        &
                       sm_levels, dim_cslayer, tsoil, cs, resp_s, resp_s_tot, &
                       npp, f_wetl, kaps_weight, dzsoil, ztot, t0_ch4,        &
                       ch4_substrate, const_tdep_cs, const_tdep_npp,          &
                       const_tdep_resps, fch4_wetl_cs, fch4_wetl_npp,         &
                       fch4_wetl_resps) 
ELSE !NOT l_ch4_tlayered
  CALL ch4_tdep(npnts, soil_pts, soil_index, tsoil_d, cs_eff, resp_s_tot,     &
                npp, f_wetl, t0_ch4, const_tdep_cs, const_tdep_npp,           &
                const_tdep_resps, fch4_wetl_cs, fch4_wetl_npp, fch4_wetl_resps)
END IF !l_ch4_tlayered

IF ( ch4_substrate == ch4_substrate_soil ) THEN
  fch4_wetl(:) = 1.0e9 * fch4_wetl_cs(:)
ELSE IF ( ch4_substrate == ch4_substrate_npp ) THEN
  fch4_wetl(:) = 1.0e9 * fch4_wetl_npp(:)
ELSE IF ( ch4_substrate == ch4_substrate_soil_resp ) THEN
  fch4_wetl(:) = 1.0e9 * fch4_wetl_resps(:)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ch4_wetl
END MODULE ch4_wetl_mod
