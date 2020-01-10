! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE p_s_parms

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module containing plant and soil variables (plus a few others)

! Note that all variables with suffix _soilt have soil tiles as their
! 2nd dimension. Cross-reference allocate_jules_arrays for more details.
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  !Used in both the UM and JULES
  bexp_soilt(:,:,:),                                                          &
    !  Exponent for soil moisture characteristic functions
    !    Clapp-Hornberger model: b is the Clapp-Hornberger exponent
    !    van Genuchten model: b=1/(n-1)  (metres)
  sathh_soilt(:,:,:),                                                         &
    !  Parameter for soil moisture characteristic functions
    !  Clapp-Hornberger model: sathh is the saturated soil water pressure (m)
    !  van Genuchten model: sathh=1/alpha
  hcap_soilt(:,:,:),                                                          &
    !  Soil heat capacity (J/K/m3)
  hcon_soilt(:,:,:),                                                          &
    !  Soil thermal conductivity (W/m/K)
  satcon_soilt(:,:,:),                                                        &
    !  Saturated hydraulic conductivity (kg/m2/s)
  smvccl_soilt(:,:,:),                                                        &
    !  Critical volumetric SMC (cubic m per cubic m of soil)
  smvcst_soilt(:,:,:),                                                        &
    !  Volumetric saturation point (m3/m3 of soil)
  smvcwt_soilt(:,:,:),                                                        &
    !  Volumetric wilting point (cubic m per cubic m of soil)
  v_close_pft(:,:,:),                                                         &
    ! Volumetric soil moisture
    ! concentration below which stomata are fully closed.
    ! (m3 H2O/m3 soil).
    ! If l_use_pft_psi=F, this variable is not used, and
    ! the soil ancil variable sm_wilt is used instead.
  v_open_pft(:,:,:),                                                          &
    ! concentration above which stomatal aperture is not limited
    ! by soil water (m3 H2O/m3 soil).
    ! If l_use_pft_psi=F, this variable is not used, and
    ! the soil ancillary variable sm_crit and
    ! the pft variable fsmc_p0 are used instead.
  clay_soilt(:,:,:)
    !  Fraction of clay (kg clay/kg soil).

  !JULES only
#if !defined(UM_JULES)
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  albsoil_soilt(:,:),                                                         &
    !  Soil albedo
  albobs_sw_gb(:),                                                            &
    !  Obs SW albedo
  albobs_vis_gb(:),                                                           &
    !  Obs VIS albedo
  albobs_nir_gb(:),                                                           &
    !  Obs NIR albedo
  catch_surft(:,:),                                                           &
    !  Surface/canopy water capacity of snow-free land tiles (kg/m2)
  catch_snow_surft(:,:),                                                      &
    !  Snow interception capacity (kg/m2)
  cosz_ij(:,:),                                                               &
    !  Cosine of the zenith angle
  infil_surft(:,:),                                                           &
    !  Maximum possible surface infiltration for tiles (kg/m2/s)
  z0_surft(:,:),                                                              &
    !  Surface roughness on tiles (m).
  z0h_bare_surft(:,:),                                                        &
    !  Surface thermal roughness on tiles before allowance for snow
    !  cover (m).
  z0m_soil_gb(:),                                                             &
    !  Bare soil roughness, for momentum (m).
  sthu_soilt(:,:,:),                                                          &
    ! Unfrozen soil moisture content of the layers as a fraction of saturation
  sthf_soilt(:,:,:),                                                          &
    !  Frozen soil moisture content of the layers as a fraction of saturation.
  sthu_min_soilt(:,:,:),                                                      &
    ! Minimum unfrozen water content for each layer. Used to normalise
    ! thaw depth calculation based on unfrozen water content fraction.
!-----------------------------------------------------------------------------
!    ECOSSE ancillaries.
!-----------------------------------------------------------------------------
    soil_ph_soilt(:,:,:)
      ! Soil pH, defined on soil layers.
#endif

END MODULE p_s_parms

