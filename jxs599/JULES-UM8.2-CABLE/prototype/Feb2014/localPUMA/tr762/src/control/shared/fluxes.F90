! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing surface fluxes.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Land

MODULE fluxes

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! anthrop_heat is required by both the UM and standalone configurations
!-----------------------------------------------------------------------------
  REAL, ALLOCATABLE ::                                                        &
    anthrop_heat(:,:)      ! Additional heat source on tiles
!                          ! used for anthropgenic urban
!                          ! heat source (W/m2)

  REAL, ALLOCATABLE ::                                                        &
    surf_ht_store(:,:)     !   Diagnostic to store values
!                          !   of C*(dT/dt) during calculation
!                          !   of energy balance

  REAL, ALLOCATABLE, SAVE ::                                                  &
    sw_sice(:,:),                                                             &
                           !  Net SW on sea ice categories
    sw_sice_rts(:,:),                                                         &
                           !  Net SW on sea ice categories on the
!                          !     radiative timestep
    alb_sice(:,:,:)        !  Albedo of sea ice categories

!-----------------------------------------------------------------------------
! Everything else is required only in standalone configuration
!-----------------------------------------------------------------------------
#if !defined(UM_JULES)
  REAL, ALLOCATABLE ::                                                        &
    ALB_TILE(:,:,:),                                                          &
      !   Albedo for surface tiles
      !     (:,:,1) direct beam visible
      !     (:,:,2) diffuse visible
      !     (:,:,3) direct beam near-IR
      !     (:,:,4) diffuse near-IR
    e_sea(:,:),                                                               &
      !   Evaporation from sea times leads fraction. Zero over land
      !                                (kg per square metre per sec)
    ECAN(:,:),                                                                &
      !   Gridbox mean evaporation from canopy/surface store (kg/m2/s)
      !     Zero over sea
    ECAN_TILE(:,:),                                                           &
      !   Canopy evaporation from for snow-free land tiles
    EI(:,:),                                                                  &
      !   Sublimation from lying snow or sea-ice (kg/m2/s)
    EI_TILE(:,:),                                                             &
      !   EI for land tiles
    ESOIL(:,:),                                                               &
      !   Surface evapotranspiration from soil moisture store (kg/m2/s)
    ESOIL_TILE(:,:),                                                          &
      !   ESOIL for snow-free land tiles
    EXT(:,:),                                                                 &
      !   Extraction of water from each soil layer (kg/m2/s)
    fqw_1(:,:),                                                               &
      !   Moisture flux between layers (kg per square metre per sec)
      !                         FQW(,1) is total water flux from surface, 'E'
    FQW_TILE(:,:),                                                            &
      !   Surface FQW for land tiles
    FQW_ICE(:,:,:),                                                           &
      !   Surface FQW for sea-ice
    fsmc(:,:),                                                                &
      !   Moisture availability factor.
    ftl_1(:,:),                                                               &
      !   FTL(,K) contains net turbulent sensible heat flux into layer K from below
      !   so FTL(,1) is the surface sensible heat, H.(W/m2)
    FTL_ICE(:,:,:),                                                           &
      !   Surface FTL for sea-ice
    ftl_tile(:,:),                                                            &
      !   Surface FTL for land tiles
    H_SEA(:,:),                                                               &
      !   Surface sensible heat flux over sea times leads fraction (W/m2)
    HF_SNOW_MELT(:),                                                          &
      !   Gridbox snowmelt heat flux (W/m2)
    LAND_ALBEDO(:,:,:),                                                       &
      !   GBM albedo
      !     (,:,1) direct beam visible
      !     (::,,2) diffuse visible
      !     (:,:,3) direct beam near-IR
      !     (:,:,4) diffuse near-IR
    LATENT_HEAT(:,:),                                                         &
      !   Surface latent heat flux, +ve upwards (Watts per sq m)
    LE_TILE(:,:),                                                             &
      !   Surface latent heat flux for land tiles
    MELT_TILE(:,:),                                                           &
      !   Snowmelt on land tiles (kg/m2/s)
    SEA_ICE_HTF(:,:,:),                                                       &
      !   Heat flux through sea-ice (W/m2, positive downwards)
    SICE_MLT_HTF(:,:,:),                                                      &
      !   Heat flux due to melting of sea-ice (Watts per sq metre)
    SNOMLT_SUB_HTF(:),                                                        &
      !   Sub-canopy snowmelt heat flux (W/m2)
    SNOMLT_SURF_HTF(:,:),                                                     &
      !   Heat flux required for surface melting of snow (W/m2)
    SNOW_MELT(:),                                                             &
      !   Snowmelt on land points (kg/m2/s)
    SNOWMELT(:,:),                                                            &
      !   Snowmelt (kg/m2/s)
    SUB_SURF_ROFF(:),                                                         &
      !   Sub-surface runoff (kg/m2/s)
    SURF_HT_FLUX(:,:),                                                        &
      !   Net downward heat flux at surface over land and sea-ice fraction of gridbox (W/m2)
    surf_htf_tile(:,:),                                                       &
      !   Surface heat flux on land tiles (W/m2)
    SURF_ROFF(:),                                                             &
      !   Surface runoff (kg/m2/s)
    RADNET_TILE(:,:),                                                         &
      !   Surface net radiation on tiles ( W/m2)
    TAUX_1(:,:),                                                              &
      !   W'ly component of surface wind stress (N/sq m)
    TAUY_1(:,:),                                                              &
      !   S'ly component of surface wind stress (N/sq m)
      !                                 On V-grid; comments as per TAUX
    TOT_TFALL(:),                                                             &
      !   Total throughfall (kg/m2/s)
    tstar(:,:),                                                               &
      !   GBM surface temperature (K)
    emis_tile(:,:),                                                           &
      !   Tile emissivity
    sw_tile(:,:)
      !   Surface net shortwave on tiles (W/m2)
#endif

END MODULE fluxes
