! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing all of the prognostic variables

MODULE prognostics

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Description
! Module containing all of the prognostic variables,
! i.e.those required to be kept from one timestep
! to another. Variables all appear in a model dump - NOT AT PRESENT!
! And some of these are not prognostics (eg smc)....

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

! Declarations:

INTEGER, ALLOCATABLE ::                                                       &
  nsnow(:,:)      !  Number of snow layers on ground on tiles

INTEGER, ALLOCATABLE ::                                                       &
  nsnow_surft(:,:)      !  Number of snow layers on ground on tiles

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  tsoil_deep_gb(:,:),                                                         &
              !  Deep soil temperatures (K)
  sice_surft(:,:,:),                                                          &
                  ! Snow layer ice mass on tiles (kg/m2)
  sliq_surft(:,:,:),                                                          &
                  ! Snow layer liquid mass on tiles (kg/m2)
  snowdepth_surft(:,:),                                                       &
                  ! Snow depth on ground on tiles (m)
  tsnow_surft(:,:,:),                                                         &
                  ! Snow layer temperature (K)
  rgrainl_surft(:,:,:),                                                       &
                  ! Snow layer grain size on tiles (microns)
  rho_snow_grnd_surft(:,:),                                                   &
                  ! Snowpack bulk density (kg/m3)

  !The status of this as a prognostic needs to be reviewed. It is in the UM
  !dump but not the JULES-standalone dump.
  rho_snow_surft(:,:,:),                                                      &
                  ! Snow layer densities (m)

  ds_surft(:,:,:),                                                            &
                  ! Snow layer thickness (m)
  wood_prod_fast_gb(:),                                                       &
!               ! Fast-turnover wood product C pool (kg m-2).
    wood_prod_med_gb(:),                                                      &
                ! Medium-turnover wood product C pool (kg m-2).
    wood_prod_slow_gb(:),                                                     &
                ! Slow-turnover wood product C pool (kg m-2).
    frac_agr_prev_gb(:),                                                      &
                ! Agricultural fraction from previous TRIFFID call
    frac_past_prev_gb(:),                                                     &
                ! Pasture fraction from previous TRIFFID call
    n_inorg_soilt_lyrs(:,:,:),                                                &
                ! Gridbox Inorganic N pool on soil levels (kg N/m2)
    n_inorg_gb(:),                                                            &
                ! Gridbox Inorganic N pool (kg N/m2)
    n_inorg_avail_pft(:,:,:),                                                 &
                ! Available inorganic N for PFTs (depends on roots) (kg N/m2)
    ns_pool_gb(:,:,:),                                                        &
                !  Soil Organic Nitrogen (kg N/m2)
!               !  If dim_cs1=1, there is a single soil C pool.
!               !  If dim_cs1=4, the pools are:
!               !  1  decomposable plant material
!               !  2  resistant plant material
!               !  3  biomass
!               !  4  humus
    triffid_co2_gb(:)
                ! Atmospheric CO2 fluxes from TRIFFID (kgC/m2/yr)
                ! exudates + wood product pool flux + harvest for adding to the
                ! atmosphere in interactive CO2 runs

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  canht_pft(:,:),                                                             &
              !  Canopy height (m)
  canopy_surft(:,:),                                                          &
              !  Surface/canopy water for snow-free land tiles (kg/m2)
  canopy_gb(:),                                                               &
              !  Gridbox canopy water content (kg/m2)
  cs_pool_soilt(:,:,:,:),                                                     &
              !  Soil carbon (kg C/m2)
!               !  If dim_cs1=1, there is a single soil C pool per layer.
!               !  If dim_cs1=4, the pools are:
!               !  1  decomposable plant material
!               !  2  resistant plant material
!               !  3  biomass
!               !  4  humus

    di_ncat_sicat(:,:,:),                                                     &
                !  "Equivalent thickness" of sea-ice catagories (m)
    k_sice_sicat(:,:,:),                                                      &
                !  Sea ice effective conductivity (2*kappai/de)
    gc_surft(:,:),                                                            &
                !  Stomatal" conductance to evaporation for land tiles(m/s)
    gs_gb(:),                                                                 &
                !  "Stomatal" conductance to evaporation (m/s)
    lai_pft(:,:),                                                             &
                !  LAI of plant functional types
    rgrain_surft(:,:),                                                        &
                !  Snow surface grain size on tiles (microns)
    smc_soilt(:,:),                                                           &
                !  Soil moisture in a layer at the surface (kg/m2).
!      Note that SMC is used twice:
!      1) SF_EXPL and SF_IMPL2 use it to return the available water in the total
!         soil column (as calculated in PHYSIOL)
!      2) HYDROL uses it to return the available water in a layer of a given
!         depth (as calculated in SOILMC)
    smcl_soilt(:,:,:),                                                        &
                !  Soil moisture content of layers on soil tiles (kg/m2)
    snow_surft(:,:),                                                          &
                !  Lying snow on tiles (kg/m2)
                !  If can_model=4,
                !    snow_surft is the snow on the canopy
                !    snow_grnd is the snow on the ground beneath canopy
                !  If can_model/=4, snow_surft is the total snow.
    snow_grnd_surft(:,:),                                                     &
                !  Snow on the ground (kg/m2)
                !  This is the snow beneath the canopy and is only
                !  used if can_model=4.
    snow_mass_ij(:,:),                                                        &
                !  Gridbox snowmass (kg/m2)
    snow_mass_sea_sicat(:,:,:),                                               &
                !  Snow on category sea-ice (kg/m2)
    soot_ij(:,:),                                                             &
                !  Snow soot content (kg/kg)
    t_soil_soilt(:,:,:),                                                      &
                !  Sub-surface temperature on layers and soil tiles (K)
    t_soil_soilt_acc(:,:,:),                                                  &
                !  Sub-surface temperature on layers and soil tiles
                !  accumulated over TRIFFID timestep (K).
                !  Gives mean temperature over TRIFFID timestep.
    ti_sicat(:,:,:),                                                          &
                !  Sea-ice surface layer
    tstar_surft(:,:),                                                         &
                !  Tile surface temperatures (K)
    tsurf_elev_surft(:,:),                                                    &
                !  Tiled land-ice bedrock subsurface temperatures (K)
    z0msea_ij(:,:),                                                           &
                !  Sea-surface roughness length for momentum (m).
    routestore(:,:)
                ! channel storage (kg). This is defined at all
!               ! points on the routing grid, both land and sea, although it is
!               ! only used at land points. I've used this full grid
!               ! approach mainly for compatability with current
!               ! code in the UM, and also because the grid is needed
!               ! for the calculated flow, so doing everything on grid
!               ! simplifies i/o as we don't need to deal with both
!               ! vector and grid i/o.

INTEGER, ALLOCATABLE ::                                                       &
  seed_rain(:)
              ! Seeding number for subdaily rainfall for use
              ! in IMOGEN or when l_daily_disagg = T

LOGICAL :: l_broadcast_soilt = .FALSE.
                      !Switch to broadcast model state around all soil tiles.
                      !Only has an effect if l_tile_soil is true.
                      !Does not do anything with ancils- this is done by
                      !l_broadcast_ancils

END MODULE prognostics
