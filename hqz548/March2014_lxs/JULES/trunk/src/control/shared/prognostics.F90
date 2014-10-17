! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing all of the prognostic variables

MODULE prognostics

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


  INTEGER, ALLOCATABLE ::                                               &
    nsnow(:,:)      !  Number of snow layers on ground on tiles

  REAL, ALLOCATABLE ::                                                  &
    sice(:,:,:),                                                        &
                    ! Snow layer ice mass on tiles (Kg/m2)
    sliq(:,:,:),                                                        &
                    ! Snow layer liquid mass on tiles (Kg/m2)
    snowdepth(:,:),                                                     &
                    ! Snow depth on ground on tiles (m)
    tsnow(:,:,:),                                                       &
                    ! Snow layer temperature (K)
    rgrainl(:,:,:),                                                     &
                    ! Snow layer grain size on tiles (microns)
    rho_snow_grnd(:,:),                                                 &
                    !   Snowpack bulk density (kg/m3)
    rho_snow(:,:,:)
                    ! Snow layer densities (m)



!-----------------------------------------------------------------------
! If we are not in UM, define everything else that is needed
!-----------------------------------------------------------------------
#if !defined(UM_JULES)

  REAL, ALLOCATABLE ::                                                  &
    canht_ft(:,:),                                                      &
                !  Canopy height (m)
    canopy(:,:),                                                        &
                !  Surface/canopy water for snow-free land tiles (kg/m2)
    canopy_gb(:),                                                       &
                !  Gridbox canopy water content (kg/m2)
    cs(:,:),                                                            &
                !  Soil carbon (kg C/m2)
!               !  If dim_cs1=1, there is a single soil C pool.
!               !  If dim_cs1=4, the pools are:
!               !  1  decomposable plant material
!               !  2  resistant plant material
!               !  3  biomass
!               !  4  humus
    di(:,:),                                                            &
                !  "Equivalent thickness" of sea-ice GBM agregate (m)
    di_ncat(:,:,:),                                                     &
                !  "Equivalent thickness" of sea-ice catagories (m)
    k_sice(:,:,:),                                                      &
                !  Sea ice effective conductivity (2*kappai/de)
    gc(:,:),                                                            &
                !  Stomatal" conductance to evaporation for land tiles(m/s)
    gs(:),                                                              &
                !  "Stomatal" conductance to evaporation (m/s)
    lai(:,:),                                                           &
                !  LAI of plant functional types
    rgrain(:,:),                                                        &
                !  Snow surface grain size on tiles (microns)
    smc(:),                                                             &
                !  Soil moisture in a layer at the surface (kg/m2).
!      Note that SMC is used twice:
!      1) SF_EXPL and SF_IMPL2 use it to return the available water in the total
!         soil column (as calculated in PHYSIOL)
!      2) HYDROL uses it to return the available water in a layer of a given
!         depth (as calculated in SOILMC)
    smcl(:,:),                                                          &
                !  Soil moisture content of layers (kg/m2)
    snow_tile(:,:),                                                     &
                !  Lying snow on tiles (kg/m2)
                !  If can_model=4,
                !    snow_tile is the snow on the canopy
                !    snow_grnd is the snow on the ground beneath canopy
                !  If can_model/=4, snow_tile is the total snow.
    snow_grnd(:,:),                                                     &
                !  Snow on the ground (kg/m2)
                !  This is the snow beneath the canopy and is only
                !  used if can_model=4.
    snow_mass(:,:),                                                     &
                !  Gridbox snowmass (kg/m2)
    snow_mass_sea_ncat(:,:,:),                                          &
                !  Snow on category sea-ice (Kg/m2)
    soot(:),                                                            &
                !  Snow soot content (kg/kg)
    t_soil(:,:),                                                        &
                !  Sub-surface temperatures (K)
    ti(:,:),                                                            &
                !  Sea-ice surface layer
    tstar_tile(:,:),                                                    &
                !  Tile surface temperatures (K)
    z0msea(:,:),                                                        &
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

#endif

END MODULE prognostics
