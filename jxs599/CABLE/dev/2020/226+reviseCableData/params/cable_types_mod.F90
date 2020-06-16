MODULE cable_types_mod
#define UM_BUILD yes

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines variable types and variables for CABLE standalone runs.
!   Based on cable_def_types_mod.F90 from the CABLE trunk.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

!calculated length of cable arrays = number of active tiles
INTEGER :: mp          ! # total no of patches/tiles

LOGICAL, ALLOCATABLE :: l_tile_pts(:,:)

!jh!elevate these to namelist status
INTEGER ::        & 
#ifdef UM_BUILD 
  mvtype=17,      & ! total # vegetation types,   from input
  mstype=9,       & ! total # soil types, needs to de defined atCompile TimeForNow
#else       
  mvtype,         & ! total # vegetation types,   from input
  mstype,         & ! total # soil types,         from input
#endif
  mland                           ! # land grid cells

!jh!elevate these to namelist status
INTEGER, PARAMETER ::                                                        &
       i_d  = KIND(9), &
#ifdef UM_BUILD 
       r_2  = KIND(1.0),&!SELECTED_REAL_KIND(12, 50), &
#else       
       r_2  = KIND(1.d0),&!SELECTED_REAL_KIND(12, 50), &
#endif
       n_tiles = 17,  & ! # possible no of different
       ncp = 3,       & ! # vegetation carbon stores
       ncs = 2,       & ! # soil carbon stores
       mf = 2,        & ! # leaves (sunlit, shaded)
!jh!this needs sorting (nrb=4[2*2])
       nrb = 3,       & ! # radiation bands
       msn = 3,       & ! max # snow layers
       swb = 2,       & ! # shortwave bands
       niter = 4,     & ! number of iterations for za/L
                                !      ms = 12          ! # soil layers
       ms = 6         ! # soil layers - standard

  INTEGER, PARAMETER :: n_ktherm = 3

!!  ! Energy and water balance variables:
!!  TYPE balances_type
!!
!!     REAL, DIMENSION(:), POINTER ::                                           &
!!          drybal,           & ! energy balance for dry canopy
!!          ebal,             & ! energy balance per time step (W/m^2)
!!          ebal_tot,         & ! cumulative energy balance (W/m^2)
!!          ebal_cncheck,     & ! energy balance consistency check (W/m^2)
!!          ebal_tot_cncheck, & ! cumulative energy balance (W/m^2)
!!          ebaltr,           & ! energy balance per time step (W/m^2)
!!          ebal_tottr,       & ! cumulative energy balance (W/m^2)
!!          evap_tot,         & ! cumulative evapotranspiration (mm/dels)
!!          osnowd0,          & ! snow depth, first time step
!!          precip_tot,       & ! cumulative precipitation (mm/dels)
!!          rnoff_tot,        & ! cumulative runoff (mm/dels)
!!          wbal,             & ! water balance per time step (mm/dels)
!!          wbal_tot,         & ! cumulative water balance (mm/dels)
!!          wbtot0,           & ! total soil water (mm), first time step
!!          wetbal,           & ! energy balance for wet canopy
!!          cansto0,          & ! canopy water storage (mm)
!!          owbtot,           & ! total soil water (mm), first time step
!!          evapc_tot,        & ! cumulative evapotranspiration (mm/dels)
!!          evaps_tot,        & ! cumulative evapotranspiration (mm/dels)
!!          rnof1_tot,        & ! cumulative runoff (mm/dels)
!!          rnof2_tot,        & ! cumulative runoff (mm/dels)
!!          snowdc_tot,       & ! cumulative runoff (mm/dels)
!!          wbal_tot1,        & ! cumulative water balance (mm/dels)
!!          delwc_tot,        & ! energy balance for wet canopy
!!          qasrf_tot,        & ! heat advected to the snow by precip.
!!          qfsrf_tot,        & ! energy of snowpack phase changes
!!          qssrf_tot, &        ! energy of snowpack phase changes
!!          Radbal, &
!!          EbalSoil, &
!!          Ebalveg, &
!!          Radbalsum
!!
!!  END TYPE balances_type
!!
!!  ! .............................................................................
!!
!!  ! Roughness variables:
!!  TYPE roughness_type
!!
!!     REAL, DIMENSION(:), POINTER ::                                           &
!!          disp,    & ! zero-plane displacement
!!          hruff,   & ! canopy height above snow level
!!          hruff_grmx,&! max ht of canopy from tiles on same grid
!!          rt0us,   & ! eq. 3.54, SCAM manual (CSIRO tech report 132)
!!          rt1usa,  & ! resistance from disp to hruf
!!          rt1usb,  & ! resist fr hruf to zruffs (zref if zref<zruffs)
!!          rt1,     & ! 1/aerodynamic conductance
!!          za_uv,   & ! level of lowest atmospheric model layer
!!          za_tq,   & ! level of lowest atmospheric model layer
!!          z0m,     & ! roughness length
!!          zref_uv, & ! Reference height for met forcing
!!          zref_tq, & ! Reference height for met forcing
!!          zruffs,  & ! SCALAR Roughness sublayer depth (ground=origin)
!!          z0soilsn,& ! roughness length of bare soil surface
!!          z0soil     ! roughness length of bare soil surface
!!
!!     ! "coexp": coefficient in exponential in-canopy wind profile
!!     ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
!!     ! canopy and roughness-sublayer U(z) at z=h
!!     REAL, DIMENSION(:), POINTER ::                                           &
!!          coexp ! Extinction coef for wind profile in canopy
!!
!!     ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
!!     REAL, DIMENSION(:), POINTER ::                                           &
!!          usuh ! Friction velocity/windspeed at canopy height
!!
!!     REAL, DIMENSION(:), POINTER ::                                           &
!!          term2, term3, term5, term6, term6a ! for aerodyn resist. calc.
!!
!!
!!
!!  END TYPE roughness_type
!!
!!  ! .............................................................................
!!
!!  ! Air variables:
!!  TYPE air_type
!!
!!     REAL, DIMENSION(:), POINTER ::                                           &
!!          rho,     & ! dry air density (kg m-3)
!!          volm,    & ! molar volume (m3 mol-1)
!!          rlam,    & ! latent heat for water (j/kg)
!!          qsat,    & ! saturation specific humidity
!!          epsi,    & ! d(qsat)/dT ((kg/kg)/K)
!!          visc,    & ! air kinematic viscosity (m2/s)
!!          psyc,    & ! psychrometric constant
!!          dsatdk,  & ! d(es)/dT (mb/K)
!!          cmolar     ! conv. from m/s to mol/m2/s
!!
!!  END TYPE air_type
!!
!!  ! Climate data:
!!  TYPE climate_type
!!
!!     INTEGER :: nyear_average = 20
!!     INTEGER :: nday_average  = 31
!!     !      INTEGER, POINTER ::                                                  &
!!     INTEGER ::                                                  &
!!          nyears, & ! number of years in climate record
!!          doy ! day of year
!!
!!     INTEGER, DIMENSION(:), POINTER ::                                   &
!!          chilldays, &   ! length of chilling period (period with T<5deg)
!!          iveg, &        ! potential vegetation type based on climatic constraints
!!          biome
!!
!!     REAL, DIMENSION(:), POINTER ::                                           &
!!          dtemp,        & ! daily temperature
!!          dmoist,        & ! daily moisture availability
!!          mtemp,       & ! mean temperature over the last 31 days
!!          qtemp,       & ! mean temperature over the last 91 days
!!          mmoist,        & ! monthly moisture availability
!!          mtemp_min,   & ! minimum monthly temperature
!!          mtemp_max,   & ! maximum monhtly temperature
!!          qtemp_max,   & ! mean temperature of the warmest quarter (so far this year)
!!          qtemp_max_last_year,   & ! mean temperature of the warmest quarter (last calendar year)
!!          mtemp_min20,   & ! minimum monthly temperature, averaged over 20 y
!!          mtemp_max20,   & ! maximum monhtly temperature, averaged over 20 y
!!          atemp_mean,  & ! annual average temperature
!!          AGDD5,       &
!!          GDD5,        & ! growing degree day sum relative to 5deg base temperature
!!          AGDD0,        & !
!!          GDD0,        & ! growing degree day sum relative to 0deg base temperature
!!          alpha_PT,    & ! ratio of annual evap to annual PT evap
!!          evap_PT,    & ! annual PT evap [mm]
!!          aevap , &       ! annual evap [mm]
!!          alpha_PT20
!!
!!     REAL, DIMENSION(:,:), POINTER ::                                   &
!!          mtemp_min_20, & ! mimimum monthly temperatures for the last 20 y
!!          mtemp_max_20, & ! maximum monthly temperatures for the last 20 y
!!          dtemp_31 , &    ! daily temperature for the last 31 days
!!          dmoist_31 , &    ! daily moisture availability for the last 31 days
!!          alpha_PT_20, &      ! priestley Taylor Coefft for last 20 y
!!          dtemp_91     ! daily temperature for the last 91 days
!!
!!  END TYPE climate_type
!!
!!  ! .............................................................................
!!
!!  ! Cumulative flux variables:
!!  TYPE sum_flux_type
!!
!!     REAL, DIMENSION(:), POINTER ::                                           &
!!          sumpn,   & ! sum of canopy photosynthesis (g C m-2)
!!          sumrp,   & ! sum of plant respiration (g C m-2)
!!          sumrpw,  & ! sum of plant respiration (g C m-2)
!!          sumrpr,  & ! sum of plant respiration (g C m-2)
!!          sumrs,   & ! sum of soil respiration (g C m-2)
!!          sumrd,   & ! sum of daytime respiration (g C m-2)
!!          dsumpn,  & ! daily sumpn
!!          dsumrp,  & ! daily sumrp
!!          dsumrs,  & ! daily sumrs
!!          dsumrd,  & ! daily sumrd
!!          sumxrp,  & ! sum plant resp. modifier
!!          sumxrs     ! sum soil resp. modifier
!!
!!  END TYPE sum_flux_type
!!
!!  ! .............................................................................
!!
!!  TYPE bgc_pool_type
!!
!!     REAL, DIMENSION(:,:), POINTER ::                                         &
!!          cplant,  & ! plant carbon (g C/m2))
!!          csoil   ! soil carbon (g C/m2)
!!
!!
!!     REAL, DIMENSION(ncp)  :: ratecp ! plant carbon rate constant (1/year)
!!
!!     REAL, DIMENSION(ncs)  :: ratecs ! soil carbon rate constant (1/year)
!!
!!  END TYPE bgc_pool_type

END MODULE cable_types_mod
