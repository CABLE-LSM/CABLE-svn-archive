MODULE cable_Duplicate_types_mod
   
   use cable_dupVars_mod
   use cable_def_types_mod, only :  alloc, mp, mvtype, mstype, mland, & 
                                    r_2, n_tiles, ncp, ncs, mf, nrb, &
                                    msn, swb, niter, ms

   IMPLICIT NONE

   PUBLIC
  
! .............................................................................
   ! Soil parameters:
   TYPE Dsoil_parameter_type 
   
      INTEGER, DIMENSION(:,:), allocatable ::                                        &
         isoilm     ! integer soil type

      REAL, DIMENSION(:,:), allocatable ::                                           &
         bch,     & ! parameter b in Campbell equation
         c3,      & ! c3 drainage coeff (fraction)
         clay,    & ! fraction of soil which is clay
         css,     & ! soil specific heat capacity [kJ/kg/K]
         hsbh,    & ! difsat * etasat (=hyds*abs(sucs)*bch)
         hyds,    & ! hydraulic conductivity @ saturation [m/s], Ksat
         i2bp3,   & ! par. one in K vis suction (=nint(bch)+2)
         ibp2,    & ! par. two in K vis suction (fn of pbch)
         rhosoil, & ! soil density [kg/m3]
         sand,    & ! fraction of soil which is sand
         sfc,     & ! vol H2O @ field capacity
         silt,    & ! fraction of soil which is silt
         ssat,    & ! vol H2O @ saturation
         sucs,    & ! suction at saturation (m)
         swilt,   & ! vol H2O @ wilting
         zse,     & ! thickness of each soil layer (1=top) in m
         zshh,    & ! distance between consecutive layer midpoints (m)
         albsoilf   ! soil reflectance
     
      REAL(r_2), DIMENSION(:,:), allocatable ::                                      &
         cnsd,    & ! thermal conductivity of dry soil [W/m/K]
         pwb_min    ! working variable (swilt/ssat)**ibp2
     
      REAL, DIMENSION(:,:,:), allocatable ::                                         &
         albsoil    ! soil reflectance (2nd dim. BP 21Oct2009)

  END TYPE Dsoil_parameter_type

! .............................................................................

   ! Soil and snow variables:
   TYPE Dsoil_snow_type 
     
     INTEGER, DIMENSION(:,:), allocatable :: isflag ! 0 => no snow 1 => snow
    
      REAL, DIMENSION(:,:), allocatable ::                                           &
         iantrct, & ! allocatable to Antarctic land points
         pudsto,  & ! puddle storage
         pudsmx,  & ! puddle storage
         cls,     & ! factor for latent heat
         dfn_dtg, & ! d(canopy%fns)/d(ssnow%tgg)
         dfh_dtg, & ! d(canopy%fhs)/d(ssnow%tgg)
         dfe_ddq, & ! d(canopy%fes)/d(dq)
         ddq_dtg, & ! d(dq)/d(ssnow%tgg)
         evapsn,  & ! snow evaporation  
         fwtop,   & ! water flux to the soil
         fwtop1,  & ! water flux to the soil
         fwtop2,  & ! water flux to the soil
         fwtop3,  & ! water flux to the soil
         osnowd,  & ! snow depth from previous time step
         potev,   & ! potential evapotranspiration
         runoff,  & ! total runoff (mm/dels)
         rnof1,   & ! surface runoff (mm/dels)
         rnof2,   & ! deep drainage (mm/dels)
         rtsoil,  & ! turbulent resistance for soil
         wbtot1,  & ! total soil water (mm)
         wbtot2,  & ! total soil water (mm)
         wb_lake, &
         sinfil,  & 
         qstss,   & 
         wetfac,  & ! surface wetness fact. at current time step
         owetfac, & ! surface wetness fact. at previous time step
         t_snwlr, & ! top snow layer depth in 3 layer snowpack
         tggav,   & ! mean soil temperature in K
         otgg,    & ! soil temperature in K
         otss,    & ! surface temperature (weighted soil, snow)
         otss_0,  & ! surface temperature (weighted soil, snow)
         tprecip, &
         tevap,   &
         trnoff,  &
         totenbal,&!
         totenbal2,&
         fland,   & ! factor for latent heat
         ifland,  & ! integer soil type
         qasrf,   & ! heat advected to the snow by precip. 
         qfsrf,   & ! energy of snowpack phase changes 
         qssrf,   & ! sublimation 
         snage,   & ! snow age
         snowd,   & ! snow depth (liquid water)
         smelt,   & ! snow melt 
         ssdnn,   & ! average snow density
         tss,     & ! surface temperature (weighted soil, snow)
         tss_p,   & ! surface temperature (weighted soil, snow)
         deltss,  & ! surface temperature (weighted soil, snow)
         owb1       ! surface temperature (weighted soil, snow)
 
      REAL, DIMENSION(:,:,:), allocatable ::                                         &
         sconds,     & !
         sdepth,     & ! snow depth
         smass,      & ! snow mass
         ssdn,       & ! snow densities
         tgg,        & ! soil temperature in K
         tggsn,      & ! snow temperature in K
         dtmlt,      & ! water flux to the soil
         albsoilsn,  & ! soil + snow reflectance
         evapfbl,    & !
         tilefrac      ! factor for latent heat
     
    
      REAL(r_2), DIMENSION(:,:), allocatable ::                                      &
         wbtot   ! total soil water (mm)
     
      REAL(r_2), DIMENSION(:,:,:), allocatable ::                                    &
         gammzz,  & ! heat capacity for each soil layer
         wb,      & ! volumetric soil moisture (solid+liq)
         wbice,   & ! soil ice
         wblf,    & !
         wbfice     !

   END TYPE Dsoil_snow_type

! .............................................................................

   ! Vegetation parameters:
   TYPE Dveg_parameter_type
     
      INTEGER, DIMENSION(:,:), allocatable ::                                        &
         iveg       ! vegetation type

      REAL, DIMENSION(:,:), allocatable ::                                           &
         canst1,  & ! max intercepted water by canopy (mm/LAI)
         dleaf,   & ! chararacteristc legnth of leaf (m)
         ejmax,   & ! max pot. electron transp rate top leaf(mol/m2/s)
         meth,    & ! method for calculation of canopy fluxes and temp.
         frac4,   & ! fraction of c4 plants
         hc,      & ! roughness height of canopy (veg - snow)
         vlai,    & ! leaf area index
         xalbnir, & 
         rp20,    & ! plant respiration coefficient at 20 C
         rpcoef,  & ! temperature coef nonleaf plant respiration (1/C)
         rs20,    & ! soil respiration at 20 C [mol m-2 s-1]
         shelrb,  & ! sheltering factor (dimensionless)
         vegcf,   & ! kdcorbin, 08/10
         tminvj,  & ! min temperature of the start of photosynthesis
         tmaxvj,  & ! max temperature of the start of photosynthesis
         vbeta,   & ! 
         vcmax,   & ! max RuBP carboxylation rate top leaf (mol/m2/s)
         xfang,   & ! leaf angle PARAMETER
         extkn,   & ! extinction coef for vertical
         vlaimax, & ! extinction coef for vertical
         wai        ! wood area index (stem+branches+twigs)

      LOGICAL, DIMENSION(:,:), allocatable ::                                        &
         deciduous ! flag used for phenology fix

      REAL, DIMENSION(:,:,:), allocatable ::                                         &
         refl,    &
         taul,    & 
         froot      ! fraction of root in each soil layer

   END TYPE Dveg_parameter_type

! .............................................................................

   ! Canopy/vegetation variables:
   TYPE Dcanopy_type
      

      REAL, DIMENSION(:,:), allocatable ::                                           &
         cansto,  & ! canopy water storage (mm)
         cduv,    & ! drag coefficient for momentum
         delwc,   & ! change in canopy water store (mm/dels)
         dewmm,   & ! dewfall (mm)
         fe,      & ! total latent heat (W/m2)
         fh,      & ! total sensible heat (W/m2)
         fpn,     & ! plant photosynthesis (g C m-2 s-1)
         frp,     & ! plant respiration (g C m-2 s-1)
         frpw,    & ! plant respiration (g C m-2 s-1)???
         frpr,    & ! plant respiration (g C m-2 s-1)???
         frs,     & ! soil respiration (g C m-2 s-1)
         fnee,    & ! net carbon flux (g C m-2 s-1)
         frday,   & ! daytime leaf resp
         fnv,     & ! net rad. avail. to canopy (W/m2)
         fev,     & ! latent hf from canopy (W/m2)
         epot,    & ! total potential evaporation 
         fnpp,    & ! npp flux
         fevw_pot,& ! potential lat heat from canopy
         gswx_T,  & ! ! stom cond for water
         cdtq,    & ! drag coefficient for momentum
         wetfac_cs,&! 
         fevw,    & ! lat heat fl wet canopy (W/m2)
         fhvw,    & ! sens heatfl from wet canopy (W/m2)
         oldcansto,&! canopy water storage (mm)
         fhv,     & ! sens heatfl from canopy (W/m2)
         fns,     & ! net rad avail to soil (W/m2)
         fhs,     & ! sensible heat flux from soil
         fhs_cor, &
         ga,      & ! ground heat flux (W/m2) ???
         ghflux,  & ! ground heat flux (W/m2) ???
         precis,  & ! throughfall to soil, after snow (mm)
         qscrn,   & ! specific humudity at screen height (g/g)
         rnet,    & ! net radiation absorbed by surface (W/m2)
         segg,    & ! latent heatfl from soil mm
         sghflux, & ! ground heat flux (W/m2) ???
         through, & ! canopy throughfall (mm)
         spill,   & ! can.storage excess after dewfall (mm)
         tscrn,   & ! air temperature at screen height (oC)
         wcint,   & ! canopy rainfall interception (mm)
         tv,      & ! vegetation temp (K)
         us,      & ! friction velocity
         uscrn,   & ! wind speed at screen height (m/s)
         vlaiw,   & ! lai adj for snow depth for calc of resistances
         rghlai,  & ! lai adj for snow depth for calc of resistances
         fwet       ! fraction of canopy wet

      REAL, DIMENSION(:,:,:), allocatable ::                                         &
         evapfbl, &
         gswx,    & ! stom cond for water
         zetar      ! stability correction

      REAL(r_2), DIMENSION(:,:), allocatable ::                                      &
         fess,    & ! latent heatfl from soil (W/m2)
         fesp,    & ! latent heatfl from soil (W/m2)
         dgdtg,   & ! derivative of gflux wrt soil temp
         fes,     & ! latent heatfl from soil (W/m2)
         fes_cor, & ! latent heatfl from soil (W/m2)
         fevc       ! dry canopy transpiration (W/m2)

   END TYPE Dcanopy_type

! .............................................................................

   ! Meterological data:
   TYPE Dmet_type
     
      INTEGER, DIMENSION(:,:), allocatable ::                                        &
         year,    & ! local time year AD 
         moy        ! local time month of year 
     
      REAL, DIMENSION(:,:), allocatable ::                                           &
         ca,      & ! CO2 concentration (mol/mol)
         doy,     & ! local time day of year = days since 0 hr 1st Jan 
         hod,     & ! local hour of day
         ofsd,    & ! downward short-wave radiation (W/m2)
         fld,     & ! downward long-wave radiation (W/m2)
         precip,  & ! rainfall (liquid+solid)(mm/dels)
         precip_sn,&! solid preipitation only (mm/dels)
         tk,      & ! surface air temperature (oK)
         tvair,   & ! within canopy air temperature (oK)
         tvrad,   & ! radiative vegetation temperature (K)
         pmb,     & ! surface air pressure (mbar)
         ua,      & ! surface wind speed (m/s)
         qv,      & ! surface specific humidity (g/g)
         qvair,   & ! within canopy specific humidity (g/g)
         da,      & ! water vap pressure deficit at ref height (Pa)
         dva,     & ! in canopy water vap pressure deficit (Pa)
         coszen     ! cos(zenith angle of sun)
     
      REAL, DIMENSION(:,:,:), allocatable ::                                         &
         fsd  ! downward short-wave radiation (W/m2)
     
   END TYPE Dmet_type

! .............................................................................

   TYPE Dbgc_pool_type
      
      REAL, DIMENSION(:,:,:), allocatable ::                                         &
         cplant,  & ! plant carbon (g C/m2))
         csoil      ! soil carbon (g C/m2)
     
      REAL, DIMENSION(ncp)  :: ratecp ! plant carbon rate constant (1/year)
      
      REAL, DIMENSION(ncs)  :: ratecs ! soil carbon rate constant (1/year)
   
   END TYPE Dbgc_pool_type

! .............................................................................

   ! Radiation variables:
   TYPE Dradiation_type
   
      REAL, DIMENSION(:,:), allocatable   ::                                         &
         transb,  & ! fraction SW beam tranmitted through canopy
         albedo_T,& ! canopy+soil albedo for VIS+NIR
         longitude,&! longitude
         workp1,  & ! absorbed short-wave radiation for soil
         workp2,  & ! absorbed short-wave radiation for soil
         workp3,  & ! absorbed short-wave radiation for soil
         extkb,   & ! beam radiation extinction coeff
         extkd2,  & ! diffuse 2D radiation extinction coeff
         extkd,   & ! diffuse radiation extinction coeff (-)
         flws,    & ! soil long-wave radiation
         latitude,& ! latitude
         lwabv,   & ! long wave absorbed by vegetation
         qssabs,  & ! absorbed short-wave radiation for soil
         transd,  & ! frac SW diffuse transmitted through canopy
         trad       !  radiative temperature (soil and veg)
     
      REAL, DIMENSION(:,:,:), allocatable  ::                                        &
         fvlai,   & ! leaf area index of big leaf
         rhocdf,  & ! canopy diffuse reflectance (-)
         rniso,   & ! sum(rad%qcan, 3) total abs by canopy (W/m2)
         scalex,  & ! scaling PARAMETER for big leaf
         albedo,  & ! canopy+soil albedo
         reffdf,  & ! effective conopy diffuse reflectance
         reffbm,  & ! effective conopy beam reflectance
         extkbm,  & ! modified k beam(6.20)(for leaf scattering)
         extkdm,  & ! modified k diffuse(6.20)(for leaf scattering)
         fbeam,   & ! beam fraction 
         cexpkbm, & ! canopy beam transmittance
         cexpkdm, & ! canopy diffuse transmittance
         rhocbm,  & ! modified canopy beam reflectance(6.21)
         gradis     ! radiative conductance
     
      REAL, DIMENSION(:,:,:,:), allocatable ::                                       &
         qcan ! absorbed radiation for canopy (W/m^2)
    
    
  END TYPE Dradiation_type

 !.............................................................................

   ! Roughness variables:
   TYPE Droughness_type
      
      REAL, DIMENSION(:,:), allocatable ::                                           &
         disp,    & ! zero-plane displacement
         hruff,   & ! canopy height above snow level
         hruff_grmx,&! max ht of canopy from tiles on same grid 
         rt0us,   & ! eq. 3.54, SCAM manual (CSIRO tech report 132)
         rt1usa,  & ! resistance from disp to hruf
         rt1usb,  & ! resist fr hruf to zruffs (zref if zref<zruffs)
         rt1,     & ! 1/aerodynamic conductance
         za_uv,   & ! level of lowest atmospheric model layer
         za_tq,   & ! level of lowest atmospheric model layer
         z0m,     & ! roughness length
         zref_uv, & ! Reference height for met forcing
         zref_tq, & ! Reference height for met forcing
         zruffs,  & ! SCALAR Roughness sublayer depth (ground=origin)
         z0soilsn,& ! roughness length of bare soil surface
         z0soil     ! roughness length of bare soil surface
      
      ! "coexp": coefficient in exponential in-canopy wind profile
      ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
      ! canopy and roughness-sublayer U(z) at z=h
      REAL, DIMENSION(:,:), allocatable ::                                           &
         coexp ! Extinction coef for wind profile in canopy
     
      ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
      REAL, DIMENSION(:,:), allocatable ::                                           &
         usuh ! Friction velocity/windspeed at canopy height
   
      REAL, DIMENSION(:,:), allocatable ::                                           &
         term2, term3, term5, term6 ! for aerodyn resist. calc.
   
   END TYPE Droughness_type

! .............................................................................

   ! Air variables:
   TYPE Dair_type
      
      REAL, DIMENSION(:,:), allocatable ::                                           &
         rho,     & ! dry air density (kg m-3)
         volm,    & ! molar volume (m3 mol-1)
         rlam,    & ! latent heat for water (j/kg)
         qsat,    & ! saturation specific humidity
         epsi,    & ! d(qsat)/dT ((kg/kg)/K)
         visc,    & ! air kinematic viscosity (m2/s)
         psyc,    & ! psychrometric constant
         dsatdk,  & ! d(es)/dT (mb/K)
         cmolar     ! conv. from m/s to mol/m2/s

   END TYPE Dair_type

! .............................................................................

   type TdupVars
      ! CABLE variables
      TYPE (Dmet_type)       :: met     ! met input variables
      TYPE (Dair_type)       :: air     ! air property variables
      TYPE (Dcanopy_type)    :: canopy  ! vegetation variables
      TYPE (Dradiation_type) :: rad     ! radiation variables
      TYPE (Droughness_type) :: rough   ! roughness varibles
      TYPE (Dsoil_snow_type) :: ssnow   ! soil and snow variables
      !TYPE (Dbalances_type)  :: bal     ! energy and water balance variables
      
      ! CABLE parameters
      TYPE (Dsoil_parameter_type) :: soil ! soil parameters	
      TYPE (Dveg_parameter_type)  :: veg  ! vegetation parameters	 
      !TYPE (Ddriver_type)    :: C         ! constants used locally  
      
      !TYPE (Dsum_flux_type)  :: sum_flux ! cumulative flux variables
      TYPE (Dbgc_pool_type)  :: bgc  ! carbon pool variables
      
      !! CASA-CNP variables 
      !TYPE (Dcasa_biome)     :: casabiome
      !TYPE (Dcasa_pool)      :: casapool
      !TYPE (Dcasa_flux)      :: casaflux
      !TYPE (Dcasa_met)       :: casamet
      !TYPE (Dcasa_balance)   :: casabal
      !TYPE (Dphen_variable)  :: phen 
 
   end type TdupVars

! .............................................................................


!   ! Energy and water balance variables:
!   TYPE balances_type 
!
!      REAL, DIMENSION(:,:), allocatable ::                                           &
!         drybal,           & ! energy balance for dry canopy
!         ebal,             & ! energy balance per time step (W/m^2)
!         ebal_tot,         & ! cumulative energy balance (W/m^2)
!         ebal_cncheck,     & ! energy balance consistency check (W/m^2)
!         ebal_tot_cncheck, & ! cumulative energy balance (W/m^2)
!         ebaltr,           & ! energy balance per time step (W/m^2)
!         ebal_tottr,       & ! cumulative energy balance (W/m^2)
!         evap_tot,         & ! cumulative evapotranspiration (mm/dels)
!         osnowd0,          & ! snow depth, first time step
!         precip_tot,       & ! cumulative precipitation (mm/dels)
!         rnoff_tot,        & ! cumulative runoff (mm/dels)
!         wbal,             & ! water balance per time step (mm/dels)
!         wbal_tot,         & ! cumulative water balance (mm/dels)
!         wbtot0,           & ! total soil water (mm), first time step
!         wetbal,           & ! energy balance for wet canopy
!         cansto0,          & ! canopy water storage (mm)
!         owbtot,           & ! total soil water (mm), first time step
!         evapc_tot,        & ! cumulative evapotranspiration (mm/dels)
!         evaps_tot,        & ! cumulative evapotranspiration (mm/dels)
!         rnof1_tot,        & ! cumulative runoff (mm/dels)
!         rnof2_tot,        & ! cumulative runoff (mm/dels)
!         snowdc_tot,       & ! cumulative runoff (mm/dels)
!         wbal_tot1,        & ! cumulative water balance (mm/dels)
!         delwc_tot,        & ! energy balance for wet canopy
!         qasrf_tot,        & ! heat advected to the snow by precip. 
!         qfsrf_tot,        & ! energy of snowpack phase changes 
!         qssrf_tot           ! energy of snowpack phase changes 
!
!   END TYPE balances_type
!
!! .............................................................................


!   ! Cumulative flux variables:
!   TYPE sum_flux_type
!     
!      REAL, DIMENSION(:), allocatable ::                                           &
!         sumpn,   & ! sum of canopy photosynthesis (g C m-2)
!         sumrp,   & ! sum of plant respiration (g C m-2)
!         sumrpw,  & ! sum of plant respiration (g C m-2)
!         sumrpr,  & ! sum of plant respiration (g C m-2)
!         sumrs,   & ! sum of soil respiration (g C m-2)
!         sumrd,   & ! sum of daytime respiration (g C m-2)
!         dsumpn,  & ! daily sumpn
!         dsumrp,  & ! daily sumrp
!         dsumrs,  & ! daily sumrs
!         dsumrd,  & ! daily sumrd
!         sumxrp,  & ! sum plant resp. modifier
!         sumxrs     ! sum soil resp. modifier
!
!   END TYPE sum_flux_type

! .............................................................................

   ! Functions for allocating these types
   ! All overloaded so code only needs to call alloc_cbm_var
   ! Alloc routines could all initialise to NaN or zero for debugging?
   ! Don't need the mp argument here as it's a module variable.
!   PUBLIC :: alloc_cbm_var
!   PRIVATE :: alloc_bgc_pool_type, dealloc_bgc_pool_type
   
!   INTERFACE alloc_cbm_var
!      MODULE PROCEDURE alloc_balances_type,                                    &
!         alloc_soil_parameter_type,                                            &
!         alloc_soil_snow_type,                                                 &
!         alloc_veg_parameter_type,                                             &
!         alloc_canopy_type,                                                    &
!         alloc_radiation_type,                                                 &
!         alloc_roughness_type,                                                 &
!         alloc_air_type,                                                       &
!         alloc_met_type,                                                       &
!         alloc_sum_flux_type,                                                  &
!         alloc_bgc_pool_type            
!   END INTERFACE
!
!   INTERFACE dealloc_cbm_var
!      MODULE PROCEDURE dealloc_balances_type,                                  &
!         dealloc_soil_parameter_type,                                          &
!         dealloc_soil_snow_type,                                               &
!         dealloc_veg_parameter_type,                                           &
!         dealloc_canopy_type,                                                  &
!         dealloc_radiation_type,                                               &
!         dealloc_roughness_type,                                               &
!         dealloc_air_type,                                                     &
!         dealloc_met_type,                                                     &
!         dealloc_sum_flux_type,                                                &
!         dealloc_bgc_pool_type            
!   END INTERFACE


CONTAINS

  
SUBROUTINE alloc_DupVars( dup, kend )
   integer :: kend
   type (TdupVars) :: dup

   call cable_dupVars( dup% bgc% cplant, mp, ncp, kend, alloc )
   call cable_dupVars( dup% bgc% csoil, mp, ncs, kend, alloc )

   call cable_dupVars( dup% soil% bch,    mp, kend, alloc )   
   call cable_dupVars( dup% soil% c3,     mp, kend, alloc )    
   call cable_dupVars( dup% soil% clay,   mp, kend, alloc )  
   call cable_dupVars( dup% soil% css,    mp, kend, alloc )   
   call cable_dupVars( dup% soil% hsbh,   mp, kend, alloc )  
   call cable_dupVars( dup% soil% hyds,   mp, kend, alloc )  
   call cable_dupVars( dup% soil% i2bp3,  mp, kend, alloc ) 
   call cable_dupVars( dup% soil% ibp2,   mp, kend, alloc )  
   call cable_dupVars( dup% soil% isoilm, mp, kend, alloc )  
   call cable_dupVars( dup% soil% rhosoil,mp, kend, alloc )  
   call cable_dupVars( dup% soil% sand,   mp, kend, alloc )   
   call cable_dupVars( dup% soil% sfc,    mp, kend, alloc )   
   call cable_dupVars( dup% soil% silt,   mp, kend, alloc )   
   call cable_dupVars( dup% soil% ssat,   mp, kend, alloc )   
   call cable_dupVars( dup% soil% sucs,   mp, kend, alloc )   
   call cable_dupVars( dup% soil% swilt,  mp, kend, alloc )  
   call cable_dupVars( dup% soil% zse,    ms, kend, alloc )    
   call cable_dupVars( dup% soil% zshh, ms+1, kend, alloc )  
   call cable_dupVars( dup% soil% cnsd,   mp, kend, alloc )  
   
   call cable_dupVars( dup% soil% albsoil,   mp, nrb, kend, alloc)   
   call cable_dupVars( dup% soil% pwb_min,   mp, kend, alloc )  
   call cable_dupVars( dup% soil% albsoilf,  mp, kend, alloc )  
   
   call cable_DupVars( dup% ssnow% iantrct,   mp, kend, alloc )
   call cable_DupVars( dup% ssnow % pudsto,   mp, kend, alloc )
   call cable_DupVars( dup% ssnow % pudsmx,   mp, kend, alloc )
   call cable_DupVars( dup% ssnow% cls,       mp, kend, alloc )     
   call cable_DupVars( dup% ssnow% dfn_dtg,   mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow% dfh_dtg,   mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow% dfe_ddq,   mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow% ddq_dtg,   mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow% evapsn,    mp, kend, alloc )  
   call cable_DupVars( dup% ssnow% fwtop,     mp, kend, alloc )   
   call cable_DupVars( dup% ssnow% fwtop1,    mp, kend, alloc )   
   call cable_DupVars( dup% ssnow% fwtop2,    mp, kend, alloc )   
   call cable_DupVars( dup% ssnow% fwtop3,    mp, kend, alloc )   
   call cable_DupVars( dup% ssnow% gammzz,    mp, ms, kend, alloc ) 
   call cable_DupVars( dup% ssnow% isflag,    mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow% osnowd,    mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow% potev,     mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow% runoff,    mp, kend, alloc )
   call cable_DupVars( dup% ssnow% rnof1,     mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow% rnof2,     mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow% rtsoil,    mp, kend, alloc )
   call cable_DupVars( dup% ssnow% sconds,    mp, msn, kend, alloc ) 
   call cable_DupVars( dup% ssnow% sdepth,    mp, msn, kend, alloc ) 
   call cable_DupVars( dup% ssnow% smass,     mp, msn, kend, alloc ) 
   call cable_DupVars( dup% ssnow% snage,     mp, kend, alloc )  
   call cable_DupVars( dup% ssnow% snowd,     mp, kend, alloc )  
   call cable_DupVars( dup% ssnow% smelt,     mp, kend, alloc )  
   call cable_DupVars( dup% ssnow% ssdn,      mp, msn, kend, alloc )  
   call cable_DupVars( dup% ssnow% ssdnn,     mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow% tgg,       mp, ms, kend, alloc )   
   call cable_DupVars( dup% ssnow% tggsn,     mp, msn, kend, alloc )  
   call cable_DupVars( dup% ssnow% tss,       mp, kend, alloc )   
   call cable_DupVars( dup% ssnow% tss_p,     mp, kend, alloc )   
   call cable_DupVars( dup% ssnow% deltss,    mp, kend, alloc )   
   call cable_DupVars( dup% ssnow% owb1,      mp, kend, alloc )   
   call cable_DupVars( dup% ssnow% wb,        mp, ms, kend, alloc )    
   call cable_DupVars( dup% ssnow% wbice,     mp, ms, kend, alloc ) 
   call cable_DupVars( dup% ssnow% wblf,      mp, ms, kend, alloc ) 
   call cable_DupVars( dup% ssnow%wbtot,      mp, kend, alloc )    
   call cable_DupVars( dup% ssnow%wbtot1,     mp, kend, alloc )    
   call cable_DupVars( dup% ssnow%wbtot2,     mp, kend, alloc )    
   call cable_DupVars( dup% ssnow%wb_lake,    mp, kend, alloc )    
   call cable_DupVars( dup% ssnow%sinfil,     mp, kend, alloc )    
   call cable_DupVars( dup% ssnow%evapfbl,    mp, ms, kend, alloc )    
   call cable_DupVars( dup% ssnow%qstss,      mp, kend, alloc )    
   call cable_DupVars( dup% ssnow%wetfac,     mp, kend, alloc )  
   call cable_DupVars( dup% ssnow%owetfac,    mp, kend, alloc )  
   call cable_DupVars( dup% ssnow%t_snwlr,    mp, kend, alloc )  
   call cable_DupVars( dup% ssnow%wbfice,     mp, ms, kend, alloc )  
   call cable_DupVars( dup% ssnow%tggav,      mp, kend, alloc )  
   call cable_DupVars( dup% ssnow%otgg,       mp, kend, alloc )   
   call cable_DupVars( dup% ssnow%otss,       mp, kend, alloc )   
   call cable_DupVars( dup% ssnow%otss_0,     mp, kend, alloc )   
   call cable_DupVars( dup% ssnow%tprecip,    mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow%tevap,      mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow%trnoff,     mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow%totenbal,   mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow%totenbal2,  mp, kend, alloc ) 
   call cable_DupVars( dup% ssnow%fland,      mp, kend, alloc )      
   call cable_DupVars( dup% ssnow%ifland,     mp, kend, alloc )  
   call cable_DupVars( dup% ssnow%tilefrac,   mp, n_tiles, kend, alloc ) 
   call cable_DupVars( dup% ssnow%qasrf,      mp, kend, alloc )  
   call cable_DupVars( dup% ssnow%qfsrf,      mp, kend, alloc )  
   call cable_DupVars( dup% ssnow%qssrf,      mp, kend, alloc )  

   call cable_DupVars( dup% ssnow % dtmlt,   mp, 3, kend, alloc )
   call cable_DupVars( dup% ssnow% albsoilsn,mp, nrb, kend, alloc ) 

   call cable_DupVars( dup% veg% canst1,   mp, kend, alloc ) 
   call cable_DupVars( dup% veg% dleaf,    mp, kend, alloc )  
   call cable_DupVars( dup% veg% ejmax,    mp, kend, alloc ) 
   call cable_DupVars( dup% veg% iveg,     mp, kend, alloc ) 
   call cable_DupVars( dup% veg% meth,     mp, kend, alloc ) 
   call cable_DupVars( dup% veg% frac4,    mp, kend, alloc )  
   call cable_DupVars( dup% veg% hc,       mp, kend, alloc )     
   call cable_DupVars( dup% veg% vlai,     mp, kend, alloc )   
   call cable_DupVars( dup% veg% xalbnir,  mp, kend, alloc ) 
   call cable_DupVars( dup% veg% rp20,     mp, kend, alloc )   
   call cable_DupVars( dup% veg% rpcoef,   mp, kend, alloc ) 
   call cable_DupVars( dup% veg% rs20,     mp, kend, alloc )   
   call cable_DupVars( dup% veg% shelrb,   mp, kend, alloc ) 
   call cable_DupVars( dup% veg% vegcf,    mp, kend, alloc )  
   call cable_DupVars( dup% veg% tminvj,   mp, kend, alloc ) 
   call cable_DupVars( dup% veg% tmaxvj,   mp, kend, alloc ) 
   call cable_DupVars( dup% veg% vbeta,    mp, kend, alloc )  
   call cable_DupVars( dup% veg% vcmax,    mp, kend, alloc )  
   call cable_DupVars( dup% veg% xfang,    mp, kend, alloc )  
   call cable_DupVars( dup% veg% extkn,     mp, kend, alloc ) 
   call cable_DupVars( dup% veg% wai,       mp, kend, alloc )   
   call cable_DupVars( dup% veg% deciduous, mp, kend, alloc ) 
   call cable_DupVars( dup% veg% froot,     mp, ms, kend, alloc ) 
   call cable_DupVars( dup% veg% refl,      mp, 2, kend, alloc ) !jhan:swb?
   call cable_DupVars( dup% veg% taul,      mp, 2, kend, alloc ) 
   call cable_DupVars( dup% veg% vlaimax,   mp, kend, alloc ) 

   call cable_DupVars( dup% canopy% fess,    mp, kend, alloc )
   call cable_DupVars( dup% canopy% fesp,    mp, kend, alloc )
   call cable_DupVars( dup% canopy% cansto,  mp, kend, alloc )  
   call cable_DupVars( dup% canopy% cduv,    mp, kend, alloc )   
   call cable_DupVars( dup% canopy% delwc,   mp, kend, alloc )  
   call cable_DupVars( dup% canopy% dewmm,   mp, kend, alloc )  
   call cable_DupVars( dup% canopy% dgdtg,   mp, kend, alloc )  
   call cable_DupVars( dup% canopy% fe,      mp, kend, alloc )      
   call cable_DupVars( dup% canopy% fh,      mp, kend, alloc )      
   call cable_DupVars( dup% canopy% fpn,     mp, kend, alloc )     
   call cable_DupVars( dup% canopy% frp,     mp, kend, alloc )     
   call cable_DupVars( dup% canopy% frpw,    mp, kend, alloc )    
   call cable_DupVars( dup% canopy% frpr,    mp, kend, alloc )    
   call cable_DupVars( dup% canopy% frs,     mp, kend, alloc )     
   call cable_DupVars( dup% canopy% fnee,    mp, kend, alloc )    
   call cable_DupVars( dup% canopy% frday,   mp, kend, alloc )   
   call cable_DupVars( dup% canopy% fnv,     mp, kend, alloc )     
   call cable_DupVars( dup% canopy% fev,     mp, kend, alloc )     
   call cable_DupVars( dup% canopy% fevc,    mp, kend, alloc )    
   call cable_DupVars( dup% canopy% fhv,     mp, kend, alloc )     
   call cable_DupVars( dup% canopy% fns,     mp, kend, alloc )     
   call cable_DupVars( dup% canopy% fhs,     mp, kend, alloc )     
   call cable_DupVars( dup% canopy% fhs_cor, mp, kend, alloc )     
   call cable_DupVars( dup% canopy% ga,      mp, kend, alloc )      
   call cable_DupVars( dup% canopy% ghflux,  mp, kend, alloc )   
   call cable_DupVars( dup% canopy% precis,  mp, kend, alloc ) 
   call cable_DupVars( dup% canopy% qscrn,   mp, kend, alloc )  
   call cable_DupVars( dup% canopy% rnet,    mp, kend, alloc )   
   call cable_DupVars( dup% canopy% segg,    mp, kend, alloc )   
   call cable_DupVars( dup% canopy% sghflux, mp, kend, alloc )  
   call cable_DupVars( dup% canopy% through, mp, kend, alloc )  
   call cable_DupVars( dup% canopy% spill,   mp, kend, alloc )  
   call cable_DupVars( dup% canopy% tscrn,   mp, kend, alloc )  
   call cable_DupVars( dup% canopy% wcint,   mp, kend, alloc )  
   call cable_DupVars( dup% canopy% tv,      mp, kend, alloc )      
   call cable_DupVars( dup% canopy% us,      mp, kend, alloc )      
   call cable_DupVars( dup% canopy% uscrn,   mp, kend, alloc )   
   call cable_DupVars( dup% canopy% rghlai,  mp, kend, alloc ) 
   call cable_DupVars( dup% canopy% vlaiw,   mp, kend, alloc ) 
   call cable_DupVars( dup% canopy% fwet,    mp, kend, alloc )   
   call cable_DupVars( dup% canopy% epot,    mp, kend, alloc )   
   call cable_DupVars( dup% canopy% fnpp,    mp, kend, alloc )   
   call cable_DupVars( dup% canopy% fevw_pot,mp, kend, alloc )  
   call cable_DupVars( dup% canopy% gswx_T,  mp, kend, alloc )  
   call cable_DupVars( dup% canopy% cdtq,    mp, kend, alloc )   
   call cable_DupVars( dup% canopy% fevw,    mp, kend, alloc )   
   call cable_DupVars( dup% canopy% fhvw,    mp, kend, alloc )   
   call cable_DupVars( dup% canopy% fes,     mp, kend, alloc )    
   call cable_DupVars( dup% canopy% fes_cor, mp, kend, alloc )    
   call cable_DupVars( dup% canopy% gswx,    mp, mf, kend, alloc )  
   call cable_DupVars( dup% canopy% evapfbl, mp, ms, kend, alloc )
   call cable_DupVars( dup% canopy% wetfac_cs,  mp, kend, alloc )  
   call cable_DupVars( dup% canopy% oldcansto,  mp, kend, alloc )  
   call cable_DupVars( dup% canopy% zetar,   mp, NITER, kend, alloc )  
 
   call cable_DupVars ( dup% met% ca,     mp, kend, alloc )
   call cable_DupVars ( dup% met% year,   mp, kend, alloc )
   call cable_DupVars ( dup% met% moy,    mp, kend, alloc )
   call cable_DupVars ( dup% met% doy,    mp, kend, alloc )
   call cable_DupVars ( dup% met% hod,    mp, kend, alloc )
   call cable_DupVars ( dup% met% ofsd,   mp, kend, alloc ) 
   call cable_DupVars ( dup% met% fld,    mp, kend, alloc )
   call cable_DupVars ( dup% met% precip, mp, kend, alloc )
   call cable_DupVars ( dup% met% tk,     mp, kend, alloc )
   call cable_DupVars ( dup% met% tvair,  mp, kend, alloc )
   call cable_DupVars ( dup% met% tvrad,  mp, kend, alloc )
   call cable_DupVars ( dup% met% pmb,    mp, kend, alloc )
   call cable_DupVars ( dup% met% ua,     mp, kend, alloc )
   call cable_DupVars ( dup% met% qv,     mp, kend, alloc )
   call cable_DupVars ( dup% met% qvair,  mp, kend, alloc )
   call cable_DupVars ( dup% met% da,     mp, kend, alloc )
   call cable_DupVars ( dup% met% dva,    mp, kend, alloc )
   call cable_DupVars ( dup% met% coszen, mp, kend, alloc )
   call cable_DupVars ( dup% met% fsd,    mp, swb, kend, alloc ) 
   call cable_DupVars ( dup% met% precip_sn, mp, kend, alloc )
   
   call cable_DupVars( dup% rad% extkb,   mp, kend, alloc )  
   call cable_DupVars( dup% rad% extkd2,  mp, kend, alloc )
   call cable_DupVars( dup% rad% extkd,   mp, kend, alloc )
   call cable_DupVars( dup% rad% flws,    mp, kend, alloc )
   call cable_DupVars( dup% rad% lwabv,   mp, kend, alloc )
   call cable_DupVars( dup% rad% qssabs,  mp, kend, alloc )
   call cable_DupVars( dup% rad% transd,  mp, kend, alloc )
   call cable_DupVars( dup% rad% trad,    mp, kend, alloc )
   call cable_DupVars( dup% rad% transb,  mp, kend, alloc )
   call cable_DupVars( dup% rad% workp1,  mp, kend, alloc )
   call cable_DupVars( dup% rad% workp2,  mp, kend, alloc )
   call cable_DupVars( dup% rad% workp3,  mp, kend, alloc )
   call cable_DupVars( dup% rad% albedo_T,   mp,   kend, alloc )
   call cable_DupVars( dup% rad% latitude,   mp,   kend, alloc )
   call cable_DupVars( dup% rad% longitude,  mp,   kend, alloc )
   call cable_DupVars( dup% rad% fvlai,   mp, mf,  kend, alloc )
   call cable_DupVars( dup% rad% rniso,   mp, mf,  kend, alloc )
   call cable_DupVars( dup% rad% scalex,  mp, mf,  kend, alloc )
   call cable_DupVars( dup% rad% gradis,  mp, mf,  kend, alloc )
   call cable_DupVars( dup% rad% albedo,  mp, nrb, kend, alloc ) 
   call cable_DupVars( dup% rad% rhocdf,  mp, nrb, kend, alloc )
   call cable_DupVars( dup% rad% reffdf,  mp, nrb, kend, alloc )
   call cable_DupVars( dup% rad% reffbm,  mp, nrb, kend, alloc )
   call cable_DupVars( dup% rad% extkbm,  mp, nrb, kend, alloc )
   call cable_DupVars( dup% rad% extkdm,  mp, nrb, kend, alloc )
   call cable_DupVars( dup% rad% cexpkbm, mp, swb, kend, alloc )
   call cable_DupVars( dup% rad% cexpkdm, mp, swb, kend, alloc )
   call cable_DupVars( dup% rad% fbeam,   mp, nrb, kend, alloc )
   call cable_DupVars( dup% rad% rhocbm,  mp, nrb, kend, alloc )
   call cable_DupVars( dup% rad% qcan,    mp, mf, nrb, kend, alloc )

   call cable_DupVars ( dup% rough% coexp,   mp, kend, alloc )
   call cable_DupVars ( dup% rough% disp,    mp, kend, alloc )
   call cable_DupVars ( dup% rough% hruff,   mp, kend, alloc )
   call cable_DupVars ( dup% rough% rt0us,   mp, kend, alloc )
   call cable_DupVars ( dup% rough% rt1usa,  mp, kend, alloc )
   call cable_DupVars ( dup% rough% rt1usb,  mp, kend, alloc )
   call cable_DupVars ( dup% rough% rt1,     mp, kend, alloc )
   call cable_DupVars ( dup% rough% term2,   mp, kend, alloc )
   call cable_DupVars ( dup% rough% term3,   mp, kend, alloc )
   call cable_DupVars ( dup% rough% term5,   mp, kend, alloc )
   call cable_DupVars ( dup% rough% term6,   mp, kend, alloc )
   call cable_DupVars ( dup% rough% usuh,    mp, kend, alloc )
   call cable_DupVars ( dup% rough% za_uv,   mp, kend, alloc )
   call cable_DupVars ( dup% rough% za_tq,   mp, kend, alloc )
   call cable_DupVars ( dup% rough% z0m,     mp, kend, alloc )
   call cable_DupVars ( dup% rough% zref_uv, mp, kend, alloc )
   call cable_DupVars ( dup% rough% zref_tq, mp, kend, alloc )
   call cable_DupVars ( dup% rough% zruffs,  mp, kend, alloc )
   call cable_DupVars ( dup% rough% z0soil,  mp, kend, alloc )
   call cable_DupVars ( dup% rough% hruff_grmx, mp, kend, alloc )
   call cable_DupVars ( dup% rough% z0soilsn,   mp, kend, alloc )

   call cable_DupVars ( dup% air% rho,    mp, kend, alloc )
   call cable_DupVars ( dup% air% volm,   mp, kend, alloc )
   call cable_DupVars ( dup% air% rlam,   mp, kend, alloc )
   call cable_DupVars ( dup% air% qsat,   mp, kend, alloc )
   call cable_DupVars ( dup% air% epsi,   mp, kend, alloc )
   call cable_DupVars ( dup% air% visc,   mp, kend, alloc )
   call cable_DupVars ( dup% air% psyc,   mp, kend, alloc )
   call cable_DupVars ( dup% air% dsatdk, mp, kend, alloc )
   call cable_DupVars ( dup% air% cmolar, mp, kend, alloc )

END SUBROUTINE alloc_DupVars 

! ------------------------------------------------------------------------------
  
subroutine Set_DupVars( ktau, dup, met, bgc, soil,veg, ssnow, canopy, &
                        rad, rough, air )

   use cable_def_types_mod, only : met_type, bgc_pool_type, soil_parameter_type, &
                                   veg_parameter_type, soil_snow_type, canopy_type, & 
                                   air_type, radiation_type, roughness_type

   use cable_dupvars_mod, only : cable_setDupVars  
   
   integer :: ktau
   type (TdupVars) :: dup
 
   ! CABLE variables
   TYPE (met_type)       :: met     ! met input variables
   TYPE (canopy_type)    :: canopy  ! vegetation variables
   TYPE (soil_snow_type) :: ssnow   ! soil and snow variables
   TYPE (air_type)       :: air     ! air property variables
   TYPE (radiation_type) :: rad     ! radiation variables
   TYPE (roughness_type) :: rough   ! roughness varibles
   
   ! CABLE parameters
   TYPE (soil_parameter_type) :: soil ! soil parameters	
   TYPE (veg_parameter_type)  :: veg  ! vegetation parameters	 
   
   TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables
 
   call cable_setdupVars( dup% bgc% cplant, bgc% cplant,   ktau, alloc )
   call cable_setdupVars( dup% bgc% csoil,  bgc% csoil,    ktau, alloc )
              
   call cable_setdupVars( dup% soil% bch,   soil% bch,     ktau, alloc )   
   call cable_setdupVars( dup% soil% c3,    soil% c3,      ktau, alloc )    
   call cable_setdupVars( dup% soil% clay,  soil% clay,    ktau, alloc )  
   call cable_setdupVars( dup% soil% css,   soil% css,     ktau, alloc )   
   call cable_setdupVars( dup% soil% hsbh,  soil% hsbh,    ktau, alloc )  
   call cable_setdupVars( dup% soil% hyds,  soil% hyds,    ktau, alloc )  
   call cable_setdupVars( dup% soil% i2bp3, soil% i2bp3,   ktau, alloc ) 
   call cable_setdupVars( dup% soil% ibp2,  soil% ibp2,    ktau, alloc )  
   call cable_setdupVars( dup% soil% sand,  soil% sand,    ktau, alloc )   
   call cable_setdupVars( dup% soil% sfc,   soil% sfc,     ktau, alloc )   
   call cable_setdupVars( dup% soil% silt,  soil% silt,    ktau, alloc )   
   call cable_setdupVars( dup% soil% ssat,  soil% ssat,    ktau, alloc )   
   call cable_setdupVars( dup% soil% sucs,  soil% sucs,    ktau, alloc )   
   call cable_setdupVars( dup% soil% swilt, soil% swilt,   ktau, alloc )  
   call cable_setdupVars( dup% soil% zse,   soil% zse,     ktau, alloc )    
   call cable_setdupVars( dup% soil% zshh, soil% zshh,     ktau, alloc )  
   call cable_setdupVars( dup% soil% cnsd, soil% cnsd,     ktau, alloc )  
   call cable_setdupVars( dup% soil% isoilm,  soil% isoilm,  ktau, alloc )  
   call cable_setdupVars( dup% soil% rhosoil, soil% rhosoil, ktau, alloc )  
   
   call cable_setdupVars( dup% soil% albsoil,  soil% albsoil,  ktau, alloc)   
   call cable_setdupVars( dup% soil% pwb_min,  soil% pwb_min,  ktau, alloc )  
   call cable_setdupVars( dup% soil% albsoilf, soil% albsoilf, ktau, alloc )  
   
   call cable_setDupVars( dup% ssnow% iantrct,   ssnow% iantrct,   ktau, alloc )
   call cable_setDupVars( dup% ssnow % pudsto,   ssnow % pudsto,   ktau, alloc )
   call cable_setDupVars( dup% ssnow % pudsmx,   ssnow % pudsmx,   ktau, alloc )
   call cable_setDupVars( dup% ssnow% cls,       ssnow% cls,       ktau, alloc )     
   call cable_setDupVars( dup% ssnow% dfn_dtg,   ssnow% dfn_dtg,   ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% dfh_dtg,   ssnow% dfh_dtg,   ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% dfe_ddq,   ssnow% dfe_ddq,   ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% ddq_dtg,   ssnow% ddq_dtg,   ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% evapsn,    ssnow% evapsn,    ktau, alloc )  
   call cable_setDupVars( dup% ssnow% fwtop,     ssnow% fwtop,     ktau, alloc )   
   call cable_setDupVars( dup% ssnow% fwtop1,    ssnow% fwtop1,    ktau, alloc )   
   call cable_setDupVars( dup% ssnow% fwtop2,    ssnow% fwtop2,    ktau, alloc )   
   call cable_setDupVars( dup% ssnow% fwtop3,    ssnow% fwtop3,    ktau, alloc )   
   call cable_setDupVars( dup% ssnow% gammzz,    ssnow% gammzz,    ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% isflag,    ssnow% isflag,    ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% osnowd,    ssnow% osnowd,    ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% potev,     ssnow% potev,     ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% runoff,    ssnow% runoff,    ktau, alloc )
   call cable_setDupVars( dup% ssnow% rnof1,     ssnow% rnof1,     ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% rnof2,     ssnow% rnof2,     ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% rtsoil,    ssnow% rtsoil,    ktau, alloc )
   call cable_setDupVars( dup% ssnow% sconds,    ssnow% sconds,    ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% sdepth,    ssnow% sdepth,    ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% smass,     ssnow% smass,     ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% snage,     ssnow% snage,     ktau, alloc )  
   call cable_setDupVars( dup% ssnow% snowd,     ssnow% snowd,     ktau, alloc )  
   call cable_setDupVars( dup% ssnow% smelt,     ssnow% smelt,     ktau, alloc )  
   call cable_setDupVars( dup% ssnow% ssdn,      ssnow% ssdn,      ktau, alloc )  
   call cable_setDupVars( dup% ssnow% ssdnn,     ssnow% ssdnn,     ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% tgg,       ssnow% tgg,       ktau, alloc )   
   call cable_setDupVars( dup% ssnow% tggsn,     ssnow% tggsn,     ktau, alloc )  
   call cable_setDupVars( dup% ssnow% tss,       ssnow% tss,       ktau, alloc )   
   call cable_setDupVars( dup% ssnow% tss_p,     ssnow% tss_p,     ktau, alloc )   
   call cable_setDupVars( dup% ssnow% deltss,    ssnow% deltss,    ktau, alloc )   
   call cable_setDupVars( dup% ssnow% owb1,      ssnow% owb1,      ktau, alloc )   
   call cable_setDupVars( dup% ssnow% wb,        ssnow% wb,        ktau, alloc )    
   call cable_setDupVars( dup% ssnow% wbice,     ssnow% wbice,     ktau, alloc ) 
   call cable_setDupVars( dup% ssnow% wblf,      ssnow% wblf,      ktau, alloc ) 
   call cable_setDupVars( dup% ssnow%wbtot,      ssnow%wbtot,      ktau, alloc )    
   call cable_setDupVars( dup% ssnow%wbtot1,     ssnow%wbtot1,     ktau, alloc )    
   call cable_setDupVars( dup% ssnow%wbtot2,     ssnow%wbtot2,     ktau, alloc )    
   call cable_setDupVars( dup% ssnow%wb_lake,    ssnow%wb_lake,    ktau, alloc )    
   call cable_setDupVars( dup% ssnow%sinfil,     ssnow%sinfil,     ktau, alloc )    
   call cable_setDupVars( dup% ssnow%evapfbl,    ssnow%evapfbl,    ktau, alloc )    
   call cable_setDupVars( dup% ssnow%qstss,      ssnow%qstss,      ktau, alloc )    
   call cable_setDupVars( dup% ssnow%wetfac,     ssnow%wetfac,     ktau, alloc )  
   call cable_setDupVars( dup% ssnow%owetfac,    ssnow%owetfac,    ktau, alloc )  
   call cable_setDupVars( dup% ssnow%t_snwlr,    ssnow%t_snwlr,    ktau, alloc )  
   call cable_setDupVars( dup% ssnow%wbfice,     ssnow%wbfice,     ktau, alloc )  
   call cable_setDupVars( dup% ssnow%tggav,      ssnow%tggav,      ktau, alloc )  
   call cable_setDupVars( dup% ssnow%otgg,       ssnow%otgg,       ktau, alloc )   
   call cable_setDupVars( dup% ssnow%otss,       ssnow%otss,       ktau, alloc )   
   call cable_setDupVars( dup% ssnow%otss_0,     ssnow%otss_0,     ktau, alloc )   
   call cable_setDupVars( dup% ssnow%tprecip,    ssnow%tprecip,    ktau, alloc ) 
   call cable_setDupVars( dup% ssnow%tevap,      ssnow%tevap,      ktau, alloc ) 
   call cable_setDupVars( dup% ssnow%trnoff,     ssnow%trnoff,     ktau, alloc ) 
   call cable_setDupVars( dup% ssnow%totenbal,   ssnow%totenbal,   ktau, alloc ) 
   call cable_setDupVars( dup% ssnow%totenbal2,  ssnow%totenbal2,  ktau, alloc ) 
   call cable_setDupVars( dup% ssnow%fland,      ssnow%fland,      ktau, alloc )      
   call cable_setDupVars( dup% ssnow%ifland,     ssnow%ifland,     ktau, alloc )  
   call cable_setDupVars( dup% ssnow%tilefrac,   ssnow%tilefrac,   ktau, alloc ) 
   call cable_setDupVars( dup% ssnow%qasrf,      ssnow%qasrf,      ktau, alloc )  
   call cable_setDupVars( dup% ssnow%qfsrf,      ssnow%qfsrf,      ktau, alloc )  
   call cable_setDupVars( dup% ssnow%qssrf,      ssnow%qssrf,      ktau, alloc )  

   call cable_setDupVars( dup% ssnow % dtmlt,    ssnow % dtmlt,     ktau, alloc )
   call cable_setDupVars( dup% ssnow% albsoilsn, ssnow% albsoilsn,  ktau, alloc ) 

   call cable_setDupVars( dup% veg% canst1,   veg% canst1,   ktau, alloc ) 
   call cable_setDupVars( dup% veg% dleaf,    veg% dleaf,    ktau, alloc )  
   call cable_setDupVars( dup% veg% ejmax,    veg% ejmax,    ktau, alloc ) 
   call cable_setDupVars( dup% veg% iveg,     veg% iveg,     ktau, alloc ) 
   call cable_setDupVars( dup% veg% meth,     veg% meth,     ktau, alloc ) 
   call cable_setDupVars( dup% veg% frac4,    veg% frac4,    ktau, alloc )  
   call cable_setDupVars( dup% veg% hc,       veg% hc,       ktau, alloc )     
   call cable_setDupVars( dup% veg% vlai,     veg% vlai,     ktau, alloc )   
   call cable_setDupVars( dup% veg% xalbnir,  veg% xalbnir,  ktau, alloc ) 
   call cable_setDupVars( dup% veg% rp20,     veg% rp20,     ktau, alloc )   
   call cable_setDupVars( dup% veg% rpcoef,   veg% rpcoef,   ktau, alloc ) 
   call cable_setDupVars( dup% veg% rs20,     veg% rs20,     ktau, alloc )   
   call cable_setDupVars( dup% veg% shelrb,   veg% shelrb,   ktau, alloc ) 
   call cable_setDupVars( dup% veg% vegcf,    veg% vegcf,    ktau, alloc )  
   call cable_setDupVars( dup% veg% tminvj,   veg% tminvj,   ktau, alloc ) 
   call cable_setDupVars( dup% veg% tmaxvj,   veg% tmaxvj,   ktau, alloc ) 
   call cable_setDupVars( dup% veg% vbeta,    veg% vbeta,    ktau, alloc )  
   call cable_setDupVars( dup% veg% vcmax,    veg% vcmax,    ktau, alloc )  
   call cable_setDupVars( dup% veg% xfang,    veg% xfang,    ktau, alloc )  
   call cable_setDupVars( dup% veg% extkn,    veg% extkn,     ktau, alloc ) 
   call cable_setDupVars( dup% veg% wai,      veg% wai,       ktau, alloc )   
   call cable_setDupVars( dup% veg% deciduous,veg% deciduous, ktau, alloc ) 
   call cable_setDupVars( dup% veg% froot,    veg% froot,     ktau, alloc ) 
   call cable_setDupVars( dup% veg% refl,     veg% refl,      ktau, alloc ) !jhan:swb?
   call cable_setDupVars( dup% veg% taul,     veg% taul,      ktau, alloc ) 
   call cable_setDupVars( dup% veg% vlaimax,  veg% vlaimax,   ktau, alloc ) 

   call cable_setDupVars( dup% canopy% fess,   canopy% fess,   ktau, alloc )
   call cable_setDupVars( dup% canopy% fesp,   canopy% fesp,   ktau, alloc )
   call cable_setDupVars( dup% canopy% cansto, canopy% cansto, ktau, alloc )  
   call cable_setDupVars( dup% canopy% cduv,   canopy% cduv,   ktau, alloc )   
   call cable_setDupVars( dup% canopy% delwc,  canopy% delwc,  ktau, alloc )  
   call cable_setDupVars( dup% canopy% dewmm,  canopy% dewmm,  ktau, alloc )  
   call cable_setDupVars( dup% canopy% dgdtg,  canopy% dgdtg,  ktau, alloc )  
   call cable_setDupVars( dup% canopy% fe,     canopy% fe,     ktau, alloc )      
   call cable_setDupVars( dup% canopy% fh,     canopy% fh,     ktau, alloc )      
   call cable_setDupVars( dup% canopy% fpn,    canopy% fpn,    ktau, alloc )     
   call cable_setDupVars( dup% canopy% frp,    canopy% frp,    ktau, alloc )     
   call cable_setDupVars( dup% canopy% frpw,   canopy% frpw,   ktau, alloc )    
   call cable_setDupVars( dup% canopy% frpr,   canopy% frpr,   ktau, alloc )    
   call cable_setDupVars( dup% canopy% frs,    canopy% frs,    ktau, alloc )     
   call cable_setDupVars( dup% canopy% fnee,   canopy% fnee,   ktau, alloc )    
   call cable_setDupVars( dup% canopy% frday,  canopy% frday,  ktau, alloc )   
   call cable_setDupVars( dup% canopy% fnv,    canopy% fnv,    ktau, alloc )     
   call cable_setDupVars( dup% canopy% fev,    canopy% fev,    ktau, alloc )     
   call cable_setDupVars( dup% canopy% fevc,   canopy% fevc,   ktau, alloc )    
   call cable_setDupVars( dup% canopy% fhv,    canopy% fhv,    ktau, alloc )     
   call cable_setDupVars( dup% canopy% fns,    canopy% fns,    ktau, alloc )     
   call cable_setDupVars( dup% canopy% fhs,    canopy% fhs,    ktau, alloc )     
   call cable_setDupVars( dup% canopy% fhs_cor,canopy% fhs_cor,ktau, alloc )     
   call cable_setDupVars( dup% canopy% ga,     canopy% ga,     ktau, alloc )      
   call cable_setDupVars( dup% canopy% ghflux, canopy% ghflux, ktau, alloc )   
   call cable_setDupVars( dup% canopy% precis, canopy% precis, ktau, alloc ) 
   call cable_setDupVars( dup% canopy% qscrn,  canopy% qscrn,  ktau, alloc )  
   call cable_setDupVars( dup% canopy% rnet,   canopy% rnet,   ktau, alloc )   
   call cable_setDupVars( dup% canopy% segg,   canopy% segg,   ktau, alloc )   
   call cable_setDupVars( dup% canopy% sghflux,canopy% sghflux,ktau, alloc )  
   call cable_setDupVars( dup% canopy% through,canopy% through,ktau, alloc )  
   call cable_setDupVars( dup% canopy% spill,  canopy% spill,  ktau, alloc )  
   call cable_setDupVars( dup% canopy% tscrn,  canopy% tscrn,  ktau, alloc )  
   call cable_setDupVars( dup% canopy% wcint,  canopy% wcint,  ktau, alloc )  
   call cable_setDupVars( dup% canopy% tv,     canopy% tv,     ktau, alloc )      
   call cable_setDupVars( dup% canopy% us,     canopy% us,     ktau, alloc )      
   call cable_setDupVars( dup% canopy% uscrn,  canopy% uscrn,  ktau, alloc )   
   call cable_setDupVars( dup% canopy% rghlai, canopy% rghlai, ktau, alloc ) 
   call cable_setDupVars( dup% canopy% vlaiw,  canopy% vlaiw,  ktau, alloc ) 
   call cable_setDupVars( dup% canopy% fwet,   canopy% fwet,   ktau, alloc )   
   call cable_setDupVars( dup% canopy% epot,   canopy% epot,   ktau, alloc )   
   call cable_setDupVars( dup% canopy% fnpp,   canopy% fnpp,   ktau, alloc )   
   call cable_setDupVars( dup% canopy% gswx_T, canopy% gswx_T, ktau, alloc )  
   call cable_setDupVars( dup% canopy% cdtq,   canopy% cdtq,   ktau, alloc )   
   call cable_setDupVars( dup% canopy% fevw,   canopy% fevw,   ktau, alloc )   
   call cable_setDupVars( dup% canopy% fhvw,   canopy% fhvw,   ktau, alloc )   
   call cable_setDupVars( dup% canopy% fes,    canopy% fes,    ktau, alloc )    
   call cable_setDupVars( dup% canopy% fes_cor,canopy% fes_cor,ktau, alloc )    
   call cable_setDupVars( dup% canopy% gswx,   canopy% gswx,   ktau, alloc )  
   call cable_setDupVars( dup% canopy% evapfbl,canopy% evapfbl,ktau, alloc )
   call cable_setDupVars( dup% canopy% wetfac_cs, canopy% wetfac_cs, ktau, alloc )  
   call cable_setDupVars( dup% canopy% oldcansto, canopy% oldcansto, ktau, alloc )  
   call cable_setDupVars( dup% canopy% zetar,     canopy% zetar,     ktau, alloc )  
   call cable_setDupVars( dup% canopy% fevw_pot,  canopy% fevw_pot,  ktau, alloc )  
 
   call cable_setDupVars ( dup% met% ca,    met% ca,     ktau, alloc )
   call cable_setDupVars ( dup% met% year,  met% year,   ktau, alloc )
   call cable_setDupVars ( dup% met% moy,   met% moy,    ktau, alloc )
   call cable_setDupVars ( dup% met% doy,   met% doy,    ktau, alloc )
   call cable_setDupVars ( dup% met% hod,   met% hod,    ktau, alloc )
   call cable_setDupVars ( dup% met% ofsd,  met% ofsd,   ktau, alloc ) 
   call cable_setDupVars ( dup% met% fld,   met% fld,    ktau, alloc )
   call cable_setDupVars ( dup% met% precip,met% precip, ktau, alloc )
   call cable_setDupVars ( dup% met% tk,    met% tk,     ktau, alloc )
   call cable_setDupVars ( dup% met% tvair, met% tvair,  ktau, alloc )
   call cable_setDupVars ( dup% met% tvrad, met% tvrad,  ktau, alloc )
   call cable_setDupVars ( dup% met% pmb,   met% pmb,    ktau, alloc )
   call cable_setDupVars ( dup% met% ua,    met% ua,     ktau, alloc )
   call cable_setDupVars ( dup% met% qv,    met% qv,     ktau, alloc )
   call cable_setDupVars ( dup% met% qvair, met% qvair,  ktau, alloc )
   call cable_setDupVars ( dup% met% da,    met% da,     ktau, alloc )
   call cable_setDupVars ( dup% met% dva,   met% dva,    ktau, alloc )
   call cable_setDupVars ( dup% met% coszen,met% coszen, ktau, alloc )
   call cable_setDupVars ( dup% met% fsd,   met% fsd,    ktau, alloc ) 
   call cable_setDupVars ( dup% met% precip_sn,  met% precip_sn,   ktau, alloc )
! ------------------------------------------------------------------------------
   call cable_setDupVars( dup% rad% extkb,   rad% extkb,   ktau, alloc )  
   call cable_setDupVars( dup% rad% extkd2,  rad% extkd2,  ktau, alloc )
   call cable_setDupVars( dup% rad% extkd,   rad% extkd,   ktau, alloc )
   call cable_setDupVars( dup% rad% flws,    rad% flws,    ktau, alloc )
   call cable_setDupVars( dup% rad% lwabv,   rad% lwabv,   ktau, alloc )
   call cable_setDupVars( dup% rad% qssabs,  rad% qssabs,  ktau, alloc )
   call cable_setDupVars( dup% rad% transd,  rad% transd,  ktau, alloc )
   call cable_setDupVars( dup% rad% trad,    rad% trad,    ktau, alloc )
   call cable_setDupVars( dup% rad% transb,  rad% transb,  ktau, alloc )
   call cable_setDupVars( dup% rad% workp1,  rad% workp1,  ktau, alloc )
   call cable_setDupVars( dup% rad% workp2,  rad% workp2,  ktau, alloc )
   call cable_setDupVars( dup% rad% workp3,  rad% workp3,  ktau, alloc )
   call cable_setDupVars( dup% rad% albedo_T,   rad% albedo_T,   ktau, alloc )
   call cable_setDupVars( dup% rad% latitude,   rad% latitude,   ktau, alloc )
   call cable_setDupVars( dup% rad% longitude,  rad% longitude,  ktau, alloc )
   call cable_setDupVars( dup% rad% fvlai,   rad% fvlai,   ktau, alloc )
   call cable_setDupVars( dup% rad% rniso,   rad% rniso,   ktau, alloc )
   call cable_setDupVars( dup% rad% scalex,  rad% scalex,  ktau, alloc )
   call cable_setDupVars( dup% rad% gradis,  rad% gradis,  ktau, alloc )
   call cable_setDupVars( dup% rad% albedo,  rad% albedo,  ktau, alloc ) 
   call cable_setDupVars( dup% rad% rhocdf,  rad% rhocdf,  ktau, alloc )
   call cable_setDupVars( dup% rad% reffdf,  rad% reffdf,  ktau, alloc )
   call cable_setDupVars( dup% rad% reffbm,  rad% reffbm,  ktau, alloc )
   call cable_setDupVars( dup% rad% extkbm,  rad% extkbm,  ktau, alloc )
   call cable_setDupVars( dup% rad% extkdm,  rad% extkdm,  ktau, alloc )
   call cable_setDupVars( dup% rad% cexpkbm, rad% cexpkbm, ktau, alloc )
   call cable_setDupVars( dup% rad% cexpkdm, rad% cexpkdm, ktau, alloc )
   call cable_setDupVars( dup% rad% fbeam,   rad% fbeam,   ktau, alloc )
   call cable_setDupVars( dup% rad% rhocbm,  rad% rhocbm,  ktau, alloc )
   call cable_setDupVars( dup% rad% qcan,    rad% qcan,    ktau, alloc )

   call cable_setDupVars ( dup% rough% coexp,   rough% coexp,   ktau, alloc )
   call cable_setDupVars ( dup% rough% disp,    rough% disp,    ktau, alloc )
   call cable_setDupVars ( dup% rough% hruff,   rough% hruff,   ktau, alloc )
   call cable_setDupVars ( dup% rough% rt0us,   rough% rt0us,   ktau, alloc )
   call cable_setDupVars ( dup% rough% rt1usa,  rough% rt1usa,  ktau, alloc )
   call cable_setDupVars ( dup% rough% rt1usb,  rough% rt1usb,  ktau, alloc )
   call cable_setDupVars ( dup% rough% rt1,     rough% rt1,     ktau, alloc )
   call cable_setDupVars ( dup% rough% term2,   rough% term2,   ktau, alloc )
   call cable_setDupVars ( dup% rough% term3,   rough% term3,   ktau, alloc )
   call cable_setDupVars ( dup% rough% term5,   rough% term5,   ktau, alloc )
   call cable_setDupVars ( dup% rough% term6,   rough% term6,   ktau, alloc )
   call cable_setDupVars ( dup% rough% usuh,    rough% usuh,    ktau, alloc )
   call cable_setDupVars ( dup% rough% za_uv,   rough% za_uv,   ktau, alloc )
   call cable_setDupVars ( dup% rough% za_tq,   rough% za_tq,   ktau, alloc )
   call cable_setDupVars ( dup% rough% z0m,     rough% z0m,     ktau, alloc )
   call cable_setDupVars ( dup% rough% zref_uv, rough% zref_uv, ktau, alloc )
   call cable_setDupVars ( dup% rough% zref_tq, rough% zref_tq, ktau, alloc )
   call cable_setDupVars ( dup% rough% zruffs,  rough% zruffs,  ktau, alloc )
   call cable_setDupVars ( dup% rough% z0soil,  rough% z0soil,  ktau, alloc )
   call cable_setDupVars ( dup% rough% hruff_grmx, rough% hruff_grmx, ktau, alloc )
   call cable_setDupVars ( dup% rough% z0soilsn,   rough% z0soilsn,   ktau, alloc )

   call cable_setDupVars ( dup% air% rho,    air% rho,    ktau, alloc )
   call cable_setDupVars ( dup% air% volm,   air% volm,   ktau, alloc )
   call cable_setDupVars ( dup% air% rlam,   air% rlam,   ktau, alloc )
   call cable_setDupVars ( dup% air% qsat,   air% qsat,   ktau, alloc )
   call cable_setDupVars ( dup% air% epsi,   air% epsi,   ktau, alloc )
   call cable_setDupVars ( dup% air% visc,   air% visc,   ktau, alloc )
   call cable_setDupVars ( dup% air% psyc,   air% psyc,   ktau, alloc )
   call cable_setDupVars ( dup% air% dsatdk, air% dsatdk, ktau, alloc )
   call cable_setDupVars ( dup% air% cmolar, air% cmolar, ktau, alloc )

! ------------------------------------------------------------------------------

END SUBROUTINE set_DupVars 

! ------------------------------------------------------------------------------

  
subroutine reSet_DupVars( ktau, dup, met, bgc, soil,veg, ssnow, canopy, &
                        rad, rough, air )

   use cable_def_types_mod, only : met_type, bgc_pool_type, soil_parameter_type, &
                                   veg_parameter_type, soil_snow_type, canopy_type, & 
                                   air_type, radiation_type, roughness_type
   use cable_dupvars_mod, only : cable_setDupVars  
   
   integer :: ktau
   type (TdupVars) :: dup
   ! CABLE variables
   TYPE (met_type)       :: met     ! met input variables
   TYPE (canopy_type)    :: canopy  ! vegetation variables
   TYPE (soil_snow_type) :: ssnow   ! soil and snow variables
   TYPE (air_type)       :: air     ! air property variables
   TYPE (radiation_type) :: rad     ! radiation variables
   TYPE (roughness_type) :: rough   ! roughness varibles
   
   ! CABLE parameters
   TYPE (soil_parameter_type) :: soil ! soil parameters	
   TYPE (veg_parameter_type)  :: veg  ! vegetation parameters	 
   
   TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables
 
   call cable_resetdupVars( dup% bgc% cplant, bgc% cplant,   ktau, alloc )
   call cable_resetdupVars( dup% bgc% csoil,  bgc% csoil,    ktau, alloc )
              
   call cable_resetdupVars( dup% soil% bch,   soil% bch,     ktau, alloc )   
   call cable_resetdupVars( dup% soil% c3,    soil% c3,      ktau, alloc )    
   call cable_resetdupVars( dup% soil% clay,  soil% clay,    ktau, alloc )  
   call cable_resetdupVars( dup% soil% css,   soil% css,     ktau, alloc )   
   call cable_resetdupVars( dup% soil% hsbh,  soil% hsbh,    ktau, alloc )  
   call cable_resetdupVars( dup% soil% hyds,  soil% hyds,    ktau, alloc )  
   call cable_resetdupVars( dup% soil% i2bp3, soil% i2bp3,   ktau, alloc ) 
   call cable_resetdupVars( dup% soil% ibp2,  soil% ibp2,    ktau, alloc )  
   call cable_resetdupVars( dup% soil% sand,  soil% sand,    ktau, alloc )   
   call cable_resetdupVars( dup% soil% sfc,   soil% sfc,     ktau, alloc )   
   call cable_resetdupVars( dup% soil% silt,  soil% silt,    ktau, alloc )   
   call cable_resetdupVars( dup% soil% ssat,  soil% ssat,    ktau, alloc )   
   call cable_resetdupVars( dup% soil% sucs,  soil% sucs,    ktau, alloc )   
   call cable_resetdupVars( dup% soil% swilt, soil% swilt,   ktau, alloc )  
   call cable_resetdupVars( dup% soil% zse,   soil% zse,     ktau, alloc )    
   call cable_resetdupVars( dup% soil% zshh, soil% zshh,     ktau, alloc )  
   call cable_resetdupVars( dup% soil% cnsd, soil% cnsd,     ktau, alloc )  
   call cable_resetdupVars( dup% soil% isoilm,  soil% isoilm,  ktau, alloc )  
   call cable_resetdupVars( dup% soil% rhosoil, soil% rhosoil, ktau, alloc )  
   
   call cable_resetdupVars( dup% soil% albsoil,  soil% albsoil,  ktau, alloc)   
   call cable_resetdupVars( dup% soil% pwb_min,  soil% pwb_min,  ktau, alloc )  
   call cable_resetdupVars( dup% soil% albsoilf, soil% albsoilf, ktau, alloc )  
   
   call cable_resetDupVars( dup% ssnow% iantrct,   ssnow% iantrct,   ktau, alloc )
   call cable_resetDupVars( dup% ssnow % pudsto,   ssnow % pudsto,   ktau, alloc )
   call cable_resetDupVars( dup% ssnow % pudsmx,   ssnow % pudsmx,   ktau, alloc )
   call cable_resetDupVars( dup% ssnow% cls,       ssnow% cls,       ktau, alloc )     
   call cable_resetDupVars( dup% ssnow% dfn_dtg,   ssnow% dfn_dtg,   ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% dfh_dtg,   ssnow% dfh_dtg,   ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% dfe_ddq,   ssnow% dfe_ddq,   ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% ddq_dtg,   ssnow% ddq_dtg,   ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% evapsn,    ssnow% evapsn,    ktau, alloc )  
   call cable_resetDupVars( dup% ssnow% fwtop,     ssnow% fwtop,     ktau, alloc )   
   call cable_resetDupVars( dup% ssnow% fwtop1,    ssnow% fwtop1,    ktau, alloc )   
   call cable_resetDupVars( dup% ssnow% fwtop2,    ssnow% fwtop2,    ktau, alloc )   
   call cable_resetDupVars( dup% ssnow% fwtop3,    ssnow% fwtop3,    ktau, alloc )   
   call cable_resetDupVars( dup% ssnow% gammzz,    ssnow% gammzz,    ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% isflag,    ssnow% isflag,    ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% osnowd,    ssnow% osnowd,    ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% potev,     ssnow% potev,     ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% runoff,    ssnow% runoff,    ktau, alloc )
   call cable_resetDupVars( dup% ssnow% rnof1,     ssnow% rnof1,     ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% rnof2,     ssnow% rnof2,     ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% rtsoil,    ssnow% rtsoil,    ktau, alloc )
   call cable_resetDupVars( dup% ssnow% sconds,    ssnow% sconds,    ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% sdepth,    ssnow% sdepth,    ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% smass,     ssnow% smass,     ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% snage,     ssnow% snage,     ktau, alloc )  
   call cable_resetDupVars( dup% ssnow% snowd,     ssnow% snowd,     ktau, alloc )  
   call cable_resetDupVars( dup% ssnow% smelt,     ssnow% smelt,     ktau, alloc )  
   call cable_resetDupVars( dup% ssnow% ssdn,      ssnow% ssdn,      ktau, alloc )  
   call cable_resetDupVars( dup% ssnow% ssdnn,     ssnow% ssdnn,     ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% tgg,       ssnow% tgg,       ktau, alloc )   
   call cable_resetDupVars( dup% ssnow% tggsn,     ssnow% tggsn,     ktau, alloc )  
   call cable_resetDupVars( dup% ssnow% tss,       ssnow% tss,       ktau, alloc )   
   call cable_resetDupVars( dup% ssnow% tss_p,     ssnow% tss_p,     ktau, alloc )   
   call cable_resetDupVars( dup% ssnow% deltss,    ssnow% deltss,    ktau, alloc )   
   call cable_resetDupVars( dup% ssnow% owb1,      ssnow% owb1,      ktau, alloc )   
   call cable_resetDupVars( dup% ssnow% wb,        ssnow% wb,        ktau, alloc )    
   call cable_resetDupVars( dup% ssnow% wbice,     ssnow% wbice,     ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow% wblf,      ssnow% wblf,      ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow%wbtot,      ssnow%wbtot,      ktau, alloc )    
   call cable_resetDupVars( dup% ssnow%wbtot1,     ssnow%wbtot1,     ktau, alloc )    
   call cable_resetDupVars( dup% ssnow%wbtot2,     ssnow%wbtot2,     ktau, alloc )    
   call cable_resetDupVars( dup% ssnow%wb_lake,    ssnow%wb_lake,    ktau, alloc )    
   call cable_resetDupVars( dup% ssnow%sinfil,     ssnow%sinfil,     ktau, alloc )    
   call cable_resetDupVars( dup% ssnow%evapfbl,    ssnow%evapfbl,    ktau, alloc )    
   call cable_resetDupVars( dup% ssnow%qstss,      ssnow%qstss,      ktau, alloc )    
   call cable_resetDupVars( dup% ssnow%wetfac,     ssnow%wetfac,     ktau, alloc )  
   call cable_resetDupVars( dup% ssnow%owetfac,    ssnow%owetfac,    ktau, alloc )  
   call cable_resetDupVars( dup% ssnow%t_snwlr,    ssnow%t_snwlr,    ktau, alloc )  
   call cable_resetDupVars( dup% ssnow%wbfice,     ssnow%wbfice,     ktau, alloc )  
   call cable_resetDupVars( dup% ssnow%tggav,      ssnow%tggav,      ktau, alloc )  
   call cable_resetDupVars( dup% ssnow%otgg,       ssnow%otgg,       ktau, alloc )   
   call cable_resetDupVars( dup% ssnow%otss,       ssnow%otss,       ktau, alloc )   
   call cable_resetDupVars( dup% ssnow%otss_0,     ssnow%otss_0,     ktau, alloc )   
   call cable_resetDupVars( dup% ssnow%tprecip,    ssnow%tprecip,    ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow%tevap,      ssnow%tevap,      ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow%trnoff,     ssnow%trnoff,     ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow%totenbal,   ssnow%totenbal,   ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow%totenbal2,  ssnow%totenbal2,  ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow%fland,      ssnow%fland,      ktau, alloc )      
   call cable_resetDupVars( dup% ssnow%ifland,     ssnow%ifland,     ktau, alloc )  
   call cable_resetDupVars( dup% ssnow%tilefrac,   ssnow%tilefrac,   ktau, alloc ) 
   call cable_resetDupVars( dup% ssnow%qasrf,      ssnow%qasrf,      ktau, alloc )  
   call cable_resetDupVars( dup% ssnow%qfsrf,      ssnow%qfsrf,      ktau, alloc )  
   call cable_resetDupVars( dup% ssnow%qssrf,      ssnow%qssrf,      ktau, alloc )  

   call cable_resetDupVars( dup% ssnow % dtmlt,    ssnow % dtmlt,     ktau, alloc )
   call cable_resetDupVars( dup% ssnow% albsoilsn, ssnow% albsoilsn,  ktau, alloc ) 

   call cable_resetDupVars( dup% veg% canst1,   veg% canst1,   ktau, alloc ) 
  
   !call check_dup_var( dup% veg% canst1(:,ktau),   veg% canst1,   ktau )
   call check_dup_var( ) 
   
   call cable_resetDupVars( dup% veg% dleaf,    veg% dleaf,    ktau, alloc )  
   call cable_resetDupVars( dup% veg% ejmax,    veg% ejmax,    ktau, alloc ) 
   call cable_resetDupVars( dup% veg% iveg,     veg% iveg,     ktau, alloc ) 
   call cable_resetDupVars( dup% veg% meth,     veg% meth,     ktau, alloc ) 
   call cable_resetDupVars( dup% veg% frac4,    veg% frac4,    ktau, alloc )  
   call cable_resetDupVars( dup% veg% hc,       veg% hc,       ktau, alloc )     
   call cable_resetDupVars( dup% veg% vlai,     veg% vlai,     ktau, alloc )   
   call cable_resetDupVars( dup% veg% xalbnir,  veg% xalbnir,  ktau, alloc ) 
   call cable_resetDupVars( dup% veg% rp20,     veg% rp20,     ktau, alloc )   
   call cable_resetDupVars( dup% veg% rpcoef,   veg% rpcoef,   ktau, alloc ) 
   call cable_resetDupVars( dup% veg% rs20,     veg% rs20,     ktau, alloc )   
   call cable_resetDupVars( dup% veg% shelrb,   veg% shelrb,   ktau, alloc ) 
   call cable_resetDupVars( dup% veg% vegcf,    veg% vegcf,    ktau, alloc )  
   call cable_resetDupVars( dup% veg% tminvj,   veg% tminvj,   ktau, alloc ) 
   call cable_resetDupVars( dup% veg% tmaxvj,   veg% tmaxvj,   ktau, alloc ) 
   call cable_resetDupVars( dup% veg% vbeta,    veg% vbeta,    ktau, alloc )  
   call cable_resetDupVars( dup% veg% vcmax,    veg% vcmax,    ktau, alloc )  
   call cable_resetDupVars( dup% veg% xfang,    veg% xfang,    ktau, alloc )  
   call cable_resetDupVars( dup% veg% extkn,    veg% extkn,     ktau, alloc ) 
   call cable_resetDupVars( dup% veg% wai,      veg% wai,       ktau, alloc )   
   call cable_resetDupVars( dup% veg% deciduous,veg% deciduous, ktau, alloc ) 
   call cable_resetDupVars( dup% veg% froot,    veg% froot,     ktau, alloc ) 
   call cable_resetDupVars( dup% veg% refl,     veg% refl,      ktau, alloc ) !jhan:swb?
   call cable_resetDupVars( dup% veg% taul,     veg% taul,      ktau, alloc ) 
   call cable_resetDupVars( dup% veg% vlaimax,  veg% vlaimax,   ktau, alloc ) 

   call cable_resetDupVars( dup% canopy% fess,   canopy% fess,   ktau, alloc )
   call cable_resetDupVars( dup% canopy% fesp,   canopy% fesp,   ktau, alloc )
   call cable_resetDupVars( dup% canopy% cansto, canopy% cansto, ktau, alloc )  
   call cable_resetDupVars( dup% canopy% cduv,   canopy% cduv,   ktau, alloc )   
   call cable_resetDupVars( dup% canopy% delwc,  canopy% delwc,  ktau, alloc )  
   call cable_resetDupVars( dup% canopy% dewmm,  canopy% dewmm,  ktau, alloc )  
   call cable_resetDupVars( dup% canopy% dgdtg,  canopy% dgdtg,  ktau, alloc )  
   call cable_resetDupVars( dup% canopy% fe,     canopy% fe,     ktau, alloc )      
   call cable_resetDupVars( dup% canopy% fh,     canopy% fh,     ktau, alloc )      
   call cable_resetDupVars( dup% canopy% fpn,    canopy% fpn,    ktau, alloc )     
   call cable_resetDupVars( dup% canopy% frp,    canopy% frp,    ktau, alloc )     
   call cable_resetDupVars( dup% canopy% frpw,   canopy% frpw,   ktau, alloc )    
   call cable_resetDupVars( dup% canopy% frpr,   canopy% frpr,   ktau, alloc )    
   call cable_resetDupVars( dup% canopy% frs,    canopy% frs,    ktau, alloc )     
   call cable_resetDupVars( dup% canopy% fnee,   canopy% fnee,   ktau, alloc )    
   call cable_resetDupVars( dup% canopy% frday,  canopy% frday,  ktau, alloc )   
   call cable_resetDupVars( dup% canopy% fnv,    canopy% fnv,    ktau, alloc )     
   call cable_resetDupVars( dup% canopy% fev,    canopy% fev,    ktau, alloc )     
   call cable_resetDupVars( dup% canopy% fevc,   canopy% fevc,   ktau, alloc )    
   call cable_resetDupVars( dup% canopy% fhv,    canopy% fhv,    ktau, alloc )     
   call cable_resetDupVars( dup% canopy% fns,    canopy% fns,    ktau, alloc )     
   call cable_resetDupVars( dup% canopy% fhs,    canopy% fhs,    ktau, alloc )     
   call cable_resetDupVars( dup% canopy% fhs_cor,canopy% fhs_cor,ktau, alloc )     
   call cable_resetDupVars( dup% canopy% ga,     canopy% ga,     ktau, alloc )      
   call cable_resetDupVars( dup% canopy% ghflux, canopy% ghflux, ktau, alloc )   
   call cable_resetDupVars( dup% canopy% precis, canopy% precis, ktau, alloc ) 
   call cable_resetDupVars( dup% canopy% qscrn,  canopy% qscrn,  ktau, alloc )  
   call cable_resetDupVars( dup% canopy% rnet,   canopy% rnet,   ktau, alloc )   
   call cable_resetDupVars( dup% canopy% segg,   canopy% segg,   ktau, alloc )   
   call cable_resetDupVars( dup% canopy% sghflux,canopy% sghflux,ktau, alloc )  
   call cable_resetDupVars( dup% canopy% through,canopy% through,ktau, alloc )  
   call cable_resetDupVars( dup% canopy% spill,  canopy% spill,  ktau, alloc )  
   call cable_resetDupVars( dup% canopy% tscrn,  canopy% tscrn,  ktau, alloc )  
   call cable_resetDupVars( dup% canopy% wcint,  canopy% wcint,  ktau, alloc )  
   call cable_resetDupVars( dup% canopy% tv,     canopy% tv,     ktau, alloc )      
   call cable_resetDupVars( dup% canopy% us,     canopy% us,     ktau, alloc )      
   call cable_resetDupVars( dup% canopy% uscrn,  canopy% uscrn,  ktau, alloc )   
   call cable_resetDupVars( dup% canopy% rghlai, canopy% rghlai, ktau, alloc ) 
   call cable_resetDupVars( dup% canopy% vlaiw,  canopy% vlaiw,  ktau, alloc ) 
   call cable_resetDupVars( dup% canopy% fwet,   canopy% fwet,   ktau, alloc )   
   call cable_resetDupVars( dup% canopy% epot,   canopy% epot,   ktau, alloc )   
   call cable_resetDupVars( dup% canopy% fnpp,   canopy% fnpp,   ktau, alloc )   
   call cable_resetDupVars( dup% canopy% gswx_T, canopy% gswx_T, ktau, alloc )  
   call cable_resetDupVars( dup% canopy% cdtq,   canopy% cdtq,   ktau, alloc )   
   call cable_resetDupVars( dup% canopy% fevw,   canopy% fevw,   ktau, alloc )   
   call cable_resetDupVars( dup% canopy% fhvw,   canopy% fhvw,   ktau, alloc )   
   call cable_resetDupVars( dup% canopy% fes,    canopy% fes,    ktau, alloc )    
   call cable_resetDupVars( dup% canopy% fes_cor,canopy% fes_cor,ktau, alloc )    
   call cable_resetDupVars( dup% canopy% gswx,   canopy% gswx,   ktau, alloc )  
   call cable_resetDupVars( dup% canopy% evapfbl,canopy% evapfbl,ktau, alloc )
   call cable_resetDupVars( dup% canopy% wetfac_cs, canopy% wetfac_cs, ktau, alloc )  
   call cable_resetDupVars( dup% canopy% oldcansto, canopy% oldcansto, ktau, alloc )  
   call cable_resetDupVars( dup% canopy% zetar,     canopy% zetar,     ktau, alloc )  
   call cable_resetDupVars( dup% canopy% fevw_pot,  canopy% fevw_pot,  ktau, alloc )  
 
   call cable_resetDupVars ( dup% met% ca,    met% ca,     ktau, alloc )
   call cable_resetDupVars ( dup% met% year,  met% year,   ktau, alloc )
   call cable_resetDupVars ( dup% met% moy,   met% moy,    ktau, alloc )
   call cable_resetDupVars ( dup% met% doy,   met% doy,    ktau, alloc )
   call cable_resetDupVars ( dup% met% hod,   met% hod,    ktau, alloc )
   call cable_resetDupVars ( dup% met% ofsd,  met% ofsd,   ktau, alloc ) 
   call cable_resetDupVars ( dup% met% fld,   met% fld,    ktau, alloc )
   call cable_resetDupVars ( dup% met% precip,met% precip, ktau, alloc )
   call cable_resetDupVars ( dup% met% tk,    met% tk,     ktau, alloc )
   call cable_resetDupVars ( dup% met% tvair, met% tvair,  ktau, alloc )
   call cable_resetDupVars ( dup% met% tvrad, met% tvrad,  ktau, alloc )
   call cable_resetDupVars ( dup% met% pmb,   met% pmb,    ktau, alloc )
   call cable_resetDupVars ( dup% met% ua,    met% ua,     ktau, alloc )
   call cable_resetDupVars ( dup% met% qv,    met% qv,     ktau, alloc )
   call cable_resetDupVars ( dup% met% qvair, met% qvair,  ktau, alloc )
   call cable_resetDupVars ( dup% met% da,    met% da,     ktau, alloc )
   call cable_resetDupVars ( dup% met% dva,   met% dva,    ktau, alloc )
   call cable_resetDupVars ( dup% met% coszen,met% coszen, ktau, alloc )
   call cable_resetDupVars ( dup% met% fsd,   met% fsd,    ktau, alloc ) 
   call cable_resetDupVars ( dup% met% precip_sn,  met% precip_sn,   ktau, alloc )
! ------------------------------------------------------------------------------
   
   call cable_resetDupVars( dup% rad% extkb,   rad% extkb,   ktau, alloc )  
   call cable_resetDupVars( dup% rad% extkd2,  rad% extkd2,  ktau, alloc )
   call cable_resetDupVars( dup% rad% extkd,   rad% extkd,   ktau, alloc )
   call cable_resetDupVars( dup% rad% flws,    rad% flws,    ktau, alloc )
   call cable_resetDupVars( dup% rad% lwabv,   rad% lwabv,   ktau, alloc )
   call cable_resetDupVars( dup% rad% qssabs,  rad% qssabs,  ktau, alloc )
   call cable_resetDupVars( dup% rad% transd,  rad% transd,  ktau, alloc )
   call cable_resetDupVars( dup% rad% trad,    rad% trad,    ktau, alloc )
   call cable_resetDupVars( dup% rad% transb,  rad% transb,  ktau, alloc )
   call cable_resetDupVars( dup% rad% workp1,  rad% workp1,  ktau, alloc )
   call cable_resetDupVars( dup% rad% workp2,  rad% workp2,  ktau, alloc )
   call cable_resetDupVars( dup% rad% workp3,  rad% workp3,  ktau, alloc )
   call cable_resetDupVars( dup% rad% albedo_T,   rad% albedo_T,   ktau, alloc )
   call cable_resetDupVars( dup% rad% latitude,   rad% latitude,   ktau, alloc )
   call cable_resetDupVars( dup% rad% longitude,  rad% longitude,  ktau, alloc )
   call cable_resetDupVars( dup% rad% fvlai,   rad% fvlai,   ktau, alloc )
   call cable_resetDupVars( dup% rad% rniso,   rad% rniso,   ktau, alloc )
   call cable_resetDupVars( dup% rad% scalex,  rad% scalex,  ktau, alloc )
   call cable_resetDupVars( dup% rad% gradis,  rad% gradis,  ktau, alloc )
   call cable_resetDupVars( dup% rad% albedo,  rad% albedo,  ktau, alloc ) 
   call cable_resetDupVars( dup% rad% rhocdf,  rad% rhocdf,  ktau, alloc )
   call cable_resetDupVars( dup% rad% reffdf,  rad% reffdf,  ktau, alloc )
   call cable_resetDupVars( dup% rad% reffbm,  rad% reffbm,  ktau, alloc )
   call cable_resetDupVars( dup% rad% extkbm,  rad% extkbm,  ktau, alloc )
   call cable_resetDupVars( dup% rad% extkdm,  rad% extkdm,  ktau, alloc )
   call cable_resetDupVars( dup% rad% cexpkbm, rad% cexpkbm, ktau, alloc )
   call cable_resetDupVars( dup% rad% cexpkdm, rad% cexpkdm, ktau, alloc )
   call cable_resetDupVars( dup% rad% fbeam,   rad% fbeam,   ktau, alloc )
   call cable_resetDupVars( dup% rad% rhocbm,  rad% rhocbm,  ktau, alloc )
   call cable_resetDupVars( dup% rad% qcan,    rad% qcan,    ktau, alloc )

   call cable_resetDupVars ( dup% rough% coexp,   rough% coexp,   ktau, alloc )
   call cable_resetDupVars ( dup% rough% disp,    rough% disp,    ktau, alloc )
   call cable_resetDupVars ( dup% rough% hruff,   rough% hruff,   ktau, alloc )
   call cable_resetDupVars ( dup% rough% rt0us,   rough% rt0us,   ktau, alloc )
   call cable_resetDupVars ( dup% rough% rt1usa,  rough% rt1usa,  ktau, alloc )
   call cable_resetDupVars ( dup% rough% rt1usb,  rough% rt1usb,  ktau, alloc )
   call cable_resetDupVars ( dup% rough% rt1,     rough% rt1,     ktau, alloc )
   call cable_resetDupVars ( dup% rough% term2,   rough% term2,   ktau, alloc )
   call cable_resetDupVars ( dup% rough% term3,   rough% term3,   ktau, alloc )
   call cable_resetDupVars ( dup% rough% term5,   rough% term5,   ktau, alloc )
   call cable_resetDupVars ( dup% rough% term6,   rough% term6,   ktau, alloc )
   call cable_resetDupVars ( dup% rough% usuh,    rough% usuh,    ktau, alloc )
   call cable_resetDupVars ( dup% rough% za_uv,   rough% za_uv,   ktau, alloc )
   call cable_resetDupVars ( dup% rough% za_tq,   rough% za_tq,   ktau, alloc )
   call cable_resetDupVars ( dup% rough% z0m,     rough% z0m,     ktau, alloc )
   call cable_resetDupVars ( dup% rough% zref_uv, rough% zref_uv, ktau, alloc )
   call cable_resetDupVars ( dup% rough% zref_tq, rough% zref_tq, ktau, alloc )
   call cable_resetDupVars ( dup% rough% zruffs,  rough% zruffs,  ktau, alloc )
   call cable_resetDupVars ( dup% rough% z0soil,  rough% z0soil,  ktau, alloc )
   call cable_resetDupVars ( dup% rough% hruff_grmx, rough% hruff_grmx, ktau, alloc )
   call cable_resetDupVars ( dup% rough% z0soilsn,   rough% z0soilsn,   ktau, alloc )

   call cable_resetDupVars ( dup% air% rho,    air% rho,    ktau, alloc )
   call cable_resetDupVars ( dup% air% volm,   air% volm,   ktau, alloc )
   call cable_resetDupVars ( dup% air% rlam,   air% rlam,   ktau, alloc )
   call cable_resetDupVars ( dup% air% qsat,   air% qsat,   ktau, alloc )
   call cable_resetDupVars ( dup% air% epsi,   air% epsi,   ktau, alloc )
   call cable_resetDupVars ( dup% air% visc,   air% visc,   ktau, alloc )
   call cable_resetDupVars ( dup% air% psyc,   air% psyc,   ktau, alloc )
   call cable_resetDupVars ( dup% air% dsatdk, air% dsatdk, ktau, alloc )
   call cable_resetDupVars ( dup% air% cmolar, air% cmolar, ktau, alloc )

! ------------------------------------------------------------------------------

END SUBROUTINE reset_DupVars 

! ------------------------------------------------------------------------------





!SUBROUTINE alloc_balances_type(var, mp)
!   
!   TYPE(balances_type) :: var
!   INTEGER :: mp
!   
!   call cable_DupVars( var% drybal, mp, alloc ) 
!   call cable_DupVars( var% ebal, mp, alloc )  
!   call cable_DupVars( var% ebal_tot, mp, alloc )
!   call cable_DupVars( var% ebaltr, mp, alloc )  
!   call cable_DupVars( var% ebal_tottr, mp, alloc )
!   call cable_DupVars( var% ebal_cncheck, mp, alloc )  
!   call cable_DupVars( var% ebal_tot_cncheck, mp, alloc )
!   call cable_DupVars( var% evap_tot, mp, alloc ) 
!   call cable_DupVars( var% osnowd0, mp, alloc ) 
!   call cable_DupVars( var% precip_tot, mp, alloc  )
!   call cable_DupVars( var% rnoff_tot, mp, alloc ) 
!   call cable_DupVars( var% wbal, mp, alloc )    
!   call cable_DupVars( var% wbal_tot, mp, alloc  )
!   call cable_DupVars( var% wbtot0, mp, alloc )  
!   call cable_DupVars( var% wetbal, mp, alloc ) 
!   call cable_DupVars( var% cansto0, mp, alloc  ) 
!   call cable_DupVars( var% evapc_tot, mp, alloc ) 
!   call cable_DupVars( var% evaps_tot, mp, alloc ) 
!   call cable_DupVars( var% rnof1_tot, mp, alloc ) 
!   call cable_DupVars( var% rnof2_tot, mp, alloc ) 
!   call cable_DupVars( var% snowdc_tot, mp, alloc )
!   call cable_DupVars( var% wbal_tot1, mp, alloc ) 
!   call cable_DupVars( var% owbtot, mp, alloc )  
!   call cable_DupVars( var% delwc_tot, mp, alloc )  
!   call cable_DupVars( var% qasrf_tot, mp, alloc ) 
!   call cable_DupVars( var% qfsrf_tot, mp, alloc )  
!   call cable_DupVars( var% qssrf_tot, mp, alloc )  
!
!END SUBROUTINE alloc_balances_type
!
!! ------------------------------------------------------------------------------
!
!SUBROUTINE alloc_soil_parameter_type(var, mp)
!   
!   TYPE(soil_parameter_type) :: var
!   INTEGER :: mp
!   
!END SUBROUTINE alloc_soil_parameter_type
! 
!! ------------------------------------------------------------------------------
!
!SUBROUTINE alloc_soil_snow_type(var, mp)
!   
!   TYPE(soil_snow_type) :: var
!   INTEGER :: mp
!  
!END SUBROUTINE alloc_soil_snow_type
!
!! ------------------------------------------------------------------------------
!   
!SUBROUTINE alloc_veg_parameter_type(var, mp)
!
!   TYPE(veg_parameter_type) :: var
!   INTEGER :: mp
!   
!END SUBROUTINE alloc_veg_parameter_type
!
!! ------------------------------------------------------------------------------
!   
!SUBROUTINE alloc_canopy_type(var, mp)
!
!   TYPE(canopy_type) :: var
!   INTEGER :: mp
!END SUBROUTINE alloc_canopy_type
!
!  
!SUBROUTINE alloc_met_type(var, mp)
!
!   TYPE(met_type) :: var
!   INTEGER :: mp
! 
!END SUBROUTINE alloc_met_type
!   
! ------------------------------------------------------------------------------

!SUBROUTINE alloc_sum_flux_type(var, mp)
!
!   TYPE(sum_flux_type) :: var
!   INTEGER :: mp
! 
!   call cable_DupVars ( var % sumpn, mp, alloc )
!   call cable_DupVars ( var % sumrp, mp, alloc )
!   call cable_DupVars ( var % sumrpw, mp, alloc )
!   call cable_DupVars ( var % sumrpr, mp, alloc )
!   call cable_DupVars ( var % sumrs, mp, alloc )
!   call cable_DupVars ( var % sumrd, mp, alloc )
!   call cable_DupVars ( var % dsumpn, mp, alloc )
!   call cable_DupVars ( var % dsumrp, mp, alloc )
!   call cable_DupVars ( var % dsumrs, mp, alloc )
!   call cable_DupVars ( var % dsumrd, mp, alloc )
!   call cable_DupVars ( var % sumxrp, mp, alloc )
!   call cable_DupVars ( var % sumxrs, mp, alloc )
!
!END SUBROUTINE alloc_sum_flux_type

! ------------------------------------------------------------------------------
! Begin deallocation routines:
!SUBROUTINE dealloc_balances_type(var)
!   
!   TYPE(balances_type) :: var
!   
!   DEALLOCATE( var% drybal ) 
!   DEALLOCATE( var% ebal )  
!   DEALLOCATE( var% ebal_tot )
!   DEALLOCATE( var% ebaltr )  
!   DEALLOCATE( var% ebal_tottr )
!   DEALLOCATE( var% ebal_cncheck )  
!   DEALLOCATE( var% ebal_tot_cncheck )
!   DEALLOCATE( var% evap_tot)
!   DEALLOCATE( var% osnowd0 )
!   DEALLOCATE( var% precip_tot )
!   DEALLOCATE( var% rnoff_tot )
!   DEALLOCATE( var% wbal )   
!   DEALLOCATE( var% wbal_tot )
!   DEALLOCATE( var% wbtot0 ) 
!   DEALLOCATE( var% wetbal )
!   DEALLOCATE( var% cansto0 ) 
!   DEALLOCATE( var% evapc_tot ) 
!   DEALLOCATE( var% evaps_tot ) 
!   DEALLOCATE( var% rnof1_tot ) 
!   DEALLOCATE( var% rnof2_tot ) 
!   DEALLOCATE( var% snowdc_tot )
!   DEALLOCATE( var% wbal_tot1 ) 
!   DEALLOCATE( var% owbtot ) 
!   DEALLOCATE( var% delwc_tot ) 
!   DEALLOCATE( var% qasrf_tot )
!   DEALLOCATE( var% qfsrf_tot ) 
!   DEALLOCATE( var% qssrf_tot ) 
!   
!!END SUBROUTINE dealloc_balances_type
!
!! ------------------------------------------------------------------------------
!
!SUBROUTINE dealloc_soil_parameter_type(var)
!  
!   TYPE(soil_parameter_type) :: var
!   
!   DEALLOCATE( var% bch )   
!   DEALLOCATE( var% c3 )    
!   DEALLOCATE( var% clay )  
!   DEALLOCATE( var% css )   
!   DEALLOCATE( var% hsbh )  
!   DEALLOCATE( var% hyds )  
!   DEALLOCATE( var% i2bp3 ) 
!   DEALLOCATE( var% ibp2 )  
!   DEALLOCATE( var% isoilm )  
!   DEALLOCATE( var% rhosoil )  
!   DEALLOCATE( var% sand )   
!   DEALLOCATE( var% sfc )   
!   DEALLOCATE( var% silt )   
!   DEALLOCATE( var% ssat )   
!   DEALLOCATE( var% sucs )   
!   DEALLOCATE( var% swilt )  
!   DEALLOCATE( var% zse )    
!   DEALLOCATE( var% zshh )  
!   DEALLOCATE( var% cnsd )  
!   DEALLOCATE( var% albsoil )  
!   DEALLOCATE( var% cnsd )  
!   DEALLOCATE( var% pwb_min)  
!   DEALLOCATE( var% albsoilf )  
!   
!END SUBROUTINE dealloc_soil_parameter_type
! 
!! ------------------------------------------------------------------------------
!
!SUBROUTINE dealloc_soil_snow_type(var)
!   
!   TYPE(soil_snow_type) :: var
!  
!   DEALLOCATE ( var % iantrct )
!   DEALLOCATE ( var % pudsto )
!   DEALLOCATE ( var % pudsmx )
!   DEALLOCATE ( var % dtmlt )
!   DEALLOCATE( var% albsoilsn ) 
!   DEALLOCATE( var% cls )     
!   DEALLOCATE( var% dfn_dtg ) 
!   DEALLOCATE( var% dfh_dtg ) 
!   DEALLOCATE( var% dfe_ddq ) 
!   DEALLOCATE( var% ddq_dtg ) 
!   DEALLOCATE( var% evapsn )  
!   DEALLOCATE( var% fwtop )   
!   DEALLOCATE( var% fwtop1 )   
!   DEALLOCATE( var% fwtop2 )   
!   DEALLOCATE( var% fwtop3 )   
!   DEALLOCATE( var% gammzz ) 
!   DEALLOCATE( var% isflag ) 
!   DEALLOCATE( var% osnowd ) 
!   DEALLOCATE( var% potev ) 
!   DEALLOCATE( var% runoff )
!   DEALLOCATE( var% rnof1 ) 
!   DEALLOCATE( var% rnof2 ) 
!   DEALLOCATE( var% rtsoil )
!   DEALLOCATE( var% sconds ) 
!   DEALLOCATE( var% sdepth ) 
!   DEALLOCATE( var% smass ) 
!   DEALLOCATE( var% snage )  
!   DEALLOCATE( var% snowd )  
!   DEALLOCATE( var% smelt )  
!   DEALLOCATE( var% ssdn ) 
!   DEALLOCATE( var% ssdnn ) 
!   DEALLOCATE( var% tgg )   
!   DEALLOCATE( var% tggsn ) 
!   DEALLOCATE( var% tss )   
!   DEALLOCATE( var% tss_p )   
!   DEALLOCATE( var% deltss )   
!   DEALLOCATE( var% owb1 )   
!   DEALLOCATE( var% wb )    
!   DEALLOCATE( var% wbice ) 
!   DEALLOCATE( var% wblf ) 
!   DEALLOCATE( var%wbtot )    
!   DEALLOCATE( var%wbtot1 )    
!   DEALLOCATE( var%wbtot2 )    
!   DEALLOCATE( var%wb_lake )    
!   DEALLOCATE( var%sinfil )    
!   DEALLOCATE( var%evapfbl)    
!   DEALLOCATE( var%qstss)    
!   DEALLOCATE( var%wetfac )  
!   DEALLOCATE( var%owetfac )  
!   DEALLOCATE( var%t_snwlr )  
!   DEALLOCATE( var%wbfice )  
!   DEALLOCATE( var%tggav )  
!   DEALLOCATE( var%otgg )   
!   DEALLOCATE( var%otss )   
!   DEALLOCATE( var%otss_0 )   
!   DEALLOCATE( var%tprecip ) 
!   DEALLOCATE( var%tevap ) 
!   DEALLOCATE( var%trnoff ) 
!   DEALLOCATE( var%totenbal ) 
!   DEALLOCATE( var%totenbal2 ) 
!   DEALLOCATE( var%fland )      
!   DEALLOCATE( var%ifland )  
!   DEALLOCATE( var%tilefrac ) 
!   DEALLOCATE( var%qasrf )  
!   DEALLOCATE( var%qfsrf )  
!   DEALLOCATE( var%qssrf )  
!   
!END SUBROUTINE dealloc_soil_snow_type
!   
!! ------------------------------------------------------------------------------
!
!SUBROUTINE dealloc_veg_parameter_type(var)
!
!   TYPE(veg_parameter_type) :: var
!
!   DEALLOCATE( var% canst1 ) 
!   DEALLOCATE( var% dleaf )  
!   DEALLOCATE( var% ejmax ) 
!   DEALLOCATE( var% iveg ) 
!   DEALLOCATE( var% meth ) 
!   DEALLOCATE( var% frac4 )  
!   DEALLOCATE( var% hc )     
!   DEALLOCATE( var% vlai )   
!   DEALLOCATE( var% xalbnir ) 
!   DEALLOCATE( var% rp20 )   
!   DEALLOCATE( var% rpcoef ) 
!   DEALLOCATE( var% rs20 )   
!   DEALLOCATE( var% shelrb ) 
!   DEALLOCATE( var% vegcf )  
!   DEALLOCATE( var% tminvj ) 
!   DEALLOCATE( var% tmaxvj ) 
!   DEALLOCATE( var% vbeta)  
!   DEALLOCATE( var% vcmax )  
!   DEALLOCATE( var% xfang )  
!   DEALLOCATE( var%extkn ) 
!   DEALLOCATE( var%wai )   
!   DEALLOCATE( var%deciduous ) 
!   DEALLOCATE( var%froot) 
!   DEALLOCATE( var%refl )
!   DEALLOCATE( var%taul ) 
!   
!END SUBROUTINE dealloc_veg_parameter_type
!   
!! ------------------------------------------------------------------------------
!
!SUBROUTINE dealloc_canopy_type(var)
!
!   TYPE(canopy_type) :: var
!
!   DEALLOCATE ( var % fess )
!   DEALLOCATE ( var % fesp )
!   DEALLOCATE( var% cansto )  
!   DEALLOCATE( var% cduv )   
!   DEALLOCATE( var% delwc )  
!   DEALLOCATE( var% dewmm )  
!   DEALLOCATE( var% dgdtg )  
!   DEALLOCATE( var% fe )      
!   DEALLOCATE( var% fh )      
!   DEALLOCATE( var% fpn )     
!   DEALLOCATE( var% frp )     
!   DEALLOCATE( var% frpw )    
!   DEALLOCATE( var% frpr )    
!   DEALLOCATE( var% frs )     
!   DEALLOCATE( var% fnee )    
!   DEALLOCATE( var% frday )   
!   DEALLOCATE( var% fnv )     
!   DEALLOCATE( var% fev )     
!   DEALLOCATE( var% fevc )    
!   DEALLOCATE( var% fhv )     
!   DEALLOCATE( var% fns )     
!   DEALLOCATE( var% fhs )     
!   DEALLOCATE( var% fhs_cor )     
!   DEALLOCATE( var% ga )      
!   DEALLOCATE( var% ghflux )   
!   DEALLOCATE( var% precis ) 
!   DEALLOCATE( var% qscrn )  
!   DEALLOCATE( var% rnet )   
!   DEALLOCATE( var% segg )   
!   DEALLOCATE( var% sghflux )  
!   DEALLOCATE( var% through )  
!   DEALLOCATE( var% spill )  
!   DEALLOCATE( var% tscrn )  
!   DEALLOCATE( var% wcint )  
!   DEALLOCATE( var% tv )      
!   DEALLOCATE( var% us )      
!   DEALLOCATE( var% uscrn )   
!   DEALLOCATE( var% rghlai ) 
!   DEALLOCATE( var% vlaiw ) 
!   DEALLOCATE( var% fwet )   
!   DEALLOCATE ( var % evapfbl )
!   DEALLOCATE( var% epot )   
!   DEALLOCATE( var% fnpp )   
!   DEALLOCATE( var% fevw_pot )  
!   DEALLOCATE( var% gswx_T )  
!   DEALLOCATE( var% cdtq )   
!   DEALLOCATE( var% wetfac_cs )  
!   DEALLOCATE( var% fevw )   
!   DEALLOCATE( var% fhvw )   
!   DEALLOCATE( var% fes )    
!   DEALLOCATE( var% fes_cor )    
!   DEALLOCATE( var% gswx )  
!   DEALLOCATE( var% oldcansto )  
!   DEALLOCATE( var% zetar )  
!
!END SUBROUTINE dealloc_canopy_type
!   
!! ------------------------------------------------------------------------------
!
!SUBROUTINE dealloc_radiation_type(var)
!   
!   TYPE(radiation_type) :: var
!         
!   DEALLOCATE( var% albedo ) 
!   DEALLOCATE( var% extkb )  
!   DEALLOCATE( var% extkd2 )
!   DEALLOCATE( var% extkd )
!   DEALLOCATE( var% flws )
!   DEALLOCATE( var% fvlai )
!   DEALLOCATE( var% latitude )
!   DEALLOCATE( var% lwabv )
!   DEALLOCATE( var% qcan )
!   DEALLOCATE( var% qssabs )
!   DEALLOCATE( var% rhocdf )
!   DEALLOCATE( var% rniso )
!   DEALLOCATE( var% scalex )
!   DEALLOCATE( var% transd )
!   DEALLOCATE( var% trad )
!   DEALLOCATE( var% reffdf )
!   DEALLOCATE( var% reffbm )
!   DEALLOCATE( var% extkbm )
!   DEALLOCATE( var% extkdm )
!   DEALLOCATE( var% fbeam )
!   DEALLOCATE( var% cexpkbm )
!   DEALLOCATE( var% cexpkdm )
!   DEALLOCATE( var% rhocbm )
!   DEALLOCATE( var% transb )
!   DEALLOCATE( var% albedo_T )
!   DEALLOCATE( var% gradis )
!   DEALLOCATE( var% longitude )
!   DEALLOCATE( var% workp1 )
!   DEALLOCATE( var% workp2 )
!   DEALLOCATE( var% workp3 )
!   
!END SUBROUTINE dealloc_radiation_type
!   
!! ------------------------------------------------------------------------------
!
!SUBROUTINE dealloc_roughness_type(var)
!   
!   TYPE(roughness_type) :: var
!   
!   DEALLOCATE ( var % coexp )
!   DEALLOCATE ( var % disp )
!   DEALLOCATE ( var % hruff )
!   DEALLOCATE ( var % hruff_grmx )
!   DEALLOCATE ( var % rt0us )
!   DEALLOCATE ( var % rt1usa )
!   DEALLOCATE ( var % rt1usb )
!   DEALLOCATE ( var % rt1 )
!   DEALLOCATE ( var % term2 )
!   DEALLOCATE ( var % term3 )
!   DEALLOCATE ( var % term5 )
!   DEALLOCATE ( var % term6 )
!   DEALLOCATE ( var % usuh )
!   DEALLOCATE ( var % za_uv )
!   DEALLOCATE ( var % za_tq )
!   DEALLOCATE ( var % z0m )
!   DEALLOCATE ( var % zref_uv )
!   DEALLOCATE ( var % zref_tq )
!   DEALLOCATE ( var % zruffs )
!   DEALLOCATE ( var % z0soilsn )
!   DEALLOCATE ( var % z0soil )
!  
!END SUBROUTINE dealloc_roughness_type
!   
!! ------------------------------------------------------------------------------
!
!SUBROUTINE dealloc_air_type(var)
!   
!   TYPE(air_type) :: var
!   
!   DEALLOCATE ( var % rho )
!   DEALLOCATE ( var % volm )
!   DEALLOCATE ( var % rlam )
!   DEALLOCATE ( var % qsat )
!   DEALLOCATE ( var % epsi )
!   DEALLOCATE ( var % visc )
!   DEALLOCATE ( var % psyc )
!   DEALLOCATE ( var % dsatdk )
!   DEALLOCATE ( var % cmolar )
!  
!END SUBROUTINE dealloc_air_type
!   
!! ------------------------------------------------------------------------------
!
!SUBROUTINE dealloc_met_type(var)
!
!   TYPE(met_type) :: var
!   
!   DEALLOCATE ( var % ca )
!   DEALLOCATE ( var % year )
!   DEALLOCATE ( var % moy )
!   DEALLOCATE ( var % doy )
!   DEALLOCATE ( var % hod )
!   DEALLOCATE ( var % fsd )
!   DEALLOCATE ( var % ofsd )
!   DEALLOCATE ( var % fld )
!   DEALLOCATE ( var % precip )
!   DEALLOCATE ( var % precip_sn )
!   DEALLOCATE ( var % tk )
!   DEALLOCATE ( var % tvair )
!   DEALLOCATE ( var % tvrad )
!   DEALLOCATE ( var % pmb )
!   DEALLOCATE ( var % ua )
!   DEALLOCATE ( var % qv )
!   DEALLOCATE ( var % qvair )
!   DEALLOCATE ( var % da )
!   DEALLOCATE ( var % dva )
!   DEALLOCATE ( var % coszen )
!
!END SUBROUTINE dealloc_met_type
!
!! ------------------------------------------------------------------------------
!
!SUBROUTINE dealloc_sum_flux_type(var)
!
!   TYPE(sum_flux_type) :: var
!  
!   DEALLOCATE ( var % sumpn )
!   DEALLOCATE ( var % sumrp )
!   DEALLOCATE ( var % sumrpw )
!   DEALLOCATE ( var % sumrpr )
!   DEALLOCATE ( var % sumrs )
!   DEALLOCATE ( var % sumrd )
!   DEALLOCATE ( var % dsumpn )
!   DEALLOCATE ( var % dsumrp )
!   DEALLOCATE ( var % dsumrs )
!   DEALLOCATE ( var % dsumrd )
!   DEALLOCATE ( var % sumxrp )
!   DEALLOCATE ( var % sumxrs )
!
!END SUBROUTINE dealloc_sum_flux_type
!
!! ------------------------------------------------------------------------------
!
!SUBROUTINE dealloc_bgc_pool_type(var)
!   
!   TYPE(bgc_pool_type) :: var
!   
!   DEALLOCATE ( var % cplant )
!   DEALLOCATE ( var % csoil )
!
!END SUBROUTINE dealloc_bgc_pool_type
  
   !subroutine check_dup_var( a, b, ktau ) 
   subroutine check_dup_var( ) 
      real, dimension(10) :: a, b
      real, dimension(10) :: d
      integer :: ktau 
      a=1.
      b=1.
      ktau=1 
      d = a-b
      if (d  > 0.0 ) then
      !if (a - b .ne. 0. ) then
         print *, "ktau, dup% veg% ,  veg% ", ktau, a,b 
      endif

   end subroutine check_dup_var

end MODULE cable_Duplicate_types_mod
