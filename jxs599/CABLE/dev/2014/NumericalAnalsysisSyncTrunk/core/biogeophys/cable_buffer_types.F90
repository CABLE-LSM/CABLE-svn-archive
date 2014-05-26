
MODULE cable_buffer_types_mod
   
   use cable_BufVars_mod
   use cable_def_types_mod, only :  alloc, mp, mvtype, mstype, mland, & 
                                    r_2, n_tiles, ncp, ncs, mf, nrb, &
                                    msn, swb, niter, ms

   IMPLICIT NONE

   PUBLIC
  
! .............................................................................
   ! Soil parameters:
   TYPE Bsoil_parameter_type 
   
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

  END TYPE Bsoil_parameter_type

! .............................................................................

   ! Soil and snow variables:
   TYPE Bsoil_snow_type 
     
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

   END TYPE Bsoil_snow_type

! .............................................................................

   ! Vegetation parameters:
   TYPE Bveg_parameter_type
     
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

   END TYPE Bveg_parameter_type

! .............................................................................

   ! Canopy/vegetation variables:
   TYPE Bcanopy_type
      

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

   END TYPE Bcanopy_type

! .............................................................................

   ! Meterological data:
   TYPE Bmet_type
     
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
     
   END TYPE Bmet_type

! .............................................................................

   TYPE Bbgc_pool_type
      
      REAL, DIMENSION(:,:,:), allocatable ::                                         &
         cplant,  & ! plant carbon (g C/m2))
         csoil      ! soil carbon (g C/m2)
     
      REAL, DIMENSION(ncp)  :: ratecp ! plant carbon rate constant (1/year)
      
      REAL, DIMENSION(ncs)  :: ratecs ! soil carbon rate constant (1/year)
   
   END TYPE Bbgc_pool_type

! .............................................................................

   ! Radiation variables:
   TYPE Bradiation_type
   
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
    
    
  END TYPE Bradiation_type

 !.............................................................................

   ! Roughness variables:
   TYPE Broughness_type
      
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
   
   END TYPE Broughness_type

! .............................................................................

   ! Air variables:
   TYPE Bair_type
      
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

   END TYPE Bair_type

! .............................................................................

   type TbufVars
      ! CABLE variables
      TYPE (Bmet_type)       :: met     ! met input variables
      TYPE (Bair_type)       :: air     ! air property variables
      TYPE (Bcanopy_type)    :: canopy  ! vegetation variables
      TYPE (Bradiation_type) :: rad     ! radiation variables
      TYPE (Broughness_type) :: rough   ! roughness varibles
      TYPE (Bsoil_snow_type) :: ssnow   ! soil and snow variables
      !TYPE (Bbalances_type)  :: bal     ! energy and water balance variables
      
      ! CABLE parameters
      TYPE (Bsoil_parameter_type) :: soil ! soil parameters	
      TYPE (Bveg_parameter_type)  :: veg  ! vegetation parameters	 
      !TYPE (Bdriver_type)    :: C         ! constants used locally  
      
      !TYPE (Bsum_flux_type)  :: sum_flux ! cumulative flux variables
      TYPE (Bbgc_pool_type)  :: bgc  ! carbon pool variables
      
      !! CASA-CNP variables 
      !TYPE (Bcasa_biome)     :: casabiome
      !TYPE (Bcasa_pool)      :: casapool
      !TYPE (Bcasa_flux)      :: casaflux
      !TYPE (Bcasa_met)       :: casamet
      !TYPE (Bcasa_balance)   :: casabal
      !TYPE (Bphen_variable)  :: phen 
 
   end type TbufVars

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

CONTAINS

  
SUBROUTINE alloc_BufVars( buf, kend )
   integer :: kend
   type (TbufVars) :: buf

   call cable_BufVars( buf% bgc% cplant, mp, ncp, kend, alloc )
   call cable_BufVars( buf% bgc% csoil, mp, ncs, kend, alloc )

   call cable_BufVars( buf% soil% bch,    mp, kend, alloc )   
   call cable_BufVars( buf% soil% c3,     mp, kend, alloc )    
   call cable_BufVars( buf% soil% clay,   mp, kend, alloc )  
   call cable_BufVars( buf% soil% css,    mp, kend, alloc )   
   call cable_BufVars( buf% soil% hsbh,   mp, kend, alloc )  
   call cable_BufVars( buf% soil% hyds,   mp, kend, alloc )  
   call cable_BufVars( buf% soil% i2bp3,  mp, kend, alloc ) 
   call cable_BufVars( buf% soil% ibp2,   mp, kend, alloc )  
   call cable_BufVars( buf% soil% isoilm, mp, kend, alloc )  
   call cable_BufVars( buf% soil% rhosoil,mp, kend, alloc )  
   call cable_BufVars( buf% soil% sand,   mp, kend, alloc )   
   call cable_BufVars( buf% soil% sfc,    mp, kend, alloc )   
   call cable_BufVars( buf% soil% silt,   mp, kend, alloc )   
   call cable_BufVars( buf% soil% ssat,   mp, kend, alloc )   
   call cable_BufVars( buf% soil% sucs,   mp, kend, alloc )   
   call cable_BufVars( buf% soil% swilt,  mp, kend, alloc )  
   call cable_BufVars( buf% soil% zse,    ms, kend, alloc )    
   call cable_BufVars( buf% soil% zshh, ms+1, kend, alloc )  
   call cable_BufVars( buf% soil% cnsd,   mp, kend, alloc )  
   
   call cable_BufVars( buf% soil% albsoil,   mp, nrb, kend, alloc)   
   call cable_BufVars( buf% soil% pwb_min,   mp, kend, alloc )  
   call cable_BufVars( buf% soil% albsoilf,  mp, kend, alloc )  
   
   call cable_BufVars( buf% ssnow% iantrct,   mp, kend, alloc )
   call cable_BufVars( buf% ssnow % pudsto,   mp, kend, alloc )
   call cable_BufVars( buf% ssnow % pudsmx,   mp, kend, alloc )
   call cable_BufVars( buf% ssnow% cls,       mp, kend, alloc )     
   call cable_BufVars( buf% ssnow% dfn_dtg,   mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow% dfh_dtg,   mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow% dfe_ddq,   mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow% ddq_dtg,   mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow% evapsn,    mp, kend, alloc )  
   call cable_BufVars( buf% ssnow% fwtop,     mp, kend, alloc )   
   call cable_BufVars( buf% ssnow% fwtop1,    mp, kend, alloc )   
   call cable_BufVars( buf% ssnow% fwtop2,    mp, kend, alloc )   
   call cable_BufVars( buf% ssnow% fwtop3,    mp, kend, alloc )   
   call cable_BufVars( buf% ssnow% gammzz,    mp, ms, kend, alloc ) 
   call cable_BufVars( buf% ssnow% isflag,    mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow% osnowd,    mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow% potev,     mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow% runoff,    mp, kend, alloc )
   call cable_BufVars( buf% ssnow% rnof1,     mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow% rnof2,     mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow% rtsoil,    mp, kend, alloc )
   call cable_BufVars( buf% ssnow% sconds,    mp, msn, kend, alloc ) 
   call cable_BufVars( buf% ssnow% sdepth,    mp, msn, kend, alloc ) 
   call cable_BufVars( buf% ssnow% smass,     mp, msn, kend, alloc ) 
   call cable_BufVars( buf% ssnow% snage,     mp, kend, alloc )  
   call cable_BufVars( buf% ssnow% snowd,     mp, kend, alloc )  
   call cable_BufVars( buf% ssnow% smelt,     mp, kend, alloc )  
   call cable_BufVars( buf% ssnow% ssdn,      mp, msn, kend, alloc )  
   call cable_BufVars( buf% ssnow% ssdnn,     mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow% tgg,       mp, ms, kend, alloc )   
   call cable_BufVars( buf% ssnow% tggsn,     mp, msn, kend, alloc )  
   call cable_BufVars( buf% ssnow% tss,       mp, kend, alloc )   
   call cable_BufVars( buf% ssnow% tss_p,     mp, kend, alloc )   
   call cable_BufVars( buf% ssnow% deltss,    mp, kend, alloc )   
   call cable_BufVars( buf% ssnow% owb1,      mp, kend, alloc )   
   call cable_BufVars( buf% ssnow% wb,        mp, ms, kend, alloc )    
   call cable_BufVars( buf% ssnow% wbice,     mp, ms, kend, alloc ) 
   call cable_BufVars( buf% ssnow% wblf,      mp, ms, kend, alloc ) 
   call cable_BufVars( buf% ssnow%wbtot,      mp, kend, alloc )    
   call cable_BufVars( buf% ssnow%wbtot1,     mp, kend, alloc )    
   call cable_BufVars( buf% ssnow%wbtot2,     mp, kend, alloc )    
   call cable_BufVars( buf% ssnow%wb_lake,    mp, kend, alloc )    
   call cable_BufVars( buf% ssnow%sinfil,     mp, kend, alloc )    
   call cable_BufVars( buf% ssnow%evapfbl,    mp, ms, kend, alloc )    
   call cable_BufVars( buf% ssnow%qstss,      mp, kend, alloc )    
   call cable_BufVars( buf% ssnow%wetfac,     mp, kend, alloc )  
   call cable_BufVars( buf% ssnow%owetfac,    mp, kend, alloc )  
   call cable_BufVars( buf% ssnow%t_snwlr,    mp, kend, alloc )  
   call cable_BufVars( buf% ssnow%wbfice,     mp, ms, kend, alloc )  
   call cable_BufVars( buf% ssnow%tggav,      mp, kend, alloc )  
   call cable_BufVars( buf% ssnow%otgg,       mp, kend, alloc )   
   call cable_BufVars( buf% ssnow%otss,       mp, kend, alloc )   
   call cable_BufVars( buf% ssnow%otss_0,     mp, kend, alloc )   
   call cable_BufVars( buf% ssnow%tprecip,    mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow%tevap,      mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow%trnoff,     mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow%totenbal,   mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow%totenbal2,  mp, kend, alloc ) 
   call cable_BufVars( buf% ssnow%fland,      mp, kend, alloc )      
   call cable_BufVars( buf% ssnow%ifland,     mp, kend, alloc )  
   call cable_BufVars( buf% ssnow%tilefrac,   mp, n_tiles, kend, alloc ) 
   call cable_BufVars( buf% ssnow%qasrf,      mp, kend, alloc )  
   call cable_BufVars( buf% ssnow%qfsrf,      mp, kend, alloc )  
   call cable_BufVars( buf% ssnow%qssrf,      mp, kend, alloc )  

   call cable_BufVars( buf% ssnow % dtmlt,   mp, 3, kend, alloc )
   call cable_BufVars( buf% ssnow% albsoilsn,mp, nrb, kend, alloc ) 

   call cable_BufVars( buf% veg% canst1,   mp, kend, alloc ) 
   call cable_BufVars( buf% veg% dleaf,    mp, kend, alloc )  
   call cable_BufVars( buf% veg% ejmax,    mp, kend, alloc ) 
   call cable_BufVars( buf% veg% iveg,     mp, kend, alloc ) 
   call cable_BufVars( buf% veg% meth,     mp, kend, alloc ) 
   call cable_BufVars( buf% veg% frac4,    mp, kend, alloc )  
   call cable_BufVars( buf% veg% hc,       mp, kend, alloc )     
   call cable_BufVars( buf% veg% vlai,     mp, kend, alloc )   
   call cable_BufVars( buf% veg% xalbnir,  mp, kend, alloc ) 
   call cable_BufVars( buf% veg% rp20,     mp, kend, alloc )   
   call cable_BufVars( buf% veg% rpcoef,   mp, kend, alloc ) 
   call cable_BufVars( buf% veg% rs20,     mp, kend, alloc )   
   call cable_BufVars( buf% veg% shelrb,   mp, kend, alloc ) 
   call cable_BufVars( buf% veg% vegcf,    mp, kend, alloc )  
   call cable_BufVars( buf% veg% tminvj,   mp, kend, alloc ) 
   call cable_BufVars( buf% veg% tmaxvj,   mp, kend, alloc ) 
   call cable_BufVars( buf% veg% vbeta,    mp, kend, alloc )  
   call cable_BufVars( buf% veg% vcmax,    mp, kend, alloc )  
   call cable_BufVars( buf% veg% xfang,    mp, kend, alloc )  
   call cable_BufVars( buf% veg% extkn,     mp, kend, alloc ) 
   call cable_BufVars( buf% veg% wai,       mp, kend, alloc )   
   call cable_BufVars( buf% veg% deciduous, mp, kend, alloc ) 
   call cable_BufVars( buf% veg% froot,     mp, ms, kend, alloc ) 
   call cable_BufVars( buf% veg% refl,      mp, 2, kend, alloc ) !jhan:swb?
   call cable_BufVars( buf% veg% taul,      mp, 2, kend, alloc ) 
   call cable_BufVars( buf% veg% vlaimax,   mp, kend, alloc ) 

   call cable_BufVars( buf% canopy% fess,    mp, kend, alloc )
   call cable_BufVars( buf% canopy% fesp,    mp, kend, alloc )
   call cable_BufVars( buf% canopy% cansto,  mp, kend, alloc )  
   call cable_BufVars( buf% canopy% cduv,    mp, kend, alloc )   
   call cable_BufVars( buf% canopy% delwc,   mp, kend, alloc )  
   call cable_BufVars( buf% canopy% dewmm,   mp, kend, alloc )  
   call cable_BufVars( buf% canopy% dgdtg,   mp, kend, alloc )  
   call cable_BufVars( buf% canopy% fe,      mp, kend, alloc )      
   call cable_BufVars( buf% canopy% fh,      mp, kend, alloc )      
   call cable_BufVars( buf% canopy% fpn,     mp, kend, alloc )     
   call cable_BufVars( buf% canopy% frp,     mp, kend, alloc )     
   call cable_BufVars( buf% canopy% frpw,    mp, kend, alloc )    
   call cable_BufVars( buf% canopy% frpr,    mp, kend, alloc )    
   call cable_BufVars( buf% canopy% frs,     mp, kend, alloc )     
   call cable_BufVars( buf% canopy% fnee,    mp, kend, alloc )    
   call cable_BufVars( buf% canopy% frday,   mp, kend, alloc )   
   call cable_BufVars( buf% canopy% fnv,     mp, kend, alloc )     
   call cable_BufVars( buf% canopy% fev,     mp, kend, alloc )     
   call cable_BufVars( buf% canopy% fevc,    mp, kend, alloc )    
   call cable_BufVars( buf% canopy% fhv,     mp, kend, alloc )     
   call cable_BufVars( buf% canopy% fns,     mp, kend, alloc )     
   call cable_BufVars( buf% canopy% fhs,     mp, kend, alloc )     
   call cable_BufVars( buf% canopy% fhs_cor, mp, kend, alloc )     
   call cable_BufVars( buf% canopy% ga,      mp, kend, alloc )      
   call cable_BufVars( buf% canopy% ghflux,  mp, kend, alloc )   
   call cable_BufVars( buf% canopy% precis,  mp, kend, alloc ) 
   call cable_BufVars( buf% canopy% qscrn,   mp, kend, alloc )  
   call cable_BufVars( buf% canopy% rnet,    mp, kend, alloc )   
   call cable_BufVars( buf% canopy% segg,    mp, kend, alloc )   
   call cable_BufVars( buf% canopy% sghflux, mp, kend, alloc )  
   call cable_BufVars( buf% canopy% through, mp, kend, alloc )  
   call cable_BufVars( buf% canopy% spill,   mp, kend, alloc )  
   call cable_BufVars( buf% canopy% tscrn,   mp, kend, alloc )  
   call cable_BufVars( buf% canopy% wcint,   mp, kend, alloc )  
   call cable_BufVars( buf% canopy% tv,      mp, kend, alloc )      
   call cable_BufVars( buf% canopy% us,      mp, kend, alloc )      
   call cable_BufVars( buf% canopy% uscrn,   mp, kend, alloc )   
   call cable_BufVars( buf% canopy% rghlai,  mp, kend, alloc ) 
   call cable_BufVars( buf% canopy% vlaiw,   mp, kend, alloc ) 
   call cable_BufVars( buf% canopy% fwet,    mp, kend, alloc )   
   call cable_BufVars( buf% canopy% epot,    mp, kend, alloc )   
   call cable_BufVars( buf% canopy% fnpp,    mp, kend, alloc )   
   call cable_BufVars( buf% canopy% fevw_pot,mp, kend, alloc )  
   call cable_BufVars( buf% canopy% gswx_T,  mp, kend, alloc )  
   call cable_BufVars( buf% canopy% cdtq,    mp, kend, alloc )   
   call cable_BufVars( buf% canopy% fevw,    mp, kend, alloc )   
   call cable_BufVars( buf% canopy% fhvw,    mp, kend, alloc )   
   call cable_BufVars( buf% canopy% fes,     mp, kend, alloc )    
   call cable_BufVars( buf% canopy% fes_cor, mp, kend, alloc )    
   call cable_BufVars( buf% canopy% gswx,    mp, mf, kend, alloc )  
   call cable_BufVars( buf% canopy% evapfbl, mp, ms, kend, alloc )
   call cable_BufVars( buf% canopy% wetfac_cs,  mp, kend, alloc )  
   call cable_BufVars( buf% canopy% oldcansto,  mp, kend, alloc )  
   call cable_BufVars( buf% canopy% zetar,   mp, NITER, kend, alloc )  
 
   call cable_BufVars ( buf% met% ca,     mp, kend, alloc )
   call cable_BufVars ( buf% met% year,   mp, kend, alloc )
   call cable_BufVars ( buf% met% moy,    mp, kend, alloc )
   call cable_BufVars ( buf% met% doy,    mp, kend, alloc )
   call cable_BufVars ( buf% met% hod,    mp, kend, alloc )
   call cable_BufVars ( buf% met% ofsd,   mp, kend, alloc ) 
   call cable_BufVars ( buf% met% fld,    mp, kend, alloc )
   call cable_BufVars ( buf% met% precip, mp, kend, alloc )
   call cable_BufVars ( buf% met% tk,     mp, kend, alloc )
   call cable_BufVars ( buf% met% tvair,  mp, kend, alloc )
   call cable_BufVars ( buf% met% tvrad,  mp, kend, alloc )
   call cable_BufVars ( buf% met% pmb,    mp, kend, alloc )
   call cable_BufVars ( buf% met% ua,     mp, kend, alloc )
   call cable_BufVars ( buf% met% qv,     mp, kend, alloc )
   call cable_BufVars ( buf% met% qvair,  mp, kend, alloc )
   call cable_BufVars ( buf% met% da,     mp, kend, alloc )
   call cable_BufVars ( buf% met% dva,    mp, kend, alloc )
   call cable_BufVars ( buf% met% coszen, mp, kend, alloc )
   call cable_BufVars ( buf% met% fsd,    mp, swb, kend, alloc ) 
   call cable_BufVars ( buf% met% precip_sn, mp, kend, alloc )
   
   call cable_BufVars( buf% rad% extkb,   mp, kend, alloc )  
   call cable_BufVars( buf% rad% extkd2,  mp, kend, alloc )
   call cable_BufVars( buf% rad% extkd,   mp, kend, alloc )
   call cable_BufVars( buf% rad% flws,    mp, kend, alloc )
   call cable_BufVars( buf% rad% lwabv,   mp, kend, alloc )
   call cable_BufVars( buf% rad% qssabs,  mp, kend, alloc )
   call cable_BufVars( buf% rad% transd,  mp, kend, alloc )
   call cable_BufVars( buf% rad% trad,    mp, kend, alloc )
   call cable_BufVars( buf% rad% transb,  mp, kend, alloc )
   call cable_BufVars( buf% rad% workp1,  mp, kend, alloc )
   call cable_BufVars( buf% rad% workp2,  mp, kend, alloc )
   call cable_BufVars( buf% rad% workp3,  mp, kend, alloc )
   call cable_BufVars( buf% rad% albedo_T,   mp,   kend, alloc )
   call cable_BufVars( buf% rad% latitude,   mp,   kend, alloc )
   call cable_BufVars( buf% rad% longitude,  mp,   kend, alloc )
   call cable_BufVars( buf% rad% fvlai,   mp, mf,  kend, alloc )
   call cable_BufVars( buf% rad% rniso,   mp, mf,  kend, alloc )
   call cable_BufVars( buf% rad% scalex,  mp, mf,  kend, alloc )
   call cable_BufVars( buf% rad% gradis,  mp, mf,  kend, alloc )
   call cable_BufVars( buf% rad% albedo,  mp, nrb, kend, alloc ) 
   call cable_BufVars( buf% rad% rhocdf,  mp, nrb, kend, alloc )
   call cable_BufVars( buf% rad% reffdf,  mp, nrb, kend, alloc )
   call cable_BufVars( buf% rad% reffbm,  mp, nrb, kend, alloc )
   call cable_BufVars( buf% rad% extkbm,  mp, nrb, kend, alloc )
   call cable_BufVars( buf% rad% extkdm,  mp, nrb, kend, alloc )
   call cable_BufVars( buf% rad% cexpkbm, mp, swb, kend, alloc )
   call cable_BufVars( buf% rad% cexpkdm, mp, swb, kend, alloc )
   call cable_BufVars( buf% rad% fbeam,   mp, nrb, kend, alloc )
   call cable_BufVars( buf% rad% rhocbm,  mp, nrb, kend, alloc )
   call cable_BufVars( buf% rad% qcan,    mp, mf, nrb, kend, alloc )

   call cable_BufVars ( buf% rough% coexp,   mp, kend, alloc )
   call cable_BufVars ( buf% rough% disp,    mp, kend, alloc )
   call cable_BufVars ( buf% rough% hruff,   mp, kend, alloc )
   call cable_BufVars ( buf% rough% rt0us,   mp, kend, alloc )
   call cable_BufVars ( buf% rough% rt1usa,  mp, kend, alloc )
   call cable_BufVars ( buf% rough% rt1usb,  mp, kend, alloc )
   call cable_BufVars ( buf% rough% rt1,     mp, kend, alloc )
   call cable_BufVars ( buf% rough% term2,   mp, kend, alloc )
   call cable_BufVars ( buf% rough% term3,   mp, kend, alloc )
   call cable_BufVars ( buf% rough% term5,   mp, kend, alloc )
   call cable_BufVars ( buf% rough% term6,   mp, kend, alloc )
   call cable_BufVars ( buf% rough% usuh,    mp, kend, alloc )
   call cable_BufVars ( buf% rough% za_uv,   mp, kend, alloc )
   call cable_BufVars ( buf% rough% za_tq,   mp, kend, alloc )
   call cable_BufVars ( buf% rough% z0m,     mp, kend, alloc )
   call cable_BufVars ( buf% rough% zref_uv, mp, kend, alloc )
   call cable_BufVars ( buf% rough% zref_tq, mp, kend, alloc )
   call cable_BufVars ( buf% rough% zruffs,  mp, kend, alloc )
   call cable_BufVars ( buf% rough% z0soil,  mp, kend, alloc )
   call cable_BufVars ( buf% rough% hruff_grmx, mp, kend, alloc )
   call cable_BufVars ( buf% rough% z0soilsn,   mp, kend, alloc )

   call cable_BufVars ( buf% air% rho,    mp, kend, alloc )
   call cable_BufVars ( buf% air% volm,   mp, kend, alloc )
   call cable_BufVars ( buf% air% rlam,   mp, kend, alloc )
   call cable_BufVars ( buf% air% qsat,   mp, kend, alloc )
   call cable_BufVars ( buf% air% epsi,   mp, kend, alloc )
   call cable_BufVars ( buf% air% visc,   mp, kend, alloc )
   call cable_BufVars ( buf% air% psyc,   mp, kend, alloc )
   call cable_BufVars ( buf% air% dsatdk, mp, kend, alloc )
   call cable_BufVars ( buf% air% cmolar, mp, kend, alloc )

END SUBROUTINE alloc_BufVars 

! ------------------------------------------------------------------------------
  
subroutine Set_BufVars( ktau, buf, met, bgc, soil,veg, ssnow, canopy, &
                        rad, rough, air )

   use cable_def_types_mod, only : met_type, bgc_pool_type, soil_parameter_type, &
                                   veg_parameter_type, soil_snow_type, canopy_type, & 
                                   air_type, radiation_type, roughness_type

   use cable_Bufvars_mod, only : cable_setBufVars  
   
   integer :: ktau
   type (TbufVars) :: buf
 
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
 
   call cable_setBufVars( buf% bgc% cplant, bgc% cplant,   ktau, alloc )
   call cable_setBufVars( buf% bgc% csoil,  bgc% csoil,    ktau, alloc )
              
   call cable_setBufVars( buf% soil% bch,   soil% bch,     ktau, alloc )   
   call cable_setBufVars( buf% soil% c3,    soil% c3,      ktau, alloc )    
   call cable_setBufVars( buf% soil% clay,  soil% clay,    ktau, alloc )  
   call cable_setBufVars( buf% soil% css,   soil% css,     ktau, alloc )   
   call cable_setBufVars( buf% soil% hsbh,  soil% hsbh,    ktau, alloc )  
   call cable_setBufVars( buf% soil% hyds,  soil% hyds,    ktau, alloc )  
   call cable_setBufVars( buf% soil% i2bp3, soil% i2bp3,   ktau, alloc ) 
   call cable_setBufVars( buf% soil% ibp2,  soil% ibp2,    ktau, alloc )  
   call cable_setBufVars( buf% soil% sand,  soil% sand,    ktau, alloc )   
   call cable_setBufVars( buf% soil% sfc,   soil% sfc,     ktau, alloc )   
   call cable_setBufVars( buf% soil% silt,  soil% silt,    ktau, alloc )   
   call cable_setBufVars( buf% soil% ssat,  soil% ssat,    ktau, alloc )   
   call cable_setBufVars( buf% soil% sucs,  soil% sucs,    ktau, alloc )   
   call cable_setBufVars( buf% soil% swilt, soil% swilt,   ktau, alloc )  
   call cable_setBufVars( buf% soil% zse,   soil% zse,     ktau, alloc )    
   call cable_setBufVars( buf% soil% zshh, soil% zshh,     ktau, alloc )  
   call cable_setBufVars( buf% soil% cnsd, soil% cnsd,     ktau, alloc )  
   call cable_setBufVars( buf% soil% isoilm,  soil% isoilm,  ktau, alloc )  
   call cable_setBufVars( buf% soil% rhosoil, soil% rhosoil, ktau, alloc )  
   
   call cable_setBufVars( buf% soil% albsoil,  soil% albsoil,  ktau, alloc)   
   call cable_setBufVars( buf% soil% pwb_min,  soil% pwb_min,  ktau, alloc )  
   call cable_setBufVars( buf% soil% albsoilf, soil% albsoilf, ktau, alloc )  
   
   call cable_setBufVars( buf% ssnow% iantrct,   ssnow% iantrct,   ktau, alloc )
   call cable_setBufVars( buf% ssnow % pudsto,   ssnow % pudsto,   ktau, alloc )
   call cable_setBufVars( buf% ssnow % pudsmx,   ssnow % pudsmx,   ktau, alloc )
   call cable_setBufVars( buf% ssnow% cls,       ssnow% cls,       ktau, alloc )     
   call cable_setBufVars( buf% ssnow% dfn_dtg,   ssnow% dfn_dtg,   ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% dfh_dtg,   ssnow% dfh_dtg,   ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% dfe_ddq,   ssnow% dfe_ddq,   ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% ddq_dtg,   ssnow% ddq_dtg,   ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% evapsn,    ssnow% evapsn,    ktau, alloc )  
   call cable_setBufVars( buf% ssnow% fwtop,     ssnow% fwtop,     ktau, alloc )   
   call cable_setBufVars( buf% ssnow% fwtop1,    ssnow% fwtop1,    ktau, alloc )   
   call cable_setBufVars( buf% ssnow% fwtop2,    ssnow% fwtop2,    ktau, alloc )   
   call cable_setBufVars( buf% ssnow% fwtop3,    ssnow% fwtop3,    ktau, alloc )   
   call cable_setBufVars( buf% ssnow% gammzz,    ssnow% gammzz,    ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% isflag,    ssnow% isflag,    ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% osnowd,    ssnow% osnowd,    ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% potev,     ssnow% potev,     ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% runoff,    ssnow% runoff,    ktau, alloc )
   call cable_setBufVars( buf% ssnow% rnof1,     ssnow% rnof1,     ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% rnof2,     ssnow% rnof2,     ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% rtsoil,    ssnow% rtsoil,    ktau, alloc )
   call cable_setBufVars( buf% ssnow% sconds,    ssnow% sconds,    ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% sdepth,    ssnow% sdepth,    ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% smass,     ssnow% smass,     ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% snage,     ssnow% snage,     ktau, alloc )  
   call cable_setBufVars( buf% ssnow% snowd,     ssnow% snowd,     ktau, alloc )  
   call cable_setBufVars( buf% ssnow% smelt,     ssnow% smelt,     ktau, alloc )  
   call cable_setBufVars( buf% ssnow% ssdn,      ssnow% ssdn,      ktau, alloc )  
   call cable_setBufVars( buf% ssnow% ssdnn,     ssnow% ssdnn,     ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% tgg,       ssnow% tgg,       ktau, alloc )   
   call cable_setBufVars( buf% ssnow% tggsn,     ssnow% tggsn,     ktau, alloc )  
   call cable_setBufVars( buf% ssnow% tss,       ssnow% tss,       ktau, alloc )   
   call cable_setBufVars( buf% ssnow% tss_p,     ssnow% tss_p,     ktau, alloc )   
   call cable_setBufVars( buf% ssnow% deltss,    ssnow% deltss,    ktau, alloc )   
   call cable_setBufVars( buf% ssnow% owb1,      ssnow% owb1,      ktau, alloc )   
   call cable_setBufVars( buf% ssnow% wb,        ssnow% wb,        ktau, alloc )    
   call cable_setBufVars( buf% ssnow% wbice,     ssnow% wbice,     ktau, alloc ) 
   call cable_setBufVars( buf% ssnow% wblf,      ssnow% wblf,      ktau, alloc ) 
   call cable_setBufVars( buf% ssnow%wbtot,      ssnow%wbtot,      ktau, alloc )    
   call cable_setBufVars( buf% ssnow%wbtot1,     ssnow%wbtot1,     ktau, alloc )    
   call cable_setBufVars( buf% ssnow%wbtot2,     ssnow%wbtot2,     ktau, alloc )    
   call cable_setBufVars( buf% ssnow%wb_lake,    ssnow%wb_lake,    ktau, alloc )    
   call cable_setBufVars( buf% ssnow%sinfil,     ssnow%sinfil,     ktau, alloc )    
   call cable_setBufVars( buf% ssnow%evapfbl,    ssnow%evapfbl,    ktau, alloc )    
   call cable_setBufVars( buf% ssnow%qstss,      ssnow%qstss,      ktau, alloc )    
   call cable_setBufVars( buf% ssnow%wetfac,     ssnow%wetfac,     ktau, alloc )  
   call cable_setBufVars( buf% ssnow%owetfac,    ssnow%owetfac,    ktau, alloc )  
   call cable_setBufVars( buf% ssnow%t_snwlr,    ssnow%t_snwlr,    ktau, alloc )  
   call cable_setBufVars( buf% ssnow%wbfice,     ssnow%wbfice,     ktau, alloc )  
   call cable_setBufVars( buf% ssnow%tggav,      ssnow%tggav,      ktau, alloc )  
   call cable_setBufVars( buf% ssnow%otgg,       ssnow%otgg,       ktau, alloc )   
   call cable_setBufVars( buf% ssnow%otss,       ssnow%otss,       ktau, alloc )   
   call cable_setBufVars( buf% ssnow%otss_0,     ssnow%otss_0,     ktau, alloc )   
   call cable_setBufVars( buf% ssnow%tprecip,    ssnow%tprecip,    ktau, alloc ) 
   call cable_setBufVars( buf% ssnow%tevap,      ssnow%tevap,      ktau, alloc ) 
   call cable_setBufVars( buf% ssnow%trnoff,     ssnow%trnoff,     ktau, alloc ) 
   call cable_setBufVars( buf% ssnow%totenbal,   ssnow%totenbal,   ktau, alloc ) 
   call cable_setBufVars( buf% ssnow%totenbal2,  ssnow%totenbal2,  ktau, alloc ) 
   call cable_setBufVars( buf% ssnow%fland,      ssnow%fland,      ktau, alloc )      
   call cable_setBufVars( buf% ssnow%ifland,     ssnow%ifland,     ktau, alloc )  
   call cable_setBufVars( buf% ssnow%tilefrac,   ssnow%tilefrac,   ktau, alloc ) 
   call cable_setBufVars( buf% ssnow%qasrf,      ssnow%qasrf,      ktau, alloc )  
   call cable_setBufVars( buf% ssnow%qfsrf,      ssnow%qfsrf,      ktau, alloc )  
   call cable_setBufVars( buf% ssnow%qssrf,      ssnow%qssrf,      ktau, alloc )  

   call cable_setBufVars( buf% ssnow % dtmlt,    ssnow % dtmlt,     ktau, alloc )
   call cable_setBufVars( buf% ssnow% albsoilsn, ssnow% albsoilsn,  ktau, alloc ) 

   call cable_setBufVars( buf% veg% canst1,   veg% canst1,   ktau, alloc ) 
   call cable_setBufVars( buf% veg% dleaf,    veg% dleaf,    ktau, alloc )  
   call cable_setBufVars( buf% veg% ejmax,    veg% ejmax,    ktau, alloc ) 
   call cable_setBufVars( buf% veg% iveg,     veg% iveg,     ktau, alloc ) 
   call cable_setBufVars( buf% veg% meth,     veg% meth,     ktau, alloc ) 
   call cable_setBufVars( buf% veg% frac4,    veg% frac4,    ktau, alloc )  
   call cable_setBufVars( buf% veg% hc,       veg% hc,       ktau, alloc )     
   call cable_setBufVars( buf% veg% vlai,     veg% vlai,     ktau, alloc )   
   call cable_setBufVars( buf% veg% xalbnir,  veg% xalbnir,  ktau, alloc ) 
   call cable_setBufVars( buf% veg% rp20,     veg% rp20,     ktau, alloc )   
   call cable_setBufVars( buf% veg% rpcoef,   veg% rpcoef,   ktau, alloc ) 
   call cable_setBufVars( buf% veg% rs20,     veg% rs20,     ktau, alloc )   
   call cable_setBufVars( buf% veg% shelrb,   veg% shelrb,   ktau, alloc ) 
   call cable_setBufVars( buf% veg% vegcf,    veg% vegcf,    ktau, alloc )  
   call cable_setBufVars( buf% veg% tminvj,   veg% tminvj,   ktau, alloc ) 
   call cable_setBufVars( buf% veg% tmaxvj,   veg% tmaxvj,   ktau, alloc ) 
   call cable_setBufVars( buf% veg% vbeta,    veg% vbeta,    ktau, alloc )  
   call cable_setBufVars( buf% veg% vcmax,    veg% vcmax,    ktau, alloc )  
   call cable_setBufVars( buf% veg% xfang,    veg% xfang,    ktau, alloc )  
   call cable_setBufVars( buf% veg% extkn,    veg% extkn,     ktau, alloc ) 
   call cable_setBufVars( buf% veg% wai,      veg% wai,       ktau, alloc )   
   call cable_setBufVars( buf% veg% deciduous,veg% deciduous, ktau, alloc ) 
   call cable_setBufVars( buf% veg% froot,    veg% froot,     ktau, alloc ) 
   call cable_setBufVars( buf% veg% refl,     veg% refl,      ktau, alloc ) !jhan:swb?
   call cable_setBufVars( buf% veg% taul,     veg% taul,      ktau, alloc ) 
   call cable_setBufVars( buf% veg% vlaimax,  veg% vlaimax,   ktau, alloc ) 

   call cable_setBufVars( buf% canopy% fess,   canopy% fess,   ktau, alloc )
   call cable_setBufVars( buf% canopy% fesp,   canopy% fesp,   ktau, alloc )
   call cable_setBufVars( buf% canopy% cansto, canopy% cansto, ktau, alloc )  
   call cable_setBufVars( buf% canopy% cduv,   canopy% cduv,   ktau, alloc )   
   call cable_setBufVars( buf% canopy% delwc,  canopy% delwc,  ktau, alloc )  
   call cable_setBufVars( buf% canopy% dewmm,  canopy% dewmm,  ktau, alloc )  
   call cable_setBufVars( buf% canopy% dgdtg,  canopy% dgdtg,  ktau, alloc )  
   call cable_setBufVars( buf% canopy% fe,     canopy% fe,     ktau, alloc )      
   call cable_setBufVars( buf% canopy% fh,     canopy% fh,     ktau, alloc )      
   call cable_setBufVars( buf% canopy% fpn,    canopy% fpn,    ktau, alloc )     
   call cable_setBufVars( buf% canopy% frp,    canopy% frp,    ktau, alloc )     
   call cable_setBufVars( buf% canopy% frpw,   canopy% frpw,   ktau, alloc )    
   call cable_setBufVars( buf% canopy% frpr,   canopy% frpr,   ktau, alloc )    
   call cable_setBufVars( buf% canopy% frs,    canopy% frs,    ktau, alloc )     
   call cable_setBufVars( buf% canopy% fnee,   canopy% fnee,   ktau, alloc )    
   call cable_setBufVars( buf% canopy% frday,  canopy% frday,  ktau, alloc )   
   call cable_setBufVars( buf% canopy% fnv,    canopy% fnv,    ktau, alloc )     
   call cable_setBufVars( buf% canopy% fev,    canopy% fev,    ktau, alloc )     
   call cable_setBufVars( buf% canopy% fevc,   canopy% fevc,   ktau, alloc )    
   call cable_setBufVars( buf% canopy% fhv,    canopy% fhv,    ktau, alloc )     
   call cable_setBufVars( buf% canopy% fns,    canopy% fns,    ktau, alloc )     
   call cable_setBufVars( buf% canopy% fhs,    canopy% fhs,    ktau, alloc )     
   call cable_setBufVars( buf% canopy% fhs_cor,canopy% fhs_cor,ktau, alloc )     
   call cable_setBufVars( buf% canopy% ga,     canopy% ga,     ktau, alloc )      
   call cable_setBufVars( buf% canopy% ghflux, canopy% ghflux, ktau, alloc )   
   call cable_setBufVars( buf% canopy% precis, canopy% precis, ktau, alloc ) 
   call cable_setBufVars( buf% canopy% qscrn,  canopy% qscrn,  ktau, alloc )  
   call cable_setBufVars( buf% canopy% rnet,   canopy% rnet,   ktau, alloc )   
   call cable_setBufVars( buf% canopy% segg,   canopy% segg,   ktau, alloc )   
   call cable_setBufVars( buf% canopy% sghflux,canopy% sghflux,ktau, alloc )  
   call cable_setBufVars( buf% canopy% through,canopy% through,ktau, alloc )  
   call cable_setBufVars( buf% canopy% spill,  canopy% spill,  ktau, alloc )  
   call cable_setBufVars( buf% canopy% tscrn,  canopy% tscrn,  ktau, alloc )  
   call cable_setBufVars( buf% canopy% wcint,  canopy% wcint,  ktau, alloc )  
   call cable_setBufVars( buf% canopy% tv,     canopy% tv,     ktau, alloc )      
   call cable_setBufVars( buf% canopy% us,     canopy% us,     ktau, alloc )      
   call cable_setBufVars( buf% canopy% uscrn,  canopy% uscrn,  ktau, alloc )   
   call cable_setBufVars( buf% canopy% rghlai, canopy% rghlai, ktau, alloc ) 
   call cable_setBufVars( buf% canopy% vlaiw,  canopy% vlaiw,  ktau, alloc ) 
   call cable_setBufVars( buf% canopy% fwet,   canopy% fwet,   ktau, alloc )   
   call cable_setBufVars( buf% canopy% epot,   canopy% epot,   ktau, alloc )   
   call cable_setBufVars( buf% canopy% fnpp,   canopy% fnpp,   ktau, alloc )   
   call cable_setBufVars( buf% canopy% gswx_T, canopy% gswx_T, ktau, alloc )  
   call cable_setBufVars( buf% canopy% cdtq,   canopy% cdtq,   ktau, alloc )   
   call cable_setBufVars( buf% canopy% fevw,   canopy% fevw,   ktau, alloc )   
   call cable_setBufVars( buf% canopy% fhvw,   canopy% fhvw,   ktau, alloc )   
   call cable_setBufVars( buf% canopy% fes,    canopy% fes,    ktau, alloc )    
   call cable_setBufVars( buf% canopy% fes_cor,canopy% fes_cor,ktau, alloc )    
   call cable_setBufVars( buf% canopy% gswx,   canopy% gswx,   ktau, alloc )  
   call cable_setBufVars( buf% canopy% evapfbl,canopy% evapfbl,ktau, alloc )
   call cable_setBufVars( buf% canopy% wetfac_cs, canopy% wetfac_cs, ktau, alloc )  
   call cable_setBufVars( buf% canopy% oldcansto, canopy% oldcansto, ktau, alloc )  
   call cable_setBufVars( buf% canopy% zetar,     canopy% zetar,     ktau, alloc )  
   call cable_setBufVars( buf% canopy% fevw_pot,  canopy% fevw_pot,  ktau, alloc )  
 
   call cable_setBufVars ( buf% met% ca,    met% ca,     ktau, alloc )
   call cable_setBufVars ( buf% met% year,  met% year,   ktau, alloc )
   call cable_setBufVars ( buf% met% moy,   met% moy,    ktau, alloc )
   call cable_setBufVars ( buf% met% doy,   met% doy,    ktau, alloc )
   call cable_setBufVars ( buf% met% hod,   met% hod,    ktau, alloc )
   call cable_setBufVars ( buf% met% ofsd,  met% ofsd,   ktau, alloc ) 
   call cable_setBufVars ( buf% met% fld,   met% fld,    ktau, alloc )
   call cable_setBufVars ( buf% met% precip,met% precip, ktau, alloc )
   call cable_setBufVars ( buf% met% tk,    met% tk,     ktau, alloc )
   call cable_setBufVars ( buf% met% tvair, met% tvair,  ktau, alloc )
   call cable_setBufVars ( buf% met% tvrad, met% tvrad,  ktau, alloc )
   call cable_setBufVars ( buf% met% pmb,   met% pmb,    ktau, alloc )
   call cable_setBufVars ( buf% met% ua,    met% ua,     ktau, alloc )
   call cable_setBufVars ( buf% met% qv,    met% qv,     ktau, alloc )
   call cable_setBufVars ( buf% met% qvair, met% qvair,  ktau, alloc )
   call cable_setBufVars ( buf% met% da,    met% da,     ktau, alloc )
   call cable_setBufVars ( buf% met% dva,   met% dva,    ktau, alloc )
   call cable_setBufVars ( buf% met% coszen,met% coszen, ktau, alloc )
   call cable_setBufVars ( buf% met% fsd,   met% fsd,    ktau, alloc ) 
   call cable_setBufVars ( buf% met% precip_sn,  met% precip_sn,   ktau, alloc )
! ------------------------------------------------------------------------------
   call cable_setBufVars( buf% rad% extkb,   rad% extkb,   ktau, alloc )  
   call cable_setBufVars( buf% rad% extkd2,  rad% extkd2,  ktau, alloc )
   call cable_setBufVars( buf% rad% extkd,   rad% extkd,   ktau, alloc )
   call cable_setBufVars( buf% rad% flws,    rad% flws,    ktau, alloc )
   call cable_setBufVars( buf% rad% lwabv,   rad% lwabv,   ktau, alloc )
   call cable_setBufVars( buf% rad% qssabs,  rad% qssabs,  ktau, alloc )
   call cable_setBufVars( buf% rad% transd,  rad% transd,  ktau, alloc )
   call cable_setBufVars( buf% rad% trad,    rad% trad,    ktau, alloc )
   call cable_setBufVars( buf% rad% transb,  rad% transb,  ktau, alloc )
   call cable_setBufVars( buf% rad% workp1,  rad% workp1,  ktau, alloc )
   call cable_setBufVars( buf% rad% workp2,  rad% workp2,  ktau, alloc )
   call cable_setBufVars( buf% rad% workp3,  rad% workp3,  ktau, alloc )
   call cable_setBufVars( buf% rad% albedo_T,   rad% albedo_T,   ktau, alloc )
   call cable_setBufVars( buf% rad% latitude,   rad% latitude,   ktau, alloc )
   call cable_setBufVars( buf% rad% longitude,  rad% longitude,  ktau, alloc )
   call cable_setBufVars( buf% rad% fvlai,   rad% fvlai,   ktau, alloc )
   call cable_setBufVars( buf% rad% rniso,   rad% rniso,   ktau, alloc )
   call cable_setBufVars( buf% rad% scalex,  rad% scalex,  ktau, alloc )
   call cable_setBufVars( buf% rad% gradis,  rad% gradis,  ktau, alloc )
   call cable_setBufVars( buf% rad% albedo,  rad% albedo,  ktau, alloc ) 
   call cable_setBufVars( buf% rad% rhocdf,  rad% rhocdf,  ktau, alloc )
   call cable_setBufVars( buf% rad% reffdf,  rad% reffdf,  ktau, alloc )
   call cable_setBufVars( buf% rad% reffbm,  rad% reffbm,  ktau, alloc )
   call cable_setBufVars( buf% rad% extkbm,  rad% extkbm,  ktau, alloc )
   call cable_setBufVars( buf% rad% extkdm,  rad% extkdm,  ktau, alloc )
   call cable_setBufVars( buf% rad% cexpkbm, rad% cexpkbm, ktau, alloc )
   call cable_setBufVars( buf% rad% cexpkdm, rad% cexpkdm, ktau, alloc )
   call cable_setBufVars( buf% rad% fbeam,   rad% fbeam,   ktau, alloc )
   call cable_setBufVars( buf% rad% rhocbm,  rad% rhocbm,  ktau, alloc )
   call cable_setBufVars( buf% rad% qcan,    rad% qcan,    ktau, alloc )

   call cable_setBufVars ( buf% rough% coexp,   rough% coexp,   ktau, alloc )
   call cable_setBufVars ( buf% rough% disp,    rough% disp,    ktau, alloc )
   call cable_setBufVars ( buf% rough% hruff,   rough% hruff,   ktau, alloc )
   call cable_setBufVars ( buf% rough% rt0us,   rough% rt0us,   ktau, alloc )
   call cable_setBufVars ( buf% rough% rt1usa,  rough% rt1usa,  ktau, alloc )
   call cable_setBufVars ( buf% rough% rt1usb,  rough% rt1usb,  ktau, alloc )
   call cable_setBufVars ( buf% rough% rt1,     rough% rt1,     ktau, alloc )
   call cable_setBufVars ( buf% rough% term2,   rough% term2,   ktau, alloc )
   call cable_setBufVars ( buf% rough% term3,   rough% term3,   ktau, alloc )
   call cable_setBufVars ( buf% rough% term5,   rough% term5,   ktau, alloc )
   call cable_setBufVars ( buf% rough% term6,   rough% term6,   ktau, alloc )
   call cable_setBufVars ( buf% rough% usuh,    rough% usuh,    ktau, alloc )
   call cable_setBufVars ( buf% rough% za_uv,   rough% za_uv,   ktau, alloc )
   call cable_setBufVars ( buf% rough% za_tq,   rough% za_tq,   ktau, alloc )
   call cable_setBufVars ( buf% rough% z0m,     rough% z0m,     ktau, alloc )
   call cable_setBufVars ( buf% rough% zref_uv, rough% zref_uv, ktau, alloc )
   call cable_setBufVars ( buf% rough% zref_tq, rough% zref_tq, ktau, alloc )
   call cable_setBufVars ( buf% rough% zruffs,  rough% zruffs,  ktau, alloc )
   call cable_setBufVars ( buf% rough% z0soil,  rough% z0soil,  ktau, alloc )
   call cable_setBufVars ( buf% rough% hruff_grmx, rough% hruff_grmx, ktau, alloc )
   call cable_setBufVars ( buf% rough% z0soilsn,   rough% z0soilsn,   ktau, alloc )

   call cable_setBufVars ( buf% air% rho,    air% rho,    ktau, alloc )
   call cable_setBufVars ( buf% air% volm,   air% volm,   ktau, alloc )
   call cable_setBufVars ( buf% air% rlam,   air% rlam,   ktau, alloc )
   call cable_setBufVars ( buf% air% qsat,   air% qsat,   ktau, alloc )
   call cable_setBufVars ( buf% air% epsi,   air% epsi,   ktau, alloc )
   call cable_setBufVars ( buf% air% visc,   air% visc,   ktau, alloc )
   call cable_setBufVars ( buf% air% psyc,   air% psyc,   ktau, alloc )
   call cable_setBufVars ( buf% air% dsatdk, air% dsatdk, ktau, alloc )
   call cable_setBufVars ( buf% air% cmolar, air% cmolar, ktau, alloc )

! ------------------------------------------------------------------------------

END SUBROUTINE set_BufVars 

! ------------------------------------------------------------------------------

  
subroutine test_BufVars( ktau, buf, met, bgc, soil,veg, ssnow, canopy, &
                        rad, rough, air )

   use cable_def_types_mod, only : met_type, bgc_pool_type, soil_parameter_type, &
                                   veg_parameter_type, soil_snow_type, canopy_type, & 
                                   air_type, radiation_type, roughness_type
   use cable_Bufvars_mod, only : cable_testBufVars  
   
   integer :: ktau
   type (TbufVars) :: buf
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
 
   call cable_testBufVars( buf% bgc% cplant, bgc% cplant,   ktau, alloc )
   call cable_testBufVars( buf% bgc% csoil,  bgc% csoil,    ktau, alloc )
              
   call cable_testBufVars( buf% soil% bch,   soil% bch,     ktau, alloc )   
   call cable_testBufVars( buf% soil% c3,    soil% c3,      ktau, alloc )    
   call cable_testBufVars( buf% soil% clay,  soil% clay,    ktau, alloc )  
   call cable_testBufVars( buf% soil% css,   soil% css,     ktau, alloc )   
   call cable_testBufVars( buf% soil% hsbh,  soil% hsbh,    ktau, alloc )  
   call cable_testBufVars( buf% soil% hyds,  soil% hyds,    ktau, alloc )  
   call cable_testBufVars( buf% soil% i2bp3, soil% i2bp3,   ktau, alloc ) 
   call cable_testBufVars( buf% soil% ibp2,  soil% ibp2,    ktau, alloc )  
   call cable_testBufVars( buf% soil% sand,  soil% sand,    ktau, alloc )   
   call cable_testBufVars( buf% soil% sfc,   soil% sfc,     ktau, alloc )   
   call cable_testBufVars( buf% soil% silt,  soil% silt,    ktau, alloc )   
   call cable_testBufVars( buf% soil% ssat,  soil% ssat,    ktau, alloc )   
   call cable_testBufVars( buf% soil% sucs,  soil% sucs,    ktau, alloc )   
   call cable_testBufVars( buf% soil% swilt, soil% swilt,   ktau, alloc )  
   call cable_testBufVars( buf% soil% zse,   soil% zse,     ktau, alloc )    
   call cable_testBufVars( buf% soil% zshh, soil% zshh,     ktau, alloc )  
   call cable_testBufVars( buf% soil% cnsd, soil% cnsd,     ktau, alloc )  
   call cable_testBufVars( buf% soil% isoilm,  soil% isoilm,  ktau, alloc )  
   call cable_testBufVars( buf% soil% rhosoil, soil% rhosoil, ktau, alloc )  
   
   call cable_testBufVars( buf% soil% albsoil,  soil% albsoil,  ktau, alloc)   
   call cable_testBufVars( buf% soil% pwb_min,  soil% pwb_min,  ktau, alloc )  
   call cable_testBufVars( buf% soil% albsoilf, soil% albsoilf, ktau, alloc )  
   
   call cable_testBufVars( buf% ssnow% iantrct,   ssnow% iantrct,   ktau, alloc )
   call cable_testBufVars( buf% ssnow % pudsto,   ssnow % pudsto,   ktau, alloc )
   call cable_testBufVars( buf% ssnow % pudsmx,   ssnow % pudsmx,   ktau, alloc )
   call cable_testBufVars( buf% ssnow% cls,       ssnow% cls,       ktau, alloc )     
   call cable_testBufVars( buf% ssnow% dfn_dtg,   ssnow% dfn_dtg,   ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% dfh_dtg,   ssnow% dfh_dtg,   ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% dfe_ddq,   ssnow% dfe_ddq,   ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% ddq_dtg,   ssnow% ddq_dtg,   ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% evapsn,    ssnow% evapsn,    ktau, alloc )  
   call cable_testBufVars( buf% ssnow% fwtop,     ssnow% fwtop,     ktau, alloc )   
   call cable_testBufVars( buf% ssnow% fwtop1,    ssnow% fwtop1,    ktau, alloc )   
   call cable_testBufVars( buf% ssnow% fwtop2,    ssnow% fwtop2,    ktau, alloc )   
   call cable_testBufVars( buf% ssnow% fwtop3,    ssnow% fwtop3,    ktau, alloc )   
   call cable_testBufVars( buf% ssnow% gammzz,    ssnow% gammzz,    ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% isflag,    ssnow% isflag,    ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% osnowd,    ssnow% osnowd,    ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% potev,     ssnow% potev,     ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% runoff,    ssnow% runoff,    ktau, alloc )
   call cable_testBufVars( buf% ssnow% rnof1,     ssnow% rnof1,     ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% rnof2,     ssnow% rnof2,     ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% rtsoil,    ssnow% rtsoil,    ktau, alloc )
   call cable_testBufVars( buf% ssnow% sconds,    ssnow% sconds,    ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% sdepth,    ssnow% sdepth,    ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% smass,     ssnow% smass,     ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% snage,     ssnow% snage,     ktau, alloc )  
   call cable_testBufVars( buf% ssnow% snowd,     ssnow% snowd,     ktau, alloc )  
   call cable_testBufVars( buf% ssnow% smelt,     ssnow% smelt,     ktau, alloc )  
   call cable_testBufVars( buf% ssnow% ssdn,      ssnow% ssdn,      ktau, alloc )  
   call cable_testBufVars( buf% ssnow% ssdnn,     ssnow% ssdnn,     ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% tgg,       ssnow% tgg,       ktau, alloc )   
   call cable_testBufVars( buf% ssnow% tggsn,     ssnow% tggsn,     ktau, alloc )  
   call cable_testBufVars( buf% ssnow% tss,       ssnow% tss,       ktau, alloc )   
   call cable_testBufVars( buf% ssnow% tss_p,     ssnow% tss_p,     ktau, alloc )   
   call cable_testBufVars( buf% ssnow% deltss,    ssnow% deltss,    ktau, alloc )   
   call cable_testBufVars( buf% ssnow% owb1,      ssnow% owb1,      ktau, alloc )   
   call cable_testBufVars( buf% ssnow% wb,        ssnow% wb,        ktau, alloc )    
   call cable_testBufVars( buf% ssnow% wbice,     ssnow% wbice,     ktau, alloc ) 
   call cable_testBufVars( buf% ssnow% wblf,      ssnow% wblf,      ktau, alloc ) 
   call cable_testBufVars( buf% ssnow%wbtot,      ssnow%wbtot,      ktau, alloc )    
   call cable_testBufVars( buf% ssnow%wbtot1,     ssnow%wbtot1,     ktau, alloc )    
   call cable_testBufVars( buf% ssnow%wbtot2,     ssnow%wbtot2,     ktau, alloc )    
   call cable_testBufVars( buf% ssnow%wb_lake,    ssnow%wb_lake,    ktau, alloc )    
   call cable_testBufVars( buf% ssnow%sinfil,     ssnow%sinfil,     ktau, alloc )    
   call cable_testBufVars( buf% ssnow%evapfbl,    ssnow%evapfbl,    ktau, alloc )    
   call cable_testBufVars( buf% ssnow%qstss,      ssnow%qstss,      ktau, alloc )    
   call cable_testBufVars( buf% ssnow%wetfac,     ssnow%wetfac,     ktau, alloc )  
   call cable_testBufVars( buf% ssnow%owetfac,    ssnow%owetfac,    ktau, alloc )  
   call cable_testBufVars( buf% ssnow%t_snwlr,    ssnow%t_snwlr,    ktau, alloc )  
   call cable_testBufVars( buf% ssnow%wbfice,     ssnow%wbfice,     ktau, alloc )  
   call cable_testBufVars( buf% ssnow%tggav,      ssnow%tggav,      ktau, alloc )  
   call cable_testBufVars( buf% ssnow%otgg,       ssnow%otgg,       ktau, alloc )   
   call cable_testBufVars( buf% ssnow%otss,       ssnow%otss,       ktau, alloc )   
   call cable_testBufVars( buf% ssnow%otss_0,     ssnow%otss_0,     ktau, alloc )   
   call cable_testBufVars( buf% ssnow%tprecip,    ssnow%tprecip,    ktau, alloc ) 
   call cable_testBufVars( buf% ssnow%tevap,      ssnow%tevap,      ktau, alloc ) 
   call cable_testBufVars( buf% ssnow%trnoff,     ssnow%trnoff,     ktau, alloc ) 
   call cable_testBufVars( buf% ssnow%totenbal,   ssnow%totenbal,   ktau, alloc ) 
   call cable_testBufVars( buf% ssnow%totenbal2,  ssnow%totenbal2,  ktau, alloc ) 
   call cable_testBufVars( buf% ssnow%fland,      ssnow%fland,      ktau, alloc )      
   call cable_testBufVars( buf% ssnow%ifland,     ssnow%ifland,     ktau, alloc )  
   call cable_testBufVars( buf% ssnow%tilefrac,   ssnow%tilefrac,   ktau, alloc ) 
   call cable_testBufVars( buf% ssnow%qasrf,      ssnow%qasrf,      ktau, alloc )  
   call cable_testBufVars( buf% ssnow%qfsrf,      ssnow%qfsrf,      ktau, alloc )  
   call cable_testBufVars( buf% ssnow%qssrf,      ssnow%qssrf,      ktau, alloc )  

   call cable_testBufVars( buf% ssnow % dtmlt,    ssnow % dtmlt,     ktau, alloc )
   call cable_testBufVars( buf% ssnow% albsoilsn, ssnow% albsoilsn,  ktau, alloc ) 

   call cable_testBufVars( buf% veg% canst1,   veg% canst1,   ktau, alloc ) 
   call cable_testBufVars( buf% veg% dleaf,    veg% dleaf,    ktau, alloc )  
   call cable_testBufVars( buf% veg% ejmax,    veg% ejmax,    ktau, alloc ) 
   call cable_testBufVars( buf% veg% iveg,     veg% iveg,     ktau, alloc ) 
   call cable_testBufVars( buf% veg% meth,     veg% meth,     ktau, alloc ) 
   call cable_testBufVars( buf% veg% frac4,    veg% frac4,    ktau, alloc )  
   call cable_testBufVars( buf% veg% hc,       veg% hc,       ktau, alloc )     
   call cable_testBufVars( buf% veg% vlai,     veg% vlai,     ktau, alloc )   
   call cable_testBufVars( buf% veg% xalbnir,  veg% xalbnir,  ktau, alloc ) 
   call cable_testBufVars( buf% veg% rp20,     veg% rp20,     ktau, alloc )   
   call cable_testBufVars( buf% veg% rpcoef,   veg% rpcoef,   ktau, alloc ) 
   call cable_testBufVars( buf% veg% rs20,     veg% rs20,     ktau, alloc )   
   call cable_testBufVars( buf% veg% shelrb,   veg% shelrb,   ktau, alloc ) 
   call cable_testBufVars( buf% veg% vegcf,    veg% vegcf,    ktau, alloc )  
   call cable_testBufVars( buf% veg% tminvj,   veg% tminvj,   ktau, alloc ) 
   call cable_testBufVars( buf% veg% tmaxvj,   veg% tmaxvj,   ktau, alloc ) 
   call cable_testBufVars( buf% veg% vbeta,    veg% vbeta,    ktau, alloc )  
   call cable_testBufVars( buf% veg% vcmax,    veg% vcmax,    ktau, alloc )  
   call cable_testBufVars( buf% veg% xfang,    veg% xfang,    ktau, alloc )  
   call cable_testBufVars( buf% veg% extkn,    veg% extkn,     ktau, alloc ) 
   call cable_testBufVars( buf% veg% wai,      veg% wai,       ktau, alloc )   
   call cable_testBufVars( buf% veg% deciduous,veg% deciduous, ktau, alloc ) 
   call cable_testBufVars( buf% veg% froot,    veg% froot,     ktau, alloc ) 
   call cable_testBufVars( buf% veg% refl,     veg% refl,      ktau, alloc ) !jhan:swb?
   call cable_testBufVars( buf% veg% taul,     veg% taul,      ktau, alloc ) 
   call cable_testBufVars( buf% veg% vlaimax,  veg% vlaimax,   ktau, alloc ) 

   call cable_testBufVars( buf% canopy% fess,   canopy% fess,   ktau, alloc )
   call cable_testBufVars( buf% canopy% fesp,   canopy% fesp,   ktau, alloc )
   call cable_testBufVars( buf% canopy% cansto, canopy% cansto, ktau, alloc )  
   call cable_testBufVars( buf% canopy% cduv,   canopy% cduv,   ktau, alloc )   
   call cable_testBufVars( buf% canopy% delwc,  canopy% delwc,  ktau, alloc )  
   call cable_testBufVars( buf% canopy% dewmm,  canopy% dewmm,  ktau, alloc )  
   call cable_testBufVars( buf% canopy% dgdtg,  canopy% dgdtg,  ktau, alloc )  
   call cable_testBufVars( buf% canopy% fe,     canopy% fe,     ktau, alloc )      
   call cable_testBufVars( buf% canopy% fh,     canopy% fh,     ktau, alloc )      
   call cable_testBufVars( buf% canopy% fpn,    canopy% fpn,    ktau, alloc )     
   call cable_testBufVars( buf% canopy% frp,    canopy% frp,    ktau, alloc )     
   call cable_testBufVars( buf% canopy% frpw,   canopy% frpw,   ktau, alloc )    
   call cable_testBufVars( buf% canopy% frpr,   canopy% frpr,   ktau, alloc )    
   call cable_testBufVars( buf% canopy% frs,    canopy% frs,    ktau, alloc )     
   call cable_testBufVars( buf% canopy% fnee,   canopy% fnee,   ktau, alloc )    
   call cable_testBufVars( buf% canopy% frday,  canopy% frday,  ktau, alloc )   
   call cable_testBufVars( buf% canopy% fnv,    canopy% fnv,    ktau, alloc )     
   call cable_testBufVars( buf% canopy% fev,    canopy% fev,    ktau, alloc )     
   call cable_testBufVars( buf% canopy% fevc,   canopy% fevc,   ktau, alloc )    
   call cable_testBufVars( buf% canopy% fhv,    canopy% fhv,    ktau, alloc )     
   call cable_testBufVars( buf% canopy% fns,    canopy% fns,    ktau, alloc )     
   call cable_testBufVars( buf% canopy% fhs,    canopy% fhs,    ktau, alloc )     
   call cable_testBufVars( buf% canopy% fhs_cor,canopy% fhs_cor,ktau, alloc )     
   call cable_testBufVars( buf% canopy% ga,     canopy% ga,     ktau, alloc )      
   call cable_testBufVars( buf% canopy% ghflux, canopy% ghflux, ktau, alloc )   
   call cable_testBufVars( buf% canopy% precis, canopy% precis, ktau, alloc ) 
   call cable_testBufVars( buf% canopy% qscrn,  canopy% qscrn,  ktau, alloc )  
   call cable_testBufVars( buf% canopy% rnet,   canopy% rnet,   ktau, alloc )   
   call cable_testBufVars( buf% canopy% segg,   canopy% segg,   ktau, alloc )   
   call cable_testBufVars( buf% canopy% sghflux,canopy% sghflux,ktau, alloc )  
   call cable_testBufVars( buf% canopy% through,canopy% through,ktau, alloc )  
   call cable_testBufVars( buf% canopy% spill,  canopy% spill,  ktau, alloc )  
   call cable_testBufVars( buf% canopy% tscrn,  canopy% tscrn,  ktau, alloc )  
   call cable_testBufVars( buf% canopy% wcint,  canopy% wcint,  ktau, alloc )  
   call cable_testBufVars( buf% canopy% tv,     canopy% tv,     ktau, alloc )      
   call cable_testBufVars( buf% canopy% us,     canopy% us,     ktau, alloc )      
   call cable_testBufVars( buf% canopy% uscrn,  canopy% uscrn,  ktau, alloc )   
   call cable_testBufVars( buf% canopy% rghlai, canopy% rghlai, ktau, alloc ) 
   call cable_testBufVars( buf% canopy% vlaiw,  canopy% vlaiw,  ktau, alloc ) 
   call cable_testBufVars( buf% canopy% fwet,   canopy% fwet,   ktau, alloc )   
   call cable_testBufVars( buf% canopy% epot,   canopy% epot,   ktau, alloc )   
   call cable_testBufVars( buf% canopy% fnpp,   canopy% fnpp,   ktau, alloc )   
   call cable_testBufVars( buf% canopy% gswx_T, canopy% gswx_T, ktau, alloc )  
   call cable_testBufVars( buf% canopy% cdtq,   canopy% cdtq,   ktau, alloc )   
   call cable_testBufVars( buf% canopy% fevw,   canopy% fevw,   ktau, alloc )   
   call cable_testBufVars( buf% canopy% fhvw,   canopy% fhvw,   ktau, alloc )   
   call cable_testBufVars( buf% canopy% fes,    canopy% fes,    ktau, alloc )    
   call cable_testBufVars( buf% canopy% fes_cor,canopy% fes_cor,ktau, alloc )    
   call cable_testBufVars( buf% canopy% gswx,   canopy% gswx,   ktau, alloc )  
   call cable_testBufVars( buf% canopy% evapfbl,canopy% evapfbl,ktau, alloc )
   call cable_testBufVars( buf% canopy% wetfac_cs, canopy% wetfac_cs, ktau, alloc )  
   call cable_testBufVars( buf% canopy% oldcansto, canopy% oldcansto, ktau, alloc )  
   call cable_testBufVars( buf% canopy% zetar,     canopy% zetar,     ktau, alloc )  
   call cable_testBufVars( buf% canopy% fevw_pot,  canopy% fevw_pot,  ktau, alloc )  
 
   call cable_testBufVars ( buf% met% ca,    met% ca,     ktau, alloc )
   call cable_testBufVars ( buf% met% year,  met% year,   ktau, alloc )
   call cable_testBufVars ( buf% met% moy,   met% moy,    ktau, alloc )
   call cable_testBufVars ( buf% met% doy,   met% doy,    ktau, alloc )
   call cable_testBufVars ( buf% met% hod,   met% hod,    ktau, alloc )
   call cable_testBufVars ( buf% met% ofsd,  met% ofsd,   ktau, alloc ) 
   call cable_testBufVars ( buf% met% fld,   met% fld,    ktau, alloc )
   call cable_testBufVars ( buf% met% precip,met% precip, ktau, alloc )
   call cable_testBufVars ( buf% met% tk,    met% tk,     ktau, alloc )
   call cable_testBufVars ( buf% met% tvair, met% tvair,  ktau, alloc )
   call cable_testBufVars ( buf% met% tvrad, met% tvrad,  ktau, alloc )
   call cable_testBufVars ( buf% met% pmb,   met% pmb,    ktau, alloc )
   call cable_testBufVars ( buf% met% ua,    met% ua,     ktau, alloc )
   call cable_testBufVars ( buf% met% qv,    met% qv,     ktau, alloc )
   call cable_testBufVars ( buf% met% qvair, met% qvair,  ktau, alloc )
   call cable_testBufVars ( buf% met% da,    met% da,     ktau, alloc )
   call cable_testBufVars ( buf% met% dva,   met% dva,    ktau, alloc )
   call cable_testBufVars ( buf% met% coszen,met% coszen, ktau, alloc )
   call cable_testBufVars ( buf% met% fsd,   met% fsd,    ktau, alloc ) 
   call cable_testBufVars ( buf% met% precip_sn,  met% precip_sn,   ktau, alloc )
! ------------------------------------------------------------------------------
   
   call cable_testBufVars( buf% rad% extkb,   rad% extkb,   ktau, alloc )  
   call cable_testBufVars( buf% rad% extkd2,  rad% extkd2,  ktau, alloc )
   call cable_testBufVars( buf% rad% extkd,   rad% extkd,   ktau, alloc )
   call cable_testBufVars( buf% rad% flws,    rad% flws,    ktau, alloc )
   call cable_testBufVars( buf% rad% lwabv,   rad% lwabv,   ktau, alloc )
   call cable_testBufVars( buf% rad% qssabs,  rad% qssabs,  ktau, alloc )
   call cable_testBufVars( buf% rad% transd,  rad% transd,  ktau, alloc )
   call cable_testBufVars( buf% rad% trad,    rad% trad,    ktau, alloc )
   call cable_testBufVars( buf% rad% transb,  rad% transb,  ktau, alloc )
   call cable_testBufVars( buf% rad% workp1,  rad% workp1,  ktau, alloc )
   call cable_testBufVars( buf% rad% workp2,  rad% workp2,  ktau, alloc )
   call cable_testBufVars( buf% rad% workp3,  rad% workp3,  ktau, alloc )
   call cable_testBufVars( buf% rad% albedo_T,   rad% albedo_T,   ktau, alloc )
   call cable_testBufVars( buf% rad% latitude,   rad% latitude,   ktau, alloc )
   call cable_testBufVars( buf% rad% longitude,  rad% longitude,  ktau, alloc )
   call cable_testBufVars( buf% rad% fvlai,   rad% fvlai,   ktau, alloc )
   call cable_testBufVars( buf% rad% rniso,   rad% rniso,   ktau, alloc )
   call cable_testBufVars( buf% rad% scalex,  rad% scalex,  ktau, alloc )
   call cable_testBufVars( buf% rad% gradis,  rad% gradis,  ktau, alloc )
   call cable_testBufVars( buf% rad% albedo,  rad% albedo,  ktau, alloc ) 
   call cable_testBufVars( buf% rad% rhocdf,  rad% rhocdf,  ktau, alloc )
   call cable_testBufVars( buf% rad% reffdf,  rad% reffdf,  ktau, alloc )
   call cable_testBufVars( buf% rad% reffbm,  rad% reffbm,  ktau, alloc )
   call cable_testBufVars( buf% rad% extkbm,  rad% extkbm,  ktau, alloc )
   call cable_testBufVars( buf% rad% extkdm,  rad% extkdm,  ktau, alloc )
   call cable_testBufVars( buf% rad% cexpkbm, rad% cexpkbm, ktau, alloc )
   call cable_testBufVars( buf% rad% cexpkdm, rad% cexpkdm, ktau, alloc )
   call cable_testBufVars( buf% rad% fbeam,   rad% fbeam,   ktau, alloc )
   call cable_testBufVars( buf% rad% rhocbm,  rad% rhocbm,  ktau, alloc )
   call cable_testBufVars( buf% rad% qcan,    rad% qcan,    ktau, alloc )

   call cable_testBufVars ( buf% rough% coexp,   rough% coexp,   ktau, alloc )
   call cable_testBufVars ( buf% rough% disp,    rough% disp,    ktau, alloc )
   call cable_testBufVars ( buf% rough% hruff,   rough% hruff,   ktau, alloc )
   call cable_testBufVars ( buf% rough% rt0us,   rough% rt0us,   ktau, alloc )
   call cable_testBufVars ( buf% rough% rt1usa,  rough% rt1usa,  ktau, alloc )
   call cable_testBufVars ( buf% rough% rt1usb,  rough% rt1usb,  ktau, alloc )
   call cable_testBufVars ( buf% rough% rt1,     rough% rt1,     ktau, alloc )
   call cable_testBufVars ( buf% rough% term2,   rough% term2,   ktau, alloc )
   call cable_testBufVars ( buf% rough% term3,   rough% term3,   ktau, alloc )
   call cable_testBufVars ( buf% rough% term5,   rough% term5,   ktau, alloc )
   call cable_testBufVars ( buf% rough% term6,   rough% term6,   ktau, alloc )
   call cable_testBufVars ( buf% rough% usuh,    rough% usuh,    ktau, alloc )
   call cable_testBufVars ( buf% rough% za_uv,   rough% za_uv,   ktau, alloc )
   call cable_testBufVars ( buf% rough% za_tq,   rough% za_tq,   ktau, alloc )
   call cable_testBufVars ( buf% rough% z0m,     rough% z0m,     ktau, alloc )
   call cable_testBufVars ( buf% rough% zref_uv, rough% zref_uv, ktau, alloc )
   call cable_testBufVars ( buf% rough% zref_tq, rough% zref_tq, ktau, alloc )
   call cable_testBufVars ( buf% rough% zruffs,  rough% zruffs,  ktau, alloc )
   call cable_testBufVars ( buf% rough% z0soil,  rough% z0soil,  ktau, alloc )
   call cable_testBufVars ( buf% rough% hruff_grmx, rough% hruff_grmx, ktau, alloc )
   call cable_testBufVars ( buf% rough% z0soilsn,   rough% z0soilsn,   ktau, alloc )

   call cable_testBufVars ( buf% air% rho,    air% rho,    ktau, alloc )
   call cable_testBufVars ( buf% air% volm,   air% volm,   ktau, alloc )
   call cable_testBufVars ( buf% air% rlam,   air% rlam,   ktau, alloc )
   call cable_testBufVars ( buf% air% qsat,   air% qsat,   ktau, alloc )
   call cable_testBufVars ( buf% air% epsi,   air% epsi,   ktau, alloc )
   call cable_testBufVars ( buf% air% visc,   air% visc,   ktau, alloc )
   call cable_testBufVars ( buf% air% psyc,   air% psyc,   ktau, alloc )
   call cable_testBufVars ( buf% air% dsatdk, air% dsatdk, ktau, alloc )
   call cable_testBufVars ( buf% air% cmolar, air% cmolar, ktau, alloc )

! ------------------------------------------------------------------------------

END SUBROUTINE test_BufVars 

! ------------------------------------------------------------------------------





!SUBROUTINE alloc_balances_type(var, mp)
!   
!   TYPE(balances_type) :: var
!   INTEGER :: mp
!   
!   call cable_BufVars( var% drybal, mp, alloc ) 
!   call cable_BufVars( var% ebal, mp, alloc )  
!   call cable_BufVars( var% ebal_tot, mp, alloc )
!   call cable_BufVars( var% ebaltr, mp, alloc )  
!   call cable_BufVars( var% ebal_tottr, mp, alloc )
!   call cable_BufVars( var% ebal_cncheck, mp, alloc )  
!   call cable_BufVars( var% ebal_tot_cncheck, mp, alloc )
!   call cable_BufVars( var% evap_tot, mp, alloc ) 
!   call cable_BufVars( var% osnowd0, mp, alloc ) 
!   call cable_BufVars( var% precip_tot, mp, alloc  )
!   call cable_BufVars( var% rnoff_tot, mp, alloc ) 
!   call cable_BufVars( var% wbal, mp, alloc )    
!   call cable_BufVars( var% wbal_tot, mp, alloc  )
!   call cable_BufVars( var% wbtot0, mp, alloc )  
!   call cable_BufVars( var% wetbal, mp, alloc ) 
!   call cable_BufVars( var% cansto0, mp, alloc  ) 
!   call cable_BufVars( var% evapc_tot, mp, alloc ) 
!   call cable_BufVars( var% evaps_tot, mp, alloc ) 
!   call cable_BufVars( var% rnof1_tot, mp, alloc ) 
!   call cable_BufVars( var% rnof2_tot, mp, alloc ) 
!   call cable_BufVars( var% snowdc_tot, mp, alloc )
!   call cable_BufVars( var% wbal_tot1, mp, alloc ) 
!   call cable_BufVars( var% owbtot, mp, alloc )  
!   call cable_BufVars( var% delwc_tot, mp, alloc )  
!   call cable_BufVars( var% qasrf_tot, mp, alloc ) 
!   call cable_BufVars( var% qfsrf_tot, mp, alloc )  
!   call cable_BufVars( var% qssrf_tot, mp, alloc )  
!
!END SUBROUTINE alloc_balances_type
!
!! ------------------------------------------------------------------------------

!SUBROUTINE alloc_sum_flux_type(var, mp)
!
!   TYPE(sum_flux_type) :: var
!   INTEGER :: mp
! 
!   call cable_BufVars ( var % sumpn, mp, alloc )
!   call cable_BufVars ( var % sumrp, mp, alloc )
!   call cable_BufVars ( var % sumrpw, mp, alloc )
!   call cable_BufVars ( var % sumrpr, mp, alloc )
!   call cable_BufVars ( var % sumrs, mp, alloc )
!   call cable_BufVars ( var % sumrd, mp, alloc )
!   call cable_BufVars ( var % dsumpn, mp, alloc )
!   call cable_BufVars ( var % dsumrp, mp, alloc )
!   call cable_BufVars ( var % dsumrs, mp, alloc )
!   call cable_BufVars ( var % dsumrd, mp, alloc )
!   call cable_BufVars ( var % sumxrp, mp, alloc )
!   call cable_BufVars ( var % sumxrs, mp, alloc )
!
!END SUBROUTINE alloc_sum_flux_type

! ------------------------------------------------------------------------------


end MODULE cable_buffer_types_mod
