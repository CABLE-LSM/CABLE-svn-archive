!jhan:re-write arg list to alloc subrs -> subr(var), use define_dimensions in subrs NOT module
MODULE define_types
  ! Contains all variables which are not subroutine-internal
  !USE define_dimensions, only : mp !soil_snow etc relies on thisinheritance 
  USE define_dimensions, only : i_d, r_1, r_2, ms, msn, mf, nrb, ncp, ncs, swb, n_tiles
  implicit none
  !MANY subrs rely on inheritance of these ATM
  private :: i_d, r_1, r_2, ms, msn, mf, nrb, ncp, ncs
  
  ! Energy and water balance variables:
  TYPE balances_type 
     REAL, DIMENSION(:), POINTER :: drybal ! energy balance for dry canopy
     REAL, DIMENSION(:), POINTER :: ebal   ! energy balance per time step (W/m^2)
     REAL, DIMENSION(:), POINTER :: ebal_tot ! cumulative energy balance (W/m^2)
     REAL, DIMENSION(:), POINTER :: ebal_cncheck ! energy balance consistency check (W/m^2)
     REAL, DIMENSION(:), POINTER :: ebal_tot_cncheck ! cumulative energy balance (W/m^2)
     REAL, DIMENSION(:), POINTER :: ebaltr   ! energy balance per time step (W/m^2)
     REAL, DIMENSION(:), POINTER :: ebal_tottr ! cumulative energy balance (W/m^2)
     REAL, DIMENSION(:), POINTER :: evap_tot ! cumulative evapotranspiration (mm/dels)
     REAL, DIMENSION(:), POINTER :: osnowd0  ! snow depth, first time step
     REAL, DIMENSION(:), POINTER :: precip_tot ! cumulative precipitation (mm/dels)
     REAL, DIMENSION(:), POINTER :: rnoff_tot  ! cumulative runoff (mm/dels)
     REAL, DIMENSION(:), POINTER :: wbal   ! water balance per time step (mm/dels)
     REAL, DIMENSION(:), POINTER :: wbal_tot ! cumulative water balance (mm/dels)
     REAL, DIMENSION(:), POINTER :: wbtot0 ! total soil water (mm), first time step
     REAL, DIMENSION(:), POINTER :: wetbal ! energy balance for wet canopy
!jhan:prev. def. in UM only
     !not used
     REAL, DIMENSION(:), POINTER :: cansto0 ! canopy water storage (mm)
!checks only
     REAL, DIMENSION(:), POINTER :: owbtot ! total soil water (mm), first time step
     REAL, DIMENSION(:), POINTER :: evapc_tot ! cumulative evapotranspiration (mm/dels)
     REAL, DIMENSION(:), POINTER :: evaps_tot ! cumulative evapotranspiration (mm/dels)
     REAL, DIMENSION(:), POINTER :: rnof1_tot  ! cumulative runoff (mm/dels)
     REAL, DIMENSION(:), POINTER :: rnof2_tot  ! cumulative runoff (mm/dels)
     REAL, DIMENSION(:), POINTER :: snowdc_tot  ! cumulative runoff (mm/dels)
     REAL, DIMENSION(:), POINTER :: wbal_tot1 ! cumulative water balance (mm/dels)
     REAL, DIMENSION(:), POINTER :: delwc_tot ! energy balance for wet canopy
     REAL, DIMENSION(:), POINTER :: qasrf_tot ! heat advected to the snow by precip. 
     REAL, DIMENSION(:), POINTER :: qfsrf_tot ! energy of snowpack phase changes 
     REAL, DIMENSION(:), POINTER :: qssrf_tot ! energy of snowpack phase changes 
!jhan:}
  END TYPE balances_type
  ! Soil parameters:
  TYPE soil_parameter_type 
     REAL, DIMENSION(:), POINTER :: bch  ! parameter b in Campbell equation
     REAL, DIMENSION(:), POINTER :: c3   ! c3 drainage coeff (fraction)
     REAL, DIMENSION(:), POINTER :: clay ! fraction of soil which is clay
     REAL, DIMENSION(:), POINTER :: css  ! soil specific heat capacity [kJ/kg/K]
     REAL, DIMENSION(:), POINTER :: hsbh  ! difsat * etasat (=hyds*abs(sucs)*bch)
     REAL, DIMENSION(:), POINTER :: hyds  ! hydraulic conductivity @ saturation [m/s], Ksat
     INTEGER, DIMENSION(:), POINTER :: i2bp3  ! par. one in K vis suction (=nint(bch)+2)
     INTEGER, DIMENSION(:), POINTER :: ibp2   ! par. two in K vis suction (fn of pbch)
     INTEGER, DIMENSION(:), POINTER :: isoilm ! integer soil type
     REAL, DIMENSION(:), POINTER :: rhosoil ! soil density [kg/m3]
     !REAL, DIMENSION(:), POINTER :: rs20  ! soil respiration at 20 C [mol m-2 s-1]
     REAL, DIMENSION(:), POINTER :: sand  ! fraction of soil which is sand
     REAL, DIMENSION(:), POINTER :: sfc   ! vol H2O @ field capacity
     REAL, DIMENSION(:), POINTER :: silt  ! fraction of soil which is silt
     REAL, DIMENSION(:), POINTER :: ssat  ! vol H2O @ saturation
     REAL, DIMENSION(:), POINTER :: sucs  ! suction at saturation (m)
     REAL, DIMENSION(:), POINTER :: swilt ! vol H2O @ wilting
     REAL, DIMENSION(:), POINTER  :: zse   ! thickness of each soil layer (1=top) in m
     REAL, DIMENSION(:), POINTER  :: zshh ! distance between consecutive layer midpoints (m)
     REAL(r_2), DIMENSION(:), POINTER :: cnsd ! thermal conductivity of dry soil [W/m/K]
     REAL(r_2), DIMENSION(:), POINTER :: pwb_min ! working variable (swilt/ssat)**ibp2
     REAL, DIMENSION(:,:), POINTER :: albsoil ! soil reflectance (2nd dim. BP 21Oct2009)
!jhan:prev. def. in UM only
     REAL, DIMENSION(:), POINTER :: albsoilf ! soil reflectance
!jhan:}
  END TYPE soil_parameter_type

  ! Soil and snow variables:
  TYPE soil_snow_type 
     REAL, DIMENSION(:), POINTER :: iantrct ! pointer to Antarctic land points
     REAL, DIMENSION(:,:), POINTER :: dtmlt   ! water flux to the soil
     REAL, DIMENSION(:), POINTER :: pudsto  ! puddle storage
     REAL, DIMENSION(:), POINTER :: pudsmx  ! puddle storage
     REAL, DIMENSION(:,:), POINTER :: albsoilsn ! soil + snow reflectance
     REAL, DIMENSION(:), POINTER :: cls     ! factor for latent heat
     REAL, DIMENSION(:), POINTER :: dfn_dtg ! d(canopy%fns)/d(ssoil%tgg)
     REAL, DIMENSION(:), POINTER :: dfh_dtg ! d(canopy%fhs)/d(ssoil%tgg)
     REAL, DIMENSION(:), POINTER :: dfe_ddq ! d(canopy%fes)/d(dq)
     REAL, DIMENSION(:), POINTER :: ddq_dtg ! d(dq)/d(ssoil%tgg)
     REAL, DIMENSION(:), POINTER :: evapsn  ! snow evaporation  
     REAL, DIMENSION(:), POINTER :: fwtop   ! water flux to the soil
 !jhan:Eva added these 3 vars. we have to clean up all o fthis anyway
        REAL, DIMENSION(:), POINTER :: fwtop1   ! water flux to the soil
     REAL, DIMENSION(:), POINTER :: fwtop2   ! water flux to the soil
     REAL, DIMENSION(:), POINTER :: fwtop3   ! water flux to the soil
     REAL(r_2), DIMENSION(:,:), POINTER :: gammzz ! heat capacity for each soil layer
     INTEGER, DIMENSION(:), POINTER :: isflag ! 0 => no snow 1 => snow
     REAL, DIMENSION(:), POINTER :: osnowd  ! snow depth from previous time step
     REAL, DIMENSION(:), POINTER :: potev   ! potential evapotranspiration
     REAL, DIMENSION(:), POINTER :: runoff  ! total runoff (mm/dels)
     REAL, DIMENSION(:), POINTER :: rnof1   ! surface runoff (mm/dels)
     REAL, DIMENSION(:), POINTER :: rnof2   ! deep drainage (mm/dels)
     REAL, DIMENSION(:), POINTER :: rtsoil  ! turbulent resistance for soil
     REAL, DIMENSION(:,:), POINTER :: sconds !
     REAL, DIMENSION(:,:), POINTER :: sdepth ! snow depth
     REAL, DIMENSION(:,:), POINTER  :: smass  ! snow mass
     REAL, DIMENSION(:), POINTER :: snage   ! snow age
     REAL, DIMENSION(:), POINTER :: snowd   ! snow depth (liquid water)
     REAL, DIMENSION(:), POINTER :: smelt   ! snow melt 
     REAL, DIMENSION(:,:), POINTER  :: ssdn ! snow densities
     REAL, DIMENSION(:), POINTER :: ssdnn   ! average snow density
     REAL, DIMENSION(:,:), POINTER :: tgg   ! soil temperature in K
     REAL, DIMENSION(:,:), POINTER  :: tggsn ! snow temperature in K
     REAL, DIMENSION(:), POINTER :: tss     ! surface temperature (weighted soil, snow)
     REAL, DIMENSION(:), POINTER :: tss_p     ! surface temperature (weighted soil, snow)
     REAL, DIMENSION(:), POINTER :: deltss     ! surface temperature (weighted soil, snow)
     REAL, DIMENSION(:), POINTER :: owb1     ! surface temperature (weighted soil, snow)
     REAL(r_2), DIMENSION(:,:), POINTER :: wb    ! volumetric soil moisture (solid+liq)
     REAL(r_2), DIMENSION(:,:), POINTER :: wbice  ! soil ice
     REAL(r_2), DIMENSION(:,:), POINTER :: wblf !
     REAL(r_2), DIMENSION(:), POINTER :: wbtot   ! total soil water (mm)
 !jhan:Eva added these 6 vars. ome are for lkaes?
     REAL, DIMENSION(:), POINTER :: wbtot1   ! total soil water (mm)
     REAL, DIMENSION(:), POINTER :: wbtot2   ! total soil water (mm)
     REAL, DIMENSION(:), POINTER :: wb_lake
     REAL, DIMENSION(:), POINTER :: sinfil 
     REAL, DIMENSION(:,:), POINTER :: evapfbl
     REAL, DIMENSION(:), POINTER :: qstss 
     REAL, DIMENSION(:), POINTER :: wetfac ! surface wetness fact. at current time step
     REAL, DIMENSION(:), POINTER :: owetfac ! surface wetness fact. at previous time step
     REAL, DIMENSION(:), POINTER :: t_snwlr ! top snow layer depth in 3 layer snowpack
     REAL(r_2), DIMENSION(:,:), POINTER :: wbfice !
!jhan:prev. def. in UM only
     REAL, DIMENSION(:), POINTER   :: tggav ! mean soil temperature in K
     REAL, DIMENSION(:), POINTER   :: otgg  ! soil temperature in K
     REAL, DIMENSION(:), POINTER :: otss     ! surface temperature (weighted soil, snow)
     REAL, DIMENSION(:), POINTER :: otss_0     ! surface temperature (weighted soil, snow)
     REAL, DIMENSION(:), POINTER :: tprecip
     REAL, DIMENSION(:), POINTER :: tevap
     REAL, DIMENSION(:), POINTER :: trnoff
     REAL, DIMENSION(:), POINTER :: totenbal
     REAL, DIMENSION(:), POINTER :: totenbal2
     REAL, DIMENSION(:), POINTER :: fland     ! factor for latent heat
     REAL, DIMENSION(:), POINTER :: ifland ! integer soil type
     REAL, DIMENSION(:,:), POINTER :: tilefrac     ! factor for latent heat
     REAL, DIMENSION(:), POINTER :: qasrf ! heat advected to the snow by precip. 
     REAL, DIMENSION(:), POINTER :: qfsrf ! energy of snowpack phase changes 
     REAL, DIMENSION(:), POINTER :: qssrf ! sublimation 
!jhan:}
  END TYPE soil_snow_type

  ! Vegetation parameters:
  TYPE veg_parameter_type
     REAL, DIMENSION(:), POINTER :: canst1 ! max intercepted water by canopy (mm/LAI)
     REAL, DIMENSION(:), POINTER :: dleaf  ! chararacteristc legnth of leaf (m)
     REAL, DIMENSION(:), POINTER :: ejmax  ! max pot. electron transp rate top leaf(mol/m2/s)
     INTEGER,DIMENSION(:), POINTER :: iveg ! vegetation type
     INTEGER,DIMENSION(:), POINTER :: meth ! method for calculation of canopy fluxes and temp.
     REAL, DIMENSION(:), POINTER :: frac4  ! fraction of c4 plants
     REAL, DIMENSION(:), POINTER :: hc     ! roughness height of canopy (veg - snow)
     REAL, DIMENSION(:), POINTER :: vlai   ! leaf area index
     REAL, DIMENSION(:), POINTER :: xalbnir 
     REAL, DIMENSION(:), POINTER :: rp20   ! plant respiration coefficient at 20 C
     REAL, DIMENSION(:), POINTER :: rpcoef ! temperature coef nonleaf plant respiration (1/C)
     REAL, DIMENSION(:), POINTER :: rs20  ! soil respiration at 20 C [mol m-2 s-1]
     REAL, DIMENSION(:), POINTER :: shelrb ! sheltering factor (dimensionless)
     REAL, DIMENSION(:), POINTER :: vegcf  ! kdcorbin, 08/10
     REAL, DIMENSION(:), POINTER :: tminvj ! min temperature of the start of photosynthesis
     REAL, DIMENSION(:), POINTER :: tmaxvj ! max temperature of the start of photosynthesis
     REAL, DIMENSION(:), POINTER :: vbeta  ! 
     REAL, DIMENSION(:), POINTER :: vcmax  ! max RuBP carboxylation rate top leaf (mol/m2/s)
     REAL, DIMENSION(:), POINTER :: xfang  ! leaf angle PARAMETER
    REAL, DIMENSION(:), POINTER :: extkn  ! extinction coef for vertical
    REAL, DIMENSION(:), POINTER :: wai   ! wood area index (stem+branches+twigs)
    LOGICAL,   DIMENSION(:), POINTER :: deciduous ! flag used for phenology fix
     REAL, DIMENSION(:,:), POINTER :: refl
     REAL, DIMENSION(:,:), POINTER :: taul 
    REAL, DIMENSION(:,:), POINTER :: froot  ! fraction of root in each soil layer
!jhan:prev. def. in UM only
     REAL, DIMENSION(:), POINTER :: vlaimax ! ???
  END TYPE veg_parameter_type

  ! Canopy/vegetation variables:
  TYPE canopy_type
       REAL(r_2), DIMENSION(:), POINTER :: fess   ! latent heatfl from soil (W/m2)
     REAL(r_2), DIMENSION(:), POINTER :: fesp   ! latent heatfl from soil (W/m2)
     REAL, DIMENSION(:), POINTER :: cansto ! canopy water storage (mm)
     REAL, DIMENSION(:), POINTER :: cduv  ! drag coefficient for momentum
     REAL, DIMENSION(:), POINTER :: delwc ! change in canopy water store (mm/dels)
     REAL, DIMENSION(:), POINTER :: dewmm ! dewfall (mm)
     REAL(r_2), DIMENSION(:), POINTER :: dgdtg ! derivative of gflux wrt soil temp
     REAL, DIMENSION(:), POINTER :: fe    ! total latent heat (W/m2)
     REAL, DIMENSION(:), POINTER :: fh    ! total sensible heat (W/m2)
     REAL, DIMENSION(:), POINTER :: fpn   ! plant photosynthesis (g C s-1)
     REAL, DIMENSION(:), POINTER :: frp   ! plant respiration (g C m-2 s-1)
     REAL, DIMENSION(:), POINTER :: frpw  ! plant respiration (g C m-2 s-1)???
     REAL, DIMENSION(:), POINTER :: frpr  ! plant respiration (g C m-2 s-1)???
     REAL, DIMENSION(:), POINTER :: frs   ! soil respiration (g C m-2 s-1)
     REAL, DIMENSION(:), POINTER :: fnee  ! net carbon flux (g C m-2 s-1)
     REAL, DIMENSION(:), POINTER :: frday ! daytime leaf resp
     REAL, DIMENSION(:), POINTER :: fnv   ! net rad. avail. to canopy (W/m2)
     REAL, DIMENSION(:), POINTER :: fev   ! latent hf from canopy (W/m2)
     REAL(r_2), DIMENSION(:), POINTER :: fevc  ! dry canopy transpiration (W/m2)
     REAL, DIMENSION(:), POINTER :: fhv   ! sens heatfl from canopy (W/m2)
     REAL, DIMENSION(:), POINTER :: fns   ! net rad avail to soil (W/m2)
     REAL, DIMENSION(:), POINTER :: fhs   ! sensible heat flux from soil
     REAL, DIMENSION(:), POINTER :: fhs_cor 
     REAL, DIMENSION(:), POINTER :: ga    ! ground heat flux (W/m2) ???
     REAL, DIMENSION(:), POINTER :: ghflux  ! ground heat flux (W/m2) ???
     REAL, DIMENSION(:), POINTER :: precis! throughfall to soil, after snow (mm)
     REAL, DIMENSION(:), POINTER :: qscrn ! specific humudity at screen height (g/g)
     REAL, DIMENSION(:), POINTER :: rnet  ! net radiation absorbed by surface (W/m2)
     REAL, DIMENSION(:), POINTER :: segg  ! latent heatfl from soil mm
     REAL, DIMENSION(:), POINTER :: sghflux ! ground heat flux (W/m2) ???
     REAL, DIMENSION(:), POINTER :: through ! canopy throughfall (mm)
     REAL, DIMENSION(:), POINTER :: spill ! can.storage excess after dewfall (mm)
     REAL, DIMENSION(:), POINTER :: tscrn ! air temperature at screen height (oC)
     REAL, DIMENSION(:), POINTER :: wcint ! canopy rainfall interception (mm)
     REAL, DIMENSION(:), POINTER :: tv    ! vegetation temp (K)
     REAL, DIMENSION(:), POINTER :: us    ! friction velocity
     REAL, DIMENSION(:), POINTER :: uscrn ! wind speed at screen height (m/s)
     REAL, DIMENSION(:), POINTER :: vlaiw  ! lai adj for snow depth for calc of resistances
     REAL, DIMENSION(:), POINTER :: rghlai  ! lai adj for snow depth for calc of resistances
     REAL, DIMENSION(:), POINTER :: fwet   ! fraction of canopy wet
     REAL, DIMENSION(:,:), POINTER :: evapfbl
!jhan:prev. def. in UM only
     REAL, DIMENSION(:), POINTER :: epot  ! total potential evaporation 
     REAL, DIMENSION(:), POINTER :: fnpp  ! npp flux
     REAL, DIMENSION(:), POINTER :: fevw_pot ! potential lat heat from canopy
     REAL, DIMENSION(:), POINTER :: gswx_T ! ! stom cond for water
     REAL, DIMENSION(:), POINTER :: gs_vs ! ! stom cond for water
     REAL, DIMENSION(:), POINTER :: cdtq  ! drag coefficient for momentum
     REAL, DIMENSION(:), POINTER :: wetfac_cs ! 
     REAL, DIMENSION(:), POINTER :: fevw  ! lat heat fl wet canopy (W/m2)
     REAL, DIMENSION(:), POINTER :: fhvw  ! sens heatfl from wet canopy (W/m2)
     REAL(r_2), DIMENSION(:), POINTER :: fes   ! latent heatfl from soil (W/m2)
     REAL(r_2), DIMENSION(:), POINTER :: fes_cor   ! latent heatfl from soil (W/m2)
     REAL, DIMENSION(:,:), POINTER :: gswx ! ! stom cond for water
     REAL, DIMENSION(:,:), POINTER :: zetar ! stability correction
     REAL, DIMENSION(:), POINTER :: oldcansto ! canopy water storage (mm)
  END TYPE canopy_type


  ! Radiation variables:
  TYPE radiation_type
     REAL, DIMENSION(:,:), POINTER :: albedo ! canopy+soil albedo
     REAL, DIMENSION(:), POINTER     :: extkb  ! beam radiation extinction coeff
     REAL, DIMENSION(:), POINTER     :: extkd2 ! diffuse 2D radiation extinction coeff
     REAL, DIMENSION(:), POINTER ::  extkd ! diffuse radiation extinction coeff (-)
     REAL, DIMENSION(:), POINTER :: flws   ! soil long-wave radiation
     REAL, DIMENSION(:,:), POINTER  :: fvlai  ! leaf area index of big leaf
     REAL, DIMENSION(:), POINTER :: latitude  ! latitude
     REAL, DIMENSION(:), POINTER :: lwabv ! long wave absorbed by vegetation
     REAL, DIMENSION(:,:,:), POINTER :: qcan ! absorbed radiation for canopy (W/m^2)
     REAL, DIMENSION(:), POINTER     :: qssabs ! absorbed short-wave radiation for soil
     REAL, DIMENSION(:,:), POINTER ::  rhocdf ! canopy diffuse reflectance (-)
     REAL, DIMENSION(:,:), POINTER  :: rniso  !  sum(rad%qcan, 3) total abs by canopy (W/m2)
     REAL, DIMENSION(:,:), POINTER  :: scalex ! scaling PARAMETER for big leaf
     REAL, DIMENSION(:), POINTER     :: transd ! frac SW diffuse transmitted through canopy
     REAL, DIMENSION(:), POINTER ::  trad  !  radiative temperature (soil and veg)
     REAL,DIMENSION(:,:), POINTER  :: reffdf  !effective conopy diffuse reflectance
     REAL,DIMENSION(:,:), POINTER  :: reffbm  !effective conopy beam reflectance
     REAL, DIMENSION(:,:), POINTER :: extkbm  !modified k beam(6.20)(for leaf scattering)
     REAL, DIMENSION(:,:), POINTER :: extkdm  !modified k diffuse(6.20)(for leaf scattering)
     REAL, DIMENSION(:,:), POINTER   :: fbeam   !beam fraction 
     REAL, DIMENSION(:,:), POINTER   :: cexpkbm ! canopy beam transmittance
     REAL, DIMENSION(:,:), POINTER   :: cexpkdm ! canopy diffuse transmittance
     REAL, DIMENSION(:,:), POINTER :: rhocbm  !modified canopy beam reflectance(6.21)
     REAL, DIMENSION(:), POINTER   :: transb  ! fraction SW beam tranmitted through canopy
!jhan:prev. def. in UM only
     REAL, DIMENSION(:), POINTER :: albedo_T ! canopy+soil albedo for VIS+NIR
     REAL, DIMENSION(:,:), POINTER  :: gradis ! radiative conductance
     REAL, DIMENSION(:), POINTER :: longitude ! longitude
     REAL, DIMENSION(:), POINTER     :: workp1 ! absorbed short-wave radiation for soil
     REAL, DIMENSION(:), POINTER     :: workp2 ! absorbed short-wave radiation for soil
     REAL, DIMENSION(:), POINTER     :: workp3 ! absorbed short-wave radiation for soil
!    REAL(r_2), DIMENSION(:,:), POINTER  :: gradis ! radiative conductance
  END TYPE radiation_type

  ! Roughness variables:
  TYPE roughness_type
     ! "coexp": coefficient in exponential in-canopy wind profile
     ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
     ! canopy and roughness-sublayer U(z) at z=h
     REAL, DIMENSION(:), POINTER   :: coexp ! Extinction coef for wind profile in canopy
     REAL, DIMENSION(:), POINTER   :: disp  ! zero-plane displacement
     REAL, DIMENSION(:), POINTER   :: hruff ! canopy height above snow level
     REAL, DIMENSION(:), POINTER :: hruff_grmx ! max ht of canopy from tiles on same grid 
     REAL, DIMENSION(:), POINTER   :: rt0us ! eq. 3.54, SCAM manual (CSIRO tech report 132)
     REAL, DIMENSION(:), POINTER   :: rt1usa ! resistance from disp to hruf
     REAL, DIMENSION(:), POINTER   :: rt1usb ! resist fr hruf to zruffs (zref if zref<zruffs)
     REAL, DIMENSION(:), POINTER   :: rt1 ! 1/aerodynamic conductance
     REAL, DIMENSION(:), POINTER   :: term2, term3, term5, term6 ! for aerodyn resist. calc.
     ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
     REAL, DIMENSION(:), POINTER   :: usuh ! Friction velocity/windspeed at canopy height
     REAL, DIMENSION(:), POINTER   :: za_uv   ! level of lowest atmospheric model layer
     REAL, DIMENSION(:), POINTER   :: za_tq   ! level of lowest atmospheric model layer
     REAL, DIMENSION(:), POINTER   :: z0m  ! roughness length
     REAL, DIMENSION(:), POINTER   :: zref_uv ! Reference height for met forcing
     REAL, DIMENSION(:), POINTER   :: zref_tq ! Reference height for met forcing
     REAL, DIMENSION(:), POINTER   :: zruffs ! SCALAR Roughness sublayer depth (ground=origin)
     REAL, DIMENSION(:), POINTER   :: z0soilsn ! roughness length of bare soil surface
     REAL, DIMENSION(:), POINTER   :: z0soil ! roughness length of bare soil surface
  END TYPE roughness_type

  ! Air variables:
  TYPE air_type
    REAL, DIMENSION(:), POINTER :: rho  ! dry air density (kg m-3)
    REAL, DIMENSION(:), POINTER :: volm ! molar volume (m3 mol-1)
    REAL, DIMENSION(:), POINTER :: rlam ! latent heat for water (j/kg)
    REAL, DIMENSION(:), POINTER :: qsat ! saturation specific humidity
    REAL, DIMENSION(:), POINTER :: epsi ! d(qsat)/dT ((kg/kg)/K)
    REAL, DIMENSION(:), POINTER :: visc ! air kinematic viscosity (m2/s)
    REAL, DIMENSION(:), POINTER :: psyc ! psychrometric constant
    REAL, DIMENSION(:), POINTER :: dsatdk ! d(es)/dT (mb/K)
    REAL, DIMENSION(:), POINTER :: cmolar ! conv. from m/s to mol/m2/s
  END TYPE air_type


  ! Meterological data:
  TYPE met_type
     REAL, DIMENSION(:), POINTER :: ca   ! CO2 concentration (mol/mol)
     INTEGER, DIMENSION(:), POINTER :: year ! local time year AD 
     INTEGER, DIMENSION(:), POINTER :: moy  ! local time month of year 
     REAL, DIMENSION(:), POINTER :: doy  ! local time day of year = days since 0 hr 1st Jan 
     REAL, DIMENSION(:), POINTER :: hod  ! local hour of day
     REAL, DIMENSION(:,:), POINTER :: fsd  ! downward short-wave radiation (W/m2)
     REAL, DIMENSION(:), POINTER :: ofsd  ! downward short-wave radiation (W/m2)
     REAL, DIMENSION(:), POINTER :: fld  ! downward long-wave radiation (W/m2)
     REAL, DIMENSION(:), POINTER :: precip  ! rainfall (liquid+solid)(mm/dels)
     REAL, DIMENSION(:), POINTER :: precip_sn ! solid preipitation only (mm/dels)
     REAL, DIMENSION(:), POINTER :: tk      ! surface air temperature (oK)
     REAL, DIMENSION(:), POINTER :: tvair   ! within canopy air temperature (oK)
     REAL, DIMENSION(:), POINTER :: tvrad   ! radiative vegetation temperature (K)
     REAL, DIMENSION(:), POINTER :: pmb     ! surface air pressure (mbar)
     REAL, DIMENSION(:), POINTER :: ua      ! surface wind speed (m/s)
     REAL, DIMENSION(:), POINTER :: qv      ! surface specific humidity (g/g)
     REAL, DIMENSION(:), POINTER :: qvair   ! within canopy specific humidity (g/g)
     REAL, DIMENSION(:), POINTER :: da      ! water vap pressure deficit at ref height (Pa)
     REAL, DIMENSION(:), POINTER :: dva     ! in canopy water vap pressure deficit (Pa)
     REAL, DIMENSION(:), POINTER :: coszen  ! cos(zenith angle of sun)
!jhan:prev. def. in UM only
     REAL, DIMENSION(:), POINTER :: tc      ! surface air temperature (oC)
  END TYPE met_type

  ! Cumulative flux variables:
  TYPE sum_flux_type
     REAL, DIMENSION(:), POINTER   :: sumpn  ! sum of canopy photosynthesis (g C m-2)
     REAL, DIMENSION(:), POINTER   :: sumrp  ! sum of plant respiration (g C m-2)
     REAL, DIMENSION(:), POINTER   :: sumrpw ! sum of plant respiration (g C m-2)
     REAL, DIMENSION(:), POINTER   :: sumrpr ! sum of plant respiration (g C m-2)
     REAL, DIMENSION(:), POINTER   :: sumrs  ! sum of soil respiration (g C m-2)
     REAL, DIMENSION(:), POINTER   :: sumrd  ! sum of daytime respiration (g C m-2)
     REAL, DIMENSION(:), POINTER   :: dsumpn ! daily sumpn
     REAL, DIMENSION(:), POINTER   :: dsumrp ! daily sumrp
     REAL, DIMENSION(:), POINTER   :: dsumrs ! daily sumrs
     REAL, DIMENSION(:), POINTER   :: dsumrd ! daily sumrd
     REAL, DIMENSION(:), POINTER   :: sumxrp ! sum plant resp. modifier
     REAL, DIMENSION(:), POINTER   :: sumxrs ! sum soil resp. modifier
  END TYPE sum_flux_type

  TYPE bgc_pool_type
     REAL, DIMENSION(:,:), POINTER :: cplant ! plant carbon (g C/m2))
     REAL, DIMENSION(:,:), POINTER :: csoil  ! soil carbon (g C/m2)
     REAL, DIMENSION(ncp)  :: ratecp ! plant carbon rate constant (1/year)
     REAL, DIMENSION(ncs)  :: ratecs ! soil carbon rate constant (1/year)
  END TYPE bgc_pool_type

  ! Functions for allocating these types
  ! All overloaded so code only needs to call alloc_cbm_var
  ! Alloc routines could all initialise to NaN or zero for debugging?
  ! Don't need the mp argument here as it's a module variable.
  PUBLIC :: alloc_cbm_var
  PRIVATE :: alloc_bgc_pool_type, dealloc_bgc_pool_type
  INTERFACE alloc_cbm_var
     MODULE PROCEDURE alloc_balances_type,            &
          alloc_soil_parameter_type,      &
          alloc_soil_snow_type,           &
          alloc_veg_parameter_type,       &
          alloc_canopy_type,              &
          alloc_radiation_type,           &
          alloc_roughness_type,           &
          alloc_air_type,                 &
          alloc_met_type,                 &
          alloc_sum_flux_type,            &
          alloc_bgc_pool_type            
  END INTERFACE
  INTERFACE dealloc_cbm_var
     MODULE PROCEDURE dealloc_balances_type,            &
          dealloc_soil_parameter_type,      &
          dealloc_soil_snow_type,           &
          dealloc_veg_parameter_type,       &
          dealloc_canopy_type,              &
          dealloc_radiation_type,           &
          dealloc_roughness_type,           &
          dealloc_air_type,                 &
          dealloc_met_type,                 &
          dealloc_sum_flux_type,            &
          dealloc_bgc_pool_type            
  END INTERFACE
CONTAINS
  
   SUBROUTINE alloc_balances_type(var, mp)
      TYPE(balances_type), INTENT(inout) :: var
      INTEGER, INTENT(in) :: mp
         allocate( var% drybal(mp) ) 
         allocate( var% ebal(mp) )  
         allocate( var% ebal_tot(mp) )
         allocate( var% ebaltr(mp) )  
         allocate( var% ebal_tottr(mp) )
         allocate( var% ebal_cncheck(mp) )  
         allocate( var% ebal_tot_cncheck(mp) )
         allocate( var% evap_tot(mp) )
         allocate( var% osnowd0(mp) )
         allocate( var% precip_tot(mp) )
         allocate( var% rnoff_tot(mp) )
         allocate( var% wbal(mp) )   
         allocate( var% wbal_tot(mp) )
         allocate( var% wbtot0(mp) ) 
         allocate( var% wetbal(mp) )
         allocate( var% cansto0(mp) ) 
         allocate( var% evapc_tot(mp) ) 
         allocate( var% evaps_tot(mp) ) 
         allocate( var% rnof1_tot(mp) ) 
         allocate( var% rnof2_tot(mp) ) 
         allocate( var% snowdc_tot(mp) )
         allocate( var% wbal_tot1(mp) ) 
         allocate( var% owbtot(mp) ) 
         allocate( var% delwc_tot(mp) ) 
         allocate( var% qasrf_tot(mp) )
         allocate( var% qfsrf_tot(mp) ) 
         allocate( var% qssrf_tot(mp) ) 
      return
   END SUBROUTINE alloc_balances_type

   SUBROUTINE alloc_soil_parameter_type(var, mp)
      TYPE(soil_parameter_type), INTENT(inout) :: var
      INTEGER, INTENT(in) :: mp
         allocate( var% bch(mp) )   
         allocate( var% c3(mp) )    
         allocate( var% clay(mp) )  
         allocate( var% css(mp) )   
         allocate( var% hsbh(mp) )  
         allocate( var% hyds(mp) )  
         allocate( var% i2bp3(mp) ) 
         allocate( var% ibp2(mp) )  
         allocate( var% isoilm(mp) )  
         allocate( var% rhosoil(mp) )  
         !allocate( var% rs20(mp) )   
         allocate( var% sand(mp) )   
         allocate( var% sfc(mp) )   
         allocate( var% silt(mp) )   
         allocate( var% ssat(mp) )   
         allocate( var% sucs(mp) )   
         allocate( var% swilt(mp) )  
         allocate( var% zse(ms) )    
         allocate( var% zshh(ms+1) )  
         allocate( var% cnsd(mp) )  
         allocate( var% albsoil(mp, nrb) )  
         allocate( var% pwb_min(mp) )  
         allocate( var% albsoilf(mp) )  
      return
   END SUBROUTINE alloc_soil_parameter_type
 
   SUBROUTINE alloc_soil_snow_type(var, mp)
      TYPE(soil_snow_type), INTENT(inout) :: var
      INTEGER, INTENT(in) :: mp
          ALLOCATE ( var % iantrct(mp) )
    ALLOCATE ( var % pudsto(mp) )
    ALLOCATE ( var % pudsmx(mp) )
    ALLOCATE ( var % dtmlt(mp,3) )
         allocate( var% albsoilsn(mp,nrb) ) 
         allocate( var% cls(mp) )     
         allocate( var% dfn_dtg(mp) ) 
         allocate( var% dfh_dtg(mp) ) 
         allocate( var% dfe_ddq(mp) ) 
         allocate( var% ddq_dtg(mp) ) 
         allocate( var% evapsn(mp) )  
         allocate( var% fwtop(mp) )   
         allocate( var% fwtop1(mp) )   
         allocate( var% fwtop2(mp) )   
         allocate( var% fwtop3(mp) )   
         allocate( var% gammzz(mp,ms) ) 
         allocate( var% isflag(mp) ) 
         allocate( var% osnowd(mp) ) 
         allocate( var% potev(mp) ) 
         allocate( var% runoff(mp) )
         allocate( var% rnof1(mp) ) 
         allocate( var% rnof2(mp) ) 
         allocate( var% rtsoil(mp) )
         allocate( var% sconds(mp,msn) ) 
         allocate( var% sdepth(mp,msn) ) 
         allocate( var% smass(mp,msn) ) 
         allocate( var% snage(mp) )  
         allocate( var% snowd(mp) )  
         allocate( var% smelt(mp) )  
         allocate( var% ssdn(mp,msn) ) 
         allocate( var% ssdnn(mp) ) 
         allocate( var% tgg(mp,ms) )   
         allocate( var% tggsn(mp,msn) ) 
         allocate( var% tss(mp) )   
         allocate( var% tss_p(mp) )   
         allocate( var% deltss(mp) )   
         allocate( var% owb1(mp) )   
         allocate( var% wb(mp,ms) )    
         allocate( var% wbice(mp,ms) ) 
         allocate( var% wblf(mp,ms) ) 
         allocate( var%wbtot(mp) )    
         allocate( var%wbtot1(mp) )    
         allocate( var%wbtot2(mp) )    
         allocate( var%wb_lake(mp) )    
         allocate( var%sinfil(mp) )    
         allocate( var%evapfbl(mp,ms) )    
         allocate( var%qstss(mp) )    
         allocate( var%wetfac(mp) )  
         allocate( var%owetfac(mp) )  
         allocate( var%t_snwlr(mp) )  
         allocate( var%wbfice(mp,ms) )  
         allocate( var%tggav(mp) )  
         allocate( var%otgg(mp) )   
         allocate( var%otss(mp) )   
         allocate( var%otss_0(mp) )   
         allocate( var%tprecip(mp) ) 
         allocate( var%tevap(mp) ) 
         allocate( var%trnoff(mp) ) 
         allocate( var%totenbal(mp) ) 
         allocate( var%totenbal2(mp) ) 
         allocate( var%fland(mp) )      
         allocate( var%ifland(mp) )  
         allocate( var%tilefrac(mp,n_tiles) ) 
         allocate( var%qasrf(mp) )  
         allocate( var%qfsrf(mp) )  
         allocate( var%qssrf(mp) )  
  END SUBROUTINE alloc_soil_snow_type
   
   SUBROUTINE alloc_veg_parameter_type(var, mp)
      TYPE(veg_parameter_type), INTENT(inout) :: var
      INTEGER, INTENT(in) :: mp
         allocate( var% canst1(mp) ) 
         allocate( var% dleaf(mp) )  
         allocate( var% ejmax(mp) ) 
         allocate( var% iveg(mp) ) 
         allocate( var% meth(mp) ) 
         allocate( var% frac4(mp) )  
         allocate( var% hc(mp) )     
         allocate( var% vlai(mp) )   
         allocate( var% xalbnir(mp) ) 
         allocate( var% rp20(mp) )   
         allocate( var% rpcoef(mp) ) 
         allocate( var% rs20(mp) )   
         allocate( var% shelrb(mp) ) 
         allocate( var% vegcf(mp) )  
         allocate( var% tminvj(mp) ) 
         allocate( var% tmaxvj(mp) ) 
         allocate( var% vbeta(mp) )  
         allocate( var% vcmax(mp) )  
         allocate( var% xfang(mp) )  
         allocate( var%extkn(mp) ) 
         allocate( var%wai(mp) )   
         allocate( var%deciduous(mp) ) 
         allocate( var%froot(mp,ms) ) 
         allocate( var%refl(mp,2) ) !jhan:swb?
         allocate( var%taul(mp,2) ) 
         allocate( var%vlaimax(mp) ) 
      return      
   end SUBROUTINE alloc_veg_parameter_type
   
   SUBROUTINE alloc_canopy_type(var, mp)
      USE physical_constants, only : niter
      TYPE(canopy_type), INTENT(inout) :: var
      INTEGER, INTENT(in) :: mp
          ALLOCATE ( var % fess(mp) )
    ALLOCATE ( var % fesp(mp) )
         allocate( var% cansto(mp) )  
         allocate( var% cduv(mp) )   
         allocate( var% delwc(mp) )  
         allocate( var% dewmm(mp) )  
         allocate( var% dgdtg(mp) )  
         allocate( var% fe(mp) )      
         allocate( var% fh(mp) )      
         allocate( var% fpn(mp) )     
         allocate( var% frp(mp) )     
         allocate( var% frpw(mp) )    
         allocate( var% frpr(mp) )    
         allocate( var% frs(mp) )     
         allocate( var% fnee(mp) )    
         allocate( var% frday(mp) )   
         allocate( var% fnv(mp) )     
         allocate( var% fev(mp) )     
         allocate( var% fevc(mp) )    
         allocate( var% fhv(mp) )     
         allocate( var% fns(mp) )     
         allocate( var% fhs(mp) )     
         allocate( var% fhs_cor(mp) )     
         allocate( var% ga(mp) )      
         allocate( var% ghflux(mp) )   
         allocate( var% precis(mp) ) 
         allocate( var% qscrn(mp) )  
         allocate( var% rnet(mp) )   
         allocate( var% segg(mp) )   
         allocate( var% sghflux(mp) )  
         allocate( var% through(mp) )  
         allocate( var% spill(mp) )  
         allocate( var% tscrn(mp) )  
         allocate( var% wcint(mp) )  
         allocate( var% tv(mp) )      
         allocate( var% us(mp) )      
         allocate( var% uscrn(mp) )   
         allocate( var% rghlai(mp) ) 
         allocate( var% vlaiw(mp) ) 
         allocate( var% fwet(mp) )   
    ALLOCATE ( var % evapfbl(mp,ms) )
         allocate( var% epot(mp) )   
         allocate( var% fnpp(mp) )   
         allocate( var% fevw_pot(mp) )  
         allocate( var% gswx_T(mp) )  
         allocate( var% cdtq(mp) )   
         allocate( var% wetfac_cs(mp) )  
         allocate( var% fevw(mp) )   
         allocate( var% fhvw(mp) )   
         allocate( var% fes(mp) )    
         allocate( var% fes_cor(mp) )    
         allocate( var% gswx(mp,mf) )  
         allocate( var% oldcansto(mp) )  
         allocate( var% zetar(mp,niter) )  
      return
   END SUBROUTINE alloc_canopy_type
   
  SUBROUTINE alloc_radiation_type(var, mp)
    TYPE(radiation_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
         allocate( var% albedo(mp,nrb) ) 
         allocate( var% extkb(mp) )  
         allocate( var% extkd2(mp) )
         allocate( var% extkd(mp) )
         allocate( var% flws(mp) )
         allocate( var% fvlai(mp,mf) )
         allocate( var% latitude(mp) )
         allocate( var% lwabv(mp) )
         allocate( var% qcan(mp,mf,nrb) )
         allocate( var% qssabs(mp) )
         allocate( var% rhocdf(mp,nrb) )
         allocate( var% rniso(mp,mf) )
         allocate( var% scalex(mp,mf) )
         allocate( var% transd(mp) )
         allocate( var% trad(mp) )
         allocate( var% reffdf(mp,nrb) )
         allocate( var% reffbm(mp,nrb) )
         allocate( var% extkbm(mp,nrb) )
         allocate( var% extkdm(mp,nrb) )
         allocate( var% cexpkbm(mp,swb) )
         allocate( var% cexpkdm(mp,swb) )
         allocate( var% fbeam(mp,nrb) )
         allocate( var% rhocbm(mp,nrb) )
         allocate( var% transb(mp) )
         allocate( var% albedo_T(mp) )
         allocate( var% gradis(mp,mf) )
         allocate( var% longitude(mp) )
         allocate( var% workp1(mp) )
         allocate( var% workp2(mp) )
         allocate( var% workp3(mp) )
      return
   END SUBROUTINE alloc_radiation_type
   
  SUBROUTINE alloc_roughness_type(var, mp)
    TYPE(roughness_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp

    ALLOCATE ( var % coexp(mp) )
    ALLOCATE ( var % disp(mp) )
    ALLOCATE ( var % hruff(mp) )
    ALLOCATE ( var % hruff_grmx(mp) )
    ALLOCATE ( var % rt0us(mp) )
    ALLOCATE ( var % rt1usa(mp) )
    ALLOCATE ( var % rt1usb(mp) )
    ALLOCATE ( var % rt1(mp) )
    ALLOCATE ( var % term2(mp) )
    ALLOCATE ( var % term3(mp) )
    ALLOCATE ( var % term5(mp) )
    ALLOCATE ( var % term6(mp) )
    ALLOCATE ( var % usuh(mp) )
    ALLOCATE ( var % za_uv(mp) )
    ALLOCATE ( var % za_tq(mp) )
    ALLOCATE ( var % z0m(mp) )
    ALLOCATE ( var % zref_uv(mp) )
    ALLOCATE ( var % zref_tq(mp) )
    ALLOCATE ( var % zruffs(mp) )
    ALLOCATE ( var % z0soilsn(mp) )
    ALLOCATE ( var % z0soil(mp) )
  END SUBROUTINE alloc_roughness_type
   
  SUBROUTINE alloc_air_type(var, mp)
    TYPE(air_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % rho(mp) )
    ALLOCATE ( var % volm(mp) )
    ALLOCATE ( var % rlam(mp) )
    ALLOCATE ( var % qsat(mp) )
    ALLOCATE ( var % epsi(mp) )
    ALLOCATE ( var % visc(mp) )
    ALLOCATE ( var % psyc(mp) )
    ALLOCATE ( var % dsatdk(mp) )
    ALLOCATE ( var % cmolar(mp) )
  END SUBROUTINE alloc_air_type
 
  
  SUBROUTINE alloc_met_type(var, mp)
    TYPE(met_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % ca(mp) )
    ALLOCATE ( var % year(mp) )
    ALLOCATE ( var % moy(mp) )
    ALLOCATE ( var % doy(mp) )
    ALLOCATE ( var % hod(mp) )
    ALLOCATE ( var % fsd(mp,swb) ) 
    ALLOCATE ( var % ofsd(mp) ) 
    ALLOCATE ( var % fld(mp) )
    ALLOCATE ( var % precip(mp) )
    ALLOCATE ( var % precip_sn(mp) )
    ALLOCATE ( var % tk(mp) )
    ALLOCATE ( var % tvair(mp) )
    ALLOCATE ( var % tvrad(mp) )
    ALLOCATE ( var % pmb(mp) )
    ALLOCATE ( var % ua(mp) )
    ALLOCATE ( var % qv(mp) )
    ALLOCATE ( var % qvair(mp) )
    ALLOCATE ( var % da(mp) )
    ALLOCATE ( var % dva(mp) )
    ALLOCATE ( var % coszen(mp) )
    !ALLOCATE ( var % tc(mp) )
  END SUBROUTINE alloc_met_type
   
  SUBROUTINE alloc_sum_flux_type(var, mp)
    TYPE(sum_flux_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % sumpn(mp) )
    ALLOCATE ( var % sumrp(mp) )
    ALLOCATE ( var % sumrpw(mp) )
    ALLOCATE ( var % sumrpr(mp) )
    ALLOCATE ( var % sumrs(mp) )
    ALLOCATE ( var % sumrd(mp) )
    ALLOCATE ( var % dsumpn(mp) )
    ALLOCATE ( var % dsumrp(mp) )
    ALLOCATE ( var % dsumrs(mp) )
    ALLOCATE ( var % dsumrd(mp) )
    ALLOCATE ( var % sumxrp(mp) )
    ALLOCATE ( var % sumxrs(mp) )
  END SUBROUTINE alloc_sum_flux_type

  SUBROUTINE alloc_bgc_pool_type(var, mp)
    TYPE(bgc_pool_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    ALLOCATE ( var % cplant(mp,ncp) )
    ALLOCATE ( var % csoil(mp,ncs) )
  END SUBROUTINE alloc_bgc_pool_type

  ! Begin deallocation routines:
   SUBROUTINE dealloc_balances_type(var, mp)
    TYPE(balances_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
         deallocate( var% drybal ) 
         deallocate( var% ebal )  
         deallocate( var% ebal_tot )
         deallocate( var% ebaltr )  
         deallocate( var% ebal_tottr )
         deallocate( var% ebal_cncheck )  
         deallocate( var% ebal_tot_cncheck )
         deallocate( var% evap_tot)
         deallocate( var% osnowd0 )
         deallocate( var% precip_tot )
         deallocate( var% rnoff_tot )
         deallocate( var% wbal )   
         deallocate( var% wbal_tot )
         deallocate( var% wbtot0 ) 
         deallocate( var% wetbal )
         deallocate( var% cansto0 ) 
         deallocate( var% evapc_tot ) 
         deallocate( var% evaps_tot ) 
         deallocate( var% rnof1_tot ) 
         deallocate( var% rnof2_tot ) 
         deallocate( var% snowdc_tot )
         deallocate( var% wbal_tot1 ) 
         deallocate( var% owbtot ) 
         deallocate( var% delwc_tot ) 
         deallocate( var% qasrf_tot )
         deallocate( var% qfsrf_tot ) 
         deallocate( var% qssrf_tot ) 
      return
  END SUBROUTINE dealloc_balances_type

  SUBROUTINE dealloc_soil_parameter_type(var, mp)
    TYPE(soil_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
         deallocate( var% bch )   
         deallocate( var% c3 )    
         deallocate( var% clay )  
         deallocate( var% css )   
         deallocate( var% hsbh )  
         deallocate( var% hyds )  
         deallocate( var% i2bp3 ) 
         deallocate( var% ibp2 )  
         deallocate( var% isoilm )  
         deallocate( var% rhosoil )  
         !deallocate( var% rs20 )   
         deallocate( var% sand )   
         deallocate( var% sfc )   
         deallocate( var% silt )   
         deallocate( var% ssat )   
         deallocate( var% sucs )   
         deallocate( var% swilt )  
         deallocate( var% zse )    
         deallocate( var% zshh )  
         deallocate( var% cnsd )  
         deallocate( var% albsoil )  
         deallocate( var% cnsd )  
         deallocate( var% pwb_min)  
         deallocate( var% albsoilf )  
      return
  END SUBROUTINE dealloc_soil_parameter_type
 
  SUBROUTINE dealloc_soil_snow_type(var, mp)
    TYPE(soil_snow_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
          deALLOCATE ( var % iantrct )
    deALLOCATE ( var % pudsto )
    deALLOCATE ( var % pudsmx )
    deALLOCATE ( var % dtmlt )

         deallocate( var% albsoilsn ) 
         deallocate( var% cls )     
         deallocate( var% dfn_dtg ) 
         deallocate( var% dfh_dtg ) 
         deallocate( var% dfe_ddq ) 
         deallocate( var% ddq_dtg ) 
         deallocate( var% evapsn )  
         deallocate( var% fwtop )   
         deallocate( var% fwtop1 )   
         deallocate( var% fwtop2 )   
         deallocate( var% fwtop3 )   
         deallocate( var% gammzz ) 
         deallocate( var% isflag ) 
         deallocate( var% osnowd ) 
         deallocate( var% potev ) 
         deallocate( var% runoff )
         deallocate( var% rnof1 ) 
         deallocate( var% rnof2 ) 
         deallocate( var% rtsoil )
         deallocate( var% sconds ) 
         deallocate( var% sdepth ) 
         deallocate( var% smass ) 
         deallocate( var% snage )  
         deallocate( var% snowd )  
         deallocate( var% smelt )  
         deallocate( var% ssdn ) 
         deallocate( var% ssdnn ) 
         deallocate( var% tgg )   
         deallocate( var% tggsn ) 
         deallocate( var% tss )   
         deallocate( var% tss_p )   
         deallocate( var% deltss )   
         deallocate( var% owb1 )   
         deallocate( var% wb )    
         deallocate( var% wbice ) 
         deallocate( var% wblf ) 
         deallocate( var%wbtot )    
         deallocate( var%wbtot1 )    
         deallocate( var%wbtot2 )    
         deallocate( var%wb_lake )    
         deallocate( var%sinfil )    
         deallocate( var%evapfbl)    
         deallocate( var%qstss)    
         deallocate( var%wetfac )  
         deallocate( var%owetfac )  
         deallocate( var%t_snwlr )  
         deallocate( var%wbfice )  
         deallocate( var%tggav )  
         deallocate( var%otgg )   
         deallocate( var%otss )   
         deallocate( var%otss_0 )   
         deallocate( var%tprecip ) 
         deallocate( var%tevap ) 
         deallocate( var%trnoff ) 
         deallocate( var%totenbal ) 
         deallocate( var%totenbal2 ) 
         deallocate( var%fland )      
         deallocate( var%ifland )  
         deallocate( var%tilefrac ) 
         deallocate( var%qasrf )  
         deallocate( var%qfsrf )  
         deallocate( var%qssrf )  
      return
   END SUBROUTINE dealloc_soil_snow_type
   
  SUBROUTINE dealloc_veg_parameter_type(var, mp)
    TYPE(veg_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
         deallocate( var% canst1 ) 
         deallocate( var% dleaf )  
         deallocate( var% ejmax ) 
         deallocate( var% iveg ) 
         deallocate( var% meth ) 
         deallocate( var% frac4 )  
         deallocate( var% hc )     
         deallocate( var% vlai )   
         deallocate( var% xalbnir ) 
         deallocate( var% rp20 )   
         deallocate( var% rpcoef ) 
         deallocate( var% rs20 )   
         deallocate( var% shelrb ) 
         deallocate( var% vegcf )  
         deallocate( var% tminvj ) 
         deallocate( var% tmaxvj ) 
         deallocate( var% vbeta)  
         deallocate( var% vcmax )  
         deallocate( var% xfang )  
         deallocate( var%extkn ) 
         deallocate( var%wai )   
         deallocate( var%deciduous ) 
         deallocate( var%froot) 
         deallocate( var%refl )
         deallocate( var%taul ) 
      return      
  END SUBROUTINE dealloc_veg_parameter_type
   
  SUBROUTINE dealloc_canopy_type(var, mp)
    TYPE(canopy_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
          deALLOCATE ( var % fess )
       deALLOCATE ( var % fesp )
         deallocate( var% cansto )  
         deallocate( var% cduv )   
         deallocate( var% delwc )  
         deallocate( var% dewmm )  
         deallocate( var% dgdtg )  
         deallocate( var% fe )      
         deallocate( var% fh )      
         deallocate( var% fpn )     
         deallocate( var% frp )     
         deallocate( var% frpw )    
         deallocate( var% frpr )    
         deallocate( var% frs )     
         deallocate( var% fnee )    
         deallocate( var% frday )   
         deallocate( var% fnv )     
         deallocate( var% fev )     
         deallocate( var% fevc )    
         deallocate( var% fhv )     
         deallocate( var% fns )     
         deallocate( var% fhs )     
         deallocate( var% fhs_cor )     
         deallocate( var% ga )      
         deallocate( var% ghflux )   
         deallocate( var% precis ) 
         deallocate( var% qscrn )  
         deallocate( var% rnet )   
         deallocate( var% segg )   
         deallocate( var% sghflux )  
         deallocate( var% through )  
         deallocate( var% spill )  
         deallocate( var% tscrn )  
         deallocate( var% wcint )  
         deallocate( var% tv )      
         deallocate( var% us )      
         deallocate( var% uscrn )   
         deallocate( var% rghlai ) 
         deallocate( var% vlaiw ) 
         deallocate( var% fwet )   
    DEALLOCATE ( var % evapfbl )
         deallocate( var% epot )   
         deallocate( var% fnpp )   
         deallocate( var% fevw_pot )  
         deallocate( var% gswx_T )  
         deallocate( var% cdtq )   
         deallocate( var% wetfac_cs )  
         deallocate( var% fevw )   
         deallocate( var% fhvw )   
         deallocate( var% fes )    
         deallocate( var% fes_cor )    
!
         deallocate( var% gswx )  
         deallocate( var% oldcansto )  
         deallocate( var% zetar )  
      return
  END SUBROUTINE dealloc_canopy_type
   
  SUBROUTINE dealloc_radiation_type(var, mp)
    TYPE(radiation_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
         deallocate( var% albedo ) 
         deallocate( var% extkb )  
         deallocate( var% extkd2 )
         deallocate( var% extkd )
         deallocate( var% flws )
         deallocate( var% fvlai )
         deallocate( var% latitude )
         deallocate( var% lwabv )
         deallocate( var% qcan )
         deallocate( var% qssabs )
         deallocate( var% rhocdf )
         deallocate( var% rniso )
         deallocate( var% scalex )
         deallocate( var% transd )
         deallocate( var% trad )
         deallocate( var% reffdf )
         deallocate( var% reffbm )
         deallocate( var% extkbm )
         deallocate( var% extkdm )
         deallocate( var% fbeam )
         deallocate( var% cexpkbm )
         deallocate( var% cexpkdm )
         deallocate( var% rhocbm )
         deallocate( var% transb )
         deallocate( var% albedo_T )
         deallocate( var% gradis )
         deallocate( var% longitude )
         deallocate( var% workp1 )
         deallocate( var% workp2 )
         deallocate( var% workp3 )
      return
  END SUBROUTINE dealloc_radiation_type
   
  SUBROUTINE dealloc_roughness_type(var, mp)
    TYPE(roughness_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % coexp )
    DEALLOCATE ( var % disp )
    DEALLOCATE ( var % hruff )
    DEALLOCATE ( var % hruff_grmx )
    DEALLOCATE ( var % rt0us )
    DEALLOCATE ( var % rt1usa )
    DEALLOCATE ( var % rt1usb )
    DEALLOCATE ( var % rt1 )
    DEALLOCATE ( var % term2 )
    DEALLOCATE ( var % term3 )
    DEALLOCATE ( var % term5 )
    DEALLOCATE ( var % term6 )
    DEALLOCATE ( var % usuh )
    DEALLOCATE ( var % za_uv )
    DEALLOCATE ( var % za_tq )
    DEALLOCATE ( var % z0m )
    DEALLOCATE ( var % zref_uv )
    DEALLOCATE ( var % zref_tq )
    DEALLOCATE ( var % zruffs )
    DEALLOCATE ( var % z0soilsn )
    DEALLOCATE ( var % z0soil )
  END SUBROUTINE dealloc_roughness_type
   
  SUBROUTINE dealloc_air_type(var, mp)
    TYPE(air_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % rho )
    DEALLOCATE ( var % volm )
    DEALLOCATE ( var % rlam )
    DEALLOCATE ( var % qsat )
    DEALLOCATE ( var % epsi )
    DEALLOCATE ( var % visc )
    DEALLOCATE ( var % psyc )
    DEALLOCATE ( var % dsatdk )
    DEALLOCATE ( var % cmolar )
  END SUBROUTINE dealloc_air_type
   
  SUBROUTINE dealloc_met_type(var, mp)
    TYPE(met_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % ca )
    DEALLOCATE ( var % year )
    DEALLOCATE ( var % moy )
    DEALLOCATE ( var % doy )
    DEALLOCATE ( var % hod )
    DEALLOCATE ( var % fsd )
    DEALLOCATE ( var % ofsd )
    DEALLOCATE ( var % fld )
    DEALLOCATE ( var % precip )
    DEALLOCATE ( var % precip_sn )
    DEALLOCATE ( var % tk )
    DEALLOCATE ( var % tvair )
    DEALLOCATE ( var % tvrad )
    DEALLOCATE ( var % pmb )
    DEALLOCATE ( var % ua )
    DEALLOCATE ( var % qv )
    DEALLOCATE ( var % qvair )
    DEALLOCATE ( var % da )
    DEALLOCATE ( var % dva )
    DEALLOCATE ( var % coszen )
  END SUBROUTINE dealloc_met_type
   
  SUBROUTINE dealloc_sum_flux_type(var, mp)
    TYPE(sum_flux_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % sumpn )
    DEALLOCATE ( var % sumrp )
    DEALLOCATE ( var % sumrpw )
    DEALLOCATE ( var % sumrpr )
    DEALLOCATE ( var % sumrs )
    DEALLOCATE ( var % sumrd )
    DEALLOCATE ( var % dsumpn )
    DEALLOCATE ( var % dsumrp )
    DEALLOCATE ( var % dsumrs )
    DEALLOCATE ( var % dsumrd )
    DEALLOCATE ( var % sumxrp )
    DEALLOCATE ( var % sumxrs )
  END SUBROUTINE dealloc_sum_flux_type

  SUBROUTINE dealloc_bgc_pool_type(var, mp)
    TYPE(bgc_pool_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: mp
    DEALLOCATE ( var % cplant )
    DEALLOCATE ( var % csoil )
  END SUBROUTINE dealloc_bgc_pool_type
  

END MODULE define_types
