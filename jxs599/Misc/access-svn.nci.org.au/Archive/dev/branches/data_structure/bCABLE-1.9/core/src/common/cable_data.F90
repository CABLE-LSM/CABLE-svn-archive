module cable_data_module 
   use define_dimensions, only : nrb 
   implicit none
   
   !public :: cable, const
   !public :: air_in, air_out 
   !private

!   integer, parameter :: i_d = KIND(9)
!   integer, parameter :: r_1  = KIND(1.0) 
!   !---at least 10-digit precision
!   integer, parameter :: r_2  = SELECTED_REAL_KIND(12, 50)
!
!   integer, parameter :: n_tiles = 17        ! # possible no of different tiles/patches
!   integer :: mp                             ! # total no of patches/tiles 
!   integer, parameter :: ncp = 3             ! # vegetation carbon stores
!   integer, parameter :: ncs = 2             ! # soil carbon stores
!   integer, parameter :: mf = 2              ! # leaves (sunlit, shaded)
!   integer, parameter :: nrb = 3             ! # radiation bands
!   integer, parameter :: msn = 3              ! max # snow layers
!   integer, parameter :: swb = 2             ! # shortwave bands 
!   INTEGER, PARAMETER :: ms = 6  ! # soil layers
!
!   integer            :: mvtype               ! total # vegetation types,   from input
!!  integer            :: mvtype               ! total # vegetation types,   from input
!   integer            :: mstype               ! total # soil types,         from input
!
!   integer :: mland                           ! # land grid cells

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   ! definition of constants used throughout model
   type physical_constants
      real :: capp   = 1004.64  ! air spec. heat (J/kg/K)
      real :: hl = 2.5104e6  ! air spec. heat (J/kg/K)
      real :: hlf = 0.335e6   ! latent heat of fusion
      real :: dheat  = 21.5E-6  ! molecular diffusivity for heat
      real :: grav   = 9.80     ! gravity acceleration (m/s2)
      real :: rgas   = 8.3143   ! universal gas const  (J/mol/K)
      real :: rmair  = 0.02897  ! molecular wt: dry air (kg/mol)
      real :: rmh2o  = 0.018016 ! molecular wt: water	(kg/mol)
      real :: sboltz = 5.67e-8  ! Stefan-Boltz. constant (W/m2/K4)
      real :: tfrz   = 273.16   ! Temp (K) corresp. to 0 C
      ! Teten coefficients
      real:: tetena = 6.106  ! ??? refs?
      real:: tetenb = 17.27
      real:: tetenc = 237.3
      ! Aerodynamic parameters, diffusivities, water density:
      real:: vonk   = 0.40 ! von Karman constant
      real:: a33    = 1.25 ! inertial sublayer sw/us
      real:: csw    = 0.50 ! canopy sw decay (Weil theory)
      real:: ctl    = 0.40 ! Wagga wheat (RDD 1992, Challenges)
      real:: apol   = 0.70 ! Polhausen coeff: single-sided plate
      real:: prandt = 0.71 ! Prandtl number: visc/diffh
      real:: schmid = 0.60 ! Schmidt number: visc/diffw
      real:: diffwc = 1.60 ! diffw/diffc = H2O/CO2 diffusivity
      real:: rhow   = 1000.0 ! liquid water density   [kg/m3]
      real:: emleaf = 1.0  ! leaf emissivity
      real:: emsoil = 1.0  ! soil emissivity
      real:: crd = 0.3     ! element drag coefficient
      real:: csd = 0.003   ! substrate drag coefficient
      !jhan:hardwire for now. note beta2 = crd/csd
      real:: beta2 = 0.3/0.003 ! ratio cr/cs
      real:: ccd   = 15.0  ! constant in d/h equation
      real:: ccw_c = 2.0   ! ccw=(zw-d)/(h-d)
      real:: usuhm = 0.3   ! (max of us/uh)
      ! Turbulence parameters:
      integer:: niter = 4  ! number of iterations for za/L
      real:: zetmul = 0.4  ! if niter=2, final zeta=zetmul*zetar(2)
      real:: zeta0  = 0.0  ! initial value of za/L
      real:: zetneg = -15.0 ! negative limit on za/L when niter>=3
      real:: zetpos = 1.0  ! positive limit on za/L when niter>=3
      real:: zdlin  = 1.0  ! height frac of d below which TL linear
      real:: umin   = 0.01
   end type physical_constants

   type math_constants
      real :: pi_c = 3.1415927
      !jhan:hardwire for now. note pi180= pi_c/180
      real :: pi180 = 3.1415927/ 180.0 ! radians / degree
      real :: two_pi = 2.0 * 3.1415927
   end type math_constants

   type other_constants
      real, DIMENSION(nrb) :: gauss_w=(/0.308,0.514,0.178/) ! Gaussian integ. weights
      real, DIMENSION(nrb) :: refl = (/ 0.07, 0.425, 0.00 /) ! YP nov2009
      real, DIMENSION(nrb) :: taul = (/ 0.07, 0.425, 0.00/)  ! leaf transmittance
      !--- jhan: can make these trigger of #defines/namelist
      real:: RAD_THRESH = 0.01 
      real:: LAI_THRESH = 0.01 
   end type other_constants

   type photosynthetic_constants
      integer:: maxiter=20 ! max # interations for leaf temperature
      real :: a1c3 = 9.0
      real :: a1c4 = 4.0
      real :: alpha3 = 0.200
      real :: alpha4  = 0.05
      real :: cfrd3  = 0.015
      real :: cfrd4  = 0.025
      real :: conkc0 = 302.e-6  !mol mol^-1
      real :: conko0 = 256.e-3  !mol mol^-1
      real :: convx3 = 1.0E-2
      real :: convx4 = 0.8
      real :: d0c3 = 1500.0
      real :: d0c4 = 1500.0
      real :: ekc = 59430.0  !J mol^-1
      real :: eko = 36000.0  !J mol^-1
      real :: gam0 = 28.0E-6  !mol mol^-1 @ 20C = 36.9 @ 25C
      real :: gam1 = 0.0509
      real :: gam2 = 0.0010
      real :: gsw03  = 0.01
      real :: gsw04  = 0.04
      real :: rgbwc  = 1.32
      real :: rgswc  = 1.57
      real :: tmaxj  = 45.0
      real :: tmaxv  = 45.0
      real :: tminj  = -5.0
      real :: tminv  = -5.0
      real :: toptj  = 20.0
      real :: toptv  = 20.0
      real :: trefk= 298.2  !reference temperature K
   end type photosynthetic_constants


   type model_constants
      type (physical_constants) :: phys
      type (math_constants) :: math
      type (other_constants) :: other
      type (photosynthetic_constants) :: photo
   end type model_constants


   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   TYPE air_auto_type
    REAL, DIMENSION(:), POINTER :: air_rho  ! dry air density (kg m-3)
    REAL, DIMENSION(:), POINTER :: air_rlam ! latent heat for water (j/kg)
    REAL, DIMENSION(:), POINTER :: air_epsi ! d(qsat)/dT ((kg/kg)/K)
    REAL, DIMENSION(:), POINTER :: air_visc ! air kinematic viscosity (m2/s)
    REAL, DIMENSION(:), POINTER :: air_psyc ! psychrometric constant
    REAL, DIMENSION(:), POINTER :: air_dsatdk ! d(es)/dT (mb/K)
    REAL, DIMENSION(:), POINTER :: air_cmolar ! conv. from m/s to mol/m2/s
   end TYPE air_auto_type

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   !type met_auto_type
   !end type met_auto_type
  !jhan:these are only those being passed to define_air 
   type met_in_type
     real, dimension(:), pointer :: met_tvair   ! within canopy air temperature (oK)
     real, dimension(:), pointer :: met_pmb     ! surface air pressure (mbar)
   end type met_in_type
  
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   
!     ! "coexp": coefficient in exponential in-canopy wind profile
!     ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
!     ! canopy and roughness-sublayer U(z) at z=h
   type rough_in_type
     REAL, DIMENSION(:), POINTER :: rough_hruff_grmx ! max ht of canopy from tiles on same grid 
     REAL, DIMENSION(:), POINTER :: rough_za_uv   ! level of lowest atmospheric model layer
     REAL, DIMENSION(:), POINTER :: rough_za_tq   ! level of lowest atmospheric model layer
   end type rough_in_type

   type rough_out_type
     REAL, DIMENSION(:), POINTER   :: rough_z0m  ! roughness length
     REAL, DIMENSION(:), POINTER   :: rough_zref_tq ! Reference height for met forcing
   end type rough_out_type
   
   type rough_auto_type
     REAL, DIMENSION(:), POINTER :: rough_hruff ! canopy height above snow level
     REAL, DIMENSION(:), POINTER   :: rough_coexp ! Extinction coef for wind profile in canopy
     REAL, DIMENSION(:), POINTER   :: rough_disp  ! zero-plane displacement
     REAL, DIMENSION(:), POINTER   :: rough_rt0us ! eq. 3.54, SCAM manual (CSIRO tech report 132)
     REAL, DIMENSION(:), POINTER   :: rough_rt1usa ! resistance from disp to hruf
     REAL, DIMENSION(:), POINTER   :: rough_rt1usb ! resist fr hruf to zruffs (zref if zref<zruffs)
     REAL, DIMENSION(:), POINTER   :: rough_rt1 ! 1/aerodynamic conductance
     REAL, DIMENSION(:), POINTER   :: rough_usuh ! Friction velocity/windspeed at canopy height
     REAL, DIMENSION(:), POINTER   :: rough_zref_uv ! Reference height for met forcing
     REAL, DIMENSION(:), POINTER   :: rough_zruffs ! SCALAR Roughness sublayer depth (ground=origin)
     REAL, DIMENSION(:), POINTER   :: rough_z0soilsn ! roughness length of bare soil surface
   end type rough_auto_type
   
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   type ssnow_in_type
     REAL, DIMENSION(:), POINTER :: ssnow_snowd   ! snow depth (liquid water)
     REAL, DIMENSION(:), POINTER :: ssnow_ssdnn   ! average snow density
   end type ssnow_in_type
      
   type ssnow_out_type
   end type ssnow_out_type
      
   type ssnow_auto_type
   end type ssnow_auto_type
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   type veg_in_type
     REAL, DIMENSION(:), POINTER :: veg_hc     ! roughness height of canopy (veg - snow)
     REAL, DIMENSION(:), POINTER :: veg_vlai   ! leaf area index
     INTEGER,DIMENSION(:), POINTER :: veg_iveg ! vegetation type
   end type veg_in_type
      
   type veg_out_type
   end type veg_out_type
   
   type veg_auto_type
   end type veg_auto_type
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   
   type canopy_in_type
   end type canopy_in_type
      
   type canopy_out_type
   end type canopy_out_type
   
   type canopy_auto_type
     REAL, DIMENSION(:), POINTER :: canopy_vlaiw  ! lai adj for snow depth for calc of resistances
     REAL, DIMENSION(:), POINTER :: canopy_rghlai  ! lai adj for snow depth for calc of resistances
   end type canopy_auto_type
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   
!   type _in_type
!   end type _in_type
!      
!   type _out_type
!   end type _out_type
!   
!   type _auto_type
!   end type _auto_type
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type model_in 
      type (met_in_type) :: met 
      type (rough_in_type) :: rough 
      type (ssnow_in_type) :: ssnow
      type (veg_in_type) :: veg
   end type model_in 

   type model_out 
      type (rough_out_type) :: rough 
   end type model_out 

   type model_auto
      type (air_auto_type) :: air 
      type (rough_auto_type) :: rough 
      type (canopy_auto_type) :: canopy
   end type model_auto

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type model_type 
      real, allocatable :: lat(:), lon(:), tile_frac(:)
      integer, allocatable :: tile(:)
      type (model_in) :: i 
      type (model_out) :: o 
      type (model_auto) :: a 
   end type model_type 


   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
      
   !instantiate types
   type (model_constants),save, target :: const
   type (model_type),save :: cable 
!  cable%a%air%met_pmb
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   

contains


   !subroutine alloc_cable_type(mp) 
   !   implicit none
   !   integer, intent(in) :: mp 
   !   
   !      type(model) :: cable
   !             
   !      allocate(cable%gridcell(mp))
   !      allocate(cable%lat(mp))
   !      allocate(cable%lon(mp))

   !   return
   !end subroutine alloc_cable_type 

end module cable_data_module 


















!===================================================================================
!cable %
!      % constants
!      % forcings
!      % fluxes
!      % state
!     
!      =================================================================================
!      constants   % math
!                  % physical  
!                  % photosynthetic
!                  % other
!
!      =================================================================================
!      forcings    % met
!                  % rad
!
!      =================================================================================
!      !this is what is returned and affects the atmosphere/climate            
!      fluxes      % canopy
!                  % co2 etc
!
!      =================================================================================
!      state       % prognostics
!                  % diagnostics
!                                              
!                  =============================================================================
!                  prognostics % canopy
!                              % soil
!                              % ssoil (change this to snow%)
!                              % albedo (which is forward)
!                  
!                  =============================================================================
!                  diagnostics % canopy
!                              % soil                 
!                              % ssoil 
!                              % met 
!                              % rad
!                              % bgc etc
!===================================================================================
!
!!jhan: ?? maybe we can have this data structure. make local pointers in "interface function" to all major components ??
!


