module cable_data_module 
   use define_dimensions, only : nrb 
   implicit none
   
   public :: const
   public :: air_in, air_out, air_const 
   private

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
   
   !type forcing 
   !end type forcing
   
   type model_fluxes
      real, pointer :: le(:)
   end type model_fluxes
   
   type state_prognostics
      real, pointer :: tsoil(:)
   end type state_prognostics
   
   type state_diagnostics
      real, pointer :: tscreen(:)
   end type state_diagnostics
   
   type model_state
      type (state_prognostics) :: prog
      type (state_diagnostics) :: diag
   end type model_state

   type canopy_type 
      real, pointer :: le(:)
   end type canopy_type 

   type component_vars 
      type (canopy_type) :: canopy 
   end type component_vars 

   type model 
      !idea here to keep track of gridcell per mp and have corresponding lat-lon map
      integer, pointer :: gridcell(:) 
      integer, pointer :: lat, lon
      !type (forcing) :: frc
      type (model_fluxes) :: flux
      type (model_state) :: state
      type (component_vars) :: vars 
   end type model 

   !subr specicific types
   type air_in
      real,dimension(:), pointer     :: &
         met_tvair, &
         met_pmb
   end type air_in

   type air_out
      real, dimension(:), pointer     :: &
         air_cmolar, &
         air_rho, &
         air_rlam, &
         air_dsatdk, &
         air_epsi, & 
         air_visc, & 
         air_psyc 
   end type air_out

   type air_const
      real, pointer     :: &
         CAPP, &      
         HL, &          
         RGAS, &      
         RMAIR, &     
         RMH2O, &     
         TFRZ, &      
         TETENA, &    
         TETENB, &    
         TETENC       
   end type air_const

   !instantiate type
   type (model_constants),save, target :: const
   

contains


   subroutine alloc_cable_type 
      implicit none
      integer :: n_gr =10, npatch = 30
      
         type(model) :: cable
                
         !allocate(cable%gridcell(n_gr))
         !allocate(cable%flux%le(npatch))
         !allocate(cable%state%tsoil(npatch))
         !allocate(cable%vars%canopy%le(npatch))

         !cable%vars%canopy%le => cable%flux%le
         !cable%vars%canopy%le = 12.1 
         !print *, 'cable%vars%canopy%le: ', cable%vars%canopy%le(1)
         !print *, 'cable%flux%le: ', cable%flux%le(1)
      return
   end subroutine alloc_cable_type 

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


