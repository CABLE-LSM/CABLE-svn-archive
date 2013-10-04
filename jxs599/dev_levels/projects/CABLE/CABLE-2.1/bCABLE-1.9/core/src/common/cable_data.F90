module cable_data_module 
   use define_dimensions, only : i_d, r_1, r_2, nrb
   implicit none
   
   !public :: cable, const
   !public :: air_in, air_out 
   !private

 

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

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   type model_constants
      type (physical_constants) :: phys
      type (math_constants) :: math
      type (other_constants) :: other
      type (photosynthetic_constants) :: photo
   end type model_constants
   

   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   
   

   
   !climate forcing to land-surface
   type cable_in
      
      real(r_1), dimension(:), pointer :: &
         met_fld, &        ! downward long-wave radiation (W/m2)
         met_coszen, &     ! cos(zenith angle of sun)
         met_precip, &     ! rainfall (liquid+solid)(mm/dels)
         met_precip_sn, &  ! solid preipitation only (mm/dels)
         met_tk, &         ! surface air temperature (oK)
         met_qv, &         ! surface specific humidity (g/g)
         met_pmb, &        ! surface air pressure (mbar)
         met_ua, &         ! surface wind speed (m/s)
         rough_za_uv, &    ! level of lowest atmospheric model layer
         rough_za_tq, &    ! level of lowest atmospheric model layer
         fland, &          ! factor for latent heat ?? do we need this here 
         met_ca            ! CO2 concentration (mol/mol)

      real(r_1), dimension(:,:), pointer ::  &
         met_fsd, &     ! downward short-wave radiation (W/m2)
         rad_fbeam      ! beam fraction 

   end type cable_in


   
   !land-surface feedback to  climate
   type cable_out
      
      real(r_1), dimension(:), pointer :: &
         snow_rnof1, &     ! surface runoff (mm/dels)
         snow_rnof2, &     ! deep drainage (mm/dels)
         rough_z0m, &      ! roughness length
         rough_zref_tq, &  ! Reference height for met forcing
         rad_trad, &       ! radiative temperature (soil and veg)
         canopy_fe, &      ! total latent heat (W/m2)
         canopy_fh, &      ! total sensible heat (W/m2)
         canopy_us, &      ! friction velocity
         canopy_cdtq, &    ! drag coefficient for momentum
         canopy_fwet, &    ! fraction of canopy wet
         canopy_rnet, &    ! net radiation absorbed by surface (W/m2)
         canopy_epot, &    ! total potential evaporation 
         canopy_through, & ! canopy throughfall (mm)
         canopy_wetfac_cs,&! 
         canopy_frs, &            ! soil respiration (g C m-2 s-1)
         canopy_fnee, &           ! net carbon flux (g C m-2 s-1)
         canopy_frday             ! daytime leaf resp
         

   end type cable_out


   
   
   !idea is that many variables passed back and fourth to/from CABLE can be kept 
   !in CABLE mem, and not reinitialized from _IN on every timestep
   !updated by CABLE each timestep. maybe initialized through host (UM), 
   !and passed back for diag., but essentially can be stored in CABLE's memory.
   type cable_mem

      integer(i_d), dimension(:), pointer :: snow_isflag ! 0 => no snow 1 => snow
      
      real(r_1), dimension(:), pointer :: &
         soil_hyds, &      ! hydraulic conductivity @ saturation [m/s], Ksat
         !this is initialized on first timestep from bexp in UM, check if spatially explicit as well
         soil_bch, &       ! parameter b in Campbell equation
         soil_ssat, &      ! vol H2O @ saturation
         soil_sucs, &      ! suction at saturation (m)
         soil_swilt, &     ! vol H2O @ wilting
         soil_sfc, &       ! vol H2O @ field capacity
         snow_snowd, &     ! snow depth (liquid water)
         snow_osnowd, &    ! snow depth from previous time step
         snow_ssdnn, &     ! average snow density
         canopy_cansto     ! canopy water storage (mm)
       
      real(r_2), dimension(:), pointer :: &
         soil_cnsd ! thermal conductivity of dry soil [W/m/K]
          
      real(r_2), dimension(:,:), pointer :: &
         snow_wb, &    ! volumetric soil moisture (solid+liq)
         snow_wbice    ! soil ice

      real(r_1), dimension(:,:), pointer :: & 
         snow_sdepth, &    ! snow depth
         snow_sconds, &    ! snow_cond in UM
         snow_smass, &     ! snow mass
         snow_ssdn, &      ! snow densities
         snow_tggsn, &     ! snow temperature in K
         snow_tgg   ! soil temperature in K
             
   end type cable_mem


   
   
   !diagnostics computed in CABLE 
   type cable_auto_bal
      
      real(r_1), dimension(:), pointer :: &
         drybal, & 
         wetbal 
   
   end type cable_auto_bal


   type cable_auto_canopy
      
      real(r_1), dimension(:), pointer ::  &
         canst1
      
      real(r_1), dimension(:,:), pointer ::  &
         zetar       ! stability correction
         
   end type cable_auto_canopy


   type model_cable_auto
      type (cable_auto_bal) :: bal
      type (cable_auto_canopy) :: canopy
   end type model_cable_auto


      
   type model_cable
      real, pointer, dimension(:) :: &
         lat, & 
         lon, &
         tile_frac

      integer, pointer, dimension(:)  :: &
         tile(:)

      type (cable_in) :: i
      type (cable_out) :: o
      type (cable_mem) :: mem
   end type model_cable

   

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   
   type canopy_in
      REAL(r_1), DIMENSION(:), POINTER :: &
         canopy_canst1 !, & !prev.veg%canst1
   end type canopy_in
   
   type canopy_out
      !REAL(r_1), DIMENSION(:), POINTER :: bal%drybal 
      !REAL(r_1), DIMENSION(:), POINTER :: wetbal 
   end type canopy_out
   
   type canopy_inout
      REAL(r_1), DIMENSION(:), POINTER :: &
         canopy_cansto!, & !prev.canopy%cansto
         
   end type canopy_inout
      
   
   !type model_canopy
   !   type (canopy_in) :: i
   !   type (canopy_out) :: o
   !   type (canopy_inout) :: io
   !end type model_canopy


   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   !instantiate types
   type (model_constants) ,save, target :: const
   type (model_cable) ,save, target :: cable
   type (model_cable_auto) ,save, target :: cable_auto
   
   
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   

contains

   
   
   
   
   subroutine allocate_cable_types( met, rad, soil, ssoil, canopy, rough, bal, veg)
      use define_dimensions, only : mp, ms, msn, swb !as well as inherited from above 
      use define_types
      implicit none
      !type (model_cable) ,intent(inout) :: cable
      type (canopy_type), intent(in)    :: canopy ! vegetation variables
      type (met_type), intent(in)       :: met  ! met input variables
      type (balances_type), intent(in)  :: bal  ! energy and water balance variables
      type (radiation_type), intent(in) :: rad  ! radiation variables
      type (roughness_type), intent(in) :: rough ! roughness varibles
      type (soil_parameter_type), intent(in) :: soil ! soil parameters	
      type (soil_snow_type), intent(in) :: ssoil ! soil and snow variables                
      type (veg_parameter_type), intent(in) :: veg 

      !jhan:all of i,o, auto can be deallocated every timestep


         allocate( cable%mem%soil_hyds(mp) )
         allocate( cable%mem%soil_bch(mp) )
         allocate( cable%mem%soil_ssat(mp) )
         allocate( cable%mem%soil_sucs(mp) )
         allocate( cable%mem%soil_swilt(mp) )
         allocate( cable%mem%soil_sfc(mp) )
         allocate( cable%mem%soil_cnsd(mp) ) 

         
         allocate( cable%mem%canopy_cansto(mp) )


         allocate( cable%mem%snow_isflag(mp) )      
         allocate( cable%mem%snow_snowd(mp) )
         allocate( cable%mem%snow_osnowd(mp) )
         allocate( cable%mem%snow_ssdnn(mp) )

         allocate( cable%mem%snow_tggsn(mp, msn) )
         allocate( cable%mem%snow_sdepth(mp, msn) )
         allocate( cable%mem%snow_smass(mp, msn) )
         allocate( cable%mem%snow_ssdn(mp, msn) )
         allocate( cable%mem%snow_sconds(mp, msn) )

         allocate( cable%mem%snow_wb(mp,ms) )
         allocate( cable%mem%snow_wbice(mp,ms) )
         allocate( cable%mem%snow_tgg(mp,ms) )   

         
         
         allocate( cable%i%met_fld(mp) )
         allocate( cable%i%met_coszen(mp) )
         allocate( cable%i%met_precip(mp) )
         allocate( cable%i%met_precip_sn(mp) )
         allocate( cable%i%met_tk(mp) )
         allocate( cable%i%met_qv(mp) )
         allocate( cable%i%met_pmb(mp) )
         allocate( cable%i%met_ua(mp) )
         allocate( cable%i%rough_za_uv(mp) )
         allocate( cable%i%rough_za_tq(mp) )
         allocate( cable%i%fland(mp) )
         allocate( cable%i%met_ca(mp) ) 

     
         allocate( cable%i%met_fsd(mp,swb) )
         allocate( cable%i%rad_fbeam(mp,nrb) )
         

         
         !out   
         allocate( cable%o%rough_z0m(mp) ) 
         allocate( cable%o%rough_zref_tq(mp) )
         allocate( cable%o%rad_trad(mp) )
         allocate( cable%o%canopy_fe(mp) )
         allocate( cable%o%canopy_fh(mp) )
         allocate( cable%o%canopy_us(mp) )
         allocate( cable%o%canopy_cdtq(mp) )
         allocate( cable%o%canopy_fwet(mp) )
         allocate( cable%o%canopy_rnet(mp) )
         allocate( cable%o%canopy_epot(mp) )
         allocate( cable%o%canopy_through(mp) )
         
         allocate( cable%o%canopy_wetfac_cs(mp) ) 
         allocate( cable%o%snow_rnof1(mp) )
         allocate( cable%o%snow_rnof2(mp) )

         allocate( cable%o%canopy_frs(mp) ) 
         allocate( cable%o%canopy_fnee(mp) )
         allocate( cable%o%canopy_frday(mp) )
  
         
         !as a starting point all of these vars are pointing to the 
         !counter parts in the old data structure. gradually replace 
         !instances in the bottom level subrs and then rm old structure

         cable%mem%soil_hyds     => soil%hyds
         cable%mem%soil_bch      => soil%bch
         cable%mem%soil_ssat     => soil%ssat
         cable%mem%soil_sucs     => soil%sucs
         cable%mem%soil_swilt    => soil%swilt
         cable%mem%soil_sfc      => soil%sfc
         cable%mem%soil_cnsd     => soil%cnsd

         cable%mem%snow_isflag   => ssoil%isflag
         cable%mem%snow_tgg      => ssoil%tgg   
         cable%mem%snow_snowd    => ssoil%snowd
         cable%mem%snow_osnowd   => ssoil%osnowd
         cable%mem%snow_ssdnn    => ssoil%ssdnn
         cable%mem%snow_sdepth   => ssoil%sdepth
         cable%mem%snow_smass    => ssoil%smass
         cable%mem%snow_ssdn     => ssoil%ssdn
         cable%mem%snow_tggsn    => ssoil%tggsn
         cable%mem%snow_sconds   => ssoil%sconds
         cable%mem%snow_wb       => ssoil%wb
         cable%mem%snow_wbice    => ssoil%wbice
         
         cable%mem%canopy_cansto => canopy%cansto
         
      
         cable%i%met_fld         => met%fld
         cable%i%met_coszen      => met%coszen
         cable%i%met_precip      => met%precip
         cable%i%met_precip_sn   => met%precip_sn
         cable%i%met_tk          => met%tk
         cable%i%met_qv          => met%qv
         cable%i%met_pmb         => met%pmb
         cable%i%met_ua          => met%ua
         cable%i%met_ca          => met%ca
         cable%i%met_fsd         => met%fsd
         
         cable%i%rad_fbeam       => rad%fbeam
         
         cable%i%rough_za_uv     => rough%za_uv
         cable%i%rough_za_tq     => rough%za_tq
        
         cable%i%fland           => ssoil%fland
         
         
         cable%o%rough_z0m       => rough%z0m
         cable%o%rough_zref_tq   => rough%zref_tq
         cable%o%rad_trad        => rad%trad
         cable%o%canopy_fe       => canopy%fe    
         cable%o%canopy_fh       => canopy%fh     
         cable%o%canopy_us       => canopy%us   
         cable%o%canopy_cdtq     => canopy%cdtq
         cable%o%canopy_fwet     => canopy%fwet 
         cable%o%canopy_rnet     => canopy%rnet
         cable%o%canopy_epot     => canopy%epot
         cable%o%canopy_through  => canopy%through
         cable%o%canopy_wetfac_cs => canopy%wetfac_cs
         
         cable%o%snow_rnof1      => ssoil%rnof1
         cable%o%snow_rnof2      => ssoil%rnof2

         cable%o%canopy_frs    => canopy%frs           
         cable%o%canopy_fnee   => canopy%fnee
         cable%o%canopy_frday  => canopy%frday
   
    
         
      return
   end subroutine allocate_cable_types 


   
  subroutine allocate_cable_auto( met, rad, soil, ssoil, canopy, rough, bal, veg)
      use define_dimensions, only : mp, ms, msn, swb !as well as inherited from above 
      use define_types
      implicit none
      !type (model_cable_auto), intent(inout) :: cable_auto
      type (canopy_type), intent(in)    :: canopy ! vegetation variables
      type (met_type), intent(in)       :: met  ! met input variables
      type (balances_type), intent(in)  :: bal  ! energy and water balance variables
      type (radiation_type), intent(in) :: rad  ! radiation variables
      type (roughness_type), intent(in) :: rough ! roughness varibles
      type (soil_parameter_type), intent(in) :: soil ! soil parameters	
      type (soil_snow_type), intent(in) :: ssoil ! soil and snow variables                
      type (veg_parameter_type), intent(in) :: veg 

   
         
         allocate( cable_auto%bal%drybal(mp) ) 
         allocate( cable_auto%bal%wetbal(mp) ) 


         allocate( cable_auto%canopy%canst1(mp) ) 

         allocate( cable_auto%canopy%zetar(mp, const%phys%niter) ) 
   
         
         cable_auto%canopy%canst1 => veg%canst1

         cable_auto%canopy%zetar => canopy%zetar
         
         cable_auto%bal%drybal   => bal%drybal 
         cable_auto%bal%wetbal   => bal%wetbal 

      return
   end subroutine allocate_cable_auto 

   

   subroutine deallocate_cable_auto() 
      use define_dimensions, only : mp!, ms, msn, swb, nrb !as well as inherited from above
      use define_types
      implicit none
      !type (model_cable_auto), intent(inout) :: cable_auto

         nullify( cable_auto%canopy%canst1 )
         nullify( cable_auto%bal%drybal ) 
         nullify( cable_auto%bal%wetbal ) 
         nullify( cable_auto%canopy%zetar ) 
      
  
      return
   end subroutine deallocate_cable_auto 

   
   

!   subroutine deallocate_cable_types() 
!      use define_dimensions, only : mp!, ms, msn, swb, nrb !as well as inherited from above
!      use define_types
!      implicit none
!      logical, save :: first_call = .true.
!      
!      if (first_call) &
!         deallocate( cable%auto%canopy%canst1 )
!
!         
!         deallocate( cable%i%met_fld )
!         deallocate( cable%i%met_coszen )
!         deallocate( cable%i%met_precip )
!         deallocate( cable%i%met_precip_sn )
!         deallocate( cable%i%met_tk )
!         deallocate( cable%i%met_qv )
!         deallocate( cable%i%met_pmb )
!         deallocate( cable%i%met_ua )
!         deallocate( cable%i%rough_za_uv )
!         deallocate( cable%i%rough_za_tq )
!         deallocate( cable%i%fland )
!         deallocate( cable%i%met_ca ) 
!
!     
!         deallocate( cable%i%met_fsd )
!         deallocate( cable%i%rad_fbeam )
!     
!
!         !out   
!         deallocate( cable%o%rough_z0m ) 
!         deallocate( cable%o%rough_zref_tq )
!         deallocate( cable%o%rad_trad )
!         deallocate( cable%o%canopy_fe )
!         deallocate( cable%o%canopy_fh )
!         deallocate( cable%o%canopy_us )
!         deallocate( cable%o%canopy_cdtq )
!         deallocate( cable%o%canopy_fwet )
!         deallocate( cable%o%canopy_rnet )
!         deallocate( cable%o%canopy_epot )
!         deallocate( cable%o%canopy_through )
!      
!         deallocate( cable%o%canopy_wetfac_cs ) 
!         deallocate( cable%o%snow_rnof1 )
!         deallocate( cable%o%snow_rnof2 )
!
!         deallocate( cable%o%canopy_frs ) 
!         deallocate( cable%o%canopy_fnee )
!         deallocate( cable%o%canopy_frday )
!  
!      
!         deallocate( cable%auto%bal%drybal ) 
!         deallocate( cable%auto%bal%wetbal ) 
!
!         deallocate( cable%auto%canopy%zetar ) 
!        
!      first_call = .false.
!      
!  
!      return
!   end subroutine deallocate_cable_types 



   subroutine allocate_canopy_types(i,o,io) 
      use define_dimensions, only : mp!, ms, msn, swb, nrb !as well as inherited from above
      use define_types
      implicit none
      !type (model_cable) ,intent(inout) :: cable
      type (canopy_in), intent(inout) :: i
      type (canopy_out), intent(inout) :: o
      type (canopy_inout), intent(inout) :: io
      logical, save :: first_call = .true.
      
      
         allocate( io%canopy_cansto(mp) )
         io%canopy_cansto => cable%mem%canopy_cansto               
         
         allocate( i%canopy_canst1(mp) )
         i%canopy_canst1 => cable_auto%canopy%canst1               

      return
   end subroutine allocate_canopy_types 



   subroutine deallocate_canopy_types( i, o, io) 
      use define_dimensions, only : mp!, ms, msn, swb, nrb !as well as inherited from above
      use define_types
      implicit none
      type (canopy_in), intent(inout) :: i
      type (canopy_out), intent(inout) :: o
      type (canopy_inout), intent(inout) :: io
         
         nullify(i%canopy_canst1)
         nullify(io%canopy_cansto)
             
      return
   end subroutine deallocate_canopy_types 





end module cable_data_module 

















