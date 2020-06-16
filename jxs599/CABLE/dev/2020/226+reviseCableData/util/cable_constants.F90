!H!s should be removed 
! cable_constants.f90
!
! Source file containing constants for CABLE
!
! Development by Ying-Ping Wang, Eva Kowalczyk, Ray Leuning
! Gab Abramowitz, Martin Dix, Harvey Davies, Mike Raupach, Matthias Cuntz
!
! bugs to bernard.pak@csiro.au
!
! This file contains modules:
!  math_constants
!  physical_constants
!  other_constants
!  photosynthetic_constants
!  spatial_heterogeneity

!H!MODULE math_constants
!H!
!H!  USE cable_def_types_mod, ONLY : i_d, r_2
!H!
!H!  IMPLICIT NONE
!H!
!H!  REAL, PARAMETER :: pi     = 3.141592653589793238462643383279502884197
!H!  REAL, PARAMETER :: pi180  = pi / 180.0 ! radians / degree
!H!  REAL, PARAMETER :: two_pi = 2.0 * pi
!H!  REAL(r_2), PARAMETER :: pi_r_2 = 3.141592653589793238462643383279502884197
!H!
!H!END MODULE math_constants
!H!
!H!!=========================================================================
!H!
!H!MODULE physical_constants
!H!
!H!  USE cable_def_types_mod, ONLY : i_d
!H!
!H!  IMPLICIT NONE
!H!
!H!  !mrd561
!H!  !constants from/for soilsnow and gw_hydro
!H!  REAL,    PARAMETER :: hl = 2.5014e6  ! air spec. heat (J/kg/K)
!H!  REAL,    PARAMETER :: hlf = 0.334e6  ! latent heat of fusion
!H!  REAL,    PARAMETER :: hls = 2.8350e6  ! latent heat of sublimation (J/kg)
!H!
!H!  REAL,    PARAMETER :: capp   = 1004.64  ! air spec. heat capacity (J/kg/K)
!H!  REAL,    PARAMETER :: dheat  = 21.5E-6  ! molecular diffusivity for heat
!H!  REAL,    PARAMETER :: grav   = 9.80     ! gravity acceleration (m/s2)
!H!  REAL,    PARAMETER :: rgas   = 8.3143   ! universal gas const  (J/mol/K)
!H!  REAL,    PARAMETER :: rmair  = 0.02897  ! molecular wt: dry air (kg/mol)
!H!  REAL,    PARAMETER :: rmh2o  = 0.018016 ! molecular wt: water       (kg/mol)
!H!  REAL,    PARAMETER :: sboltz = 5.67e-8  ! Stefan-Boltz. constant (W/m2/K4)
!H!  REAL,    PARAMETER :: tfrz   = 273.16   ! Temp (K) corresp. to 0 C
!H!  ! Teten coefficients
!H!  REAL,    PARAMETER :: tetena = 6.106    ! ??? refs?
!H!  REAL,    PARAMETER :: tetenb = 17.27
!H!  REAL,    PARAMETER :: tetenc = 237.3
!H!  !mrd561 the parameters for sat above ice
!H!  REAL,    PARAMETER :: tetena_ice = 6.1078  ! ??? refs?
!H!  REAL,    PARAMETER :: tetenb_ice = 21.875 
!H!  REAL,    PARAMETER :: tetenc_ice = 265.5 
!H!  ! Aerodynamic parameters, diffusivities, water density:
!H!  REAL,    PARAMETER :: vonk   = 0.40     ! von Karman constant
!H!  REAL,    PARAMETER :: a33    = 1.25     ! inertial sublayer sw/us
!H!  REAL,    PARAMETER :: csw    = 0.50     ! canopy sw decay (Weil theory)
!H!  REAL,    PARAMETER :: ctl    = 0.40     ! Wagga wheat (RDD 1992, Challenges)
!H!  REAL,    PARAMETER :: apol   = 0.70     ! Polhausen coeff: single-sided plate
!H!  REAL,    PARAMETER :: prandt = 0.71     ! Prandtl number: visc/diffh
!H!  REAL,    PARAMETER :: schmid = 0.60     ! Schmidt number: visc/diffw
!H!  REAL,    PARAMETER :: diffwc = 1.60     ! diffw/diffc = H2O/CO2 diffusivity
!H!  REAL,    PARAMETER :: rhow   = 1000.0   ! liquid water density   [kg/m3]
!H!  REAL,    PARAMETER :: emleaf = 1.0      ! leaf emissivity
!H!  REAL,    PARAMETER :: emsoil = 1.0      ! soil emissivity
!H!  REAL,    PARAMETER :: cr     = 0.3      ! element drag coefficient
!H!  REAL,    PARAMETER :: cs     = 0.003    ! substrate drag coefficient
!H!  REAL,    PARAMETER :: beta   = cr/cs    ! ratio cr/cs
!H!  REAL,    PARAMETER :: ccd    = 15.0     ! constant in d/h equation
!H!  REAL,    PARAMETER :: ccw    = 2.0      ! ccw=(zw-d)/(h-d)
!H!  REAL,    PARAMETER :: usuhm  = 0.3      ! (max of us/uh)
!H!  ! Turbulence  parameters:
!H!  INTEGER(i_d), PARAMETER :: niter  = 10       ! number of iterations for za/L
!H!  REAL,    PARAMETER :: zetmul = 0.4      ! if niter=2, final zeta=zetmul*zetar(2)
!H!  REAL,    PARAMETER :: zeta0  = 0.0      ! initial value of za/L
!H!  REAL,    PARAMETER :: zetneg = -10.0    ! negative limit on za/L when niter>=3
!H!  REAL,    PARAMETER :: zetpos = 0.5      ! positive limit on za/L when niter>=3
!H!  REAL,    PARAMETER :: zdlin  = 1.0      ! height frac of d below which TL linear
!H!  REAL,    PARAMETER :: umin   = 0.1
!H!
!H!END MODULE physical_constants
!H!
!H!!=========================================================================
!H!
!H!MODULE other_constants
!H!
!H!  USE cable_def_types_mod, ONLY : i_d, nrb,r_2
!H!
!H!  IMPLICIT NONE
!H!
!H!  REAL,    PARAMETER, DIMENSION(nrb) :: gauss_w =(/0.308,0.514,0.178/)   ! Gaussian integ. weights
!H!  ! values in refl and taul are slightly modified since Oct07 and Mar08 (YP)
!H!  ! leaf reflectance
!H!  REAL,    PARAMETER, DIMENSION(nrb) :: refl    = (/ 0.1, 0.425, 0.02 /) ! mar08
!H!  ! leaf transmittance
!H!  REAL,    PARAMETER, DIMENSION(nrb) :: taul    = (/ 0.1, 0.425, 0.02 /) ! mar08
!H!  INTEGER(i_d), PARAMETER                 :: istemp  = 4                      ! soil temp:     1,2,3,4 = FR,kf,mrr,mrrkf
!H!  INTEGER(i_d), PARAMETER                 :: ismois  = 2                      ! soil moist:  1,2,3     = MP84,NP89,Richards
!H!  INTEGER(i_d), PARAMETER                 :: isinf   = 2                      ! soil infilt: 1,2       = MP84, FC96
!H!  INTEGER(i_d), PARAMETER                 :: isevap  = 2                      ! soil evap: 1,2,3 = alfa,beta,threshold
!H!  INTEGER(i_d), PARAMETER                 :: itherm  = 1                      ! VW or KGK algorithm for hconds,rkapps
!H!  INTEGER(i_d), PARAMETER                 :: irktem  = 5                      ! RK steps in soil temp schemes
!H!  INTEGER(i_d), PARAMETER                 :: irkmoi  = 5                      ! RK steps in soil moisture schemes
!H!  ! soil water  parameters:
!H!  REAL,    PARAMETER                 :: etarct  = 0.7                    ! rel soil moisture for finding zst1,zst2
!H!  REAL,    PARAMETER                 :: dbde    = 1.3333                 ! d(beta)/d(etar): 4/3, 8 for D78,KGK91
!H!  INTEGER(i_d), PARAMETER                 :: istsw   = 1                      !
!H!  INTEGER(i_d), PARAMETER                 :: iresp   = 0                      ! unscaled (iresp=0) or scaled (iresp=1) respiration
!H!
!H!END MODULE other_constants
!H!
!H!
!H!!=========================================================================
!H!
!H!MODULE photosynthetic_constants
!H!
!H!  USE cable_def_types_mod, ONLY : i_d
!H!
!H!  IMPLICIT NONE
!H!
!H!  INTEGER(i_d), PARAMETER :: maxiter         = 20       ! max # interations for leaf temperature
!H!  ! a1c3 is defined inside cable_canopy.f90
!H!  REAL,    PARAMETER :: a1c4_default    = 4.0
!H!  REAL,    PARAMETER :: a1c3_default    = 9.0
!H!  REAL,    PARAMETER :: d0c3_hawkesbury = 5000.0   ! Empirical coef for vpd sensitivity of stomata
!H!  REAL,    PARAMETER :: d0c3_default    = 1500.0   ! Empirical coef for vpd sensitivity of stomata
!H!  REAL,    PARAMETER :: d0c4_default    = 1500.0   ! Empirical coef for vpd sensitivity of stomata
!H!  REAL,    PARAMETER :: alpha3          = 0.2
!H!  REAL,    PARAMETER :: alpha4          = 0.05
!H!  REAL,    PARAMETER :: cfrd3           = 0.015
!H!  REAL,    PARAMETER :: cfrd4           = 0.025
!H!  REAL,    PARAMETER :: conkc0          = 302.0E-6 !mol mol^-1
!H!  REAL,    PARAMETER :: conko0          = 256.0E-3 !mol mol^-1
!H!  REAL,    PARAMETER :: convx3          = 0.01
!H!  REAL,    PARAMETER :: convx4          = 0.8
!H!  REAL,    PARAMETER :: ekc             = 59430.0  !J mol^-1
!H!  REAL,    PARAMETER :: eko             = 36000.0  !J mol^-1
!H!  REAL,    PARAMETER :: gam0            = 28.0E-6  !mol mol^-1 @ 20C = 36.9 @ 25C
!H!  REAL,    PARAMETER :: gam1            = 0.0509
!H!  REAL,    PARAMETER :: gam2            = 0.0010
!H!  REAL,    PARAMETER :: gsw03           = 0.01
!H!  REAL,    PARAMETER :: gsw04           = 0.04
!H!  REAL,    PARAMETER :: rgbwc           = 1.32
!H!  REAL,    PARAMETER :: rgswc           = 1.57
!H!  REAL,    PARAMETER :: tmaxj           = 45.0
!H!  REAL,    PARAMETER :: tmaxv           = 45.0
!H!  REAL,    PARAMETER :: tminj           = -5.0
!H!  REAL,    PARAMETER :: tminv           = -5.0
!H!  REAL,    PARAMETER :: toptj           = 20.0
!H!  REAL,    PARAMETER :: toptv           = 20.0
!H!  REAL,    PARAMETER :: trefk           = 298.2    !reference temperature K
!H!
!H!END MODULE photosynthetic_constants
