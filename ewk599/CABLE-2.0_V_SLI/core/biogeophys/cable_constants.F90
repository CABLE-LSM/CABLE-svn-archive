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

MODULE math_constants

  USE cable_def_types_mod, ONLY : i_d, r_1, r_2

  IMPLICIT NONE

  REAL(r_1), PARAMETER :: pi     = 3.141592653589793238462643383279502884197
  REAL(r_1), PARAMETER :: pi180  = pi / 180.0 ! radians / degree
  REAL(r_1), PARAMETER :: two_pi = 2.0 * pi
  REAL(r_2), PARAMETER :: pi_r_2 = 3.141592653589793238462643383279502884197

END MODULE math_constants

!=========================================================================

MODULE physical_constants

  USE cable_def_types_mod, ONLY : i_d, r_1

  IMPLICIT NONE

  REAL(r_1),    PARAMETER :: capp   = 1004.64  ! air spec. heat capacity (J/kg/K)
  REAL(r_1),    PARAMETER :: dheat  = 21.5E-6  ! molecular diffusivity for heat
  REAL(r_1),    PARAMETER :: grav   = 9.80     ! gravity acceleration (m/s2)
  REAL(r_1),    PARAMETER :: rgas   = 8.3143   ! universal gas const  (J/mol/K)
  REAL(r_1),    PARAMETER :: rmair  = 0.02897  ! molecular wt: dry air (kg/mol)
  REAL(r_1),    PARAMETER :: rmh2o  = 0.018016 ! molecular wt: water       (kg/mol)
  REAL(r_1),    PARAMETER :: sboltz = 5.67e-8  ! Stefan-Boltz. constant (W/m2/K4)
  REAL(r_1),    PARAMETER :: tfrz   = 273.16   ! Temp (K) corresp. to 0 C
  ! Teten coefficients
  REAL(r_1),    PARAMETER :: tetena = 6.106    ! ??? refs?
  REAL(r_1),    PARAMETER :: tetenb = 17.27
  REAL(r_1),    PARAMETER :: tetenc = 237.3
  ! Aerodynamic parameters, diffusivities, water density:
  REAL(r_1),    PARAMETER :: vonk   = 0.40     ! von Karman constant
  REAL(r_1),    PARAMETER :: a33    = 1.25     ! inertial sublayer sw/us
  REAL(r_1),    PARAMETER :: csw    = 0.50     ! canopy sw decay (Weil theory)
  REAL(r_1),    PARAMETER :: ctl    = 0.40     ! Wagga wheat (RDD 1992, Challenges)
  REAL(r_1),    PARAMETER :: apol   = 0.70     ! Polhausen coeff: single-sided plate
  REAL(r_1),    PARAMETER :: prandt = 0.71     ! Prandtl number: visc/diffh
  REAL(r_1),    PARAMETER :: schmid = 0.60     ! Schmidt number: visc/diffw
  REAL(r_1),    PARAMETER :: diffwc = 1.60     ! diffw/diffc = H2O/CO2 diffusivity
  REAL(r_1),    PARAMETER :: rhow   = 1000.0   ! liquid water density   [kg/m3]
  REAL(r_1),    PARAMETER :: emleaf = 1.0      ! leaf emissivity
  REAL(r_1),    PARAMETER :: emsoil = 1.0      ! soil emissivity
  REAL(r_1),    PARAMETER :: cr     = 0.3      ! element drag coefficient
  REAL(r_1),    PARAMETER :: cs     = 0.003    ! substrate drag coefficient
  REAL(r_1),    PARAMETER :: beta   = cr/cs    ! ratio cr/cs
  REAL(r_1),    PARAMETER :: ccd    = 15.0     ! constant in d/h equation
  REAL(r_1),    PARAMETER :: ccw    = 2.0      ! ccw=(zw-d)/(h-d)
  REAL(r_1),    PARAMETER :: usuhm  = 0.3      ! (max of us/uh)
  ! Turbulence  parameters:
  INTEGER(i_d), PARAMETER :: niter  = 10       ! number of iterations for za/L
  REAL(r_1),    PARAMETER :: zetmul = 0.4      ! if niter=2, final zeta=zetmul*zetar(2)
  REAL(r_1),    PARAMETER :: zeta0  = 0.0      ! initial value of za/L
  REAL(r_1),    PARAMETER :: zetneg = -10.0    ! negative limit on za/L when niter>=3
  REAL(r_1),    PARAMETER :: zetpos = 0.5      ! positive limit on za/L when niter>=3
  REAL(r_1),    PARAMETER :: zdlin  = 1.0      ! height frac of d below which TL linear
  REAL(r_1),    PARAMETER :: umin   = 1.0

END MODULE physical_constants

!=========================================================================

MODULE other_constants

  USE cable_def_types_mod, ONLY : i_d, r_1, nrb

  IMPLICIT NONE

  REAL(r_1),    PARAMETER, DIMENSION(nrb) :: gauss_w =(/0.308,0.514,0.178/)   ! Gaussian integ. weights
  ! values in refl and taul are slightly modified since Oct07 and Mar08 (YP)
  ! leaf reflectance
  REAL(r_1),    PARAMETER, DIMENSION(nrb) :: refl    = (/ 0.1, 0.425, 0.02 /) ! mar08
  ! leaf transmittance
  REAL(r_1),    PARAMETER, DIMENSION(nrb) :: taul    = (/ 0.1, 0.425, 0.02 /) ! mar08
  INTEGER(i_d), PARAMETER                 :: istemp  = 4                      ! soil temp:     1,2,3,4 = FR,kf,mrr,mrrkf
  INTEGER(i_d), PARAMETER                 :: ismois  = 2                      ! soil moist:  1,2,3     = MP84,NP89,Richards
  INTEGER(i_d), PARAMETER                 :: isinf   = 2                      ! soil infilt: 1,2       = MP84, FC96
  INTEGER(i_d), PARAMETER                 :: isevap  = 2                      ! soil evap: 1,2,3 = alfa,beta,threshold
  INTEGER(i_d), PARAMETER                 :: itherm  = 1                      ! VW or KGK algorithm for hconds,rkapps
  INTEGER(i_d), PARAMETER                 :: irktem  = 5                      ! RK steps in soil temp schemes
  INTEGER(i_d), PARAMETER                 :: irkmoi  = 5                      ! RK steps in soil moisture schemes
  ! soil water  parameters:
  REAL(r_1),    PARAMETER                 :: etarct  = 0.7                    ! rel soil moisture for finding zst1,zst2
  REAL(r_1),    PARAMETER                 :: dbde    = 1.3333                 ! d(beta)/d(etar): 4/3, 8 for D78,KGK91
  INTEGER(i_d), PARAMETER                 :: istsw   = 1                      !
  INTEGER(i_d), PARAMETER                 :: iresp   = 0                      ! unscaled (iresp=0) or scaled (iresp=1) respiration

END MODULE other_constants

!=========================================================================

MODULE photosynthetic_constants

  USE cable_def_types_mod, ONLY : i_d, r_1

  IMPLICIT NONE

  INTEGER(i_d), PARAMETER :: maxiter         = 20       ! max # interations for leaf temperature
  ! a1c3 is defined inside cable_canopy.f90
  REAL(r_1),    PARAMETER :: a1c4_default    = 4.0
  REAL(r_1),    PARAMETER :: a1c3_default    = 9.0
  REAL(r_1),    PARAMETER :: d0c3_hawkesbury = 5000.0   ! Empirical coef for vpd sensitivity of stomata
  REAL(r_1),    PARAMETER :: d0c3_default    = 1500.0   ! Empirical coef for vpd sensitivity of stomata
  REAL(r_1),    PARAMETER :: d0c4_default    = 1500.0   ! Empirical coef for vpd sensitivity of stomata
  REAL(r_1),    PARAMETER :: alpha3          = 0.2
  REAL(r_1),    PARAMETER :: alpha4          = 0.05
  REAL(r_1),    PARAMETER :: cfrd3           = 0.015
  REAL(r_1),    PARAMETER :: cfrd4           = 0.025
  REAL(r_1),    PARAMETER :: conkc0          = 302.0E-6 !mol mol^-1
  REAL(r_1),    PARAMETER :: conko0          = 256.0E-3 !mol mol^-1
  REAL(r_1),    PARAMETER :: convx3          = 0.01
  REAL(r_1),    PARAMETER :: convx4          = 0.8
  REAL(r_1),    PARAMETER :: ekc             = 59430.0  !J mol^-1
  REAL(r_1),    PARAMETER :: eko             = 36000.0  !J mol^-1
  REAL(r_1),    PARAMETER :: gam0            = 28.0E-6  !mol mol^-1 @ 20C = 36.9 @ 25C
  REAL(r_1),    PARAMETER :: gam1            = 0.0509
  REAL(r_1),    PARAMETER :: gam2            = 0.0010
  REAL(r_1),    PARAMETER :: gsw03           = 0.01
  REAL(r_1),    PARAMETER :: gsw04           = 0.04
  REAL(r_1),    PARAMETER :: rgbwc           = 1.32
  REAL(r_1),    PARAMETER :: rgswc           = 1.57
  REAL(r_1),    PARAMETER :: tmaxj           = 45.0
  REAL(r_1),    PARAMETER :: tmaxv           = 45.0
  REAL(r_1),    PARAMETER :: tminj           = -5.0
  REAL(r_1),    PARAMETER :: tminv           = -5.0
  REAL(r_1),    PARAMETER :: toptj           = 20.0
  REAL(r_1),    PARAMETER :: toptv           = 20.0
  REAL(r_1),    PARAMETER :: trefk           = 298.2    !reference temperature K

END MODULE photosynthetic_constants
