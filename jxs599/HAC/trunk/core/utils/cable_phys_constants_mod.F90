MODULE cable_phys_constants_mod

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
! Description:
!   CABLE physical constants
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

  REAL,    PARAMETER :: tfrz   = 273.16   ! Temp (K) corresp. to 0 C

! Will be added to as required.
  !constants from/for soilsnow and gw_hydro
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
!H!  REAL,    PARAMETER :: umin   = 1.0


END MODULE cable_phys_constants_mod
