!jhan:rename
module physical_constants
   use define_dimensions, ONLY : i_d, r_1
   implicit none
   private i_d, r_1
   
   real, PARAMETER :: capp   = 1004.64  ! air spec. heat (J/kg/K)
   real, PARAMETER :: hl = 2.5104e6  ! air spec. heat (J/kg/K)
  REAL, PARAMETER :: hlf = 0.335e6   ! latent heat of fusion
   real, PARAMETER :: dheat  = 21.5E-6  ! molecular diffusivity for heat
   real, PARAMETER :: grav   = 9.80     ! gravity acceleration (m/s2)
   real, PARAMETER :: rgas   = 8.3143   ! universal gas const  (J/mol/K)
   real, PARAMETER :: rmair  = 0.02897  ! molecular wt: dry air (kg/mol)
   real, PARAMETER :: rmh2o  = 0.018016 ! molecular wt: water	(kg/mol)
   real, PARAMETER :: sboltz = 5.67e-8  ! Stefan-Boltz. constant (W/m2/K4)
   real, PARAMETER :: tfrz   = 273.16   ! Temp (K) corresp. to 0 C
   ! Teten coefficients
   real, PARAMETER :: tetena = 6.106  ! ??? refs?
   real, PARAMETER :: tetenb = 17.27
   real, PARAMETER :: tetenc = 237.3
   ! Aerodynamic parameters, diffusivities, water density:
   real, PARAMETER :: vonk   = 0.40 ! von Karman constant
   real, PARAMETER :: a33    = 1.25 ! inertial sublayer sw/us
   real, PARAMETER :: csw    = 0.50 ! canopy sw decay (Weil theory)
   real, PARAMETER :: ctl    = 0.40 ! Wagga wheat (RDD 1992, Challenges)
   real, PARAMETER :: apol   = 0.70 ! Polhausen coeff: single-sided plate
   real, PARAMETER :: prandt = 0.71 ! Prandtl number: visc/diffh
   real, PARAMETER :: schmid = 0.60 ! Schmidt number: visc/diffw
   real, PARAMETER :: diffwc = 1.60 ! diffw/diffc = H2O/CO2 diffusivity
   real, PARAMETER :: rhow   = 1000.0 ! liquid water density   [kg/m3]
   real, PARAMETER :: emleaf = 1.0  ! leaf emissivity
   real, PARAMETER :: emsoil = 1.0  ! soil emissivity
   real, PARAMETER :: crd = 0.3     ! element drag coefficient
   real, PARAMETER :: csd = 0.003   ! substrate drag coefficient
!#  ifdef ONLINE_UM
!   real, PARAMETER :: cr = 0.3     ! element drag coefficient
!   real, PARAMETER :: cs = 0.003   ! substrate drag coefficient
!#  endif
   real, PARAMETER :: beta2 = crd/csd ! ratio cr/cs
   real, PARAMETER :: ccd   = 15.0  ! constant in d/h equation
   real, PARAMETER :: ccw_c = 2.0   ! ccw=(zw-d)/(h-d)
   real, PARAMETER :: usuhm = 0.3   ! (max of us/uh)
   ! Turbulence parameters:
   integer, PARAMETER :: niter = 4  ! number of iterations for za/L
   real, PARAMETER :: zetmul = 0.4  ! if niter=2, final zeta=zetmul*zetar(2)
   real, PARAMETER :: zeta0  = 0.0  ! initial value of za/L
!#  ifdef ONLINE_UM
   real, PARAMETER :: zetneg = -15.0 ! negative limit on za/L when niter>=3
   real, PARAMETER :: zetpos = 1.0  ! positive limit on za/L when niter>=3
!#  else
!   real, PARAMETER :: zetneg = -10.0 ! negative limit on za/L when niter>=3
!   real, PARAMETER :: zetpos = 0.5  ! positive limit on za/L when niter>=3
!#  endif
   real, PARAMETER :: zdlin  = 1.0  ! height frac of d below which TL linear
!#  ifdef ONLINE_UM
   real, PARAMETER :: umin   = 0.01
!#  else
!   real, PARAMETER :: umin   = 1.0
!#  endif

end module physical_constants
