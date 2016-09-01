MODULE Constants
!-------------------------------------------------------------------------------
! * Physical constants and parameters
! * MRR, 2004
!-------------------------------------------------------------------------------
USE TypeDef
implicit none
save
real(sp),parameter:: Pi    = 3.141592653589793238462643383279502884197
real(sp),parameter:: TwoPi = 6.283185307179586476925286766559005768394
real(sp),parameter:: Rgas      = 8.3143   ! universal gas constant[J/mol/K]
real(sp),parameter:: RMA       = 0.02897  ! molecular wt dry air [kg/mol]
real(sp),parameter:: RMW       = 0.018016 ! molecular wt of water [kg/mol]
real(sp),parameter:: RMC       = 0.012000 ! atomic wt of C [kg/mol]
real(sp),parameter:: SBoltz    = 5.67e-8  ! Stefan-Boltzmann constant [W/m2/K4]
real(sp),parameter:: Rlat      = 44140.0  ! lat heat evap H2O at 20C [J/molW]
real(sp),parameter:: Capp      = 29.09    ! isobaric spec heat air [J/molA/K]
real(sp),parameter:: RhoW      = 55506.0  ! liquid water density [molW/m3]
real(sp),parameter:: ViscW     = 1.002    ! water viscosity, 25 C [Pa s]
real(sp),parameter:: Grav      = 9.81     ! gravity acceleration [m/s2]
real(sp),parameter:: VonK      = 0.4      ! von Karman constant [-]
real(sp),parameter:: SecDay    = 86400.0  ! seconds/day [-]
real(sp),parameter:: QMJSolar  = 2.0      ! (molPAR)/(MJsolar) [molQ/MJ]
real(sp),parameter:: SolCMJday = 118.4    ! solar constant [MJ/m2/d]
real(sp),parameter:: SolCWm2   = 1370.0   ! solar constant [W/m2]
real(sp),parameter:: REarth    = 6.37e6   ! average earth radius [m]
REAL(sp), PARAMETER :: rlam =  2.5104e6  ! latent heat of evaporation (J kg-1)
! edit vh 13/02/08
END MODULE Constants

