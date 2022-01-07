MODULE cable_other_constants_mod

USE grid_constants_mod_cbl, ONLY : nrb, nsl, nsCs, nvCs
USE grid_constants_mod_cbl, ONLY : msn =>  nsnl

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
! Description:
!   Other CABLE constants
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

INTEGER, PARAMETER ::                                                          &
  swb = 2,           & ! 2 shortwave bands (initial division - visible /
                       ! near infrared)
  n_sw_bands = 4,    & ! total number of shortwave radiation bands
                       ! (above divided into direct / diffuse)
  mf = 2,            & ! types of leaves (sunlit / shaded)
  r_2  = SELECTED_REAL_KIND(12, 50), &!this will be removed
                       ! double precision real dimension
  niter = 4,         & ! number of iterations for za/L
  n_assim_rates = 3, & ! Rubisco, RuBP and Sink-limited rates of photosynthesis
  n_soiltypes = 9      ! number of soil types

REAL, PARAMETER ::                                                             &
  max_snow_depth = 50000.0,  & ! maximum depth of lying snow on tiles (kg/m2)
  init_snow_rho1l = 140.0      ! Initial value for snow mean density
! Gaussian integ. weights
!REAL, PARAMETER :: gauss_w(nrb)=(/0.308,0.514,0.178 /) ! F90 
REAL, PARAMETER :: gauss_w(nrb)=[0.308,0.514,0.178 ]    ! F03
REAL, PARAMETER :: rad_thresh = 0.001
                        ! minimum zenithal angle for downward SW radiation
REAL, PARAMETER :: lai_thresh = 0.001
                        ! threshold for minimum significant LAI

! minimum (cosine)zenith angle of sun signalling sunrise 
REAL, PARAMETER :: coszen_tols = 1.0e-4

REAL, PARAMETER :: z0surf_min = 1.e-7 ! min. roughness of bare soil surface

!H!REAL, PARAMETER :: z0snow_min = 1.e-7 ! min. roughness of bare snow surface

  !Meng! leaf reflectance/transmittance - comes from pft_params
  !Meng!REAL,    PARAMETER, DIMENSION(nrb) :: refl    = (/ 0.1, 0.425, 0.02 /) ! mar08
  !Meng!REAL,    PARAMETER, DIMENSION(nrb) :: taul    = (/ 0.1, 0.425, 0.02 /) ! mar08

  ! Are these all Meng's parameters: rm KIND dependence
  INTEGER, PARAMETER                 :: istemp  = 4                      ! soil temp:     1,2,3,4 = FR,kf,mrr,mrrkf
  INTEGER, PARAMETER                 :: ismois  = 2                      ! soil moist:  1,2,3     = MP84,NP89,Richards
  INTEGER, PARAMETER                 :: isinf   = 2                      ! soil infilt: 1,2       = MP84, FC96
  INTEGER, PARAMETER                 :: isevap  = 2                      ! soil evap: 1,2,3 = alfa,beta,threshold
  INTEGER, PARAMETER                 :: itherm  = 1                      ! VW or KGK algorithm for hconds,rkapps
  INTEGER, PARAMETER                 :: irktem  = 5                      ! RK steps in soil temp schemes
  INTEGER, PARAMETER                 :: irkmoi  = 5                      ! RK steps in soil moisture schemes
  ! soil water  parameters:
  REAL,    PARAMETER                 :: etarct  = 0.7                    ! rel soil moisture for finding zst1,zst2
  REAL,    PARAMETER                 :: dbde    = 1.3333                 ! d(beta)/d(etar): 4/3, 8 for D78,KGK91
  INTEGER, PARAMETER                 :: istsw   = 1                      !
  INTEGER, PARAMETER                 :: iresp   = 0                      ! unscaled (iresp=0) or scaled (iresp=1) respiration

END MODULE cable_other_constants_mod
