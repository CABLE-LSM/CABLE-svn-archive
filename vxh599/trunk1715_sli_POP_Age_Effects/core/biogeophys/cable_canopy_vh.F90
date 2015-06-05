! cable_canopy_vh.f90
!********************************************************************************
! Edited by VH and MC 16/10/09: Moved incanopy T to top of iteration loop;
! moved stability calc to top of iteration loop;
! removed restriction on fevc due to extractible water;
! introduced tolerances for qvair and Tvair to reduce number of iterations
! merge wet canopy calc into leaf T iteration loop
! recalc wetbal and drybal
! fwsoil modified to max(alpha2)
! initialised all local variables
!
! 6/6/2012: modified to convert input screen t and q to values at ref height (for BIOS2)
!
! 24/6/2012: convert to a single loop, with Tvair and qvair updated once using fluxes from prev. timestep.
!
! MC 2012: make all carbon double precision and water single precision
!          works with all double precision: included interface routines
!
! VH March 2014: energy closure checks using both sli and soilsnow.
! Note small changes to cable_radiation.F90 required for energy closure.
! If soil model is sli, set fwsoil to canopy%fwsoil, which is calculated in
! getrex (subroutine called from cable_sli_main).

! n.b. No testing for over-extraction of water by soil-snow. Suggest using getrex (in cable_sil_utils)
! to caclulate both fwsoil and root extraction, as implemented for sli. This avoids over-extraction.
! Currently wbal is non-zero when cable_user%soil_struc=='default'. Above suggestion should solve this.

!********************************************************************************

! Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach,
! Ray Leuning, Vanessa Haverd, et al. at CSIRO Marine and Atmospheric Research.
!
! Fortran-95 coding by Harvey Davies, Gab Abramowitz, Martin Dix, and Matthias Cuntz
! bugs to gabsun@gmail.com.
!
! This file contains module canopy_module and subroutine define_canopy.
! The functions included are:
!   qsatf,
!   ej3
!   ej4x,
!   xvcmxt4,
!   xvcmxt3,
!   xejmxt3,
!   psim,
!   psis,
!   rplant, and
!   rsoil
!
MODULE canopy_vh_module

  USE photosynthetic_constants, ONLY: trefk, alpha3, alpha4, convx3, convx4, &
       a1c3_default, a1c4_default, d0c3_default, d0c4_default, d0c3_hawkesbury, &
       gsw03, gsw04, conkc0, conko0, ekc, eko, gam0, gam1, gam2, cfrd3, cfrd4, &
       rgswc, rgbwc, maxiter
  USE physical_constants,       ONLY: tfrz, rmh2o, rmair, tetena, tetenb, tetenc, &
       rgas, zeta0, zetpos, zetneg, vonk, grav, capp, rhow, apol, prandt, emleaf, &
       emsoil, sboltz, umin, dheat, niter
  USE cable_radiation_module,         ONLY: radiation
  USE cable_air_module,               ONLY: define_air
  USE cable_def_types_mod,             ONLY: air_type, balances_type, bgc_pool_type, canopy_type, &
       met_type, radiation_type, roughness_type, ms, mp, mf, &
       soil_parameter_type, soil_snow_type, veg_parameter_type, i_d, r_2
  USE cable_common_module,  ONLY: cable_user
  USE cable_data_module, ONLY : icanopy_type, point2constants







  IMPLICIT NONE

  PRIVATE

  PUBLIC :: define_canopy_vh

  INTERFACE qsatf
     MODULE PROCEDURE qsatf_r_1
     MODULE PROCEDURE qsatf_r_2
  END INTERFACE qsatf
  INTERFACE ej3x
     MODULE PROCEDURE ej3x_r_1
     MODULE PROCEDURE ej3x_r_2
  END INTERFACE ej3x
  INTERFACE ej4x
     MODULE PROCEDURE ej4x_r_1
     MODULE PROCEDURE ej4x_r_2
  END INTERFACE ej4x
  INTERFACE xvcmxt4
     MODULE PROCEDURE xvcmxt4_r_1
     MODULE PROCEDURE xvcmxt4_r_2
  END INTERFACE xvcmxt4
  INTERFACE xvcmxt3
     MODULE PROCEDURE xvcmxt3_r_1
     MODULE PROCEDURE xvcmxt3_r_2
  END INTERFACE xvcmxt3
  INTERFACE xejmxt3
     MODULE PROCEDURE xejmxt3_r_1
     MODULE PROCEDURE xejmxt3_r_2
  END INTERFACE xejmxt3
  INTERFACE psim
     MODULE PROCEDURE psim_r_1
     MODULE PROCEDURE psim_r_2
  END INTERFACE psim
  INTERFACE psis
     MODULE PROCEDURE psis_r_1
     MODULE PROCEDURE psis_r_2
  END INTERFACE psis
  INTERFACE rplant
     MODULE PROCEDURE rplant_r_1
     MODULE PROCEDURE rplant_r_2
  END INTERFACE rplant
  INTERFACE rsoil
     MODULE PROCEDURE rsoil_r_1
     MODULE PROCEDURE rsoil_r_2
  END INTERFACE rsoil

  TYPE( icanopy_type ) :: C
  REAL, PARAMETER :: xvccoef      = 1.17461 !derived parameter
  REAL, PARAMETER :: EHaVc        = 73637.0  !J/mol (Leuning 2002)
  REAL, PARAMETER :: EHdVc        = 149252.0 !J/mol (Leuning 2002)
  REAL, PARAMETER :: EntropVc     = 486.0  !J/mol/K (Leuning 2002)
  !   xvccoef=1.0+exp(EHdjx)/(Rconst*TrefK))
  REAL, PARAMETER :: EHaJx        = 50300.0  !J/mol (Leuning 2002)
  REAL, PARAMETER :: EHdJx        = 152044.0    !J/mol (Leuning 2002)
  REAL, PARAMETER :: EntropJx     = 495.0     !J/mol/K (Leuning 2002)
  REAL, PARAMETER :: xjxcoef      = 1.16715    !derived parameter
  REAL(r_2), PARAMETER :: tfrz_r_2     = real(tfrz,r_2)
  REAL(r_2), PARAMETER :: trefk_r_2    = real(trefk,r_2)
  REAL(r_2), PARAMETER :: rmh2o_r_2    = real(rmh2o,r_2)
  REAL(r_2), PARAMETER :: rmair_r_2    = real(rmair,r_2)
  REAL(r_2), PARAMETER :: tetena_r_2   = real(tetena,r_2)
  REAL(r_2), PARAMETER :: tetenb_r_2   = real(tetenb,r_2)
  REAL(r_2), PARAMETER :: tetenc_r_2   = real(tetenc,r_2)
  REAL(r_2), PARAMETER :: alpha3_r_2   = real(alpha3,r_2)
  REAL(r_2), PARAMETER :: convx3_r_2   = real(convx3,r_2)
  REAL(r_2), PARAMETER :: alpha4_r_2   = real(alpha4,r_2)
  REAL(r_2), PARAMETER :: convx4_r_2   = real(convx4,r_2)
  REAL(r_2), PARAMETER :: xvccoef_r_2  = real(xvccoef,r_2)
  REAL(r_2), PARAMETER :: ehavc_r_2    = real(ehavc,r_2)
  REAL(r_2), PARAMETER :: rgas_r_2     = real(rgas,r_2)
  REAL(r_2), PARAMETER :: entropvc_r_2 = real(entropvc,r_2)
  REAL(r_2), PARAMETER :: ehdvc_r_2    = real(ehdvc,r_2)
  REAL(r_2), PARAMETER :: xjxcoef_r_2  = real(xjxcoef,r_2)
  REAL(r_2), PARAMETER :: ehajx_r_2    = real(ehajx,r_2)
  REAL(r_2), PARAMETER :: entropjx_r_2 = real(entropjx,r_2)
  REAL(r_2), PARAMETER :: ehdjx_r_2    = real(ehdjx,r_2)

CONTAINS

  !--------------------------------------------------------------------------
  SUBROUTINE define_canopy_vh(ktau,bal,rad,rough,air,met,dels,ssnow,soil, &
       veg,canopy)


    TYPE (balances_type),        INTENT(INOUT) :: bal
    TYPE (radiation_type),       INTENT(INOUT) :: rad
    TYPE (roughness_type),       INTENT(INOUT) :: rough
    TYPE (air_type),             INTENT(INOUT) :: air
    TYPE (met_type),             INTENT(INOUT) :: met
    REAL,                        INTENT(IN)    :: dels ! integration time setp (s)
    TYPE (soil_snow_type),       INTENT(INOUT) :: ssnow
    TYPE (soil_parameter_type),  INTENT(INOUT) :: soil
    TYPE (veg_parameter_type),   INTENT(INOUT) :: veg
    TYPE (canopy_type),          INTENT(INOUT) :: canopy
    INTEGER,                     INTENT(IN)    :: ktau ! integration step number


    REAL(r_2), DIMENSION(mp,mf)   :: a1c3 ! Spatially varying a1c3
    REAL(r_2), DIMENSION(mp,mf)   :: a1c4 ! Spatially varying a1c3
    REAL(r_2), DIMENSION(mp,mf)   :: abs_deltlf ! ABS(deltlf)
    REAL(r_2), DIMENSION(mp,mf,3) :: ancj ! soln to quad eqn
    REAL(r_2), DIMENSION(mp,mf)   :: anx ! net photos. prev iteration
    REAL(r_2), DIMENSION(mp,mf)   :: an_y ! net photosynthesis soln
    REAL(r_2), DIMENSION(mp,mf)   :: ca2   ! 2D CO2 concentration
    REAL(r_2), DIMENSION(mp)      :: cansat ! max canopy intercept. (mm)
    REAL(r_2), DIMENSION(mp,mf,3) :: ci ! intercellular CO2 conc.
    REAL(r_2), PARAMETER                :: co2cp3 = 0.0 ! CO2 compensation pt C3
    REAL(r_2), DIMENSION(mp,mf,3) :: coef0 ! CO2 comp. pt coeff 1
    REAL(r_2), DIMENSION(mp,mf,3) :: coef1 ! " 2
    REAL(r_2), DIMENSION(mp,mf,3) :: coef2 ! " 3
    REAL(r_2), DIMENSION(mp,mf)   :: conkct ! Michaelis Menton const.
    REAL(r_2), DIMENSION(mp,mf)   :: conkot ! Michaelis Menton const.
    REAL(r_2), DIMENSION(mp)      :: cscale ! scaling between 2 hawkesbury elev co2 vals
    REAL(r_2), DIMENSION(mp,mf)   :: csx ! leaf surface CO2 concentration
    REAL(r_2), DIMENSION(mp,mf,3) :: cx  ! "d_{3}" in Wang and Leuning, 1998, appendix E
    REAL(r_2), DIMENSION(mp,mf)   :: da2 ! 2D sat vap pres deficit
    REAL(r_2), DIMENSION(mp,mf,3) :: delcx ! discriminant  in quadratic in eq. E7 Wang and Leuning, 1998
    REAL(r_2), DIMENSION(mp,mf)   :: deltlf ! deltlfy of prev iter.
    REAL(r_2), DIMENSION(mp,mf)   :: deltlfy ! del temp successive iteration

    REAL(r_2), DIMENSION(mp,mf)    :: dsatdk2      ! 2D dsatdk
    REAL(r_2), DIMENSION(mp,mf)    :: dsx ! leaf surface vpd
    REAL(r_2), DIMENSION(mp,mf)    :: ecx ! lat. hflux big leaf
    REAL(r_2), DIMENSION(mp,mf)    :: ejmax2 ! jmax of big leaf
    REAL(r_2), DIMENSION(mp,mf)    :: ejmxt3 ! jmax big leaf C3 plants
    REAL(r_2), DIMENSION(mp,mf)    :: ecy ! lat heat fl dry big leaf
    REAL(r_2), DIMENSION(mp,mf)    :: d0c3 ! Empirical coef for vpd sensitivity of stomata
    REAL(r_2), DIMENSION(mp,mf)    :: d0c4 ! Empirical coef for vpd sensitivity of stomata
    REAL(r_2), DIMENSION(mp,mf)    :: frac42       ! 2D frac4
    REAL(r_2), DIMENSION(mp,mf)    :: fwsoil2 ! soil water modifier of stom. cond.
    REAL(r_2), DIMENSION(mp,mf)    :: gbhf ! freeConvectionBndLayerConductance mol/m2/s
    REAL(r_2), DIMENSION(mp,mf)    :: gbhu ! forcedConvectionBoundaryLayerConductance
    REAL(r_2), DIMENSION(mp)       :: gbvtop ! bnd layer cond. top leaf
    REAL(r_2), DIMENSION(mp,mf)    :: gras ! Grashof coeff
    REAL(r_2), DIMENSION(mp,mf)    :: gswmin ! min stomatal conductance
    REAL(r_2), DIMENSION(mp,mf)    :: gswx ! stom cond for water
    REAL(r_2), DIMENSION(mp,mf)    :: gw  ! cond for water for a dry canopy
    REAL(r_2), DIMENSION(mp,mf)    :: gh  ! cond for heat for a dry canopy
    REAL(r_2), DIMENSION(mp,mf)    :: ghr ! dry canopy cond for heat & thermal radiat'n
    REAL(r_2), DIMENSION(mp,mf)    :: ghwet  ! cond for heat for a wet canopy
    REAL(r_2), DIMENSION(mp,mf)    :: gbw ! cond for water transfer (leaf boundary layer)
    REAL(r_2), DIMENSION(mp,mf)    :: hcx ! sens heat fl big leaf prev iteration
    REAL(r_2), DIMENSION(mp,mf)    :: hcy ! veg. sens heat
    INTEGER(i_d)                         :: iter ! iteration #
    INTEGER(i_d)                         :: iterplus !
    INTEGER(i_d)                         :: k            ! iteration count
    INTEGER(i_d)                         :: kk           ! iteration count
    REAL(r_2)                            :: rk ! k in real
    REAL(r_2), DIMENSION(mp,mf)    :: psycst ! modified pych. constant
    REAL(r_2), DIMENSION(mp,mf)    :: rdx ! daytime leaf resp rate, prev iteration
    REAL(r_2), DIMENSION(mp,mf)    :: rdy ! daytime leaf resp rate
    REAL(r_2), DIMENSION(mp,mf)    :: rnx ! net rad prev timestep
    REAL(r_2), DIMENSION(mp,mf)    :: rny ! net rad
    REAL(r_2), DIMENSION(mp)       :: ortsoil ! turbulent resistance, prev time step
    REAL(r_2), DIMENSION(mp,mf)    :: tair2 ! 2D tair
    REAL(r_2), DIMENSION(mp,mf)    :: tvair2 ! 2D tair
    REAL(r_2), DIMENSION(mp,mf)    :: tdiff ! leaf air temp diff.
    REAL(r_2), DIMENSION(mp,mf)    :: tlfx ! leaf temp prev. iteration
    REAL(r_2), DIMENSION(mp,mf)    :: tlfxx ! leaf temperature of current iteration
    REAL(r_2), DIMENSION(mp,mf)    :: tlfy ! leaf temp
    REAL(r_2), DIMENSION(mp,mf)    :: vcmax2 ! vcmax big leaf
    REAL(r_2), DIMENSION(mp,mf)    :: vcmxt3 ! vcmax big leaf C3
    REAL(r_2), DIMENSION(mp,mf)    :: vcmxt4 ! vcmax big leaf C4
    REAL(r_2), DIMENSION(mp,mf,2)  :: vx3 ! carboxylation C3 plants
    REAL(r_2), DIMENSION(mp,mf,2)  :: vx4 ! carboxylation C4 plants
    REAL(r_2), DIMENSION(mp,mf)    :: xdleaf2      ! 2D dleaf
    REAL(r_2), DIMENSION(mp,mf)    :: xleuning ! leuning stomatal coeff
    REAL(r_2), DIMENSION(mp,mf)    :: temp ! vcmax big leaf C3
    REAL(r_2), DIMENSION(mp,mf)    :: fwet ! fraction wet canopy
    REAL(r_2), DIMENSION(mp,mf)    :: Ecansto  ! supply limited evap from wet leaves
    ! Bonan,LSM version 1.0, p106)
    REAL, DIMENSION(mp)       :: oldcansto ! prev t step canopy storage
    REAL(r_2), DIMENSION(mp)  :: cc ! limitation term for canopy interception per timestep
    REAL, DIMENSION(mp,niter) :: zetar, zetash ! stability correction
    REAL(r_2), PARAMETER      :: jtomol = 4.6e-6 ! Conversion from Joule to Mol for light
    REAL(r_2), PARAMETER      :: effc4 = 4000.0  !Vc=effc4*Ci*Vcmax (see
    REAL, DIMENSION(mp)       :: fwsoil ! soil water modifier of stom. cond.
    REAL, DIMENSION(mp)       :: rt0 ! turbulent resistance
    REAL, DIMENSION(mp)       :: rt1usc ! eq. 3.53, SCAM manual, 1997
    REAL, DIMENSION(mp)       :: denom ! denominator in calculating screen temperature, humidity etc
    REAL, DIMENSION(mp)       :: tstar !
    REAL, DIMENSION(mp)       :: zscrn !
    REAL, DIMENSION(mp)       :: qstar !
    REAL, DIMENSION(mp)       :: rsts  !
    REAL, DIMENSION(mp)       :: qsurf !
    REAL, DIMENSION(mp)       :: qtgnet !
    REAL, DIMENSION(mp,ms)    :: tmp2d1, tmp2d2
    REAL, DIMENSION(mp)       :: phenps ! Leaf phenology influence on vcmax and jmax
    REAL, DIMENSION(mp)       :: poolcoef1 ! leaf carbon turnover rate * leaf pool size
    REAL, DIMENSION(mp)       :: poolcoef1w ! wood carbon turnover rate * wood pool size
    REAL, DIMENSION(mp)       :: poolcoef1r ! root carbon turnover rate * root pool size
    REAL, DIMENSION(mp)            :: rbw ! leaf boundary layer resistance for water
    REAL, DIMENSION(mp)            :: rrbw ! recipr. leaf boundary layer resistance for water
    REAL, DIMENSION(mp)            :: rsw ! stomatal resistance for water
    REAL, DIMENSION(mp)            :: rrsw ! recipr. stomatal resistance for water
    REAL, DIMENSION(mp)            :: dmah ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL, DIMENSION(mp)            :: dmbh ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL, DIMENSION(mp)            :: dmch ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL, DIMENSION(mp)            :: dmae ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL, DIMENSION(mp)            :: dmbe ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL, DIMENSION(mp)            :: dmce ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL, DIMENSION(mp)       :: tss4 ! soil/snow temperature**4
    REAL, DIMENSION(mp)       :: sss ! variable for Penman-Monteith evap for soil
    REAL, DIMENSION(mp)       :: cc1 ! variable for Penman-Monteith evap for soil
    REAL, DIMENSION(mp)       :: cc2 ! variable for Penman-Monteith evap for soil
    REAL, DIMENSION(mp)       :: qstvair ! sat spec hunidity at leaf temperature
    REAL, DIMENSION(mp)       :: xx ! delta-type function for sparse canopy limit, p20 SCAM manual
    !%% changes by Ashok Luhar (low wind speed)
    REAL, PARAMETER                 :: alpha1 = 4.0
    REAL, PARAMETER                 :: beta1  = 0.5
    REAL, PARAMETER                 :: gamma1 = 0.3
    REAL, DIMENSION(mp)       :: zeta1
    REAL, DIMENSION(mp)       :: zeta2
    !**************************************************************************************** ! vh 17/07/09
    REAL, DIMENSION(mp,ms)    :: alpha2a_root, alpha2_root, delta_root
    !**************************************************************************************** ! vh 17/07/09
    LOGICAL, DIMENSION(mp,mf)      :: Flag_fwet
    LOGICAL, DIMENSION(mp)         :: mdb_mask ! Needs to be set
    REAL, DIMENSION(mp)       :: tmp1d
    REAL, DIMENSION(mp)       :: zstar, rL, phist, csw, psihat,rt0bus

    ! assign local ptrs to constants defined in cable_data_module
    CALL point2constants(C)

    a1c3         = 0.0_r_2
    a1c4         = 0.0_r_2
    abs_deltlf   = 0.0_r_2
    ancj         = 0.0_r_2
    anx          = 0.0_r_2
    an_y         = 0.0_r_2
    ca2          = 0.0_r_2
    cansat       = 0.0_r_2
    ci           = 0.0_r_2
    coef0        = 0.0_r_2
    coef1        = 0.0_r_2
    coef2        = 0.0_r_2
    conkct       = 0.0_r_2
    conkot       = 0.0_r_2
    cscale       = 0.0_r_2
    csx          = 0.0_r_2
    cx           = 0.0_r_2
    da2          = 0.0_r_2
    delcx        = 0.0_r_2
    deltlf       = 0.0_r_2
    deltlfy      = 0.0_r_2
    dsatdk2      = 0.0_r_2
    dsx          = 0.0_r_2
    ecx          = 0.0_r_2
    ejmax2       = 0.0_r_2
    ejmxt3       = 0.0_r_2
    ecy          = 0.0_r_2
    d0c3         = 0.0_r_2
    d0c4         = 0.0_r_2
    frac42       = 0.0_r_2
    fwsoil2      = 0.0_r_2
    gbhf         = 0.0_r_2
    gbhu         = 0.0_r_2
    gbvtop       = 0.0_r_2
    gras         = 0.0_r_2
    gswmin       = 0.0_r_2
    gswx         = 0.0_r_2
    gw           = 0.0_r_2
    gh           = 0.0_r_2
    ghr          = 0.0_r_2
    ghwet        = 0.0_r_2
    gbw          = 0.0_r_2
    hcx          = 0.0_r_2
    hcy          = 0.0_r_2
    iter         = 1_i_d
    iterplus     = 0_i_d
    k            = 0_i_d
    kk           = 0_i_d
    rk           = 0.0_r_2
    psycst       = 0.0_r_2
    rdx          = 0.0_r_2
    rdy          = 0.0_r_2
    rnx          = 0.0_r_2
    rny          = 0.0_r_2
    ortsoil      = 0.0_r_2
    tair2        = 0.0_r_2
    tvair2       = 0.0_r_2
    tdiff        = 0.0_r_2
    tlfx         = 0.0_r_2
    tlfxx        = 0.0_r_2
    tlfy         = 0.0_r_2
    vcmax2       = 0.0_r_2
    vcmxt3       = 0.0_r_2
    vcmxt4       = 0.0_r_2
    vx3          = 0.0_r_2
    vx4          = 0.0_r_2
    xdleaf2      = 0.0_r_2
    xleuning     = 0.0_r_2
    temp         = 0.0_r_2
    fwet         = 0.0_r_2
    Ecansto      = 0.0_r_2
    oldcansto    = 0.0
    cc           = 0.0_r_2
    zetar        = 0.0
    fwsoil       = 0.0
    rt0          = 0.0
    rt1usc       = 0.0
    denom        = 0.0
    tstar        = 0.0
    zscrn        = 0.0
    qstar        = 0.0
    rsts         = 0.0
    qsurf        = 0.0
    qtgnet       = 0.0
    tmp2d2       = 0.0
    tmp2d1       = 0.0
    phenps       = 0.0
    poolcoef1    = 0.0
    poolcoef1w   = 0.0
    poolcoef1r   = 0.0
    tss4         = 0.0
    sss          = 0.0
    cc1          = 0.0
    cc2          = 0.0
    qstvair      = 0.0
    xx           = 0.0
    zeta1        = 0.0
    zeta2        = 0.0
    alpha2a_root = 0.0
    alpha2_root  = 0.0
    delta_root   = 0.0
    Flag_fwet    = .false.
    mdb_mask     = .false.
    tmp1d        = 0.0

    a1c3 = a1c3_default
    a1c4 = a1c3_default
    d0c3 = d0c3_default
    d0c4 = d0c4_default

    a1c3(:,1) = veg%a1c3
    a1c3(:,2) = a1c3(:,1)
    a1c4(:,1) = a1c4_default ! set in photosynthetic constants

    a1c4(:,2) = a1c4(:,1)
    d0c3(:,1) = veg%d0c3
    d0c3(:,2) = d0c3(:,1)
    d0c4(:,1) = veg%d0c3
    d0c4(:,2) = d0c4(:,1)
    canopy%fevw = 0.0
    canopy%delwc = 0.0



    ! Set surface water vapour pressure deficit:
    met%da = (qsatf(met%tk-tfrz,met%pmb) - met%qv )*rmair/rmh2o*met%pmb*100.0

    ! Lai and Katul formulation for root efficiency function  vh 17/07/09
    alpha2a_root = real(max(ssnow%wb-soil%swilt_vec, 0.001_r_2)/(soil%ssat_vec))
    tmp2d1 = REAL(ssnow%wb -soil%swilt_vec)
    tmp2d2 = SPREAD(real(veg%gamma),2,ms)/tmp2d1*log(alpha2a_root)
    WHERE ((tmp2d1>1.e-8) .and. (tmp2d2 > -10.0))
       alpha2_root = exp(tmp2d2)
    ELSEWHERE
       alpha2_root = 0.0
    ENDWHERE

    WHERE (veg%froot>0.0)
       delta_root = 1.0
    ELSEWHERE
       delta_root = 0.0
    ENDWHERE

    fwsoil  = maxval(alpha2_root*delta_root, 2)
    fwsoil  = max(0.0, fwsoil)
    IF (cable_user%soil_struc=='sli') fwsoil = real(canopy%fwsoil)

    fwsoil2 = real(SPREAD(fwsoil,2,mf),r_2)


    ! BATS-type canopy saturation proportional to LAI:
    cansat = real(veg%canst1*canopy%vlaiw,r_2)
    ! Leaf phenology influence on vcmax and jmax
    ! rml 22/10/07 only apply to deciduous types
    WHERE (veg%deciduous)
       phenps = max(1.0e-4, MIN(1.,1. - ( (veg%tmaxvj - real(ssnow%tgg(:,4))+tfrz)/ &
            (veg%tmaxvj - veg%tminvj) )**2 ) )
       WHERE ( ssnow%tgg(:,4) < real(veg%tminvj + tfrz,r_2) ) phenps = 0.0
       WHERE ( ssnow%tgg(:,4) > real(veg%tmaxvj + tfrz,r_2) ) phenps = 1.0
    ELSEWHERE
       phenps = 1.0
    END WHERE

    ! Set previous time step canopy water storage:
    oldcansto = canopy%cansto
    ! to avoid excessive direct canopy evaporation, rainfall rate is limited,
    ! hence canopy interception is limited (EK nov2007, snow scheme)
    ! modified further by timestep requirement to avoid canopy temperature
    ! oscillations (EAK aug08)
    cc1 = MIN(met%precip-met%precip_sn, 4.*MIN(dels,1800.)/60./1440.)
    ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
    canopy%wcint = MERGE(MIN(MAX(real(cansat) - canopy%cansto,0.0),cc1), 0.0, cc1>0.0)  ! EK nov2007, snow scheme
    !         cc > 0.0  .AND. met%tk > tfrz)

    ! Delete line below in case it affects snow sites (precip_s) (EK Jul08)
    !    canopy%through = MIN(met%precip,MAX(0.0, met%precip - canopy%wcint))

    ssnow%wetfac = MAX(0.0, MIN(1.0, &
         (REAL(ssnow%wb(:,1)) - soil%swilt) / (soil%sfc - soil%swilt)))
    ! owetfac introduced to reduce sharp changes in dry regions,
    ! especially in offline runs where there may be discrepancies between
    ! timing of precip and temperature change (EAK apr2009)
    !ssnow%wetfac = 0.5*(ssnow%wetfac + ssnow%owetfac)
    ! Temporay fixer for accounting the reduction of soil evap due to freezing
    WHERE ( ssnow%wbice(:,1) > 0.0_r_2 ) ! Prevents divide by zero at glaciated
       ! points where wb and wbice=0.
       ssnow%wetfac = ssnow%wetfac &
            * (1.0 - REAL(ssnow%wbice(:,1)/ssnow%wb(:,1)))**2 !! only used if soilsnow
    END WHERE
    !ssnow%wetfac = 1.0 ! test vh!
    zetar(:,1) = zeta0 ! stability correction terms
    zetar(:,2) = zetpos + 1.0
    xdleaf2 = real(SPREAD(veg%dleaf,2,mf),r_2) ! characteristic leaf length
    ca2 = real(SPREAD(met%ca,2,mf),r_2)        ! CO2 concentration
    csx = ca2                     ! initialise leaf surface CO2 concentration
    da2 = real(SPREAD(met%da,2,mf),r_2) ! water vapour pressure deficit
    dsx = da2                     ! init. leaf surface vpd
    tair2  = real(SPREAD(met%tvair-tfrz,2,mf),r_2) ! air temp in C
    ejmax2 = real(SPREAD(veg%ejmax*phenps,2,mf),r_2) !max. pot. electr transp. rate top leaf(mol/m2s)
    vcmax2 = real(SPREAD(veg%vcmax*phenps,2,mf),r_2) !max. RuBP carboxylsation rate top leaf(mol/m2s)
    tlfx = tair2  ! initialise leaf temp iteration memory variable
    tlfy = tair2  ! initialise current leaf temp
    frac42 = real(SPREAD(veg%frac4,2,mf),r_2) ! frac C4 plants
    ! weight min stomatal conductance by C3 an C4 plant fractions
    rdy  = 0.0_r_2       ! init daytime leaf respiration rate
    rdx  = 0.0_r_2       ! init daytime leaf respiration rate
    an_y = 0.0_r_2       ! init current estimate net photos.
    gswx = 1e-15_r_2     ! default stomatal conuctance
    gbhf = 1e-3_r_2      ! default free convection boundary layer conductance
    gbhu = 1e-3_r_2      ! default forced convection boundary layer conductance
    anx  = 0.0_r_2       ! init net photos. iteration memory variable
    ancj = 0.0_r_2
    !    psycst = SPREAD(air%psyc, 2, mf) ! modified psyc. constant
    ! add on by ypw 1-oct-2002
    gw    = 1.0e-3_r_2 ! default values of conductance
    gh    = 1.0e-3_r_2
    ghr   = 1.0e-3_r_2
    ghwet = 1.0e-3_r_2
    ! Initialise in-canopy temperatures and humidity:
    met%tvair = met%tk
    met%tvrad = met%tk
    met%qvair = met%qv
    ortsoil   = ssnow%rtsoil
    IF (cable_user%soil_struc=='default') then
       ssnow%tss =  real((1-ssnow%isflag))*ssnow%tgg(:,1) + real(ssnow%isflag)*ssnow%tggsn(:,1)
    elseif (cable_user%soil_struc=='sli') then
       ssnow%tss = real(ssnow%Tsurface) + tfrz
    endif
    tss4      = ssnow%tss**4

    ! Calculate fraction of canopy which is wet:
    canopy%fwet = MAX(0.0,MIN(1.0,real(0.8_r_2*canopy%cansto/MAX(cansat,0.01_r_2))))
    fwet = real(SPREAD(canopy%fwet,2,mf),r_2)
    iter = 0


    CALL define_air(met, air)
    dsatdk2 = real(SPREAD(air%dsatdk,2,mf),r_2)! deriv of sat vap pressure wrt temp
    CALL radiation(ssnow, veg, air, met, rad, canopy)
    where (rad%fvlai > 1.e-5)
       ! supply limited evap from wet leaves
       Ecansto = fwet*SPREAD(air%rlam,2,mf)*real(rhow/(dels*1.0e3)*rad%fvlai/SPREAD(canopy%vlaiw,2,mf),r_2)
       !Ecansto = SPREAD(canopy%fwet,2,mf)*rhow*SPREAD(air%rlam,2,mf)/(dels*1.0e3) * &
       !    SPREAD(veg%canst1,2,mf) *rad%fvlai/SPREAD(canopy%vlaiw,2,mf) ! supply limited evap from wet leaves
    elsewhere
       Ecansto = 0.0_r_2
    endwhere

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! SPECIAL for BIOS2: set screen t and humidity as met inputs
    !   canopy%tscrn=   met%tk
    ! canopy%qscrn  = met%qv
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DO WHILE (iter<3) ! iterate over in-canopy T and q
       iter = iter+1
       canopy%cansto = oldcansto
       ! Add canopy interception to canopy storage term:
       canopy%cansto = canopy%cansto + canopy%wcint
       ! Define canopy throughfall (100% of precip if temp < 0C, see above):
       canopy%through = met%precip_sn + MIN( met%precip - met%precip_sn , &
            MAX(0.0, met%precip - met%precip_sn - real(canopy%wcint)) )  ! EK nov2007

       ssnow%otss(:) = ssnow%tss(:)
       tss4 = ssnow%tss**4
       CALL define_air(met, air)
       dsatdk2 = SPREAD(air%dsatdk, 2, mf)

       ! monin-obukhov stability parameter zetar=zref/l

       if (iter==1) then
          zetar(:,iter) = zeta0
          zetash(:,iter) = zeta0
          if (ktau.eq.1) then
             canopy%fh = 0.0
             canopy%fe = 0.0
          endif
       elseif (iter ==2) then
          zetar(:,iter) = -(vonk*grav*rough%zref_tq*(canopy%fh+0.07*canopy%fe))/ &
               max( (air%rho*capp*met%tk*canopy%us**3), 1.e-12)
          where (canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
             zetash(:,iter) = -(vonk*grav*(0.1*rough%hruff)*(canopy%fhs+0.07*real(canopy%fes)))/ &
                  max( (air%rho*capp*met%tk*(canopy%us*rough%term6a)**3), 1.e-12)
          elsewhere
             zetash(:,iter) = zetash(:,iter-1)
          end where
       elseif (iter ==3) then
          zetar(:,iter) = 0.4 * zetar(:,iter-1)
          zetash(:,iter) = 0.4 * zetash(:,iter-1)

          !          zetar(:,iter) = -(vonk*grav*rough%zref_tq*(canopy%fh+0.07*canopy%fe))/ &
          !               max( (air%rho*capp*met%tk*canopy%us**3), 1.e-12)
          !
          !          zetash(:,iter) = -(vonk*grav*(0.1*rough%hruff)*(canopy%fhs+0.07*canopy%fes))/ &
          !               max( (air%rho*capp*met%tk*(canopy%us*rough%term6a)**3), 1.e-12)
       endif
       ! write(*,*) "after zetar 1",  ktau, iter, zetar(:,iter)
       !        constrain zeta to zetpos and zetneg (set in param0)
       zetar(:,iter) = min(zetpos,zetar(:,iter))         ! zetar too +
       zetar(:,iter) = max(zetneg,zetar(:,iter))         ! zetar too -
       !  write(*,*) "after zetar 2",  ktau, iter, zetar(:,iter)
       !if (ktau>100) then
       ! stop
       !  endif
       gswmin = real(rad%scalex*gsw03,r_2) * (1.-frac42) + real(rad%scalex*gsw04,r_2) * frac42
       ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
       ! resistances rt0, rt1 (elements of dispersion matrix):
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
       canopy%us = MAX(1.e-6, vonk * MAX(met%ua,umin) / ( LOG(rough%zref_tq / rough%z0m) - &
            psim(zetar(:,iter)) + psim(zetar(:,iter) * rough%z0m / rough%zref_tq) ))
       ! write(*,*) ktau, iter, zetar(:,iter), rough%zref_tq, rough%disp
       WHERE (canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
          rt0bus = (LOG(0.1*rough%hruff/rough%z0soilsn) - psis(zetash(:,iter)) + &
               psis(zetash(:,iter)*rough%z0soilsn/(0.1*rough%hruff))) / &
               vonk/rough%term6a
       end where

       !%%change by Ashok Luhar - low wind formulation
       where (zetar(:,iter) > 0.7)
          zeta1=zetar(:,iter) * rough%z0m / rough%zref_tq
          canopy%us = MAX(1.e-6, vonk * MAX(met%ua,umin) / ( &
               alpha1* ((zetar(:,iter)**beta1* (1.0+gamma1*zetar(:,iter)**(1.0-beta1))) &
               - (zeta1**beta1*(1.0+gamma1*zeta1**(1.0-beta1))))))
       endwhere
   
       !%%
       ! Turbulent aerodynamic resistance from roughness sublayer depth to reference height,
       ! x=1 if zref+disp>zruffs, 0 otherwise: thus rt1usc = 0 if zref+disp<zruffs
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
       xx = 0.5 + sign(0.5,rough%zref_tq+rough%disp-rough%zruffs)
       !              correction  by Ian Harman to the 2nd psis term
       rt1usc = xx * (LOG(rough%zref_tq/MAX(rough%zruffs-rough%disp, rough%z0soilsn)) &
            - psis( zetar(:,iter) ) &
            + psis( zetar(:,iter)*(MAX(rough%zruffs-rough%disp,rough%z0soilsn))/rough%zref_tq ) &
            )/vonk


       !      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! vh 14/04/14
       ! for stable conditions, update rough%rt0us & rough%rt1usa by replacing C%CSW by
       ! csw = cd/2* (U(hc)/ust)**2 according to Eqs 15 & 19 from notes by Ian Harman (9-9-2011)
       WHERE (canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
          zstar = rough%disp + 1.5*(veg%hc - rough%disp)
          ! write(*,*) zstar, rough%disp, veg%hc,rough%disp
          psihat = log((zstar - rough%disp)/ (veg%hc - rough%disp)) + &
               (veg%hc - zstar)/(zstar - rough%disp)
          rL = -(vonk*grav*(zstar - rough%disp)*(canopy%fh))/ &  ! 1/Monin-Obokov Length
               max( (air%rho*capp*met%tk*canopy%us**3), 1.e-12)
          phist = 1 + 5.0*(zstar - rough%disp)*rL

          where ((zetar(:,iter) .gt. 1.e-6).and.(.not.( canopy%vlaiw .LT. 0.01 .OR.                                          &
               rough%hruff .LT. rough%z0soilsn )))! stable conditions

             csw = min(0.3*((log((veg%hc-rough%disp)/rough%z0m) + phist*psihat - &
                  psim(zetar(:,iter)*(veg%hc-rough%disp)/(rough%zref_tq-rough%disp))+ &
                  psim(zetar(:,iter)*rough%z0m/(rough%zref_tq-rough%disp)))/0.4)**2/2., 3.0)* c%csw

             rough%term2  = EXP( 2 * CSW * canopy%rghlai *                          &
                  ( 1 - rough%disp / rough%hruff ) )
             rough%term3  = C%A33**2 * C%CTL * 2 * CSW * canopy%rghlai
             rough%term5  = MAX( ( 2. / 3. ) * rough%hruff / rough%disp, 1.0 )
             rough%term6 =  EXP( 3. * rough%coexp * ( rough%disp / rough%hruff -1. ) )

             !! vh ! Haverd et al., Biogeosciences 10, 2011-2040, 2013
             rough%rt0us  = log(rough%disp / rough%z0soilsn) * &
                  EXP(2. * C%CSW * canopy%rghlai) * rough%disp &
                  / rough%hruff / (c%a33 ** 2 * c%ctl)

             ! vh ! Modify rt0us to be resistance between shear height = 0.1h and disp
             ! use this form when including addtional resistance from z0soil to 0.1hc (done in cable_canopy_vh)
             rough%rt0us  = log(rough%disp/(0.1 * rough%hruff)) * &
                  EXP(2. * C%CSW * canopy%rghlai) * rough%disp &
                  / rough%hruff / (c%a33 ** 2 * c%ctl)

             rough%rt1usa = rough%term5 * ( rough%term2 - 1.0 ) / rough%term3
          elsewhere
             csw = c%csw
          endwhere
       endwhere
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       ! rt0 = turbulent resistance from soil to canopy:
       ! use this when accounting for separate resistance from z0soil to 0.1hc
       WHERE (canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
          rt0 = (rough%rt0us+rt0bus) / canopy%us
       elsewhere
          rt0 = (rough%rt0us) / canopy%us
       endwhere


       ! rt0 = (rough%rt0us) / canopy%us
       ! Aerodynamic resistance (sum 3 height integrals)/us
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
       rough%rt1 = max(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)
       WHERE (ssnow%snowd > 0.1)
          ssnow%wetfac = 1.0
       END WHERE
       ! change canopy%vlaiw requirement to 0.01 for conformity (BP may 2009)
       WHERE (canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
          ssnow%rtsoil = rt0
       ELSEWHERE
          ssnow%rtsoil = rt0 + rough%rt1
       END WHERE
       ssnow%rtsoil = max(25.,ssnow%rtsoil)
       !WHERE (ssnow%rtsoil>2.0_r_2*ortsoil .OR. ssnow%rtsoil<0.5_r_2*ortsoil)
       !   ssnow%rtsoil = MAX(25._r_2, 0.5_r_2*(ssnow%rtsoil + ortsoil))
       !END WHERE

       ! Vegetation boundary-layer conductance (mol/m2/s)
       ! prandt = kinematic viscosity/molecular diffusivity
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
       WHERE (canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
          gbvtop = real(air%cmolar*apol * air%visc / prandt / veg%dleaf *       &
               sqrt(canopy%us / MIN(rough%usuh, 0.2) * veg%dleaf / air%visc) * prandt**(0.3333333) / veg%shelrb, r_2)
          ! Forced convection boundary layer conductance (see Wang & Leuning 1998, AFM):
          !                                gbhu corrected by F.Cruz & A.Pitman on 13/03/07
          tmp1d = 0.5*rough%coexp+rad%extkb
          where (tmp1d < 40.0) ! MC: avoid floating underflow
             gbhu(:,1) = gbvtop*(1.0_r_2-EXP(real(-canopy%vlaiw*tmp1d,r_2))) / real(tmp1d,r_2)
          elsewhere
             gbhu(:,1) = gbvtop* 1.0_r_2 / real(tmp1d,r_2)
          endwhere
          ! MC include max because can go <0 in some weird cirumstances
          gbhu(:,2) = max(0.0_r_2, (2.0_r_2/real(rough%coexp,r_2))*gbvtop* &
               (1.0_r_2-EXP(real(-0.5*rough%coexp*canopy%vlaiw,r_2)))-gbhu(:,1) )
       ENDWHERE

       ! special for BIOS2, where input Ta is screen temperature
       !********************************************************************************************
       !      tstar = - (canopy%fh) / ( air%rho*capp*canopy%us)
       !      qstar = - (canopy%fe) / ( air%rho*air%rlam *canopy%us)
       !      zscrn = max(rough%z0m,1.8-rough%disp)
       !      denom = ( log(rough%zref/zscrn)- psis(zetar(:,iter)) + &
       !           psis(zetar(:,iter) * zscrn / rough%zref) ) /vonk

       !      where (zetar(:,iter) > 0.7)
       !         zeta2=zetar(:,iter) * zscrn / rough%zref
       !         denom =alpha1* ((zetar(:,iter)**beta1* &
       !              (1.0+gamma1*zetar(:,iter)**(1.0-beta1)))  &
       !              - (zeta2*beta1*(1.0+gamma1*zeta2**(1.0-beta1)))) /vonk
       !      endwhere

       ! Calculate temperature, humidity at ref height:
       !       DO i = 1, max_vegpatches ! over all tiles
       !          met%tk(i:mp:max_vegpatches)=canopy%tscrn(1:mp:max_vegpatches) + &
       !              tstar(1:mp:max_vegpatches) * denom(1:mp:max_vegpatches)
       !          met%qv(i:mp:max_vegpatches)=canopy%qscrn(1:mp:max_vegpatches) + &
       !              qstar(1:mp:max_vegpatches) * denom(1:mp:max_vegpatches)
       !       ENDDO
       !**********************************************************************************************

       WHERE (veg%meth > 0 .and. canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
          !      use the dispersion matrix (DM) to find the air temperature and specific humidity
          !      (Raupach, Finkele and Zhang 1997, pp 17)
          ! leaf boundary layer resistance for water
          rbw = air%cmolar/real(sum(gbhu+gbhf,2))    ! gbhf initially set to 1e-3
          rrbw = real(sum(gbhu+gbhf,2))/air%cmolar  ! MJT
          ! leaf stomatal resistance for water
          rsw = air%cmolar/real(sum(gswx,2))
          rrsw = real(sum(gswx,2))/air%cmolar ! MJT
          ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmah = (rt0+rough%rt1)*((1.+air%epsi)*rrsw +rrbw) &
               + air%epsi * (rt0*rough%rt1)*(rrbw*rrsw)
          ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmbh = (-air%rlam/capp)*(rt0*rough%rt1)*(rrbw*rrsw)
          ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132

          dmch = ((1.+air%epsi)*rrsw +rrbw)*rt0*rough%rt1* &
               (canopy%fhv + canopy%fhs)/(air%rho*capp)
          ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmae = (-air%epsi*capp/air%rlam)*(rt0*rough%rt1)*(rrbw*rrsw)
          ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmbe = (rt0+ssnow%wetfac*rough%rt1)*((1.+air%epsi)*rrsw +rrbw)+(rt0*rough%rt1)*(rrbw*rrsw)
          ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmce = ((1.+air%epsi)*rrsw +rrbw)*rt0*rough%rt1*(canopy%fev + real(canopy%fes))/ &
               (air%rho*air%rlam)

          where (veg%vlai>0.1)
             met%tvair = met%tk  + (dmbe*dmch-dmbh*dmce)/(dmah*dmbe-dmae*dmbh+1.0e-12)
          elsewhere
             met%tvair = met%tk
          endwhere


          where (abs(met%tvair-met%tk)>15.0)  ! in case of crazy deviations of incaopy air T: shouldn't need this!
             met%tvair = met%tk
          endwhere

          ! Within canopy specific humidity:
          where (veg%vlai>0.1)
             met%qvair = met%qv   + (dmah*dmce-dmae*dmch)/(dmah*dmbe-dmae*dmbh+1.0e-12)
          elsewhere
             met%qvair = met%qv
          endwhere
          met%qvair = max(0.0,met%qvair)
       END WHERE

       ! write(*,*) "LNF T", ktau, iter, canopy%fh, met%tvair- met%tk
       ! write(*,*) "LNF q", ktau, iter, canopy%fe, met%qvair- met%qv
       !   write(*,*) "u*", canopy%us, rt0, rough%rt1, ssnow%rtsoil
       !if (ktau.le.3001) then
       !    write(75,"(2i5, 100f16.6)") ktau, iter, met%ua, met%tvair-met%tk, canopy%us, zetar(:,iter), canopy%fh, &
       !    ssnow%rtsoil, canopy%fhs, csw
       !endif


       ! Saturated specific humidity in canopy:
       qstvair = qsatf((met%tvair-tfrz),met%pmb)
       met%qvair = min(qstvair,met%qvair)          ! avoid -ve dva
       ! Saturated vapour pressure deficit in canopy:
       met%dva = (qstvair - met%qvair) *  rmair/rmh2o * met%pmb * 100.
       ! 2 dim Within canopy air temperature in degrees C:
       tvair2 = real(SPREAD(met%tvair-tfrz, 2, mf),r_2)  ! N.B. tvair2 and tair2 ned to be equal because
       !the temperature difference (Ta-Tl) need to be
       !the same for both the sensible heat flux and non-isothermal radiation flux at the leaf
       tair2 = tvair2
       ! Set radiative temperature as within canopy air temp:
       met%tvrad = met%tvair
       ! call radiation here: longwave isothermal radiation absorption
       ! and radiation conductance depends on tvrad
       CALL radiation(ssnow, veg, air, met, rad, canopy)

       CALL define_air(met, air)
       dsatdk2 = SPREAD(air%dsatdk, 2, mf)

       hcx = 0.0_r_2             ! init sens heat iteration memory variable
       ecx = real(rad%rniso,r_2) ! init lat heat iteration memory variable
       rnx = real(rad%rniso,r_2) ! init net rad iteration memory variable
       rny = real(rad%rniso,r_2) ! init current estimate net rad
       hcy = 0.0_r_2             ! init current estimate lat heat
       ecy = rny - hcy           ! init current estimate lat heat
       abs_deltlf = 999.0
       deltlfy    = 999.0
       ! Initialise, over each gridpoint, sunlit and shaded leaves:
       DO k=1, mp
          DO kk=1, mf
             IF (rad%fvlai(k,kk) <= 1.0e-7) THEN
                abs_deltlf(k,kk) = 0.0
                hcx(k,kk)  = 0.0_r_2   ! intialise
                ecx(k,kk)  = 0.0_r_2   ! intialise
                anx(k,kk)  = 0.0_r_2   ! intialise
                rnx(k,kk)  = 0.0_r_2   ! intialise
                rny(k,kk)  = rnx(k,kk) ! store initial values
                hcy(k,kk)  = hcx(k,kk) ! store initial values
                ecy(k,kk)  = ecx(k,kk) ! store initial values
                rdy(k,kk)  = rdx(k,kk) ! store initial values
                an_y(k,kk) = anx(k,kk) ! store initial values
             END IF
          ENDDO
       ENDDO
       deltlfy   = abs_deltlf
       k         = 0
       Flag_fwet = .false.

       DO WHILE ((ANY(abs_deltlf > 0.1) )  .AND.  k < maxiter)
          k = k + 1
          ! Where vegetation and no convergence...

          WHERE (rad%fvlai > 1e-5 .and. abs_deltlf > 0.1 .or. Flag_fwet)
             ! Grashof number (Leuning et al, 1995) eq E4:
             gras = max(1.0e-6_r_2,1.595E8_r_2*ABS(tlfx-tair2)*(xdleaf2**3))
             ! See Appendix E in (Leuning et al, 1995):
             gbhf = real(rad%fvlai*SPREAD(air%cmolar,2,mf)*0.5*dheat,r_2)*(gras**0.25)/xdleaf2
             ! Conductance for heat:
             gh    = 2.0_r_2 * (gbhu + gbhf)
             ghwet = 2.0_r_2 * (gbhu + gbhf) ! NB changed wet leaf bl conductance to include forced component VH 16/10/09
             gh    = fwet*ghwet + (1.0_r_2-fwet)*gh
             ! Conductance for heat and longwave radiation:
             ghr  = real(rad%gradis,r_2) + gh
             temp =  xvcmxt3(tlfx+tfrz_r_2)
             !  Leuning 2002 (P C & E) equation for temperature response
             !  used for Vcmax for C3 plants:
             vcmxt3 = (1.0_r_2-frac42)*vcmax2*real(rad%scalex,r_2) *temp

             temp=  xvcmxt4(tlfx)
             ! Temperature of Vcmax for C4 plants (Collatz et al 1989):
             vcmxt4 = frac42 * vcmax2 * real(rad%scalex,r_2) * temp
             temp= xejmxt3(tlfx+tfrz_r_2)
             !  Leuning 2002 (P C & E) equation for temperature response
             !  used for Jmax for C3 plants:
             ejmxt3 = (1.0_r_2-frac42) * ejmax2 * real(rad%scalex,r_2) * temp
             ! Difference between leaf temperature and reference temperature:
             tdiff  =  tlfx + tfrz_r_2 - trefk_r_2
             ! Michaelis menten constant of Rubisco for CO2:
             conkct = conkc0*EXP(real(ekc/(rgas*trefk),r_2) *(1.0_r_2-trefk_r_2/(tlfx+tfrz_r_2)))
             ! Michaelis menten constant of Rubisco for oxygen:
             conkot = conko0*EXP(real(eko/(rgas*trefk),r_2) *(1.0_r_2-trefk_r_2/(tlfx+tfrz_r_2)))
             ! "d_{3}" in Wang and Leuning, 1998, appendix E:
             cx(:,:,1) = conkct*(1.0_r_2+0.21_r_2/conkot)
             cx(:,:,2) = 2.0_r_2* gam0*(1.0_r_2+gam1*tdiff + gam2*tdiff*tdiff) !gamma*
             ! All equations below in appendix E in Wang and Leuning 1998 are
             ! for calculating anx, csx and gswx for Rubisco limited, RuBP limited,
             ! sink limited
             vx3(:,:,1) = vcmxt3
             vx4(:,:,1) = vcmxt4
             temp = rad%qcan(:,:,1)*jtomol*(1.0_r_2-frac42)
             vx3(:,:,2) = ej3x(temp,ejmxt3)
             temp = frac42*rad%qcan(:,:,1)*jtomol
             vx4(:,:,2) = ej4x(temp,vcmxt4)
             rdx = cfrd3*vcmxt3+cfrd4*vcmxt4 !*fwsoil2
             xleuning = (1.0_r_2-frac42)*a1c3/(1.0_r_2+dsx/d0c3) +frac42*a1c4/(1.0_r_2+dsx/d0c4)
             xleuning = xleuning * fwsoil2 / (csx-co2cp3)
             ! Rubisco limited:
             coef2(:,:,1) = gswmin*fwsoil2/rgswc+xleuning *(vx3(:,:,1)-(rdx-vx4(:,:,1)))
             coef1(:,:,1) = (1.0_r_2-csx*xleuning) *(vx3(:,:,1)+vx4(:,:,1)-rdx)     &
                  +(gswmin*fwsoil2/rgswc)*(cx(:,:,1)-csx) -xleuning*(vx3(:,:,1)*cx(:,:,2)*0.5_r_2 &
                  +cx(:,:,1)*(rdx-vx4(:,:,1)))
             coef0(:,:,1) = -(1.0_r_2-csx*xleuning) *(vx3(:,:,1)*cx(:,:,2)*0.5_r_2  &
                  +cx(:,:,1)*(rdx-vx4(:,:,1))) -(gswmin*fwsoil2/rgswc)*cx(:,:,1)*csx
             ! Discriminant in quadratic in eq. E7 Wang and Leuning, 1998
             delcx(:,:,1) = coef1(:,:,1)**2 -4.0_r_2*coef0(:,:,1)*coef2(:,:,1)
             where (coef2(:,:,1)>1.e-8_r_2)
                ci(:,:,1) = (-coef1(:,:,1)+SQRT(MAX(0.0_r_2,delcx(:,:,1)))) /(2.0_r_2*coef2(:,:,1))
             elsewhere
                ci(:,:,1) = cx(:,:,2)*0.5_r_2
             endwhere
             ancj(:,:,1) = vx3(:,:,1)*(ci(:,:,1)-cx(:,:,2)*0.5_r_2) &
                  / (ci(:,:,1) + cx(:,:,1)) + vx4(:,:,1) - rdx
             ! RuBP limited:
             coef2(:,:,2) = real(gswmin*fwsoil2/rgswc,r_2)+xleuning *(vx3(:,:,2)-(rdx-vx4(:,:,2)))
             coef1(:,:,2) = (1.0_r_2-csx*xleuning) *(vx3(:,:,2)+vx4(:,:,2)-rdx)     &
                  +(gswmin*fwsoil2/rgswc)*(cx(:,:,2)-csx) -xleuning*(vx3(:,:,2)*cx(:,:,2)*0.5_r_2 &
                  +cx(:,:,2)*(rdx-vx4(:,:,2)))
             coef0(:,:,2) = -(1.0_r_2-csx*xleuning) *(vx3(:,:,2)*cx(:,:,2)*0.5_r_2  &
                  +cx(:,:,2)*(rdx-vx4(:,:,2))) -(gswmin*fwsoil2/rgswc)*cx(:,:,2)*csx
             delcx(:,:,2) = coef1(:,:,2)**2 -4.0_r_2*coef0(:,:,2)*coef2(:,:,2)
             where (coef2(:,:,2)>1.e-8_r_2)
                ci(:,:,2) = (-coef1(:,:,2)+SQRT(MAX(0.0_r_2,delcx(:,:,2)))) /(2.0_r_2*coef2(:,:,2))
                !ci(:,:,2) = MAX(0.0_r_2,ci(:,:,2))
             elsewhere
                ci(:,:,2) = cx(:,:,2)*0.5_r_2
             endwhere

             ancj(:,:,2) = vx3(:,:,2)*(ci(:,:,2)-cx(:,:,2)*0.5_r_2) &
                  /(ci(:,:,2)+cx(:,:,2)) +vx4(:,:,2)-rdx
             ! Sink limited:
             coef2(:,:,3) = xleuning
             coef1(:,:,3) = gswmin*fwsoil2/rgswc + xleuning * (rdx - 0.5_r_2*vcmxt3)  +  &
                  effc4 * vcmxt4 - xleuning * csx * effc4 *vcmxt4
             coef0(:,:,3) = -(gswmin*fwsoil2/rgswc)*csx *effc4*vcmxt4 + &
                  (rdx -0.5_r_2*vcmxt3)*gswmin*fwsoil2/rgswc
             delcx(:,:,3) = coef1(:,:,3)**2 -4.0_r_2*coef0(:,:,3)*coef2(:,:,3)

             where (coef2(:,:,3) > 1.e-8_r_2)
                ancj(:,:,3)  = (-coef1(:,:,3)+SQRT(MAX(0.0_r_2,delcx(:,:,3)))) &
                     /(2.0_r_2*coef2(:,:,3))
             elsewhere
                ancj(:,:,3)  = 0.0_r_2
             endwhere

             anx = MIN(ancj(:,:,1),ancj(:,:,2),ancj(:,:,3))
             csx  = ca2 - anx * (gbhu + gbhf) / rgbwc
             gswx = gswmin*fwsoil2 + MAX(0.0_r_2,rgswc*xleuning*anx)
             ! Recalculate conductance for water:
             gbw  = 1.075_r_2*(gbhu+gbhf)
             where (gswx > 1e-15_r_2)
                gw = 1.0_r_2/(1.0/gswx + 1.0/gbw)
             elsewhere
                gw = 0.0_r_2
             endwhere
             ! corrected vapour conductance for influence of wet part of leaf.
             ! Still need to reduce carbon conducatnce by factor of (1-fwet)
             gw = (1.0_r_2-fwet)*gw + fwet*gbw
             ! Modified psychrometric constant (Monteith and Unsworth, 1990)
             where (gw.gt.1e-15_r_2)
                psycst = real(SPREAD(air%psyc,2,mf),r_2) * ghr/gw
                ! Update canopy latent heat flux:
                ecx = (dsatdk2*real(rad%rniso,r_2)+real(capp*rmair,r_2)*da2*ghr) /(dsatdk2+psycst)
             elsewhere
                psycst = 0.0_r_2
                ecx = 0.0_r_2
             endwhere

             ! Store leaf temperature:
             tlfxx = tlfx
             ! Update canopy sensible heat flux:
             hcx = (real(rad%rniso,r_2)-ecx)*gh/ghr
             ! Update leaf temperature:
             tlfx=tair2+hcx/(real(capp*rmair,r_2)*gh)

             ! Update net radiation for canopy:
             rnx = real(rad%rniso,r_2)-real(capp*rmair*rad%gradis,r_2)*(tlfx-tair2)
             ! Update leaf surface vapour pressure deficit:
             ! dsx = ecx*100.0* SPREAD(met%pmb, 2, mf) /(gswx*rmh2o*SPREAD(air%rlam, 2, mf))
             dsx = da2 + dsatdk2 * (tlfx-tair2)
             ! Store change in leaf temperature between successive iterations:
             deltlf = tlfxx-tlfx
             abs_deltlf = ABS(deltlf)
          END WHERE



          ! Where leaf temp change b/w iterations is significant, and difference is
          ! smaller than the previous iteration, store results:
          !MC-Guess: update y-values after END DO, i.e. move 30 lines down
          WHERE (abs_deltlf > 0.1_r_2 .AND. abs_deltlf < ABS(deltlfy) )
             deltlfy = deltlf
             tlfy    = tlfx
             rny     = rnx
             hcy     = hcx
             ecy     = ecx
             rdy     = rdx
             an_y    = anx
          END WHERE

          rk = real(k,r_2)
          WHERE (abs_deltlf > 0.1_r_2)
             !        after four iteration, take the mean value of current and previous estimates
             !        as the next estimate of leaf temperature, to avoid oscillation
             tlfx = (0.5_r_2*(MAX(0._r_2,rk-5._r_2)/(rk-4.9999_r_2))) *tlfxx + &
                  (1.0_r_2- (0.5_r_2*(MAX(0._r_2,rk-5._r_2)/(rk-4.9999_r_2))))*tlfx
             !MC-Guess: update tlfy here?
             !          i.e. tlfy = tlfx
          END WHERE
          IF (k==1) THEN
             !        taken the first iterated estimates as the defaults
             tlfy = tlfx
             rny  = rnx
             hcy  = hcx
             ecy  = ecx
             rdy  = rdx
             an_y = anx
          END IF
       END DO  ! DO WHILE (ANY(abs_deltlf > 0.1)        .AND.  k < maxiter)

       tlfy    = tlfx
       rny     = rnx
       hcy     = hcx
       ecy     = ecx
       rdy     = rdx
       an_y    = anx
       canopy%fev = real(sum(ecy,2))
       canopy%fhv = real(sum(hcy,2))
       canopy%fnv = real(sum(rny,2))

       an_y = (an_y+rdy)*(1.0_r_2-fwet) - rdy    ! only allow gross photosynthesis on dry part of leaf

       ! evaulate canopy%fes, canopy%fhs for use in dispersion matrix calc
       rad%lwabv = (capp*rmair*real(tlfy(:,1) - tvair2(:,1))*rad%gradis(:,1) &
            +capp*rmair*real(tlfy(:,2) - tvair2(:,2))*rad%gradis(:,2)) ! non-isothermal emitted long-wave radiation
       WHERE (canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
          canopy%tv = (rad%lwabv / (2.0*(1.0-rad%transd)*sboltz*emleaf)+met%tvrad**4)**0.25
       ELSEWHERE ! sparse canopy
          canopy%tv = met%tvair
       END WHERE

       canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*emleaf* &
            sboltz*canopy%tv**4 - emsoil*sboltz* tss4

       ! Calculate radiative/skin temperature:
       rad%trad = ( (1.-rad%transd)*canopy%tv**4 &
            + rad%transd * ssnow%tss**4 )**0.25

       IF (cable_user%soil_struc=='default') THEN
          ! Penman-Monteith formula
          sss = air%dsatdk
          cc1 = sss/(sss+air%psyc )
          cc2 = air%psyc /(sss+air%psyc )

          ssnow%potev = cc1 * (canopy%fns - canopy%ga) + cc2 * air%rho  &
               * real(air%rlam) * (qsatf((met%tk-tfrz),met%pmb) - met%qv) / ssnow%rtsoil
          ! Soil latent heat:
          canopy%fes= real(ssnow%wetfac,r_2) * ssnow%potev
          WHERE (ssnow%snowd < 0.1_r_2 .AND. canopy%fes > 0.0)
             ! Reduce for wilting point limitation:
             canopy%fes= MIN( canopy%fes, &
                  MAX(0.0_r_2, (ssnow%wb(:,1)-real(soil%swilt,r_2))) * &
                  soil%zse(1)*1000.0*real(air%rlam)/dels)
             ! Reduce for soil ice limitation:
             canopy%fes = MIN(canopy%fes, (ssnow%wb(:,1)-ssnow%wbice(:,1)) &
                  * soil%zse(1) * 1000. * real(air%rlam) / dels)
          END WHERE
          ssnow%cls = 1.
          WHERE (ssnow%snowd >= 0.1)
             ssnow%cls = 1.1335
             canopy%fes= MIN(ssnow%wetfac * real(ssnow%potev),real(ssnow%snowd*air%rlam*ssnow%cls)/dels)
          END WHERE
          ! Calculate soil sensible heat:
          canopy%fhs = air%rho*capp*(ssnow%tss - met%tk) /ssnow%rtsoil
          ! Calculate ground heat flux:
          canopy%ga = canopy%fns - canopy%fhs - real(canopy%fes)*ssnow%cls

       ELSEIF (cable_user%soil_struc=='sli') THEN

          CALL sli_main(ktau,dels,veg,soil,ssnow,met,canopy,air,rad,1)


       END IF

       ! Calculate total latent heat:
       canopy%fe = canopy%fev + real(canopy%fes)
       ! Calculate total sensible heat:
       canopy%fh = real(canopy%fhv) + canopy%fhs
       ! write(*,*) "end do iter" , iter, canopy%fe, canopy%fh
    END DO      ! do iter = 1, niter

    where ((fwet*ecy) > Ecansto)  ! move this inside Tleaf loop?
       fwet = Ecansto/ecy
       Flag_fwet = .true.
    elsewhere (ecy < 0.0_r_2) ! dew formation
       fwet = 1.0_r_2
       Flag_fwet = .true.
    elsewhere
       Flag_fwet = .false.
    endwhere

    where (sum(ecy,2) > 0.0_r_2)  ! evaporation
       canopy%fevw = min(real(sum(ecy*fwet*gbw/max(gw,1.e-6_r_2),2)), &
            max(0.0,canopy%cansto)*air%rlam/dels)
       canopy%fwet = REAL(canopy%fevw/sum(ecy,2))
    elsewhere ! condensation
       canopy%fevw = real(sum(ecy*fwet*gbw/max(gw,1.e-6_r_2),2))
       canopy%fwet = 1.0
    endwhere

    where (sum(gswx,2)>1e-15_r_2)
       canopy%fevc = canopy%fev - canopy%fevw
    elsewhere
       canopy%fevc = 0.0_r_2
    endwhere

    where (canopy%fevc<0.0_r_2)  ! negative values of fevc due to precision
       canopy%fevc = 0.0_r_2
       canopy%fevw = canopy%fev
    endwhere
    canopy%fhvw = canopy%fhv*canopy%fwet



    if (1==0) then  ! comment out calc of variables at screen height vh 15/07/09:
       ! this was causing model to crash. Update with CABLE 2.0 formulations?
       ! screen temp., windspeed and relative humidity at 1.8m
       tstar = - canopy%fh / ( air%rho*capp*canopy%us)
       qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us)
       zscrn = max(rough%z0m,1.8-rough%disp)
       !    denom = ( log(rough%zref/zscrn)- psim(zetar(:,iterplus)) + &
       !         psim(zetar(:,iterplus) * zscrn / rough%zref) ) /vonk
       denom = ( log(rough%zref_tq/zscrn)- psis(zetar(:,iterplus)) + &
            psis(zetar(:,iterplus) * zscrn / rough%zref_tq) ) /vonk

       !%% change by Ashok Luhar
       where (zetar(:,iterplus) > 0.7)
          zeta2=zetar(:,iterplus) * zscrn / rough%zref_tq
          denom =alpha1* ((zetar(:,iterplus)**beta1* &
               (1.0+gamma1*zetar(:,iterplus)**(1.0-beta1)))  &
               - (zeta2*beta1*(1.0+gamma1*zeta2**(1.0-beta1)))) /vonk
       endwhere
       !%%

       ! Calculate screen temperature:
       canopy%tscrn = met%tk-tfrz - tstar * denom
       rsts = qsatf(canopy%tscrn, met%pmb)
       qtgnet = rsts * ssnow%wetfac - met%qv
       canopy%cduv = canopy%us * canopy%us / (MAX(met%ua,umin))**2 ! EK jun08
       !    canopy%cduv = canopy%us * canopy%us / max(met%ua,umin)
       WHERE (qtgnet > 0.0)
          qsurf = rsts * ssnow%wetfac
       ELSEWHERE
          qsurf = 0.1*rsts*ssnow%wetfac + 0.9*met%qv
       END WHERE
       canopy%qscrn = qsurf + qstar * denom
       canopy%uscrn = max(0.0, max(met%ua,umin) - canopy%us * denom )    ! at present incorrect
    endif              ! comment out calc of variables at screen height vh 15/07/09

    canopy%frday = 12.0 * real(sum(rdy, 2))
    canopy%fpn   = -12.0 * real(sum(an_y, 2))
    ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
    canopy%dewmm = - REAL((min(0.0,canopy%fevw) + min(0.0_r_2,canopy%fevc))/air%rlam) *  &
         dels * 1.0e3 / rhow
    ! Add dewfall to canopy water storage:
    canopy%cansto = canopy%cansto + canopy%dewmm
    ! Modify canopy water storage for evaporation:
    canopy%cansto =  max(canopy%cansto-max(0.0,canopy%fevw)/air%rlam*1.0e3*dels/rhow, 0.0)
    ! Calculate canopy water storage excess:
    canopy%spill= real(max(0.0_r_2,min(0.2_r_2*canopy%cansto,max(0.0_r_2,canopy%cansto-cansat))))
    ! Move excess canopy water to throughfall:
    canopy%through = canopy%through + canopy%spill
    ! Initialise 'throughfall to soil' as 'throughfall from canopy'; snow may absorb
    canopy%precis = canopy%through
    ! Update canopy storage term:
    canopy%cansto = canopy%cansto - canopy%spill

    ! Calculate the total change in canopy water store (mm/dels):
    !canopy%delwc = canopy%cansto-oldcansto
    ! calculate dgdtg, derivative of ga
    ssnow%dfn_dtg = (-1.)*4.*emsoil*sboltz*tss4/ssnow%tss            ! d(canopy%fns)/d(ssnow%tgg)
    ssnow%dfh_dtg = air%rho*capp/ssnow%rtsoil      ! d(canopy%fhs)/d(ssnow%tgg)
    ssnow%dfe_ddq = ssnow%wetfac*air%rho/ssnow%rtsoil*air%rlam  ! d(canopy%fes)/d(dq)
    ssnow%ddq_dtg = (rmh2o/rmair)/met%pmb *tetena*tetenb*tetenc &
         /((tetenc+ssnow%tss-tfrz)**2)*exp(tetenb*(ssnow%tss-tfrz)/(tetenc+ssnow%tss-tfrz))
    ssnow%cls = 1.0
    WHERE (ssnow%snowd >= 0.1) ssnow%cls = 1.1335
    canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg &
         - ssnow%cls * ssnow%dfe_ddq * ssnow%ddq_dtg


    ! owetfac will need to be outside define_canopy
    ! because the UM driver will call define_canopy twice
    ssnow%owetfac = ssnow%wetfac



    ! Calculate the total change in canopy water store (mm/dels):
    canopy%delwc = canopy%cansto-oldcansto
    !******************************************************************************************
    canopy%gw     = gw   ! edit vh 6/7/09
    canopy%ancj   = ancj*(-12.0) ! edit vh 7/7/09
    canopy%gswx   = real(gswx) ! edit vh 7/7/09
    canopy%tlfy   = tlfy ! edit vh 7/7/09
    canopy%ecy    = ecy ! edit vh 7/7/09
    canopy%ecx    = ecx ! edit vh 7/7/09
    canopy%ci     = ci ! edit vh 7/7/09
    !******************************************************************************************

  END SUBROUTINE define_canopy_vh

  !--------------------------------------------------------------------------
  ELEMENTAL FUNCTION qsatf_r_1(tair,pmb)
    ! MRR, 1987
    ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
    ! HUMIDITY (KG/KG) FROM TETEN FORMULA
    REAL, INTENT(IN) :: tair ! air temperature (C)
    REAL, INTENT(IN) :: pmb  ! pressure PMB (mb)
    REAL           :: qsatf_r_1    ! result; sat sp humidity

    qsatf_r_1 = (rmh2o/rmair) * (tetena*EXP(tetenb*tair/(tetenc+tair))) / pmb

  END FUNCTION qsatf_r_1

  ELEMENTAL FUNCTION qsatf_r_2(tair,pmb)
    ! MRR, 1987
    ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
    ! HUMIDITY (KG/KG) FROM TETEN FORMULA
    REAL(r_2), INTENT(IN) :: tair ! air temperature (C)
    REAL(r_2), INTENT(IN) :: pmb  ! pressure PMB (mb)
    REAL(r_2)           :: qsatf_r_2    ! result; sat sp humidity

    qsatf_r_2 = rmh2o_r_2/rmair_r_2 * (tetena_r_2*EXP(tetenb_r_2*tair/(tetenc_r_2+tair))) / pmb

  END FUNCTION qsatf_r_2

  !---------------------------------------------------------
  ELEMENTAL FUNCTION ej3x_r_1(parx,x)

    REAL, INTENT(IN)     :: parx
    REAL, INTENT(IN)     :: x
    REAL                 :: ej3x_r_1

    ej3x_r_1 = max(0.0, &
         0.25*((alpha3*parx+x-sqrt((alpha3*parx+x)**2 - &
         4.0*convx3*alpha3*parx*x)) /(2.0*convx3)) )

  END FUNCTION ej3x_r_1

  ELEMENTAL FUNCTION ej3x_r_2(parx,x)

    REAL(r_2), INTENT(IN)     :: parx
    REAL(r_2), INTENT(IN)     :: x
    REAL(r_2)                 :: ej3x_r_2

    ej3x_r_2 = max(0.0_r_2, &
         0.25_r_2*((alpha3_r_2*parx+x-sqrt((alpha3_r_2*parx+x)**2 - &
         4.0_r_2*convx3_r_2*alpha3_r_2*parx*x)) /(2.0_r_2*convx3_r_2)) )

  END FUNCTION ej3x_r_2

  !---------------------------------------------------------
  ELEMENTAL FUNCTION ej4x_r_1(parx,x)

    REAL, INTENT(IN)     :: parx
    REAL, INTENT(IN)     :: x
    REAL                 :: ej4x_r_1

    ej4x_r_1 = max(0.0, &
         (alpha4*parx+x-sqrt((alpha4*parx+x)**2 - &
         4.0*convx4*alpha4*parx*x))/(2.0*convx4))

  END FUNCTION ej4x_r_1

  ELEMENTAL FUNCTION ej4x_r_2(parx,x)

    REAL(r_2), INTENT(IN)     :: parx
    REAL(r_2), INTENT(IN)     :: x
    REAL(r_2)                 :: ej4x_r_2

    ej4x_r_2 = max(0.0_r_2, &
         (alpha4_r_2*parx+x-sqrt((alpha4_r_2*parx+x)**2 - &
         4.0_r_2*convx4_r_2*alpha4_r_2*parx*x))/(2.0_r_2*convx4_r_2))

  END FUNCTION ej4x_r_2

  !---------------------------------------------------------
  ! Explicit array dimensions as temporary work around for NEC inlining problem
  FUNCTION xvcmxt4_r_1(x)

    REAL, PARAMETER      :: q10c4 = 2.0
    REAL, DIMENSION(mp,mf), INTENT(IN)     :: x
    REAL, DIMENSION(mp,mf)                 :: xvcmxt4_r_1

    xvcmxt4_r_1 = q10c4 ** (0.1 * x - 2.5) / &
         ((1.0 + exp(0.3 * (13.0 - x))) * (1.0 + exp(0.3 * (x - 36.0))))

  END FUNCTION xvcmxt4_r_1

  FUNCTION xvcmxt4_r_2(x)

    REAL(r_2), PARAMETER      :: q10c4 = 2.0
    REAL(r_2), DIMENSION(mp,mf), INTENT(IN)     :: x
    REAL(r_2), DIMENSION(mp,mf)                 :: xvcmxt4_r_2

    xvcmxt4_r_2 = q10c4 ** (0.1_r_2 * x - 2.5_r_2) / &
         ((1.0_r_2 + exp(0.3_r_2 * (13.0_r_2 - x))) * (1.0_r_2 + exp(0.3_r_2 * (x - 36.0_r_2))))

  END FUNCTION xvcmxt4_r_2

  !---------------------------------------------------------
  FUNCTION xvcmxt3_r_1(x)
    !  leuning 2002 (p c & e) equation for temperature response
    !  used for vcmax for c3 plants
    REAL, DIMENSION(mp,mf), INTENT(IN)     :: x
    REAL, DIMENSION(mp,mf)         :: xvcnum
    REAL, DIMENSION(mp,mf)         :: xvcden
    REAL, DIMENSION(mp,mf)         :: xvcmxt3_r_1

    xvcnum  = xvccoef*exp((ehavc/(rgas*trefk))*(1.-trefk/x))
    xvcden  = 1.0+exp((entropvc*x-ehdvc)/(rgas*x))
    xvcmxt3_r_1 = max(0.0,xvcnum/xvcden)

  END FUNCTION xvcmxt3_r_1

  FUNCTION xvcmxt3_r_2(x)
    !  leuning 2002 (p c & e) equation for temperature response
    !  used for vcmax for c3 plants
    REAL(r_2), DIMENSION(mp,mf), INTENT(IN)     :: x
    REAL(r_2), DIMENSION(mp,mf)         :: xvcnum
    REAL(r_2), DIMENSION(mp,mf)         :: xvcden
    REAL(r_2), DIMENSION(mp,mf)         :: xvcmxt3_r_2

    xvcnum  = xvccoef_r_2*exp((ehavc_r_2/(rgas_r_2*trefk_r_2))*(1.0_r_2-trefk_r_2/x))
    xvcden  = 1.0_r_2+exp((entropvc_r_2*x-ehdvc_r_2)/(rgas_r_2*x))
    xvcmxt3_r_2 = max(0.0_r_2,xvcnum/xvcden)

  END FUNCTION xvcmxt3_r_2

  !---------------------------------------------------------
  FUNCTION xejmxt3_r_1(x)
    !  leuning 2002 (p c & e) equation for temperature response
    !  used for jmax for c3 plants
    REAL, DIMENSION(mp,mf), INTENT(IN)     :: x
    REAL, DIMENSION(mp,mf)         :: xjxnum
    REAL, DIMENSION(mp,mf)         :: xjxden
    REAL, DIMENSION(mp,mf)         :: xejmxt3_r_1

    xjxnum  = xjxcoef*exp((ehajx/(rgas*trefk))*(1.-trefk/x))
    xjxden  = 1.0+exp((entropjx*x-ehdjx)/(rgas*x))
    xejmxt3_r_1 = max(0.0, xjxnum/xjxden)

  END FUNCTION xejmxt3_r_1

  FUNCTION xejmxt3_r_2(x)
    !  leuning 2002 (p c & e) equation for temperature response
    !  used for jmax for c3 plants
    REAL(r_2), DIMENSION(mp,mf), INTENT(IN)     :: x
    REAL(r_2), DIMENSION(mp,mf)         :: xjxnum
    REAL(r_2), DIMENSION(mp,mf)         :: xjxden
    REAL(r_2), DIMENSION(mp,mf)         :: xejmxt3_r_2

    xjxnum  = xjxcoef_r_2*exp((ehajx_r_2/(rgas_r_2*trefk_r_2))*(1.0_r_2-trefk_r_2/x))
    xjxden  = 1.0_r_2+exp((entropjx_r_2*x-ehdjx_r_2)/(rgas_r_2*x))
    xejmxt3_r_2 = max(0.0_r_2, xjxnum/xjxden)

  END FUNCTION xejmxt3_r_2

  !---------------------------------------------------------
  ELEMENTAL FUNCTION psim_r_1(zeta)
    ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
    ! computes integrated stability function psim(z/l) (z/l=zeta)
    ! for momentum, using the businger-dyer form for unstable cases
    ! and the webb form for stable cases. see paulson (1970).
    USE math_constants

    REAL, INTENT(IN)     :: zeta
    REAL                 :: x
    REAL, PARAMETER      :: gu = 16.0
    REAL, PARAMETER      :: gs = 5.0
    REAL                 :: z
    REAL                 :: stable
    REAL                 :: unstable
    REAL                 :: psim_r_1

    !      x = (1.0 + gu*abs(zeta))**0.25
    !      r = merge(log((1.0+x*x)*(1.0+x)**2/8.0) - 2.0*atan(x) &
    !           + pi*0.5, -gs*zeta, zeta < 0.0)
    z        = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable
    stable   = -gs*zeta
    x        = (1.0 + gu*abs(zeta))**0.25
    unstable = log((1.0+x*x)*(1.0+x)**2/8.0) - 2.0*atan(x) + pi*0.5
    psim_r_1 = z*stable + (1.0-z)*unstable

  END FUNCTION psim_r_1

  ELEMENTAL FUNCTION psim_r_2(zeta)
    ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
    ! computes integrated stability function psim(z/l) (z/l=zeta)
    ! for momentum, using the businger-dyer form for unstable cases
    ! and the webb form for stable cases. see paulson (1970).
    USE math_constants

    REAL(r_2), INTENT(IN)     :: zeta
    REAL(r_2)                 :: x
    REAL(r_2), PARAMETER      :: gu = 16.0
    REAL(r_2), PARAMETER      :: gs = 5.0
    REAL(r_2)                 :: z
    REAL(r_2)                 :: stable
    REAL(r_2)                 :: unstable
    REAL(r_2)                 :: psim_r_2

    !      x = (1.0 + gu*abs(zeta))**0.25
    !      r = merge(log((1.0+x*x)*(1.0+x)**2/8.0) - 2.0*atan(x) &
    !           + pi*0.5, -gs*zeta, zeta < 0.0)
    z        = 0.5_r_2 + sign(0.5_r_2,zeta)    ! z=1 in stable, 0 in unstable
    stable   = -gs*zeta
    x        = (1.0_r_2 + gu*abs(zeta))**0.25_r_2
    unstable = log((1.0_r_2+x*x)*(1.0_r_2+x)**2/8.0_r_2) - 2.0_r_2*atan(x) + pi_r_2*0.5_r_2
    psim_r_2     = z*stable + (1.0_r_2-z)*unstable

  END FUNCTION psim_r_2

  !---------------------------------------------------------
  ELEMENTAL FUNCTION psis_r_1(zeta)
    ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
    ! computes integrated stability function psis(z/l) (z/l=zeta)
    ! for scalars, using the businger-dyer form for unstable cases
    ! and the webb form for stable cases. see paulson (1970).
    REAL, INTENT(IN)     :: zeta
    REAL, PARAMETER      :: gu = 16.0
    REAL, PARAMETER      :: gs = 5.0
    REAL                 :: z
    REAL                 :: y
    REAL                 :: stable
    REAL                 :: unstable
    REAL                 :: psis_r_1

    !      r = merge(2.0 * log((1.0 + sqrt(1.0 + gu * abs(zeta))) * 0.5), &
    !           - gs * zeta, zeta < 0.0)
    z        = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable
    stable   = -gs*zeta
    y        = (1.0 + gu*abs(zeta))**0.5
    unstable = 2.0 * log((1.0+y)*0.5)
    psis_r_1 = z*stable + (1.0-z)*unstable

  END FUNCTION psis_r_1

  ELEMENTAL FUNCTION psis_r_2(zeta)
    ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
    ! computes integrated stability function psis(z/l) (z/l=zeta)
    ! for scalars, using the businger-dyer form for unstable cases
    ! and the webb form for stable cases. see paulson (1970).
    REAL(r_2), INTENT(IN)     :: zeta
    REAL(r_2), PARAMETER      :: gu = 16.0_r_2
    REAL(r_2), PARAMETER      :: gs = 5.0_r_2
    REAL(r_2)                 :: z
    REAL(r_2)                 :: y
    REAL(r_2)                 :: stable
    REAL(r_2)                 :: unstable
    REAL(r_2)                 :: psis_r_2

    !      r = merge(2.0 * log((1.0 + sqrt(1.0 + gu * abs(zeta))) * 0.5), &
    !           - gs * zeta, zeta < 0.0)
    z        = 0.5_r_2 + sign(0.5_r_2,zeta)    ! z=1 in stable, 0 in unstable
    stable   = -gs*zeta
    y        = (1.0_r_2 + gu*abs(zeta))**0.5_r_2
    unstable = 2.0_r_2 * log((1.0_r_2+y)*0.5_r_2)
    psis_r_2     = z*stable + (1.0_r_2-z)*unstable

  END FUNCTION psis_r_2

  !---------------------------------------------------------
  ELEMENTAL FUNCTION rplant_r_1(rpconst, rpcoef, tair)

    REAL, INTENT(IN)     :: rpconst
    REAL, INTENT(IN)     :: rpcoef
    REAL, INTENT(IN)     :: tair
    REAL                 :: rplant_r_1

    rplant_r_1 = rpconst * exp(rpcoef * tair)

  END FUNCTION rplant_r_1

  ELEMENTAL FUNCTION rplant_r_2(rpconst, rpcoef, tair)

    REAL(r_2), INTENT(IN)     :: rpconst
    REAL(r_2), INTENT(IN)     :: rpcoef
    REAL(r_2), INTENT(IN)     :: tair
    REAL(r_2)                 :: rplant_r_2

    rplant_r_2 = rpconst * exp(rpcoef * tair)

  END FUNCTION rplant_r_2

  !---------------------------------------------------------
  ELEMENTAL FUNCTION rsoil_r_1(rsconst, avgwrs, avgtrs)

    REAL, INTENT(IN)     :: rsconst
    REAL, INTENT(IN)     :: avgwrs
    REAL, INTENT(IN)     :: avgtrs
    REAL                 :: rsoil_r_1

    rsoil_r_1 = rsconst * min(1.0, max(0.0, min( &
         -0.0178+0.2883*avgwrs+5.0176*avgwrs*avgwrs-4.5128*avgwrs*avgwrs*avgwrs, &
         0.3320+22.6726*exp(-5.8184*avgwrs)))) &
         * min(1.0, max(0.0, min( 0.0104*(avgtrs**1.3053), 5.5956-0.1189*avgtrs)))

  END FUNCTION rsoil_r_1

  ELEMENTAL FUNCTION rsoil_r_2(rsconst, avgwrs, avgtrs)

    REAL(r_2), INTENT(IN)     :: rsconst
    REAL(r_2), INTENT(IN)     :: avgwrs
    REAL(r_2), INTENT(IN)     :: avgtrs
    REAL(r_2)                 :: rsoil_r_2

    rsoil_r_2 = rsconst * min(1.0_r_2, max(0.0_r_2, min( &
         -0.0178_r_2+0.2883_r_2*avgwrs+5.0176_r_2*avgwrs*avgwrs-4.5128_r_2*avgwrs*avgwrs*avgwrs, &
         0.3320_r_2+22.6726_r_2*exp(-5.8184_r_2*avgwrs)))) &
         * min(1.0_r_2, max(0.0_r_2, min( 0.0104_r_2*(avgtrs**1.3053_r_2), 5.5956_r_2-0.1189_r_2*avgtrs)))

  END FUNCTION rsoil_r_2

END MODULE canopy_vh_module

