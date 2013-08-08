!
! Original code by P.J. Ross 2001-2007  See: Ross, P.J. (2003)
!     Modeling soil water and solute transport - fast, simplified numerical solutions. Agron. J. 95:1352-1361.
!
! Modified by V. Haverd 2008: matrix expanded by factor of 2 to allow solution of coupled heat/moisture fluxes
! Explicit calculation of heat/moisture fluxes at surface
! Switchable option for litter
! Isotope subroutine (V. Haverd and M. Cuntz, September 2008)
!
! Frozen Soil (V. Haverd, M. Cuntz Sep 2010 - Jan 2011)
! Pond lumped with top soil layer (V. Haverd Jan 2011)
! Snow included in layers -2:0 (V. Haverd Jan 2011)
! Convergence test for dthetaldT (V.Haverd Feb 2011)
! Include heat advection by liquid water flux (V.Haverd Feb 2011)

! Code cleaned and isotope routine updated (V. Haverd, M. Cuntz Aug-Sep 2012)
! lots of changes, especially relating to ponding and freezing together
! Outstanding problems with surface runoff and advection (set h0max large and turn off advection for now)
!
! Sep 17 2012 (V.Haverd)
! remove separate solution for maxpond (h0 > h0max)
! Instead, convert excess pond to runoff at t=tfin, and correct energy stores for loss of pond
! Advection and surface runoff now functioning
!
! Dec 31 2012 (V.Haverd)
! bug fixes to improve energy conservation associated with change in ice,
! especially correction at time of pond disappearance (removal of negative pond)
!
! Jan 5 2013 (V.Haverd)
! Remove iteration over frozen soil within "(while iok==0)" loop
! For frozen soil, at the time of updating variables, evaluate new T and ice content,
! consistent with J0 + deltaJ, where deltaJ is the change in energy obtained from the matrix solution (=LHS_h*dt)
!
! Jan 7-8 2013 (V.Haverd)
! Calls to soil_snow albedo and density
! Improve iteration convergence at time of updating new T and ice content in frozen soil
! Bug fix relating to moisture conservation in snow
!
! Jan 14-15 2013 (V. Haverd)
! Modified definition of Sliq for frozen soil (in hyofS): now defined relative to (theta_sat - theta_r - theta_ice)
! Allow for snow-pack initalisation in the absence of pond
!
! Feb 7 2013 (V.Haverd)
! Revised formulation of dphidT for frozen soil (see hyofS): fixes a few remaining -ve T spikes in frozen soil
!
! Feb 9 2013 (V.Haverd)
! Moved snow layer adjustments to subroutine
!
! Feb 24 2013 (V. Haverd)
! Snow structure expanded to accommodate 3 layers
!
! March 30 (V. Haverd)

MODULE sli_solve

  USE cable_def_types_mod, ONLY: r_2, i_d
  USE sli_numbers,         ONLY: &
       experiment, &
       zero, one, two, half, thousand, e3, e5, &  ! numbers
       Tzero, rlambda, lambdaf, lambdas, Dva, rhocp, rhow, gf, hmin, & ! parameters
       csice , cswat, &
       snmin, nsnow_max , fsnowliq_max, &
       params, vars_aquifer, vars_met, vars, vars_snow, solve_type, & ! types
       dSfac, h0min, Smax, h0max, dSmax, dSmaxr, dtmax, dSmax, dSmaxr, & ! numerical limits
       dtmax, dtmin, dTsoilmax, dTLmax, nsteps_ice_max, tol_dthetaldT, &
       hbot, botbc, count_sparse, count_hyofS, &
       MW, RGAS ! boundary condition

  USE sli_utils,           ONLY: &
       x, Sofh, hyofh, hyofS, litter_props, massman_sparse, tri, &
       aquifer_props, Tfrz, thetalmax, Tthetalmax, dthetalmaxdTh, &
       getfluxes_vp, getheatfluxes, flux, sol, &
       csat, slope_csat, potential_evap, tri, setsol, zerovars, &
       esat_ice, slope_esat_ice, Tfrozen, rtbis_Tfrozen, GTFrozen, &
       JSoilLayer

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: solve ! solution routine

  INTEGER(i_d), DIMENSION(:), ALLOCATABLE :: nless, n_noconverge ! global counters

  ! Definitions of public entities and private parameters (see above for default
  ! values):
  ! botbc    - bottom boundary condn for water; "constant head", "free drainage",
  !    "seepage", or "zero flux". Constant head means that matric head h
  !    is specified. Free drainage means zero gradient of matric head,
  !    i.e. unit hydraulic gradient. Seepage means zero flux when the
  !    matric head is below zero and an upper limit of zero for the head.
  ! h0max    - max pond depth allowed before runoff.
  ! hbot    - matric head at bottom of profile when botbc set to "constant head".
  ! dSmax    - max change in S (the "effective saturation") of any unsaturated
  !    layer to aim for each time step; controls time step size.
  ! dSmaxr   - maximum negative relative change in S each time step. This
  !    parameter helps avoid very small or negative S.
  ! dtmax    - max time step allowed.
  ! dsmmax   - max solute change per time step (see dSmax); user should set this
  !    according to solute units used. Units for different solutes can be
  !    scaled by the user (e.g. to an expected max of around 1.0).
  ! dSfac    - a change in S of up to dSfac*dSmax is accepted.
  ! dpmaxr   - relative change in matric flux potential (MFP) phi that is
  !    accepted for convergence when finding head h at soil interfaces.
  ! h0min    - min (negative) value for surface pond when it empties.
  ! Smax    - max value for layer saturation to allow some overshoot.
  ! solve    - sub to call to solve RE
  !

CONTAINS

  !**********************************************************************************************************************

  SUBROUTINE solve(ts, tfin, irec, mp, qprec,qprec_snow, n, dx, h0, S,thetai, Jsensible, Tsoil, evap, evap_pot, runoff, &
       infil, drainage, discharge, qh, nsteps, vmet, vlit,vsnow, csoil, kth, phi, T0, rh0, Tsurface, rhsurface, Hcum, lEcum, &
       Gcum, Qadvcum, Jcol_sensible, Jcol_latent_S, Jcol_latent_T, deltaice_cum_T, deltaice_cum_S, dxL, zdelta, &
       SL, TL, plit, par, qex, wex, heads, FS,  &
       ciso, cisoice, cisos, ciso0, cisoL, cprec, cali, qali, qiso_in, qiso_out, qiso_evap_cum, qiso_trans_cum, qiso_liq_adv, &
       qiso_vap_adv, qiso_liq_diff, qiso_vap_diff, qvsig, qlsig, qvTsig, qvh, deltaTa, lE_old, &
       dolitter, doisotopologue, dosepts, docondition, doadvection)


    IMPLICIT NONE

    REAL(r_2),                             INTENT(IN)              :: ts, tfin
    INTEGER(i_d),                          INTENT(IN)              :: irec, mp
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN)              :: qprec
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN)              :: qprec_snow
    INTEGER(i_d),                          INTENT(IN)              :: n
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(IN)              :: dx
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: h0
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT)           :: S
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT)             :: thetai
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT)             :: Jsensible
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT)           :: Tsoil
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)           :: evap, evap_pot, runoff, infil
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)             :: drainage, discharge
    REAL(r_2),      DIMENSION(1:mp,0:n),   INTENT(OUT)             :: qh
    INTEGER(i_d),   DIMENSION(1:mp),       INTENT(OUT)           :: nsteps
    TYPE(vars_met), DIMENSION(1:mp),       INTENT(INOUT)           :: vmet
    TYPE(vars),     DIMENSION(1:mp),       INTENT(INOUT)           :: vlit
    TYPE(vars_snow),     DIMENSION(1:mp),       INTENT(INOUT)           :: vsnow
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: T0, rh0, Tsurface, rhsurface
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)             :: Hcum, lEcum, Gcum,Qadvcum, deltaice_cum_T, deltaice_cum_S
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)             :: Jcol_sensible, Jcol_latent_S, Jcol_latent_T
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT)             :: csoil, kth
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT)             :: phi
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN)              :: dxL
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: zdelta, SL, Tl
    TYPE(params),   DIMENSION(1:mp),       INTENT(IN)              :: plit
    TYPE(params),   DIMENSION(1:mp,1:n),   INTENT(IN)              :: par
    REAL(r_2),       DIMENSION(1:mp,1:n),  INTENT(IN), OPTIONAL    :: qex
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT), OPTIONAL :: wex
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT),   OPTIONAL :: heads
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(IN),    OPTIONAL :: FS
    REAL(r_2),      DIMENSION(1:mp,0:n),   INTENT(INOUT), OPTIONAL :: ciso
    REAL(r_2),      DIMENSION(1:mp,0:n),   INTENT(INOUT), OPTIONAL :: cisoice
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT), OPTIONAL :: cisos, ciso0, cisoL
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: cprec
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: cali
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: qali
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT),   OPTIONAL :: qiso_in, qiso_out
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT),   OPTIONAL :: qiso_evap_cum, qiso_trans_cum
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT),   OPTIONAL :: qiso_liq_adv, qiso_vap_adv
    REAL(r_2),      DIMENSION(1:mp,1:n-1), INTENT(OUT),   OPTIONAL :: qiso_liq_diff, qiso_vap_diff
    REAL(r_2),      DIMENSION(1:mp,0:n),   INTENT(OUT),   OPTIONAL :: qvsig, qlsig, qvTsig, qvh
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT), OPTIONAL :: deltaTa
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: lE_old
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: dolitter       ! 0: no; 1: normal; 2: resistance
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: doisotopologue ! 0: no isotope; 1: HDO; 2: H218O
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: dosepts        ! 0: normal; 1: uncouple T & S
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: docondition    ! 0: no cond., 1: columns, 2: lines, 3: both
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: doadvection       ! 0: off; 1: onn

    ! Solves the RE and, optionally, the ADE from time ts to tfin.
    ! Definitions of arguments:
    ! Required args:
    ! ts   - start time (h).
    ! tfin   - finish time.
    ! qprec   - precipitation (or water input) rate (fluxes are in cm/h).
    ! qevap   - potl evaporation rate from soil surface.
    ! n    - no. of soil layers.
    ! nsol   - no. of solutes.
    ! dx(1:n) - layer thicknesses.
    ! h0   - surface head, equal to depth of surface pond.
    ! S(1:n)  - degree of saturation ("effective satn") of layers.
    ! evap   - cumulative evaporation from soil surface (cm, not initialised).
    ! runoff  - cumulative runoff.
    ! infil   - cumulative net infiltration (time integral of flux across surface).
    ! drn   - cumulative net drainage (time integral of flux across bottom).
    ! nsteps  - cumulative no. of time steps for RE soln.
    ! Optional args:
    ! heads(1:n)   - matric heads h of layers at finish.
    ! qexsub    - subroutine to get layer water extraction rates (cm/h) by
    !     plants. Note that there is no solute extraction and osmotic
    !     effects due to solute are ignored. Arguments:
    !     qex(1:n) - layer extraction rates; qexh(1:n) - partial
    !     derivs of qex wrt h.
    ! wex(1:n)    - cumulative water extraction from layers.
    ! cin(1:nsol)   - solute concns in water input (user's units/cc).
    ! c0(1:nsol)   - solute concns in surface pond.
    ! sm(1:n,1:nsol)  - solute (mass) concns in layers.
    ! soff(1:nsol)   - cumulative solute runoff (user's units).
    ! sinfil(1:nsol)  - cumulative solute infiltration.
    ! sdrn(1:nsol)   - cumulative solute drainage.
    ! nssteps(1:nsol) - cumulative no. of time steps for ADE soln.
    ! isosub    - subroutine to get adsorbed solute (units/g soil) from concn
    !     in soil water according to chosen isotherm code.
    !     Arguments: iso - 2 character code; c - concn in soil water;
    !     p(:) - isotherm parameters; f - adsorbed mass/g soil;
    !     fc - deriv of f wrt c (slope of isotherm curve). Note that
    !     linear adsorption does not require a sub, and other types
    !     are available in sub isosub.

    REAL(r_2),    DIMENSION(1:mp)       :: precip, qevap
    REAL(r_2),    DIMENSION(1:mp)       :: qL, qhL, qybL, qTbL, qhTbL, qhybL, rexcol, wcol, ql0, qv0
    LOGICAL,      DIMENSION(1:mp)       :: again, getq0,getqn,init
    LOGICAL,      DIMENSION(1:mp,1:n)   :: again_ice
    INTEGER(i_d), DIMENSION(1:mp)       :: ih0, iok, itmp, ns, nsat, nsatlast, nsteps0
    REAL(r_2),    DIMENSION(1:mp)       :: accel, dmax, dt, dwinfil, dwoff, fac, Khmin1, Kmin1, phimin1, phip
    REAL(r_2),    DIMENSION(1:mp)       :: qpme, rsig, rsigdt, sig, t
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: Sbot, Tbot
    REAL(r_2),    DIMENSION(1:mp,1:n-1) :: dz
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: hint, phimin, Khmin, qexd
    REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   :: aa, bb, cc, dd, ee, ff, gg, dy
    REAL(r_2),    DIMENSION(1:mp,-nsnow_max+1:n)   :: aah, bbh, cch, ddh, eeh, ffh, ggh, de
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: q, qya, qyb, qTa, qTb,qhya, qhyb, qhTa, qhTb
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: qadv, qadvya, qadvyb, qadvTa, qadvTb
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: qcond, qcondya, qcondyb, qcondTa, qcondTb
    REAL(r_2),    DIMENSION(1:mp,0:nsnow_max)   :: q_snow, qya_snow, qyb_snow, qTa_snow, qTb_snow, qhya_snow, qhyb_snow
    REAL(r_2),    DIMENSION(1:mp,0:nsnow_max)   :: qh_snow, qhTa_snow, qhTb_snow
    REAL(r_2),    DIMENSION(1:mp,0:nsnow_max)   :: qadv_snow, qadvyb_snow, qadvTb_snow

    TYPE(vars)                          :: vtmp
    TYPE(vars),   DIMENSION(1:mp,1:n)   :: var
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: qsig, qhsig, qadvsig
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: qliq, qv, qvT, qlya, qlyb, qvya, qvyb, qlTb, qvTa, qvTb
    REAL(r_2),    DIMENSION(1:mp,0:nsnow_max)   :: qsig_snow, qhsig_snow, qadvsig_snow
    TYPE(vars),   DIMENSION(1:mp,1:n)   :: vcall
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: deltaS, dTsoil
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: tmp2d1, tmp2d2
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: S0, Sliq0, Sliq, deltaSliq, cv0, deltacv
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: Sliqice0, Sliqice, deltaSliqice
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: Sice0, Sice, deltaSice
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: tmp_thetasat, tmp_thetar, tmp_tortuosity
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: delthetai, dthetaldT, thetal
    INTEGER(i_d), DIMENSION(1:mp,1:n)   :: isave, nsteps_ice, imelt


    TYPE(vars),         DIMENSION(1:mp) :: vtop, vbot
    TYPE(vars_aquifer), DIMENSION(1:mp) :: v_aquifer
    REAL(r_2),          DIMENSION(1:mp) :: qd, dwcol, dwdrainage, drn,inlit, dwinlit, drexcol, dwdischarge
    REAL(r_2),          DIMENSION(1:mp) :: dJcol_latent_S, dJcol_latent_T, dJcol_sensible
    REAL(r_2),          DIMENSION(1:mp,1:n):: deltaJ_latent_S, deltaJ_latent_T, deltaJ_sensible_S, deltaJ_sensible_T
    REAL(r_2),          DIMENSION(1:mp) :: qevapsig
    REAL(r_2),          DIMENSION(1:mp) :: qrunoff
    REAL(r_2),          DIMENSION(1:mp) :: tmp1d1, tmp1d2, tmp1d3,  tmp1d4
    REAL(r_2),          DIMENSION(1:mp) :: deltah0
    REAL(r_2),          DIMENSION(1:mp) :: Tsurface_pot, Epot, Hpot, Gpot, dEdrha, dEdTa, dEdTsoil, dGdTa, dGdTsoil
    REAL(r_2),          DIMENSION(1:mp) :: SL0, deltaSL, cvL0, SLliq0, deltacvL, SLliq, deltaSLliq
    REAL(r_2),          DIMENSION(1:mp) :: qiso_evap, qiso_trans
    REAL(r_2),          DIMENSION(1:mp) :: lE0, G0, E_vap, E_liq, dE_vapdT1
    REAL(r_2),          DIMENSION(1:mp) :: Tfreezing, dT0
    REAL(r_2),          DIMENSION(1:mp) :: dtdT
    REAL(r_2),          DIMENSION(1:mp,0:n) :: LHS, RHS,LHS_h, RHS_h
    REAL(r_2),          DIMENSION(1:mp,1:nsnow_max) :: LHS_snow, RHS_snow,LHS_h_snow, RHS_h_snow
    INTEGER(i_d),       DIMENSION(1:mp) :: pondcase

    INTEGER(i_d),       DIMENSION(1:mp) :: nns, iflux
    LOGICAL      :: litter
    INTEGER(i_d) :: i, j, k, kk, condition
    INTEGER(i_d) :: littercase, isotopologue, advection, septs ! switches
    REAL(r_2)    :: ztmp, c2, theta
    REAL(r_2)    :: dTqwdTa, dTqwdTb, Tqw, keff
    REAL(r_2),          DIMENSION(1:mp) :: cp, cpeff, hice, deltahice,h0_0, hice_0, h0_tmp, hice_tmp
    REAL(r_2),          DIMENSION(1:mp,nsnow_max) :: qmelt, hsnow
    REAL(r_2),          DIMENSION(1:mp,nsnow_max) :: delta_snowcol, delta_snowT, delta_snowliq, dTsnow
    REAL(r_2),          DIMENSION(1:mp) :: melt ! cumulative loss of snow pack as melt water
    REAL(r_2),      DIMENSION(1:mp,1:n)       :: thetai_0, J0
    REAL(r_2) :: tmp1, tmp2
    REAL(r_2),          DIMENSION(1:mp,1:n) :: iqex, thetal_max
    REAL(r_2),          DIMENSION(1:mp)     :: icali
    INTEGER(i_d),       DIMENSION(1:mp) :: nfac1, nfac2, nfac3, nfac4, nfac5, &
         nfac6, nfac7, nfac8, nfac9, nfac10
    REAL(r_2),          DIMENSION(1:mp)     :: v1, v2, lnSliq, kmin, ipi


    !open (unit=7, file="Test.out", status="replace", position="rewind")
    ! The derived types params and vars hold soil water parameters and variables.
    ! Parameter names often end in e, which loosely denotes "air entry", i.e.,
    ! values at h=he. While values of water content th and hydraulic conductivity K
    ! at h=he are equal to those at saturation, the derivs wrt S are nonzero. The
    ! MFP phi for h>he is given by phi=phie+Ke*(h-he). The saturation status of a
    ! layer is stored as 0 or 1 in isat since S may be >1 (because of previous
    ! overshoot) when a layer desaturates. Fluxes at the beginning of a time step
    ! and their partial derivs wrt S or phi of upper and lower layers or boundaries
    ! are stored in q, qya and qyb.

    ! set switches
    if (present(dolitter)) then
       littercase = dolitter
    else
       littercase = 0
    endif
    if (present(doadvection)) then
       advection = doadvection
    else
       advection = 0
    endif
    if (littercase > 2) then
       write(*,*) 'dolitter not in [0-2]: ', littercase
       stop
    endif

    if (present(doisotopologue)) then
       isotopologue = doisotopologue
    else
       isotopologue = 0
    endif
    if (isotopologue > 2) then
       write(*,*) 'doisotopologue not in [0-2]: ', isotopologue
       stop
    endif
    if (isotopologue /= 0 .and. (.not. present(ciso))) then
       write(*,*) 'doisotopologue /= 0 but no ciso present.'
       stop
    endif

    if (present(dosepts)) then
       septs = dosepts
    else
       septs = 0
    endif
    if (septs > 1) then
       write(*,*) 'dosepts not in [0-1]: ', septs
       stop
    endif

    if (present(docondition)) then
       condition = docondition
    else
       condition = 0
    endif
    if (condition < 0 .or. condition > 3) then
       write(*,*) 'docondition not in [0-3]: ', condition
       stop
    endif

    if (present(qex)) then
       iqex = qex
    else
       iqex = zero
    endif

    if (present(qali) .and. present(cali)) then
       where (qali>zero)
          icali = cali
       elsewhere
          icali = zero
       endwhere
    else
       icali = zero
    endif

    ! global counters
    if (.not. allocated(nless)) allocate(nless(mp))
    nless(:) = 0
    if (.not. allocated(n_noconverge)) allocate(n_noconverge(mp))
    n_noconverge(:) = 0

    ! set solve_type for numerical derivatives
    if (.not. allocated(sol)) call setsol(mp)

    ! initialise cumulative variables
    wcol(:)      = zero
    Jcol_sensible(:) = zero
    Jcol_latent_S(:) = zero
    Jcol_latent_T(:) = zero
    deltaice_cum_T(:) = zero
    deltaice_cum_S(:) = zero
    deltaJ_sensible_S(:,:) = zero
    deltaJ_sensible_T(:,:) = zero
    deltaJ_latent_S(:,:) = zero
    deltaJ_latent_T(:,:) = zero
    drainage(:)  = zero
    discharge(:) = zero
    infil(:)     = zero
    inlit(:)     = zero
    dwinlit(:)   = zero
    evap(:)      = zero
    evap_pot(:)  = zero
    runoff(:)    = zero
    melt(:) = zero
    rexcol(:)    = zero
    Hcum(:)      = zero
    Gcum(:)      = zero
    lEcum(:)     = zero
    Qadvcum(:)  = zero
    wex(:,:)          = zero
    precip(:)         = zero
    drn(:)            = zero
    if (isotopologue /= 0) then
       qiso_evap_cum(:)  = zero
       qiso_trans_cum(:) = zero
    endif
    deltah0(:) = zero
    ! zero var-structure that contains all the hydrological variables
    vtmp = zerovars()
    vtmp%h       = one
    vtmp%lambdav = rlambda
    vtmp%lambdaf = lambdaf
    ! zero vars at the bottom and top of the soil column
    vtop = spread(vtmp,1,mp)
    vbot = spread(vtmp,1,mp)
    vlit = spread(vtmp,1,mp)
    ! Vanessa: try this with limited stacksize, otherwise the double-loop
    var = spread(vtop,2,n)
    ! do i=1, mp
    !    do k=1, n
    !       var(i,k) = vtmp
    !    end do
    ! end do
    hint(:,:)   = zero
    phimin(:,:) = zero
    q(:,:)      = zero
    qya(:,:)    = zero
    qyb(:,:)    = zero
    qTa(:,:)    = zero
    qTb(:,:)    = zero
    qhya(:,:)   = zero
    qhyb(:,:)   = zero
    qhTa(:,:)   = zero
    qhTb(:,:)   = zero
    qadvyb(:,:) = zero
    qadvya(:,:) = zero
    qcondyb(:,:) = zero
    qcondya(:,:) = zero
    qcondTb(:,:) = zero
    qcondTa(:,:) = zero
    qya_snow(:,:)  = zero
    qyb_snow(:,:)  = zero
    qTa_snow(:,:)  = zero
    qTb_snow(:,:)  = zero
    qhya_snow(:,:) = zero
    qhyb_snow(:,:) = zero
    qhTa_snow(:,:) = zero
    qhTb_snow(:,:) = zero
    aa(:,:)     = zero
    aah(:,:)    = zero
    bb(:,:)     = zero
    bbh(:,:)    = zero
    cc(:,:)     = zero
    cch(:,:)    = zero
    dd(:,:)     = zero
    ddh(:,:)    = zero
    ee(:,:)     = zero
    eeh(:,:)    = zero
    ff(:,:)     = zero
    ffh(:,:)    = zero
    gg(:,:)     = zero
    ggh(:,:)    = zero
    dy(:,:)     = zero
    dTsoil(:,:) = zero
    de(:,:)            = zero
    keff               = zero
    dTqwdTa            = zero
    dTqwdTb            = zero
    Tqw                = zero
    ztmp               = zero
    c2                 = zero
    theta              = zero
    cp                 = zero
    cpeff              = zero
    hice               = zero
    deltahice          = zero
    h0_0               = zero
    hice_0             = zero
    h0_tmp             = zero
    hice_tmp           = zero
    qmelt(:,:)         = zero
    hsnow(:,:)         = zero
    delta_snowcol(:,:) = zero
    delta_snowT(:,:)   = zero
    delta_snowliq(:,:) = zero
    dTsnow(:,:)        = zero
    melt(:)            = zero
    thetai_0(:,:)      = zero
    J0(:,:)            = zero
    dT0(:)             = zero
    thetal_max         = zero
    nsteps             = 0

    ! initialise snow Ebal diagnostics
    vsnow(:)%Qadv_rain = zero
    vsnow(:)%Qadv_snow = zero
    do i=1, nsnow_max
       vsnow(:)%Qadv_vap(i)       = zero
       vsnow(:)%Qcond_net(i)      = zero
       vsnow(:)%deltaJlatent(i)   = zero
       vsnow(:)%deltaJsensible(i) = zero
       vsnow(:)%Qadv_transfer(i)  = zero
       vsnow(:)%Qadv_melt(i)      = zero
       vsnow(:)%FluxDivergence(i) = zero
    enddo

    litter = .false.
    if (littercase == 1) litter=.true. ! full litter model

    qexd(:,:) = zero
    phip(:)   = zero !max(par(:,1)%phie-par(:,1)%he*par(:,1)%Ke, 1.00001_r_2*par(:,1)%phie) ! phi at h=0

    ! get K, Kh and phi at hmin (hmin is smallest h, stored in hy-props)
    call hyofh(spread(hmin,1,mp), par(:,1), Kmin1(:), Khmin1(:), phimin1(:))

    dz(:,:) = half*(dx(:,1:n-1)+dx(:,2:n)) ! flow paths

    !----- set up for boundary conditions
    getq0(:) = .false.
    getqn(:) = .false.
    if (botbc == "constant head") then ! h at bottom bdry specified
       getqn(:)  = .true.
       tmp1d1(:)  = hbot
       ! for hbot < he
       Sbot(:,:) = spread(Sofh(tmp1d1,par(:,n)),1,n)

       Tbot(:,:) = spread(Tsoil(:,n),1,n)
       count_hyofS = count_hyofS + 1
       ! Debug for mp=1: remove elemental from hyofS and do loop instead of next line
       call hyofS(Sbot, Tbot, par, vcall)
       ! do i=1, n
       !    call hyofS(Sbot(1,i), Tbot(1,i), par(1,i), vcall(1,i))
       ! end do
       ! End debug hyofS
       ! for hbot >= he
       vtmp = zerovars()
       vtmp%isat    = 1
       vtmp%h       = hbot
       vtmp%rh      = one
       vtmp%lambdav = rlambda
       vtmp%lambdaf = lambdaf
       vbot = spread(vtmp,1,mp)
       vbot(:)%phi = (hbot-par(:,n)%he)*par(:,n)%Ke+var(:,n)%phie
       vbot(:)%K   = par(:,n)%Ke
       where (par(:,n)%he > hbot)
          vbot(:)      = vcall(:,n)
          vbot(:)%isat = 0
       endwhere
    end if
    !----- end set up for boundary conditions

    !----- initialise
    nfac1  = 0
    nfac2  = 0
    nfac3  = 0
    nfac4  = 0
    nfac5  = 0
    nfac6  = 0
    nfac7  = 0
    nfac8  = 0
    nfac9  = 0
    nfac10 = 0
    t(:)       = ts
    nsteps0(:) = nsteps
    nsat(:)    = 0
    qd(:)      = zero
    ! initialise saturated regions
    var(:,:)%isat  = 0
    where (S(:,:) >= one)
       var(:,:)%K    = par(:,:)%Ke
       var(:,:)%isat = 1
    endwhere
    !    Sliq = S
    !    where ((Tsoil<Tfrz(S,par%he,one/par%lam)) .and. (var%isat == 1))
    !       thetal_max    = thetalmax(Tsoil,S,par%he,one/par%lam,par%thre,par%the)
    !       Sliq          = (thetal_max - (par%the-par%thre))/par%thre
    !       var(:,:)%phie = par%phie*exp(-log(Sliq)/par%lam)*exp(par%eta*log(Sliq))
    !       var(:,:)%phi  = par%phie*exp(-log(Sliq)/par%lam)*exp(par%eta*log(Sliq))
    !    endwhere
    !    where(h0(:)>zero)
    !       hice(:) = h0(:)*(S(:,1)-Sliq(:,1))
    !       var(:,1)%phi = max((var(:,1)%phie -var(:,1)%he*var(:,1)%Ksat), &
    !            (one+e5)*var(:,1)%phie)+(h0(:)-hice(:))*var(:,1)%Ksat
    !    endwhere

    vlit(:)%isat = 0
    where (SL(:) >= one) vlit(:)%isat = 1

    ! ! initialise acquifer
    ! v_aquifer(:)%zsoil  = sum(dx(:,:),2)
    ! v_aquifer(:)%zdelta = zdelta(:)
    ! call aquifer_props(v_aquifer(:))

    ! initialise litter
    if (littercase == 1 .or. littercase == 2) then
       call litter_props(Sl(:), Tl(:), vlit(:), plit(:), h0(:))
    endif
    ! Add resistance through litter for simple litter model
    if (littercase == 2) then
       vmet(:)%rbw = vmet(:)%rbw + dxL(:)/vlit(:)%Dv
       ztmp        = one/rhocp
       vmet(:)%rbh = vmet(:)%rbh + dxL(:)/(vlit(:)%kth*ztmp)
    endif

    lE0(:) = lE_old(:) ! used for initial guess of litter temperature
    !----- end initialise

    !----- solve until tfin
    init(:) = .true. ! flag to initialise h at soil interfaces
    do kk=1, mp
       do while (t(kk) < tfin)

          !----- take next time step
          iflux(kk)=1
          again(kk)  = .true. ! flag for recalcn of fluxes (default=false)
          imelt = 0 ! initialise imelt (==1 at onset of melting)
          qmelt(kk,:) = zero
          vsnow(kk)%nsnow_last = vsnow(kk)%nsnow ! for detecting onset and disappearance of dedicated snow layer
          do while (again(kk)) ! sometimes need twice to adjust phi at satn

             nsatlast(kk) = nsat(kk) ! for detecting onset of profile saturation
             nsat(kk)     = sum(var(kk,:)%isat,1) ! no. of sat layers
             sig(kk)      = half
             if (nsat(kk) /= 0) sig(kk) = one ! time weighting sigma
             rsig(kk)     = one/sig(kk)

             ! update variables
             if (iflux(kk)==1) then
                ! Calc flux matric potentials (and derivatives) from S
                ! this set var-structure
                isave(kk,:) = var(kk,:)%isat
                count_hyofS = count_hyofS + 1 ! counter for calls to hyofS
                ! Debug for mp=1: remove elemental from hyofS and do loop instead of next line
                call hyofS(S(kk,:), Tsoil(kk,:), par(kk,:), var(kk,:)) ! for layers where S<1
                ! do i=1, n
                !    call hyofS(S(kk,i), Tsoil(kk,i), par(kk,i), var(kk,i))
                ! end do
                ! End debug hyofS
                cp(kk) = real(1-var(kk,1)%iice,r_2)*cswat*rhow & ! heat capacity of pond
                     + real(var(kk,1)%iice,r_2)*rhow* &
                     ((one-var(kk,1)%thetai/par(kk,1)%thre)*cswat + (var(kk,1)%thetai/par(kk,1)%thre)*csice)
                cpeff(kk) = cp(kk) + rhow*lambdaf*var(kk,1)%dthetaldT/par(kk,1)%thre
                phip(kk) = max(var(kk,1)%phie-var(kk,1)%he*var(kk,1)%Ksat, (one+e5)*var(kk,1)%phie) !at onset of ponding
                var(kk,:)%isat  = isave(kk,:)
                thetai(kk,:)    = var(kk,:)%thetai ! ice profile
                thetal(kk,:)    = var(kk,:)%thetal ! liq water profile
                thetai_0(kk,:)  = thetai(kk,:) ! initial ice profile
                dthetaldT(kk,:) = var(kk,:)%dthetaldT
                hice(kk)   = h0(kk)*var(kk,1)%thetai/par(kk,1)%thre
                hice_0(kk) = hice(kk)
                h0_0(kk)   = h0(kk)
                ! sensible heat stored in top layer + pond
                Jsensible(kk,1)   = (var(kk,1)%csoil* dx(kk,1)+h0(kk)*cp(kk))*(Tsoil(kk,1))
                ! sensible heat stored in soil column (2:n)
                Jsensible(kk,2:n) = var(kk,2:n)%csoil*(Tsoil(kk,2:n))* dx(kk,2:n)

                Sice0(kk,1:n)    = var(kk,1:n)%thetai/par(kk,1:n)%thre
                Sliqice0(kk,1:n) = (S(kk,1:n) - var(kk,1:n)%cv)/(one-var(kk,1:n)%cv)
                Sliq0(kk,1:n)    = Sliqice0(kk,1:n) - Sice0(kk,1:n)
                Sice0(kk,1) = Sice0(kk,1)  + hice_0(kk)/(dx(kk,1)*par(kk,1)%thre)
                Sliq0(kk,1) = Sliq0(kk,1) + (h0_0(kk)-hice_0(kk))/(dx(kk,1)*par(kk,1)%thre) ! add pond component to Sliq(kk,1)

                CALL snow_adjust(irec,mp,n,kk,ns,h0,h0_0,hice,hice_0,thetai,dx,vsnow,var,par,S,Tsoil, &
                     Jcol_latent_S, Jcol_latent_T, Jcol_sensible,deltaJ_sensible_S,qmelt)
             endif ! iflux==1

             ! initialise litter vars
             if (iflux(kk)==1 .and. (littercase==1 .or. littercase == 2)) then
                call litter_props(Sl(kk), Tl(kk), vlit(kk), plit(kk), h0(kk))
             endif

             ! ! phi is solution var at satn, so h calc from phi where S>=1 - done in hyofS above for S<1
             ! where (S(kk,:) >= one) & ! (special for frozen soil)
             !      var(kk,:)%h = var(kk,:)%he + (var(kk,:)%phi-var(kk,:)%phie)/var(kk,:)%Ksat

             !----- get fluxes and derivs
             ! get surface condition
             ! ns==1 if no pond or full pond, i.e. do not solve for change in pond height
             ! ns==0 then change for pond height
             if ((var(kk,1)%phi <= phip(kk) .and. h0(kk) <= zero .and. nsat(kk) < n).or.(var(kk,1)%isat==0)) then ! no ponding
                ns(kk)    = 1 ! start index for eqns
             else ! ponding
                ns(kk)    = 0
                TL(kk)    = vmet(kk)%Ta ! initialise pond T
                vtop(kk)%isat    = 1
                vtop(kk)%h       = h0(kk)
                vtop(kk)%phi     = max((var(kk,1)%phie -var(kk,1)%he*var(kk,1)%Ksat), &
                     (one+e5)*var(kk,1)%phie)+(h0(kk)-hice(kk))*var(kk,1)%Ksat
                vtop(kk)%K       = var(kk,1)%Ksat
                vtop(kk)%rh      = one
                vtop(kk)%lambdav = rlambda
                vtop(kk)%lambdaf = lambdaf
                ! var(kk,1)%phi = var(kk,1)%phie
                ! calculates phi1,eff (pond + top soil layer)
                call flux(par(kk,1), vtop(kk), var(kk,1), half*dx(kk,1), &
                     q(kk,0), qya(kk,0), qyb(kk,0), qTa(kk,0), qTb(kk,0))
             endif

             if (Sl(kk) >= one) then
                vlit(kk)%isat = 1
                vlit(kk)%h    = plit(kk)%he
             endif

             ! get bottom boundary condn
             if (botbc=="seepage") then
                vtmp = zerovars()
                vbot = spread(vtmp,1,mp)
                if (var(kk,n)%h > -half*gf*dx(kk,n)) then
                   getqn(kk) = .true.
                   vbot(kk)%isat    = 1
                   vbot(kk)%phi     = (zero-var(kk,n)%he)*var(kk,n)%Ksat+var(kk,n)%phie ! (special for frozen soil)
                   vbot(kk)%K       = var(kk,n)%Ksat
                   vbot(kk)%rh      = one
                   vbot(kk)%lambdav = rlambda
                   vbot(kk)%lambdaf = lambdaf
                else
                   getqn(kk) = .false.
                endif
             end if

             ! Fluxes and derivatives at the surfaces (air/litter, litter/soil, and similar)
             ! three cases
10           pondcase(kk) = 0
             if (ns(kk)==1 .and. (.not. litter) ) pondcase(kk) = 1 ! no litter and no pond
             if (ns(kk)==1 .and.        litter ) pondcase(kk) = 2 ! litter and no pond
             if (ns(kk)==0 .and.vsnow(kk)%nsnow==0) pondcase(kk) = 3 ! pond with and without litter
             if (vsnow(kk)%nsnow>0) pondcase(kk) = 4
             sol(kk)%Ks = var(kk,1)%Ksat
             select case (pondcase(kk))
             case (1) ! no litter and no pond
                ! structure that can be seen by rh0_sol_0d etc. which is used in rtbis (bisection)
                ! Look in paper Appendix A.6: upper boundary condition: no ponding or litter
                ! it just equates fluxes towards and away from the interface; 1 for water, 1 for heat
                ! have to solve numerically because of non-linearity in phi(h)
                sol(kk)%T1      = Tsoil(kk,1)
                sol(kk)%Ta      = vmet(kk)%Ta
                sol(kk)%cva     = vmet(kk)%cva*thousand
                sol(kk)%Rnet    = vmet(kk)%Rn
                sol(kk)%hr1     = var(kk,1)%rh
                sol(kk)%hra     = vmet(kk)%rha
                sol(kk)%Dv      = var(kk,1)%Dv
                sol(kk)%gv      = one/vmet(kk)%rbw
                sol(kk)%gh      = one/vmet(kk)%rrc * rhocp
                sol(kk)%Dh      = var(kk,1)%kth
                sol(kk)%dz      = half*dx(kk,1)
                sol(kk)%phie    = var(kk,1)%phie
                sol(kk)%he      = var(kk,1)%he
                sol(kk)%K1      = var(kk,1)%K
                sol(kk)%eta     = par(kk,1)%eta
                sol(kk)%lambda  = par(kk,1)%lam
                sol(kk)%lambdav = var(kk,1)%lambdav

                call potential_evap(vmet(kk)%Rn, vmet(kk)%rrc, vmet(kk)%rbw, vmet(kk)%Ta, vmet(kk)%rha, &
                     Tsoil(kk,1), var(kk,1)%kth, half*dx(kk,1)+h0(kk),sol(kk)%lambdav, Tsurface_pot(kk), Epot(kk), Hpot(kk), &
                     Gpot(kk), dEdrha(kk), dEdTa(kk), dEdTsoil(kk), dGdTa(kk), dGdTsoil(kk))
                if (sol(kk)%Dv > 1.e-12_r_2) then
                   E_vap(kk) = (sol(kk)%hr1*csat(sol(kk)%T1)-sol(kk)%cva)/(1/sol(kk)%gv + sol(kk)%dz/sol(kk)%Dv)*sol(kk)%lambdav
                   dE_vapdT1(kk) = (sol(kk)%hr1*slope_csat(sol(kk)%T1))/(1/sol(kk)%gv + sol(kk)%dz/sol(kk)%Dv)*sol(kk)%lambdav
                else
                   E_vap(kk) = zero
                   dE_vapdT1(kk) = zero
                endif
                call hyofh(hmin, par(kk,1), kmin(kk), Khmin(kk,1), phimin(kk,1)) ! get phi at hmin
                E_liq(kk) = ((var(kk,1)%phi-phimin(kk,1))/sol(kk)%dz-sol(kk)%K1)*thousand*sol(kk)%lambdav
                lE0(kk) = min(Epot(kk),E_vap(kk)+ E_liq(kk)) ! analytic approximation (See Haverd et al. 2012, Appxx)
                qevap(kk) = lE0(kk)/(thousand*sol(kk)%lambdav)
                q(kk,0)  = qprec(kk) + qprec_snow(kk) - qevap(kk)

                T0(kk) = (-sol(kk)%dz*lE0(kk) + sol(kk)%dz*sol(kk)%Rnet + &
                     sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta) &
                     /(sol(kk)%Dh + sol(kk)%dz*sol(kk)%gh)

                G0(kk)    = sol(kk)%Dh/sol(kk)%dz*(T0(kk)-sol(kk)%T1)
                Tsurface(kk) = T0(kk)
                qcond(kk,0) = G0(kk)
                qh(kk,0) = qcond(kk,0)

                ! derivatives
                ! q refers to moisture in numerator
                ! qh refers to heat in numerator
                ! y refers to moisture in denominator
                ! T refers to temperature in denominator
                ! a refers to the layer above
                ! b refers to the layer below

                ! initialise derivatives to zero
                qya(kk,0)  = zero
                qTa(kk,0)  = zero
                qhya(kk,0) = zero
                qhyb(kk,0) = zero
                qhTa(kk,0) = zero

                qybL(kk) = zero
                qL(kk)   = zero
                qTbL(kk) = zero

                ! liquid and vapour fluxes

                qlya(kk,0) = zero
                qlyb(kk,0) = zero

                qvya(kk,0) = qya(kk,0)
                qvyb(kk,0) = qyb(kk,0)

                if (Epot(kk)<=(E_vap(kk)+ E_liq(kk))) then ! potential evap independent of S(1), dependent on T1
                   qyb(kk,0) = zero
                   qTb(kk,0) = -dEdTsoil(kk)/(thousand*sol(kk)%lambdav)
                else ! supply limited
                   qTb(kk,0) = -dE_vapdT1(kk)/(thousand*sol(kk)%lambdav)
                   qyb(kk,0) = -(var(kk,1)%phiS/sol(kk)%dz - var(kk,1)%KS)  !!vh!! include vapour component??
                endif

                if (Epot(kk)<=(E_vap(kk)+ E_liq(kk))) then ! potential evap independent of S(1), dependent on T1
                   qliq(kk,0) = -Epot(kk)/(thousand*sol(kk)%lambdav)
                   qv(kk,0) = zero
                   qlyb(kk,0) = zero
                   qvyb(kk,0) = qyb(kk,0)
                   qlTb(kk,0) = qTb(kk,0)
                   qvTb(kk,0) = zero
                else ! supply limited
                   qliq(kk,0) = -E_liq(kk)/(thousand*sol(kk)%lambdav)
                   qv(kk,0) = -E_vap(kk)/(thousand*sol(kk)%lambdav)
                   qlyb(kk,0) = -(var(kk,1)%phiS/sol(kk)%dz - var(kk,1)%KS)
                   qvyb(kk,0) = zero
                   qlTb(kk,0) = zero
                   qvTb(kk,0) = qTb(kk,0)
                endif
                qvya(kk,0) = qya(kk,0)
                qlya(kk,0) = zero

                ! end of partial derivative evaluation

                ! advective component of heat flux
                qadv(kk,0) = rhow*cswat*qprec(kk)*(vmet(kk)%Ta) + rhow*csice*qprec_snow(kk)*(min(vmet(kk)%Ta,zero)) &
                     - rhow*qprec_snow(kk)*lambdaf

                Tqw  = merge(vmet(kk)%Ta, Tsoil(kk,1), (-qevap(kk))>zero)
                dTqwdTb = merge(zero, one, (-qevap(kk))>zero)
                qadv(kk,0) = qadv(kk,0) + rhow*cswat*Tqw*(-qevap(kk))
                qadvTb(kk,0) = dTqwdTb + rhow*cswat*Tqw*qTb(kk,0)
                qadvyb(kk,0) =  rhow*cswat*qyb(kk,0)*Tqw
                qadvya(kk,0) =  zero
                qadvTa(kk,0) =  zero
                qh(kk,0) = qadv(kk,0) + qcond(kk,0)
                qhyb(kk,0) = qcondyb(kk,0) +  qadvyb(kk,0)
                qhTb(kk,0) = qcondTb(kk,0) + qadvTb(kk,0)

             case (2) ! litter and no pond
                !!vh!!
                ! need to recode using analytic approx

             case (3) ! pond with and without litter
                deltaTa(kk) = zero
                call potential_evap(vmet(kk)%Rn, vmet(kk)%rrc, vmet(kk)%rbw, vmet(kk)%Ta, vmet(kk)%rha, &
                     Tsoil(kk,1), var(kk,1)%kth, half*dx(kk,1)+h0(kk),var(kk,1)%lambdav, &
                     Tsurface_pot(kk), Epot(kk), Hpot(kk), &
                     Gpot(kk), dEdrha(kk), dEdTa(kk), dEdTsoil(kk), dGdTa(kk), dGdTsoil(kk))
                Tsurface(kk)  = Tsurface_pot(kk)
                T0(kk)        = Tsurface_pot(kk)
                rh0(kk)       = one
                rhsurface(kk) = one
                G0(kk)        = Gpot(kk)
                lE0(kk)       = Epot(kk)
                qevap(kk)     = Epot(kk)/(thousand*var(kk,1)%lambdav)
                qpme(kk)        = qprec_snow(kk) + qprec(kk) - qevap(kk)
                q(kk,0) = qprec_snow(kk)+ qprec(kk) - qevap(kk)
                if (ns(kk)==0) then
                   ! pond included in top soil layer: calculate phi1(eff)
                   var(kk,1)%phi = max((var(kk,1)%phie -var(kk,1)%he*var(kk,1)%Ksat), &
                        (one+e5)*var(kk,1)%phie)+(h0(kk)-hice(kk))*var(kk,1)%Ksat
                   q(kk,0) = qpme(kk)
                   qya(kk,0) = zero
                   qyb(kk,0) = zero
                   qTa(kk,0) = zero
                   qTb(kk,0) = zero
                endif

                qh(kk,0)   = G0(kk)
                qcond(kk,0)   = G0(kk)
                qh(kk,0)   = qcond(kk,0)
                qhya(kk,0) = zero
                qhyb(kk,0) = zero
                qhTa(kk,0) = zero
                qcondTb(kk,0)  = dGdTsoil(kk)
                qhTb(kk,0) = qcondTb(kk,0)

                qL(kk)    = zero
                qybL(kk)  = zero
                qTbL(kk)  = zero
                qhL(kk)   = zero
                qhTbL(kk) = zero
                qhybL(kk) = zero

                ql0(kk)   = -qevap(kk)
                qv0(kk)   = zero
                if (litter) then
                   qL(kk)  = qpme(kk)
                   qhL(kk) = G0(kk)
                endif

                ! liquid and vapour fluxes
                qyb(kk,0) = zero
                qTb(kk,0) = -dEdTsoil(kk)/(thousand*var(kk,1)%lambdav)
                qlya(kk,0) = zero
                qvya(kk,0) = zero

                qliq(kk,0) = -Epot(kk)/(thousand*var(kk,1)%lambdav)
                qv(kk,0) = zero
                qlyb(kk,0) = zero
                qvyb(kk,0) = qyb(kk,0)
                qlTb(kk,0) = qTb(kk,0)
                qvTb(kk,0) = zero

                ! advective component of heat flux

                qadv(kk,0) = rhow*cswat*qprec(kk)*(vmet(kk)%Ta) + rhow*csice*qprec_snow(kk)*(min(vmet(kk)%Ta,zero)) &
                     - rhow*qprec_snow(kk)*lambdaf

                Tqw  = merge(vmet(kk)%Ta, Tsoil(kk,1), (-qevap(kk))>zero)
                dTqwdTb = merge(zero, one, (-qevap(kk))>zero)
                qadv(kk,0) = qadv(kk,0) + rhow*cswat*Tqw*(-qevap(kk))
                qadvTb(kk,0) = dTqwdTb + rhow*cswat*Tqw*qTb(kk,0)
                qadvyb(kk,0) =  rhow*cswat*qyb(kk,0)*Tqw
                qadvya(kk,0) =  zero
                qadvTa(kk,0) =  zero
                qh(kk,0) = qadv(kk,0) + qcond(kk,0)
                qhyb(kk,0) = qcondyb(kk,0) +  qadvyb(kk,0)
                qhTb(kk,0) = qcondTb(kk,0) + qadvTb(kk,0)

             case(4) !dedicated snow layer
                ! SEB at snow/air interface
                if (vsnow(kk)%hliq(vsnow(kk)%nsnow)>zero) then
                   sol(kk)%lambdav = lambdaf
                else
                   sol(kk)%lambdav = lambdas
                endif
                call potential_evap(vmet(kk)%Rn, vmet(kk)%rrc, vmet(kk)%rbw, vmet(kk)%Ta, vmet(kk)%rha, &
                     vsnow(kk)%tsn(vsnow(kk)%nsnow), vsnow(kk)%kth(vsnow(kk)%nsnow), half*vsnow(kk)%depth(vsnow(kk)%nsnow), &
                     sol(kk)%lambdav, Tsurface_pot(kk), Epot(kk), Hpot(kk), &
                     Gpot(kk), dEdrha(kk), dEdTa(kk), dEdTsoil(kk), dGdTa(kk), dGdTsoil(kk))

                qevap(kk) = Epot(kk)/(thousand*sol(kk)%lambdav)
                Tsurface(kk) = Tsurface_pot(kk)

                ! moisture flux at air/snow interface
                q_snow(kk,0) = qprec_snow(kk)+qprec(kk)-qevap(kk)
                qyb_snow(kk,0) = zero
                qTb_snow(kk,0) = -dEdTsoil(kk)/(thousand*sol(kk)%lambdav)

                ! conductive heat flux at air/snow interface
                qh_snow(kk,0) = Gpot(kk)
                qhTb_snow(kk,0) = dGdTsoil(kk)
                if (vsnow(kk)%hliq(vsnow(kk)%nsnow)>zero) then
                   qhTb_snow(kk,0) = zero
                   qTb_snow(kk,0) = zero
                endif

                ! advective heat flux at air/snow interface
                qadv_snow(kk,0) = rhow*(qprec_snow(kk))*(csice*(min(vmet(kk)%Ta,zero))-lambdaf) + &
                     rhow*(qprec(kk))*cswat*(vmet(kk)%Ta)
                Tqw  = merge(vmet(kk)%Ta, vsnow(kk)%tsn(vsnow(kk)%nsnow), -qevap(kk)>zero)
                dTqwdTb = merge(zero,one, -qevap(kk)>zero)
                if (vsnow(kk)%hliq(vsnow(kk)%nsnow)>zero) then
                   qadv_snow(kk,0) = qadv_snow(kk,0) + rhow*(-qevap(kk))*cswat*Tqw
                   qadvTb_snow(kk,0) = zero
                else
                   ! qadv_snow(kk,0) = qadv_snow(kk,0) + rhow*(-qevap(kk))*(csice*Tqw-lambdaf)
                   ! qadvTb_snow(kk,0) = rhow*csice*-qevap(kk)*dTqwdTb  +  rhow*csice*Tqw*qTb_snow(kk,0)
                   qadv_snow(kk,0) = qadv_snow(kk,0) + rhow*(-qevap(kk))*cswat*Tqw
                   qadvTb_snow(kk,0) = rhow*cswat*(-qevap(kk))*dTqwdTb  +  rhow*cswat*Tqw*qTb_snow(kk,0)

                endif
                ! total heat flux at air/snow interface
                qh_snow(kk,0) = qh_snow(kk,0) + qadv_snow(kk,0)
                qhTb_snow(kk,0) = qhTb_snow(kk,0) + qadvTb_snow(kk,0)
                qadvyb_snow(kk,0) = zero
                qhyb_snow(kk,0) = zero

                ! vapour flux at soil/snow interface
                if (var(kk,1)%isat==1) then
                   keff = vsnow(kk)%kE(1)/(thousand*lambdaf)/vsnow(kk)%depth(1)/2_r_2
                endif
                if (var(kk,1)%isat==0) then
                   keff = 2_r_2*((vsnow(kk)%kE(1)/(thousand*lambdaf))*(var(kk,1)%kE/thousand/var(kk,1)%lambdav))/ &
                        ((vsnow(kk)%kE(1)/(thousand*lambdaf))*dx(kk,1)+(var(kk,1)%kE/thousand/var(kk,1)%lambdav)* &
                        vsnow(kk)%depth(1))
                endif
                q(kk,0) = keff*(vsnow(kk)%tsn(1)-Tsoil(kk,1))
                qTa(kk,0) = keff
                qTb(kk,0) =-keff
                qya(kk,0) = zero
                qyb(kk,0) = zero
                qv(kk,0)   = q(kk,0)
                qvyb(kk,0) = qyb(kk,0)
                qvTb(kk,0) = qTb(kk,0)
                qliq(kk,0) = zero
                qlyb(kk,0) = zero
                qlTb(kk,0) = zero
                if (vsnow(kk)%hliq(1)>zero) then
                   qcondTa(kk,0) = zero
                   qhTa(kk,0) = zero
                   qTa(kk,0) = zero
                endif

                ! conductive heat flux at snow/soil interface  ! check this!
                keff = 2_r_2*(vsnow(kk)%kth(1)*var(kk,1)%kth)/(vsnow(kk)%kth(1)*dx(kk,1)+var(kk,1)%kth*vsnow(kk)%depth(1))
                qcond(kk,0) = keff*(vsnow(kk)%tsn(1)-Tsoil(kk,1))
                if (vsnow(kk)%hliq(1)>zero) then
                   qcondTa(kk,0) = zero
                else
                   qcondTa(kk,0) = keff
                endif
                qcondTb(kk,0) = -keff
                qh(kk,0) = qcond(kk,0)
                qhTa(kk,0) = qcondTa(kk,0)
                qhTb(kk,0) = qcondTb(kk,0)

                ! advective heat flux at snow/soil interface
                Tqw  = merge(vsnow(kk)%tsn(1), Tsoil(kk,1), q(kk,0)>zero)
                dTqwdTb = merge(zero,one, q(kk,0)>zero)
                dTqwdTa = merge(one,zero, q(kk,0)>zero)
                if (vsnow(kk)%hliq(1)>zero) then
                   qadv(kk,0) = rhow*q(kk,0)*cswat*Tqw
                   qadvTa(kk,0) = zero
                else
                   !qadv(kk,0) = rhow*q(kk,0)*(csice*Tqw-lambdaf)
                   !qadvTa(kk,0) = rhow*csice*q(kk,0)*dTqwdTa  +  rhow*csice*Tqw*qTa(kk,0)
                   qadv(kk,0) = rhow*q(kk,0)*cswat*Tqw
                   qadvTa(kk,0) = rhow*cswat*q(kk,0)*dTqwdTa  +  rhow*cswat*Tqw*qTa(kk,0)
                endif
                qadvTb(kk,0) = rhow*cswat*q(kk,0)*dTqwdTb + rhow*cswat*Tqw*qTb(kk,0)

                if (ns(kk)==0) then
                   ! pond included in top soil layer: calculate phi1(eff)
                   var(kk,1)%phi = max((var(kk,1)%phie -var(kk,1)%he*var(kk,1)%Ksat), &
                        (one+e5)*var(kk,1)%phie)+(h0(kk)-hice(kk))*var(kk,1)%Ksat
                endif

             case default
                write(*,*) "solve: illegal pond case."
                stop
             end select ! pondcase
             ! finished all the surfaces

             ! get moisture fluxes and derivatives (at time t=0, i.e. q0 etc.)

             call getfluxes_vp(n, ns(kk), dx(kk,1:n), vtop(kk), vbot(kk), par(kk,1:n), var(kk,1:n), & ! moisture fluxes
                  hint(kk,1:n), phimin(kk,1:n), q(kk,0:n), qya(kk,0:n), qyb(kk,0:n), qTa(kk,0:n), qTb(kk,0:n), &
                  qliq(kk,0:n), qlya(kk,0:n), qlyb(kk,0:n), qv(kk,0:n), qvT(kk,0:n), qvh(kk,0:n), qvya(kk,0:n), &
                  qvyb(kk,0:n), &
                  iflux(kk), init(kk), getq0(kk), getqn(kk), Tsoil(kk,1:n), T0(kk), nsat(kk), nsatlast(kk))
             qTa(kk,n) = zero
             qTb(kk,n) = zero
             qvTa(kk,1:n) = qTa(kk,1:n)
             qvTb(kk,1:n) = qTb(kk,1:n)
             qlTb(kk,1:n) = zero
             ! get  fluxes heat and derivatives (at time t=0, i.e. q0 etc.)
             !             call getheatfluxes(ns(kk:kk), h0(kk:kk), dx(kk:kk,1:n), dxL(kk:kk), &
             !                  qh(kk:kk,0:n), qhya(kk:kk,0:n), qhyb(kk:kk,0:n), qhTa(kk:kk,0:n), qhTb(kk:kk,0:n), &
             !                  var(kk:kk,1:n), vlit(kk:kk), Tsoil(kk:kk,1:n), TL(kk:kk), T0(kk:kk), litter, &
             !                  q(kk:kk,0:n), qya(kk:kk,0:n), qyb(kk:kk,0:n), qTa(kk:kk,0:n), qTb(kk:kk,0:n), &
             !                  qadv(kk:kk,0:n),qadvya(kk:kk,0:n), qadvyb(kk:kk,0:n), qadvTa(kk:kk,0:n), qadvTb(kk:kk,0:n), &
             !                  advection) ! heat fluxes
             call getheatfluxes(n, ns(kk), h0(kk), dx(kk,1:n), dxL(kk), &
                  qh(kk,0:n), qhya(kk,0:n), qhyb(kk,0:n), qhTa(kk,0:n), qhTb(kk,0:n), &
                  var(kk,1:n), vlit(kk), Tsoil(kk,1:n), TL(kk), T0(kk), litter, &
                  q(kk,0:n), qya(kk,0:n), qyb(kk,0:n), qTa(kk,0:n), qTb(kk,0:n), &
                  qadv(kk,0:n),qadvya(kk,0:n), qadvyb(kk,0:n), qadvTa(kk,0:n), qadvTb(kk,0:n), &
                  advection) ! heat fluxes

             ! get heat and vapour fluxes and derivatives in snow-pack
             !MC to be included with more snow layers !!!
             ! if (vsnow(kk)%nsnow>1) then
             !    ! vapour flux at interface between top snow layer and layer beneath
             !    keff = 2_r_2*((vsnow(kk)%kE(2)/(thousand*lambdaf))*(vsnow(kk)%kE(2)/(thousand*lambdaf))/ &
             !         ((vsnow(kk)%kE(2)/(thousand*lambdaf))*vsnow(kk)%depth(1)+(vsnow(kk)%kE(1)/ &
             !         (thousand*lambdaf))*vsnow(kk)%depth(2)) )
             !    q_snow(kk,1) = keff*(vsnow(kk)%tsn(2)-vsnow(kk)%tsn(1))
             !    qTa_snow(kk,1) = merge(keff,zero,vsnow(kk)%hliq(2)<=zero)
             !    qTb_snow(kk,1) = merge(-keff,zero,vsnow(kk)%hliq(1)<=zero)
             !    qya_snow(kk,1) = zero
             !    qyb_snow(kk,1) = zero
             !    ! conductive heat flux at interface between top snow layer and layer beneath
             !    keff = 2_r_2*(vsnow(kk)%kth(2)*vsnow(kk)%kth(1))/ &
             !         (vsnow(kk)%kth(2)*vsnow(kk)%depth(1)+vsnow(kk)%kth(1)*vsnow(kk)%depth(2))  ! check this!
             !    qcond_snow(kk,1) = keff*(vsnow(kk)%tsn(2)-vsnow(kk)%tsn(2))
             !    qcondTa_snow(kk,1) = merge(zero,keff,vsnow(kk)%hliq(2)>zero)
             !    qcondTb_snow(kk,1) = merge(zero,-keff,vsnow(kk)%hliq(1)>zero)

             !    qh_snow(kk,1) = qcond_snow(kk,1)
             !    qhTa_snow(kk,1) = qcondTa_snow(kk,1)
             !    qhTb_snow(kk,1) = qcondTb_snow(kk,1)

             !    ! advective heat flux at snow/soil interface
             !    Tqw  = merge(vsnow(kk)%tsn(2), vsnow(kk)%tsn(1), q_snow(kk,1)>zero)
             !    dTqwdTb = merge(zero,one, q_snow(kk,1)>zero)
             !    dTqwdTa = merge(one,zero, q_snow(kk,1)>zero)
             !    if (vsnow(kk)%hliq(1)>zero) then
             !       qadv(kk,0) = rhow*q(kk,0)*cswat*Tqw
             !       qadvTa(kk,0) = zero
             !    else
             !       qadv(kk,0) = rhow*q(kk,0)*(csice*Tqw-lambdaf)
             !       qadvTa(kk,0) = rhow*csice*q(kk,0)*dTqwdTa  +  rhow*csice*Tqw*qTa(kk,0)
             !    endif
             !    qadvTb(kk,0) = rhow*cswat*q(kk,0)*dTqwdTb + rhow*cswat*Tqw*qTb(kk,0)
             ! endif

             if (ns(kk)==0) then ! pond included in top soil layer
                ! change qya(1) from dq/dphi (returned by getfluxes) to dq/dh
                qya(kk,1) = var(kk,1)%Ksat*qya(kk,1)
                if (advection==1) then
                   qhya(kk,1) = qhya(kk,1) - qadvya(kk,1)
                   Tqw  = merge(Tsoil(kk,1), Tsoil(kk,2), q(kk,1)>zero)
                   qadvya(kk,1) =  rhow*cswat*qya(kk,1)*Tqw  ! apply corrected qya(kk,1) to qadvya(kk,1)
                   qhya(kk,1) = qhya(kk,1) + qadvya(kk,1)
                endif
             endif

             ! adjust for bottom boundary condition
             if (botbc=="zero flux") then
                qliq(kk,n) = zero
                qv(kk,n)   = zero
                q(kk,n)    = zero
                qya(kk,n)  = zero
                qlya(kk,n) = zero
                qvya(kk,n) = zero
             endif

             ! specify mositure flux at bottom of soil column (heat flux set to zero)
             if (botbc /= "constant head") then
                select case (botbc)
                case ("zero flux")
                   q(kk,n)   = zero
                   qya(kk,n) = zero
                case ("free drainage")
                   q(kk,n) = gf*var(kk,n)%K
                   if (var(kk,n)%isat == 0) then
                      qya(kk,n) = gf*var(kk,n)%KS
                   else
                      qya(kk,n) = zero
                   end if
                case ("seepage")
                   if (var(kk,n)%h <= -half*gf*dx(kk,n)) then
                      q(kk,n)   = zero
                      qya(kk,n) = zero
                   end if
                case default
                   write(*,*) "solve: illegal bottom boundary condition."
                   stop
                end select
             end if
             if (present(qali)) then
                if (qali(kk)>zero) then
                   q(kk,n)   = -qali(kk)
                   qya(kk,n) = zero
                end if
             endif
             if (experiment==7 .or. experiment==8) then
                q(kk,n)   = q(kk,0)
                qya(kk,n) = zero
             endif
             if (experiment==8) qh(kk,n) = G0(kk)



             ! adjust lower heat flux for advection
             if (advection==1) then
                qadv(kk,n) = rhow*cswat*(Tsoil(kk,n))*q(kk,n)
                qadvya(kk,n) = rhow*cswat*(Tsoil(kk,n))*qya(kk,n)
                qadvTa(kk,n)= rhow*cswat*q(kk,n)
                qh(kk,n) = qh(kk,n) + qadv(kk,n)
                qhya(kk,n) = qhya(kk,n) + qadvya(kk,n)
                qhTa(kk,n) = qhTa(kk,n) + qadvTa(kk,n)
             else
                qadv(:,:)   = zero
                qadvya(:,:) = zero
                qadvyb(:,:) = zero
                qadvTa(:,:) = zero
                qadvTb(:,:) = zero
             endif

             qexd(kk,1:n) = zero ! time derivative for root extraction (assumed fixed at input value)
             again(kk)  = .false. ! flag for recalcn of fluxes (default=false)
             !----- end get fluxes and derivs

             !----- first estimate of time step dt before the calculation
             !      gets revised after the calculation
             dmax(kk)     = zero
             tmp2d1(kk,:) = zero
             tmp2d2(kk,:) = zero !  temp storage
             ! estimate rate of change of moisture storage [m/s]
             where (var(kk,1:n)%isat==0) tmp2d1(kk,1:n) = &
                  abs(q(kk,1:n)-q(kk,0:n-1)-iqex(kk,1:n))/(par(kk,1:n)%thre*dx(kk,1:n))
             ! estimate rate of change of temperature [K/s]
             tmp2d2(kk,1:n) = abs(qh(kk,1:n)-qh(kk,0:n-1))/(var(kk,1:n)%csoileff*dx(kk,1:n))
             if (advection==1) then
                tmp2d2(kk,1:n) = abs((qh(kk,1:n)-qadv(kk,1:n))-(qh(kk,0:n-1)-qadv(kk,0:n-1)))/(var(kk,1:n)%csoileff*dx(kk,1:n))
             endif
             if (litter .and. ns(kk)==1 ) then ! litter , no pond
                if (vlit(kk)%isat==0) then ! estimate rate of change of moisture storage [m/s]
                   tmp2d1(kk,0) = abs(q(kk,0) - qL(kk))/(plit(kk)%thre*dxL(kk))
                endif
                tmp2d2(kk,0) = abs(qh(kk,0) - qhL(kk))/(vlit(kk)%csoileff*dxL(kk)) ! estimate rate of change of heat storage [K/s]
                tmp1d3(kk) = (dTLmax-abs(deltaTa(kk))) / tmp2d2(kk,0)
             else
                tmp2d1(kk,0) = zero
                tmp2d2(kk,0) = zero
                tmp1d3(kk)   = dtmax
             endif

             dmax(kk) = maxval(tmp2d1(kk,:),1) ! max derivative |dS/dt|

             if (abs(minval(tmp2d1(kk,:),1)) > maxval(tmp2d1(kk,:),1)) then
                write(*,*) 'Should not be here (01)'
                stop
                dmax(kk) = minval(tmp2d1(kk,:),1)
             endif
             tmp1d1(kk) = maxval(tmp2d2(kk,0:n),1) ! max derivative |dTsoil/dt|
             if (dmax(kk) > zero) then
                dt(kk) = min(dSmax/dmax(kk), dTLmax/tmp1d1(kk), tmp1d3(kk)) ! constrained either by moisture or temp
                ! if pond, overwrite dt
                ! if (h0(kk)>zero .and. (q(kk,1)-qpme(kk))*dt(kk)>h0(kk).and.ns(kk)==0) &
                !      dt(kk) = (h0(kk)-half*h0min)/(q(kk,1)-qpme(kk))
             else ! steady state flow
                if (qpme(kk)>=q(kk,n)) then ! if saturated soil columnn and more precip then drainige -> finish
                   dt(kk) = tfin-t(kk) ! step to finish
                else ! otherwise adjust dt because change of pond height
                   dt(kk) = -(h0(kk)-half*h0min)/(qpme(kk)-q(kk,n))
                end if
                dt(kk) = min(dt(kk), dTLmax/tmp1d1(kk), tmp1d3(kk)) ! constrained by  temp
             end if
             if (dt(kk)>dtmax) dt(kk) = dtmax ! user's limit

             ! if initial step, improve phi where S>=1
             ! might be that you get better derivatives, especially at sat/non-sat interfaces
             if (nsteps(kk)==nsteps0(kk) .and. nsat(kk)>0 .and. iflux(kk)==1) then
                again(kk) = .true.
                dt(kk)    = dtmin
             end if
             ! if fully saturated but was not fully saturated before, adjust even within each time step iteration
             if (nsat(kk)==n .and. nsatlast(kk)<n .and. iflux(kk)==1) then
                again(kk) = .true. ! profile has just become saturated so adjust phi values
                dt(kk)    = dtmin
             end if
             ! sprint to the end
             if (t(kk)+1.1_r_2*dt(kk)>tfin) then ! step to finish
                dt(kk) = tfin-t(kk)
                t(kk)  = tfin
             else
                t(kk) = t(kk)+dt(kk) ! tentative update
                if (again(kk)) t(kk) = t(kk)-dt(kk)
             end if

             !----- end estimate time step dt

             !----- get and solve eqns
             rsigdt(kk) = one/(sig(kk)*dt(kk))
             if (.not. again(kk))  then

                dwoff(kk) = max(h0(kk)*(one-var(kk,1)%thetai/par(kk,1)%thre)-h0max,zero)
                dwoff(kk) = min(dwoff(kk),max((q(kk,0)-qprec_snow(kk))*dt(kk),zero))
                qrunoff(kk) = dwoff(kk)/dt(kk)

             else
                qrunoff(kk) = zero
             endif

             ! aa, bb, cc and dd hold coeffs and rhs of linear equation eqn set
             if (septs == 1) then ! uncoupling of T and S
                qTa(kk,:)   = zero
                qTb(kk,:)   = zero
                qhya(kk,:)  = zero
                qhyb(kk,:)  = zero
                qTbL(kk)    = zero
                qhybL(kk)   = zero
             endif

             iok(kk)            = 0 ! flag for time step test
             itmp(kk)           = 0 ! counter to abort if not getting solution
             again_ice(kk,:)      = .false.
             nsteps_ice(kk,1:n) = 0
             h0_tmp(kk)=h0(kk)
             do while (iok(kk)==0) ! keep reducing time step until all ok
                itmp(kk)  = itmp(kk) + 1
                accel(kk) = one - 0.05_r_2*real(min(10,max(0,itmp(kk)-4)),r_2) ! acceleration [0.5,1], start with 1
                if (itmp(kk) > 1000) then
                   write(*,*) "solve: too many iterations of equation solution"
                   write(*,*) " irec, kk, S"
                   write(*,*) irec, kk, S(kk,:)
                   write(*,*) " irec, kk, Tsoil"
                   write(*,*) irec, kk, Tsoil(kk,:)
                   write(*,*) " irec, kk, qex"
                   write(*,*) irec, kk, iqex(kk,:)
                   write(*,*) " irec, kk, h0, hsnow, hsnowliq"
                   write(*,*) irec, kk, h0(kk), vsnow(kk)%hsnow, vsnow(kk)%hliq
                   stop
                end if

                ! prelim estimate of new snow layer depth for use in energy cons eq'n             !
                if ((vsnow(kk)%nsnow>0))  then
                   do i =1,vsnow(kk)%nsnow
                      hsnow(kk,i) = vsnow(kk)%hsnow(i) +(-q(kk,0) + q_snow(kk,0))*dt(kk)
                   enddo
                endif

                aa(kk,1:n)   =  qya(kk,0:n-1)
                ee(kk,0:n-1) = -qyb(kk,0:n-1)
                bb(kk,1:n)   =  qTa(kk,0:n-1)
                ff(kk,0:n-1) = -qTb(kk,0:n-1)
                gg(kk,1:n) = -(q(kk,0:n-1)-q(kk,1:n)-iqex(kk,1:n))*rsig(kk)
                gg(kk,1) = gg(kk,1)+qrunoff(kk)*rsig(kk)

                aah(kk,1:n)   =  qhya(kk,0:n-1)
                eeh(kk,0:n-1) = -qhyb(kk,0:n-1)
                bbh(kk,1:n)   =  qhTa(kk,0:n-1)
                ffh(kk,0:n-1) = -qhTb(kk,0:n-1)

                ggh(kk,1:n) =  -(qh(kk,0:n-1)-qh(kk,1:n))*rsig(kk) + &
                     rhow*(lambdaf+(cswat-csice)*(Tsoil(kk,1:n)))* &
                     thetai(kk,1:n)*dx(kk,1:n)*rsigdt(kk)*real(imelt(kk,1:n),r_2)

                ggh(kk,1) = ggh(kk,1) + &
                     rhow*(lambdaf+(cswat-csice)*(Tsoil(kk,1)))* &
                     thetai(kk,1)*h0(kk)/par(kk,1)%thre*rsigdt(kk)*real(imelt(kk,1),r_2)

                if (advection==1) then
                   ggh(kk,1:n) = ggh(kk,1:n) + iqex(kk,1:n)*rsig(kk)*(Tsoil(kk,1:n))*cswat*rhow
                   ggh(kk,1) = ggh(kk,1) + qrunoff(kk)*rsig(kk)*(Tsoil(kk,1))*cswat*rhow
                endif

                if (litter) then ! full litter model: litter in zeroth layer
                   ! only use zeroth layer for litter (pond included in layer 1)
                   cc(kk,0) = -qya(kk,0) - rsigdt(kk)*plit(kk)%thre*dxL(kk) + qybL(kk)
                   gg(kk,0)  = -(qprec(kk)-qevap(kk)-q(kk,0))*rsig(kk)
                   ggh(kk,0) = -(G0(kk)-qh(kk,0))*rsig(kk)
                   dd(kk,0)  = -qTa(kk,0)
                endif
                aa(kk,0)  = qya(kk,0)
                aah(kk,0) = zero
                bbh(kk,0) = zero

                where (var(kk,1:n)%isat==0) ! unsaturated layers
                   cc(kk,1:n) = qyb(kk,0:n-1) - qya(kk,1:n) - par(kk,1:n)%thre*dx(kk,1:n)*rsigdt(kk) - qexd(kk,1:n)
                elsewhere ! saturated layers
                   cc(kk,1:n) = qyb(kk,0:n-1) - qya(kk,1:n) - qexd(kk,1:n)
                endwhere

                if (ns(kk)<1) then ! pond included in top soil layer, solving for change in pond height
                   cc(kk,1) = -qya(kk,1)-rsigdt(kk) -qexd(kk,1)
                endif

                cch(kk,1:n) = qhyb(kk,0:n-1)-qhya(kk,1:n) +  real(1-imelt(kk,1:n),r_2)* &
                     real(var(kk,1:n)%iice,r_2)*real(1-var(kk,1:n)%isat,r_2)*rhow*lambdaf*par(kk,1:n)%thre*dx(kk,1:n)*rsigdt(kk)

                if (ns(kk)==0) then ! change in pond height (top layer saturated)
                   cch(kk,1) = cch(kk,1) + real(1-imelt(kk,1),r_2)* &
                        real(var(kk,1)%iice,r_2)*rhow*lambdaf*rsigdt(kk)*(var(kk,1)%thetai/par(kk,1)%thre)
                endif

                ! modification to cch for advection
                if (advection==1) then
                   if (ns(kk)==0) then  ! changing pond included in top soil layer
                      cch(kk,1) = cch(kk,1) -rhow*cswat*(Tsoil(kk,1))*rsigdt(kk)*real(1-var(kk,1)%iice,r_2) &
                           -rhow*csice*(Tsoil(kk,1))*rsigdt(kk)*real(var(kk,1)%iice,r_2) &
                           *(var(kk,1)%thetai/par(kk,1)%thre)*real(1-imelt(kk,1),r_2) &
                           -rhow*cswat*(Tsoil(kk,1))*rsigdt(kk)*real(var(kk,1)%iice,r_2)* &
                           (one-(var(kk,1)%thetai/par(kk,1)%thre))*real(1-imelt(kk,1),r_2) &
                           -rhow*(Tsoil(kk,1))*rsigdt(kk)*real(imelt(kk,1),r_2)*cswat
                   else ! no pond
                      cch(kk,1) = cch(kk,1) -rhow*(Tsoil(kk,1))*par(kk,1)%thre*dx(kk,1)*rsigdt(kk)* &
                           real(1-var(kk,1)%isat,r_2)*(cswat*real(1-var(kk,1)%iice,r_2)+cswat*real(imelt(kk,1),r_2) &
                           +csice*real(1-imelt(kk,1),r_2)*real(var(kk,1)%iice,r_2))
                   endif

                   cch(kk,2:n) = cch(kk,2:n) -rhow*(Tsoil(kk,2:n))*par(kk,2:n)%thre*dx(kk,2:n)*rsigdt(kk)* &
                        real(1-var(kk,2:n)%isat,r_2)*(cswat*real(1-var(kk,2:n)%iice,r_2)+cswat*real(imelt(kk,2:n),r_2) &
                        +csice*real(1-imelt(kk,2:n),r_2)*real(var(kk,2:n)%iice,r_2))
                endif

                dd(kk,1:n)  = qTb(kk,0:n-1)-qTa(kk,1:n)

                ddh(kk,1:n) = qhTb(kk,0:n-1)-qhTa(kk,1:n) - &
                                ! Only apply latent heat component of heat capacity to total deltaT if soil remains frozen
                     var(kk,1:n)%csoileff*dx(kk,1:n)*rsigdt(kk)*real(1-imelt(kk,1:n),r_2) - &
                     var(kk,1:n)%csoil*dx(kk,1:n)*rsigdt(kk)*real(imelt(kk,1:n),r_2) - &
                     (cswat-csice)*dx(kk,1:n)*var(kk,1:n)%dthetaldt*rhow*(Tsoil(kk,1:n))* &
                     rsigdt(kk)*real(1-imelt(kk,1:n),r_2)*real(var(kk,1:n)%iice,r_2)
                ! modification to ddh(kk,1) for pond
                ddh(kk,1) = ddh(kk,1) - cpeff(kk)*h0_tmp(kk)*rsigdt(kk)*real(1-imelt(kk,1),r_2) - &
                     cp(kk)*h0(kk)*rsigdt(kk)*real(imelt(kk,1),r_2) - &
                     (cswat-csice)*h0(kk)/par(kk,1)%thre*var(kk,1)%dthetaldt*rhow*(Tsoil(kk,1))* &
                     rsigdt(kk)*real(1-imelt(kk,1),r_2)*real(var(kk,1)%iice,r_2)

                if (advection==1) then
                   ddh(kk,1:n) = ddh(kk,1:n) - iqex(kk,1:n)*cswat*rhow
                   ddh(kk,1) = ddh(kk,1) - qrunoff(kk)*cswat*rhow
                endif

                ! modification of matrix to incorporate single snow layer
                if (vsnow(kk)%nsnow==1) then
                   cc(kk,0) = -rsigdt(kk)
                   dd(kk,0) =  qTb_snow(kk,0)-qTa(kk,0)
                   ee(kk,0) =   -qyb(kk,0)
                   ff(kk,0) =   -qTb(kk,0)
                   gg(kk,0) = -(q_snow(kk,0)-q(kk,0))*rsig(kk)
                   if (vsnow(kk)%hliq(1)>zero) then ! liquid phase present, solve for change in liq content
                      ddh(kk,0) =  - rhow*((vsnow(kk)%tsn(1))*(cswat-csice)+lambdaf)*rsigdt(kk)
                   else ! solid phase only; solve for change in snow t
                      ddh(kk,0) = qhTb_snow(kk,0) - qhTa(kk,0) - rhow*csice*hsnow(kk,1)*rsigdt(kk)
                   endif
                   cch(kk,0) = qhyb_snow(kk,0)-qhya(kk,0) - rhow*(csice*(vsnow(kk)%tsn(1))-lambdaf)*rsigdt(kk)
                   eeh(kk,0) = -qhyb(kk,0)
                   ffh(kk,0) = -qhTb(kk,0)
                   ggh(kk,0) = -(qh_snow(kk,0)-qh(kk,0))*rsig(kk)
                endif
                if (litter .and. ns(kk)==1) then ! litter and no pond
                   ! watch for deltaTa
                   cc(kk,0)  = qybL(kk) -qya(kk,0) -rsigdt(kk)*plit(kk)%thre*dxL(kk)
                   dd(kk,0)  = qTbL(kk) - qTa(kk,0)
                   ee(kk,0)  = -qyb(kk,0)
                   ff(kk,0)  = -qTb(kk,0)
                   gg(kk,0)  = (q(kk,0) - qL(kk))*rsig(kk) + deltaTa(kk)*(qTbL(kk)-qTa(kk,0))
                   gg(kk,1)  = -(q(kk,0)-q(kk,1)-iqex(kk,1))*rsig(kk) + deltaTa(kk)*qTa(kk,0)
                   cch(kk,0) = qhybL(kk) - qhya(kk,0)
                   ddh(kk,0) = -qhTa(kk,0)-vlit(kk)%csoil*dxL(kk)*rsigdt(kk)+qhTbL(kk)
                   eeh(kk,0) = -qhyb(kk,0)
                   ffh(kk,0) = -qhTb(kk,0)
                   ggh(kk,0) = (qh(kk,0)-qhL(kk))*rsig(kk) + deltaTa(kk)*(qhTbL(kk) - qhTa(kk,0))
                   ggh(kk,1) = -(qh(kk,0)-qh(kk,1))*rsig(kk) + deltaTa(kk)*qhTa(kk,0)
                endif

                ! litter and pond !!vh!! need to check this now that pond is lumped with top soil layer
                if (litter .and. ns(kk)==0) then
                   ddh(kk,0) = -qhTa(kk,0)-cswat*rhow*h0(kk)*rsigdt(kk) -vlit(kk)%csoil*dxL(kk)*rsigdt(kk)
                endif

                if (septs == 1) then ! uncoupled of T and S
                   ! Solve two tridiagonal matrices instead of 1 block matrix
                   ! Could try two different time steps for moisture and temperature
                   if ((ns(kk) == 1 ) .and. (.not. litter)) then ! no pond , no litter
                      nns(kk) = 1
                   else
                      nns(kk) = 0
                   endif
                   !nn(kk) = n-nns(kk) + 1
                   !call dgtsv(nn, 1, aa(nns+1:n), cc(nns:n), ee(nns:n-1), gg(nns:n), nn, info)
                   !dy(nns:n) = gg(nns:n)
                   call tri(nns(kk), n, aa(kk,0:n), cc(kk,0:n), ee(kk,0:n), gg(kk,0:n), dy(kk,0:n))
                   if (nns(kk) == 1) dy(kk,0) = zero
                   !call dgtsv(nn, 1, bbh(nns+1:n), ddh(nns:n), ffh(nns:n-1), ggh(nns:n), nn, info)
                   !dTsoil(nns:n) = ggh(nns:n)
                   call tri(nns(kk), n, bbh(kk,0:n), ddh(kk,0:n), ffh(kk,0:n), gg(kk,0:n), dTsoil(kk,0:n))
                   if (nns(kk) == 1) de(kk,0) = zero
                   if (nns(kk)==0 .and. h0(kk)<e3 .and. (.not. litter)) de(kk,0) = zero
                else ! coupled of T and S
                   nns(kk) = 1  ! pond included in top soil layer
                   if (vsnow(kk)%nsnow>0) then
                      nns(kk) = vsnow(kk)%nsnow-1
                   endif

                   count_sparse = count_sparse + 1
                   call massman_sparse(aa(kk,nns(kk)+1:n), aah(kk,nns(kk)+1:n), bb(kk,nns(kk)+1:n), bbh(kk,nns(kk)+1:n), &
                        cc(kk,nns(kk):n), cch(kk,nns(kk):n), dd(kk,nns(kk):n), ddh(kk,nns(kk):n), ee(kk,nns(kk):n-1), &
                        eeh(kk,nns(kk):n-1), &
                        ff(kk,nns(kk):n-1), ffh(kk,nns(kk):n-1), gg(kk,nns(kk):n), ggh(kk,nns(kk):n), &
                        dy(kk,nns(kk):n), de(kk,nns(kk):n),condition=condition)

                   dTsoil(kk,1:n) = de(kk,1:n)
                   where (vsnow(kk)%hliq(:)>zero)
                      dTsnow(kk,:) = zero
                   elsewhere
                      dTsnow(kk,:) = de(kk,lbound(de,2):0)
                   endwhere

                   ! evaluate soil fluxes at sigma of time step
                   qsig(kk,0)  = q(kk,0)+sig(kk)*qyb(kk,0)*dy(kk,1) + sig(kk)*qya(kk,0)*dy(kk,0) + &
                        sig(kk)*qTb(kk,0)*dTsoil(kk,1) + sig(kk)*qTa(kk,0)*de(kk,0)

                   qhsig(kk,0) = qh(kk,0) + sig(kk)*qhyb(kk,0)*dy(kk,1) + sig(kk)*qhya(kk,0)*dy(kk,0) + &
                        sig(kk)*qhTb(kk,0)*dTsoil(kk,1) + sig(kk)*qhTa(kk,0)*de(kk,0)

                   qadvsig(kk,0) = qadv(kk,0) + sig(kk)*qadvyb(kk,0)*dy(kk,1) &
                        + sig(kk)*qadvya(kk,0)*dy(kk,0) + sig(kk)*qadvTb(kk,0)*dTsoil(kk,1) + sig(kk)*qadvTa(kk,0)*de(kk,0)

                   qsig(kk,1:n-1) = q(kk,1:n-1) + sig(kk)*(qya(kk,1:n-1)*dy(kk,1:n-1)+qyb(kk,1:n-1)*dy(kk,2:n) &
                        +qTa(kk,1:n-1)*dTsoil(kk,1:n-1)+qTb(kk,1:n-1)*dTsoil(kk,2:n))
                   qsig(kk,n)     = q(kk,n) + sig(kk)*(qya(kk,n)*dy(kk,n)+qTa(kk,n)*dTsoil(kk,n))

                   qhsig(kk,1:n-1) = qh(kk,1:n-1) + sig(kk)*(qhya(kk,1:n-1)*dy(kk,1:n-1)+qhyb(kk,1:n-1)*dy(kk,2:n) &
                        +qhTa(kk,1:n-1)*dTsoil(kk,1:n-1)+qhTb(kk,1:n-1)*dTsoil(kk,2:n))
                   qhsig(kk,n) = qh(kk,n) + sig(kk)*(qhya(kk,n)*dy(kk,n)+qhTa(kk,n)*dTsoil(kk,n))

                   qadvsig(kk,1:n-1) = qadv(kk,1:n-1) + sig(kk)*(qadvya(kk,1:n-1)*dy(kk,1:n-1)+qadvyb(kk,1:n-1)*dy(kk,2:n) &
                        +qadvTa(kk,1:n-1)*dTsoil(kk,1:n-1)+qadvTb(kk,1:n-1)*dTsoil(kk,2:n))
                   qadvsig(kk,n) = qadv(kk,n) + sig(kk)*(qadvya(kk,n)*dy(kk,n)+qadvTa(kk,n)*dTsoil(kk,n))

                   !LHS_h(kk,1) = (dx(kk,1)*var(kk,1)%csoileff+h0(kk)*cpeff(kk))*dTsoil(kk,1)/dt(kk) - &
                   !     dx(kk,1)*dy(kk,1)/dt(kk)*par(kk,1)%thre*var(kk,1)%iice*rhow*lambdaf*real(1-var(kk,1)%isat,r_2) - &
                   !     dy(kk,1)/dt(kk)*var(kk,1)%iice*rhow*lambdaf*real(1-ns(kk),r_2)*(var(kk,1)%thetai/par(kk,1)%thre) + &
                   !     var(kk,1)%iice*var(kk,1)%dthetaldT*(Tsoil(kk,1)+Tzero)*rhow*(cswat-csice)*dx(kk,1)*dTsoil(kk,1)/dt(kk) + &
                   !     var(kk,1)%iice*var(kk,1)%dthetaldT*(Tsoil(kk,1)+Tzero)*rhow*(cswat-csice)*h0(kk)/par(kk,1)%thre &
                   !             *dTsoil(kk,1)/dt(kk)

                   !LHS_h(kk,2:n) = (dx(kk,2:n)*var(kk,2:n)%csoileff)*dTsoil(kk,2:n)/dt(kk) - &
                   !     real(var(kk,2:n)%iice,r_2)*real(1-var(kk,2:n)%isat,r_2)*dx(kk,2:n)*dy(kk,2:n)/dt(kk) * &
                   !     par(kk,2:n)%thre*rhow*lambdaf + &
                   !     real(var(kk,2:n)%iice,r_2)*(Tsoil(kk,2:n)+Tzero)*rhow*dx(kk,2:n)*(cswat-csice) * &
                   !     var(kk,2:n)%dthetaldT*dTsoil(kk,2:n)/dt(kk)

                   LHS_h(kk,1) = (dx(kk,1)*var(kk,1)%csoileff+h0_tmp(kk)*cpeff(kk))*dTsoil(kk,1)/dt(kk)* &
                        real(1-imelt(kk,1),r_2) + &
                        (dx(kk,1)*var(kk,1)%csoil+h0(kk)*cp(kk))*dTsoil(kk,1)/dt(kk)*real(imelt(kk,1),r_2) - &
                        dy(kk,1)*par(kk,1)%thre*dx(kk,1)*rhow*lambdaf/dt(kk)*var(kk,1)%iice* &
                        real(1-var(kk,1)%isat,r_2)*real(1-imelt(kk,1),r_2) - &
                        dy(kk,1)*(var(kk,1)%thetai/par(kk,1)%thre)*rhow*lambdaf/dt(kk)*real(1-ns(kk),r_2)* &
                        var(kk,1)%iice*real(1-imelt(kk,1),r_2) + &
                        dx(kk,1)*rhow*(lambdaf+(cswat-csice)*(Tsoil(kk,1)))*thetai(kk,1)/dt(kk)*real(imelt(kk,1),r_2) + &
                        h0(kk)/par(kk,1)%thre*rhow*(lambdaf+(cswat-csice)*(Tsoil(kk,1)))*thetai(kk,1)/dt(kk)* &
                        real(imelt(kk,1),r_2) + &
                        var(kk,1)%dthetaldT*(Tsoil(kk,1))*rhow*(cswat-csice)*dx(kk,1)*dTsoil(kk,1)/dt(kk)* &
                        var(kk,1)%iice*real(1-imelt(kk,1),r_2) + &
                        var(kk,1)%dthetaldT*(Tsoil(kk,1))*rhow*(cswat-csice)*h0(kk)/par(kk,1)%thre* &
                        dTsoil(kk,1)/dt(kk)*var(kk,1)%iice*real(1-imelt(kk,1),r_2)

                   LHS_h(kk,2:n) = (dx(kk,2:n)*var(kk,2:n)%csoileff)*dTsoil(kk,2:n)/dt(kk)*real(1-imelt(kk,2:n),r_2) + &
                        (dx(kk,2:n)*var(kk,2:n)%csoil)*dTsoil(kk,2:n)/dt(kk)*real(imelt(kk,2:n),r_2) - &
                        dy(kk,2:n)*par(kk,2:n)%thre*dx(kk,2:n)*rhow*lambdaf/dt(kk)*var(kk,2:n)%iice* &
                        real(1-var(kk,2:n)%isat,r_2)*real(1-imelt(kk,2:n),r_2) + &
                        dx(kk,2:n)*rhow*(lambdaf+(cswat-csice)*(Tsoil(kk,2:n)))*thetai(kk,2:n)/dt(kk)* &
                        real(imelt(kk,2:n),r_2) + &
                        var(kk,2:n)%dthetaldT*(Tsoil(kk,2:n))*rhow*(cswat-csice)*dx(kk,2:n)* &
                        dTsoil(kk,2:n)/dt(kk)*var(kk,2:n)%iice*real(1-imelt(kk,2:n),r_2)

                   RHS_h(kk,1:n) = qhsig(kk,0:n-1) -qhsig(kk,1:n)

                   if (advection==1) then
                      LHS_h(kk,1:n) = LHS_h(kk,1:n) &
                           + real(1-var(kk,1:n)%iice,r_2) * real(1-var(kk,1:n)%isat,r_2) * &
                           dx(kk,1:n) * dy(kk,1:n)/dt(kk) * par(kk,1:n)%thre * (Tsoil(kk,1:n)) * rhow * cswat &
                           + real(imelt(kk,1:n),r_2) * real(1-var(kk,1:n)%isat,r_2) * &
                           dx(kk,1:n) * dy(kk,1:n)/dt(kk) * par(kk,1:n)%thre * (Tsoil(kk,1:n)) * rhow * cswat &
                           + real(var(kk,1:n)%iice,r_2) * real(1-imelt(kk,1:n),r_2) * real(1-var(kk,1:n)%isat,r_2) * &
                           dx(kk,1:n) * dy(kk,1:n)/dt(kk) * par(kk,1:n)%thre * rhow * csice * (Tsoil(kk,1:n))

                      LHS_h(kk,1) = LHS_h(kk,1) &
                           + real(1-ns(kk),r_2) * real(1-var(kk,1)%iice,r_2) * &
                           dy(kk,1)/dt(kk) * rhow * cswat * (Tsoil(kk,1)) &
                           + real(1-ns(kk),r_2) * real(imelt(kk,1),r_2) * &
                           dy(kk,1)/dt(kk) * rhow * cswat * (Tsoil(kk,1)) &
                           + real(1-ns(kk),r_2) * real(var(kk,1)%iice,r_2) * real(1-imelt(kk,1),r_2) * &
                           dy(kk,1)/dt(kk) * rhow * csice * (Tsoil(kk,1)) * (var(kk,1)%thetai/par(kk,1)%thre) &
                           + real(1-ns(kk),r_2) * real(var(kk,1)%iice,r_2) * real(1-imelt(kk,1),r_2) * &
                           dy(kk,1)/dt(kk) * rhow * cswat * (Tsoil(kk,1)) * (one-(var(kk,1)%thetai/par(kk,1)%thre))

                      RHS_h(kk,1:n) = RHS_h(kk,1:n) - iqex(kk,1:n)*cswat*rhow*(Tsoil(kk,1:n)+ sig(kk)*dTsoil(kk,1:n))
                      RHS_h(kk,1) = RHS_h(kk,1) - qrunoff(kk)*cswat*rhow*(Tsoil(kk,1) + sig(kk)*dTsoil(kk,1))
                   endif

                   tmp2d1(kk,1:n) = (LHS_h(kk,1:n)-RHS_h(kk,1:n)) !*dt(kk)
                   if (any(abs(tmp2d1(kk,1:n))>1.e-6_r_2)) then
                      write(*,*) "energy imbalance after matrix solution"
                      write(*,*) "LHS: ", LHS_h(kk,1:n)
                      write(*,*) "RHS: ", RHS_h(kk,1:n)
                      write(*,*) "LHS-RHS: ", LHS_h(kk,1:n)-RHS_h(kk,1:n)
                      stop
                   endif

                   ! check mass balance on top soil layer
                   if (ns(kk)==0) then  ! pond included in top soil layer
                      LHS(kk,1) = dy(kk,1)/dt(kk)
                   else
                      LHS(kk,1) = dy(kk,1)*par(kk,1)%thre*dx(kk,1)*real(-var(kk,1)%isat+1)/dt(kk)
                   endif
                   RHS(kk,1) = qsig(kk,0) - qsig(kk,1)- iqex(kk,1) - qrunoff(kk)

                   ! snow layer
                   if (vsnow(kk)%nsnow>0) then
                      qsig_snow(kk,0) = q_snow(kk,0) + sig(kk)*qyb_snow(kk,0)*dy(kk,0) + &
                           sig(kk)*qTb_snow(kk,0)*de(kk,0)

                      qhsig_snow(kk,0) = qh_snow(kk,0) + sig(kk)*qhyb_snow(kk,0)*dy(kk,0) + &
                           sig(kk)*qhTb_snow(kk,0)*de(kk,0)

                      qadvsig_snow(kk,0) = qadv_snow(kk,0) + sig(kk)*qadvyb_snow(kk,0)*dy(kk,0) + &
                           sig(kk)*qadvTb_snow(kk,0)*de(kk,0)

                      RHS_snow(kk,1) = qsig_snow(kk,0) - qsig(kk,0)         ! rate of change of snow water col
                      LHS_snow(kk,1) = dy(kk,0)/dt(kk)

                      RHS_h_snow(kk,1) = qhsig_snow(kk,0) - qhsig(kk,0)
                      if (vsnow(kk)%hliq(1)>zero) then ! liquid phase present, solve for change in liq content
                         LHS_h_snow(kk,1) = (dy(kk,0)*rhow*(csice*(vsnow(kk)%tsn(1))-lambdaf) + &
                              de(kk,0)*rhow*((cswat-csice)*(vsnow(kk)%tsn(1))+lambdaf))/dt(kk)
                      else
                         LHS_h_snow(kk,1) = (dy(kk,0)*rhow*(csice*(vsnow(kk)%tsn(1))-lambdaf) + &
                              de(kk,0)*csice*rhow*hsnow(kk,1))/dt(kk)
                      endif
                   endif ! snow layer

                endif ! septs==1

                tmp1d1 = nless
                ! dy contains dS or, for sat layers, dphi values
                iok(kk) = 1
                fac(kk) = one

                if (.not. again(kk)) then

                   ! check if time step ok, if not then set fac to make it less
                   iok(kk) = 1
                   do i=1, n
                      if (var(kk,i)%isat==0) then ! check change in S in initially unsaturated layers
                         if (abs(dy(kk,i)) > dSfac*dSmax) then
                            fac(kk) = max(half,accel(kk)*abs(dSmax/dy(kk,i)))
                            nfac1(kk) = nfac1(kk)+1
                            iok(kk) = 0
                            exit
                         end if
                         if (-dy(kk,i) > dSmaxr*S(kk,i)) then ! Check relative moisture change
                            fac(kk) = max(half,accel(kk)*dSmaxr*S(kk,i)/(-dSfac*dy(kk,i)))
                            nfac2(kk) = nfac2(kk)+1
                            iok(kk) = 0
                            exit
                         end if
                         if (S(kk,i) < one .and. S(kk,i)+dy(kk,i) > Smax) then ! Check for oversaturating
                            fac(kk) = accel(kk)*(half*(one+Smax)-S(kk,i))/dy(kk,i)
                            nfac3(kk) = nfac3(kk)+1
                            !write(*,*) 'incrementing nfac3 due to layer', i
                            iok(kk) = 0
                            exit
                         end if
                         if (S(kk,i) >= one .and. dy(kk,i) > half*(Smax-one)) then ! Check for limit at oversaturation
                            fac(kk) = 0.25_r_2*(Smax-one)/dy(kk,i)
                            nfac4(kk) = nfac4(kk)+1
                            iok(kk) = 0
                            exit
                         end if
                      end if
                   end do

                   ! Check absolute soil temperature change in frozen soil layers where updated Tsoil exceeds freezing point
                   if (iok(kk)==1) then
                      do i=1,n
                         if(var(kk,i)%iice==1) then
                            Tfreezing(kk) = Tfrz(S(kk,i)+dy(kk,i)*real(1-var(kk,i)%isat),par(kk,i)%he,one/par(kk,i)%lam)

                            if ((Tsoil(kk,i)+dTsoil(kk,i))>Tfreezing(kk)) then

                               theta = (S(kk,i)+dy(kk,i)*real(1-var(kk,i)%isat))*(par(kk,i)%thre) + &
                                    (par(kk,i)%the - par(kk,i)%thre)
                               deltaJ_latent_T(kk,i) = thetai(kk,i)*dx(kk,i)*rhow*lambdaf
                               tmp1d1(kk) = Tsoil(kk,i)
                               deltaJ_sensible_S(kk,i) = (tmp1d1(kk))*(rhow*cswat*(theta*dx(kk,i)) + &
                                    par(kk,i)%rho*par(kk,i)%css*dx(kk,i)) - &
                                    (tmp1d1(kk))*(rhow*cswat*(var(kk,i)%thetal*dx(kk,i)) + &
                                    par(kk,i)%rho*par(kk,i)%css*dx(kk,i)) - &
                                    (tmp1d1(kk))*(rhow*csice*(thetai(kk,i)*dx(kk,i)))

                               if ((i==1).and.(h0(kk)>zero)) then
                                  deltaJ_latent_T(kk,i) = deltaJ_latent_T(kk,i) + &
                                       thetai(kk,i)*h0(kk)/par(kk,i)%thre*rhow*lambdaf
                                  deltaJ_sensible_S(kk,i) = deltaJ_sensible_S(kk,i) + &
                                       (tmp1d1(kk))*(rhow*cswat*(h0(kk)+dy(kk,1)*real(1-ns(kk)))) - &
                                       (tmp1d1(kk))*(rhow*cswat*(h0(kk)* &
                                       (one-thetai(kk,1)/par(kk,1)%thre))) - &
                                       (tmp1d1(kk))*(rhow*csice*(h0(kk)* &
                                       thetai(kk,1)/par(kk,1)%thre))

                                  tmp1d1(kk) = (LHS_h(kk,i)*dt(kk) - (deltaJ_latent_T(kk,i) + deltaJ_sensible_S(kk,i)))/ &
                                       (rhow*cswat*(theta*dx(kk,1)+(h0(kk)+dy(kk,1) * &
                                       real(1-ns(kk))))+par(kk,1)%rho*par(kk,1)%css*dx(kk,1))
                               else
                                  tmp1d1(kk) = (LHS_h(kk,i)*dt(kk) - (deltaJ_latent_T(kk,i) + deltaJ_sensible_S(kk,i)))/ &
                                       (dx(kk,i)*(cswat*theta*rhow + par(kk,i)%css*par(kk,i)%rho))
                               endif
                               if (tmp1d1(kk) > dTsoilmax) then ! Check absolute soil temperature change
                                  fac(kk) = min(fac(kk),max(half,accel(kk)*abs(dTsoilmax/tmp1d1(kk))))
                                  iok(kk) = 0
                                  exit
                               end if
                            endif
                         endif
                      enddo
                   endif

                   do i=1, n
                      if (abs(dTsoil(kk,i)) > dTsoilmax) then ! Check absolute soil temperature change
                         fac(kk) = min(fac(kk),max(half,accel(kk)*abs(dTsoilmax/dTsoil(kk,i))))
                         nfac5(kk) = nfac5(kk)+1

                         iok(kk) = 0
                         exit
                      end if
                   enddo

                   if (vsnow(kk)%nsnow>0) then
                      do i=1, vsnow(kk)%nsnow
                         if ((vsnow(kk)%hsnow(i) + dy(kk,-i+1))<zero) then
                            fac(kk) = 0.99*vsnow(kk)%hsnow(i)/abs(dy(kk,-i+1))
                            nfac6(kk) = nfac6(kk)+1
                            iok(kk) =0
                            exit
                         endif

                         if (vsnow(kk)%hliq(i)>zero) then
                            if ((vsnow(kk)%hliq(i) + de(kk,-i+1))>(dy(kk,-i+1) + vsnow(kk)%hsnow(i))) then
                               fac(kk) = 0.5
                               nfac7(kk) = nfac7(kk)+1
                               iok(kk) =0
                               exit
                            endif

                            if ((vsnow(kk)%hliq(i) + de(kk,-i+1))<-0.01*vsnow(kk)%hsnow(i)) then
                               fac(kk) = (-0.009*vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))/de(kk,-i+1)
                               nfac8(kk) = nfac8(kk)+1
                               iok(kk) =0
                               exit
                            endif
                         endif
                      enddo
                   endif

                   if (litter .and. ns(kk)==0) then ! litter and pond
                      if (abs(de(kk,0)) > dTLmax) then ! Check absolute litter temperature change
                         fac(kk) = min(fac(kk),max(half,accel(kk)*abs(dTLmax/de(kk,0))))
                         iok(kk) = 0
                      end if
                   endif

                   if (litter .and. ns(kk)==1) then ! litter, no pond
                      if (-dy(kk,0) > SL(kk)*dSmaxr) then ! Check relative litter moisture change
                         fac(kk) = min(fac(kk),max(half,accel(kk)*dSmaxr*SL(kk)/(-dSfac*dy(kk,0))))
                         iok(kk) = 0
                      endif

                      if (abs(deltaTa(kk)) > zero) then
                         tmp1d1(kk) = de(kk,0) - deltaTa(kk)
                         if (abs(tmp1d1(kk)) > dTLmax) then ! Check absolute litter temperature change
                            fac(kk) = min(fac(kk),max(half,accel(kk)*dTLmax/abs(tmp1d1(kk))))
                            iok(kk) = 0
                            if (itmp(kk) > 10) then
                               iok(kk) = 1
                               t(kk)   = t(kk)-dt(kk)
                               n_noconverge(kk) = n_noconverge(kk) + 1
                               goto 10
                            endif
                         endif
                      else
                         if (abs(de(kk,0)) > dTLmax) then ! Check absolute litter temperature change
                            fac(kk) = min(fac(kk),max(half,accel(kk)*abs(dTLmax/de(kk,0))))
                            iok(kk) = 0
                         end if
                      endif
                   endif ! litter, no pond

                   ! pond decreasing or runoff starts
                   if (iok(kk)==1 .and. ns(kk)<1 .and. h0(kk)+dy(kk,1)<h0min) then ! pond going
                      fac(kk) = -(h0(kk)-half*h0min)/dy(kk,1)
                      nfac9(kk) = nfac9(kk) +1
                      iok(kk) = 0
                   end if

                   ! J0+deltaJ<J(thetai=theta) i.e. too much energy extracted
                   J0(kk,1) = rhow*cswat*(Tsoil(kk,1))*dx(kk,1)*var(kk,1)%thetal + &
                        rhow*dx(kk,1)*thetai(kk,1)*(csice*(Tsoil(kk,1))-lambdaf) + &
                        par(kk,1)%rho*par(kk,1)%css*dx(kk,1)*(Tsoil(kk,1)) + &
                        (h0(kk)-hice(kk))*cswat*rhow*(Tsoil(kk,1)) + &
                        hice(kk)*rhow*(csice*(Tsoil(kk,1))-lambdaf)
                   tmp1d2(kk) = J0(kk,1) + LHS_h(kk,1)*dt(kk)
                   tmp1d3(kk) = Tsoil(kk,1)+dTsoil(kk,1)
                   tmp1d4(kk) = h0(kk) + dy(kk,1)*real(one-ns(kk))
                   theta      = (S(kk,1)+ dy(kk,1)*real(one-var(kk,1)%isat))*(par(kk,1)%thre) + (par(kk,1)%the - par(kk,1)%thre)
                   c2 = rhow*((tmp1d3(kk))*csice-lambdaf)*(dx(kk,1)*theta+tmp1d4(kk)*theta/par(kk,1)%thre) + &
                        rhow*(tmp1d3(kk))*cswat*tmp1d4(kk)*(one-theta/par(kk,1)%thre) + &
                        dx(kk,1)*(tmp1d3(kk))*par(kk,1)%rho*par(kk,1)%css
                   if(((tmp1d2(kk))-c2)<zero) then
                      fac(kk) = 0.5
                      iok(kk) = 0
                      nfac10(kk) = nfac10(kk)+1
                   endif

                   if (fac(kk) < one) then
                      again_ice(kk,1:n) = .false.  ! reset all again_ice if calc is to be repeated with smaller time-step
                      imelt(kk,:) = 0
                      var(kk,:)%thetai = thetai(kk,:)
                      var(kk,:)%dthetaldT = dthetaldT(kk,:)
                      do i=1,n
                         theta         = S(kk,i)*(par(kk,i)%thre) + (par(kk,i)%the - par(kk,i)%thre)
                         var(kk,i)%csoil = par(kk,i)%css*par(kk,i)%rho + rhow*cswat*(theta-var(kk,i)%thetai) + &
                              rhow*csice*var(kk,i)%thetai
                         if ((i==1) .and. (ns(kk)==0)) then
                            cp(kk) = real(1-var(kk,1)%iice,r_2)*cswat*rhow & ! heat capacity of pond
                                 + real(var(kk,1)%iice,r_2)*rhow* &
                                 ((one-var(kk,1)%thetai/par(kk,1)%thre)*cswat + (var(kk,1)%thetai/par(kk,1)%thre)*csice)
                         endif
                      enddo
                      var(kk,:)%csoileff = var(kk,:)%csoil + rhow*lambdaf*var(kk,:)%dthetaldT*real(var(kk,:)%iice,r_2)
                      cpeff(kk)= cp(kk) + rhow*lambdaf*var(kk,1)%dthetaldT/par(kk,1)%thre*real(var(kk,1)%iice,r_2)
                      h0_tmp(kk) = h0(kk)

                   endif

                   ! reduce time step
                   if (iok(kk)==0) then
                      t(kk)      = t(kk)-dt(kk)
                      ! dt(kk)     = max(fac(kk)*dt(kk), dtmin)
                      dt(kk)     = fac(kk)*dt(kk)
                      t(kk)      = t(kk)+dt(kk)
                      rsigdt(kk) = one/(sig(kk)*dt(kk))
                      nless(kk)  = nless(kk) + 1 ! count step size reductions
                      !temporary update to h0 for use in energy cons eqn
                      !h0_tmp(kk) = h0(kk) + (qsig(kk,0)-qsig(kk,1))*dt(kk)

                      !write(339,"(i8,i8,i8,14f15.6)") irec, nsteps(kk),nless(kk), dt(kk), fac
                   end if
                   if (var(kk,1)%isat/=0 .and. iflux(kk)==1 .and. var(kk,1)%phi<phip(kk) .and. &
                        var(kk,1)%phi+dy(kk,1)>phip(kk)) then
                      ! incipient (onset of) ponding - adjust state of saturated regions
                      t(kk)      = t(kk)-dt(kk)
                      dt(kk)     = dtmin
                      rsigdt(kk) = one/(sig(kk)*dt(kk))
                      again(kk)  = .true.
                      iok(kk)    = 0
                   end if
                end if
                nsteps(kk)        = nsteps(kk) + 1
                !                if (irec.eq.10) then
                !                    write(*,*) 'writing diags'
                !
                !                     write(345,"(11i8,1500e16.6)") nsteps, nfac1, nfac2, nfac3, nfac4, nfac5, nfac6, nfac7, nfac8, nfac9, nfac10, &
                !                     q, qsig, qH, qhsig, dy(kk,1:n), de(kk,1:n), dTsoil, S,thetai, Tsoil, real(var(kk,1)%iice), &
                !                     real(var(kk,1)%isat), &
                !                      h0(kk), real(iok(kk)), var(kk,1)%phie, var(kk,1)%phi, phip(kk)
                !
                !                     if (nsteps(kk).gt.100) STOP
                !
                !                 endif
             end do ! while (iok==0) ----- end get and solve eqns

             !----- update unknowns
             ! i.e. update state variables to end of time step
             !      cumulate some surface fluxes
             ih0(kk) = 0
             if (.not. again(kk)) then
                ! get surface fluxes at sigma of the time step
                dwoff(kk)  = zero
                deltah0(kk) = zero
                precip(kk) = precip(kk) + (qprec(kk)+qprec_snow(kk))*dt(kk)
                dwcol(kk) = sum(par(kk,1:n)%thre*dx(kk,1:n)*dy(kk,1:n)*(abs(var(kk,1:n)%isat-1)),1) + &
                     real(1-ns(kk),r_2)*dy(kk,1)

                if (vsnow(kk)%nsnow==1) dwcol(kk) = dwcol(kk) + dy(kk,0)

                ! absolute heat stored in soil column before update
                do j=1,n
                   if (j==1) then
                      J0(kk,j) = rhow*cswat*(Tsoil(kk,j))*dx(kk,j)*var(kk,j)%thetal + &
                           rhow*dx(kk,j)*thetai(kk,j)*(csice*(Tsoil(kk,j))-lambdaf) + &
                           par(kk,j)%rho*par(kk,j)%css*dx(kk,j)*(Tsoil(kk,j)) + &
                           (h0(kk)-hice(kk))*cswat*rhow*(Tsoil(kk,1)) + &
                           hice(kk)*rhow*(csice*(Tsoil(kk,1))-lambdaf)
                   else ! (j>1)
                      J0(kk,j) = rhow*cswat*(Tsoil(kk,j))*dx(kk,j)*var(kk,j)%thetal + &
                           rhow*dx(kk,j)*thetai(kk,j)*(csice*(Tsoil(kk,j))-lambdaf) + &
                           par(kk,j)%rho*par(kk,j)%css*dx(kk,j)*(Tsoil(kk,j))
                   endif
                end do ! absolute heat stored in soil column

                ! change in heat stored in soil column
                do j=1,n
                   if (j==1) then
                      ! pond included in top soil layer
                      deltaJ_latent_S(kk,j) = -dx(kk,j)*dy(kk,j)*par(kk,j)%thre*real(var(kk,j)%iice,r_2)* &
                           real(1-var(kk,j)%isat,r_2)*rhow*lambdaf - &
                           real(1-ns(kk),r_2)*dy(kk,1)*real(var(kk,j)%iice,r_2)*rhow*lambdaf*(var(kk,1)%thetai/par(kk,1)%thre)

                      deltaJ_latent_T(kk,j) = var(kk,j)%dthetaldT*dx(kk,j)*dTsoil(kk,j)*rhow*lambdaf* &
                           real(var(kk,j)%iice,r_2) + &
                           h0(kk)/par(kk,1)%thre*var(kk,j)%dthetaldT*dTsoil(kk,j)*rhow*lambdaf*real(var(kk,j)%iice,r_2)

                      deltaJ_sensible_T(kk,j) = var(kk,j)%csoil*dx(kk,j)*dTsoil(kk,j) + &
                           cp(kk)*h0(kk)*dTsoil(kk,j) + &
                           var(kk,1)%iice*var(kk,j)%dthetaldT*(Tsoil(kk,j))*rhow*(cswat-csice)*dx(kk,j)*dTsoil(kk,j) + &
                           var(kk,j)%iice*var(kk,j)%dthetaldT*(Tsoil(kk,j))*rhow*(cswat-csice)* &
                           h0(kk)/par(kk,j)%thre*dTsoil(kk,j)

                      if (advection==1) then
                         deltaJ_sensible_S(kk,j) = rhow*cswat*(Tsoil(kk,j))*dx(kk,j)*par(kk,j)%thre*dy(kk,j)* &
                              real(1-var(kk,j)%isat,r_2)*real(1-var(kk,j)%iice,r_2)  + &
                              rhow*csice*(Tsoil(kk,j))*dx(kk,j)*par(kk,j)%thre*dy(kk,j)* &
                              real(1-var(kk,j)%isat,r_2)*real(var(kk,j)%iice,r_2) + &
                              real(1-ns(kk),r_2)*real(1-var(kk,1)%iice,r_2)* &
                              dy(kk,1)*rhow*cswat*(Tsoil(kk,1))  + &
                              real(1-ns(kk),r_2)*real(var(kk,1)%iice,r_2)* &
                              dy(kk,1)*rhow*csice*(Tsoil(kk,1))*(var(kk,1)%thetai/par(kk,1)%thre) + &
                              real(1-ns(kk),r_2)*real(var(kk,1)%iice,r_2)* &
                              dy(kk,1)*rhow*cswat*(Tsoil(kk,1))*(one-(var(kk,1)%thetai/par(kk,1)%thre))
                      endif
                   else ! (j>1)
                      deltaJ_latent_S(kk,j) = -dx(kk,j)*dy(kk,j)*par(kk,j)%thre*real(var(kk,j)%iice,r_2)* &
                           real(1-var(kk,j)%isat,r_2)*rhow*lambdaf
                      deltaJ_sensible_T(kk,j) = var(kk,j)%csoil*dx(kk,j)*dTsoil(kk,j) + &
                           var(kk,1)%iice*var(kk,j)%dthetaldT*(Tsoil(kk,j))*rhow*(cswat-csice)*dx(kk,j)*dTsoil(kk,j)

                      if (advection==1) then
                         deltaJ_sensible_S(kk,j) = rhow*cswat*(Tsoil(kk,j))*dx(kk,j)*dy(kk,j)*par(kk,j)%thre* &
                              real(1-var(kk,j)%isat,r_2)*real(1-var(kk,j)%iice,r_2) &
                              + rhow*csice*(Tsoil(kk,j))*dx(kk,j)*dy(kk,j)*par(kk,j)%thre &
                              *real(1-var(kk,j)%isat,r_2)*real(var(kk,j)%iice,r_2)
                      endif

                      deltaJ_latent_T(kk,j) = var(kk,j)%dthetaldT*dx(kk,j)*dTsoil(kk,j)*rhow*lambdaf*real(var(kk,j)%iice,r_2)
                   endif
                end do

                if (ns(kk)<1) then ! change in pond height
                   h0(kk) = h0(kk) + dy(kk,1)
                   deltah0(kk) = dy(kk,1)
                   if (h0(kk)<zero .and. dy(kk,1)<zero) then
                      ih0(kk)=1 ! pond gone
                   endif
                endif

                ! cumulate evaporation from top of soil column or top of litter/pond or top of snow pack
                select case (pondcase(kk))
                case(1)  ! no pond or litter
                   evap(kk)     = evap(kk) +qevap(kk)*dt(kk) - sig(kk)*(qyb(kk,0)*dy(kk,1)+qTb(kk,0)*DTsoil(kk,1))*dt(kk)
                   qevapsig(kk) = qevap(kk) - sig(kk)*(qyb(kk,0)*dy(kk,1)+qTb(kk,0)*DTsoil(kk,1))
                   Gcum(kk) = Gcum(kk)+(qhsig(kk,0)-qadvsig(kk,0))*dt(kk)
                   Qadvcum(kk) = Qadvcum(kk) + qadvsig(kk,0)*dt(kk) - qadvsig(kk,n)*dt(kk)
                   lEcum(kk)    = lEcum(kk) + qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk)
                   Hcum(kk) = Hcum(kk) + (vmet(kk)%Rn*dt(kk)-(qhsig(kk,0)-qadvsig(kk,0))*dt(kk)- &
                        qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk))
                   dwinfil(kk) = (qprec(kk)+qprec_snow(kk)-qevapsig(kk))*dt(kk)
                case(2)  ! litter, no pond
                   evap(kk)     = evap(kk) - (qL(kk)-qprec(kk) + sig(kk)*(qybL(kk)*dy(kk,0)+qTbL(kk)*(de(kk,0) &
                        -deltaTa(kk))))*dt(kk)
                   qevapsig(kk) = - (qL(kk)-qprec(kk) + sig(kk)*(qybL(kk)*dy(kk,0)+qTbL(kk)*(de(kk,0)-deltaTa(kk))))
                   Gcum(kk)  = Gcum(kk)+qhsig(kk,0)*dt(kk)
                   lEcum(kk) = lEcum(kk) + qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk)
                   Hcum(kk)  = Hcum(kk) + (vmet(kk)%Rn*dt(kk)-qhsig(kk,0)*dt(kk) - &
                        qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk))
                   dwinlit(kk)    = (qL(kk) + sig(kk)*(qybL(kk)*dy(kk,0)+qTbL(kk)*(de(kk,0)-deltaTa(kk))))*dt(kk)
                   dwinfil(kk) = qsig(kk,0)*dt(kk)
                case(3)! pond
                   evap(kk)     = evap(kk) +qevap(kk)*dt(kk) - sig(kk)*(qyb(kk,0)*dy(kk,1)+qTb(kk,0)*DTsoil(kk,1))*dt(kk)
                   qevapsig(kk) = qevap(kk) - sig(kk)*(qyb(kk,0)*dy(kk,1)+qTb(kk,0)*DTsoil(kk,1))
                   dwinfil(kk)  = (q(kk,0)+sig(kk)*(qya(kk,0)*dy(kk,0)+qyb(kk,0)*dy(kk,1)+qTb(kk,0)*dTsoil(kk,1)))*dt(kk) - &
                        deltah0(kk) ! pond included in top soil layer

                   Gcum(kk) = Gcum(kk)+(qhsig(kk,0)-qadvsig(kk,0))*dt(kk)
                   Qadvcum(kk) = Qadvcum(kk) + qadvsig(kk,0)*dt(kk) - qadvsig(kk,n)*dt(kk)
                   lEcum(kk)    = lEcum(kk) + qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk)
                   Hcum(kk) = Hcum(kk) + (vmet(kk)%Rn*dt(kk)-(qhsig(kk,0)-qadvsig(kk,0))*dt(kk)- &
                        qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk))
                case(4) ! dedicated snow layer
                   evap(kk)     = evap(kk) +qevap(kk)*dt(kk) - sig(kk)*(qyb_snow(kk,0)*dy(kk,0)+qTb_snow(kk,0)*de(kk,0))*dt(kk)
                   qevapsig(kk) = qevap(kk) - sig(kk)*(qyb_snow(kk,0)*dy(kk,0)+qTb_snow(kk,0)*de(kk,0))
                   Gcum(kk) = Gcum(kk)+(qhsig_snow(kk,0)-qadvsig_snow(kk,0))*dt(kk)
                   Qadvcum(kk) = Qadvcum(kk) + qadvsig_snow(kk,0)*dt(kk) - qadvsig(kk,n)*dt(kk)
                   lEcum(kk)    = lEcum(kk) + qevapsig(kk)*thousand*var(kk,1)%lambdaf*dt(kk)
                   Hcum(kk) = Hcum(kk) + (vmet(kk)%Rn*dt(kk)-(qhsig_snow(kk,0)-qadvsig_snow(kk,0))*dt(kk)- &
                        qevapsig(kk)*thousand*var(kk,1)%lambdaf*dt(kk))
                   dwinfil(kk) = (qprec(kk)+qprec_snow(kk)-qevapsig(kk))*dt(kk)
                end select ! pondcase

                dwdrainage(kk)     = q(kk,n)*dt(kk) +sig(kk)*dt(kk)*(qya(kk,n)*dy(kk,n)+qTa(kk,n)*dTsoil(kk,n))
                if (botbc=="aquifer" .and. v_aquifer(kk)%isat==0) then
                   dwdischarge(kk) = v_aquifer(kk)%discharge*dt(kk)
                else
                   dwdischarge(kk) = dwdrainage(kk)
                endif

                drexcol(kk) = sum(iqex(kk,:),1)*dt(kk)

                ! cumulative fluxes
                infil(kk)     = infil(kk) + dwinfil(kk)
                inlit(kk)     = inlit(kk) + dwinlit(kk)
                wcol(kk)      = wcol(kk) + dwcol(kk)
                drainage(kk)  = drainage(kk) + dwdrainage(kk)
                discharge(kk) = discharge(kk) + dwdischarge(kk)
                rexcol(kk)    = rexcol(kk) + drexcol(kk)
                Qadvcum(kk) = Qadvcum(kk) - sum(iqex(kk,:)*(Tsoil(kk,:)),1)*dt(kk)*rhow*cswat - &
                     qrunoff(kk)*(Tsoil(kk,1))*dt(kk)*rhow*cswat

                ! evaluate soil fluxes at sigma of time step for use in isotope routine

                qsig(kk,0)  = q(kk,0)  + sig(kk)*(qyb(kk,0)*dy(kk,1)  + qya(kk,0)*dy(kk,0)  + qTb(kk,0)*dTsoil(kk,1) &
                     + qTa(kk,0)*de(kk,0))
                qhsig(kk,0) = qh(kk,0) + sig(kk)*(qhyb(kk,0)*dy(kk,1) + qhya(kk,0)*dy(kk,0) + qhTb(kk,0)*dTsoil(kk,1) &
                     + qhTa(kk,0)*de(kk,0))
                qvsig(kk,0) = qv(kk,0)+sig(kk)*qvyb(kk,0)*dy(kk,1)+sig(kk)*qvTb(kk,0)*dTsoil(kk,1)
                qlsig(kk,0) = qliq(kk,0)+sig(kk)*qlyb(kk,0)*dy(kk,1)+sig(kk)*qlTb(kk,0)*dTsoil(kk,1)

                qsig(kk,1:n-1)   = q(kk,1:n-1) + sig(kk)*(qya(kk,1:n-1)*dy(kk,1:n-1) + qyb(kk,1:n-1)*dy(kk,2:n) &
                     + qTa(kk,1:n-1)*dTsoil(kk,1:n-1) + qTb(kk,1:n-1)*dTsoil(kk,2:n))
                qsig(kk,n)       = q(kk,n) + sig(kk)*(qya(kk,n)*dy(kk,n) + qTa(kk,n)*dTsoil(kk,n))
                qhsig(kk,1:n-1)  = qh(kk,1:n-1) + sig(kk)*(qhya(kk,1:n-1)*dy(kk,1:n-1) + qhyb(kk,1:n-1)*dy(kk,2:n) &
                     + qhTa(kk,1:n-1)*dTsoil(kk,1:n-1) + qhTb(kk,1:n-1)*dTsoil(kk,2:n))
                qhsig(kk,n)      = qh(kk,n) + sig(kk)*(qhya(kk,n)*dy(kk,n) + qhTa(kk,n)*dTsoil(kk,n))
                qvsig(kk,1:n-1)  = qv(kk,1:n-1) + sig(kk)*(qvya(kk,1:n-1)*dy(kk,1:n-1) + qvyb(kk,1:n-1)*dy(kk,2:n) &
                     + qvTa(kk,1:n-1)*dTsoil(kk,1:n-1) + qvTb(kk,1:n-1)*dTsoil(kk,2:n))
                qvsig(kk,n)      = zero

                qvTsig(kk,0)     = qvT(kk,0) + sig(kk)*qvTb(kk,0)*dTsoil(kk,1)
                qvTsig(kk,1:n-1) = qvT(kk,1:n-1) + sig(kk)*(qvTa(kk,1:n-1)*dTsoil(kk,1:n-1) + qvTb(kk,1:n-1)*dTsoil(kk,2:n))
                qvTsig(kk,n)     = zero
                qlsig(kk,1:n-1)  = qsig(kk,1:n-1) - qvsig(kk,1:n-1)
                qlsig(kk,n)      = qsig(kk,n)

                if (botbc=="aquifer") then ! update aquifer props
                   if (v_aquifer(kk)%isat==0) then
                      if (v_aquifer(kk)%WA+(dwdrainage(kk)-dwdischarge(kk))*dt(kk) > &
                           (v_aquifer(kk)%zzero-v_aquifer(kk)%zsoil)*v_aquifer(kk)%Sy) then
                         v_aquifer(kk)%isat = 1  ! aquifer saturated
                         S(kk,n) = S(kk,n) + (-(v_aquifer(kk)%WA+(dwdrainage(kk)-dwdischarge(kk))*dt(kk)) &
                              +(v_aquifer(kk)%zzero-v_aquifer(kk)%zsoil)*v_aquifer(kk)%Sy)/(dx(kk,n)*par(kk,n)%thre)
                      endif
                      v_aquifer(kk)%WA        = min(v_aquifer(kk)%WA+(dwdrainage(kk)-dwdischarge(kk))*dt(kk), &
                           (v_aquifer(kk)%zzero-v_aquifer(kk)%zsoil)*v_aquifer(kk)%Sy)
                      v_aquifer(kk)%zdelta    = v_aquifer(kk)%zzero - v_aquifer(kk)%Wa/v_aquifer(kk)%Sy
                      ! new discharge rate
                      v_aquifer(kk)%discharge = v_aquifer(kk)%Rsmax*exp(-v_aquifer(kk)%f*v_aquifer(kk)%zdelta)
                   elseif (v_aquifer(kk)%isat==1) then
                      ! check for desat of aquifer
                      if (dwdrainage(kk) < v_aquifer(kk)%Rsmax*exp(-v_aquifer(kk)%f*v_aquifer(kk)%zsoil)) then
                         v_aquifer(kk)%isat      = 0
                         v_aquifer(kk)%zdelta    = v_aquifer(kk)%zsoil
                         ! new discharge rate
                         v_aquifer(kk)%discharge = v_aquifer(kk)%Rsmax*exp(-v_aquifer(kk)%f*v_aquifer(kk)%zdelta)
                      endif
                   endif
                   zdelta(kk) = v_aquifer(kk)%zdelta
                endif       ! end update aquifer props

                if (botbc=="constant head") then
                   drn(kk) = drn(kk)+(q(kk,n)+sig(kk)*qya(kk,n)*dy(kk,n))*dt(kk)
                else
                   drn(kk) = drn(kk)+(q(kk,n)+sig(kk)*qya(kk,n)*dy(kk,n))*dt(kk)
                end if

                if (present(wex)) then
                   wex(kk,1:n) = wex(kk,1:n)+(iqex(kk,1:n)+spread(sig(kk),1,n)*qexd(kk,1:n)*dy(kk,1:n))*spread(dt(kk),1,n)
                end if

                if (litter .and. ns(kk)==1 ) then ! adjust litter moisture content if no ponding
                   SL(kk)      = SL(kk) + dy(kk,0)
                   deltaSL(kk) = dy(kk,0)
                else
                   deltaSL(kk) = zero
                endif

                TL(kk) = TL(kk) + de(kk,0)   ! Tlitter assumed the same as pond T
                if (SL(kk) >= one) vlit(kk)%isat = 1   ! check for litter saturation

                csoil(kk,:) = var(kk,1:n)%csoil
                kth(kk,:) = var(kk,1:n)%kth
                phi(kk,:) = var(kk,1:n)%phi

                ! update snow layer
                if (vsnow(kk)%nsnow>0) then
                   vsnow(kk)%Qadv_snow = vsnow(kk)%Qadv_snow + rhow*(qprec_snow(kk))* &
                        (csice*(min(vmet(kk)%Ta,zero))-lambdaf)*dt(kk)
                   vsnow(kk)%Qadv_rain = vsnow(kk)%Qadv_rain + rhow*(qprec(kk))*cswat*(vmet(kk)%Ta)*dt(kk)

                   ! update heat flux components
                   Tqw  = merge(vmet(kk)%Ta, vsnow(kk)%tsn(vsnow(kk)%nsnow), -qevap(kk)>zero)
                   !if (vsnow(kk)%hliq(1)>zero) then
                   vsnow(kk)%Qadv_vap(vsnow(kk)%nsnow) = vsnow(kk)%Qadv_vap(vsnow(kk)%nsnow) + &
                        rhow*(-qevap(kk))*cswat*Tqw*dt(kk) - &
                        qadvsig(kk,0)*dt(kk)
                   !       else
                   !          vsnow(kk)%Qadv_vap(vsnow(kk)%nsnow) = vsnow(kk)%Qadv_vap(vsnow(kk)%nsnow) + &
                   !                                                rhow*(-qevap(kk))*(csice*Tqw-lambdaf)*dt(kk)
                   !       endif
                   vsnow(kk)%Qcond_net(1) = vsnow(kk)%Qcond_net(1) + (qhsig_snow(kk,0) -  qhsig(kk,0))*dt(kk) - &
                        (qadvsig_snow(kk,0) -  qadvsig(kk,0))*dt(kk)

                   do i=1,vsnow(kk)%nsnow
                      if (vsnow(kk)%hliq(i)>zero) then

                         if ((vsnow(kk)%hliq(i) + de(kk,-i+1))>(dy(kk,-i+1) + vsnow(kk)%hsnow(i))) then
                            ! account here for new hliq exceeding new hsnow
                            tmp1d1(kk) = vsnow(kk)%hliq(i)*cswat*rhow*zero + &
                                 (vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow*(zero*csice-lambdaf)
                            tmp1d2(kk) = ((tmp1d1(kk) + LHS_h_snow(kk,i)*dt(kk)))/ &
                                 ((dy(kk,-i+1) + vsnow(kk)%hsnow(i))*cswat*rhow)
                            vsnow(kk)%tsn(i) = tmp1d2(kk) - zero

                            vsnow(kk)%deltaJlatent(i) = vsnow(kk)%deltaJlatent(i) + &
                                 lambdaf*(vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow
                            vsnow(kk)%deltaJsensible(i) = vsnow(kk)%deltaJsensible(i) + LHS_h_snow(kk,i)*dt(kk) - &
                                 lambdaf*(vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i))*rhow
                            vsnow(kk)%hliq(i) = vsnow(kk)%hsnow(i) + dy(kk,-i+1)
                         elseif ((vsnow(kk)%hliq(i)+de(kk,-i+1))<zero) then
                            ! account here for new hliq less than zero
                            tmp1d1(kk) = (csice*(vsnow(kk)%tsn(i))-lambdaf)*rhow*(vsnow(kk)%hsnow(i)-vsnow(kk)%hliq(i)) + &
                                 cswat*(vsnow(kk)%tsn(i))*rhow*vsnow(kk)%hliq(i) + LHS_h_snow(kk,i)*dt(kk)
                            delta_snowT(kk,i)  = (tmp1d1(kk)/rhow/(dy(kk,-i+1) + vsnow(kk)%hsnow(i))+lambdaf)/csice - zero
                            vsnow(kk)%tsn(i) = delta_snowT(kk,i)

                            vsnow(kk)%deltaJlatent(i) = vsnow(kk)%deltaJlatent(i) - lambdaf*dy(kk,-i+1)*rhow + &
                                 lambdaf*vsnow(kk)%hliq(i)*rhow
                            vsnow(kk)%deltaJsensible(i) = vsnow(kk)%deltaJsensible(i) + LHS_h_snow(kk,i)*dt(kk) - &
                                 (- lambdaf*dy(kk,-i+1)*rhow + lambdaf*vsnow(kk)%hliq(i)*rhow)
                            vsnow(kk)%hliq(i) = zero
                         else
                            vsnow(kk)%hliq(i) = vsnow(kk)%hliq(i) + de(kk,-i+1)
                            delta_snowliq(kk,i) = de(kk,-i+1)
                            delta_snowT(kk,i) = zero
                            vsnow(kk)%deltaJlatent(i) = vsnow(kk)%deltaJlatent(i) - &
                                 lambdaf*dy(kk,-i+1)*rhow + lambdaf*de(kk,-i+1)*rhow
                            vsnow(kk)%deltaJsensible(i) = vsnow(kk)%deltaJsensible(i) + &
                                 dy(kk,-i+1)*rhow*csice*(vsnow(kk)%tsn(i)) + &
                                 de(kk,-i+1)*rhow*(cswat-csice)*(vsnow(kk)%tsn(i))
                         endif
                      else !hliq = zero
                         if ((vsnow(kk)%tsn(i) + de(kk,-i+1))<zero) then
                            vsnow(kk)%deltaJlatent(i) = vsnow(kk)%deltaJlatent(i) - lambdaf*dy(kk,-i+1)*rhow
                            vsnow(kk)%deltaJsensible(i) = vsnow(kk)%deltaJsensible(i) + &
                                 dy(kk,-i+1)*rhow*csice*(vsnow(kk)%tsn(i)) + &
                                 de(kk,-i+1)*csice*rhow*hsnow(kk,1)
                            vsnow(kk)%tsn(i) = vsnow(kk)%tsn(i) + de(kk,-i+1)
                            delta_snowT(kk,i) = de(kk,-i+1)
                            delta_snowliq(kk,i) = zero
                            vsnow(kk)%hliq(i) = zero
                         else ! vsnow(kk)%tsn + de(kk,0))>zero
                            ! LHS_h_snow(kk,1) = (dy(kk,0)*rhow*(csice*(vsnow(kk)%tsn(1))-lambdaf) + &
                            !                         de(kk,0)*rhow*((cswat-csice)*(vsnow(kk)%tsn(1))+lambdaf))/dt(kk)

                            ! else
                            !  LHS_h_snow(kk,1) = (dy(kk,0)*rhow*(csice*(vsnow(kk)%tsn(1))-lambdaf) + &
                            !                          de(kk,0)*csice*rhow*hsnow(kk,1))/dt(kk)
                            delta_snowT(kk,i) = zero - vsnow(kk)%tsn(i)
                            delta_snowliq(kk,i) = (LHS_h_snow(kk,1)*dt(kk) - &
                                 delta_snowT(kk,vsnow(kk)%nsnow)*csice*hsnow(kk,i) *rhow - &
                                 dy(kk,-i+1)*rhow*(csice*(vsnow(kk)%tsn(i))-lambdaf))/ &
                                 (rhow*(zero*(cswat-csice)+lambdaf))

                            vsnow(kk)%hliq(i) = vsnow(kk)%hliq(i) + delta_snowliq(kk,vsnow(kk)%nsnow)
                            vsnow(kk)%deltaJlatent(i) = vsnow(kk)%deltaJlatent(i) &
                                 - lambdaf*dy(kk,-i+1)*rhow + &
                                 lambdaf*rhow*delta_snowliq(kk,vsnow(kk)%nsnow)
                            vsnow(kk)%deltaJsensible(i) = vsnow(kk)%deltaJsensible(i)  &
                                 + dy(kk,-i+1)*rhow*csice*(vsnow(kk)%tsn(i)) + &
                                 delta_snowT(kk,i)*csice*rhow*hsnow(kk,i)
                            vsnow(kk)%tsn(i) = vsnow(kk)%tsn(i) + delta_snowT(kk,i)

                         endif ! vsnow(kk)%tsn + de(kk,0))>zero
                      endif ! hliq = zero
                      vsnow(kk)%hsnow(i) = dy(kk,-i+1) + vsnow(kk)%hsnow(i)
                      delta_snowcol(kk,i) = dy(kk,-i+1)
                      vsnow(kk)%depth(i) = vsnow(kk)%hsnow(i)/(vsnow(kk)%dens(i)/rhow)
                      ! if (i==vsnow(kk)%nsnow) then
                      !    vsnow(kk)%FluxDivergence(i) = vsnow(kk)%Qadv_rain + vsnow(kk)%Qadv_snow + &
                      !                                  vsnow(kk)%Qadv_vap(i) + &
                      !                                + vsnow(kk)%Qadv_melt(i) + vsnow(kk)%Qcond_net(i)  +  &
                      !                                  vsnow(kk)%Qadv_transfer(i)
                      ! else
                      !    vsnow(kk)%FluxDivergence(i) = vsnow(kk)%Qadv_vap(i) + &
                      !                                + vsnow(kk)%Qadv_melt(i) + vsnow(kk)%Qcond_net(i)  + &
                      !                                  vsnow(kk)%Qadv_transfer(i)
                      ! endif

                      ! if (abs(vsnow(kk)%FluxDivergence(i)- &
                      !     (vsnow(kk)%deltaJsensible(i)+vsnow(kk)%deltaJlatent(i)))>1e3) then
                      !    write(*,*)
                      ! endif
                   enddo ! loop over snow layers
                   vsnow(kk)%wcol = sum(vsnow(kk)%hsnow(:))
                endif ! end update snowlayers
                vsnow(kk)%FluxDivergence(1) = vsnow(kk)%Qadv_rain + vsnow(kk)%Qadv_snow + vsnow(kk)%Qadv_vap(1) + &
                     vsnow(kk)%Qadv_melt(1) + vsnow(kk)%Qcond_net(1)  +  vsnow(kk)%Qadv_transfer(1)
                ! heat stored in snowpack
                where (vsnow(kk)%hsnow(:).gt.zero)
                   vsnow(kk)%Jsensible(:) = (vsnow(kk)%hsnow(:)-vsnow(kk)%hliq(:))*rhow*csice*(vsnow(kk)%Tsn(:)) + &
                        vsnow(kk)%hliq(:)*rhow*cswat*(vsnow(kk)%Tsn(:))

                   vsnow(kk)%Jlatent(:) = (vsnow(kk)%hsnow(:)-vsnow(kk)%hliq(:))*rhow*(-lambdaf)
                elsewhere
                   vsnow(kk)%Jsensible(:) = zero
                   vsnow(kk)%Jlatent(:) = zero
                end where
             endif ! .not. again

             ! update variables (S,T) to end of time step
             ! update variables to sig for use in isotope routine
             do i=1, n
                if (.not.again(kk)) Tsoil(kk,i)  = Tsoil(kk,i) + dTsoil(kk,i)
                if (var(kk,i)%isat==0) then
                   if (.not.again(kk)) then
                      deltaS(kk,i) = dy(kk,i)  !required for isotope subroutine
                      S(kk,i)      = S(kk,i)+dy(kk,i)
                      if (S(kk,i)>one .and. dy(kk,i)>zero) then ! onset of saturation of layer
                         var(kk,i)%isat = 1
                         var(kk,i)%K    = var(kk,i)%Ksat
                         var(kk,i)%phi  = var(kk,i)%phie
                      endif
                   endif
                else
                   deltaS(kk,i)  = zero    !required for isotope subroutine
                   if (i==1) then
                      var(kk,i)%phi = var(kk,i)%phi + dy(kk,i)*real(ns(kk),r_2) ! pond included in top soil layer
                   else
                      var(kk,i)%phi = var(kk,i)%phi + dy(kk,i)
                   endif

                   ! var(kk,i)%phi = zero from Ross
                   ! VH thinks it is o.k. because hyofS is called later again
                   ! MC think that we should probably get rid of -he*Ksat in phip etc.
                   !   because we merged the pond with the first layer.
                   if (i==1 .and. ih0(kk)/=0 .and. var(kk,i)%phi>=var(kk,i)%phie) var(kk,i)%phi = zero ! pond gone
                   if (var(kk,i)%phi<var(kk,i)%phie .and. var(kk,i)%iice==0) then ! desaturation of layer
                      var(kk,i)%isat = 0
                      var(kk,i)%K    = var(kk,i)%Ksat
                      var(kk,i)%phi  = var(kk,i)%phie
                      var(kk,i)%KS   = par(kk,i)%KSe
                      var(kk,i)%phiS = par(kk,i)%phiSe
                   elseif (var(kk,i)%phi<var(kk,i)%phie .and. var(kk,i)%iice==1) then
                      var(kk,i)%isat = 0
                      var(kk,i)%K    = var(kk,i)%Ksat
                      var(kk,i)%phi  = var(kk,i)%phie
                      var(kk,i)%KS   = var(kk,i)%KS
                      var(kk,i)%phiS = var(kk,i)%phiS
                   endif
                end if

                if (.not. again(kk)) then
                   if (h0(kk)<zero .and. var(kk,1)%isat==0 .and. (i==1)) then ! start negative pond correction
                      infil(kk)     = infil(kk)+h0(kk)
                      ! negative pond converted to loss of soil moisture in top layer
                      S(kk,1) = S(kk,1) + h0(kk)/(dx(kk,1)*par(kk,1)%thre)
                      deltaS(kk,1) = h0(kk)/(dx(kk,1)*par(kk,1)%thre)
                      deltah0(kk) = -(h0(kk)-deltah0(kk)) ! whole pond lost (new change = - original pond height)
                      hice(kk) = zero
                      h0(kk) = zero  ! zero pond remaining
                   endif


                   ! corrections for onset of freezing and melting
                   ! calculate freezing point temperature at new S
                   Tfreezing(kk) = Tfrz(S(kk,i),par(kk,i)%he,one/par(kk,i)%lam)
                   ! check for onset of freezing and adjust (increase) temperature to account for latent heat
                   ! release by ice formation

                   if (Tsoil(kk,i) < Tfreezing(kk) .and. var(kk,i)%iice==0) then  ! start correction for onset of freezing
                      dtdT(kk)      = dthetalmaxdTh(Tfreezing(kk), S(kk,i), par(kk,i)%he, one/par(kk,i)%lam, &
                           par(kk,i)%thre,par(kk,i)%the)
                      tmp1d3(kk) = 1000.
                      tmp1d2(kk) = Tsoil(kk,i) - dTsoil(kk,i)
                      k = 1
                      if ((i==1).and.(h0(kk)- deltah0(kk))>zero) then
                         h0(kk) = h0(kk) - deltah0(kk)  ! calculate stored heat using h0 at t-dt

                         do while ((k<nsteps_ice_max).and.(tmp1d3(kk)>tol_dthetaldT))
                            tmp1d1(kk)  = LHS_h(kk,i)*dt(kk)-deltaJ_sensible_S(kk,i) + tmp1d2(kk)* &
                                 (var(kk,i)%csoil*dx(kk,i)+cp(kk)*h0(kk)) + &
                                 Tfreezing(kk)*rhow*(lambdaf+(tmp1d2(kk))*(cswat-csice)) * &
                                 dtdT(kk)*(dx(kk,i)+h0(kk)/par(kk,i)%thre)
                            tmp1d1(kk)  = tmp1d1(kk) /( (var(kk,i)%csoil*dx(kk,i)+cp(kk)*h0(kk)) + &
                                 rhow*(lambdaf+(tmp1d2(kk))*(cswat-csice)) * &
                                 dtdT(kk)*(dx(kk,i)+h0(kk)/par(kk,i)%thre))
                            dtdT(kk) = (thetalmax(tmp1d1(kk), S(kk,i), par(kk,i)%he, one/par(kk,i)%lam, &
                                 par(kk,i)%thre, par(kk,i)%the) &
                                 - (S(kk,i)*par(kk,i)%thre + (par(kk,i)%the-par(kk,i)%thre)))/ &
                                 (tmp1d1(kk) - Tfreezing(kk))
                            tmp1d3(kk) = abs((tmp1d1(kk)- Tfreezing(kk))*dtdT(kk))
                            k = k + 1
                         enddo  ! end of do while
                         dTsoil(kk,i) = tmp1d1(kk) - tmp1d2(kk)
                         Tsoil(kk,i) = tmp1d1(kk)
                         deltaJ_sensible_T(kk,i) = (var(kk,i)%csoil*dx(kk,i)+ cp(kk)*h0(kk))*dTsoil(kk,i) + &
                              (Tsoil(kk,i)-Tfreezing(kk))*rhow*dtdT(kk)*(dx(kk,1)+h0(kk)/par(kk,1)%thre)* &
                              ((cswat-csice)*(tmp1d2(kk)))
                         deltaJ_latent_T(kk,i) = (Tsoil(kk,i)-Tfreezing(kk))*rhow*dtdT(kk)*lambdaf* &
                              (dx(kk,1)+h0(kk)/par(kk,1)%thre)
                         h0(kk) = h0(kk) + deltah0(kk)
                      else
                         do while ((k<nsteps_ice_max).and.(tmp1d3(kk)>tol_dthetaldT))
                            ! zeroth estimate of corrected temperature
                            tmp1d1(kk)  = LHS_h(kk,i)*dt(kk)-deltaJ_sensible_S(kk,i) + tmp1d2(kk)* &
                                 (var(kk,i)%csoil*dx(kk,i)) + Tfreezing(kk)*rhow*(lambdaf+(tmp1d2(kk))*(cswat-csice)) * &
                                 dtdT(kk)*dx(kk,i)
                            tmp1d1(kk)  = tmp1d1(kk) /( (var(kk,i)%csoil*dx(kk,i)) + &
                                 rhow*(lambdaf+(tmp1d2(kk))*(cswat-csice)) * dtdT(kk)*dx(kk,i))
                            dtdT(kk) =  (thetalmax(tmp1d1(kk), S(kk,i), par(kk,i)%he, one/par(kk,i)%lam, par(kk,i)%thre, &
                                 par(kk,i)%the) &
                                 - (S(kk,i)*par(kk,i)%thre + (par(kk,i)%the-par(kk,i)%thre)))/ &
                                 (tmp1d1(kk) - Tfreezing(kk))
                            tmp1d3(kk) = abs((tmp1d1(kk)- Tfreezing(kk))*dtdT(kk))
                            k = k + 1
                         enddo  ! end of do while
                         dTsoil(kk,i) = tmp1d1(kk)- (Tsoil(kk,i)-dTsoil(kk,i))
                         Tsoil(kk,i) = tmp1d1(kk)
                         deltaJ_sensible_T(kk,i) = (var(kk,i)%csoil*dx(kk,i))*dTsoil(kk,i) + &
                              (Tsoil(kk,i)-Tfreezing(kk))*rhow*dtdT(kk)*dx(kk,i)* &
                              ((cswat-csice)*(tmp1d2(kk)))
                         deltaJ_latent_T(kk,i) = (Tsoil(kk,i)-Tfreezing(kk))*rhow*dtdT(kk)*lambdaf*dx(kk,i)
                      endif
                   endif  ! end correction for onset of freezing

                   if (var(kk,i)%iice==1) then
                      theta         = S(kk,i)*(par(kk,i)%thre) + (par(kk,i)%the - par(kk,i)%thre)
                      tmp1d1(kk) = Tsoil(kk,i) - dTsoil(kk,i)
                      ! energy for complete melting
                      if ((i==1).and.(h0(kk)>zero)) then


                         tmp1d3(kk) = thetai(kk,i)*dx(kk,i)*rhow*lambdaf + &
                              thetai(kk,i)*(h0(kk)- deltah0(kk))/par(kk,i)%thre*rhow*lambdaf

                         tmp1d2(kk) = (Tfreezing(kk))*(rhow*cswat*(theta*dx(kk,1)+h0(kk))+ &
                              par(kk,1)%rho*par(kk,1)%css*dx(kk,1)) - &
                              (tmp1d1(kk))*(rhow*cswat*(var(kk,i)%thetal*dx(kk,1)+(h0(kk)- deltah0(kk))* &
                              (one-thetai(kk,1)/par(kk,1)%thre))+par(kk,1)%rho*par(kk,1)%css*dx(kk,1)) - &
                              (tmp1d1(kk))*(rhow*csice*(thetai(kk,i)*dx(kk,1)+(h0(kk)- deltah0(kk))* &
                              thetai(kk,1)/par(kk,1)%thre))
                      else
                         tmp1d3(kk) = thetai(kk,i)*dx(kk,i)*rhow*lambdaf

                         tmp1d2(kk) = (Tfreezing(kk))*(rhow*cswat*(theta*dx(kk,i))+par(kk,i)%rho*par(kk,i)%css*dx(kk,i)) - &
                              (tmp1d1(kk))*(rhow*cswat*(var(kk,i)%thetal*dx(kk,i))+par(kk,i)%rho*par(kk,i)%css*dx(kk,i)) - &
                              (tmp1d1(kk))*(rhow*csice*(thetai(kk,i)*dx(kk,i)))
                      endif
                   endif
                   ! check for onset of thawing and decrease temperature to account for latent heat required to melt ice
                   ! if ((var(kk,i)%iice==1).and.((LHS_h(kk,i)*dt(kk)-tmp1d2(kk)-tmp1d3(kk))>=zero)) then
                   if ((var(kk,i)%iice==1).and.((J0(kk,i) + LHS_h(kk,i)*dt(kk)).ge.JSoilLayer(Tfreezing(kk), &
                        dx(kk,i), theta,par(kk,i)%css, par(kk,i)%rho, &
                        merge(h0(kk),zero,i==1), par(kk,i)%thre, par(kk,i)%the, &
                        par(kk,i)%he, one/par(kk,i)%lam))) then
                      ! start correction for onset of thawing
                      ! Correct for "overthawing". This extra energy is used to heat the soil even more
                      ! The correction comes from equating the energy balances before and after correction.
                      tmp1d1(kk) = Tsoil(kk,i) - dTsoil(kk,i)

                      if ((i==1) .and. (h0(kk)>zero)) then
                         deltaJ_latent_T(kk,i)=      tmp1d3(kk)
                         deltaJ_latent_S(kk,i) = zero
                         h0(kk) = h0(kk) - deltah0(kk)

                         deltaJ_sensible_S(kk,i) = (tmp1d1(kk))*(rhow*cswat*(theta*dx(kk,1)+(h0(kk)+deltah0(kk)))+ &
                              par(kk,1)%rho*par(kk,1)%css*dx(kk,1)) - &
                              (tmp1d1(kk))*(rhow*cswat*(var(kk,i)%thetal*dx(kk,1)+(h0(kk))* &
                              (one-thetai(kk,1)/par(kk,1)%thre))+par(kk,1)%rho*par(kk,1)%css*dx(kk,1)) - &
                              (tmp1d1(kk))*(rhow*csice*(thetai(kk,i)*dx(kk,1)+(h0(kk))* &
                              thetai(kk,1)/par(kk,1)%thre))

                         dTsoil(kk,i) = (LHS_h(kk,i)*dt(kk) - (deltaJ_latent_T(kk,i) + deltaJ_sensible_S(kk,i)))/ &
                              (rhow*cswat*(theta*dx(kk,1)+(h0(kk)+deltah0(kk)))+par(kk,1)%rho*par(kk,1)%css*dx(kk,1))

                         deltaJ_sensible_T(kk,i) = dTsoil(kk,i)* &
                              (rhow*cswat*(theta*dx(kk,1)+(h0(kk)+deltah0(kk)))+par(kk,1)%rho*par(kk,1)%css*dx(kk,1))

                         h0(kk) = h0(kk) + deltah0(kk)
                      else
                         deltaJ_latent_T(kk,i) = thetai(kk,i)*dx(kk,i)*rhow*lambdaf
                         deltaJ_latent_S(kk,i) = zero

                         theta         = S(kk,i)*(par(kk,i)%thre) + (par(kk,i)%the - par(kk,i)%thre)
                         var(kk,i)%csoil = cswat*theta*rhow + par(kk,i)%css*par(kk,i)%rho

                         deltaJ_sensible_S(kk,i) = (tmp1d1(kk))*(rhow*cswat*(theta*dx(kk,i))+ &
                              par(kk,i)%rho*par(kk,i)%css*dx(kk,i)) - &
                              (tmp1d1(kk))*(rhow*cswat*(var(kk,i)%thetal*dx(kk,i))+par(kk,i)%rho*par(kk,i)%css*dx(kk,i)) - &
                              (tmp1d1(kk))*(rhow*csice*(thetai(kk,i)*dx(kk,i)))
                         dTsoil(kk,i) = (LHS_h(kk,i)*dt(kk) - (deltaJ_latent_T(kk,i) + deltaJ_sensible_S(kk,i)))/ &
                              (dx(kk,i)*var(kk,i)%csoil)
                         deltaJ_sensible_T(kk,i) = var(kk,i)%csoil*dTsoil(kk,i)*dx(kk,i)
                      endif
                      Tsoil(kk,i) = tmp1d1(kk) + dTsoil(kk,i)
                      ! correct thetai and T in frozen soil
                      !elseif ((var(kk,i)%iice==1).and.((LHS_h(kk,i)*dt(kk)-tmp1d2(kk)-tmp1d3(kk))<zero)) then
                   elseif ((var(kk,i)%iice==1).and.((J0(kk,i) + LHS_h(kk,i)*dt(kk))<JSoilLayer(Tfreezing(kk), &
                        dx(kk,i), theta,par(kk,i)%css, par(kk,i)%rho, &
                        merge(h0(kk),zero,i==1), par(kk,i)%thre, par(kk,i)%the, &
                        par(kk,i)%he, one/par(kk,i)%lam))) then


                      tmp1d2(kk) = J0(kk,i) + LHS_h(kk,i)*dt(kk) ! total energy in  soil layer
                      !check there is a zero
                      tmp1 = GTfrozen(real(Tsoil(kk,i)-50., r_2), tmp1d2(kk), dx(kk,i), theta,par(kk,i)%css, par(kk,i)%rho, &
                           merge(h0(kk),zero,i==1), par(kk,i)%thre, par(kk,i)%the, par(kk,i)%he, one/par(kk,i)%lam)
                      tmp2 = GTFrozen(Tfreezing(kk), tmp1d2(kk), dx(kk,i), theta,par(kk,i)%css, par(kk,i)%rho, &
                           merge(h0(kk),zero,i==1), par(kk,i)%thre, par(kk,i)%the, par(kk,i)%he, one/par(kk,i)%lam)
                      ! there is a zero in between
                      if ((tmp1*tmp2) < zero) then
                         tmp1d3(kk) = rtbis_Tfrozen(tmp1d2(kk), dx(kk,i), theta,par(kk,i)%css, par(kk,i)%rho, &
                              merge(h0(kk),zero,i==1), par(kk,i)%thre, par(kk,i)%the, par(kk,i)%he, one/par(kk,i)%lam, &
                              real(Tsoil(kk,i)-50., r_2), Tfreezing(kk), real(0.0001,r_2))
                         tmp1d4(kk) = thetalmax(tmp1d3(kk), min(S(kk,i),one), par(kk,i)%he, one/par(kk,i)%lam, &
                              par(kk,i)%thre, par(kk,i)%the) ! liquid content at inital Tsoil
                      else
                         write(*,*) "Found no solution for Tfrozen 1. Stop. ", kk, i
                         write(*,*) nsteps(kk), S(kk,i), Tsoil(kk,i), h0(kk), tmp1, tmp2, tmp1d2(kk), theta
                         stop
                      endif
                      var(kk,i)%thetal = tmp1d4(kk)
                      var(kk,i)%thetai = theta - tmp1d4(kk)
                      if (i==1) then
                         hice_tmp(kk) = hice(kk)
                         hice(kk) = h0(kk)*var(kk,1)%thetai/par(kk,1)%thre

                         ! correct total energy stored in pond + soil
                         deltaJ_latent_S(kk,i) = - rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
                              dx(kk,1)*(var(kk,1)%thetai-thetai(kk,1)))

                         deltaJ_latent_T(kk,i) = zero

                         Jcol_sensible(kk) = Jcol_sensible(kk) - deltaJ_sensible_T(kk,i)

                         deltaJ_sensible_T(kk,i) = &
                              dx(kk,1)*tmp1d3(kk)*par(kk,1)%rho*par(kk,1)%css + &
                              dx(kk,1)*var(kk,1)%thetal*cswat*rhow*tmp1d3(kk) + &
                              (h0(kk)-hice(kk))*cswat*rhow*tmp1d3(kk) + &
                              var(kk,1)%thetai*dx(kk,1)*rhow*(csice*(tmp1d3(kk))-lambdaf) + &
                              hice(kk)*rhow*(csice*(tmp1d3(kk))-lambdaf)- &
                              J0(kk,1) - &
                              deltaJ_latent_S(kk,1)

                         Jcol_sensible(kk) = Jcol_sensible(kk) + deltaJ_sensible_T(kk,i)

                         deltaJ_sensible_S(kk,i) = 0.

                      else
                         ! correct total energy stored in  soil
                         deltaJ_latent_S(kk,i) = - rhow*lambdaf*( + &
                              dx(kk,i)*(var(kk,i)%thetai-thetai(kk,i)))

                         deltaJ_latent_T(kk,i) = zero

                         Jcol_sensible(kk) = Jcol_sensible(kk) - deltaJ_sensible_T(kk,i)

                         deltaJ_sensible_T(kk,i) = &
                              dx(kk,i)*tmp1d3(kk)*par(kk,i)%rho*par(kk,i)%css + &
                              dx(kk,i)*var(kk,i)%thetal*cswat*rhow*tmp1d3(kk) + &
                              var(kk,i)%thetai*dx(kk,i)*rhow*(csice*(tmp1d3(kk))-lambdaf) + &
                              J0(kk,i) - &
                              deltaJ_latent_S(kk,i)

                         Jcol_sensible(kk) = Jcol_sensible(kk) + deltaJ_sensible_T(kk,i)

                         deltaJ_sensible_S(kk,i) = 0.
                      endif
                      thetai(kk,i) = var(kk,i)%thetai
                      ! ! debug
                      ! if ((irec.eq.10).and.(nsteps(kk).eq.32).and.(i.eq.1)) then
                      !    write(*,*) 'writing diags'


                      !    write(*,*) Tsoil(kk,1) , tmp1d3(kk) , dTsoil(kk,1), (tmp1d3(kk)-(Tsoil(kk,1)-dTsoil(kk,1)))

                      !    !if (nsteps(kk).gt.100) STOP

                      ! endif
                      Tsoil(kk,i) = tmp1d3(kk)
                   endif ! soil remains frozen

                endif ! if .not.again

             end do ! i=1, n => update variables S,T

             if (.not. again(kk)) then
                cv0(kk,1:n)   = var(kk,1:n)%cv  ! save cv before updating
                isave(kk,1:n) = var(kk,1:n)%iice ! save isat before updating
                ! update variables (particularly thetaice)
                ! Debug for mp=1: remove elemental from hyofS and do loop instead of next line
                call hyofS(S(kk,1:n), Tsoil(kk,1:n), par(kk,1:n), var(kk,1:n))
                ! do i=1, n
                !    call hyofS(S(kk,i), Tsoil(kk,i), par(kk,i), var(kk,i))
                ! end do
                ! End debug hyofS
                var(kk,1:n)%iice = isave(kk,1:n)
             endif ! if .not.again

             iflux(kk) = iflux(kk) + 1

          end do ! do while (again(kk)) ! iflux loop

          !nsteps(kk)        = nsteps(kk) + 1
          runoff(kk) = runoff(kk) + qrunoff(kk)*dt(kk)

          S0(kk,1:n)    = S(kk,1:n) - deltaS(kk,1:n)
          Sliqice0(kk,1:n) = (S0(kk,1:n) - cv0(kk,1:n))/(one-cv0(kk,1:n))
          deltacv(kk,1:n)   = var(kk,1:n)%cv - cv0(kk,1:n)
          Sliqice(kk,1:n)      = (S(kk,1:n) - var(kk,1:n)%cv)/(one-var(kk,1:n)%cv)
          deltaSliqice(kk,1:n) = Sliqice(kk,1:n) - Sliqice0(kk,1:n)

          delthetai(kk,1:n) = (var(kk,1:n)%thetai-thetai_0(kk,1:n))

          hice(kk) = h0(kk)*var(kk,1)%thetai/par(kk,1)%thre
          deltahice(kk) =  hice(kk) - hice_0(kk)
          deltahice(kk) =  hice(kk) - h0_0(kk)*thetai_0(kk,1)/par(kk,1)%thre

          ! change in snow pack (lumped with top layer)
          if (vsnow(kk)%nsnow==0.) then
             if (var(kk,1)%iice==1) then
                vsnow(kk)%wcol = vsnow(kk)%wcol + &
                     min(qprec_snow(kk)*dt(kk), max(zero,(deltahice(kk)+dx(kk,1)*delthetai(kk,1)))) + & ! accumulation
                     max(-vsnow(kk)%wcol,min(zero,(deltahice(kk)+dx(kk,1)*delthetai(kk,1)))) ! melting

                if (h0(kk)>zero) then
                   vsnow(kk)%wcol = min(vsnow(kk)%wcol, hice(kk))
                else
                   vsnow(kk)%wcol = min(vsnow(kk)%wcol, dx(kk,1)*var(kk,1)%thetai)
                endif
             else
                vsnow(kk)%wcol = zero
             endif

             vsnow(kk)%tsn = merge(Tsoil(kk,1),zero,vsnow(kk)%wcol>zero)
             if (vsnow(kk)%wcol>zero) then
                vsnow(kk)%depth(1) = vsnow(kk)%wcol/(vsnow(kk)%dens(1)/rhow)
             else
                vsnow(kk)%depth(1) = zero
             endif

          endif

          Sice(kk,1:n) = var(kk,1:n)%thetai/par(kk,1:n)%thre
          deltaSice(kk,1:n) = delthetai(kk,1:n)/par(kk,1:n)%thre

          Sliq(kk,1:n) = Sliqice(kk,1:n) - Sice(kk,1:n)
          deltaSliq(kk,1:n) = deltaSliqice(kk,1:n) - deltaSice(kk,1:n)

          Sice(kk,1) = Sice(kk,1)  + hice(kk)/(dx(kk,1)*par(kk,1)%thre)
          deltaSice(kk,1) = deltaSice(kk,1) + deltahice(kk)/(dx(kk,1)*par(kk,1)%thre)

          Sliq(kk,1) = Sliq(kk,1) + (h0(kk)-hice(kk))/(dx(kk,1)*par(kk,1)%thre) ! add pond component to Sliq(kk,1)
          deltaSliq(kk,1) = deltaSliq(kk,1) + (deltah0(kk)-deltahice(kk))/(dx(kk,1)*par(kk,1)%thre)

          thetai_0(kk,1:n)    = var(kk,1:n)%thetai
          thetai(kk,1:n) = var(kk,1:n)%thetai
          dthetaldt(kk,1:n) = var(kk,1:n)%dthetaldT
          h0_0(kk) = h0(kk)
          hice_0(kk) = hice(kk)
          Sliq0(kk,1:n) = Sliq(kk,1:n)
          Sice0(kk,1:n)  = Sice(kk,1:n)
          Sliqice0(kk,1:n) = Sliqice(kk,1:n)
          if (littercase==1) then  !!!! vh needs attention !!!!
             SL0(kk)    = SL(kk) - deltaSL(kk)
             cvL0(kk)   = vlit(kk)%cv
             SLliq0(kk) = (SL0(kk) - cvL0(kk))/(one-cvL0(kk))
             call litter_props(Sl(kk), Tl(kk), vlit(kk), plit(kk), h0(kk))
             deltacvL(kk)   = vlit(kk)%cv - cvL0(kk)
             SLliq(kk)      = (SL(kk) - vlit(kk)%cv)/(one-vlit(kk)%cv)
             deltaSLliq(kk) = SLliq(kk) - SLliq0(kk)
          endif

          Jsensible(kk,1) = (var(kk,1)%csoil* dx(kk,1)+h0(kk)*cswat)*(Tsoil(kk,1))
          Jsensible(kk,2:n) = var(kk,2:n)%csoil*(Tsoil(kk,2:n))* dx(kk,2:n)

          ! change in heat stored in soil column
          dJcol_latent_S(kk) = sum(deltaJ_latent_S(kk,1:n))
          Jcol_latent_S(kk) = Jcol_latent_S(kk) + dJcol_latent_S(kk)
          dJcol_latent_T(kk) = sum(deltaJ_latent_T(kk,1:n))
          Jcol_latent_T(kk) = Jcol_latent_T(kk) + dJcol_latent_T(kk)
          dJcol_sensible(kk) = sum(deltaJ_sensible_T(kk,1:n)) + sum(deltaJ_sensible_S(kk,1:n))
          Jcol_sensible(kk) = Jcol_sensible(kk) + dJcol_sensible(kk)

          deltaice_cum_S(kk) = -Jcol_latent_S(kk)/rhow/lambdaf
          deltaice_cum_T(kk) = -Jcol_latent_T(kk)/rhow/lambdaf

          deltaTa(kk)       = zero
          init(kk)          = .false.

          if (isotopologue /= 0) then

             tmp_thetasat(kk,1:n)   = par(kk,1:n)%thre
             tmp_tortuosity(kk,1:n) = par(kk,1:n)%tortuosity
             tmp_thetar(kk,1:n)     = par(kk,1:n)%the - par(kk,1:n)%thre
             ql0(kk) = qlsig(kk,0)
             qv0(kk) = qvsig(kk,0)

             call isotope_vap(isotopologue, litter, n, &
                  1, dx(kk,1:n), dz(kk,1:n-1), sig(kk), dt(kk), dxL(kk), &
                  Tsoil(kk,1:n), dTsoil(kk,1:n), dT0(kk), Sliqice(kk,1:n), deltaSliqice(kk,1:n), &
                  Sliq(kk,1:n), deltaSliq(kk,1:n),Sice(kk,1:n), deltaSice(kk,1:n), &
                  Tsurface(kk), TL(kk), T0(kk), h0(kk), deltah0(kk), SLliq(kk), deltaSLliq(kk), &
                  qsig(kk,0:n), qlsig(kk,0:n), qvsig(kk,0:n), &
                  qprec(kk), qevapsig(kk), qrunoff(kk), iqex(kk,1:n), qd(kk), &
                  var(kk,1:n), tmp_thetasat(kk,1:n), tmp_thetar(kk,1:n), tmp_tortuosity(kk,1:n), &
                  deltacv(kk,1:n), vmet(kk)%ra, vmet(kk)%rs, vlit(kk), vmet(kk)%cva, vmet(kk)%civa, &
                  plit(kk)%the, deltacvL(kk), &
                  cprec(kk), icali(kk), &
                  ql0(kk), qv0(kk), &
                  ciso(kk,0:n),cisoice(kk,0:n), ciso0(kk), cisoL(kk), cisos(kk), &
                  qiso_in(kk), qiso_out(kk), qiso_evap(kk), qiso_trans(kk), &
                  qiso_liq_adv(kk,1:n), qiso_vap_adv(kk,1:n), qiso_liq_diff(kk,1:n-1), qiso_vap_diff(kk,1:n-1))

             qiso_evap_cum(kk)  = qiso_evap_cum(kk)  + qiso_evap(kk)*dt(kk)
             qiso_trans_cum(kk) = qiso_trans_cum(kk) + qiso_trans(kk)*dt(kk)

          endif ! isotopologue/=0

       end do ! while (t<tfin)

    end do ! kk=1, mp

    ! get heads if required
    if (present(heads)) then
       isave    = var%isat
       var%isat = isave
       heads    = var%h
       where (S(:,:) >= one) heads(:,:) = par(:,:)%he + (var(:,:)%phi-par(:,:)%phie)/par(:,:)%Ke
    end if

  END SUBROUTINE solve

  !**********************************************************************************************************************

  ! SUBROUTINE solute(ti,tf,thi,thf,win,cin,n,ns,dx,jt,dsmmax,sm,sdrn,nssteps,c, &
  !      isosub)

  !   USE sli_utils, ONLY: dis, isotype, bd, isopar

  !   IMPLICIT NONE

  !   INTEGER(i_d),INTENT(IN)::n,ns,jt(n)
  !   REAL(r_2),INTENT(IN)::ti,tf,thi(n),thf(n),win,cin(ns),dx(n),dsmmax
  !   INTEGER(i_d),INTENT(INOUT)::nssteps(ns)
  !   REAL(r_2),INTENT(INOUT)::sm(n,ns),sdrn(ns),c(n,ns)
  !   OPTIONAL::isosub
  !   INTERFACE
  !      SUBROUTINE isosub(iso,c,p,f,fc)
  !        USE sli_numbers, ONLY: r_2
  !        CHARACTER(LEN=2),INTENT(IN)::iso
  !        REAL(r_2),INTENT(IN)::c
  !        REAL(r_2),DIMENSION(:),INTENT(INOUT)::p
  !        REAL(r_2),INTENT(OUT)::f, fc
  !      END SUBROUTINE isosub
  !   END INTERFACE
  !   ! Solves the ADE from time ti to tf. Diffusion of solute ignored - dispersion
  !   ! coeff = dispersivity * abs(pore water velocity).
  !   ! Definitions of arguments:
  !   ! Required args:
  !   ! ti   - start time (h).
  !   ! tf   - finish time.
  !   ! thi(1:n)  - initial layer water contents.
  !   ! thf(1:n)  - final layer water contents.
  !   ! win   - water in at top of profile.
  !   ! cin(1:ns)  - solute concn in win.
  !   ! n    - no. of soil layers.
  !   ! ns   - no. of solutes.
  !   ! dx(1:n)  - layer thicknesses.
  !   ! jt(1:n)  - layer soil type nos.
  !   ! dsmmax(1:ns) - max change in sm of any layer to aim for each time step;
  !   !      controls time step size.
  !   ! sm(1:n,1:ns) - layer masses of solute per cc.
  !   ! sdrn(1:ns) - cumulative solute drainage.
  !   ! nssteps(1:ns) - cumulative no. of time steps for ADE soln.
  !   ! Optional args:
  !   ! isosub  - subroutine to get adsorbed solute (units/g soil) from concn
  !   !      in soil water according to chosen isotherm code.
  !   !      Arguments: iso - 2 character code; c - concn in soil water;
  !   !      p(:) - isotherm parameters; f - adsorbed mass/g soil;
  !   !      fc - deriv of f wrt c (slope of isotherm curve).
  !   INTEGER(i_d),PARAMETER::itmax=20 ! max iterations for finding c from sm
  !   REAL(r_2),PARAMETER::eps=0.00001 ! for stopping
  !   INTEGER(i_d)::i,it,j,k
  !   REAL(r_2)::dc,dm,dmax,dt,dz(n-1),f,fc,r,rsig,rsigdt,sig,sigdt,t,tfin,th,v1,v2
  !   REAL(r_2),DIMENSION(n-1)::coef1,coef2
  !   REAL(r_2),DIMENSION(n)::csm,tht
  !   REAL(r_2),DIMENSION(0:n)::aa,bb,cc,dd,dy,q,qw,qya,qyb
  !   INTEGER(i_d) :: info

  !   sig=half
  !   rsig=one/sig
  !   tfin=tf
  !   dz=half*(dx(1:n-1)+dx(2:n))
  !   !get average water fluxes
  !   r=one/(tf-ti)
  !   qw(0)=r*win
  !   tht=r*(thf-thi)
  !   do i=1,n
  !      qw(i)=qw(i-1)-dx(i)*tht(i)
  !   end do
  !   !get constant coefficients
  !   do i=1,n-1
  !      v1=half*qw(i)
  !      v2=half*(dis(jt(i))+dis(jt(i+1)))*abs(qw(i))/dz(i)
  !      coef1(i)=v1+v2
  !      coef2(i)=v1-v2
  !   end do
  !   do j=1,ns
  !      t=ti
  !      if (qw(0)>zero) then
  !         q(0)=qw(0)*cin(j)
  !      else
  !         q(0)=zero
  !      end if
  !      qyb(0)=zero
  !      do while (t<tfin)
  !         ! get fluxes
  !         do i=1,n
  !            ! get c and csm=dc/dsm (with theta constant)
  !            k=jt(i)
  !            th=thi(i)+(t-ti)*tht(i)
  !            if (isotype(k,j)=="no" .or. sm(i,j)<zero) then ! handle sm<0 here
  !               csm(i)=one/th
  !               c(i,j)=csm(i)*sm(i,j)
  !            else if (isotype(k,j)=="li") then
  !               csm(i)=one/(th+bd(k)*isopar(k,j)%p(1))
  !               c(i,j)=csm(i)*sm(i,j)
  !            else
  !               do it=1,itmax ! get c from sm using Newton's method and bisection
  !                  if (c(i,j)<zero) c(i,j)=zero ! c and sm are >=0
  !                  call isosub(isotype(k,j), c(i,j), isopar(k,j)%p(:), f, fc)
  !                  csm(i)=one/(th+bd(k)*fc)
  !                  dm=sm(i,j)-(bd(k)*f+th*c(i,j))
  !                  dc=dm*csm(i)
  !                  if (sm(i,j)>=zero .and. c(i,j)+dc<zero) then
  !                     c(i,j)=half*c(i,j)
  !                  else
  !                     c(i,j)=c(i,j)+dc
  !                  end if
  !                  if (abs(dm)<eps*(sm(i,j)+10.0_r_2*dsmmax)) exit
  !                  if (it==itmax) then
  !                     write(*,*) "solute: too many iterations getting c"
  !                     stop
  !                  end if
  !               end do
  !            end if
  !         end do
  !         q(1:n-1)=coef1*c(1:n-1,j)+coef2*c(2:n,j)
  !         qya(1:n-1)=coef1*csm(1:n-1)
  !         qyb(1:n-1)=coef2*csm(2:n)
  !         q(n)=qw(n)*c(n,j)
  !         qya(n)=qw(n)*csm(n)
  !         ! get time step
  !         dmax=maxval(abs(q(1:n)-q(0:n-1))/dx)
  !         if (dmax==zero) then
  !            dt=tfin-t
  !         elseif (dmax<zero) then
  !            write(*,*) "solute: errors in fluxes prevent continuation"
  !            stop
  !         else
  !            dt=dsmmax/dmax
  !         end if
  !         if (t+1.1_r_2*dt>tfin) then
  !            dt=tfin-t
  !            t=tfin
  !         else
  !            t=t+dt
  !         end if
  !         sigdt=sig*dt
  !         rsigdt=one/sigdt
  !         ! adjust q for change in theta
  !         q(1:n-1)=q(1:n-1)-sigdt*(qya(1:n-1)*tht(1:n-1)*c(1:n-1,j)+ &
  !              qyb(1:n-1)*tht(2:n)*c(2:n,j))
  !         q(n)=q(n)-sigdt*qya(n)*tht(n)*c(n,j)
  !         ! get and solve eqns
  !         aa(2:n)=qya(1:n-1)
  !         cc(1:n-1)=-qyb(1:n-1)
  !         bb(1:n)=qyb(0:n-1)-qya(1:n)-dx*rsigdt
  !         dd(1:n)=-(q(0:n-1)-q(1:n))*rsig
  !         call dgtsv(n, 1, aa(2:n), bb(1:n), cc(1:n-1), dd(1:n), n, info)
  !         if (info > 0) then
  !            write(*,*) 'solute: singular matrix (01).'
  !            stop
  !         endif
  !         dy(1:n) = dd(1:n)

  !         ! update unknowns
  !         sdrn(j)=sdrn(j)+(q(n)+sig*qya(n)*dy(n))*dt
  !         sm(:,j)=sm(:,j)+dy(1:n)
  !         nssteps(j)=nssteps(j)+1
  !      end do
  !   end do
  ! END SUBROUTINE solute
  !**********************************************************************************************************************

  SUBROUTINE snow_adjust(irec,mp,n,kk,ns,h0,h0_0,hice,hice_0,thetai,dx,vsnow,var,par,S,Tsoil, &
       Jcol_latent_S, Jcol_latent_T, Jcol_sensible,deltaJ_sensible_S,qmelt)

    INTEGER(i_d),                  INTENT(IN)    :: irec            ! # of grid-cells
    INTEGER(i_d),                  INTENT(IN)    :: mp            ! # of grid-cells
    INTEGER(i_d),                  INTENT(IN)    :: n            ! # of soil layers
    INTEGER(i_d),                  INTENT(IN)    :: kk            ! grid-cell reference
    INTEGER(i_d),   DIMENSION(mp),  INTENT(INOUT)    :: ns           ! pond (0), np ond (1)
    REAL(r_2),   DIMENSION(mp,1:n),   INTENT(IN)    :: dx           ! soil depths
    REAL(r_2),   DIMENSION(mp,1:n),   INTENT(INOUT)    :: Tsoil       ! soil temperatures soil
    REAL(r_2),   DIMENSION(mp,1:n),   INTENT(INOUT)    :: S       ! soil temperatures soil
    REAL(r_2),   DIMENSION(mp),   INTENT(INOUT)    :: h0 ,h0_0, hice , hice_0      ! pond
    REAL(r_2),          DIMENSION(1:mp,1:n),   INTENT(INOUT) :: thetai
    REAL(r_2),          DIMENSION(1:mp),   INTENT(INOUT) :: Jcol_latent_S, Jcol_latent_T, Jcol_sensible
    REAL(r_2),          DIMENSION(1:mp,1:n),   INTENT(INOUT) :: deltaJ_sensible_S
    REAL(r_2),          DIMENSION(1:mp,nsnow_max),   INTENT(INOUT) :: qmelt
    TYPE(vars),   DIMENSION(1:mp,1:n),   INTENT(INOUT)   :: var
    TYPE(params),   DIMENSION(1:mp,1:n),   INTENT(IN)              :: par
    TYPE(vars_snow),     DIMENSION(1:mp),       INTENT(INOUT)           :: vsnow
    REAL(r_2),          DIMENSION(1:mp) :: tmp1d1, tmp1d2, tmp1d3,  tmp1d4
    REAL(r_2),          DIMENSION(1:mp) :: h0_tmp, hice_tmp
    REAL(r_2) :: theta, tmp1, tmp2 ,Tfreezing(1:mp), Jsoil
    INTEGER(i_d) :: i ! counter

    tmp1d1(kk) = h0(kk)+dx(kk,1)*(var(kk,1)%thetai+var(kk,1)%thetal) ! total moisture content of top soil layer + pond
    vsnow(kk)%melt = zero
    ! no dedicated snow pack if solid part of cumulated snow is less than min thresshold
    ! also, don't initialise dedicated snowpack if this would depleted water in top soil layer to less than 1 mm
    if (((vsnow(kk)%wcol-sum(vsnow(kk)%hliq(:)))<snmin*(vsnow(kk)%dens(1)/rhow)).or. &
         ((vsnow(kk)%nsnow==0).and.(tmp1d1(kk)-vsnow(kk)%wcol)<0.001)) then
       vsnow(kk)%nsnow = 0
       !vsnow(kk)%depth = vsnow(kk)%wcol/(vsnow(kk)%dens/rhow)
       theta         = S(kk,1)*(par(kk,1)%thre) + (par(kk,1)%the - par(kk,1)%thre)
       if (vsnow(kk)%hsnow(1)>zero)  then ! termination of dedicated snow layer
          ! total energy in old snow layer
          tmp1d1(kk) = (vsnow(kk)%tsn(1))*rhow*vsnow(kk)%hliq(1)*cswat + &
               (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*((vsnow(kk)%tsn(1))*csice-lambdaf)
          vsnow(kk)%Qadv_transfer(1) = vsnow(kk)%Qadv_transfer(1) - tmp1d1(kk)
          vsnow(kk)%deltaJlatent(1) = vsnow(kk)%deltaJlatent(1) +lambdaf*(vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow
          vsnow(kk)%deltaJsensible(1) = vsnow(kk)%deltaJsensible(1) -(vsnow(kk)%tsn(1))*rhow*vsnow(kk)%hliq(1)*cswat- &
               (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*(vsnow(kk)%tsn(1))*csice
          ! total energy in old top soil layer
          tmp1d2(kk) = var(kk,1)%csoil*dx(kk,1)*(Tsoil(kk,1)) -lambdaf*dx(kk,1)*var(kk,1)%thetai + &
               (h0(kk)-hice_0(kk))*cswat*rhow*(Tsoil(kk,1)) + &
               hice_0(kk)*rhow*(csice*(Tsoil(kk,1))-lambdaf)
          tmp1d3(kk) = Tsoil(kk,1)
          !calculate new thetal, consistent with total energy and new pond height
          h0_tmp(kk) = h0(kk)
          if (h0(kk)>zero) then
             h0(kk) = h0(kk) + vsnow(kk)%hsnow(1)
          elseif ((theta+vsnow(kk)%hsnow(1)/dx(kk,1))<par(kk,1)%the) then
             h0(kk)=zero
             theta = theta+vsnow(kk)%hsnow(1)/dx(kk,1)
             S(kk,1) = (theta - (par(kk,1)%the - par(kk,1)%thre) )/(par(kk,1)%thre)
          elseif ((theta+vsnow(kk)%hsnow(1)/dx(kk,1))>par(kk,1)%the) then
             h0(kk) = ((theta+vsnow(kk)%hsnow(1)/dx(kk,1))-par(kk,1)%the)*dx(kk,1)
             theta = par(kk,1)%the
             S(kk,1) = one
             var(kk,1)%isat = 1
             ns(kk) = 0
          endif

          Tfreezing(kk) = Tfrz(S(kk,1),par(kk,1)%he,one/par(kk,1)%lam)

          ! check if total energy in old snow layer and top soil (+pond) is enough for complete melting
          tmp1d3(kk) = var(kk,1)%thetai*dx(kk,1)*rhow*lambdaf + &
               var(kk,1)%thetai*h0_tmp(kk)/par(kk,1)%thre*rhow*lambdaf
          tmp1d4(kk) = (Tfreezing(kk))*(rhow*cswat*(theta*dx(kk,1)+h0(kk))+par(kk,1)%rho*par(kk,1)%css*dx(kk,1))

          if ((tmp1d1(kk)+tmp1d2(kk)-tmp1d4(kk))>zero) then !  complete melting

             Jcol_latent_T(kk) =     Jcol_latent_T(kk) + tmp1d3(kk)
             deltaJ_sensible_S(kk,1) = zero
             Tsoil(kk,1) = var(kk,1)%Tfrz + (tmp1d1(kk)+tmp1d2(kk) -tmp1d4(kk))/ &
                  (rhow*cswat*(theta*dx(kk,1)+h0(kk))+par(kk,1)%rho*par(kk,1)%css*dx(kk,1))

             var(kk,1)%iice = 0
             var(kk,1)%thetai = zero
             var(kk,1)%thetal = theta
             hice(kk) = zero
             vsnow(kk)%hsnow(1) = zero
             vsnow(kk)%hliq(1)= zero
          else  ! soil remains frozen

             Jsoil = tmp1d1(kk)+tmp1d2(kk)! total energy in  soil layer
             !check there is a zero
             tmp1 = GTfrozen(real(Tsoil(kk,1)-50., r_2), Jsoil, dx(kk,1), theta,par(kk,1)%css, par(kk,1)%rho, &
                  h0(kk), par(kk,1)%thre, par(kk,1)%the, par(kk,1)%he, one/par(kk,1)%lam)

             tmp2 = GTFrozen(Tfreezing(kk), Jsoil, dx(kk,1), theta,par(kk,1)%css, par(kk,1)%rho, &
                  h0(kk), par(kk,1)%thre, par(kk,1)%the, par(kk,1)%he, one/par(kk,1)%lam)

             ! there is a zero in between
             if ((tmp1*tmp2) < zero) then
                tmp1d3(kk) = rtbis_Tfrozen(tmp1d2(kk), dx(kk,1), theta,par(kk,1)%css, par(kk,1)%rho, &
                     h0(kk), par(kk,1)%thre, par(kk,1)%the, par(kk,1)%he, one/par(kk,1)%lam, &
                     real(Tsoil(kk,1)-50., r_2), Tfreezing(kk), real(0.0001,r_2))

                tmp1d4(kk) = thetalmax(tmp1d3(kk), min(S(kk,1),one), par(kk,1)%he, one/par(kk,1)%lam, &
                     par(kk,1)%thre, par(kk,1)%the) ! liquid content at new Tsoil
             else
                write(*,*) "Found no solution for Tfrozen 2. Stop. " , kk
                stop
             endif


             hice_tmp(kk) = hice(kk)
             hice(kk) = h0(kk)*var(kk,1)%thetai/par(kk,1)%thre
             vsnow(kk)%hsnow(1) = zero
             vsnow(kk)%hliq(1) = zero
             var(kk,1)%thetal = tmp1d4(kk)
             var(kk,1)%thetai = theta - tmp1d4(kk)

             ! correct total energy stored in pond + soil
             Jcol_latent_S(kk) = Jcol_latent_S(kk)  - rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
                  dx(kk,1)*(var(kk,1)%thetai-thetai(kk,1)))

             !            deltaJ_latent_S(kk,1) = - rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
             !                 dx(kk,1)*(var(kk,1)%thetai-thetai(kk,1)))

             !             Jcol_sensible(kk) = Jcol_sensible(kk) + &
             !                 dx(kk,1)*tmp1d3(kk)*par(kk,1)%rho*par(kk,1)%css + &
             !                 dx(kk,1)*var(kk,1)%thetal*cswat*rhow*tmp1d3(kk) + &
             !                 (h0(kk)-hice(kk))*cswat*rhow*tmp1d3(kk) + &
             !                 var(kk,1)%thetai*dx(kk,1)*rhow*(csice*(tmp1d3(kk))-lambdaf) + &
             !                 hice(kk)*rhow*(csice*(tmp1d3(kk))-lambdaf)- &
             !                  tmp1d2(kk) - &
             !                - rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
             !                 dx(kk,1)*(var(kk,1)%thetai-thetai(kk,1)))

             Jcol_sensible(kk) = Jcol_sensible(kk) + &
                  JSoilLayer(tmp1d3(kk), &
                  dx(kk,1), theta,par(kk,1)%css, par(kk,1)%rho, &
                  h0(kk), par(kk,1)%thre, par(kk,1)%the, &
                  par(kk,1)%he, one/par(kk,1)%lam) - &
                  tmp1d2(kk) - &
                  (-rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
                  dx(kk,1)*(var(kk,1)%thetai-thetai(kk,1))))




             !             vsnow(kk)%hsnow(1) = zero
             !             vsnow(kk)%hliq(1) = zero
             !             var(kk,1)%thetal = tmp1d4(kk)
             !             var(kk,1)%thetai = theta - tmp1d4(kk)
             !             hice(kk) = h0(kk)*var(kk,1)%thetai/par(kk,1)%thre

             !             ! correct total energy stored in pond + soil
             !             Jcol_latent_S(kk) = Jcol_latent_S(kk)- rhow*lambdaf*((hice(kk)-hice_0(kk)) + &
             !                  dx(kk,1)*(var(kk,1)%thetai-thetai(kk,1)))
             !
             !             Jcol_sensible(kk) = Jcol_sensible(kk) - &
             !                  var(kk,1)%csoil*dx(kk,1)*(Tsoil(kk,1)) - &
             !                  (h0_0(kk)-hice_0(kk))*cswat*rhow*(Tsoil(kk,1)) - &
             !                  hice_0(kk)*rhow*(csice*(Tsoil(kk,1))) + &
             !                  dx(kk,1)*(tmp1d3(kk))*par(kk,1)%rho*par(kk,1)%css + &
             !                  (h0(kk)-hice(kk)+var(kk,1)%thetal*dx(kk,1))*cswat*rhow*(tmp1d3(kk)) + &
             !                  (hice_0(kk)+var(kk,1)%thetai)*rhow*(csice*(tmp1d3(kk)))

             !thetai(kk,1) = var(kk,1)%thetai
             Tsoil(kk,1) = tmp1d3(kk)
             if (var(kk,1)%thetai>zero) then
                var(kk,1)%iice = 1
             endif
          endif

       endif ! termination of dedicated snow layer

    else
       vsnow(kk)%nsnow = 1
       if  (vsnow(kk)%nsnow_last==0) then ! snow layer initialisation (transfer pond water to dedicated snow layer)
          vsnow(kk)%hsnow(1) = vsnow(kk)%wcol

          if (h0(kk)>vsnow(kk)%wcol) then ! extract new snow layer from pond
             ! total energy in new snow layer
             tmp1d1(kk) = (Tsoil(kk,1))*rhow*vsnow(kk)%hsnow(1)*(cswat*(h0(kk)-hice(kk))/h0(kk) + &
                  csice*hice(kk)/h0(kk)) - rhow*lambdaf*vsnow(kk)%hsnow(1)*hice(kk)/h0(kk)
             ! correct total energy stored in pond + soil
             Jcol_latent_S(kk) = Jcol_latent_S(kk) + rhow*lambdaf*vsnow(kk)%hsnow(1)*hice(kk)/h0(kk)
             Jcol_sensible(kk) = Jcol_sensible(kk) - (Tsoil(kk,1))*rhow*vsnow(kk)%hsnow(1) &
                  *(cswat*(h0(kk)-hice(kk))/h0(kk) +  csice*hice(kk)/h0(kk))

             ! total energy in snowpack totally frozen at 0degC
             tmp1d2(kk) = rhow*vsnow(kk)%hsnow(1)*(csice*zero - lambdaf)
             if (tmp1d1(kk)<=tmp1d2(kk)) then
                ! alll snow water frozen
                vsnow(kk)%hliq(1)=zero
                vsnow(kk)%tsn(1) = (tmp1d1(kk)+rhow*lambdaf*vsnow(kk)%hsnow(1))/(csice*rhow*vsnow(kk)%hsnow(1))
             else
                ! liquid snow water
                vsnow(kk)%hliq(1)=(tmp1d1(kk)-vsnow(kk)%hsnow(1)*rhow*(zero*csice-lambdaf))/ &
                     (rhow*(zero*cswat-zero*csice+lambdaf))
                vsnow(kk)%tsn(1) = zero
             endif

             h0(kk) = h0(kk) - vsnow(kk)%wcol
             hice(kk) = h0(kk)*var(kk,1)%thetai/par(kk,1)%thre

          else  ! extract new snow layer from soil ice + pond
             ! total energy in new snow layer (component extracted from soil ice)
             tmp1d1(kk) = (Tsoil(kk,1))*rhow*(vsnow(kk)%hsnow(1)-h0(kk))*csice - &
                  rhow*lambdaf*(vsnow(kk)%hsnow(1)-h0(kk))
             h0_tmp(kk) = h0(kk)
             hice_tmp(kk) = hice(kk)
             if (h0(kk)>zero) then
                tmp1d1(kk) = tmp1d1(kk) + (Tsoil(kk,1))*rhow*(cswat*(h0(kk)-hice(kk)) + &
                     csice*hice(kk)) - rhow*lambdaf*hice(kk)
             endif
             ! total energy in snowpack totally frozen at 0degC
             tmp1d2(kk) = rhow*vsnow(kk)%hsnow(1)*(csice*zero - lambdaf)
             if (tmp1d1(kk)<=tmp1d2(kk)) then
                ! alll snow water frozen
                vsnow(kk)%hliq(1)=zero
                vsnow(kk)%tsn(1) = (tmp1d1(kk)+rhow*lambdaf*vsnow(kk)%hsnow(1))/(csice*rhow*vsnow(kk)%hsnow(1))
             else
                ! liquid snow water
                vsnow(kk)%hliq(1)=(tmp1d1(kk)-vsnow(kk)%hsnow(1)*rhow*(zero*csice-lambdaf))/ &
                     (rhow*(zero*cswat-zero*csice+lambdaf))
                vsnow(kk)%tsn(1) = zero
             endif

             ! correct soil moisture
             S(kk,1) = S(kk,1) - (vsnow(kk)%hsnow(1)-h0(kk))/dx(kk,1)/par(kk,1)%thre
             Tfreezing(kk) = Tfrz(S(kk,1),par(kk,1)%he,one/par(kk,1)%lam)
             if (S(kk,1)<zero) then
                write(*,*) "error: over-extraction of soil water during snow pack init"
                write(*,*) "S(kk,1), snow-col, deltaS"
                write(*,*) S(kk,1) , vsnow(kk)%hsnow(1), - (vsnow(kk)%hsnow(1)-h0(kk))/dx(kk,1)/par(kk,1)%thre
             endif
             h0(kk) = zero
             hice(kk) = zero
             ! correct total energy stored in pond + soil
             ! correct soil temperature
             theta         = S(kk,1)*(par(kk,1)%thre) + (par(kk,1)%the - par(kk,1)%thre)
             ! total energy added to top soil layer
             tmp1d1(kk) = - tmp1d1(kk)
             ! total energy in old top soil layer
             tmp1d2(kk) = var(kk,1)%csoil*dx(kk,1)*(Tsoil(kk,1)) -lambdaf*dx(kk,1)*var(kk,1)%thetai + &
                  (h0_tmp(kk)-hice_tmp(kk))*cswat*rhow*(Tsoil(kk,1)) + &
                  hice_tmp(kk)*rhow*(csice*(Tsoil(kk,1))-lambdaf)
             ! calculate energy in new top soil layer
             if (var(kk,1)%iice==0) then
                var(kk,1)%csoil = theta*rhow*cswat+par(kk,1)%rho*par(kk,1)%css
                tmp1d3(kk) = (tmp1d1(kk)+tmp1d2(kk))/((theta*rhow*cswat+par(kk,1)%rho*par(kk,1)%css)*dx(kk,1)+ &
                     h0(kk)*rhow*cswat)
                Jcol_sensible(kk) = Jcol_sensible(kk) - &
                     var(kk,1)%csoil*dx(kk,1)*(Tsoil(kk,1)) - &
                     (h0_tmp(kk))*cswat*rhow*(Tsoil(kk,1)) + &
                     (theta*rhow*cswat+par(kk,1)%rho*par(kk,1)%css)*dx(kk,1)*(tmp1d3(kk)) + &
                     (h0_tmp(kk))*cswat*rhow*(tmp1d3(kk))

                var(kk,1)%csoil = theta*rhow*cswat+par(kk,1)%rho*par(kk,1)%css
                Tsoil(kk,1) = tmp1d3(kk)

             else
                ! check if total energy in melt water and top soil (+pond) is enough for complete melting
                tmp1d3(kk) = var(kk,1)%thetai*dx(kk,1)*rhow*lambdaf + &
                     var(kk,1)%thetai*h0_tmp(kk)/par(kk,1)%thre*rhow*lambdaf

                tmp1d4(kk) = (Tfreezing(kk))*(rhow*cswat*(theta*dx(kk,1)+h0(kk))+par(kk,1)%rho*par(kk,1)%css*dx(kk,1))
                if ((tmp1d1(kk)+tmp1d2(kk)-tmp1d4(kk))>zero) then
                   !  complete melting
                   Jcol_latent_T(kk) =     Jcol_latent_T(kk) + tmp1d3(kk)
                   deltaJ_sensible_S(kk,1) = zero
                   Tsoil(kk,1) = var(kk,1)%Tfrz + (tmp1d1(kk)+tmp1d2(kk) -tmp1d4(kk))/ &
                        (rhow*cswat*(theta*dx(kk,1)+h0(kk))+par(kk,1)%rho*par(kk,1)%css*dx(kk,1))

                   var(kk,1)%iice = 0
                   var(kk,1)%thetai = zero
                   thetai(kk,1) = var(kk,1)%thetai
                   var(kk,1)%thetal = theta
                else

                   Jsoil = tmp1d1(kk)+tmp1d2(kk)! total energy in  soil layer
                   !check there is a zero
                   tmp1 = GTfrozen(real(Tsoil(kk,1)-50., r_2), Jsoil, dx(kk,1), theta,par(kk,1)%css, par(kk,1)%rho, &
                        h0(kk), par(kk,1)%thre, par(kk,1)%the, par(kk,1)%he, one/par(kk,1)%lam)

                   tmp2 = GTFrozen(Tfreezing(kk), Jsoil, dx(kk,1), theta,par(kk,1)%css, par(kk,1)%rho, &
                        h0(kk), par(kk,1)%thre, par(kk,1)%the, par(kk,1)%he, one/par(kk,1)%lam)

                   ! there is a zero in between
                   if ((tmp1*tmp2) < zero) then
                      tmp1d3(kk) = rtbis_Tfrozen(tmp1d2(kk), dx(kk,1), theta,par(kk,1)%css, par(kk,1)%rho, &
                           h0(kk), par(kk,1)%thre, par(kk,1)%the, par(kk,1)%he, one/par(kk,1)%lam, &
                           real(Tsoil(kk,1)-50., r_2), Tfreezing(kk), real(0.0001,r_2))

                      tmp1d4(kk) = thetalmax(tmp1d3(kk), min(S(kk,1),one), par(kk,1)%he, one/par(kk,1)%lam, &
                           par(kk,1)%thre, par(kk,1)%the) ! liquid content at new Tsoil
                   else
                      write(*,*) "Found no solution for Tfrozen 3. Stop.", kk
                      stop
                   endif

                   var(kk,1)%thetal = tmp1d4(kk)
                   var(kk,1)%thetai = theta - tmp1d4(kk)
                   ! correct total energy stored in pond + soil
                   Jcol_latent_S(kk) = Jcol_latent_S(kk)- rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
                        dx(kk,1)*(var(kk,1)%thetai-thetai(kk,1)))

                   Jcol_sensible(kk) = Jcol_sensible(kk) - &
                        var(kk,1)%csoil*dx(kk,1)*(Tsoil(kk,1)) - &
                        (h0_tmp(kk)-hice_tmp(kk))*cswat*rhow*(Tsoil(kk,1)) - &
                        hice_tmp(kk)*rhow*(csice*(Tsoil(kk,1))) + &
                        dx(kk,1)*(tmp1d3(kk))*par(kk,1)%rho*par(kk,1)%css + &
                        (h0(kk)-hice(kk)+var(kk,1)%thetal*dx(kk,1))*cswat*rhow*(tmp1d3(kk)) + &
                        (hice_tmp(kk)+var(kk,1)%thetai)*rhow*(csice*(tmp1d3(kk)))

                   thetai(kk,1) = var(kk,1)%thetai
                   Tsoil(kk,1) = tmp1d3(kk)

                endif ! incomplete melting

             endif ! iice=1

             vsnow(kk)%deltaJsensible(1) = vsnow(kk)%deltaJsensible(1) + &
                  (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*csice*(vsnow(kk)%tsn(1)) + &
                  vsnow(kk)%hliq(1)*rhow*cswat*(vsnow(kk)%tsn(1))

             vsnow(kk)%deltaJlatent(1) = vsnow(kk)%deltaJlatent(1) + &
                  (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*(-lambdaf)

             vsnow(kk)%Qadv_transfer(1) = vsnow(kk)%Qadv_transfer(1) + &
                  (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*csice*(vsnow(kk)%tsn(1)) + &
                  vsnow(kk)%hliq(1)*rhow*cswat*(vsnow(kk)%tsn(1)) + &
                  (vsnow(kk)%hsnow(1)-vsnow(kk)%hliq(1))*rhow*(-lambdaf)

          endif ! extract new snow layer from soil ice

       endif        ! snow layer initialisation

       !add snow-melt to top soil + pond
       if ((vsnow(kk)%hliq(1) - vsnow(kk)%hsnow(1)*fsnowliq_max)>zero) then
          h0_tmp(kk) = h0(kk)
          hice_tmp(kk) = hice(kk)
          qmelt(kk,1) = max((vsnow(kk)%hliq(1) - vsnow(kk)%hsnow(1)*fsnowliq_max),zero)
          qmelt(kk,1) = min(qmelt(kk,1), max(0.9*(vsnow(kk)%hsnow(1)-snmin*(vsnow(kk)%dens(1)/rhow)),zero))
          vsnow(kk)%melt(1) = vsnow(kk)%melt(1) + qmelt(kk,1)
          vsnow(kk)%hliq(1) = vsnow(kk)%hliq(1) - qmelt(kk,1)
          vsnow(kk)%hsnow(1) = vsnow(kk)%hsnow(1) - qmelt(kk,1)
          ! total energy in melt-water
          tmp1d1(kk) = (vsnow(kk)%tsn(1))*rhow*qmelt(kk,1)*cswat
          vsnow(kk)%Qadv_melt(1) = vsnow(kk)%Qadv_melt(1)-tmp1d1(kk)
          vsnow(kk)%deltaJsensible(1) = vsnow(kk)%deltaJsensible(1) - tmp1d1(kk)
          ! total energy in old top soil layer
          !          tmp1d2(kk) = var(kk,1)%csoil*dx(kk,1)*(Tsoil(kk,1)) -lambdaf*dx(kk,1)*var(kk,1)%thetai + &
          !               (h0_tmp(kk)-hice_tmp(kk))*cswat*rhow*(Tsoil(kk,1)) + &
          !               hice_tmp(kk)*rhow*(csice*(Tsoil(kk,1))-lambdaf)

          tmp1d2(kk) = dx(kk,1)*par(kk,1)%rho*par(kk,1)%css*Tsoil(kk,1) + &
               dx(kk,1)*rhow*cswat*Tsoil(kk,1)*var(kk,1)%thetal + &
               dx(kk,1)*rhow*var(kk,1)%thetai*(csice*Tsoil(kk,1) - lambdaf) + &
               (h0_tmp(kk)-hice_tmp(kk))*cswat*rhow*(Tsoil(kk,1)) + &
               hice_tmp(kk)*rhow*(csice*(Tsoil(kk,1))-lambdaf)
          !calculate new thetal, consistent with total energy and new pond height
          theta         = S(kk,1)*(par(kk,1)%thre) + (par(kk,1)%the - par(kk,1)%thre)
          if (h0(kk)>zero) then
             h0(kk) = h0(kk) +  qmelt(kk,1)
          elseif ((theta+qmelt(kk,1)/dx(kk,1))<par(kk,1)%the) then
             h0(kk)=zero
             theta = theta+qmelt(kk,1)/dx(kk,1)
             S(kk,1) = (theta - (par(kk,1)%the - par(kk,1)%thre) )/(par(kk,1)%thre)
          elseif ((theta+qmelt(kk,1)/dx(kk,1))>par(kk,1)%the) then
             h0(kk) = ((theta+qmelt(kk,1)/dx(kk,1))-par(kk,1)%the)*dx(kk,1)
             theta = par(kk,1)%the
             S(kk,1) = one
             var(kk,1)%isat = 1
             ns(kk) = 0
          endif
          Tfreezing(kk) = Tfrz(S(kk,1),par(kk,1)%he,one/par(kk,1)%lam)

          ! calculate energy in new top soil layer
          if (var(kk,1)%iice==0) then
             var(kk,1)%csoil = theta*rhow*cswat+par(kk,1)%rho*par(kk,1)%css
             tmp1d3(kk) = (tmp1d1(kk)+tmp1d2(kk))/((theta*rhow*cswat+par(kk,1)%rho*par(kk,1)%css)*dx(kk,1)+ &
                  h0(kk)*rhow*cswat)
             Jcol_sensible(kk) = Jcol_sensible(kk) - &
                  var(kk,1)%csoil*dx(kk,1)*(Tsoil(kk,1)) - &
                  (h0_tmp(kk))*cswat*rhow*(Tsoil(kk,1)) + &
                  (theta*rhow*cswat+par(kk,1)%rho*par(kk,1)%css)*dx(kk,1)*(tmp1d3(kk)) + &
                  (h0_tmp(kk))*cswat*rhow*(tmp1d3(kk))

             var(kk,1)%csoil = theta*rhow*cswat+par(kk,1)%rho*par(kk,1)%css
             Tsoil(kk,1) = tmp1d3(kk)

          else
             ! check if total energy in melt water and top soil (+pond) is enough for complete melting
             tmp1d3(kk) = var(kk,1)%thetai*dx(kk,1)*rhow*lambdaf + &
                  var(kk,1)%thetai*h0_tmp(kk)/par(kk,1)%thre*rhow*lambdaf
             tmp1d4(kk) = (Tfreezing(kk))*(rhow*cswat*(theta*dx(kk,1)+h0(kk))+par(kk,1)%rho*par(kk,1)%css*dx(kk,1))
             if ((tmp1d1(kk)+tmp1d2(kk)-tmp1d4(kk))>zero) then
                !  complete melting
                Jcol_latent_T(kk) =     Jcol_latent_T(kk) + tmp1d3(kk)
                deltaJ_sensible_S(kk,1) = zero
                Tsoil(kk,1) = var(kk,1)%Tfrz + (tmp1d1(kk)+tmp1d2(kk) -tmp1d4(kk))/ &
                     (rhow*cswat*(theta*dx(kk,1)+h0(kk))+par(kk,1)%rho*par(kk,1)%css*dx(kk,1))

                var(kk,1)%iice = 0
                var(kk,1)%thetai = zero
                var(kk,1)%thetal = theta
             else

                Jsoil = tmp1d1(kk)+tmp1d2(kk)! total energy in  soil layer
                !check there is a zero
                tmp1 = GTfrozen(real(Tsoil(kk,1)-50., r_2), Jsoil, dx(kk,1), theta,par(kk,1)%css, par(kk,1)%rho, &
                     h0(kk), par(kk,1)%thre, par(kk,1)%the, par(kk,1)%he, one/par(kk,1)%lam)

                tmp2 = GTFrozen(Tfreezing(kk), Jsoil, dx(kk,1), theta,par(kk,1)%css, par(kk,1)%rho, &
                     h0(kk), par(kk,1)%thre, par(kk,1)%the, par(kk,1)%he, one/par(kk,1)%lam)

                ! there is a zero in between
                if ((tmp1*tmp2) < zero) then
                   tmp1d3(kk) = rtbis_Tfrozen(tmp1d2(kk), dx(kk,1), theta,par(kk,1)%css, par(kk,1)%rho, &
                        h0(kk), par(kk,1)%thre, par(kk,1)%the, par(kk,1)%he, one/par(kk,1)%lam, &
                        real(Tsoil(kk,1)-50., r_2), Tfreezing(kk), real(0.0001,r_2))

                   tmp1d4(kk) = thetalmax(tmp1d3(kk), min(S(kk,1),one), par(kk,1)%he, one/par(kk,1)%lam, &
                        par(kk,1)%thre, par(kk,1)%the) ! liquid content at new Tsoil
                else
                   write(*,*) "Found no solution for Tfrozen 4. Stop.", irec, qmelt(kk,1), h0(kk)
                   stop
                endif

                var(kk,1)%thetal = tmp1d4(kk)
                var(kk,1)%thetai = theta - tmp1d4(kk)
                hice_tmp(kk) = hice(kk)
                hice = h0(kk)*var(kk,1)%thetai/par(kk,1)%thre

                ! correct total energy stored in pond + soil
                Jcol_latent_S(kk) = Jcol_latent_S(kk)- rhow*lambdaf*((hice(kk)-hice_tmp(kk)) + &
                     dx(kk,1)*(var(kk,1)%thetai-thetai(kk,1)))

                Jcol_sensible(kk) = Jcol_sensible(kk) - &
                     var(kk,1)%csoil*dx(kk,1)*(Tsoil(kk,1)) - &
                     (h0_tmp(kk)-hice_tmp(kk))*cswat*rhow*(Tsoil(kk,1)) - &
                     hice_tmp(kk)*rhow*(csice*(Tsoil(kk,1))) + &
                     dx(kk,1)*(tmp1d3(kk))*par(kk,1)%rho*par(kk,1)%css + &
                     (h0(kk)-hice(kk)+var(kk,1)%thetal*dx(kk,1))*cswat*rhow*(tmp1d3(kk)) + &
                     (hice_tmp(kk)+var(kk,1)%thetai)*rhow*(csice*(tmp1d3(kk)))

                thetai(kk,1) = var(kk,1)%thetai
                Tsoil(kk,1) = tmp1d3(kk)

             endif ! incomplete melting
             if (S(kk,1)<zero) then
                write(*,*) "error: over-extraction of soil water during snow pack init"
                write(*,*) "S(kk,1), snow-col, qmelt"
                write(*,*) S(kk,1) , vsnow(kk)%hsnow(1), qmelt(kk,1)
             endif
          endif ! iice=1

       endif ! melt water

       do i=1,vsnow(kk)%nsnow
          vsnow(kk)%Dv(i) = Dva*((vsnow(kk)%tsn(i)+Tzero)/Tzero)**1.88_r_2 ! m2 s-1
          vsnow(kk)%sl(i) = slope_esat_ice(vsnow(kk)%tsn(i)) * Mw/thousand/Rgas/(vsnow(kk)%tsn(i)+Tzero)
          vsnow(kk)%kE(i)     = vsnow(kk)%Dv(i)*vsnow(kk)%sl(i)*thousand*lambdaf
          vsnow(kk)%kth(i) = vsnow(kk)%kE(i) + vsnow(kk)%kH(i)
          vsnow(kk)%cv(i) = esat_ice(vsnow(kk)%tsn(i))*Mw/thousand/Rgas/(vsnow(kk)%tsn(i)+Tzero) ! m3 m-3
          vsnow(kk)%depth(i) = vsnow(kk)%hsnow(i)/(vsnow(kk)%dens(i)/rhow)
       enddo

    endif ! dedicated snow layer

  END SUBROUTINE snow_adjust

  !**********************************************************************************************************************

  SUBROUTINE isotope_vap(isotopologue, litter, n, & ! scalar in
       ns, dx, deltaz, sig, dt, dxL, &     ! in soil
       Tsoil0, dTsoil, dT0, Sliqice, deltaSliqice, & ! in variables
       Sliq, deltaSliq,Sice, deltaSice, &
       Ts, TL, T0, h0new, dh0, SLliq, deltaSLliq, &
       qsig, qlsig, qvsig, qprec, qevap, qrunoff, qex, qd, & ! in fluxes
       var, thetasat, thetar, tortuosity, deltacv, ram, rbh, vlit, cva, civa, & ! in parameter
       thetasatL, deltacvL, &
       cprec, cali, & ! in iso
       ql0, qv0, & ! inout water
       ciso,cisoice, ciso0, cisoL, cisos, & ! inout iso
       qiso_in, qiso_out, qiso_evap, qiso_trans, qiso_liq_adv, qiso_vap_adv, qiso_liq_diff, qiso_vap_diff) ! out iso

    IMPLICIT NONE

    INTEGER(i_d),                  INTENT(IN)    :: isotopologue ! which isotope
    LOGICAL,                       INTENT(IN)    :: litter       ! litter or not
    INTEGER(i_d),                  INTENT(IN)    :: n            ! # of soil layers
    INTEGER(i_d),                  INTENT(IN)    :: ns           ! pond or litter (0), nothing (1)
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: dx           ! soil depths
    REAL(r_2),   DIMENSION(1:n-1), INTENT(IN)    :: deltaz       ! soil layer thickness
    REAL(r_2),                     INTENT(IN)    :: sig          ! implicit/explicit time steping constant
    REAL(r_2),                     INTENT(IN)    :: dt           ! time step
    REAL(r_2),                     INTENT(IN)    :: dxL          ! litter layer thickness
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: Tsoil0       ! soil temperatures soil
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: dTsoil       ! soil temperatures change
    REAL(r_2),                     INTENT(IN)    :: dT0          ! soil temperatures change
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: Sliqice         ! soil saturation
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: deltaSliqice    ! soil sturation change
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: Sliq         ! soil saturation
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: deltaSliq    ! soil sturation change
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: Sice         ! soil saturation
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: deltaSice    ! soil sturation change
    REAL(r_2),                     INTENT(IN)    :: Ts           ! surface temperature
    REAL(r_2),                     INTENT(IN)    :: TL           ! litter/pond temperature
    REAL(r_2),                     INTENT(IN)    :: T0           ! air temperature
    REAL(r_2),                     INTENT(IN)    :: h0new        ! pond height
    REAL(r_2),                     INTENT(IN)    :: dh0          ! pond height change
    REAL(r_2),                     INTENT(IN)    :: SLliq        ! litter saturation
    REAL(r_2),                     INTENT(IN)    :: deltaSLliq   ! litter saturation change
    REAL(r_2),   DIMENSION(0:n),   INTENT(IN)    :: qsig         ! water flux
    REAL(r_2),   DIMENSION(0:n),   INTENT(IN)    :: qlsig        ! liquid water flux
    REAL(r_2),   DIMENSION(0:n),   INTENT(IN)    :: qvsig        ! vapour flux
    REAL(r_2),                     INTENT(IN)    :: qprec        ! precip
    REAL(r_2),                     INTENT(IN)    :: qevap        ! evaporation
    REAL(r_2),                     INTENT(IN)    :: qrunoff      ! runoff
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: qex          ! root extraction
    REAL(r_2),                     INTENT(IN)    :: qd           ! litter to soil drainage
    TYPE(vars),  DIMENSION(1:n),   INTENT(IN)    :: var          ! soil variables
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: thetasat     ! saturation moisture
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: thetar       ! residual moisture
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: tortuosity   ! soil tortuosity
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: deltacv      !
    REAL(r_2),                     INTENT(IN)    :: ram          ! aerodynamic resistance
    REAL(r_2),                     INTENT(IN)    :: rbh          ! boundary layer resistance
    TYPE(vars),                    INTENT(IN)    :: vlit         ! litter variables
    REAL(r_2),                     INTENT(IN)    :: cva          !
    REAL(r_2),                     INTENT(IN)    :: civa         !
    REAL(r_2),                     INTENT(IN)    :: thetasatL    ! litter saturation moisture
    REAL(r_2),                     INTENT(IN)    :: deltacvL     !
    REAL(r_2),                     INTENT(IN)    :: cprec        ! iso conc in precip
    REAL(r_2),                     INTENT(IN)    :: cali         ! iso conc in alimenation water (from below)
    REAL(r_2),                     INTENT(INOUT) :: ql0          ! liquid flux into soil
    REAL(r_2),                     INTENT(INOUT) :: qv0          ! vapour flux into soil
    REAL(r_2),   DIMENSION(0:n),   INTENT(INOUT) :: ciso         ! iso conc in soil
    REAL(r_2),   DIMENSION(0:n),   INTENT(INOUT) :: cisoice      ! iso conc in soil
    REAL(r_2),                     INTENT(INOUT) :: ciso0        ! iso conc in air
    REAL(r_2),                     INTENT(INOUT) :: cisoL        ! iso conc in litter/pond
    REAL(r_2),                     INTENT(INOUT) :: cisos        ! iso conc on surface
    REAL(r_2),                     INTENT(OUT)   :: qiso_in      ! iso flux into soil
    REAL(r_2),                     INTENT(OUT)   :: qiso_out     ! iso flux out of soil
    REAL(r_2),                     INTENT(OUT)   :: qiso_evap    ! iso flux of evaporation
    REAL(r_2),                     INTENT(OUT)   :: qiso_trans   ! iso flux of transpiration
    REAL(r_2),   DIMENSION(1:n),   INTENT(OUT)   :: qiso_liq_adv ! liquid iso flux in soil due to advection
    REAL(r_2),   DIMENSION(1:n),   INTENT(OUT)   :: qiso_vap_adv ! vapour iso flux in soil due to advection
    REAL(r_2),   DIMENSION(1:n-1), INTENT(OUT)   :: qiso_liq_diff ! liquid iso flux in soil due to diffusion
    REAL(r_2),   DIMENSION(1:n-1), INTENT(OUT)   :: qiso_vap_diff ! vapour iso flux in soil due to diffusion

    ! Local variables

    REAL(r_2), DIMENSION(0:n)   :: aa, bb, cc, dd, dc,  LHS, RHS
    REAL(r_2), DIMENSION(0:n)   :: alphaplus, dalphaplusdT, deltaT
    ! diffusivities (liquid, liq-vap, coefft for D, surface H2O vapour, surface minor isotopologue vapour)
    REAL(r_2), DIMENSION(0:n)   :: Dl, Dv
    REAL(r_2)                   :: Dvs, Divs
    REAL(r_2)                   :: patm, nk, alphak, alphak_vdiff, alphak_ldiff
    REAL(r_2)                   :: cevapin, cevapout, qevapin, qevapout, dcevapoutdciso
    ! concentrations of advective fluxes and corresponding partial derivs wrt ciso
    REAL(r_2), DIMENSION(0:n)   :: cql, dcqldca, dcqldcb, cqv, dcqvdca, dcqvdcb
    REAL(r_2), DIMENSION(0:n-1) :: wcql, wcqv
    REAL(r_2), DIMENSION(0:n)   :: beta, deltabeta, betaqv, dbetaqv
    REAL(r_2), DIMENSION(0:n)   :: Dlmean, Dvmean, wl, wv
    REAL(r_2)                   :: coefA, coefB, coefC
    REAL(r_2), DIMENSION(0:n)   :: Seff, deltaSeff, S, Tsoil, cvsig, Sliqsig, Sicesig
    REAL(r_2), DIMENSION(0:n)   :: thetaice, deltathetaice, dcice
    REAL(r_2)                   :: h0
    INTEGER(i_d)                :: ns_ciso, j
    REAL(r_2)                   :: num, den, cvL, cv1
    REAL(r_2)                   :: alphaplus_s, alphaplus_0
    REAL(r_2)                   :: alphaplus_liqice
    REAL(r_2)                   :: deltaz0, qevapL, cv0,cvs, qevapoutL, qevapinL
    REAL(r_2)                   :: cevapinL, cevapoutL, dcevapoutdcisoL, dcevapindcisoL
    REAL(r_2), DIMENSION(0:n)   :: Dveff
    REAL(r_2)                   :: w1
    INTEGER(i_d), PARAMETER     :: formulation= 2
    REAL(r_2),    PARAMETER     :: qmin=1.e-15
    REAL(r_2), DIMENSION(1:n)   :: ifreeze         ! ==1 if ice accumulating (deltathetice>0)
    REAL(r_2), DIMENSION(1:n)   :: kfreeze, kfreeze2         ! combination of freezing variables (used when ifreeze=1)

    dcice = zero
    aa(0) = zero
    bb(0) = zero
    cc(0) = zero
    dd(0) = zero
    ns_ciso = ns ! index of top layer (0 (pond or litter) or 1 (soil))
    if ( litter) ns_ciso = 0
    if (litter .and. (ns==1)) ciso(0) = cisoL
    deltaz0 = half*dxL + half*dx(1)
    patm    = one
    deltaT(1:n) = dTsoil(1:n)
    deltaT(0)   = dT0

    Tsoil(1:n) = Tsoil0(1:n) + (sig-one)*deltaT(1:n)
    Tsoil(0)   = TL + (sig-one)*deltaT(0)         ! litter or pond temperature

    ! equilibrium fractionation factors at the surface and in the soil
    coefA = zero
    coefB = zero
    coefC = zero
    if (isotopologue==1) then
       coefA = 24844.0_r_2
       coefB = -76.248_r_2
       coefC = 0.052612_r_2
       alphaplus_liqice = 1.0212_r_2
    endif
    if (isotopologue==2) then
       coefA = 1137.0_r_2
       coefB = -0.4156_r_2
       coefC = -0.0020667_r_2
       alphaplus_liqice = 1.00291_r_2
    endif
    alphaplus_s       = one/exp(coefA/((Ts+Tzero)**2)+coefB/(Ts+Tzero)+coefC)       ! at soil or litter surface
    alphaplus_0       = one/exp(coefA/((T0+Tzero)**2)+coefB/(T0+Tzero)+coefC)       ! at soil/litter interface
    alphaplus(0:n)    = one/exp(coefA/((Tsoil(0:n)+Tzero)**2)+coefB/(Tsoil(0:n)+Tzero)+coefC)
    dalphaplusdT(0:n) = (two*coefA/(Tsoil(0:n)+Tzero)**3 + coefB/(Tsoil(0:n)+Tzero)**2) &
         / exp(coefA/((Tsoil(0:n)+Tzero)**2)+coefB/(Tsoil(0:n)+Tzero)+coefC)
    if (experiment==1 .or. experiment==2.or. experiment==9.or. experiment==10) then
       alphaplus_s    = one
       alphaplus      = one
       dalphaplusdT   = zero
       alphaplus_0    = one
       alphaplus_liqice = one
    endif

    !beta = cv/ cl
    beta(0:n) = alphaplus(0:n)      ! dimensionless
    if (ns==0 ) beta(0) = zero
    ! delta beta
    deltabeta(0:n) =  dalphaplusdT(0:n)*deltaT(0:n)
    if (ns==0) deltabeta(0) = zero
    beta = beta + sig*deltabeta         !beta_sig

    ! adjust S and h and cv for sig of time-step
    if (litter) then
       S(0)         = SLliq + deltaSLliq*(sig-one)
       cvsig(0)     =  vlit%cv + deltacvL*(sig-one)
       Seff(0)      = (S(0) + sig*deltaSLliq) + (cvsig(0)*beta(0) + sig*(beta(0)*deltacvL + cvsig(0)*deltabeta(0))) &
            - ( cvsig(0)*S(0)*beta(0) + sig*(S(0)*beta(0)*deltacvL &
            + cvsig(0)*beta(0)*deltaSLliq + cvsig(0)*S(0)*deltabeta(0)))

       deltaSeff(0) = deltaSLliq + beta(0)*deltacvL + cvsig(0)*deltabeta(0) &
            - (S(0)*beta(0)*deltacvL + cvsig(0)*beta(0)*deltaSLliq + cvsig(0)*S(0)*deltabeta(0))
       Sliqsig(0) = zero
       Sicesig(0) = zero
       thetaice(0) = zero
       deltathetaice(0) = zero
    else
       S(0)         = one
       cvsig(0)     = zero
       Seff(0)      = one
       deltaSeff(0) = one
       Sliqsig(0) = zero
       Sicesig(0) = zero
       thetaice(0) = zero
       deltathetaice(0) = zero
    endif

    S(1:n)     = Sliqice(1:n) + deltaSliqice(1:n)*(sig-one)      ! set S to Ssig. N.B. S = Sliqice
    Sliqsig(1:n)     = Sliq(1:n) + deltaSliq(1:n)*(sig-one)
    Sicesig(1:n)     = Sice(1:n) + deltaSice(1:n)*(sig-one)
    thetaice(1:n) = Sice(1:n) * (thetasat(1:n)-thetar(1:n))
    deltathetaice(1:n) = deltaSice(1:n) * (thetasat(1:n)-thetar(1:n))
    ifreeze(1:n) = merge(one,zero,deltathetaice(1:n)>zero)

    where (deltathetaice(1:n)>zero)
       kfreeze(1:n) = alphaplus_liqice*deltaSice(1:n)/Sice(1:n)
       kfreeze2(1:n) = deltaSice(1:n)/Sice(1:n)*(alphaplus_liqice*ciso(1:n)-cisoice(1:n))
    elsewhere
       kfreeze(1:n) = zero
       kfreeze2(1:n) = zero
    endwhere

    if ((Sice(1)<=zero) .and. (deltaSice(1)>zero)) then
       write(*,*) "deltaSice inconistent with Sice"
    endif

    h0         = h0new +dh0*(sig-one)    ! set h0 to h0sig
    cvsig(1:n) = var(1:n)%cv + deltacv(1:n)*(sig-one)

    where (S(1:n)<one)
       Seff(1:n) = Sliqsig(1:n) + cvsig(1:n)*beta(1:n)- cvsig(1:n)*S(1:n)*beta(1:n) &
            + thetar(1:n)/thetasat(1:n) !+ &


       deltaSeff(1:n) = deltaSliq(1:n) + beta(1:n)*deltacv(1:n) + cvsig(1:n)*deltabeta(1:n) &
            - (S(1:n)*beta(1:n)*deltacv(1:n) + cvsig(1:n)*beta(1:n)*deltaSliqice(1:n) + cvsig(1:n)*S(1:n)*deltabeta(1:n)) !+ &
    elsewhere
       Seff(1:n) = Sliqsig(1:n) + thetar(1:n)/thetasat(1:n)

       deltaSeff(1:n) = deltaSliq(1:n)
    endwhere

    ! diffusional fractionation factor in air
    alphak_vdiff =  one
    if (isotopologue==1) alphak_vdiff = one / 1.0251_r_2   ! HDO diffusivity in air (Merlivat 1978)
    if (isotopologue==2) alphak_vdiff = one / 1.0285_r_2   ! H218O diffusivity in air (Merlivat 1978)
    if (experiment >= 1 .and. experiment <= 5) alphak_vdiff = one
    if (experiment == 9 .or.experiment == 10) alphak_vdiff = one

    ! kinetic fractionation factor at the surface
    Dvs  = Dva*1.e5_r_2/patm*((Ts+Tzero)/Tzero)**1.88_r_2 ! vapour diffuxivity of water in air (m2s-1)
    Divs =  Dvs
    if (isotopologue==1) Divs = Dvs * alphak_vdiff  ! HDO diffusivity in air
    if (isotopologue==2) Divs = Dvs * alphak_vdiff
    nk = ((thetasat(1)*min(S(1),one)-zero)*half + (thetasat(1)*(one-min(S(1),one))))/(thetasat(1)-zero)
    if (experiment==7 .or. experiment==8) nk = one
    alphak = one/((Dvs/Divs)**nk)
    !alphak = one/(((Dvs/Divs)**nk + ram/rbh)/(one + ram/rbh)) ! kinetic fractionation factor (< 1)
    if (experiment >= 1 .and. experiment <= 5) alphak = one
    if (experiment == 9 .or.experiment == 10) alphak = one

    ! liquid diffusivity in the pond and soil
    alphak_ldiff = one
    if (isotopologue==1) alphak_ldiff = one / 1.013_r_2
    if (isotopologue==2) alphak_ldiff = one / 1.026_r_2
    !molecular diffusion of HDO in normal liquid water (m2s-1)
    Dl(1:n) = tortuosity(1:n) * alphak_ldiff * 1.0e-7_r_2*exp(-577.0_r_2/((Tsoil(1:n)+Tzero)-145._r_2))
    !molecular diffusion of HDO in normal liquid water (m2s-1) (pond)
    Dl(0)   = alphak_ldiff * 1.0e-7_r_2*exp(-577.0_r_2/((T0+Tzero)-145._r_2))
    Dl(1:n) = Dl(1:n) * (min(S(1:n),one) * (thetasat(1:n)-thetar(1:n)) + thetar(1:n))
    if (experiment >= 1 .and. experiment <= 4) Dl = zero
    if (experiment == 9 .or.experiment == 10) Dl = zero

    ! vapour diffusivity in the soil
    Dv(1:n) = var%Dv * alphak_vdiff  ! isotope diffusivity in soil air spaces
    if (experiment >= 1 .and. experiment <= 5) Dv(1:n) = var%Dv
    if (experiment == 9 .or.experiment == 10) Dv(1:n) = var%Dv
    Dv(0)   = zero
    Dv(1:n) = Dv(1:n) * cvsig(1:n) !* thetasat(1:n) * (1. - S(1:n))

    do j=1,n-1
       if (abs(cvsig(j)-cvsig(j+1)) > 1.e-8) then
          Dveff(j) = qvsig(j)/(cvsig(j)-cvsig(j+1))*deltaz(j) * alphak_vdiff
          if (experiment >= 1 .and. experiment <= 5) then
             Dveff(j) = qvsig(j)/(cvsig(j)-cvsig(j+1))*deltaz(j)
          endif
          if (experiment == 9 .or.experiment == 10) then
             Dveff(j) = qvsig(j)/(cvsig(j)-cvsig(j+1))*deltaz(j)
          endif
       else
          Dveff(j) = var(j)%Dv*alphak_vdiff
          if (experiment >= 1 .and. experiment <= 5) then
             Dveff(j) = var(j)%Dv
          endif
          if (experiment == 9 .or.experiment == 10) then
             Dveff(j) = var(j)%Dv
          endif

       endif

       if (Dveff(j) <= zero) then
          Dveff(j) = var(j)%Dv*alphak_vdiff
          if (experiment >= 1 .and. experiment <= 5) then
             Dveff(j) = var(j)%Dv
          endif
          if (experiment == 9 .or.experiment == 10) then
             Dveff(j) = var(j)%Dv
          endif
       endif
    enddo

    if (litter .and. (ns==1)) then
       Dv(0) = vlit%Dv*alphak_vdiff
       if (experiment >= 1 .and. experiment <= 5) Dv(0) = vlit%Dv
       if (experiment == 9 .or.experiment == 10) Dv(0) = vlit%Dv
       Dv(0) = Dv(0) * cvsig(0)
    endif

    ! upper boundary condition
    cvs = cva + qevap*(ram+rbh) ! concentration of water vapour at soil/air interface (m3 (H2O liq)/ m3 (air))

    ! denominator can be zero, check
    if ((cvs-(var(1)%cv-deltacv(1))) /= zero) then
       Dveff(0) = qv0/(cvs-(var(1)%cv-deltacv(1)))*(half*dx(1))
    else
       ! What should it be in this case?
       Dveff(0) = var(1)%Dv
    endif

    if (var(1)%Dv /= zero) then
       cv1 = -qv0*(half*dx(1))/var(1)%Dv + cvs
    else
       cv1 = var(1)%cv
    endif

    ! do as before
    cvL = zero
    if (litter) then
       Dveff(0) = -qevap/(cvs-vlit%cv)*(half*dxL)
       if (vlit%Dv /= zero) then
          cvL = qevap*(half*dxL)/vlit%Dv + cvs
       else
          cvL = vlit%cv
       endif
    endif

    if (ql0 > zero) then
       w1 = zero
    else
       w1 = one
    endif

    if (litter .and. (ns==1)) then
       num            = alphak_vdiff*vlit%Dv/(half*dxL)*cvL*alphaplus(0)*ciso(0)  +civa*alphak/(ram+rbh)
       den            = alphak*cvs*alphaplus_s/(ram+rbh) + alphaplus(0)*alphak_vdiff*cvs*vlit%Dv/(half*dxL)
       cisos          = num/den
       dcevapoutdciso = alphak*alphaplus_s*alphak_vdiff*vlit%Dv/(half*dxL)*cvL*alphaplus(0)/den

       qevapL = -(qsig(0) - qd)                     ! vapour flux from top soil layer to litter
       cv0    = qevapL*(half*dxL)/vlit%Dv + (vlit%cv -deltacvL)   ! vapour conc at soil/liter interface
       num    = alphak_vdiff*vlit%Dv/(half*dxL)*(vlit%cv -deltacvL)*alphaplus(0)*ciso(0) &
            +alphak_vdiff*var(1)%Dv/(half*dx(1))*(var(1)%cv-deltacv(1))*alphaplus(1)*ciso(1) &
            - ql0*ciso(1)*w1  + Dl(1)*ciso(1)/(half*dx(1))

       den             = alphak_vdiff*vlit%Dv/(half*dxL)*cv0*alphaplus(0)  + ql0*(one-w1) + &
            alphaplus(1)*alphak_vdiff*cv0*var(1)%Dv/(half*dx(1)) +Dl(1)/(half*dx(1))
       ciso0           = num/den
       dcevapoutdcisoL = alphak*alphaplus_0/den*alphak_vdiff*vlit%Dv/(half*dxL)*vlit%cv*alphaplus(0)*ciso(0)


    elseif ((.not. litter) .and. (ns==1)) then
       Dveff(0) = var(1)%Dv
       num            = alphak_vdiff*Dveff(0)/(half*dx(1))*cv1*alphaplus(1)*ciso(1) &
            - ql0*ciso(1)*w1 +civa*alphak/(ram+rbh) + Dl(1)*ciso(1)/(half*dx(1))
       den            = alphak*cvs*alphaplus_s/(ram+rbh) + ql0*(one-w1) + &
            alphaplus(1)*alphak_vdiff*cvs*Dveff(0)/(half*dx(1)) +Dl(1)/(half*dx(1))
       cisos          = num/den
       dcevapoutdciso = alphak*alphaplus_s
       dcevapoutdciso = dcevapoutdciso*(alphak_vdiff*Dveff(0)/(half*dx(1))*var(1)%cv*alphaplus(1) - &
            ql0*w1  + Dl(1)/(half*dx(1)))/den
    else
       cisos          = ciso(0)
       num            = zero
       den            = zero
       dcevapoutdciso = alphak*alphaplus_s
    endif

    qevapin  = cva/(ram+rbh)
    qevapout = cvs/(ram+rbh)
    cevapin  = civa/cva * alphak
    cevapout = alphak*alphaplus_s * cisos

    qevapL          = zero
    qevapoutL       = zero
    qevapinL        = zero
    cevapinL        = zero
    cevapoutL       = zero
    dcevapoutdcisoL = zero
    dcevapindcisoL  = zero
    if (litter  .and. (ns==1)) then
       qevapL          = -(qsig(0) - qd)                        ! vapour flux from top soil layer to litter
       cv0             = qevapL*(half*dxL)/vlit%Dv + vlit%cv    ! vapour conc at soil/liter interface
       qevapoutL       = cv0*vlit%Dv/(half*dxL)
       qevapinL        = vlit%cv*vlit%Dv/(half*dxL)
       cevapinL        = alphak*alphaplus(0)*ciso(0)
       cevapoutL       = alphak*alphaplus_0*ciso(1)
       dcevapoutdcisoL = alphak*alphaplus_0*ciso0
       dcevapindcisoL  = alphak*alphaplus(0)
    endif

    ! concentrations of advective fluxes and corresponding partial derivs wrt ciso
    select case (formulation)
    case (1)
       cql(0:n-1)     = merge(ciso(0:n-1), ciso(1:n), qlsig(0:n-1)>zero)
       dcqldca(0:n-1) = merge(one,  zero, qlsig(0:n-1)>zero)
       dcqldcb(0:n-1) = merge(zero, one,  qlsig(0:n-1)>zero)
       cql(n)         = ciso(n)
       dcqldca(n)     = one
       dcqldcb(n)     = zero

       cqv(0:n-1)     = merge(ciso(0:n-1), ciso(1:n), qvsig(0:n-1)>zero)
       dcqvdca(0:n-1) = merge(one,  zero, qvsig(0:n-1)>zero)
       dcqvdcb(0:n-1) = merge(zero, one,  qvsig(0:n-1)>zero)
       cqv(n)         = ciso(n)
       dcqvdca(n)     = one
       dcqvdcb(n)     = zero

       betaqv(0:n-1)  = merge(beta(0:n-1), beta(1:n), qvsig(0:n-1)>0) * alphak_vdiff
       betaqv(n)      = beta(n) * alphak_vdiff
       dbetaqv(1:n-1) = (beta(1:n-1) - beta(2:n))/deltaz(1:n-1)
       dbetaqv(0)     = zero
       dbetaqv(n)     = zero
    case (2)
       cql(0:n-1)     = half*(ciso(0:n-1)+ciso(1:n))
       dcqldca(0:n-1) = half
       dcqldcb(0:n-1) = half
       cql(n)         = ciso(n)
       dcqldca(n)     = one
       dcqldcb(n)     = zero

       cqv(0:n-1)     = half*(ciso(0:n-1)*beta(0:n-1)+ciso(1:n)*beta(1:n))
       dcqvdca(0:n-1) = half * beta(0:n-1)
       dcqvdcb(0:n-1) = half * beta(1:n)
       cqv(n)         = ciso(n)*beta(n)
       dcqvdca(n)     = one*beta(n)
       dcqvdcb(n)     = zero

       betaqv(0:n-1)  = alphak_vdiff
       betaqv(n)      = alphak_vdiff
       dbetaqv(1:n-1) = (beta(1:n-1) - beta(2:n))/deltaz(1:n-1)
       dbetaqv(0)     = zero
       dbetaqv(n)     = zero
    case (3)
       wcql(0) = one
       where (qlsig(1:n-1) > qmin)
          wcql(1:n-1) = var(1:n-1)%K*(1.-var(1:n-1)%h)/(var(1:n-1)%K*(1.-var(1:n-1)%h) +var(2:n)%K*var(2:n)%h)
          wcql(1:n-1) = dx(1:n-1)
       elsewhere
          wcql(1:n-1) = half
       endwhere
       cql(0:n-1)     = (wcql(0:n-1)*ciso(0:n-1)+(one-wcql(0:n-1))*ciso(1:n))
       dcqldca(0:n-1) = wcql(0:n-1)
       dcqldcb(0:n-1) = (one-wcql(0:n-1))
       cql(n)         = ciso(n)
       dcqldca(n)     = one
       dcqldcb(n)     = zero

       where (qvsig(1:n-1) > qmin)
          wcqv(1:n-1) = (var(1:n-1)%Kv*var(1:n-1)%h - var(1:n-1)%kE/rhow/rlambda*Tsoil(1:n-1))/ &
               (var(1:n-1)%Kv*var(1:n-1)%h - var(1:n-1)%kE/rhow/rlambda*Tsoil(1:n-1) &
               + var(2:n)%Kv*var(2:n)%h - var(2:n)%kE/rhow/rlambda*Tsoil(2:n))
       elsewhere
          wcqv(1:n-1) = half
       endwhere

       cqv(0)         = zero
       dcqvdca(0)     = zero
       dcqvdcb(0)     = zero
       cqv(1:n-1)     = (wcqv(1:n-1)*ciso(1:n-1)*beta(1:n-1)+(one-wcqv(1:n-1))*ciso(2:n)*beta(2:n))
       dcqvdca(1:n-1) = wcqv(1:n-1)* beta(1:n-1)
       dcqvdcb(1:n-1) = (one-wcqv(1:n-1)) * beta(2:n)
       cqv(n)         = zero
       dcqvdca(n)     = zero
       dcqvdcb(n)     = zero

       betaqv(0:n-1)  = alphak_vdiff
       betaqv(n)      = alphak_vdiff
       dbetaqv(1:n-1) = (beta(1:n-1) - beta(2:n))/deltaz(1:n-1)
       dbetaqv(0)     = zero
       dbetaqv(n)     = zero
    case default
       write(*,*) "isotope_vap: illegal formulation [1-3]: ", formulation
       stop
    end select

    ! mean diffusivities
    wl(0)         = h0
    wl(1:n)       =  dx(1:n)
    Dlmean(0:n-1) = (Dl(0:n-1)*wl(0:n-1) + Dl(1:n)*wl(1:n))/(wl(0:n-1)+wl(1:n))
    DLmean(n)     = zero
    where (Dv(1:n) > 1.e-16_r_2)
       wv(1:n) = dx(1:n)*Dv(1:n)
    elsewhere
       wv(1:n) = dx(1:n)
    endwhere
    wv(1:n)       = dx(1:n)
    wv(0)         = dxL
    Dvmean(0)     = zero
    Dvmean(1:n-1) = (Dv(1:n-1)*wv(1:n-1) + Dv(2:n)*wv(2:n))/(wv(1:n-1)+wv(2:n))
    Dvmean(n)     = zero
    if (litter .and. (ns==1)) then
       Dvmean(0) = (Dv(0)*wv(0) + Dv(1)*wv(1))/(wv(0)+wv(1))
    endif

    ! coefficients of tridiagonal matrix
    aa(1) = zero

    aa(2:n) = qlsig(1:n-1)*dcqldca(1:n-1) &
         + qvsig(1:n-1)*betaqv(1:n-1)*dcqvdca(1:n-1) &
         + dbetaqv(1:n-1)*dcqvdca(1:n-1)*Dvmean(1:n-1) &
         + Dlmean(1:n-1)/deltaz(1:n-1) &
         + Dvmean(1:n-1)/deltaz(1:n-1)*beta(1:n-1)

    bb(1)  = -(Seff(1)+kfreeze(1)*Sicesig(1)) * &
         thetasat(1)*dx(1)/sig/dt &                         !!!vh!!! NB thetasat should be (thetasat-thetar)??
         - qevapout*dcevapoutdciso &
         - qlsig(1)*dcqldca(1) - qvsig(1)*betaqv(1)*dcqvdca(1) &
         - dbetaqv(1)*dcqvdca(1)*Dvmean(1) &
         - Dlmean(1)/deltaz(1) &
         - Dvmean(1)/deltaz(1)*beta(1) &
         - (qex(1)+qrunoff)

    bb(2:n-1) = -(Seff(2:n-1)+kfreeze(2:n-1)*Sicesig(2:n-1)) * &
         thetasat(2:n-1)*dx(2:n-1)/sig/dt &
         + qlsig(1:n-2)*dcqldcb(1:n-2)  +qvsig(1:n-2)*betaqv(1:n-2)*dcqvdcb(1:n-2) &
         + dbetaqv(1:n-2)*dcqvdcb(1:n-2)* Dvmean(1:n-2) &
         - qlsig(2:n-1)*dcqldca(2:n-1) - qvsig(2:n-1)*betaqv(2:n-1)*dcqvdca(2:n-1) &
         - dbetaqv(2:n-1)*dcqvdca(2:n-1)*Dvmean(2:n-1) &
         - (Dlmean(1:n-2)/deltaz(1:n-2) + Dlmean(2:n-1)/deltaz(2:n-1) ) &
         - (Dvmean(1:n-2)/deltaz(1:n-2) + Dvmean(2:n-1)/deltaz(2:n-1))*beta(2:n-1) &
         - qex(2:n-1)

    bb(n) = -(Seff(n)+kfreeze(n)*Sicesig(n)) * &
         thetasat(n)*dx(n)/sig/dt &
         + qlsig(n-1)*dcqldcb(n-1)  +qvsig(n-1)*betaqv(n-1)*dcqvdcb(n-1) &
         + dbetaqv(n-1)*dcqvdcb(n-1) *Dvmean(n-1) &
         - qlsig(n)*dcqldca(n) &
         - Dlmean(n-1)/deltaz(n-1) &
         - Dvmean(n-1)/deltaz(n-1)*beta(n) &
         - qex(n)

    cc(1:n-1) = -qlsig(1:n-1)*dcqldcb(1:n-1)  - qvsig(1:n-1)*betaqv(1:n-1)*dcqvdcb(1:n-1) &
         - dbetaqv(1:n-1)*dcqvdcb(1:n-1)*Dvmean(1:n-1) &
         + Dlmean(1:n-1)/deltaz(1:n-1) &
         + Dvmean(1:n-1)/deltaz(1:n-1)*beta(2:n)

    cc(n) = zero

    dd(1) = thetasat(1)*dx(1)/sig/dt* &
         (ciso(1)*deltaSeff(1) + cisoice(1)*deltaSice(1) + &
         kfreeze2(1)*Sicesig(1)) &
         - qprec*cprec/sig - qevapin*cevapin/sig + qevapout*cevapout/sig &
         + qlsig(1)*cql(1)/sig + qvsig(1)*betaqv(1)*cqv(1)/sig &
         + dbetaqv(1)*cqv(1)*Dvmean(1)/sig &
         + Dlmean(1)/deltaz(1)/sig*(ciso(1) - ciso(2)) &
         + Dvmean(1)/deltaz(1)/sig*(ciso(1)*beta(1) - ciso(2)*beta(2) ) &
         + (qex(1)+qrunoff)*ciso(1)/sig

    dd(2:n-1) = thetasat(2:n-1)*dx(2:n-1)/sig/dt* &
         (ciso(2:n-1)*deltaSeff(2:n-1) + cisoice(2:n-1)*deltaSice(2:n-1) + &
         kfreeze2(2:n-1)*Sicesig(2:n-1)) &
         - qlsig(1:n-2)*cql(1:n-2)/sig &
         - qvsig(1:n-2)*betaqv(1:n-2)*cqv(1:n-2)/sig &
         - dbetaqv(1:n-2)*cqv(1:n-2)*Dvmean(1:n-2)/sig &
         + qlsig(2:n-1)*cql(2:n-1)/sig &
         + qvsig(2:n-1)*betaqv(2:n-1)*cqv(2:n-1)/sig &
         + dbetaqv(2:n-1)*cqv(2:n-1)*Dvmean(2:n-1)/sig &
         - Dlmean(1:n-2)/deltaz(1:n-2)/sig*(ciso(1:n-2) - ciso(2:n-1)) &
         - Dvmean(1:n-2)/deltaz(1:n-2)/sig*(ciso(1:n-2)*beta(1:n-2) - ciso(2:n-1)*beta(2:n-1)) &
         + Dlmean(2:n-1)/deltaz(2:n-1)/sig*(ciso(2:n-1) - ciso(3:n)) &
         + Dvmean(2:n-1)/deltaz(2:n-1)/sig*(ciso(2:n-1)*beta(2:n-1) - ciso(3:n)*beta(3:n)) &
         + qex(2:n-1)*ciso(2:n-1)/sig

    dd(n) = thetasat(n)*dx(n)/sig/dt* &
         (ciso(n)*deltaSeff(n) + cisoice(n)*deltaSice(n) + &
         kfreeze2(n)*Sicesig(n)) &
         - qlsig(n-1)*cql(n-1)/sig &
         - qvsig(n-1)*betaqv(n-1)*cqv(n-1)/sig &
         - dbetaqv(n-1)*cqv(n-1)*Dvmean(n-1)/sig &
         + qlsig(n)*cql(n)/sig &
         - Dlmean(n-1)/deltaz(n-1)/sig*(ciso(n-1) - ciso(n)) &
         - Dvmean(n-1)/deltaz(n-1)/sig*(ciso(n-1)*beta(n-1) - ciso(n)*beta(n) ) &
         + qex(n)*ciso(n)/sig

    if (cali>zero .or. experiment==7 .or. experiment==8) then
       bb(n) = bb(n) + qlsig(n)*dcqldca(n)
       dd(n) = dd(n) + qlsig(n)*(cali-cql(n))/sig
    endif

    if (litter) ns_ciso = 0

    if (ns_ciso==0) then
       aa(0) = zero

       bb(0) = -qevapout*dcevapoutdciso -qsig(0) - h0/sig/dt &
            - Dlmean(0)/(half*h0 + half*dx(1))

       cc(0) = Dlmean(0)/(half*h0 + half*dx(1))

       dd(0) = ciso(0)*dh0/dt/sig - cprec*qprec/sig + qevapout*cevapout/sig - qevapin*cevapin/sig + &
            qsig(0)*ciso(0)/sig + Dlmean(0)/(half*h0 + half*dx(1))/sig*(ciso(0) - ciso(1))

       aa(1) = qsig(0) + Dlmean(0)/(h0 + half*dx(1))

       bb(1)  = - Seff(1)*thetasat(1)*dx(1)/sig/dt &
            - qlsig(1)*dcqldca(1) - qvsig(1)*betaqv(1)*dcqvdca(1) &
            - (Dlmean(0)/(half*h0 + half*dx(1)) + Dlmean(1)/deltaz(1)) &
            - ( Dvmean(1)/deltaz(1))*beta(1) &
            - qex(1)

       cc(1) = -qlsig(1)*dcqldcb(1)  - qvsig(1)*betaqv(1)*dcqvdcb(1) &
            + Dlmean(1)/deltaz(1) &
            + Dvmean(1)/deltaz(1)*beta(2)

       dd(1) = thetasat(1)*dx(1)/sig/dt*ciso(1)*deltaSeff(1) &
            - qsig(0)*ciso(0)/sig &
            + qlsig(1)*cql(1)/sig + qvsig(1)*betaqv(1)*cqv(1)/sig &
            - Dlmean(0)/(half*h0 + half*dx(1))/sig*(ciso(0) - ciso(1)*S(1)*thetasat(1)) &
            + Dlmean(1)/deltaz(1)/sig*(ciso(1) - ciso(2)) &
            + Dvmean(1)/deltaz(1)/sig*(ciso(1)*beta(1) - ciso(2)*beta(2)) &  !! NB !! needs attention, esp S= Sliq or Sliqice??
            + qex(1)*ciso(1)/sig

       if (litter .and. (ns==1)) then               ! litter and no ponding
          bb(0) = -qevapout*dcevapoutdciso -qevapinL*dcevapindcisoL &
               - Seff(0)*thetasatL*dxL/sig/dt  - ( Dvmean(0)/deltaz0*beta(0) )

          cc(0) = qevapoutL*dcevapoutdcisoL + Dvmean(0)/deltaz0*beta(1)

          dd(0) = thetasatL*dxL/sig/dt*ciso(0)*deltaSeff(0) &
               - cprec*qprec/sig + qevapout*cevapout/sig - qevapin*cevapin/sig &
               - qevapoutL*cevapoutL/sig  + qevapinL*cevapinL/sig + qd*cprec/sig &
               + Dvmean(0)/deltaz0*(ciso(0)*beta(0) - ciso(1)*beta(1))/sig

          aa(1) = qevapinL*dcevapindcisoL + Dvmean(0)/deltaz0*beta(0)

          bb(1) = -Seff(1)*thetasat(1)*dx(1)/sig/dt &
               - qevapoutL*dcevapoutdcisoL &
               - qlsig(1)*dcqldca(1) - qvsig(1)*betaqv(1)*dcqvdca(1) &
               - Dvmean(0)/deltaz0*beta(1) &
               - Dlmean(1)/deltaz(1) &
               - (Dvmean(1)/deltaz(1))*beta(1) &
               - (qex(1)+qrunoff)

          cc(1) = -qlsig(1)*dcqldcb(1)  -qvsig(1)*betaqv(1)*dcqvdcb(1) &
               + Dlmean(1)/deltaz(1) &
               + Dvmean(1)/deltaz(1)*beta(2)

          dd(1) = thetasat(1)*dx(1)/sig/dt*ciso(1)*deltaSeff(1) &
               - qd*cprec/sig + qevapoutL*cevapoutL/sig -qevapinL*cevapinL/sig &
               + qlsig(1)*cql(1)/sig + qvsig(1)*betaqv(1)*cqv(1)/sig &
               - Dvmean(0)/(deltaz0)/sig*(ciso(0)*beta(0) - ciso(1)*beta(1)) &
               + Dlmean(1)/deltaz(1)/sig*(ciso(1) - ciso(2)) &
               + Dvmean(1)/deltaz(1)/sig*(ciso(1)*beta(1) - ciso(2)*beta(2)) &
               + (qex(1)+qrunoff)*ciso(1)/sig
       endif                      ! end litter

    endif ! ns_ciso = 0

    dc(0) = zero

    call tri(ns_ciso,n,aa,bb,cc,dd,dc)

    dcice(1:n) = sig*dc(1:n)*kfreeze(1:n) + kfreeze2(1:n)

    ! check for mass balance
    if (1 == 1) then
       LHS (1:n) = thetasat(1:n)*dx(1:n)/dt*(dc(1:n)*Seff(1:n) + ciso(1:n)*deltaSeff(1:n)) + &
            thetasat(1:n)*dx(1:n)/dt*(dcice(1:n)*Sice(1:n) + cisoice(1:n)*deltaSice(1:n))

       RHS(2:n-1) = qlsig(1:n-2)*(cql(1:n-2) + sig*dc(1:n-2)*dcqldca(1:n-2) + sig*dc(2:n-1)*dcqldcb(1:n-2)) &
            + qvsig(1:n-2)*betaqv(1:n-2)*(cqv(1:n-2) + sig*dc(1:n-2)*dcqvdca(1:n-2) &
            + sig*dc(2:n-1)*dcqvdcb(1:n-2)) &
            -qlsig(2:n-1)*(cql(2:n-1) + sig*dc(2:n-1)*dcqldca(2:n-1) + sig*dc(3:n)*dcqldcb(2:n-1) ) &
            - qvsig(2:n-1)*betaqv(2:n-1)*(cqv(2:n-1) + sig*dc(2:n-1)*dcqvdca(2:n-1) + sig*dc(3:n)*dcqvdcb(2:n-1)) &
            + Dlmean(1:n-2)*((ciso(1:n-2) + sig*dc(1:n-2))-(ciso(2:n-1)+sig*dc(2:n-1)))/deltaz(1:n-2) &
            - Dlmean(2:n-1)*((ciso(2:n-1) + sig*dc(2:n-1))-(ciso(3:n)+sig*dc(3:n)))/deltaz(2:n-1) &
            + Dvmean(1:n-2)*((ciso(1:n-2) + sig*dc(1:n-2))*beta(1:n-2) - (ciso(2:n-1) &
            + sig*dc(2:n-1))*beta(2:n-1))/deltaz(1:n-2) &
            - Dvmean(2:n-1)*((ciso(2:n-1) + sig*dc(2:n-1))*beta(2:n-1) - (ciso(3:n) &
            + sig*dc(3:n))*beta(3:n))/deltaz(2:n-1) &
            - qex(2:n-1)*(ciso(2:n-1) + sig*dc(2:n-1))

       RHS(1) = qprec*cprec - qevapout*(cevapout + sig*dc(1)*dcevapoutdciso) +qevapin*cevapin &
            -qlsig(1)*(cql(1) + sig*dc(1)*dcqldca(1) + sig*dc(2)*dcqldcb(1)) &
            - qvsig(1)*betaqv(1)*(cqv(1) + sig*dc(1)*dcqvdca(1) + sig*dc(2)*dcqvdcb(1)) &
            - Dlmean(1)*((ciso(1) + sig*dc(1))-(ciso(2)+sig*dc(2)))/deltaz(1) &
            - Dvmean(1)*((ciso(1) + sig*dc(1))*beta(1) - (ciso(2) + sig*dc(2))*beta(2))/deltaz(1) &
            - (qex(1)+qrunoff)*(ciso(1) + sig*dc(1))

       RHS(n) =  qlsig(n-1)*(cql(n-1) + sig*dc(n-1)*dcqldca(n-1) + sig*dc(n)*dcqldcb(n-1)) &
            + qvsig(n-1)*betaqv(n-1)*(cqv(n-1) + sig*dc(n-1)*dcqvdca(n-1) + sig*dc(n)*dcqvdcb(n-1)) &
            -qsig(n)*(ciso(n) + sig*dc(n)) &
            + Dlmean(n-1)*((ciso(n-1) + sig*dc(n-1))-(ciso(n)+sig*dc(n)))/deltaz(n-1) &
            + Dvmean(n-1)*((ciso(n-1) + sig*dc(n-1))*beta(n-1) - (ciso(n) + sig*dc(n))*beta(n))/deltaz(n-1) &
            -qex(n)*(ciso(n)+sig*dc(n))

       if (cali>zero .or. experiment==7 .or. experiment==8) then
          RHS(n) = RHS(n) - qsig(n) * (cali - (ciso(n)+sig*dc(n)))
       endif

       !MC if (any(abs(LHS(1:n)-RHS(1:n))>epsilon(one))) then
       if (any(abs(LHS(1:n)-RHS(1:n))>1.e-11_r_2)) then
          write(*,*) 'Max of abs(LHS-RHS): ', maxval(abs(LHS(1:n)-RHS(1:n)))
       endif

       if (ns_ciso==0) then
          LHS(0) = (ciso(0)*dh0 + h0*dc(0))/dt

          RHS(0) = qprec*cprec - qevapout*(cevapout + sig*dc(0)*dcevapoutdciso) &
               + qevapin*cevapin - qsig(0)*(ciso(0) + sig*dc(0)) &
               - Dlmean(0)*((ciso(0) + sig*dc(0))-(ciso(1)+sig*dc(1)))/(half*h0 + half*dx(1))

          RHS(1) = qsig(0)*(ciso(0)+sig*dc(0)) &
               -qlsig(1)*(cql(1) + sig*dc(1)*dcqldca(1) + sig*dc(2)*dcqldcb(1)) &
               - qvsig(1)*betaqv(1)*(cqv(1) + sig*dc(1)*dcqvdca(1) + sig*dc(2)*dcqvdcb(1)) &
               + Dlmean(0)*((ciso(0) + sig*dc(0))-(ciso(1)+sig*dc(1)))/(half*h0 + half*dx(1)) &
               - Dlmean(1)*((ciso(1) + sig*dc(1))-(ciso(2)+sig*dc(2)))/deltaz(1) &
               - Dvmean(1)*((ciso(1) + sig*dc(1))*beta(1) - (ciso(2) + sig*dc(2))*beta(2))/deltaz(1) &
               - qex(1)*(ciso(1) + sig*dc(1))

          if (litter .and. (ns==1)) then               ! litter and no ponding
             RHS(0) = thetasatL*dxL/dt*(ciso(0)*deltaSeff(0)+ dc(0)*Seff(0) )

             LHS(0) = qprec*cprec-qd*cprec - qevapout*(cevapout+sig*dcevapoutdciso*dc(0)) &
                  + qevapin*cevapin &
                  + qevapoutL*(cevapoutL+sig*dcevapoutdcisoL*dc(1)) &
                  -qevapinL*(cevapinL + sig*dcevapindcisoL *dc(0)) &
                  -Dvmean(0)/(deltaz0)*(beta(0)*(ciso(0)+dc(0)) - beta(1)*(ciso(1)+sig*dc(1)) )

             RHS(1) = qd*cprec - qevapoutL*(cevapoutL+ sig*dcevapoutdcisoL*dc(1)) &
                  +qevapinL*(cevapinL + sig*dcevapindcisoL*dc(0)) &
                  -qlsig(1)*(cql(1) + sig*dc(1)*dcqldca(1) + sig*dc(2)*dcqldcb(1)) &
                  - qvsig(1)*betaqv(1)*(cqv(1) + sig*dc(1)*dcqvdca(1) + sig*dc(2)*dcqvdcb(1)) &
                  + Dvmean(0)*(beta(0)*(ciso(0) + sig*dc(0))-beta(1)*(ciso(1)+sig*dc(1)))/(deltaz0) &
                  - Dlmean(1)*((ciso(1) + sig*dc(1))-(ciso(2)+sig*dc(2)))/deltaz(1) &
                  - Dvmean(1)*((ciso(1) + sig*dc(1))*beta(1) - (ciso(2) + sig*dc(2))*beta(2))/deltaz(1) &
                  - (qex(1)+qrunoff)*(ciso(1) + sig*dc(1))
          endif                      ! end litter and no ponding
       end if ! ns_ciso==0

    endif ! 1==1

    ! isotopic fluxes
    qiso_in    = qprec*cprec +qevapin*cevapin
    qiso_out   = qevapout*(cevapout + sig*dc(1)*dcevapoutdciso)
    qiso_evap  = qevapout*(cevapout + sig*dc(1)*dcevapoutdciso) - qevapin*cevapin
    qiso_trans = sum(qex(1:n)*(ciso(1:n)+sig*dc(1:n)),1)

    qiso_liq_diff(1:n-1) = Dlmean(1:n-1)*((ciso(1:n-1) + sig*dc(1:n-1))-(ciso(2:n)+sig*dc(2:n)))/deltaz(1:n-1)
    qiso_vap_diff(1:n-1) = Dvmean(1:n-1)*((ciso(1:n-1) + sig*dc(1:n-1))*beta(1:n-1) &
         - (ciso(2:n) + sig*dc(2:n))*beta(2:n))/deltaz(1 :n-1)
    qiso_liq_adv(1:n-1)  = qlsig(1:n-1)*(cql(1:n-1) + sig*dc(1:n-1)*dcqldca(1:n-1) + sig*dc(2:n)*dcqldcb(1:n-1))
    qiso_liq_adv(n)      = qlsig(n)*(cql(n) + sig*dc(n))
    qiso_vap_adv(1:n-1)  = qvsig(1:n-1) * betaqv(1:n-1) &
         * (cqv(1:n-1) + sig*dc(1:n-1)*dcqvdca(1:n-1) + sig*dc(2:n)*dcqvdcb(1:n-1))
    qiso_vap_adv(n)      = zero

    ciso = ciso + dc
    cisoice(1:n) = cisoice(1:n) + dcice(1:n)

    if (litter .and. (ns==1)) cisoL = ciso(0)

  END SUBROUTINE isotope_vap

  !**********************************************************************************************************************

END MODULE sli_solve
