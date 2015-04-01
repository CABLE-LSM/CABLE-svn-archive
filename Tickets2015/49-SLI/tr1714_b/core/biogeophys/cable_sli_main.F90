SUBROUTINE sli_main(ktau, dt, veg, soil, ssoil, met, canopy, air, rad, SEB_only)

  ! Main subroutine for Soil-litter-iso soil model
  ! Vanessa Haverd, CSIRO Marine and Atmospheric Research
  ! and Matthias Cuntz, UFZ - Helmholtz Centre for Environmental Research, 2010
  ! Modified to operate for multiple veg tiles but a single soil column, March 2011
  ! Rewritten for same number of soil columns as veg tiles May 2012
  USE cable_def_types_mod,       ONLY: veg_parameter_type, soil_parameter_type, soil_snow_type, met_type, &
       canopy_type, air_type, radiation_type, ms, mp, r_2, i_d
  !USE physical_constants, ONLY: tfrz
  USE sli_numbers,        ONLY:  zero, half, one, four, thousand, & ! numbers
       Tzero, &                                       ! variables
       vars_met, vars, params, vars_snow, &                                  ! types
       MW, snmin, Rgas, Lambdas, lambdaf, csice, cswat, rhow, nsnow_max, e5
  USE sli_utils,          ONLY: x, dx, par, setpar, setx, plit, dxL, setlitterpar, esat, &
       esat_ice, slope_esat_ice, thetalmax, Tfrz,  hyofS, SEB
  USE sli_roots,          ONLY: setroots, getrex
  USE sli_solve,          ONLY: solve
  USE cable_soil_snow_module, ONLY: snowdensity



  IMPLICIT NONE


  REAL,                      INTENT(IN)    :: dt
  TYPE(veg_parameter_type),  INTENT(INOUT) :: veg     ! all r_1
  TYPE(soil_parameter_type), INTENT(INOUT) :: soil    ! all r_1
  TYPE(soil_snow_type),      INTENT(INOUT) :: ssoil   ! r_1, r_2 desaster
  TYPE(met_type),            INTENT(INOUT) :: met     ! all r_1
  TYPE(canopy_type),         INTENT(INOUT) :: canopy  ! all r_1
  TYPE(air_type),            INTENT(INOUT) :: air     ! all r_1
  TYPE (radiation_type),       INTENT(IN) :: rad
  INTEGER,      INTENT(IN) :: ktau ! integration step number
  INTEGER,      INTENT(IN) :: SEB_only ! integration step number


  REAL(r_2),    PARAMETER :: emsoil=0.97
  REAL(r_2),    PARAMETER :: rhocp=1.1822e3
  REAL(r_2),    PARAMETER :: Dva = 2.17e-5
  INTEGER(i_d) :: i, j, k, kk, setroot
  REAL(r_2)    :: ti, tf
  TYPE(vars_met), DIMENSION(1:mp)      :: vmet ! Meteorology above soil
  TYPE(vars),     DIMENSION(1:mp)      :: vlit
  TYPE(vars),     DIMENSION(1:mp,1:ms)      :: var
  TYPE(vars_snow),     DIMENSION(1:mp) :: vsnow
  INTEGER(i_d),   DIMENSION(1:mp)      :: nsteps
  REAL(r_2),      DIMENSION(1:mp,1:ms) :: Tsoil, S, thetai, Jsensible
  REAL(r_2),      DIMENSION(1:mp)      :: SL, TL, T0
  REAL(r_2),      DIMENSION(1:mp)      :: drn, evap, infil, qprec,qprec_snow, qprec_tot, runoff, runoff_sat
  REAL(r_2),      DIMENSION(1:mp)      :: win, wp, wpi, h0,deltah0, h0old,hsnowold, discharge
  REAL(r_2),      DIMENSION(1:mp)      :: ip, ipi ! volumetric ice content of profile (final and initial)
  REAL(r_2),      DIMENSION(1:mp,1:ms) ::   wex, csoil, qex, kth, phi, thetal_max, Sliq, Ksat
  REAL(r_2),      DIMENSION(1:mp,1:ms) ::  FS
  REAL(r_2),      DIMENSION(1:mp)      :: rh0, rhsurface
  ! surface temperature (top of top soil layer or top of litter layer)
  REAL(r_2),      DIMENSION(1:mp)      :: Tsurface
  REAL(r_2),      DIMENSION(1:mp)      :: gr, grc
  REAL(r_2),      DIMENSION(1:mp)      :: Etrans
  REAL(r_2),      DIMENSION(1:mp)      :: gamma
  REAL(r_2),      DIMENSION(1:mp)      :: G0, H, lE
  REAL(r_2),      DIMENSION(1:mp,-nsnow_max:ms) :: qh, qvsig, qlsig, qvTsig, qvh
  REAL(r_2),      DIMENSION(1:mp)      :: deltaTa, lE_old !, SA, SB, wpAi, wpBi, wpA, wpB
  REAL(r_2),      DIMENSION(1:mp)      :: evap_pot, deltaice_cum_T, deltaice_cum_S, zdelta
  REAL(r_2),      DIMENSION(1:mp)      ::fws
  REAL(r_2),      DIMENSION(1:mp)      :: Qadvcum,Jcol_sensible,Jcol_latent_S,Jcol_latent_T
  REAL(r_2),      DIMENSION(1:mp)      :: tmp1d1, deltaEsnow
  REAL(r_2),      DIMENSION(1:mp)      :: hice, phie
  REAL(r_2)                           ::  tmp1d1a, tmp1d2, tmp1d3, tmp1d4, &
                                           tmp1d5, tmp1d6, tmp1d7, tmp1d8, tmp1d9,tmp1d10, tmp1d11, &
                                           tmp1d12,tmp1d13, tmp1d14, tmp1d15, tmp1d16

  REAL(r_2), PARAMETER :: alpha = 0.1 ! anistropy param for lateral flow
  REAL(r_2), PARAMETER :: fsat_max = 2.0 ! exponent for vetical profile of Ksat
  !REAL(r_2),      DIMENSION(1:mp)      :: qb ! topmodel baseflow
  REAL(r_2),      DIMENSION(1:mp)      :: fsat ! topmodel saturated area
  INTEGER(i_d), DIMENSION(1:mp)      :: index
  ! Model switches
  INTEGER(i_d), PARAMETER :: litter       = 0 ! which litter model
  ! 0: no litter
  ! 1: full litter
  ! 2: litter resistance
  INTEGER(i_d), PARAMETER :: advection    = 1 ! heat advection by water
  INTEGER(i_d), PARAMETER :: isotopologue = 0 ! which isotope
  ! 0: no isotope calculations
  ! 1: HDO
  ! 2: H218O
  ! 3: HDO & H218O
  ! 0: normal run
  INTEGER(i_d), PARAMETER :: septs        = 0 ! coupled or uncoupled energy and water calculation
  ! 0: coupled calc
  ! 1: uncoupled energy (T) and moisture (S)
  INTEGER(i_d), PARAMETER :: condition    = 3 ! condition matrix before solving
  ! 0: no conditioning
  ! 1: condition columns
  ! 2: condition lines
  ! 3: condition first lines then columns
  INTEGER, PARAMETER :: dosnow       = 1 ! implement snow model
  LOGICAL, SAVE :: first = .true.
  INTEGER(i_d), SAVE  :: counter

  ! initialise cumulative variables
  Jcol_sensible(:) = zero
  Jcol_latent_S(:) = zero
  Jcol_latent_T(:) = zero
  deltaice_cum_T(:) = zero
  deltaice_cum_S(:) = zero
  Jsensible(:,:) = zero
  drn(:)  = zero
  discharge(:) = zero
  infil(:)     = zero
  evap(:)      = zero
  evap_pot(:)  = zero
  runoff(:)    = zero
  qh = zero
  H(:)      = zero
  G0(:)      = zero
  lE(:)     = zero
  csoil = zero
  kth = zero
  Qadvcum(:)  = zero
  wex(:,:)          = zero
  qlsig(:,:)        = zero
  qvsig(:,:)        = zero
  qvh(:,:)        = zero
  qvtsig(:,:)        = zero
  thetal_max   = zero
  Sliq         = zero
  Ksat         = zero
  phie         = zero
  phi          = zero
  hice         = zero

  ! output files for testing purposes
  if (first) then
     open (unit=332,file="vh08.out",status="replace",position="rewind")
     open (unit=334,file="S.out",status="replace",position="rewind")
     open (unit=336,file="Tsoil.out",status="replace",position="rewind")
     open (unit=335,file="SEB.out",status="replace",position="rewind")
     !   !open(unit=37, file = "c:\soil_model\cable_met_test.inp",status="replace",position="rewind")
     open (unit=337,file="soil_log.out",status="replace",position="rewind")
     open(unit=338, file="thetai.out", status="replace", position="rewind")
     open(unit=339, file="latlong.out",status="replace", position="rewind")
     write(339,"(20000f8.2)") rad%latitude
     write(339,"(20000f8.2)") rad%longitude
     open(unit=340, file="snow.out", status="replace", position="rewind")
     open(unit=345, file="diags.out",status="replace", position="rewind")
	 counter = 0
  endif

  fsat = 0.0
  zdelta = 0.0
  counter = counter + 1

  ! Save soil / snow surface temperature from last time step:
  ssoil%otss(:) = ssoil%tss(:)

  ! set layer thicknesses
  if (.not. allocated(x)) call setx(mp, ms, soil)

  ! Set root density distribution (leave in for sli offline)
  setroot = 1  ! reset rooting depths
  if (setroot == 1) then
     call setroots(x(:,:)*100.0_r_2, real(veg%F10(:),r_2), real(veg%ZR(:),r_2)*100.0_r_2, FS(:,:))
  else
     FS(:,:) = real(veg%froot(:,:),r_2)
  endif

  ! set required soil hydraulic params
  if (.not. allocated(par)) then
     index=(/(i,i=1,mp,1)/)
     call setpar(mp, ms, x(:,:)-half*dx(:,:), soil,index)
  endif

  ! If we want solutes:
!!$     if (.not. allocated(bd)) allocate(bd(soil%nhorizons(k)))

  ! Litter parameters:
  if (.not. allocated(plit)) then
     index=(/(i,i=1,mp,1)/)
     call setlitterpar(mp, soil,index)
  endif

  ! Met data above soil:
 
  vmet(:)%Ta    = real(met%Tvair-273.16,r_2)
  vmet(:)%Da    = real(met%dva,r_2)
  write(55,*) ssoil%rtsoil
  vmet(:)%rbh   = ssoil%rtsoil
  vmet(:)%rbw   = vmet(:)%rbh
  vmet(:)%rha   = max(min((esat(vmet(:)%Ta)-vmet(:)%Da)/esat(vmet(:)%Ta),one),0.1_r_2)
  vmet(:)%cva   = vmet(:)%rha * esat(vmet(:)%Ta)*0.018_r_2/thousand/8.314_r_2/(vmet(:)%Ta+Tzero) ! m3 H2O (liq) m-3 (air)
  vmet(:)%phiva = Dva * vmet(:)%cva
  vmet(:)%Rn    = canopy%fns
  Etrans(:)     = max(canopy%fevc/air%rlam/thousand, zero) ! m s-1
  h0(:) = ssoil%h0
 

  ! Initialisations:
  if (first) then
     rh0(:)         = vmet(:)%rha
     rhsurface(:)   = rh0(:) ! initialise rel. humidity at surface
     Tsurface(:)    = vmet(:)%Ta
     T0(:)          = Tsurface(:)
     deltaTa(:)     = zero
     lE_old(:)      = zero
     ssoil%h0(:)    = zero
     ssoil%thetai(:,:) = zero
     ! do not seem to be initialised elsewhere
     ssoil%snowd(:) = 0.0
     ssoil%smass(:,:) = 0.0
     zdelta(:)      = x(:,ms)
     ssoil%gammzz = zero
     !ssoil%tgg(:,6)=5.0+Tfrz

     ssoil%ssdn(:,:) = 120_r_2 ! snow density kg m-3
     ssoil%sconds(:,:) = 0.06_r_2    ! snow thermal cond (W m-2 K-1)
     ssoil%S  = min(one,ssoil%S)
     ssoil%Ta_daily(:,:) = spread(vmet(:)%Ta,2,100)

     !h0 = zero
     do kk=1, mp
        vsnow(kk)%hsnow(:) = zero
        vsnow(kk)%depth(:) = zero
		vsnow(kk)%totdepth = zero
        vsnow(kk)%wcol     = zero
        vsnow(kk)%dens(:)  = 120._r_2 ! snow density kg m-3
        vsnow(kk)%tsn(:)   = zero
        vsnow(kk)%kH(:)    = 0.16_r_2 ! snow thermal cond (W m-2 K-1)
        vsnow(kk)%Dv(:)    = Dva      ! m2 s-1
        vsnow(kk)%sl(:)    = zero
        vsnow(kk)%kE(:)    = zero
        vsnow(kk)%kth(:)   = vsnow(kk)%kH
        vsnow(kk)%cv(:)    = zero
        vsnow(kk)%hliq(:)  = zero
        vsnow(kk)%melt(:)  = zero
        vsnow(kk)%nsnow    = 0
        vsnow(kk)%nsnow_last    = 0
        ssoil%cls(kk) = one
     enddo
     first = .false.
  else
     rhsurface(:)   = ssoil%rhsurface
     Tsurface(:)    = ssoil%Tsurface
     T0(:)          = ssoil%Tsurface
     rh0(:)         = ssoil%rh0
     deltaTa(:)     = zero
     lE_old(:)      = ssoil%lE
     zdelta(:)      = ssoil%zdelta
  endif

  SL(:)      = 0.5   ! degree of litter saturation
  Tsoil(:,:) = real(ssoil%tgg-273.16,r_2)
  TL(:)      = Tsoil(:,1) ! litter T

  thetai = ssoil%thetai
  ssoil%smelt = zero
  ssoil%cls = 1.0
  S(:,:)         = ssoil%S                ! degree of soil saturation


  ! ----------------------------------------------------------------
  ! Iinitialise phi where it is (frozen and saturated) and where (pond >zero)

  where (Tsoil<Tfrz(S,par%he,one/par%lam))
     thetal_max    = thetalmax(Tsoil,S,par%he,one/par%lam,par%thre,par%the)
     Sliq          = (thetal_max - (par%the-par%thre))/par%thre
     Ksat          = par%Ke*exp(par%eta*log(Sliq))
  elsewhere
     Sliq          = S
     Ksat          = par%Ke
  endwhere

  where ((Tsoil<Tfrz(S,par%he,one/par%lam)) .and. (S>=1))
     phi  = par%phie*exp(-log(Sliq)/par%lam)*exp(par%eta*log(Sliq))
  endwhere

  where (h0(:)>zero)
     hice(:)  = h0(:)*(S(:,1)-Sliq(:,1))
     phie(:)  = par(:,1)%phie*exp(-log(Sliq(:,1))/par(:,1)%lam)*exp(par(:,1)%eta*log(Sliq(:,1)))
     !phi(:,1) = max((phie -par(:,1)%he*Ksat(:,1)), (one+e5)*phie)+(h0(:)-hice(:))*Ksat(:,1)
     phi(:,1) = (one+e5)*phie + (h0(:)-hice(:))*Ksat(:,1)
  endwhere
  var(:,:)%phi = phi


  ! ----------------------------------------------------------------
  do kk=1, mp
     vsnow(kk)%nsnow = ssoil%nsnow(kk)
     if (ssoil%sdepth(kk,1).gt.zero)  then ! snow
        ! define variables associated with snow
        ssoil%isflag(kk) = 1
        vsnow(kk)%hsnow(:) = ssoil%smass(kk,1:nsnow_max)/thousand
        vsnow(kk)%depth(:) = ssoil%sdepth(kk,1:nsnow_max)  ! depth of snow pack (m)
        vsnow(kk)%dens(:) =ssoil%ssdn(kk,1:nsnow_max)
		where (vsnow(kk)%dens(:).le.200.)
		   vsnow(kk)%fsnowliq_max(:) = 0.03
		elsewhere
           vsnow(kk)%fsnowliq_max(:) = 0.03 + (0.1 - 0.03)*(vsnow(kk)%dens(:)-200.)/vsnow(kk)%dens(:)
		endwhere
        vsnow(kk)%fsnowliq_max(:) = 0.1
        vsnow(kk)%wcol= ssoil%snowd(kk)/thousand

        vsnow(kk)%tsn(:) = ssoil%tggsn(kk,1:nsnow_max) - Tzero
        vsnow(kk)%kH(:) =  ssoil%sconds(kk,1:nsnow_max)
        vsnow(kk)%Dv(:) = Dva*(ssoil%tggsn(kk,1:nsnow_max)/Tzero)**1.88_r_2 ! m2 s-1
        vsnow(kk)%sl(:) = slope_esat_ice(vsnow(kk)%tsn(:)) * Mw/thousand/Rgas/(vsnow(kk)%tsn(:)+Tzero)
        vsnow(kk)%kE(:)     = vsnow(kk)%Dv(:)*vsnow(kk)%sl(:)*thousand*lambdaf
        vsnow(kk)%kth(:) = vsnow(kk)%kE(:) + vsnow(kk)%kH(:)
        vsnow(kk)%cv(:) = esat_ice(vsnow(kk)%tsn(:))*Mw/thousand/Rgas/(vsnow(kk)%tsn(:)+Tzero) ! m3 m-3
        vsnow(kk)%hliq(:) = ssoil%snowliq(kk,1:nsnow_max)/thousand ! amount of liq snow water
        vsnow(kk)%melt(:) = zero ! amount of melted snow leaving each snowlayer (mm/dt)


     else
        ssoil%isflag(kk) = 0
        vsnow(kk)%hsnow(:) = zero
        vsnow(kk)%depth(:) = zero
        vsnow(kk)%wcol = zero
        vsnow(kk)%dens(:) = 120_r_2 ! snow density kg m-3
        vsnow(kk)%tsn(:) = zero
        vsnow(kk)%kH(:) = 0.16_r_2    ! snow thermal cond (W m-2 K-1)
        vsnow(kk)%Dv(:) = Dva*(Tzero/Tzero)**1.88_r_2 ! m2 s-1
        vsnow(kk)%sl(:) = zero
        vsnow(kk)%kE(:)     = zero
        vsnow(kk)%kth(:) = vsnow(kk)%kH
        !vsnow(kk)%kth(:) = 0.07 ! test ! vh
        vsnow(kk)%cv(:) = zero
        vsnow(kk)%hliq(:) = zero
        vsnow(kk)%melt(:) = zero
        vsnow(kk)%Jlatent(:) = zero
        vsnow(kk)%Jsensible(:) = zero
		vsnow(kk)%J = zero
		vsnow(kk)%nsnow = 0
		vsnow(kk)%fsnowliq_max(:) = 0.03
     endif
     ssoil%osnowd(kk) = ssoil%snowd(kk)
  enddo

  deltaTa(:)  = zero
  lE_old(:)   = ssoil%lE
  gamma(:)    = real(veg%gamma(:),r_2)

  gr(:)       = four * emsoil * (vmet(:)%Ta+273.16)**3 *5.67e-8_r_2 ! radiation conductance Wm-2K-1
  grc(:)      = one/vmet(:)%rbh + gr(:)/rhocp
  vmet(:)%rrc = one/grc(:)                            ! resistance to radiative and convective heat transfer
  qprec(:)    = (canopy%through-met%precip_sn)/thousand/dt              ! liq precip rate (m s-1)
  qprec_snow(:) = (met%precip_sn)/thousand/dt

  ! re-calculate qprec_snow and qprec based on total precip and air T (ref Jin et al. Table II, Hyd Proc, 1999
  qprec_tot = qprec + qprec_snow
!  where (vmet(:)%Ta.gt.2.5)
!     qprec_snow = zero
!     qprec = qprec_tot
!  elsewhere (vmet(:)%Ta.le.2.5.and.vmet(:)%Ta.gt.2.0)
!     qprec_snow = 0.6*qprec_tot
!     qprec = qprec_tot - qprec_snow
!  elsewhere (vmet(:)%Ta.le.2.0.and.vmet(:)%Ta.gt.0)
!     qprec_snow = (1. - (54.62 - 0.2 *(vmet(:)%Ta + 273.16)))*qprec_tot
!     qprec = qprec_tot - qprec_snow
!  elsewhere (vmet(:)%Ta.le.0)
!     qprec = zero
!     qprec_snow = qprec_tot
!  endwhere


 ! write(*,*) ktau, vmet(:)%Ta, qprec_tot, qprec, qprec_snow





  h0old(:)    = ssoil%h0 ! pond height

      do kk=1, mp
         hsnowold(kk) = sum(vsnow(kk)%hsnow(1:nsnow_max))

      enddo


  ! Heat balance variables

  ! Water balance variables:
  ipi(:)   = sum(ssoil%thetai*dx(:,:),2)  + h0(:)*ssoil%thetai(:,1)/par(:,1)%thre       ! ice in profile initially
  ! water in profile initially
  wpi(:)   = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*dx(:,:),2)  + plit(:)%thre*SL(:)*dxL(:)


  ! saturated fraction
  !zdelta(:) = x(:,ms)-wpi(:)
  ! zdelta = max((sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr))*dx(:,:),2) - wpi(:)),zero)
  ! fsat(:) = min(par(:,ms)%fsatmax * exp(-zdelta),1.0)
  !runoff_sat = fsat(:)*qprec(:)*thousand*dt
  ! qprec = qprec*(1-fsat)  ! precip available for infiltration

  runoff_sat = zero

  nsteps    = 0
  ti        = zero
  tf        = dt ! initial and final times

  win(:)    = zero ! water input (total precip)
  evap(:)   = zero
  runoff(:) = zero
  infil(:)  = zero
  drn(:)    = zero

  do kk=1, mp
     call getrex(ssoil%S(kk,:), ssoil%rex(kk,:), fws(kk), FS(kk,:), soil%ssat_vec(kk,:), &
          soil%swilt_vec(kk,:), real(Etrans(kk),r_2), gamma(kk), dx(kk,:), real(dt,r_2))
  enddo

  do k=1,ms
     qex(:,k)= ssoil%rex(:,k)
  enddo

  ! calculate base flow (topmodel)
  !qb = alpha * par(:,ms)%Ke*exp(-(par(:,ms)%zeta*zdelta))
  qex(:,ms) =  qex(:,ms) !+qb
  !  write(*,*) vmet, SEB_only
 ! stop
 
  if (SEB_only.eq.1) then

     do kk=1, mp
     !call hyofS(S(kk,:), Tsoil(kk,:), par(kk,:), var(kk,:)) 
      do i=1, 1
                    call hyofS(S(kk,i), Tsoil(kk,i), par(kk,i), var(kk,i))
      end do
     CALL SEB(ms, par(kk,:), vmet(kk), vsnow(kk), var(kk,:), qprec(kk), qprec_snow(kk),  1, dx(kk,:), &
                      h0(kk), hice(kk), S(kk,:), Tsoil(kk,:), &
                      Tsurface(kk), G0(kk), lE(kk),  &
                      tmp1d1a, tmp1d2, tmp1d3, tmp1d4, &
                      tmp1d5, tmp1d6, tmp1d7, tmp1d8, tmp1d9,tmp1d10, tmp1d11, &
                      tmp1d12,tmp1d13, tmp1d14, tmp1d15, tmp1d16,ktau)
     enddo
       canopy%ga       = real(G0(:))
       canopy%fes      = real(lE(:))
       canopy%fhs = canopy%fns - canopy%ga - canopy%fes
       ssoil%tss       = real(Tsurface(:)) + 273.16
       !write(*,"(a10, 100f16.6)") "sli_main SEB", canopy%fhs, canopy%fes, ssoil%rtsoil, vmet(1)%Ta
  else 
    
     call solve(ti, tf, ktau, mp, qprec(:),qprec_snow(:), ms,  dx(:,:), &
           h0(:), S(:,:), thetai(:,:), Jsensible(:,:), Tsoil(:,:), evap(:), &
           evap_pot(:), runoff(:), infil(:), drn(:), discharge(:), qh(:,:), &
           nsteps, vmet(:), vlit(:),vsnow(:), var(:,:),csoil(:,:), kth(:,:), phi(:,:), T0(:), rh0(:), Tsurface(:), &
           rhsurface(:), H(:), lE(:), G0(:),Qadvcum(:),Jcol_sensible(:), &
           Jcol_latent_S(:),Jcol_latent_T(:), deltaice_cum_T(:), &
           deltaice_cum_S(:), dxL(:), zdelta(:), SL(:), TL(:), &
           plit(:), par(:,:), qex=qex(:,:), &
           wex=wex(:,:), FS=FS(:,:), qvsig=qvsig(:,:), qlsig=qlsig(:,:), qvTsig=qvTsig(:,:), qvh=qvh(:,:), &
           deltaTa=deltaTa(:), lE_old=lE_old(:), &
           dolitter=litter, doisotopologue=isotopologue, dosepts=septs, docondition=condition, &
           doadvection=advection)

  H(:)      = H(:)/(tf-ti)
  lE(:)     = lE(:)/(tf-ti)
  G0(:)     = G0(:)/(tf-ti)
  Jcol_latent_S(:) = Jcol_latent_S(:)/(tf-ti)
  Jcol_latent_T(:) = Jcol_latent_T(:)/(tf-ti)
  Jcol_sensible(:) = Jcol_sensible(:)/(tf-ti)
  Qadvcum(:) = Qadvcum(:)/(tf-ti)


  do kk=1, mp
     tmp1d1(kk) = (sum(vsnow(kk)%Jsensible(:)) + sum(vsnow(kk)%Jlatent(:)))
     ! heat stored in snowpack
     where (vsnow(kk)%hsnow(:).gt.zero)
        vsnow(kk)%Jsensible(:) = (vsnow(kk)%hsnow(:)-vsnow(kk)%hliq(:))*rhow*csice*(vsnow(kk)%Tsn(:)) + &
             vsnow(kk)%hliq(:)*rhow*cswat*(vsnow(kk)%Tsn(:))

        vsnow(kk)%Jlatent(:) = (vsnow(kk)%hsnow(:)-vsnow(kk)%hliq(:))*rhow*(-lambdaf)
     elsewhere
        vsnow(kk)%Jsensible(:) = zero
        vsnow(kk)%Jlatent(:) = zero
     endwhere
     deltaEsnow(kk) = sum(vsnow(kk)%Jsensible(:)) + sum(vsnow(kk)%Jlatent(:)) - &
          tmp1d1(kk)



     deltah0(kk) = h0(kk)-h0old(kk)+sum(vsnow(kk)%hsnow(1:vsnow(kk)%nsnow))-hsnowold(kk)
  enddo
  ssoil%thetai = thetai(:,:)
  ip(:)  = sum(ssoil%thetai*dx(:,:),2)   + h0(:)*ssoil%thetai(:,1)/par(:,1)%thre   ! ice in profile at tf
  ! water at tf
  wp(:)  = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*dx(:,:),2) + plit(:)%thre*SL(:)*dxL(:)
  win(:) = win(:) + (qprec(:)+qprec_snow(:))*(tf-ti)



  if (1 == 1) then

     k=1

     write(332,"(i8,i8,16e16.6)") ktau,nsteps(k),wp(k)-wpi(k),infil(k)-drn(k),runoff(k),&
          win(k)-(wp(k)-wpi(k)+deltah0(k)+runoff(k)+evap(k)+drn(k))-Etrans(k)*dt,wp(k),evap(k),evap_pot(k),infil(k),&
          drn(k),h0(k),Etrans(k)*dt,discharge(k), fws(k), (ip(k)-ipi(k)), fsat(k)
     write(334,"(100f15.6)") S(k,:)
     write(336,"(100f15.6)") Tsoil(k,:)
     write(335,"(100f20.6)") vmet(k)%Ta, T0(k), rh0(k), H(k), lE(k), &
          G0(k),Jcol_sensible(k),Jcol_latent_S(k), Jcol_latent_T(k), &
          vmet(k)%Rn, TL(k), SL(k), deltaice_cum_T(k), &
          deltaice_cum_S(k), rhsurface(k), Tsurface(k), vmet(k)%rha, &
           Qadvcum(k), sum((Jsensible(k,:)-ssoil%gammzz(k,:)),1)
     write(338,"(100f18.6)") thetai(k,:)


  endif

  tmp1d1 = sum((Jsensible(1,:)-ssoil%gammzz(1,:)),1)-Jcol_sensible(1)*3600

  ! Update variables for output:
  ssoil%tss       = real(Tsurface(:)) + 273.16
  ssoil%tgg       = real(Tsoil(:,:)) + 273.16
  ssoil%wb        = real(S(:,:)*(par(:,:)%thr+(par(:,:)%the-par(:,:)%thr)))
  ssoil%wbice     = thetai
  ssoil%wbtot     = real(wp(:)*thousand)
  canopy%ga       = real(G0(:))
  canopy%fhs      = real(H(:))
  canopy%fes      = real(lE(:))
  ssoil%hflux     = qh(:,:)
  ssoil%rnof1     = real(runoff(:)*thousand)/dt + runoff_sat
  ssoil%rnof2     = real(drn(:)*thousand)/dt !+ qb
  ssoil%runoff = ssoil%rnof1 + ssoil%rnof2
  ssoil%zdelta    = zdelta(:)
  ssoil%S         = S(:,:)
  ssoil%SL        = SL(:)
  ssoil%TL        = TL(:)
  ssoil%delwcol   = (wp(:)-wpi(:)+deltah0(:))*thousand  ! includes cange in snow pack
  ssoil%Tsurface  = Tsurface(:)
  ssoil%rh0       = rh0(:)
  ssoil%rhsurface = rhsurface(:)
  ssoil%lE        = lE(:)
  ssoil%evap      = evap(:)*thousand
  ssoil%rex       = wex(:,:)*thousand
  ssoil%kth       = kth(:,:)
  ssoil%nsteps = real(nsteps)
  canopy%fwsoil = real(fws(:))

 ! write(*,"(a10, 100f16.6)") "sli_main full", canopy%fhs, canopy%fes, ssoil%rtsoil, vmet(1)%Ta

  if (litter==0) then
     ssoil%rlitt  = zero
  else
     ssoil%rlitt  = dxL(:)/vlit(:)%Dv
  endif
  ssoil%gammzz    = Jsensible(:,:)
  ssoil%h0        = h0(:)

  ! update CABLE snow variables
  do kk=1, mp
     ssoil%snowd(kk) = vsnow(kk)%wcol*thousand ! amount of snow  (mm liq water eq)
     ssoil%smass(kk,1:nsnow_max) = vsnow(kk)%hsnow(:)*thousand ! amount of snow in dedicated snow pack (mm liq water eq)

     ssoil%sdepth(kk,1:nsnow_max) = vsnow(kk)%depth(:) ! depth of snow pack (m)
     ssoil%ssdn(kk,1:nsnow_max) = vsnow(kk)%dens(:) ! density of snow (kg m-3)
     ssoil%tggsn(kk,1:nsnow_max) = vsnow(kk)%tsn(:) + Tzero  ! abs T of snowpack
     ssoil%sconds(kk,1:nsnow_max) = vsnow(kk)%kH(:) ! thermal conductivty of snowpack
     ssoil%snowliq(kk,1:nsnow_max) = vsnow(kk)%hliq(:)*thousand ! amount of liq snow water

     ssoil%smelt(kk) = vsnow(kk)%Qmelt*thousand/dt ! amount of melted snow leaving bottom of snow pack (mm/dt)
     ssoil%nsnow(kk) = vsnow(kk)%nsnow
     if (sum(ssoil%sdepth(kk,1:nsnow_max)).gt.zero) then
        ssoil%ssdnn(kk) = ssoil%snowd(kk)/sum(ssoil%sdepth(kk,1:nsnow_max))
     endif
  enddo


  ssoil%isflag = 0


  ! snow output
  if (1 == 1) then
     k = 1
     write(340,"(100e16.6)") sum(vsnow(k)%hsnow(1:vsnow(k)%nsnow)), vsnow(k)%tsn(1),sum(vsnow(k)%hliq(1:vsnow(k)%nsnow)), &
	      qprec_snow(k)*dt, vsnow(k)%Qmelt, qprec(k)*dt, &
          vsnow(k)%Qevap,vsnow(k)%Qvap,ssoil%albsoilsn(k,1), ssoil%albsoilsn(k,2), ssoil%sconds(k,1), &
		  vsnow(k)%dens(1),sum(vsnow(k)%depth(1:vsnow(k)%nsnow)), vsnow(k)%J, &
		  vsnow(k)%MoistureFluxDivergence, vsnow(k)%FluxDivergence, vsnow(k)%dens(nsnow_max), vsnow(k)%tsn(nsnow_max)
  endif



  canopy%ofes(:) = canopy%fes(:)
  ! Update total latent heat to reflect updated soil component:
  canopy%fe(:) = real(canopy%fev(:)) + canopy%fes(:)
  ! Update total sensible heat to reflect updated soil component:
  canopy%fh(:) = real(canopy%fhv(:)) + canopy%fhs(:)

 endif ! SEB only

END SUBROUTINE sli_main
