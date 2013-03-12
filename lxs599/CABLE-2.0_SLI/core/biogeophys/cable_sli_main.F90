SUBROUTINE sli_main(irec, dt, veg, soil, ssoil, met, canopy, air)

  ! Main subroutine for Soil-litter-iso soil model
  ! Vanessa Haverd, CSIRO Marine and Atmospheric Research
  ! and Matthias Cuntz, UFZ - Helmholtz Centre for Environmental Research, 2010
  ! Modified to operate for multiple veg tiles but a single soil column, March 2011
  ! Rewritten for same number of soil columns as veg tiles May 2012
  USE cable_def_types_mod,       ONLY: veg_parameter_type, soil_parameter_type, soil_snow_type, met_type, &
       canopy_type, air_type, ms, mp, r_1, r_2, i_d 
  USE physical_constants, ONLY: tfrz
  USE sli_numbers,        ONLY: zero, half, one, four, thousand, & ! numbers
       Tzero, &                                       ! variables
       vars_met, vars, params, vars_snow, &                                  ! types
       count_sparse, count_hyofS, &
	   MW, snmin, Rgas, Lambdas, lambdaf, csice, cswat, rhow
  USE sli_utils,          ONLY: x, dx, par, setpar, setx, plit, dxL, setlitterpar, esat, &
                                esat_ice, slope_esat_ice
  USE sli_roots,          ONLY: setroots, getrex
  USE sli_solve,          ONLY: solve
  USE cable_soil_snow_module, ONLY: snowdensity

  IMPLICIT NONE

  INTEGER(i_d),              INTENT(IN)    :: irec
  REAL(r_1),                 INTENT(IN)    :: dt
  TYPE(veg_parameter_type),  INTENT(INOUT) :: veg     ! all r_1
  TYPE(soil_parameter_type), INTENT(INOUT) :: soil    ! all r_1
  TYPE(soil_snow_type),      INTENT(INOUT) :: ssoil   ! r_1, r_2 desaster
  TYPE(met_type),            INTENT(INOUT) :: met     ! all r_1
  TYPE(canopy_type),         INTENT(INOUT) :: canopy  ! all r_1
  TYPE(air_type),            INTENT(INOUT) :: air     ! all r_1

  REAL(r_2),    PARAMETER :: emsoil=0.97
  REAL(r_2),    PARAMETER :: rhocp=1.1822e3
  REAL(r_2),    PARAMETER :: Dva = 2.17e-5
  INTEGER(i_d) :: i, j, k, kk, setroot
  REAL(r_2)    :: ti, tf
  TYPE(vars_met), DIMENSION(1:mp)      :: vmet ! Meteorology above soil
  TYPE(vars),     DIMENSION(1:mp)      :: vlit
  TYPE(vars_snow),     DIMENSION(1:mp) :: vsnow
  INTEGER(i_d),   DIMENSION(1:mp)      :: nsteps
  REAL(r_2),      DIMENSION(1:mp,1:ms) :: Tsoil, S, thetai, Jsensible
  REAL(r_2),      DIMENSION(1:mp)      :: SL, TL, T0
  REAL(r_2),      DIMENSION(1:mp)      :: drn, evap, infil, qprec,qprec_snow, runoff, runoff_sat
  REAL(r_2),      DIMENSION(1:mp)      :: win, wp, wpi, h0,deltah0, h0old,hsnowold, discharge
  REAL(r_2),      DIMENSION(1:mp)      :: ip, ipi ! volumetric ice content of profile (final and initial)
  REAL(r_2),      DIMENSION(1:mp,1:ms) ::   wex, kth, qex
  REAL(r_2),      DIMENSION(1:mp,1:ms) ::  FS
  REAL(r_2),      DIMENSION(1:mp)      :: rh0, rhsurface
  ! surface temperature (top of top soil layer or top of litter layer)
  REAL(r_2),      DIMENSION(1:mp)      :: Tsurface
  REAL(r_2),      DIMENSION(1:mp)      :: gr, grc
  REAL(r_2),      DIMENSION(1:mp)      :: Etrans
  REAL(r_2),      DIMENSION(1:mp)      :: gamma
  REAL(r_2),      DIMENSION(1:mp)      :: G0, H, lE
  REAL(r_2),      DIMENSION(1:mp,0:ms) :: qh, qvsig, qlsig, qvTsig, qvh
  REAL(r_2),      DIMENSION(1:mp)      :: deltaTa, lE_old, SA, SB, wpAi, wpBi, wpA, wpB
  REAL(r_2),      DIMENSION(1:mp)      :: evap_pot, deltaice_cum_T, deltaice_cum_S, zdelta
  REAL(r_2),      DIMENSION(1:mp)      ::fws
  REAL(r_2),      DIMENSION(1:mp)      :: Qadvcum,Jcol_sensible,Jcol_latent_S,Jcol_latent_T
  REAL(r_2),      DIMENSION(1:mp)      :: tmp1d1
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
  INTEGER(i_d), PARAMETER :: testcase     = 0 ! isotopic test cases 1-8
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
  LOGICAL, SAVE :: first = .false.


    ! output files for testing purposes
  if (first) then
     open (unit=332,file="vh08.out",status="replace",position="rewind")
     open (unit=334,file="S.out",status="replace",position="rewind")
     open (unit=336,file="Tsoil.out",status="replace",position="rewind")
     open (unit=335,file="SEB.out",status="replace",position="rewind")
     open(unit=338, file="thetai.out", status="replace", position="rewind")
     open(unit=339, file="test_dt.out",status="replace", position="rewind")
	 open(unit=340, file="snow.out", status="replace", position="rewind")
	 first = .false.
  endif

  fsat = 0.0
  zdelta = 0.0

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
     soil%swilt_vec(:,:) = merge(spread(real(soil%swiltB(:),r_2),2,ms), &
          soil%swilt_vec(:,:), par(:,:)%isbhorizon==one)
     soil%ssat_vec(:,:) = merge(spread(real(soil%ssatB(:),r_2),2,ms) &
          ,soil%ssat_vec(:,:), par(:,:)%isbhorizon==one)
     soil%sfc_vec(:,:) = merge(spread(real(soil%sfcB(:),r_2),2,ms) &
          ,soil%sfc_vec(:,:), par(:,:)%isbhorizon==one)
  endif

  ! If we want solutes:
!!$     if (.not. allocated(bd)) allocate(bd(soil%nhorizons(k)))

  ! Litter parameters:
  if (.not. allocated(plit)) then
	 index=(/(i,i=1,mp,1)/)
     call setlitterpar(mp, soil,index)
  endif

  ! Met data above soil:
  vmet(:)%Ta    = real(met%Tvair-tfrz,r_2) 
  vmet(:)%Da    = real(met%dva,r_2)
  vmet(:)%rbh   = ssoil%rtsoil
  vmet(:)%rbw   = vmet(:)%rbh
  vmet(:)%rha   = max(min((esat(vmet(:)%Ta)-vmet(:)%Da)/esat(vmet(:)%Ta),one),0.1_r_2)
  vmet(:)%cva   = vmet(:)%rha * esat(vmet(:)%Ta)*0.018_r_2/thousand/8.314_r_2/(vmet(:)%Ta+Tzero) ! m3 H2O (liq) m-3 (air)
  vmet(:)%phiva = Dva * vmet(:)%cva
  vmet(:)%Rn    = canopy%fns
  Etrans(:)     = canopy%fevc/air%rlam/thousand ! m s-1
  h0(:) = ssoil%h0

  ! Initialisations:
  if (irec == 1) then
     rh0(:)         = vmet(:)%rha
     rhsurface(:)   = rh0(:) ! initialise rel. humidity at surface
     Tsurface(:)    = vmet(:)%Ta
     T0(:)          = Tsurface(:)
     deltaTa(:)     = zero
     lE_old(:)      = zero
     ssoil%h0(:)    = zero
     ssoil%thetai(:,:) = zero
     met%tk_old(:)     = met%tk(:)
	 met%qv_old(:)     = met%qv(:)
     ! do not seem to be initialised elsewhere
     ssoil%snowd(:) = 0.0_r_1
	 ssoil%smass(:,:) = 0.0_r_1
     zdelta(:)      = x(:,ms)
     ssoil%gammzz = zero
     !ssoil%tgg(:,6)=5.0+Tfrz
     count_sparse = 0
     count_hyofS = 0

     ssoil%ssdn(:,:) = 120_r_2 ! snow density kg m-3
     ssoil%sconds(:,:) = 0.06_r_2    ! snow thermal cond (W m-2 K-1)
	 ssoil%S  = min(one,ssoil%S)

  else
     rhsurface(:)   = ssoil%rhsurface
     Tsurface(:)    = ssoil%Tsurface
     T0(:)          = ssoil%Tsurface
     rh0(:)         = ssoil%rh0
     deltaTa(:)     = zero
     lE_old(:)      = ssoil%lE
     zdelta(:)      = ssoil%zdelta
  endif

  SL(:)       = zero   ! degree of litter saturation
  TL(:)       = zero  ! litter T
  Tsoil(:,:)  = real(ssoil%tgg-tfrz,r_2)
  S(:,:)         = ssoil%S                ! degree of soil saturation
  thetai = ssoil%thetai
  ssoil%smelt = zero

  do kk=1, mp
     if (ssoil%sdepth(kk,1).gt.zero)  then ! snow
        ! define variables associated with snow
	    ssoil%isflag(kk) = 1
		vsnow(kk)%hsnow(:) = ssoil%smass(kk,:)/thousand 
		vsnow(kk)%depth(:) = ssoil%sdepth(kk,:)  ! depth of snow pack (m)
	    vsnow(kk)%dens(:) = 120_r_2 ! snow density kg m-3
		vsnow(kk)%dens(:) =ssoil%ssdn(kk,:)
        vsnow(kk)%wcol= ssoil%snowd(kk)/thousand 

	    vsnow(kk)%tsn(:) = ssoil%tggsn(kk,:) - Tzero
	    !vsnow(kk.:)%kH = 0.06_r_2    ! snow thermal cond (W m-2 K-1)
		vsnow(kk)%kH(:) =  ssoil%sconds(kk,:) 
	    vsnow(kk)%Dv(:) = Dva*(ssoil%tggsn(kk,:)/Tzero)**1.88_r_2 ! m2 s-1
	    vsnow(kk)%sl(:) = slope_esat_ice(vsnow(kk)%tsn(:)) * Mw/thousand/Rgas/(vsnow(kk)%tsn(:)+Tzero)
	    vsnow(kk)%kE(:)     = vsnow(kk)%Dv(:)*vsnow(kk)%sl(:)*thousand*lambdaf
	    vsnow(kk)%kth(:) = vsnow(kk)%kE(:) + vsnow(kk)%kH(:)
        vsnow(kk)%cv(:) = esat_ice(vsnow(kk)%tsn(:))*Mw/thousand/Rgas/(vsnow(kk)%tsn(:)+Tzero) ! m3 m-3
        vsnow(kk)%hliq(:) = ssoil%snowliq(kk,:)/thousand ! amount of liq snow water
        vsnow(kk)%melt(:) = zero ! amount of melted snow leaving bottom of snow pack (mm/dt)
		! heat stored in snowpack
	    where (vsnow(kk)%hsnow(:).gt.zero)  
	       vsnow(kk)%Jsensible(:) = (vsnow(kk)%hsnow(:)-vsnow(kk)%hliq(:))*rhow*csice*(vsnow(kk)%Tsn(:)+Tzero) + &
	            vsnow(kk)%hliq(:)*rhow*cswat*(vsnow(kk)%Tsn(:)+Tzero)

	       vsnow(kk)%Jlatent(:) = (vsnow(kk)%hsnow(:)-vsnow(kk)%hliq(:))*rhow*(-lambdaf)
	    endwhere
		if (vsnow(kk)%hsnow(3).gt.zero) then
		  vsnow(kk)%nsnow = 3
		elseif (vsnow(kk)%hsnow(2).gt.zero) then
		  vsnow(kk)%nsnow = 2
		elseif (vsnow(kk)%hsnow(1).gt.zero) then
		  vsnow(kk)%nsnow = 1
		else
		  vsnow(kk)%nsnow = 0
		endif
		   
 
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
	    vsnow(kk)%cv(:) = zero
	    vsnow(kk)%hliq(:) = zero
	    vsnow(kk)%melt(:) = zero
		
      endif
	 ssoil%osnowd = ssoil%snowd(kk)
   enddo  

  deltaTa(:)  = zero
  lE_old(:)   = ssoil%lE
  gamma(:)    = real(veg%gamma(:),r_2)

  gr(:)       = four * emsoil * (vmet(:)%Ta+tfrz)**3 *5.67e-8_r_2 ! radiation conductance Wm-2K-1
  grc(:)      = one/vmet(:)%rbh + gr(:)/rhocp
  vmet(:)%rrc = one/grc(:)                            ! resistance to radiative and convective heat transfer
  qprec(:)    = (canopy%through-met%precip_sn)/thousand/dt              ! liq precip rate (m s-1)
  qPrec_snow(:) = (met%precip_sn)/thousand/dt 
  
 where (vmet%Ta<zero)

!qPrec_snow(:) = merge(one/thousand/dt,zero,irec<100) 
!qPrec_snow(:) = zero
!qPrec_snow(:) =qprec_snow(:)*10.


  endwhere
  h0old(:)    = ssoil%h0 ! pond height
  do kk=1, mp
     !hsnowold(kk) = sum(ssoil%smass(kk,:)/thousand) ! water in decicated snow layer(s)
     hsnowold(kk) = ssoil%smass(kk,1)/thousand
  enddo

  ! Heat balance variables

  ! Water balance variables:
  ipi(:)   = sum(ssoil%thetai*dx(:,:),2)  + h0(:)*ssoil%thetai(:,1)/par(:,1)%thre       ! ice in profile initially
  ! water in profile initially
  wpi(:)   = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*dx(:,:),2)  + plit(:)%thre*SL(:)*dxL(:)
  wpAi(:)  = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*par(:,:)%isahorizon*dx(:,:),2) + &
       plit(:)%thre*SL(:)*dxL(:) ! water in profile initially
  wpBi(:)  = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*par(:,:)%isbhorizon*dx(:,:),2)

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
          soil%swilt_vec(kk,:), real(canopy%fevc(kk)/air%rlam(kk),r_2)/thousand , gamma(kk))
  enddo

  do k=1,ms
     qex(:,k)= ssoil%rex(:,k)
  enddo

  ! calculate base flow (topmodel)
  !qb = alpha * par(:,ms)%Ke*exp(-(par(:,ms)%zeta*zdelta))
  qex(:,ms) =  qex(:,ms) !+qb

  if (irec.eq.1660) then ! for debugging
      write (*,*) irec, vsnow(1)%hsnow(:), vsnow(1)%wcol
   endif

  call solve(ti, tf, irec, mp, qprec(:),qprec_snow(:), ms,  dx(:,:), &
       h0(:), S(:,:), thetai(:,:), Jsensible(:,:), Tsoil(:,:), evap(:), &
       evap_pot(:), runoff(:), infil(:), drn(:), discharge(:), qh(:,:), &
       nsteps, vmet(:), vlit(:),vsnow(:),kth(:,:),T0(:), rh0(:), Tsurface(:), &
       rhsurface(:), H(:), lE(:), G0(:),Qadvcum(:),Jcol_sensible(:), &
       Jcol_latent_S(:),Jcol_latent_T(:), deltaice_cum_T(:), &
       deltaice_cum_S(:), dxL(:), zdelta(:), SL(:), TL(:), &
       plit(:), par(:,:), qex=qex(:,:), &
       wex=wex(:,:), FS=FS(:,:), qvsig=qvsig(:,:), qlsig=qlsig(:,:), qvTsig=qvTsig(:,:), qvh=qvh(:,:), &
       deltaTa=deltaTa(:), lE_old=lE_old(:), &
       dolitter=litter, doisotopologue=isotopologue, dotestcase=testcase, dosepts=septs, docondition=condition, &
       doadvection=advection)

  if (count_sparse.lt.count_hyofS) then
     write(*,*) 'count_sparse.lt.count_hyofS: ', irec, count_sparse, count_hyofS
  endif

  H(:)      = H(:)/(tf-ti)
  lE(:)     = lE(:)/(tf-ti)
  G0(:)     = G0(:)/(tf-ti)
  Jcol_latent_S(:) = Jcol_latent_S(:)/(tf-ti)
  Jcol_latent_T(:) = Jcol_latent_T(:)/(tf-ti)
  Jcol_sensible(:) = Jcol_sensible(:)/(tf-ti)
  Qadvcum(:) = Qadvcum(:)/(tf-ti)

  do kk=1, mp
     deltah0(kk) = h0(kk)-h0old(kk)+sum(vsnow(kk)%hsnow(:))-hsnowold(kk)
  enddo
  ssoil%thetai = thetai(:,:)
  ip(:)  = sum(ssoil%thetai*dx(:,:),2)   + h0(:)*ssoil%thetai(:,1)/par(:,1)%thre   ! ice in profile at tf
  ! water at tf
  wp(:)  = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*dx(:,:),2) + plit(:)%thre*SL(:)*dxL(:) 
  win(:) = win(:) + (qprec(:)+qprec_snow(:))*(tf-ti)
  wpA(:) = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*par(:,:)%isahorizon*dx(:,:),2) + &
       plit(:)%thre*SL(:)*dxL(:) ! water in profile initially
  wpB(:) = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*par(:,:)%isbhorizon*dx(:,:),2)
  SA(:)  = (wpA(:)/sum(dx(:,:)*par(:,:)%isahorizon,2) - par(:,1)%thr)/(par(:,1)%the - par(:,1)%thr)
  SB(:)  = (wpB(:)/sum(dx(:,:)*par(:,:)%isbhorizon,2) - par(:,ms)%thr)/(par(:,ms)%the - par(:,ms)%thr)
  where (all(par(:,:)%isbhorizon==zero,2)) SB(:) = zero


  if (1 == 1) then
     k=1
     write(332,"(i8,i8,16e16.6)") irec,nsteps(k),wp(k)-wpi(k),infil(k)-drn(k),runoff(k),&
          win(k)-(wp(k)-wpi(k)+deltah0(k)+runoff(k)+evap(k)+drn(k))-Etrans(k)*dt,wp(k),evap(k),evap_pot(k),infil(k),&
          drn(k),h0(k),Etrans(k)*dt,discharge(k), fws(k), (ip(k)-ipi(k)), fsat(k)
     write(334,"(100f15.6)") S(k,:)
     write(336,"(100f15.6)") Tsoil(k,:)
     write(335,"(100f20.6)") vmet(k)%Ta, T0(k), rh0(k), H(k), lE(k), &
          G0(k),Jcol_sensible(k),Jcol_latent_S(k), Jcol_latent_T(k), &
          vmet(k)%Rn, TL(k), SL(k), deltaice_cum_T(k), &
          deltaice_cum_S(k), rhsurface(k), Tsurface(k), vmet(k)%rha, &
          kth(k,:), Qadvcum(k), sum((Jsensible(k,:)-ssoil%gammzz(k,:)),1)
     write(338,"(100f18.6)") thetai(k,:)
  endif

  tmp1d1 = sum((Jsensible(1,:)-ssoil%gammzz(1,:)),1)-Jcol_sensible(1)*3600

  ! Update variables for output:
  ssoil%tss       = real(Tsurface(:),r_1) + tfrz
  ssoil%tgg     = real(Tsoil(:,:),r_1) + tfrz
  ssoil%wb      = real(S(:,:)*(par(:,:)%thr+(par(:,:)%the-par(:,:)%thr)),r_1)
  ssoil%wbtot     = real(wp(:)*thousand,r_1)
  canopy%ga       = real(G0(:),r_1)
  canopy%fhs      = real(H(:),r_1)
  canopy%fes      = real(lE(:),r_1)
  ssoil%hflux   = qh(:,:)
  ssoil%rnof1     = real(runoff(:)*thousand,r_1) + runoff_sat
  ssoil%rnof2     = real(discharge(:)*thousand,r_1) !+ qb
  ssoil%zdelta   = zdelta(:)
  ssoil%S       = S(:,:)
  ssoil%SL        = SL(:)
  ssoil%TL        = TL(:)
  ssoil%delwcol   = (wp(:)-wpi(:)+deltah0(:))*thousand 
  ssoil%Tsurface  = Tsurface(:)
  ssoil%rh0      = rh0(:)
  ssoil%rhsurface = rhsurface(:)
  ssoil%lE       = lE(:)
  ssoil%evap      = evap(:)*thousand
  ssoil%SA        = SA(:)
  ssoil%SB        = SB(:)
  ssoil%delwcolA  = (wpA(:)-wpAi(:))*thousand
  ssoil%delwcolB  = (wpB(:)-wpBi(:))*thousand
  ssoil%rex       = wex(:,:)*thousand
  if (litter==0) then
     ssoil%rlitt  = zero
  else
     ssoil%rlitt  = dxL(:)/vlit(:)%Dv
  endif
  ssoil%rexA      = sum(wex(:,:)*par(:,:)%isahorizon,2)*thousand
  ssoil%rexB      = sum(wex(:,:)*par(:,:)%isbhorizon,2)*thousand

  ssoil%kth = kth(:,:)
  ssoil%gammzz = Jsensible(:,:)
  ssoil%h0 = h0(:)
  do j=1, ms
     where(par(:,j)%isahorizon==one) ssoil%leachAB = (qvsig(:,j)+qlsig(:,j))*thousand
  enddo

  ! update CABLE snow variables
  do kk=1, mp
  ssoil%snowd(kk) = vsnow(kk)%wcol*thousand ! amount of snow  (mm liq water eq)
  ssoil%smass(kk,:) = vsnow(kk)%hsnow(:)*thousand ! amount of snow in dedicated snow pack (mm liq water eq)

  ssoil%sdepth(kk,:) = vsnow(kk)%depth(:) ! depth of snow pack (m)
  ssoil%ssdn(kk,:) = vsnow(kk)%dens(:) ! density of snow (kg m-3)
  ssoil%tggsn(kk,:) = vsnow(kk)%tsn(:) + Tzero  ! abs T of snowpack
  ssoil%sconds(kk,:) = vsnow(kk)%kH(:) ! thermal conductivty of snowpack
  ssoil%snowliq(kk,:) = vsnow(kk)%hliq(:)*thousand ! amount of liq snow water
  
  ssoil%smelt(kk) = sum(vsnow(kk)%melt(:)*thousand) ! amount of melted snow leaving bottom of snow pack (mm/dt)
  enddo
  !where (ssoil%sdepth(:,1).gt.zero)  
  !   ssoil%isflag = 1
  !elsewhere
  !   ssoil%isflag = 0
  !endwhere

  ssoil%isflag = 0
  !CALL snow_albedo(dt, irec, met, ssoil, soil )
  CALL snowdensity (dt, ssoil, soil)

  !ssoil%sconds = 0.07 ! test vh

  ! snow output
  if (1 == 1) then
    k = 1
        write(340,"(100e16.6)") vsnow(k)%wcol, vsnow(k)%tsn(1), vsnow(k)%hliq(1), qprec_snow(k)*dt, vsnow(k)%melt(1), qprec(k)*dt, &
	 ssoil%albsoilsn(k,1), ssoil%albsoilsn(k,2), ssoil%sconds(k,1), ssoil%ssdnn(k),vsnow(k)%Qadv_vap(1), vsnow(k)%Qadv_rain, &
	 vsnow(k)%Qadv_snow, vsnow(k)%Qcond_net(1),vsnow(k)%deltaJlatent(1),vsnow(k)%deltaJsensible(1),vsnow(k)%Qadv_melt(1), &
	 vsnow(k)%Qadv_transfer(1),vsnow(k)%FluxDivergence(1),vsnow(k)%Jlatent(1), vsnow(k)%Jsensible(1)
  endif

  canopy%ofes(:) = canopy%fes(:)
  ! Update total latent heat to reflect updated soil component:
  canopy%fe(:) = real(canopy%fev(:),r_1) + canopy%fes(:)
  ! Update total sensible heat to reflect updated soil component:
  canopy%fh(:) = real(canopy%fhv(:),r_1) + canopy%fhs(:)

  ! store air temperature
  met%tk_old(:)  = met%tk(:)
  met%qv_old(:)  = met%qv(:)


END SUBROUTINE sli_main
