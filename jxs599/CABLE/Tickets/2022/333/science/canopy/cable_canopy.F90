MODULE cable_canopy_module

! physical constants
USE cable_phys_constants_mod, ONLY : CTFRZ   => TFRZ
USE cable_phys_constants_mod, ONLY : CRMAIR  => RMAIR
USE cable_phys_constants_mod, ONLY : CRGAS   => RGAS
USE cable_phys_constants_mod, ONLY : CDHEAT  => DHEAT
USE cable_phys_constants_mod, ONLY : CZETNEG => ZETNEG
USE cable_phys_constants_mod, ONLY : CZETMUL => ZETMUL
USE cable_phys_constants_mod, ONLY : CZETPOS => ZETPOS
USE cable_phys_constants_mod, ONLY : CGRAV   => GRAV
USE cable_phys_constants_mod, ONLY : CUMIN   => UMIN
USE cable_phys_constants_mod, ONLY : CRHOW   => RHOW
USE cable_phys_constants_mod, ONLY : CCTL    => CTL
USE cable_phys_constants_mod, ONLY : CCSW    => CSW
USE cable_phys_constants_mod, ONLY : CEMLEAF => EMLEAF
USE cable_phys_constants_mod, ONLY : CEMSOIL => EMSOIL
USE cable_phys_constants_mod, ONLY : CSBOLTZ => SBOLTZ
USE cable_phys_constants_mod, ONLY : CPRANDT => PRANDT
USE cable_phys_constants_mod, ONLY : CCAPP   => CAPP
USE cable_phys_constants_mod, ONLY : CRMH2O  => RMH2O
USE cable_phys_constants_mod, ONLY : CAPOL   => APOL
USE cable_phys_constants_mod, ONLY : CA33    => A33
USE cable_phys_constants_mod, ONLY : CVONK   => VONK
USE cable_phys_constants_mod, ONLY : CZETA0  => ZETA0
USE cable_phys_constants_mod, ONLY : CTETENA     => TETENA
USE cable_phys_constants_mod, ONLY : CTETENB     => TETENB
USE cable_phys_constants_mod, ONLY : CTETENC     => TETENC
USE cable_phys_constants_mod, ONLY : CTETENA_ICE => TETENA_ICE
USE cable_phys_constants_mod, ONLY : CTETENB_ICE => TETENB_ICE
USE cable_phys_constants_mod, ONLY : CTETENC_ICE => TETENC_ICE
! photosynthetic constants
USE cable_photo_constants_mod, ONLY : CRGSWC => RGSWC
USE cable_photo_constants_mod, ONLY : CGAM0  => GAM0
USE cable_photo_constants_mod, ONLY : CGAM2  => GAM2
USE cable_photo_constants_mod, ONLY : CRGBWC => RGBWC
USE cable_photo_constants_mod, ONLY : CGAM1  => GAM1
USE cable_photo_constants_mod, ONLY : CTREFK => TREFK
USE cable_photo_constants_mod, ONLY : CMAXITER  => MAXITER ! only integer here
! maths & other constants
USE cable_math_constants_mod,  ONLY : CPI_C  => PI
USE cable_other_constants_mod, ONLY : CLAI_THRESH  => LAI_THRESH

  IMPLICIT NONE

  PUBLIC define_canopy
  PRIVATE

CONTAINS


SUBROUTINE define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy,climate, sunlit_veg_mask, reducedLAIdue2snow )
    USE cable_def_types_mod
   USE cbl_radiation_module, ONLY : radiation
    USE cable_air_module
    USE cable_common_module
    USE cable_roughness_module

USE cable_canopy_module_subrs_module, ONLY:  Surf_wetness_fact, dryLeaf,       &
                                             photosynthesis, fwsoil_calc_std,  &
                                             fwsoil_calc_non_linear,           &
                                             fwsoil_calc_Lai_Ktaul,            &
                                             fwsoil_calc_sli,  getrex_1d

USE cbl_friction_vel_module, ONLY: comp_friction_vel, psim

    TYPE (balances_type), INTENT(INOUT)  :: bal
    TYPE (radiation_type), INTENT(INOUT) :: rad
    TYPE (roughness_type), INTENT(INOUT) :: rough
    TYPE (air_type), INTENT(INOUT)       :: air
    TYPE (met_type), INTENT(INOUT)       :: met
    TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
    TYPE (canopy_type), INTENT(INOUT)    :: canopy
    TYPE (climate_type), INTENT(IN)    :: climate

    TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)    :: veg

REAL :: reducedLAIdue2snow(mp)
logical :: sunlit_veg_mask(mp) 
    REAL, INTENT(IN)               :: dels ! integration time setp (s)
    INTEGER  ::                                                                 &
         iter,  & ! iteration #
         iterplus !

    REAL, DIMENSION(mp) ::                                                      &
         rt0,           & ! turbulent resistance
         ortsoil,       & ! turb. resist. prev t-step
         rt1usc,        & ! eq. 3.53, SCAM manual, 1997
         tstar,         & !
         zscrn,         & !
         qstar,         & !
         rsts,          & !
         qsurf,         & !
         qtgnet,        & !
         tss4,          & ! soil/snow temperature**4
         qstvair,       & ! sat spec humidity at leaf temperature
         xx,            & ! delta-type func 4 sparse canopy limit, p20 SCAM manual
         r_sc,          & !
         zscl,          & !
         pwet,          & !
         dq,            & ! sat sp
         dq_unsat,      & ! spec hum diff including rh at srf
         xx1,           & !
         sum_rad_rniso, & !
         sum_rad_gradis,& !
         rttsoil,       & ! REV_CORR working variable for sensitivity terms
         rhlitt,        & ! REV_CORR working variables for litter resistances
         relitt,        & !
         alpm1,         & ! REV_CORR working variables for Or scheme
         beta2,         & ! beta_div_alpm1 = beta2/alpm1 (goes to zero without
         beta_div_alpm    ! division when no canopy)

    ! temporary buffers to simplify equations
    REAL, DIMENSION(mp) ::                                                      &
         ftemp,z_eff,psim_arg, psim_1, psim_2, rlower_limit,                      &
         term1, term2, term3, term5
    ! arguments for potential_evap (sli)
    REAL(r_2), DIMENSION(mp) ::  Rn, rbh, rbw, Ta, rha,Ts, &
         Tsoil, Epot, Hpot, Gpot, &
         kth, dz,lambdav, &
         dEdrha, dEdTa, dEdTsoil, dGdTa, dGdTsoil
    REAL, DIMENSION(mp) :: qsat

    REAL, DIMENSION(:), POINTER ::                                              &
         cansat,        & ! max canopy intercept. (mm)
         dsx,           & ! leaf surface vpd
         fwsoil,        & ! soil water modifier of stom. cond
         tlfx,          & ! leaf temp prev. iter (K)
         tlfy             ! leaf temp (K)

    REAL(r_2), DIMENSION(mp) ::                                                 &
         gbvtop                   ! bnd layer cond. top leaf

    REAL(r_2), DIMENSION(:), POINTER ::                                         &
         ecy,           & ! lat heat fl dry big leaf
         hcy,           & ! veg. sens heat
         rny,           & ! net rad
         ghwet             ! cond for heat for a wet canopy

    REAL(r_2), DIMENSION(:,:), POINTER ::                                       &
         gbhu,          & ! forcedConvectionBndryLayerCond
         gbhf,          & ! freeConvectionBndryLayerCond
         csx              ! leaf surface CO2 concentration

    REAL  :: rt_min
    REAL, DIMENSION(mp)       :: zstar, rL, phist, csw, psihat,rt0bus

    INTEGER :: j

    INTEGER, SAVE :: call_number =0

    ! END header

    call_number = call_number + 1

    !H! vNot sure that this is appropriate for JULES standalone - HaC
    !H!IF( .NOT. cable_runtime%um)                                                 &
         canopy%cansto =  canopy%oldcansto

    ALLOCATE( cansat(mp), gbhu(mp,mf))
    ALLOCATE( dsx(mp), fwsoil(mp), tlfx(mp), tlfy(mp) )
    ALLOCATE( ecy(mp), hcy(mp), rny(mp))
    ALLOCATE( gbhf(mp,mf), csx(mp,mf))
    ALLOCATE( ghwet(mp))

    ! BATS-type canopy saturation proportional to LAI:
    cansat = veg%canst1 * canopy%vlaiw

    !---compute surface wetness factor, update cansto, through
    CALL surf_wetness_fact( cansat, canopy, ssnow,veg,met, soil, dels )

    canopy%fevw_pot = 0.0
    canopy%gswx = 1e-3     ! default stomatal conuctance
    gbhf = 1e-3     ! default free convection boundary layer conductance
    gbhu = 1e-3     ! default forced convection boundary layer conductance
    ssnow%evapfbl = 0.0
    ssnow%rex = 0.0
    ! Initialise in-canopy temperatures and humidity:
    csx = SPREAD(met%ca, 2, mf) ! initialise leaf surface CO2 concentration
    met%tvair = met%tk
    met%qvair = met%qv
    canopy%tv = met%tvair
    canopy%fwsoil = 1.0

    CALL define_air (met, air)

    CALL qsatfjh(qstvair,met%tvair-Ctfrz,met%pmb)

    met%dva = (qstvair - met%qvair) *  Crmair/Crmh2o * met%pmb * 100.0
    dsx = met%dva     ! init. leaf surface vpd
    dsx= MAX(dsx,0.0)
    tlfx = met%tk  ! initialise leaf temp iteration memory variable (K)
    tlfy = met%tk  ! initialise current leaf temp (K)

    ortsoil = ssnow%rtsoil
    IF (cable_user%soil_struc=='sli') THEN
       ssnow%tss = REAL(ssnow%Tsurface) + Ctfrz
    ELSE
       ssnow%tss =  REAL((1-ssnow%isflag))*ssnow%tgg(:,1) +                    &
            REAL(ssnow%isflag)*ssnow%tggsn(:,1)
    ENDIF
    tss4 = ssnow%tss**4
    canopy%fes = 0.
    canopy%fess = 0.
    canopy%fesp = 0.
    ssnow%potev = 0.
    canopy%fevw_pot = 0.

    !L_REV_CORR - initialise sensitivity/ACCESS correction terms
    !NB %fes_cor is NOT initialised to zero at this point
    canopy%fhs_cor = 0.0
    canopy%fns_cor = 0.0
    canopy%ga_cor = 0.0
    canopy%fes_cor = 0.0

    !L_REV_CORR - new working variables
    rttsoil = 0.
    rhlitt = 0.
    relitt = 0.
    alpm1  = 0.
    beta2  = 0.

CALL radiation( ssnow, veg, air, met, rad, canopy, sunlit_veg_mask, &
  !constants
  clai_thresh, Csboltz, Cemsoil, Cemleaf, Ccapp &
)

    canopy%zetar(:,1) = CZETA0 ! stability correction terms
    canopy%zetar(:,2) = CZETPOS + 1
    canopy%zetash(:,1) = CZETA0 ! stability correction terms
    canopy%zetash(:,2) = CZETPOS + 1


    DO iter = 1, NITER

       ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
       ! resistances rt0, rt1 (elements of dispersion matrix):
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
       CALL comp_friction_vel(canopy%us, iter, mp, CVONK, CUMIN, CPI_C,      &
                             canopy%zetar, rough%zref_uv, rough%zref_tq,     &
                             rough%z0m, met%ua )

       ! E.Kowalczyk 2014
       IF (cable_user%l_new_roughness_soil)                                     &
        CALL ruff_resist( veg, rough, ssnow, canopy, veg%vlai, veg%hc, canopy%vlaiw )


       ! Turbulent aerodynamic resistance from roughness sublayer depth
       ! to reference height, x=1 if zref+disp>zruffs,
       ! 0 otherwise: thus rt1usc = 0 if zref+disp<zruffs
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
       xx = 0.5 + SIGN(0.5,rough%zref_tq+rough%disp-rough%zruffs)

       ! correction  by Ian Harman to the 2nd psis term
       rt1usc = xx * ( LOG( rough%zref_tq/MAX( rough%zruffs-rough%disp,         &
            rough%z0soilsn ) )               &
            - psis( canopy%zetar(:,iter) )                                  &
            + psis( canopy%zetar(:,iter) * ( MAX( rough%zruffs-rough%disp,  &
            rough%z0soilsn ) )            &
            / rough%zref_tq ) ) / CVONK
       rt_min = 5.

       !! vh_js !!
       IF (cable_user%soil_struc=='sli') THEN
          ! for stable conditions, update rough%rt0us & rough%rt1usa by replacing CCSW by
          ! csw = cd/2* (U(hc)/ust)**2 according to Eqs 15 & 19 from notes by Ian Harman (9-9-2011)
          WHERE (canopy%vlaiw > CLAI_thresh .AND. rough%hruff > rough%z0soilsn)
             rt0bus = (LOG(0.1*rough%hruff/rough%z0soilsn) - psis(canopy%zetash(:,iter)) + &
                  psis(canopy%zetash(:,iter)*rough%z0soilsn/(0.1*rough%hruff))) / &
                  Cvonk/rough%term6a

             zstar = rough%disp + 1.5*(veg%hc - rough%disp)

             psihat = LOG((zstar - rough%disp)/ (veg%hc - rough%disp)) + &
                  (veg%hc - zstar)/(zstar - rough%disp)
             rL = -(Cvonk*Cgrav*(zstar - rough%disp)*(canopy%fh))/ &  ! 1/Monin-Obokov Length
                  MAX( (air%rho*Ccapp*met%tk*canopy%us**3), 1.e-12)
             phist = 1 + 5.0*(zstar - rough%disp)*rL

             WHERE (canopy%zetar(:,iter) .GT. 1.e-6)! stable conditions
             csw = MIN(0.3*((LOG((veg%hc-rough%disp)/rough%z0m) + phist*psihat - &
             psim( canopy%zetar(:,iter)*(veg%hc-rough%disp) /                  &
                    (rough%zref_tq-rough%disp), mp, CPI_C )+ &
             psim(canopy%zetar(:,iter)*rough%z0m /                             &
                  (rough%zref_tq-rough%disp), mp, CPI_C ))                     &
                  / 0.4)**2/2., 3.0)* Ccsw

             ELSEWHERE
                csw = Ccsw
             endwhere

             rough%term2  = EXP( 2. * CSW * canopy%rghlai * &
                  ( 1 - rough%disp / rough%hruff ) )
             rough%term3  = CA33**2 * CCTL * 2. * CSW * canopy%rghlai
             rough%term5  = MAX( ( 2. / 3. ) * rough%hruff / rough%disp, 1.0 )
             rough%term6 =  EXP( 3. * rough%coexp * ( rough%disp / rough%hruff -1. ) )

             rough%rt0us  = LOG(rough%disp/(0.1 * rough%hruff)) * &
                  EXP(2. * CCSW * canopy%rghlai) * rough%disp &
                  / rough%hruff / (Ca33 ** 2 * Cctl)

             rough%rt1usa = rough%term5 * ( rough%term2 - 1.0 ) / rough%term3
             rt0 = MAX(rt_min, rough%rt0us+rt0bus) / canopy%us
          ELSEWHERE
             rt0 = MAX(rt_min,rough%rt0us) / canopy%us
          ENDWHERE

       ELSE ! NOT sli
          rt0 = MAX(rt_min,rough%rt0us / canopy%us)

          IF (cable_user%litter) THEN
             ! Mathews (2006), A process-based model of offine fuel moisture,
             !                 International Journal of Wildland Fire 15,155-168
             ! assuming here u=1.0 ms-1, bulk litter density 63.5 kgm-3
             canopy%kthLitt = 0.3_r_2 ! ~ 0.2992125984251969 = 0.2+0.14*0.045*1000.0/63.5
             canopy%DvLitt = 3.1415841138194147e-05_r_2 ! = 2.17e-5*exp(1.0*2.6)*exp(-0.5*(2.08+(1.0*2.38)))
          ENDIF

       ENDIF

!       ! Aerodynamic resistance (sum 3 height integrals)/us
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
       rough%rt1 = MAX(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)

       DO j=1,mp

          IF(canopy%vlaiw(j) > CLAI_THRESH) THEN
             ssnow%rtsoil(j) = rt0(j)
          ELSE
             ssnow%rtsoil(j) = rt0(j) + rough%rt1(j)
          ENDIF

       ENDDO


       ssnow%rtsoil = MAX(rt_min,ssnow%rtsoil)

       DO j=1,mp

          IF( ssnow%rtsoil(j) > 2.*ortsoil(j) .OR.                              &
               ssnow%rtsoil(j) < 0.5*ortsoil(j) ) THEN

             ssnow%rtsoil(j) = MAX(rt_min,0.5*(ssnow%rtsoil(j) + ortsoil(j)))

          ENDIF

       ENDDO


       IF (cable_user%or_evap) THEN
        write(6,*) "GW or ORevepis not an option right now"
        !H!          call or_soil_evap_resistance(soil,air,met,canopy,ssnow,veg,rough)
       END IF

       ! Vegetation boundary-layer conductance (mol/m2/s)
       ! Cprandt = kinematic viscosity/molecular diffusivity
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
       DO j=1,mp

          IF(canopy%vlaiw(j) > CLAI_THRESH) THEN
             gbvtop(j) = air%cmolar(j)*CAPOL * air%visc(j) / Cprandt /        &
                  veg%dleaf(j) * (canopy%us(j) / MAX(rough%usuh(j),1.e-6)&
                  * veg%dleaf(j) / air%visc(j) )**0.5                    &
                  * Cprandt**(1.0/3.0) / veg%shelrb(j)
             gbvtop(j) = MAX (0.05_r_2,gbvtop(j) )      ! for testing (BP aug2010)

             ! Forced convection boundary layer conductance
             ! (see Wang & Leuning 1998, AFM):


             !vh! inserted 'min' to avoid floating underflow
             gbhu(j,1) = gbvtop(j)*(1.0-EXP(-MIN(canopy%vlaiw(j)                    &
                  *(0.5*rough%coexp(j)+rad%extkb(j) ),20.0))) /            &
                  (rad%extkb(j)+0.5*rough%coexp(j))

             gbhu(j,2) = (2.0/rough%coexp(j))*gbvtop(j)*  &
                  (1.0-EXP(-MIN(0.5*rough%coexp(j)*canopy%vlaiw(j),20.0))) &
                  - gbhu(j,1)
          ENDIF

       ENDDO


       rny = SUM(rad%rniso,2) ! init current estimate net rad
       hcy = 0.0              ! init current estimate lat heat
       ecy = rny - hcy        ! init current estimate lat heat

       sum_rad_rniso = SUM(rad%rniso,2)

       CALL dryLeaf( dels, rad, rough, air, met,                                &
            veg, canopy, soil, ssnow, dsx,                             &
            fwsoil, tlfx, tlfy, ecy, hcy,                              &
            rny, gbhu, gbhf, csx, cansat,                              &
            ghwet,  iter,climate )


       CALL wetLeaf( dels, rad, rough, air, met,                                &
            veg, canopy, cansat, tlfy,                                 &
            gbhu, gbhf, ghwet )


       ! Calculate latent heat from vegetation:
       ! Calculate sensible heat from vegetation:
       ! Calculate net rad absorbed by canopy:
       canopy%fev = REAL(canopy%fevc + canopy%fevw)
       ftemp = (1.0 - canopy%fwet) *  REAL(hcy) + canopy%fhvw
       canopy%fhv = REAL(ftemp)
       ftemp= (1.0-canopy%fwet)*REAL(rny)+canopy%fevw+canopy%fhvw
       canopy%fnv = REAL(ftemp)

       ! canopy rad. temperature calc from long-wave rad. balance
       sum_rad_gradis = SUM(rad%gradis,2)

       DO j=1,mp

          IF ( canopy%vlaiw(j) > CLAI_THRESH .AND.                             &
               rough%hruff(j) > rough%z0soilsn(j) ) THEN

             rad%lwabv(j) = CCAPP * Crmair * ( tlfy(j) - met%tk(j) ) *        &
                  sum_rad_gradis(j)
             !! vh_js !!

             IF (  (rad%lwabv(j) / (2.0*(1.0-rad%transd(j))            &
                  * CSBOLTZ*CEMLEAF)+met%tvrad(j)**4) .GT. 0.0) THEN

                canopy%tv(j) = (rad%lwabv(j) / (2.0*(1.0-rad%transd(j))            &
                     * CSBOLTZ*CEMLEAF)+met%tvrad(j)**4)**0.25
             ELSE
                canopy%tv(j) = met%tvrad(j)
             ENDIF


          ELSE! sparse canopy

             canopy%tv(j) = met%tvrad(j)

          ENDIF

       ENDDO



       ! Calculate net rad to soil:
       canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*CEMLEAF* &
            CSBOLTZ*canopy%tv**4 - CEMSOIL*CSBOLTZ* tss4



       ! Saturation specific humidity at soil/snow surface temperature:
       CALL qsatfjh(ssnow%qstss,ssnow%tss-Ctfrz,met%pmb)

      if (cable_user%gw_model .OR.  cable_user%or_evap) & 
      write(6,*) "GW or ORevepis not an option right now"
      !H!        call pore_space_relative_humidity(ssnow,soil,veg)

       IF (cable_user%soil_struc=='default') THEN

          !REV_CORR - single location for calculating litter resistances
          !can go earlier in the code if needed
          IF (cable_user%litter) THEN
             rhlitt = REAL((1-ssnow%isflag))*veg%clitt*0.003/ &
                  canopy%kthLitt/(air%rho*CCAPP)
             relitt = REAL((1-ssnow%isflag))*veg%clitt*0.003/ &
                  canopy%DvLitt
          ENDIF

          !latent heat flux density of water (W/m2) - soil componenet of 
          !potential evapotranspiration 
          IF(cable_user%ssnow_POTEV== "P-M") THEN

             !--- uses %ga from previous timestep
             ssnow%potev = Penman_Monteith(canopy%ga)

          ELSE !by default assumes Humidity Deficit Method

             ! Humidity deficit
             ! INH: I think this should be - met%qvair
             dq = ssnow%qstss - met%qv
             dq_unsat = ssnow%rh_srf*ssnow%qstss - met%qv
             ssnow%potev =  Humidity_deficit_method(dq, dq_unsat,ssnow%qstss)

          ENDIF

          ! Soil latent heat:
          CALL latent_heat_flux()

          ! Calculate soil sensible heat:
          ! INH: I think this should be - met%tvair
          !canopy%fhs = air%rho*CCAPP*(ssnow%tss - met%tk) /ssnow%rtsoil
          IF (cable_user%gw_model .OR. cable_user%or_evap) THEN
             canopy%fhs =  air%rho*CCAPP*(ssnow%tss - met%tk) / &
                  (ssnow%rtsoil + ssnow%rt_qh_sublayer)
             !note if or_evap and litter are true then litter resistance is
             !incluyded above in ssnow%rt_qh_sublayer
          ELSEIF (cable_user%litter) THEN
             !! vh_js !! account for additional litter resistance to sensible heat transfer
             !! INH simplifying code using rhlitt
             canopy%fhs =  air%rho*CCAPP*(ssnow%tss - met%tk) / &
                                !(ssnow%rtsoil + real((1-ssnow%isflag))*veg%clitt*0.003/canopy%kthLitt/(air%rho*CCAPP))
                  (ssnow%rtsoil + rhlitt)
          ELSE
             canopy%fhs = air%rho*CCAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
          ENDIF

       ELSE

write(6,*) "SLI is not an option right now"
          ! SLI SEB to get canopy%fhs, canopy%fess, canopy%ga
          ! (Based on old Tsoil, new canopy%tv, new canopy%fns)
          !H!CALL sli_main(1,dels,veg,soil,ssnow,met,canopy,air,rad,1)

       ENDIF

       CALL within_canopy( gbhu, gbhf, rt0, rhlitt, relitt )

       ! Saturation specific humidity at soil/snow surface temperature:
       CALL qsatfjh(ssnow%qstss,ssnow%tss-Ctfrz,met%pmb)

       IF (cable_user%soil_struc=='default') THEN

          IF(cable_user%ssnow_POTEV== "P-M") THEN

             !--- uses %ga from previous timestep
             ssnow%potev = Penman_Monteith(canopy%ga)

          ELSE !by default assumes Humidity Deficit Method

             ! Humidity deficit
             dq = ssnow%qstss - met%qvair
             dq_unsat = ssnow%rh_srf*ssnow%qstss - met%qvair
             ssnow%potev =  Humidity_deficit_method(dq, dq_unsat,ssnow%qstss)

          ENDIF

          ! Soil latent heat:
          CALL latent_heat_flux()

          ! Soil sensible heat:
          !canopy%fhs = air%rho*CCAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
          IF (cable_user%gw_model .OR. cable_user%or_evap) THEN
             canopy%fhs =  air%rho*CCAPP*(ssnow%tss - met%tvair) / &
                  (ssnow%rtsoil + REAL(ssnow%rt_qh_sublayer))

          ELSEIF (cable_user%litter) THEN
             !! vh_js !! account for additional litter resistance to sensible heat transfer
             !! INH simplifying code using rhlitt
             canopy%fhs =  air%rho*CCAPP*(ssnow%tss - met%tvair) / &
                                !(ssnow%rtsoil +  real((1-ssnow%isflag))*veg%clitt*0.003/canopy%kthLitt/(air%rho*CCAPP))
                  (ssnow%rtsoil + rhlitt)
          ELSE
             canopy%fhs = air%rho*CCAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
             

          ENDIF

          !! Ticket #90 ssnow%cls factor should be retained: required for energy balance
          !! INH: %cls factor included in %fes already - do not include here
          canopy%ga = canopy%fns-canopy%fhs-canopy%fes !*ssnow%cls
       ELSE

write(6,*) "SLI is not an option right now"
          ! SLI SEB to get canopy%fhs, canopy%fess, canopy%ga
          ! (Based on old Tsoil, new canopy%tv, new canopy%fns)
          !H!CALL sli_main(1,dels,veg,soil,ssnow,met,canopy,air,rad,1)

       ENDIF

       ! Set total latent heat:
       canopy%fe = canopy%fev + canopy%fes

       ! Set total sensible heat:
       canopy%fh = canopy%fhv + canopy%fhs

       !---diagnostic purposes
       DO j=1,mp

          IF (ssnow%potev(j) .GE. 0.) THEN
             ssnow%potev(j) = MAX(0.00001,ssnow%potev(j))
          ELSE
             ssnow%potev(j) = MIN(-0.0002,ssnow%potev(j))
          ENDIF

          IF (canopy%fevw_pot(j) .GE. 0.) THEN
             canopy%fevw_pot(j) = MAX(0.000001,canopy%fevw_pot(j))
          ELSE
             canopy%fevw_pot(j) = MIN(-0.002,canopy%fevw_pot(j))
          ENDIF

       ENDDO


       canopy%rnet = canopy%fnv + canopy%fns

       !upwards flux density of water (kg/m2/s) - potential evapotranspiration 
       !radiation weighted soil and canopy contributions
       !Note: If PM routine corrected then match changes here
       canopy%epot = ((1.-rad%transd)*canopy%fevw_pot +                         &
            rad%transd*ssnow%potev*ssnow%cls) * dels/air%rlam



       canopy%rniso = SUM(rad%rniso,2) + rad%qssabs + rad%transd*met%fld + &
            (1.0-rad%transd)*CEMLEAF* &
            CSBOLTZ*met%tvrad**4 - CEMSOIL*CSBOLTZ*met%tvrad**4

       rlower_limit = canopy%epot * air%rlam / dels
       WHERE (rlower_limit == 0 ) rlower_limit = 1.e-7 !prevent from 0. by adding 1.e-7 (W/m2)


       canopy%wetfac_cs = MAX(0., MIN(1.0,canopy%fe / rlower_limit ))

       DO j=1,mp

          IF ( canopy%wetfac_cs(j) .LE. 0. )                                    &
               canopy%wetfac_cs(j) = MAX( 0., MIN( 1.,                            &
               MAX( canopy%fev(j)/canopy%fevw_pot(j),       &
               REAL(canopy%fes(j))/ssnow%potev(j) ) ) )

       ENDDO

       CALL update_zetar()
       !!880!CALL update_zetar(mp, sunlit_veg_mask)

    END DO           ! do iter = 1, NITER



    canopy%cduv = canopy%us * canopy%us / (MAX(met%ua,CUMIN))**2

    !---diagnostic purposes
    canopy%gswx_T = rad%fvlai(:,1)/MAX( CLAI_THRESH, canopy%vlaiw(:) )         &
         * canopy%gswx(:,1) + rad%fvlai(:,2) / MAX(CLAI_THRESH,     &
         canopy%vlaiw(:))*canopy%gswx(:,2)

    ! The surface conductance below is required by dust scheme; it is composed from canopy and soil conductances
    canopy%gswx_T = (1.-rad%transd)*MAX(1.e-06,canopy%gswx_T ) +  &   !contribution from  canopy conductance
         rad%transd*(.01*ssnow%wb(:,1)/soil%sfc)**2 ! + soil conductance; this part is done as in Moses
    WHERE ( soil%isoilm == 9 ) canopy%gswx_T = 1.e6   ! this is a value taken from Moses for ice points

    canopy%cdtq = canopy%cduv * &
                    ( LOG( rough%zref_uv / rough%z0m) -              &
                      psim( canopy%zetar(:,NITER) * rough%zref_uv/rough%zref_tq, mp, CPI_C )   &
                     + psim( canopy%zetar(:,NITER) * rough%z0m    /rough%zref_tq, mp, CPI_C ) & ! new term from Ian Harman
                    )                                          &
                    / ( LOG( rough%zref_tq /(0.1*rough%z0m) )                   &
                      - psis( canopy%zetar(:,NITER) )                                  &
                      + psis(canopy%zetar(:,NITER)*0.1*rough%z0m/rough%zref_tq ) ) ! n

    !INH - the screen level calculations should be split off into a new subroutine -------

    ! Calculate screen temperature: 1) original method from SCAM
    ! screen temp., windspeed and relative humidity at 1.5m
    ! screen temp., windspeed and relative humidity at 2.0m
    ! cls factor included in qstar
    tstar = - canopy%fh / ( air%rho*CCAPP*canopy%us)
    qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us * ssnow%cls)
    zscrn = MAX(rough%z0m,2.0-rough%disp)
    ftemp = ( LOG(rough%zref_tq/zscrn)- psis(canopy%zetar(:,iterplus)) +       &
         psis(canopy%zetar(:,iterplus) * zscrn / rough%zref_tq) ) /CVONK

    ! Calculate screen temperature:
    canopy%tscrn = met%tk - Ctfrz - tstar * ftemp

    ! Calculate radiative/skin temperature;
    ! at this stage old soil temperature is used
    ! calculation of screen temepratures for LAI > 0.1 . Method by Ian Harman

    term1=0.
    term2=0.
    term5=0.
    term3 = 0. ! Work around for Intel compiler problem with nested whres
    r_sc = 0.  ! sum of resistance from ground to screen level
    zscl = MAX(rough%z0soilsn,2.0)

    ! assume screen temp of bareground if all these conditions are not met
    DO j=1,mp

       IF ( canopy%vlaiw(j) > CLAI_THRESH .AND. rough%hruff(j) > 0.01) THEN

          IF ( rough%disp(j)  > 0.0 ) THEN

             term1(j) = EXP(2*CCSW*canopy%rghlai(j)*(1-zscl(j)/rough%hruff(j)))
             term2(j) = EXP(2*CCSW*canopy%rghlai(j) *                          &
                  (1-rough%disp(j)/rough%hruff(j)))
             term5(j) = MAX(2./3.*rough%hruff(j)/rough%disp(j), 1.)

          ENDIF

          term3(j) = CA33**2*CCTL*2*CCSW*canopy%rghlai(j)

          IF( zscl(j) < rough%disp(j) ) THEN

             !Ticket #154
             !r_sc(j) = term5(j) * LOG(zscl(j)/rough%z0soilsn(j)) *              &
             !     ( EXP(2*CCSW*canopy%rghlai(j)) - term1(j) ) / term3(j)
             r_sc(j) = term5(j) * LOG(zscl(j)/rough%z0soilsn(j)) *              &
                  ( EXP(2*CCSW*canopy%rghlai(j)) - term2(j) ) / term3(j)
             r_sc(j) = r_sc(j) + term5(j) * LOG(rough%disp(j)/rough%z0soilsn(j)) *  &
                  ( EXP(2*CCSW*canopy%rghlai(j)) - term1(j) ) / term3(j)

          ELSEIF( rough%disp(j) <= zscl(j) .AND.                                &
               zscl(j) < rough%hruff(j) ) THEN

             r_sc(j) = rough%rt0us(j) + term5(j) * ( term2(j) - term1(j) ) /    &
                  term3(j)

          ELSEIF( rough%hruff(j) <= zscl(j) .AND.                               &
               zscl(j) <  rough%zruffs(j) ) THEN

             r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + term5(j) *            &
                  ( zscl(j) - rough%hruff(j) ) /                           &
                  ( CA33**2 * CCTL * rough%hruff(j) )


          ELSEIF( zscl(j) >= rough%zruffs(j) ) THEN
             !Ticket #67 - Modify order of operations to avoid potential error
             r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j) +     &
                  ( LOG( (zscl(j) - rough%disp(j)) /                       &
                  MAX( rough%zruffs(j)-rough%disp(j),                      &
                  rough%z0soilsn(j) ) ) - psis( (zscl(j)-rough%disp(j))    &
                                !Ticket #67 - change order of operations to avoid /0
                                !        / (rough%zref_tq(j)/canopy%zetar(j,iterplus) ) )        &
                  * canopy%zetar(j,iterplus)/rough%zref_tq(j) )            &
                  + psis( (rough%zruffs(j) - rough%disp(j) )               &
                                !        / (rough%zref_tq(j)/canopy%zetar(j,iterplus ) ) ) )     &
                  * canopy%zetar(j,iterplus)/rough%zref_tq(j) ) )          &
                  / CVONK

          ENDIF

          !extensions for litter and Or evaporation model
          IF (cable_user%litter) THEN
             canopy%tscrn(j) = ssnow%tss(j) + (met%tk(j) - ssnow%tss(j)) *     &
                  MIN(1., ( (r_sc(j)+rhlitt(j)*canopy%us(j))  / MAX( 1.,          &
                  rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)              &
                  + rt1usc(j) + rhlitt(j)*canopy%us(j) )) ) - Ctfrz
          ELSEIF (cable_user%or_evap .OR. cable_user%gw_model) THEN
             canopy%tscrn(j) = ssnow%tss(j) + (met%tk(j) - ssnow%tss(j)) *     &
                  MIN(1., ( (ssnow%rt_qh_sublayer(j)*canopy%us(j) + r_sc(j) ) /   &
                  MAX( 1., rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)     &
                  + rt1usc(j) + ssnow%rt_qh_sublayer(j)*canopy%us(j) )) ) - Ctfrz
          ELSE
             canopy%tscrn(j) = ssnow%tss(j) + (met%tk(j) - ssnow%tss(j)) *      &
                  MIN(1., (r_sc(j) / MAX( 1.,                            &
                  rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)   &
                  + rt1usc(j))) )  - Ctfrz
          ENDIF

       ENDIF

    ENDDO


    !screen level humdity - this is only approximate --------------------------
    CALL qsatfjh(rsts,canopy%tscrn,met%pmb)

    qtgnet = rsts * ssnow%wetfac - met%qv

    DO j=1,mp

       IF (qtgnet(j) .GT. 0. ) THEN
          qsurf(j) = rsts(j) * ssnow%wetfac(j)
       ELSE
          qsurf(j) = 0.1*rsts(j)*ssnow%wetfac(j) + 0.9*met%qv(j)
       ENDIF

       canopy%qscrn(j) = met%qv(j) - qstar(j) * ftemp(j)

       IF( canopy%vlaiw(j) >CLAI_THRESH .AND. rough%hruff(j) > 0.01) THEN
          IF (cable_user%litter) THEN

          !extensions for litter and Or model
             canopy%qscrn(j) = qsurf(j) + (met%qv(j) - qsurf(j)) *             &
                  MIN(1., ( ( r_sc(j)+relitt(j)*canopy%us(j) ) / MAX( 1.,         &
                  rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)            &
                  + rt1usc(j) + relitt(j)*canopy%us(j) )) )

          ELSEIF (cable_user%or_evap .OR. cable_user%gw_model) THEN
             !using alpm1 as a dumy variable
             alpm1(j) = REAL(&
                  ssnow%satfrac(j)/(REAL(ssnow%rtsoil(j),r_2)+&
                  ssnow%rtevap_sat(j)) &
                  + (1.0-ssnow%satfrac(j))/(REAL(ssnow%rtsoil(j),r_2)+ ssnow%rtevap_unsat(j)) &
                  )

             canopy%qscrn(j) = qsurf(j) + (met%qv(j) - qsurf(j)) *             &
                  MIN(1., ( (r_sc(j) + canopy%us(j)/alpm1(j) ) / MAX( 1.,         &
                  rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)              &
                  + rt1usc(j) + canopy%us(j)/alpm1(j) )) )

          ELSE
             canopy%qscrn(j) = qsurf(j) + (met%qv(j) - qsurf(j)) *             &
                  MIN(1., (r_sc(j) / MAX( 1.,                                     &
                  rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)              &
                  + rt1usc(j))) )
          ENDIF

       ENDIF

    ENDDO

    ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
    canopy%dewmm = - (MIN(0.0,canopy%fevw) + MIN(0.0_r_2,canopy%fevc)) * &
         dels * 1.0e3 / (CRHOW*air%rlam)

    ! Add dewfall to canopy water storage:
    canopy%cansto = canopy%cansto + canopy%dewmm

    ! Modify canopy water storage for evaporation:
    canopy%cansto = MAX(canopy%cansto-MAX(0.0,REAL(canopy%fevw))*dels &
         *1.0e3/(CRHOW*air%rlam), 0.0)

    ! Calculate canopy water storage excess:
    canopy%spill=MAX(0.0, canopy%cansto-cansat)

    ! Move excess canopy water to throughfall:
    ! %through is /dels in UM app. (unpacked in hyd driver) for STASH output
    canopy%through = canopy%through + canopy%spill

    ! Initialise 'throughfall to soil' as 'throughfall from canopy';
    ! snow may absorb
    canopy%precis = MAX(0.,canopy%through)

    ! Update canopy storage term:
    canopy%cansto=canopy%cansto - canopy%spill

    ! Calculate the total change in canopy water store (mm/dels):
    canopy%delwc = canopy%cansto-canopy%oldcansto

    ! calculate dgdtg, derivative of ghflux 3 instances
    ! d(canopy%fns)/d(ssnow%tgg)
    ! d(canopy%fhs)/d(ssnow%tgg)
    ! d(canopy%fes)/d(dq)
    !IF (cable_user%soil_struc=='default') THEN
       ssnow%dfn_dtg = (-1.)*4.*CEMSOIL*CSBOLTZ*tss4/ssnow%tss

       !INH: REV_CORR revised sensitivity terms working variable
       rttsoil = ssnow%rtsoil
       IF (cable_user%L_REV_CORR) THEN
          WHERE (canopy%vlaiw > CLAI_THRESH)
            !if %vlaiw<=%LAI_THRESH then %rt1 already added to %rtsoil
            rttsoil = rttsoil + rough%rt1
          ENDWHERE
       ENDIF

IF (cable_user%gw_model .or. cable_user%or_evap) THEN

  ssnow%dfh_dtg = air%rho*CCAPP/(ssnow%rtsoil+ real(ssnow%rt_qh_sublayer))
  
  !! INH simplifying code for legibility
  !ssnow%dfe_ddq = real(ssnow%satfrac)*air%rho*air%rlam*ssnow%cls/ &
  !     (ssnow%rtsoil+ real(ssnow%rtevap_sat))  +
  !     (1.0-real(ssnow%satfrac))*real(ssnow%rh_srf)*&
  !      air%rho*air%rlam*ssnow%cls/ (ssnow%rtsoil+
  !      real(ssnow%rtevap_unsat) )
   ssnow%dfe_ddq = real(ssnow%satfrac)/(ssnow%rtsoil+ real(ssnow%rtevap_sat))  &
              + (1.0-real(ssnow%satfrac))*real(ssnow%rh_srf)                   &
                   / (ssnow%rtsoil+ real(ssnow%rtevap_unsat) )

  !mrd561 fixes.  Do same thing as INH but has been tested.
  IF (cable_user%L_REV_CORR) THEN
    alpm1  = real(ssnow%satfrac/(real(ssnow%rtsoil,r_2)+ ssnow%rtevap_sat) +     &
                  (1.0-ssnow%satfrac) / (real(ssnow%rtsoil,r_2)+ ssnow%rtevap_unsat ) )
    beta2 = real(ssnow%satfrac/(real(ssnow%rtsoil,r_2)+ ssnow%rtevap_sat) +     &
                 (1.0-ssnow%satfrac) * ssnow%rh_srf                  &
                 / (real(ssnow%rtsoil,r_2)+ ssnow%rtevap_unsat ) )
    WHERE (canopy%vlaiw > CLAI_THRESH)
      alpm1 = alpm1 + 1._r_2/real(rough%rt1,r_2)
      beta_div_alpm  = beta2 / alpm1  !might need limit here
      rttsoil = ssnow%rtsoil + rough%rt1
    ELSEWHERE!if there is no canopy then qa should not change
      beta_div_alpm=0.0  !do not divide by aplm1 prevent issues
      rttsoil = ssnow%rtsoil 
    ENDWHERE
    ssnow%dfh_dtg = air%rho*CCAPP/(rttsoil +               & 
                              real(ssnow%rt_qh_sublayer))
    ssnow%dfe_ddq = real(ssnow%satfrac*(1.0-real(beta_div_alpm,r_2)) /        & 
                    (real(ssnow%rtsoil,r_2)+ ssnow%rtevap_sat) +           &
                    (1.0-ssnow%satfrac)* (ssnow%rh_srf - real(beta_div_alpm,r_2)) /    &
                    (real(ssnow%rtsoil,r_2)+ ssnow%rtevap_unsat ) )

  ELSE ! IF (cable_user%L_REV_CORR) THEN
    ssnow%dfh_dtg = air%rho*CCAPP/(ssnow%rtsoil+ real(ssnow%rt_qh_sublayer))
    ssnow%dfe_ddq = real(ssnow%satfrac)/(ssnow%rtsoil+ real(ssnow%rtevap_sat))  &
                         + (1.0-real(ssnow%satfrac))*real(ssnow%rh_srf)                   &
                              / (ssnow%rtsoil+ real(ssnow%rtevap_unsat) )
  ENDIF ! IF (cable_user%L_REV_CORR) THEN
                
  !cls applies for both REV_CORR false and true          
  ssnow%dfe_ddq = ssnow%dfe_ddq*air%rho*air%rlam*ssnow%cls
  
  !REV_CORR: factor %wetfac needed for potev>0. and gw_model &/or snow cover
  !NB %wetfac=1. if or_evap
  IF (cable_user%L_REV_CORR) THEN
    WHERE (ssnow%potev >= 0.)
      ssnow%dfe_ddq = ssnow%dfe_ddq*ssnow%wetfac
    ENDWHERE       
  ENDIF


ELSEIF (cable_user%litter) THEN ! IF (cable_user%gw_model .or. cable_user%or_evap) THEN
  !!vh_js!! INH simplifying code for legibility and REV_CORR
  !ssnow%dfh_dtg = air%rho*CCAPP/(ssnow%rtsoil+ &
  !     real((1-ssnow%isflag))*veg%clitt*0.003/canopy%kthLitt/(air%rho*CCAPP))
  !ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/ &
  !     (ssnow%rtsoil+ real((1-ssnow%isflag))*veg%clitt*0.003/canopy%DvLitt)
  
  !recalculated - probably not needed 
  rhlitt = real((1-ssnow%isflag))*veg%clitt*0.003/canopy%kthLitt/(air%rho*CCAPP)
  relitt = real((1-ssnow%isflag))*veg%clitt*0.003/canopy%DvLitt
  
  !incorporates REV_CORR changes
  ssnow%dfh_dtg = air%rho*CCAPP/(rttsoil+rhlitt)
  ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/(rttsoil+relitt)

  !REV_CORR: factor ssnow%wetfac is not applied if dew/frost i.e. potev<0
  IF (cable_user%L_REV_CORR) THEN
     WHERE (ssnow%potev < 0.)
         ssnow%dfe_ddq = air%rho*air%rlam*ssnow%cls/(rttsoil+relitt)
     ENDWHERE       
  ENDIF

ELSE ! i.e. NOT (%gw_model .or. %or_evap or SLI)
  !ssnow%dfh_dtg = air%rho*CCAPP/ssnow%rtsoil
  !ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/ssnow%rtsoil
  
  !incorporates REV_CORR changes
  ssnow%dfh_dtg = air%rho*CCAPP/rttsoil
  ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/rttsoil
  
  !REV_CORR: factor ssnow%wetfac is not applied if dew/frost i.e. potev<0
  IF (cable_user%L_REV_CORR) THEN
     WHERE (ssnow%potev < 0.)
         ssnow%dfe_ddq = air%rho*air%rlam*ssnow%cls/(rttsoil+relitt)
     ENDWHERE       
  ENDIF      

ENDIF ! IF (cable_user%gw_model .or. cable_user%or_evap) THEN

ssnow%ddq_dtg = (Crmh2o/Crmair) /met%pmb * CTETENA*CTETENB * CTETENC   &
     / ( ( CTETENC + ssnow%tss-Ctfrz )**2 )*EXP( CTETENB *       &
     ( ssnow%tss-Ctfrz ) / ( CTETENC + ssnow%tss-Ctfrz ) )

!canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%dfe_ddq *    &
!     ssnow%ddq_dtg

!INH: REV_CORR Rewritten for flexibility
ssnow%dfe_dtg = ssnow%dfe_ddq * ssnow%ddq_dtg
canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%dfe_dtg

bal%drybal = REAL(ecy+hcy) - SUM(rad%rniso,2)                               &
     + CCAPP*Crmair*(tlfy-met%tk)*SUM(rad%gradis,2)  ! YP nov2009

bal%wetbal = canopy%fevw + canopy%fhvw - SUM(rad%rniso,2) * canopy%fwet      &
     + CCAPP*Crmair * (tlfy-met%tk) * SUM(rad%gradis,2) *          &
     canopy%fwet  ! YP nov2009

DEALLOCATE(cansat,gbhu)
DEALLOCATE(dsx, fwsoil, tlfx, tlfy)
DEALLOCATE(ecy, hcy, rny)
DEALLOCATE(gbhf, csx)
DEALLOCATE(ghwet)

RETURN

  CONTAINS

    FUNCTION Penman_Monteith( ground_H_flux ) RESULT(ssnowpotev)
      USE cable_def_types_mod, ONLY : mp
      REAL, INTENT(IN), DIMENSION(mp)  :: ground_H_flux
      REAL, DIMENSION(MP)  ::                                                     &
           ssnowpotev,      & ! returned result of function
           sss,             & ! var for Penman-Monteith soil evap
           cc1,             & ! var for Penman-Monteith soil evap
           cc2,             & ! var for Penman-Monteith soil evap
           qsatfvar           !
      INTEGER :: j

      ! Penman-Monteith formula
      sss=air%dsatdk
      cc1=sss/(sss+air%psyc )
      cc2=air%psyc /(sss+air%psyc )

      CALL qsatfjh(qsatfvar,met%tvair-Ctfrz,met%pmb)

      !INH 10-1-2017 - this P-M implementation is incorrect over snow.
      !variable ssnowpotev is actually the latent heat flux associated with
      !potential evaporation.
      !Needs to be addressed/simplified at a later date - involves changes
      !to HDM method and latent_heat_flux() and elsewhere

      IF (cable_user%litter) THEN
         !! vh_js !!
         ssnowpotev = cc1 * (canopy%fns - ground_H_flux) + &
              cc2 * air%rho * air%rlam*(qsatfvar - met%qvair)/ &
              (ssnow%rtsoil+ REAL((1-ssnow%isflag))*veg%clitt*0.003/canopy%DvLitt)
      ELSE
         ssnowpotev = cc1 * (canopy%fns - ground_H_flux) + &
              cc2 * air%rho * air%rlam*(qsatfvar  - met%qvair)/ssnow%rtsoil
      ENDIF

    END FUNCTION Penman_Monteith


    ! ------------------------------------------------------------------------------
    ! method alternative to P-M formula above
    FUNCTION humidity_deficit_method(dq,dqu,qstss ) RESULT(ssnowpotev)

      USE cable_def_types_mod, ONLY : mp

      REAL, DIMENSION(mp) ::                                                      &
           ssnowpotev,    & !
           dq,            & ! sat spec hum diff.
           dqu,           & ! sat spec hum diff.
           qstss             !dummy var for compilation

      INTEGER :: j
      REAL, DIMENSION(mp) :: q_air

      q_air = qstss - dq

      DO j=1,mp
         !if(ssnow%snowd(j) > 1.0) dq(j) = max( -0.1e-3, dq(j))
         IF( ssnow%snowd(j)>1.0 .OR. ssnow%tgg(j,1).EQ.Ctfrz)   THEN
            dq(j) = MAX( -0.1e-3, dq(j))
            dqu(j) = MAX( -0.1e-3, dqu(j))
         END IF

         IF (dq(j) .LE. 0.0 .AND. dqu(j) .LT. dq(j)) THEN
            dqu(j) = dq(j)
         END IF

         IF (dq(j) .GE. 0.0 .AND. dqu(j) .LT. 0.0) THEN
            dqu(j) = 0.0
         ENDIF
      ENDDO

      IF (cable_user%or_evap .or. cable_user%gw_model) then

write(6,*) "GW or ORevepis not an option right now"
!H!        IF (cable_user%or_evap) THEN
!H!          do j=1,mp
!H!       
!H!             if (veg%iveg(j) .lt. 16 .and. ssnow%snowd(j) .lt. 1e-7) THEN
!H!       
!H!                if (dq(j) .le. 0.0) THEN
!H!                   ssnow%rtevap_sat(j) = min(rtevap_max,canopy%sublayer_dz(j)/rt_Dff)
!H!                end if
!H!       
!H!                if (dqu(j) .le. 0.0) THEN
!H!                   ssnow%rtevap_unsat(j) = min(rtevap_max,canopy%sublayer_dz(j)/rt_Dff)
!H!                end if
!H!       
!H!             end if
!H!
!H!          end do
!H!
!H!        END IF

         ssnowpotev = air%rho * air%rlam * ( &
              REAL(ssnow%satfrac) * dq /(ssnow%rtsoil + REAL(ssnow%rtevap_sat)) + &
              (1.0 - REAL(ssnow%satfrac))* dqu/( &
              ssnow%rtsoil + REAL(ssnow%rtevap_unsat)) )

      ELSEIF (cable_user%litter) THEN
         !! vh_js !!
         ssnowpotev =air%rho * air%rlam * dq /(ssnow%rtsoil + &
              REAL((1-ssnow%isflag))*veg%clitt*0.003/canopy%DvLitt)
      ELSE
         ssnowpotev =air%rho * air%rlam * dq /ssnow%rtsoil
      ENDIF

    END FUNCTION Humidity_deficit_method

    ! ------------------------------------------------------------------------------

    SUBROUTINE Latent_heat_flux()

      USE cable_common_module
      USE cable_def_types_mod, ONLY : mp

      REAL, DIMENSION(mp) ::                                                      &
           frescale,  flower_limit, fupper_limit

      INTEGER :: j

      !Ticket 137 - adjustments made so that one of four cases occurs
      !             i) evaporation from/dew onto surfaces with no snow, T>Tfrz
      !            ii) sublimation from/frost onto surfaces with snow cover
      !           iii) evaporation of liquid water if no snow but frozen soil
      !            iv) deposition of frost onto frozen soils if no snow cover
      !
      !IMPORTANTLY the value of %cls set here is used to control whether
      !water fluxes are from the snow pack or soil column in _soilsnow

      ! Soil latent heat:
      WHERE (ssnow%potev < 0. ) ssnow%wetfac(:) = 1.0
      canopy%fess= ssnow%wetfac * ssnow%potev
      !removed below, set wetfac to 1.0 is potev < 0.0
      !WHERE (ssnow%potev < 0. ) canopy%fess = ssnow%potev

      ! Reduce soil evap due to presence of puddle
      pwet = MAX(0.,MIN(0.2,ssnow%pudsto/MAX(1.,ssnow%pudsmx)))
      canopy%fess = canopy%fess * (1.-pwet)

      frescale = soil%zse(1) * 1000. * air%rlam / dels

      DO j=1,mp

         IF(ssnow%snowd(j) < 0.1 .AND. canopy%fess(j) .GT. 0. ) THEN

            IF (.NOT.cable_user%l_new_reduce_soilevp) THEN
               flower_limit(j) = REAL(ssnow%wb(j,1))-soil%swilt(j)/2.0
            ELSE
               ! E.Kowalczyk 2014 - reduces the soil evaporation
               flower_limit(j) = REAL(ssnow%wb(j,1))-soil%swilt(j)
            ENDIF
            fupper_limit(j) = MAX( 0.,                                        &
                 flower_limit(j) * frescale(j)                       &
                 - ssnow%evapfbl(j,1)*air%rlam(j)/dels)

            canopy%fess(j) = MIN(canopy%fess(j), REAL(fupper_limit(j),r_2))

            !Ticket 137 - case iii)
            !evaporation from frozen soils needs to respect the assumption that
            !ice fraction of soil moisture cannot exceed frozen_limit=0.85
            !see soilsnow: if frozen_limit changes need to be consistent
            fupper_limit(j) = REAL(ssnow%wb(j,1)-ssnow%wbice(j,1)/0.85)*frescale(j)
            fupper_limit(j) = MAX(REAL(fupper_limit(j),r_2),0.)

            canopy%fess(j) = MIN(canopy%fess(j), REAL(fupper_limit(j),r_2))

         END IF

         ssnow%cls(j)=1.

         !Ticket 137 - case ii) deposition of frost onto snow
         ! case of sublimation of snow overwrites later
         IF (ssnow%snowd(j) >=0.1 ) THEN
            ssnow%cls(j) = 1.1335
            canopy%fess(j) = ssnow%cls(j)*ssnow%potev(j)
         ENDIF

         !Ticket 137 - case iv) deposition of frost onto frozen soil, no snow
         IF (ssnow%snowd(j) < 0.1 .AND. ssnow%potev(j) < 0. .AND. &
              ssnow%tss(j)<CTFRZ) THEN
            ssnow%cls(j)=1.1335
            canopy%fess(j) = ssnow%cls(j)*ssnow%potev(j)
         ENDIF

         !Ticket 137 - case ii) sublimation of snow
         IF (ssnow%snowd(j) >= 0.1 .AND. ssnow%potev(j) > 0.) THEN

            ssnow%cls(j) = 1.1335

            !INH - if changes to PM routine then matching changes here
            canopy%fess(j) = MIN( (ssnow%wetfac(j)*ssnow%potev(j))*ssnow%cls(j), &
                 ssnow%snowd(j)/dels*air%rlam(j)*ssnow%cls(j))

         ENDIF

      ENDDO

      ! Evaporation from soil puddle
      canopy%fesp = MIN(ssnow%pudsto/dels*air%rlam,MAX(pwet*ssnow%potev,0.))
      canopy%fes = canopy%fess + canopy%fesp

    END SUBROUTINE latent_heat_flux

    ! -----------------------------------------------------------------------------

    SUBROUTINE within_canopy( gbhu, gbhf, rt0, rhlitt, relitt )

      USE cable_def_types_mod, ONLY : mp, r_2

      REAL(r_2), INTENT(IN), DIMENSION(:,:) ::                                    &
           gbhu,    &  ! forcedConvectionBndryLayerCond
           gbhf        ! freeConvectionBndryLayerCond

      REAL, INTENT(IN), DIMENSION(mp) ::                                       &
           rt0,     &  ! resistance from ground to canopy air
           rhlitt,  &  ! additional litter resistance for heat (=0 by default)
           relitt      ! additional litter resistance for water

      REAL, DIMENSION(mp) ::                                                      &
           rrsw,             & ! recipr. stomatal resistance for water
           rrbw,             & ! recipr. leaf boundary layer resistance for water
           dmah,             & ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmbh,             & ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmch,             & ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmae,             & ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmbe,             & ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
           dmce                ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132

      REAL  :: lower_limit, upper_limit
      REAL, DIMENSION(mp) :: fix_eqn,fix_eqn2

      INTEGER :: j

      !INH: rhlitt=relitt=0. if litter resistance not active but case included
      !dmah through to dmce are not A_{H} through C_{E} as per Eqn 3.40
      !in SCAM documentation but rt0*((1+esp)/rs + 1/rb)*A_{H} etc.
      !
      !changes from v1.4 for %cls package, litter and Or hydrology

      rrbw = SUM(gbhu+gbhf,2)/air%cmolar  ! MJT

      ! leaf stomatal resistance for water
      rrsw = SUM(canopy%gswx,2)/air%cmolar ! MJT

      IF (cable_user%or_evap) THEN
         fix_eqn(:) = rt0(:)*(REAL(ssnow%satfrac(:))/(rt0(:)+REAL(ssnow%rtevap_sat(:))) + &
              (1-REAL(ssnow%satfrac(:)))/(rt0(:)+REAL(ssnow%rtevap_unsat(:))))
         !lakes/ice rtevap=0 and wetfac is .ne. 1
         fix_eqn(:) = ssnow%wetfac(:) * fix_eqn(:)*ssnow%cls(:)   !INH correction. & M.Dekker +d wetfac

         fix_eqn2(:) = rt0(:) / (rt0(:) + REAL(ssnow%rt_qh_sublayer) )

      ELSE  !with INH corrections for litter and cls

         fix_eqn(:) = ssnow%cls(:)*rt0(:)/(rt0(:)+relitt(:))
         WHERE (ssnow%potev>0.) fix_eqn(:)=fix_eqn(:)*ssnow%wetfac(:)
         fix_eqn2(:) = rt0(:)/(rt0(:)+rhlitt(:))

      END IF

      DO j=1,mp

         IF(veg%meth(j) > 0 .AND. canopy%vlaiw(j) > CLAI_THRESH .AND.              &
              rough%hruff(j) > rough%z0soilsn(j) ) THEN

            !   use the dispersion matrix (DM) to find the air temperature
            !   and specific humidity
            !   (Raupach, Finkele and Zhang 1997, pp 17)
            ! leaf boundary layer resistance for water
            ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmah(j) = (rt0(j)+fix_eqn2(j)*rough%rt1(j))*((1.+air%epsi(j))*rrsw(j) + rrbw(j))  &
                 + air%epsi(j) * (rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))

            ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmbh(j) = (-air%rlam(j)/CCAPP)*(rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))

            ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmch(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)*   &
                 (canopy%fhv(j) + canopy%fhs(j))/(air%rho(j)*CCAPP)

            ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmae(j) = (-air%epsi(j)*CCAPP/air%rlam(j))*(rt0(j)*rough%rt1(j)) *   &
                 (rrbw(j)*rrsw(j))

            ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            dmbe(j) = ( rt0(j) + fix_eqn(j) * rough%rt1(j) ) *               &
                 ( (1.+air%epsi(j) ) * rrsw(j) + rrbw(j) ) +                 &
                 ( rt0(j) * rough%rt1(j) ) * ( rrbw(j) * rrsw(j) )

            ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
            ! INH: includes modifications for %cls
            dmce(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)*   &
                 (canopy%fev(j) + canopy%fes(j)/ssnow%cls(j)) /                   &
                 (air%rho(j)*air%rlam(j))

            ! Within canopy air temperature:
            met%tvair(j) = met%tk(j) + ( dmbe(j) * dmch(j) - dmbh(j) * dmce(j) )  &
                 / (dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)

            !---set limits for comparisson
            lower_limit =  MIN( ssnow%tss(j), met%tk(j)) - 5.0
            upper_limit =  MAX( ssnow%tss(j), met%tk(j)) + 5.0

            !--- tvair within these limits
            met%tvair(j) = MAX(met%tvair(j) , lower_limit)
            met%tvair(j) = MIN(met%tvair(j) , upper_limit)

            ! recalculate using canopy within temperature
            met%qvair(j) = met%qv(j) + (dmah(j)*dmce(j)-dmae(j)*dmch(j)) /        &
                 ( dmah(j)*dmbe(j)-dmae(j)*dmbh(j)+1.0e-12)
            met%qvair(j) = MAX(0.0,met%qvair(j))

            !---set limits for comparisson

            lower_limit =  MIN( ssnow%qstss(j), met%qv(j))
            upper_limit =  MAX( ssnow%qstss(j), met%qv(j))

            !--- qvair within these limits
            met%qvair(j) =  MAX(met%qvair(j),lower_limit)
            met%qvair(j) =  MIN(met%qvair(j), upper_limit)

            ! Saturated specific humidity in canopy:
            CALL qsatfjh2(qstvair(j),met%tvair(j)-Ctfrz,met%pmb(j))

            ! Saturated vapour pressure deficit in canopy:
            met%dva(j) = ( qstvair(j) - met%qvair(j) ) *  Crmair/CRMH2o         &
                 * met%pmb(j) * 100.
         ENDIF

      ENDDO

    END SUBROUTINE within_canopy

    ! -----------------------------------------------------------------------------

    SUBROUTINE update_zetar()

      INTEGER :: j

      ! monin-obukhov stability parameter zetar=zref/l
      ! recompute zetar for the next iteration, except on last iteration
      IF (iter < NITER) THEN ! dont compute zetar on the last iter

         iterplus = MAX(iter+1,2)
         canopy%zetar(:,iterplus) = -( CVONK * CGRAV * rough%zref_tq *              &
              ( canopy%fh + 0.07 * canopy%fe ) ) /          &
              ( air%rho * CCAPP * met%tk * canopy%us**3 )

         ! stability parameter at shear height: needed for Harman in-canopy stability correction
         IF (cable_user%soil_struc=='sli') THEN
            WHERE (canopy%vlaiw > CLAI_THRESH .AND. rough%hruff > rough%z0soilsn)
               canopy%zetash(:,iterplus) = -(Cvonk*Cgrav*(0.1*rough%hruff)*(canopy%fhs+0.07*REAL(canopy%fes)))/ &
                    MAX( (air%rho*Ccapp*met%tk*(canopy%us*rough%term6a)**3), 1.e-12)
            ELSEWHERE
               canopy%zetash(:,iterplus) = canopy%zetash(:,iter)
            ENDWHERE
         ENDIF

         ! case NITER=2: final zetar=CZETmul*zetar(2) (compute only when iter=1)
         IF (NITER == 2) THEN

            canopy%zetar(:,2) = CZETmul * canopy%zetar(:,2)

            ! stability parameter at shear height: needed for Harman in-canopy stability correction
            IF (cable_user%soil_struc=='sli') THEN
               canopy%zetash(:,2) = CZETmul * canopy%zetash(:,2)
            ENDIF


            DO j=1,mp
               IF ( (met%fsd(j,1)+met%fsd(j,2))  ==  0.0 ) &
                    canopy%zetar(j,2) = 0.5 * canopy%zetar(j,2)
            ENDDO

         END IF

         ! constrain zeta to CZETPOS and CZETNEG (set in param0)

         ! zetar too +
         canopy%zetar(:,iterplus) = MIN(CZETPOS,canopy%zetar(:,iterplus))
         !jhan: to get past rigorous build - however (:,i) cant be compared
         !if ( canopy%zetash(:,iterplus) .NE. CZETPOS ) &
         IF (cable_user%soil_struc=='sli') &
              canopy%zetash(:,iterplus) = MIN(CZETPOS,canopy%zetash(:,iterplus))

         ! zetar too -
         canopy%zetar(:,iterplus) = MAX(CZETNEG,canopy%zetar(:,iterplus))
         !if ( canopy%zetash(:,iterplus) .NE. CZETNEG ) &
         IF (cable_user%soil_struc=='sli') &
              canopy%zetash(:,iterplus) = MAX(CZETNEG,canopy%zetash(:,iterplus))

      END IF ! (iter < NITER)

    END SUBROUTINE update_zetar


    FUNCTION qsatf(j,tair,pmb) RESULT(r)
      ! MRR, 1987
      ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
      ! HUMIDITY (KG/KG) FROM TETEN FORMULA

      REAL, INTENT(IN) ::                                                         &
           tair,         & ! air temperature (C)
           pmb             ! pressure PMB (mb)

      INTEGER, INTENT(IN) :: j

      REAL           :: r    ! result; sat sp humidity

      r = (CRMH2o/Crmair) * (CTETENA*EXP(CTETENB*tair/(CTETENC+tair))) / pmb

    END FUNCTION qsatf

    ! -----------------------------------------------------------------------------

    SUBROUTINE qsatfjh(var,tair,pmb)
      USE cable_def_types_mod, ONLY : mp
      REAL, INTENT(IN), DIMENSION(mp) ::                                          &
           tair,                        & ! air temperature (C)
           pmb                            ! pressure PMB (mb)

      REAL, INTENT(OUT), DIMENSION(mp) ::                                         &
           var                            ! result; sat sp humidity

      INTEGER :: j

      DO j=1,mp

         var(j) = (CRMH2o/Crmair) * (CTETENA*EXP(CTETENB*tair(j)/(CTETENC+tair(j))))    &
              / pmb(j)
      ENDDO

    END SUBROUTINE qsatfjh

    ! -----------------------------------------------------------------------------

    SUBROUTINE qsatfjh2(var,tair,pmb)

      REAL, INTENT(IN) ::                                                         &
           tair,         & ! air temperature (C)
           pmb             ! pressure PMB (mb)

      REAL, INTENT(OUT) ::                                                        &
           var             ! result; sat sp humidity

      var = (CRMH2o/Crmair) * (CTETENA*EXP(CTETENB*tair/(CTETENC+tair))) / pmb

    END SUBROUTINE qsatfjh2

    ELEMENTAL FUNCTION psis(zeta) RESULT(r)

      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psis(z/l) (z/l=zeta)
      ! for scalars, using the businger-dyer form for unstable cases
      ! and the webb form for stable cases. see paulson (1970).

      REAL, INTENT(IN)     :: zeta

      REAL, PARAMETER      ::                                                     &
           gu = 16.0,        & !
           gs = 5.0,         & !
           a = 1.0,          & !
           b = 0.667,        & !
           c = 5.0,          & !
           d = 0.35

      REAL                 ::                                                     &
           r,                & !
           stzeta,           & !
           ustzeta,          & !
           z,                & !
           y,                & !
           stable,           & !
           unstable

      z      = 0.5 + SIGN(0.5,zeta)    ! z=1 in stable, 0 in unstable

      ! Beljaars and Holtslag (1991) for stable
      stzeta = MAX(0.,zeta)
      stable = -(1.+2./3.*a*stzeta)**(3./2.) -  &
           b*(stzeta-c/d)*EXP(-d*stzeta) - b*c/d + 1.
      y      = (1.0 + gu*ABS(zeta))**0.5
      unstable = 2.0 * alog((1+y)*0.5)
      r   = z*stable + (1.0-z)*unstable

    END FUNCTION psis

    ! -----------------------------------------------------------------------------

    ELEMENTAL FUNCTION rplant(rpconst, rpcoef, tair) RESULT(z)
      REAL, INTENT(IN)     :: rpconst
      REAL, INTENT(IN)     :: rpcoef
      REAL, INTENT(IN)     :: tair
      REAL                 :: z
      z = rpconst * EXP(rpcoef * tair)
    END FUNCTION rplant

    ! -----------------------------------------------------------------------------

    SUBROUTINE wetLeaf( dels, rad, rough, air, met, veg, canopy, cansat, tlfy,     &
         gbhu, gbhf, ghwet )

      USE cable_def_types_mod

      TYPE (radiation_type), INTENT(INOUT) :: rad
      TYPE (roughness_type), INTENT(INOUT) :: rough
      TYPE (air_type),       INTENT(INOUT) :: air
      TYPE (met_type),       INTENT(INOUT) :: met
      TYPE (canopy_type),    INTENT(INOUT) :: canopy

      TYPE (veg_parameter_type), INTENT(INOUT)    :: veg

      REAL,INTENT(IN), DIMENSION(:) ::                                            &
           tlfy,          & ! leaf temp (K) - assCUMINg the temperature of
                                ! wet leaf is equal that of dry leaf ="tlfy"
           cansat           ! max canopy intercept. (mm)

      REAL(r_2), INTENT(IN), DIMENSION(:,:) ::                                    &
           gbhu,          & ! forcedConvectionBndryLayerCond
           gbhf             ! freeConvectionBndryLayerCond

      REAL(r_2), INTENT(OUT), DIMENSION(:) ::                                     &
           ghwet            ! cond for heat for a wet canopy

      REAL, INTENT(IN)     :: dels ! integration time step (s)

      ! local variables
      REAL, DIMENSION(mp) ::                                                      &
           ccfevw,        & ! limitation term for
           gwwet,         & ! cond for water for a wet canopy
           ghrwet           ! wet canopy cond: heat & thermal rad

      !i sums, terms of convenience/readability
      REAL, DIMENSION(mp) ::                                                      &
           sum_gbh, sum_rad_rniso, sum_rad_gradis, xx1

      INTEGER :: j

      ! END header

      ghwet = 1.0e-3
      gwwet = 1.0e-3
      ghrwet= 1.0e-3
      canopy%fevw = 0.0
      canopy%fhvw = 0.0
      sum_gbh = SUM((gbhu+gbhf),2)
      sum_rad_rniso = SUM(rad%rniso,2)
      sum_rad_gradis = SUM(rad%gradis,2)

      DO j=1,mp

         IF(canopy%vlaiw(j) > CLAI_THRESH) THEN

            ! VEG SENSIBLE & LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
            ! calculate total thermal resistance, rthv in s/m
            ghwet(j) = 2.0   * sum_gbh(j)
            gwwet(j) = 1.075 * sum_gbh(j)
            ghrwet(j) = sum_rad_gradis(j) + ghwet(j)

            ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
            ! to avoid excessive evaporation:
            ccfevw(j) = MIN(canopy%cansto(j) * air%rlam(j) / dels, &
                 2.0 / (1440.0 / (dels/60.0)) * air%rlam(j) )

            canopy%fevw(j) = MIN( canopy%fwet(j) * ( air%dsatdk(j) *              &
                 ( sum_rad_rniso(j)- CCAPP*Crmair*( met%tvair(j)     &
                 - met%tk(j) ) * sum_rad_gradis(j) )                   &
                 + CCAPP * Crmair * met%dva(j) * ghrwet(j) )         &
                 / ( air%dsatdk(j)+air%psyc(j)*ghrwet(j) / gwwet(j) )  &
                 , ccfevw(j) )

            !upwards flux density of water (kg/m2/s) - canopy componenet of 
            !potential evapotranspiration 
            canopy%fevw_pot(j) = ( air%dsatdk(j)* (sum_rad_rniso(j) -             &
                 CCAPP * Crmair * ( met%tvair(j) - met%tk(j) )  &
                 *sum_rad_gradis(j) )                             &
                 + CCAPP * Crmair * met%dva(j) * ghrwet(j))     &
                 / (air%dsatdk(j)+air%psyc(j)*ghrwet(j)/gwwet(j) )

            ! calculate sens heat from wet canopy:
            canopy%fhvw(j) = canopy%fwet(j) * ( sum_rad_rniso(j) -CCAPP * Crmair&
                 * ( tlfy(j) - met%tk(j) ) * sum_rad_gradis(j) )      &
                 - canopy%fevw(j)

         ENDIF

      ENDDO

    END SUBROUTINE wetLeaf

END SUBROUTINE define_canopy

END MODULE cable_canopy_module
