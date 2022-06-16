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

USE cable_photo_constants_mod, ONLY : CGSW03 => GSW03
USE cable_photo_constants_mod, ONLY : CGSW04 => GSW04
USE cable_photo_constants_mod, ONLY : CCONKC0=> CONKC0
USE cable_photo_constants_mod, ONLY : CCONKO0=> CONKO0
USE cable_photo_constants_mod, ONLY : CEKC   => EKC
USE cable_photo_constants_mod, ONLY : CEKO   => EKO
USE cable_photo_constants_mod, ONLY : CD0C3   => D0C3
USE cable_photo_constants_mod, ONLY : CD0C4   => D0C4
USE cable_photo_constants_mod, ONLY : CA1C3   => A1C3
USE cable_photo_constants_mod, ONLY : CA1C4   => A1C4
USE cable_photo_constants_mod, ONLY : CCFRD3 => CFRD3
USE cable_photo_constants_mod, ONLY : CCFRD4 => CFRD4
USE cable_photo_constants_mod, ONLY : CALPHA3 => ALPHA3  
USE cable_photo_constants_mod, ONLY : CALPHA4 => ALPHA4  
USE cable_photo_constants_mod, ONLY : CCONVX3 => CONVX3
USE cable_photo_constants_mod, ONLY : CCONVX4 => CONVX4

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
USE cable_climate_type_mod, ONLY : climate_type

USE cbl_friction_vel_module,  ONLY : comp_friction_vel, psim, psis
USE cbl_pot_evap_snow_module, ONLY : Penman_Monteith, Humidity_deficit_method
USE cbl_qsat_module,          ONLY : qsatfjh,  qsatfjh2
USE cbl_zetar_module,         ONLY : update_zetar
USE cable_latent_heat_module, ONLY : latent_heat_flux
USE cable_wetleaf_module,     ONLY : wetleaf 
USE cbl_dryLeaf_module,       ONLY : dryLeaf
USE cbl_SurfaceWetness_module, ONLY : Surf_wetness_fact
USE cable_within_canopy_module, ONLY : within_canopy

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
      relitt           !
   
   ! temporary buffers to simplify equations
   REAL, DIMENSION(mp) ::                                                      &
      ftemp,rlower_limit,                      &
      term1, term2, term3, term5 

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

   INTEGER :: j
   
   INTEGER, SAVE :: call_number =0
   
real :: tv_denom(mp)
real :: tv_denom_transd(mp)
real :: tv_frac(mp)
real :: tv_tk(mp)
real :: tv(mp)

   ! END header
   
   call_number = call_number + 1
           
   IF( .NOT. cable_runtime%um)                                                 &
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

   ! Initialise in-canopy temperatures and humidity:
   csx = SPREAD(met%ca, 2, mf) ! initialise leaf surface CO2 concentration
   met%tvair = met%tk
   met%qvair = met%qv
   canopy%tv = met%tvair

   CALL define_air (met, air)
   
   !ESM15!CALL qsatfjh(qstvair,met%tvair-Ctfrz,met%pmb)
CALL qsatfjh(mp, qstvair, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC, met%tvair-CTfrz,met%pmb)

   met%dva = (qstvair - met%qvair) *  Crmair/Crmh2o * met%pmb * 100.0
   dsx = met%dva     ! init. leaf surface vpd
   
   tlfx = met%tk  ! initialise leaf temp iteration memory variable (K)
   tlfy = met%tk  ! initialise current leaf temp (K)
   
   ortsoil = ssnow%rtsoil
   ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)
   tss4 = ssnow%tss**4
   canopy%fes = 0.
   canopy%fess = 0.
   canopy%fesp = 0.
   ssnow%potev = 0.
   canopy%fevw_pot = 0.

CALL radiation( ssnow, veg, air, met, rad, canopy, sunlit_veg_mask, &
  !constants
  clai_thresh, Csboltz, Cemsoil, Cemleaf, Ccapp &
)

   canopy%zetar(:,1) = CZETA0 ! stability correction terms
   canopy%zetar(:,2) = CZETPOS + 1 


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
      xx = 0.5 + sign(0.5,rough%zref_tq+rough%disp-rough%zruffs)
      
      ! correction  by Ian Harman to the 2nd psis term
      rt1usc = xx * ( LOG( rough%zref_tq/MAX( rough%zruffs-rough%disp,         &
                                              rough%z0soilsn ) )               &
               - psis( canopy%zetar(:,iter) )                                  &
               + psis( canopy%zetar(:,iter) * ( MAX( rough%zruffs-rough%disp,  &
                                                rough%z0soilsn ) )            &
                       / rough%zref_tq ) ) / CVONK
      
      rt_min = 5.      
      rt0 = max(rt_min,rough%rt0us / canopy%us)
      
      ! Aerodynamic resistance (sum 3 height integrals)/us
      ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
      rough%rt1 = MAX(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)
      
      DO j=1,mp
     
         IF(canopy%vlaiw(j) > CLAI_THRESH) THEN
            ssnow%rtsoil(j) = rt0(j)
         ELSE
            ssnow%rtsoil(j) = rt0(j) + rough%rt1(j)
         ENDif 
     
      ENDDO 
      
      ssnow%rtsoil = max(rt_min,ssnow%rtsoil)   
      
      DO j=1,mp
      
         IF( ssnow%rtsoil(j) > 2.*ortsoil(j) .OR.                              &
             ssnow%rtsoil(j) < 0.5*ortsoil(j) ) THEN
              
            ssnow%rtsoil(j) = MAX(rt_min,0.5*(ssnow%rtsoil(j) + ortsoil(j)))
        
        ENDIF    
      
      ENDDO 
   
      ! Vegetation boundary-layer conductance (mol/m2/s)
      ! Cprandt = kinematic viscosity/molecular diffusivity
      ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
      DO j=1,mp

         IF(canopy%vlaiw(j) > CLAI_THRESH) THEN
            gbvtop(j) = air%cmolar(j)*CAPOL * air%visc(j) / Cprandt /        &
                        veg%dleaf(j) * (canopy%us(j) / MAX(rough%usuh(j),1.e-6)&
                        * veg%dleaf(j) / air%visc(j) )**0.5                    &
                        * Cprandt**(1.0/3.0) / veg%shelrb(j)
            gbvtop(j) = MAX (0.05,gbvtop(j) )      ! for testing (BP aug2010)
            
            ! Forced convection boundary layer conductance                     
            ! (see Wang & Leuning 1998, AFM):
            gbhu(j,1) = gbvtop(j)*(1.0-EXP(-canopy%vlaiw(j)                    &
                        *(0.5*rough%coexp(j)+rad%extkb(j) ))) /                &
                        (rad%extkb(j)+0.5*rough%coexp(j))
            
            gbhu(j,2) = (2.0/rough%coexp(j))*gbvtop(j)*  &
                        (1.0-EXP(-0.5*rough%coexp(j)*canopy%vlaiw(j)))         &
                        - gbhu(j,1)
         ENDIF 
      
      ENDDO 
     
      rny = SUM(rad%rniso,2) ! init current estimate net rad
      hcy = 0.0              ! init current estimate lat heat
      ecy = rny - hcy        ! init current estimate lat heat

sum_rad_rniso = SUM(rad%rniso,2)
sum_rad_gradis = SUM(rad%gradis,2)

      CALL dryLeaf( dels, rad, rough, air, met,                                &
                    veg, canopy, soil, ssnow, dsx,                             &
                    fwsoil, tlfx, tlfy, ecy, hcy,                              &
                    rny, gbhu, gbhf, csx, cansat,                              &
                    ghwet,  iter, climate )
     
      CALL wetLeaf( dels,                                 &
                    cansat, tlfy,                                 &
                    gbhu, gbhf, ghwet, &
                    mp, CLAI_thresh, CCAPP, CRmair, & 
                    reducedLAIdue2snow, sum_rad_rniso, sum_rad_gradis, & 
                    canopy%fevw, canopy%fevw_pot, canopy%fhvw, &
                    canopy%fwet, canopy%cansto, air%rlam, air%dsatdk, &
                    met%tvair, met%tk, met%dva, air%psyc )

     
      ! Calculate latent heat from vegetation:
      ! Calculate sensible heat from vegetation:
      ! Calculate net rad absorbed by canopy:
      canopy%fev = REAL(canopy%fevc + canopy%fevw)
      ftemp = (1.0 - canopy%fwet) *  REAL(hcy) + canopy%fhvw
      canopy%fhv = real(ftemp) 
      ftemp= (1.0-canopy%fwet)*REAL(rny)+canopy%fevw+canopy%fhvw
      canopy%fnv = real(ftemp)

      ! canopy rad. temperature calc from long-wave rad. balance
      sum_rad_gradis = SUM(rad%gradis,2)

      DO j=1,mp

         IF ( canopy%vlaiw(j) > CLAI_THRESH .AND.                             &
              rough%hruff(j) > rough%z0soilsn(j) ) THEN

            rad%lwabv(j) = CCAPP * Crmair * ( tlfy(j) - met%tk(j) ) *        &
                           sum_rad_gradis(j) 

!esm1.5 using LAI_thresh=0.001 formulation for canopy%tv fails 
!esm1.5            canopy%tv(j) = (rad%lwabv(j) / (2.0*(1.0-rad%transd(j))            &
!esm1.5                     * CSBOLTZ*CEMLEAF)+met%tk(j)**4)**0.25

tv_denom_transd(j) = ( 1.0-rad%transd(j) )
tv_denom(j) = 2.0*(tv_denom_transd(j)) * CSBOLTZ*CEMLEAF
tv_frac(j) = rad%lwabv(j) / tv_denom(j)
tv_tk(j) = met%tk(j) **4
tv(j) = ( tv_frac(j) + tv_tk(j) )
!effectivel MAX here does IF( |tv_frac| >  tv_tk ) tv_frac = -1.0 * tv_tk
tv(j) = MAX( 0.0, tv(j) )
tv(j) = tv(j) ** .25
canopy%tv(j) = tv(j)

         
         ELSE! sparse canopy
         
           canopy%tv(j) = met%tk(j)
         
         ENDIF
          
      ENDDO 
     

      ! Calculate net rad to soil:
      canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*CEMLEAF* &
            CSBOLTZ*canopy%tv**4 - CEMSOIL*CSBOLTZ* tss4


      ! Saturation specific humidity at soil/snow surface temperature:
      !ESM15!call qsatfjh(ssnow%qstss,ssnow%tss-Ctfrz,met%pmb)
       CALL qsatfjh(mp, ssnow%qstss, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC,ssnow%tss-CTfrz,met%pmb)

      IF(cable_user%ssnow_POTEV== "P-M") THEN
         
         !--- uses %ga from previous timestep    
         ssnow%potev =  Penman_Monteith( mp, Ctfrz, CRMH2o, Crmair, CTETENA, CTETENB,         &
                 CTETENC, air%dsatdk, air%psyc, air%rho, air%rlam,             & 
                 met%tvair, met%pmb, met%qvair,                       &
                 canopy%ga, canopy%fns, ssnow%rtsoil, ssnow%isflag )
      ELSE !by default assumes Humidity Deficit Method
      
         dq = ssnow%qstss - met%qv
         ssnow%potev = Humidity_deficit_method( mp, Ctfrz,  &
                                 air%rho,air%rlam,           & 
                      dq, ssnow%rtsoil,ssnow%snowd, ssnow%tgg )
          
      ENDIF

      ! Soil latent heat:

      CALL Latent_heat_flux( mp, CTFRZ, dels, soil%zse(1), soil%swilt,           &
                             cable_user%l_new_reduce_soilevp, pwet, air%rlam,  &
                             ssnow%snowd, ssnow%wb(:,1), ssnow%wbice(:,1),             &
                             ssnow%pudsto, ssnow%pudsmx, ssnow%potev,          &
                             ssnow%wetfac, ssnow%evapfbl(:,1), ssnow%cls,          & 
                             ssnow%tss, canopy%fes, canopy%fess, canopy%fesp  )

      ! Calculate soil sensible heat:
      canopy%fhs = air%rho*CCAPP*(ssnow%tss - met%tk) /ssnow%rtsoil


  CALL within_canopy( mp, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC,           &
                      CLAI_thresh, CCAPP, CTFRZ, rad, rough, air, met, veg,    &
                      canopy, ssnow, gbhu, gbhf, qstvair, rt0, rhlitt, relitt )

      ! Saturation specific humidity at soil/snow surface temperature:
      !ESM15!call qsatfjh(ssnow%qstss,ssnow%tss-Ctfrz,met%pmb)
       CALL qsatfjh(mp, ssnow%qstss, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC,ssnow%tss-Ctfrz,met%pmb)

      IF(cable_user%ssnow_POTEV== "P-M") THEN
         
         !--- uses %ga from previous timestep    
         ssnow%potev =  Penman_Monteith( mp, Ctfrz, CRMH2o, Crmair, CTETENA, CTETENB,         &
                 CTETENC, air%dsatdk, air%psyc, air%rho, air%rlam,             & 
                 met%tvair, met%pmb, met%qvair,                       &
                 canopy%ga, canopy%fns, ssnow%rtsoil, ssnow%isflag )
      
      ELSE !by default assumes Humidity Deficit Method
      
         dq = ssnow%qstss - met%qvair
         ssnow%potev = Humidity_deficit_method( mp, Ctfrz,  &
                                 air%rho,air%rlam,           & 
                      dq, ssnow%rtsoil,ssnow%snowd, ssnow%tgg )
          
      ENDIF

         
      ! Soil latent heat:
      CALL Latent_heat_flux( mp, CTFRZ, dels, soil%zse(1), soil%swilt,           &
                             cable_user%l_new_reduce_soilevp, pwet, air%rlam,  &
                             ssnow%snowd, ssnow%wb(:,1), ssnow%wbice(:,1),             &
                             ssnow%pudsto, ssnow%pudsmx, ssnow%potev,          &
                             ssnow%wetfac, ssnow%evapfbl(:,1), ssnow%cls,          & 
                             ssnow%tss, canopy%fes, canopy%fess, canopy%fesp  )


      ! Soil sensible heat:
      canopy%fhs = air%rho*CCAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
      !canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssnow%cls
      canopy%ga = canopy%fns-canopy%fhs-canopy%fes
      
      ! Set total latent heat:
      canopy%fe = canopy%fev + canopy%fes
      
      ! Set total sensible heat:
      canopy%fh = canopy%fhv + canopy%fhs
      
      !---diagnostic purposes
      DO j=1,mp
      
         IF (ssnow%potev(j) .GE. 0.) THEN
            ssnow%potev(j) = max(0.00001,ssnow%potev(j))
         ELSE
            ssnow%potev(j) = min(-0.0002,ssnow%potev(j))
         ENDIF
         
         IF (canopy%fevw_pot(j) .ge. 0.) then
            canopy%fevw_pot(j) = max(0.000001,canopy%fevw_pot(j))
         ELSE
            canopy%fevw_pot(j) = min(-0.002,canopy%fevw_pot(j))
         ENDIF

      ENDDO 


      canopy%rnet = canopy%fnv + canopy%fns  
      canopy%epot = ((1.-rad%transd)*canopy%fevw_pot +                         &
                    rad%transd*ssnow%potev) * dels/air%rlam  

      rlower_limit = canopy%epot * air%rlam / dels  
      where (rlower_limit == 0 ) rlower_limit = 1.e-7 !prevent from 0. by adding 1.e-7 (W/m2)

      
      canopy%wetfac_cs = max(0., min(1.0,canopy%fe / rlower_limit ))
      
      DO j=1,mp

         IF ( canopy%wetfac_cs(j) .LE. 0. )                                    &
            canopy%wetfac_cs(j) = MAX( 0., MIN( 1.,                            &
                                  MAX( canopy%fev(j)/canopy%fevw_pot(j),       &
                                  canopy%fes(j)/ssnow%potev(j) ) ) )
      
      ENDDO 

  CALL update_zetar( mp, NITER, canopy%zetar, iter, nrb, CVONK, CGRAV, CCAPP,  &
                     CLAI_THRESH, CZETmul, CZETPOS, CZETNEG,          &
                     cable_user%soil_struc, air%rho, met%tk,  met%fsd, &
                     rough%zref_tq, rough%hruff, rough%term6a, rough%z0soilsn,   &
                     canopy%vlaiw, canopy%zetash,  canopy%us, &
                     canopy%fh, canopy%fe, canopy%fhs, REAL(canopy%fes) ) 

   END DO           ! do iter = 1, NITER


   canopy%cduv = canopy%us * canopy%us / (max(met%ua,CUMIN))**2

   !---diagnostic purposes
   canopy%gswx_T = rad%fvlai(:,1)/MAX( CLAI_THRESH, canopy%vlaiw(:) )         & 
                   * canopy%gswx(:,1) + rad%fvlai(:,2) / MAX(CLAI_THRESH,     &
                   canopy%vlaiw(:))*canopy%gswx(:,2)

    ! The surface conductance below is required by dust scheme; it is composed from canopy and soil conductances
    canopy%gswx_T = (1.-rad%transd)*max(1.e-06,canopy%gswx_T ) +  &   !contribution from  canopy conductance
                  rad%transd*(.01*ssnow%wb(:,1)/soil%sfc)**2 ! + soil conductance; this part is done as in Moses
    where ( soil%isoilm == 9 ) canopy%gswx_T = 1.e6   ! this is a value taken from Moses for ice points

    canopy%cdtq = canopy%cduv *( LOG( rough%zref_uv / rough%z0m) -              &
               psim( canopy%zetar(:,NITER) * rough%zref_uv/rough%zref_tq, mp, CPI_C ) &
             + psim( canopy%zetar(:,NITER) * rough%z0m/rough%zref_tq, mp, CPI_C ) & ! new term from Ian Harman
                 ) / ( LOG( rough%zref_tq /(0.1*rough%z0m) )                   &
               - psis( canopy%zetar(:,NITER))                                  &
               + psis(canopy%zetar(:,NITER)*0.1*rough%z0m/rough%zref_tq) ) ! n


   ! Calculate screen temperature: 1) original method from SCAM
   ! screen temp., windspeed and relative humidity at 1.5m
   ! screen temp., windspeed and relative humidity at 2.0m
    tstar = - canopy%fh / ( air%rho*CCAPP*canopy%us)
    qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us)
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
   r_sc = 0.
   zscl = MAX(rough%z0soilsn,2.0)

   ! assume screen temp of bareground if all these conditions are not met
   DO j=1,mp
      
      IF ( canopy%vlaiw(j) > CLAI_THRESH .and. rough%hruff(j) > 0.01) THEN
      
         IF ( rough%disp(j)  > 0.0 ) then
     
            term1(j) = EXP(2*CCSW*canopy%rghlai(j)*(1-zscl(j)/rough%hruff(j)))
            term2(j) = EXP(2*CCSW*canopy%rghlai(j) *                          &
                       (1-rough%disp(j)/rough%hruff(j)))
            term5(j) = MAX(2./3.*rough%hruff(j)/rough%disp(j), 1.)
         
         ENDIF
        
         term3(j) = CA33**2*CCTL*2*CCSW*canopy%rghlai(j)

         IF( zscl(j) < rough%disp(j) ) THEN

            r_sc(j) = term5(j) * LOG(zscl(j)/rough%z0soilsn(j)) *              &
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
            
            r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j) +     &
                      ( LOG( (zscl(j) - rough%disp(j)) /                       &
                      MAX( rough%zruffs(j)-rough%disp(j),                      &
                      rough%z0soilsn(j) ) ) - psis1( (zscl(j)-rough%disp(j))   &
                      / (rough%zref_tq(j)/canopy%zetar(j,iterplus) ) )         &
                      + psis1( (rough%zruffs(j) - rough%disp(j) )              &
                      / (rough%zref_tq(j)/canopy%zetar(j,iterplus ) ) ) )      &
                      / CVONK

         ENDIF

        canopy%tscrn(j) = ssnow%tss(j) + (met%tk(j) - ssnow%tss(j)) *          &
                          MIN(1.,r_sc(j) / MAX( 1.,                            &
                          rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)   &
                          + rt1usc(j))) - Ctfrz 
      ENDIF  

   ENDDO  
  
   !ESM15!CALL qsatfjh(rsts,canopy%tscrn,met%pmb)
    CALL  qsatfjh(mp, rsts, CRMH2o, Crmair, CTETENA, CTETENB, CTETENC, canopy%tscrn,met%pmb)
     
   qtgnet = rsts * ssnow%wetfac - met%qv
   
   DO j=1,mp
      
      IF (qtgnet(j) .GT. 0. ) THEN
         qsurf(j) = rsts(j) * ssnow%wetfac(j)
      ELSE
         qsurf(j) = 0.1*rsts(j)*ssnow%wetfac(j) + 0.9*met%qv(j)
      ENDIF

      canopy%qscrn(j) = met%qv(j) - qstar(j) * ftemp(j)

      IF( canopy%vlaiw(j) >CLAI_THRESH .and. rough%hruff(j) > 0.01)           &

            canopy%qscrn(j) = qsurf(j) + (met%qv(j) - qsurf(j)) * MIN( 1.,     &
                              r_sc(j) / MAX( 1., rough%rt0us(j) +              &
                              rough%rt1usa(j) + rough%rt1usb(j) + rt1usc(j) ) )

   ENDDO 


   ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
   canopy%dewmm = - (min(0.0,canopy%fevw) + min(0.0_r_2,canopy%fevc)) * &
        dels * 1.0e3 / (CRHOW*air%rlam)

   ! Add dewfall to canopy water storage:
   canopy%cansto = canopy%cansto + canopy%dewmm
   
   ! Modify canopy water storage for evaporation:
   canopy%cansto = MAX(canopy%cansto-MAX(0.0,REAL(canopy%fevw))*dels &
     *1.0e3/(CRHOW*air%rlam), 0.0)

   ! Calculate canopy water storage excess:
   canopy%spill=max(0.0, canopy%cansto-cansat)

   ! Move excess canopy water to throughfall:
   canopy%through = canopy%through + canopy%spill
   
   ! Initialise 'throughfall to soil' as 'throughfall from canopy'; 
   ! snow may absorb
   canopy%precis = max(0.,canopy%through)

   ! this change of units does not affect next timestep as canopy%through is
   ! re-calc in surf_wetness_fact routine
   canopy%through = canopy%through / dels   ! change units for stash output
   
   ! Update canopy storage term:
   canopy%cansto=canopy%cansto - canopy%spill
   
   ! Calculate the total change in canopy water store (mm/dels):
   canopy%delwc = canopy%cansto-canopy%oldcansto
   
   ! calculate dgdtg, derivative of ghflux 3 instances
   ! d(canopy%fns)/d(ssnow%tgg)
   ! d(canopy%fhs)/d(ssnow%tgg)
   ! d(canopy%fes)/d(dq)
   ssnow%dfn_dtg = (-1.)*4.*CEMSOIL*CSBOLTZ*tss4/ssnow%tss  
   ssnow%dfh_dtg = air%rho*CCAPP/ssnow%rtsoil      
   ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/ssnow%rtsoil  
  
   ssnow%ddq_dtg = (Crmh2o/Crmair) /met%pmb * CTETENA*CTETENB * CTETENC   &
                   / ( ( CTETENC + ssnow%tss-Ctfrz )**2 )*EXP( CTETENB *       &
                   ( ssnow%tss-Ctfrz ) / ( CTETENC + ssnow%tss-Ctfrz ) )
   canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%dfe_ddq *    &
                  ssnow%ddq_dtg

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

CONTAINS

!psis1 function obsoleted in trunk
! -----------------------------------------------------------------------------

FUNCTION psis1(zeta) RESULT(r)
   ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
   ! computes integrated stability function psis(z/l) (z/l=zeta)
   ! for scalars, using the businger-dyer form for unstable cases
   ! and the webb form for stable cases. see paulson (1970).
   REAL, INTENT(IN)     :: zeta
   
   REAL, PARAMETER      :: gu = 16.0
   REAL, PARAMETER      :: gs = 5.0
   REAL, PARAMETER      :: a = 1.0
   REAL, PARAMETER      :: b = 0.667
   REAL, PARAMETER      :: c = 5.0
   REAL, PARAMETER      :: d = 0.35
 
   REAL                 :: r
   REAL                 :: stable
   REAL                 :: unstable
   REAL                 :: stzeta
 
   REAL                 :: z
   REAL                 :: y
   !REAL                 :: stable
   !REAL                 :: unstable
 
   z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable 
   
   ! Beljaars and Holtslag (1991) for stable
   stzeta = max(0.,zeta)
   stable = -(1.+2./3.*a*stzeta)**(3./2.) -  &
             b*(stzeta-c/d)*exp(-d*stzeta) - b*c/d + 1.
 
   y      = (1.0 + gu*abs(zeta))**0.5
   unstable = 2.0 * alog((1+y)*0.5)
   r   = z*stable + (1.0-z)*unstable

END FUNCTION psis1

END SUBROUTINE define_canopy

END MODULE cable_canopy_module

