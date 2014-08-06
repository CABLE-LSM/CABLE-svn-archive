!break canopy down into modules in canopy directory.
!first photosynthesis module
MODULE cable_canopy_module
   
   USE cable_data_module, ONLY : icanopy_type, point2constants 
   USE cable_canopy_dryleaf_mod, only : dryleaf
   
   IMPLICIT NONE
   
   PUBLIC define_canopy
   PRIVATE
   
   TYPE( icanopy_type ) :: C
  
     
CONTAINS
 

SUBROUTINE define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy)
   USE cable_def_types_mod
   USE cable_radiation_module
   USE cable_air_module
   USE cable_common_module   
   USE cable_roughness_module

   TYPE (balances_type), INTENT(INOUT)  :: bal
   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (air_type), INTENT(INOUT)       :: air
   TYPE (met_type), INTENT(INOUT)       :: met
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (canopy_type), INTENT(INOUT)    :: canopy

   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
   
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
      xx1,           & !
      sum_rad_rniso, & ! 
      sum_rad_gradis   ! 
   
   ! temporary buffers to simplify equations
   REAL, DIMENSION(mp) ::                                                      &
      ftemp,z_eff,psim_arg, psim_1, psim_2, rlower_limit,                      &
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
   
   ! END header
   
   call_number = call_number + 1
           
   ! assign local ptrs to constants defined in cable_data_module
   CALL point2constants(C)    

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
   
   CALL qsatfjh(qstvair,met%tvair-C%tfrz,met%pmb)

   met%dva = (qstvair - met%qvair) *  C%rmair/C%rmh2o * met%pmb * 100.0
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

   CALL radiation( ssnow, veg, air, met, rad, canopy )

   canopy%zetar(:,1) = C%ZETA0 ! stability correction terms
   canopy%zetar(:,2) = C%ZETPOS + 1 


   DO iter = 1, NITER

      ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
      ! resistances rt0, rt1 (elements of dispersion matrix):
      ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
      CALL comp_friction_vel()

      ! E.Kowalczyk 2014
      IF (cable_user%l_new_roughness_soil)                                     &
         CALL ruff_resist(veg, rough, ssnow, canopy)

      
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
                       / rough%zref_tq ) ) / C%VONK
      
      rt_min = 5.      
      rt0 = max(rt_min,rough%rt0us / canopy%us)
      
      ! Aerodynamic resistance (sum 3 height integrals)/us
      ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
      rough%rt1 = MAX(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)
      
      DO j=1,mp
     
         IF(canopy%vlaiw(j) > C%LAI_THRESH) THEN
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
      ! C%prandt = kinematic viscosity/molecular diffusivity
      ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
      DO j=1,mp

         IF(canopy%vlaiw(j) > C%LAI_THRESH) THEN
            gbvtop(j) = air%cmolar(j)*C%APOL * air%visc(j) / C%prandt /        &
                        veg%dleaf(j) * (canopy%us(j) / MAX(rough%usuh(j),1.e-6)&
                        * veg%dleaf(j) / air%visc(j) )**0.5                    &
                        * C%prandt**(1.0/3.0) / veg%shelrb(j)
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

      CALL dryLeaf( dels, rad, rough, air, met,                                &
                    veg, canopy, soil, ssnow, dsx,                             &
                    fwsoil, tlfx, tlfy, ecy, hcy,                              &
                    rny, gbhu, gbhf, csx, cansat,                              &
                    ghwet,  iter )
     
      CALL wetLeaf( dels, rad, rough, air, met,                                &
                    veg, canopy, cansat, tlfy,                                 &
                    gbhu, gbhf, ghwet )

     
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

         IF ( canopy%vlaiw(j) > C%LAI_THRESH .AND.                             &
              rough%hruff(j) > rough%z0soilsn(j) ) THEN

            rad%lwabv(j) = C%CAPP * C%rmair * ( tlfy(j) - met%tk(j) ) *        &
                           sum_rad_gradis(j) 

            canopy%tv(j) = (rad%lwabv(j) / (2.0*(1.0-rad%transd(j))            &
                           * C%SBOLTZ*C%EMLEAF)+met%tk(j)**4)**0.25
         
         ELSE! sparse canopy
         
           canopy%tv(j) = met%tk(j)
         
         ENDIF
          
      ENDDO 
     

      ! Calculate net rad to soil:
      canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*C%EMLEAF* &
            C%SBOLTZ*canopy%tv**4 - C%EMSOIL*C%SBOLTZ* tss4


      ! Saturation specific humidity at soil/snow surface temperature:
      call qsatfjh(ssnow%qstss,ssnow%tss-C%tfrz,met%pmb)

      IF(cable_user%ssnow_POTEV== "P-M") THEN
         
         !--- uses %ga from previous timestep    
         ssnow%potev =  Penman_Monteith(canopy%ga) 
      
      ELSE !by default assumes Humidity Deficit Method
      
         dq = ssnow%qstss - met%qv
         ssnow%potev =  Humidity_deficit_method(dq,ssnow%qstss )
          
      ENDIF

      ! Soil latent heat:
      CALL latent_heat_flux()

      ! Calculate soil sensible heat:
      canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tk) /ssnow%rtsoil


      CALL within_canopy( gbhu, gbhf )

      ! Saturation specific humidity at soil/snow surface temperature:
      call qsatfjh(ssnow%qstss,ssnow%tss-C%tfrz,met%pmb)

      IF(cable_user%ssnow_POTEV== "P-M") THEN
         
         !--- uses %ga from previous timestep    
         ssnow%potev =  Penman_Monteith(canopy%ga) 
      
      ELSE !by default assumes Humidity Deficit Method
      
         dq = ssnow%qstss - met%qvair
         ssnow%potev =  Humidity_deficit_method(dq,ssnow%qstss )
          
      ENDIF

         
      ! Soil latent heat:
      CALL latent_heat_flux()

      ! Soil sensible heat:
      canopy%fhs = air%rho*C%CAPP*(ssnow%tss - met%tvair) /ssnow%rtsoil
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

      CALL update_zetar()

   END DO           ! do iter = 1, NITER


   canopy%cduv = canopy%us * canopy%us / (max(met%ua,C%UMIN))**2

   !---diagnostic purposes
   canopy%gswx_T = rad%fvlai(:,1)/MAX( C%LAI_THRESH, canopy%vlaiw(:) )         & 
                   * canopy%gswx(:,1) + rad%fvlai(:,2) / MAX(C%LAI_THRESH,     &
                   canopy%vlaiw(:))*canopy%gswx(:,2)

    ! The surface conductance below is required by dust scheme; it is composed from canopy and soil conductances
    canopy%gswx_T = (1.-rad%transd)*max(1.e-06,canopy%gswx_T ) +  &   !contribution from  canopy conductance
                  rad%transd*(.01*ssnow%wb(:,1)/soil%sfc)**2 ! + soil conductance; this part is done as in Moses
    where ( soil%isoilm == 9 ) canopy%gswx_T = 1.e6   ! this is a value taken from Moses for ice points

    canopy%cdtq = canopy%cduv *( LOG( rough%zref_uv / rough%z0m) -              &
                 psim( canopy%zetar(:,NITER) * rough%zref_uv/rough%zref_tq )   &
               + psim( canopy%zetar(:,NITER) * rough%z0m/rough%zref_tq ) & ! new term from Ian Harman
                 ) / ( LOG( rough%zref_tq /(0.1*rough%z0m) )                   &
               - psis( canopy%zetar(:,NITER))                                  &
               + psis(canopy%zetar(:,NITER)*0.1*rough%z0m/rough%zref_tq) ) ! n


   ! Calculate screen temperature: 1) original method from SCAM
   ! screen temp., windspeed and relative humidity at 1.5m
   ! screen temp., windspeed and relative humidity at 2.0m
    tstar = - canopy%fh / ( air%rho*C%CAPP*canopy%us)
    qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us)
    zscrn = MAX(rough%z0m,2.0-rough%disp)
    ftemp = ( LOG(rough%zref_tq/zscrn)- psis(canopy%zetar(:,iterplus)) +       &
            psis(canopy%zetar(:,iterplus) * zscrn / rough%zref_tq) ) /C%VONK

    ! Calculate screen temperature:
    canopy%tscrn = met%tk - C%tfrz - tstar * ftemp

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
      
      IF ( canopy%vlaiw(j) > C%LAI_THRESH .and. rough%hruff(j) > 0.01) THEN
      
         IF ( rough%disp(j)  > 0.0 ) then
     
            term1(j) = EXP(2*C%CSW*canopy%rghlai(j)*(1-zscl(j)/rough%hruff(j)))
            term2(j) = EXP(2*C%CSW*canopy%rghlai(j) *                          &
                       (1-rough%disp(j)/rough%hruff(j)))
            term5(j) = MAX(2./3.*rough%hruff(j)/rough%disp(j), 1.)
         
         ENDIF
        
         term3(j) = C%A33**2*C%CTL*2*C%CSW*canopy%rghlai(j)

         IF( zscl(j) < rough%disp(j) ) THEN

            r_sc(j) = term5(j) * LOG(zscl(j)/rough%z0soilsn(j)) *              &
                      ( EXP(2*C%CSW*canopy%rghlai(j)) - term1(j) ) / term3(j)

         ELSEIF( rough%disp(j) <= zscl(j) .AND.                                &
                 zscl(j) < rough%hruff(j) ) THEN
                 
            r_sc(j) = rough%rt0us(j) + term5(j) * ( term2(j) - term1(j) ) /    &
                      term3(j)

         ELSEIF( rough%hruff(j) <= zscl(j) .AND.                               &
                 zscl(j) <  rough%zruffs(j) ) THEN

            r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + term5(j) *            &
                      ( zscl(j) - rough%hruff(j) ) /                           &
                      ( C%A33**2 * C%CTL * rough%hruff(j) )


         ELSEIF( zscl(j) >= rough%zruffs(j) ) THEN
            
            r_sc(j) = rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j) +     &
                      ( LOG( (zscl(j) - rough%disp(j)) /                       &
                      MAX( rough%zruffs(j)-rough%disp(j),                      &
                      rough%z0soilsn(j) ) ) - psis1( (zscl(j)-rough%disp(j))   &
                      / (rough%zref_tq(j)/canopy%zetar(j,iterplus) ) )         &
                      + psis1( (rough%zruffs(j) - rough%disp(j) )              &
                      / (rough%zref_tq(j)/canopy%zetar(j,iterplus ) ) ) )      &
                      / C%VONK

         ENDIF

        canopy%tscrn(j) = ssnow%tss(j) + (met%tk(j) - ssnow%tss(j)) *          &
                          MIN(1.,r_sc(j) / MAX( 1.,                            &
                          rough%rt0us(j) + rough%rt1usa(j) + rough%rt1usb(j)   &
                          + rt1usc(j))) - C%tfrz 
      ENDIF  

   ENDDO  
  
   CALL qsatfjh(rsts,canopy%tscrn,met%pmb)
     
   qtgnet = rsts * ssnow%wetfac - met%qv
   
   DO j=1,mp
      
      IF (qtgnet(j) .GT. 0. ) THEN
         qsurf(j) = rsts(j) * ssnow%wetfac(j)
      ELSE
         qsurf(j) = 0.1*rsts(j)*ssnow%wetfac(j) + 0.9*met%qv(j)
      ENDIF

      canopy%qscrn(j) = met%qv(j) - qstar(j) * ftemp(j)

      IF( canopy%vlaiw(j) >C%LAI_THRESH .and. rough%hruff(j) > 0.01)           &

            canopy%qscrn(j) = qsurf(j) + (met%qv(j) - qsurf(j)) * MIN( 1.,     &
                              r_sc(j) / MAX( 1., rough%rt0us(j) +              &
                              rough%rt1usa(j) + rough%rt1usb(j) + rt1usc(j) ) )

   ENDDO 


   ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
   canopy%dewmm = - (min(0.0,canopy%fevw) + min(0.0_r_2,canopy%fevc)) * &
        dels * 1.0e3 / (C%RHOW*air%rlam)

   ! Add dewfall to canopy water storage:
   canopy%cansto = canopy%cansto + canopy%dewmm
   
   ! Modify canopy water storage for evaporation:
   canopy%cansto = MAX(canopy%cansto-MAX(0.0,REAL(canopy%fevw))*dels &
     *1.0e3/(C%RHOW*air%rlam), 0.0)

   ! Calculate canopy water storage excess:
   canopy%spill=max(0.0, canopy%cansto-cansat)

   ! Move excess canopy water to throughfall:
   ! %through is /dels in UM app. (unpacked in hyd driver) for STASH output  
   canopy%through = canopy%through + canopy%spill
   
   ! Initialise 'throughfall to soil' as 'throughfall from canopy'; 
   ! snow may absorb
   canopy%precis = max(0.,canopy%through)

   ! Update canopy storage term:
   canopy%cansto=canopy%cansto - canopy%spill
   
   ! Calculate the total change in canopy water store (mm/dels):
   canopy%delwc = canopy%cansto-canopy%oldcansto
   
   ! calculate dgdtg, derivative of ghflux 3 instances
   ! d(canopy%fns)/d(ssnow%tgg)
   ! d(canopy%fhs)/d(ssnow%tgg)
   ! d(canopy%fes)/d(dq)
   ssnow%dfn_dtg = (-1.)*4.*C%EMSOIL*C%SBOLTZ*tss4/ssnow%tss  
   ssnow%dfh_dtg = air%rho*C%CAPP/ssnow%rtsoil      
   ssnow%dfe_ddq = ssnow%wetfac*air%rho*air%rlam*ssnow%cls/ssnow%rtsoil  
  
   ssnow%ddq_dtg = (C%rmh2o/C%rmair) /met%pmb * C%TETENA*C%TETENB * C%TETENC   &
                   / ( ( C%TETENC + ssnow%tss-C%tfrz )**2 )*EXP( C%TETENB *       &
                   ( ssnow%tss-C%tfrz ) / ( C%TETENC + ssnow%tss-C%tfrz ) )
   canopy%dgdtg = ssnow%dfn_dtg - ssnow%dfh_dtg - ssnow%dfe_ddq *    &
                  ssnow%ddq_dtg

   bal%drybal = REAL(ecy+hcy) - SUM(rad%rniso,2)                               &
                + C%CAPP*C%rmair*(tlfy-met%tk)*SUM(rad%gradis,2)  ! YP nov2009

   bal%wetbal = canopy%fevw + canopy%fhvw - SUM(rad%rniso,2) * canopy%fwet      &
                + C%CAPP*C%rmair * (tlfy-met%tk) * SUM(rad%gradis,2) *          &
                canopy%fwet  ! YP nov2009

   DEALLOCATE(cansat,gbhu)
   DEALLOCATE(dsx, fwsoil, tlfx, tlfy)
   DEALLOCATE(ecy, hcy, rny)
   DEALLOCATE(gbhf, csx)
   DEALLOCATE(ghwet)

CONTAINS

! ------------------------------------------------------------------------------

SUBROUTINE comp_friction_vel()
   USE cable_def_types_mod, only : mp
   REAL, DIMENSION(mp)  :: lower_limit, rescale

   psim_1 = psim(canopy%zetar(:,iter)) 
   
      rescale = C%VONK * MAX(met%ua,C%UMIN)
      z_eff = rough%zref_uv / rough%z0m
   
   psim_arg = canopy%zetar(:,iter) / z_eff 
   !---fix for compiler limitation. bitwise reproducable whilst we  
   !---we know it to 11th decimal. psim_arg typically of a few 
   !psim_arg = nint(psim_arg * 1.e11)*1.e-11

   psim_2 = psim( psim_arg )
          
   lower_limit = rescale / ( LOG(z_eff) - psim_1 + psim_2 )

   canopy%us = MAX(1.e-6, lower_limit )

END SUBROUTINE comp_friction_vel

! ------------------------------------------------------------------------------

FUNCTION Penman_Monteith( ground_H_flux ) RESULT(ssnowpotev)
   USE cable_def_types_mod, only : mp
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
   
   CALL qsatfjh(qsatfvar,met%tvair-C%tfrz,met%pmb)

   ssnowpotev = cc1 * (canopy%fns - ground_H_flux) + &
   cc2 * air%rho * air%rlam*(qsatfvar  - met%qvair)/ssnow%rtsoil
 
END FUNCTION Penman_Monteith


! ------------------------------------------------------------------------------
! method alternative to P-M formula above
FUNCTION humidity_deficit_method(dq,qstss ) RESULT(ssnowpotev)

   USE cable_def_types_mod, only : mp
   
   REAL, DIMENSION(mp) ::                                                      &
      ssnowpotev,    & ! 
      dq,            & ! sat spec hum diff.
      qstss             !dummy var for compilation
       
   INTEGER :: j
   
   DO j=1,mp
      !if(ssnow%snowd(j) > 1.0) dq(j) = max( -0.1e-3, dq(j))
      IF( ssnow%snowd(j)>1.0 .OR. ssnow%tgg(j,1).EQ.C%tfrz)                      &
         dq(j) = max( -0.1e-3, dq(j))
   ENDDO 
   
   ssnowpotev =air%rho * air%rlam * dq /ssnow%rtsoil
   
END FUNCTION Humidity_deficit_method

! ------------------------------------------------------------------------------
 
SUBROUTINE Latent_heat_flux() 

   USE cable_common_module
   USE cable_def_types_mod, only : mp

   REAL, DIMENSION(mp) ::                                                      &
      frescale,  flower_limit, fupper_limit

   INTEGER :: j
   
   ! Soil latent heat:
   canopy%fess= ssnow%wetfac * ssnow%potev
   WHERE (ssnow%potev < 0. ) canopy%fess = ssnow%potev
   
   ! Reduce soil evap due to presence of puddle
   pwet = max(0.,min(0.2,ssnow%pudsto/max(1.,ssnow%pudsmx)))
   canopy%fess = canopy%fess * (1.-pwet)

   frescale = soil%zse(1) * 1000. * air%rlam / dels         

   DO j=1,mp
      
      IF(ssnow%snowd(j) < 0.1 .AND. canopy%fess(j) .GT. 0. ) THEN

        IF (.not.cable_user%l_new_reduce_soilevp) THEN
         flower_limit(j) = REAL(ssnow%wb(j,1))-soil%swilt(j)/2.0
        ELSE
         ! E.Kowalczyk 2014 - reduces the soil evaporation
         flower_limit(j) = REAL(ssnow%wb(j,1))-soil%swilt(j)
        ENDIF
         fupper_limit(j) = MAX( 0._r_2,                                        &
                           flower_limit(j) * frescale(j)                       &
                           - ssnow%evapfbl(j,1)*air%rlam(j)/dels)

         canopy%fess(j) = MIN(canopy%fess(j), fupper_limit(j))
         
         fupper_limit(j) = REAL(ssnow%wb(j,1)-ssnow%wbice(j,1)) * frescale(j)

         canopy%fess(j) = min(canopy%fess(j), fupper_limit(j))

      END IF

      ssnow%cls(j)=1.
      
      IF (ssnow%snowd(j) >= 0.1 .and. ssnow%potev(j) > 0.) THEN

         ssnow%cls(j) = 1.1335
         canopy%fess(j) = MIN( (ssnow%wetfac(j)*ssnow%potev(j))*ssnow%cls(j), &
                          ssnow%snowd(j)/dels*air%rlam(j)*ssnow%cls(j))
      
      ENDIF

   ENDDO 
   
   ! Evaporation form soil puddle
   canopy%fesp = min(ssnow%pudsto/dels*air%rlam,max(pwet*ssnow%potev,0.))
   canopy%fes = canopy%fess + canopy%fesp

END SUBROUTINE latent_heat_flux

! -----------------------------------------------------------------------------

SUBROUTINE within_canopy( gbhu, gbhf )

   USE cable_def_types_mod, only : mp, r_2

   REAL(r_2), INTENT(IN), DIMENSION(:,:) ::                                    &
      gbhu,    &  ! forcedConvectionBndryLayerCond
      gbhf        ! freeConvectionBndryLayerCond
      
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
 
   INTEGER :: j
   
   rrbw = sum(gbhu+gbhf,2)/air%cmolar  ! MJT 
   
   ! leaf stomatal resistance for water
   rrsw = sum(canopy%gswx,2)/air%cmolar ! MJT
   
   DO j=1,mp
   
      IF(veg%meth(j) > 0 .AND. canopy%vlaiw(j) > C%LAI_THRESH .AND.              &
         rough%hruff(j) > rough%z0soilsn(j) ) THEN

         !   use the dispersion matrix (DM) to find the air temperature 
         !   and specific humidity 
         !   (Raupach, Finkele and Zhang 1997, pp 17)
         ! leaf boundary layer resistance for water
         ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmah(j) = (rt0(j)+rough%rt1(j))*((1.+air%epsi(j))*rrsw(j) + rrbw(j))  &
                   + air%epsi(j) * (rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))
         
         ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmbh(j) = (-air%rlam(j)/C%CAPP)*(rt0(j)*rough%rt1(j))*(rrbw(j)*rrsw(j))
         
         ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmch(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)*   &
                   (canopy%fhv(j) + canopy%fhs(j))/(air%rho(j)*C%CAPP)

         ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmae(j) = (-air%epsi(j)*C%CAPP/air%rlam(j))*(rt0(j)*rough%rt1(j)) *   &
                   (rrbw(j)*rrsw(j))

         ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmbe(j) = ( rt0(j) + ssnow%wetfac(j) * rough%rt1(j) ) *               &
                   ( (1.+air%epsi(j) ) * rrsw(j) + rrbw(j) ) +                 &
                   ( rt0(j) * rough%rt1(j) ) * ( rrbw(j) * rrsw(j) )

         ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmce(j) = ((1.+air%epsi(j))*rrsw(j) + rrbw(j))*rt0(j)*rough%rt1(j)*   &
                   (canopy%fev(j) + canopy%fes(j))/(air%rho(j)*air%rlam(j))
      
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
         CALL qsatfjh2(qstvair(j),met%tvair(j)-C%tfrz,met%pmb(j))
         
         ! Saturated vapour pressure deficit in canopy:
         met%dva(j) = ( qstvair(j) - met%qvair(j) ) *  C%rmair/C%RMH2o         &
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
      canopy%zetar(:,iterplus) = -( C%VONK * C%GRAV * rough%zref_tq *              &
                                 ( canopy%fh + 0.07 * canopy%fe ) ) /          &
                                 ( air%rho * C%CAPP * met%tk * canopy%us**3 )

      ! case NITER=2: final zetar=C%ZETmul*zetar(2) (compute only when iter=1)
      IF (NITER == 2) THEN
    
         canopy%zetar(:,2) = C%ZETmul * canopy%zetar(:,2)
    
         DO j=1,mp
            IF ( (met%fsd(j,1)+met%fsd(j,2))  ==  0.0 ) &
               canopy%zetar(j,2) = 0.5 * canopy%zetar(j,2)
         ENDDO

      END IF

      ! constrain zeta to C%ZETPOS and C%ZETNEG (set in param0)
      
      ! zetar too +
      canopy%zetar(:,iterplus) = MIN(C%ZETPOS,canopy%zetar(:,iterplus))        
      
      ! zetar too -
      canopy%zetar(:,iterplus) = MAX(C%ZETNEG,canopy%zetar(:,iterplus))        
    
   END IF ! (iter < NITER)
      
END SUBROUTINE update_zetar

! -----------------------------------------------------------------------------

FUNCTION qsatf(j,tair,pmb) RESULT(r)
  ! MRR, 1987
  ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
  ! HUMIDITY (KG/KG) FROM TETEN FORMULA
   
   REAL, INTENT(IN) ::                                                         &
      tair,         & ! air temperature (C)
      pmb             ! pressure PMB (mb)
  
   INTEGER, INTENT(IN) :: j 
  
   REAL           :: r    ! result; sat sp humidity
  
   r = (C%RMH2o/C%rmair) * (C%TETENA*EXP(C%TETENB*tair/(C%TETENC+tair))) / pmb
END FUNCTION qsatf

! -----------------------------------------------------------------------------

SUBROUTINE qsatfjh(var,tair,pmb) 
   USE cable_def_types_mod, only : mp
   REAL, INTENT(IN), DIMENSION(mp) ::                                          &
      tair,                        & ! air temperature (C)
      pmb                            ! pressure PMB (mb)
  
   REAL, INTENT(OUT), DIMENSION(mp) ::                                         &
      var                            ! result; sat sp humidity
  
  INTEGER :: j
  
  DO j=1,mp

     var(j) = (C%RMH2o/C%rmair) * (C%TETENA*EXP(C%TETENB*tair(j)/(C%TETENC+tair(j))))    &
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
      
      var = (C%RMH2o/C%rmair) * (C%TETENA*EXP(C%TETENB*tair/(C%TETENC+tair))) / pmb

END SUBROUTINE qsatfjh2

! -----------------------------------------------------------------------------

FUNCTION psim(zeta) RESULT(r)
   USE cable_def_types_mod, only : mp
   ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
   ! computes integrated stability function psim(z/l) (z/l=zeta)
   ! for momentum, using the businger-dyer form for unstable cases
   ! and the Beljaars and Holtslag (1991) form for stable cases.
   
   
   REAL, INTENT(IN), DIMENSION(mp) ::  zeta       !

   ! function result   
   REAL, DIMENSION(mp) :: r
   
   REAL, DIMENSION(mp) ::                                                      &
      x,       & !
      z,       & !
      stable,  & !
      unstable   !
      
   REAL, PARAMETER ::                                                          &
      gu = 16.0,  & !
      gs = 5.0
  
   REAL, PARAMETER ::                                                          &
      a = 1.0,       & !
      b = 0.667,     & !
      xc = 5.0,       & !
      d = 0.35         !

   z = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable
   
   ! Beljaars and Holtslag (1991) for stable
   stable = -a*zeta - b*(zeta - xc/d)*exp( -d*zeta) - b*xc/d
   x      = (1.0 + gu*abs(zeta))**0.25
   unstable = ALOG((1.0+x*x)*(1.0+x)**2/8) - 2.0*atan(x) + C%PI_C*0.5
   r = z*stable + (1.0-z)*unstable
END FUNCTION psim

! -----------------------------------------------------------------------------

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
   
   z      = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable 

   ! Beljaars and Holtslag (1991) for stable
   stzeta = MAX(0.,zeta)
   stable = -(1.+2./3.*a*stzeta)**(3./2.) -  &
             b*(stzeta-c/d)*exp(-d*stzeta) - b*c/d + 1.
   y      = (1.0 + gu*abs(zeta))**0.5
   unstable = 2.0 * alog((1+y)*0.5)
   r   = z*stable + (1.0-z)*unstable

END FUNCTION psis

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

! -----------------------------------------------------------------------------

ELEMENTAL FUNCTION rplant(rpconst, rpcoef, tair) result(z)
   REAL, INTENT(IN)     :: rpconst
   REAL, INTENT(IN)     :: rpcoef
   REAL, INTENT(IN)     :: tair
   REAL                 :: z
   z = rpconst * exp(rpcoef * tair)
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
      tlfy,          & ! leaf temp (K) - assC%UMINg the temperature of 
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

      IF(canopy%vlaiw(j) > C%LAI_THRESH) THEN

         ! VEG SENSIBLE & LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
         ! calculate total thermal resistance, rthv in s/m
         ghwet(j) = 2.0   * sum_gbh(j) 
         gwwet(j) = 1.075 * sum_gbh(j) 
         ghrwet(j) = sum_rad_gradis(j) + ghwet(j)
         
         ! Calculate fraction of canopy which is wet:
         canopy%fwet(j) = MAX( 0.0, MIN( 1.0,                                  &
                          0.8 * canopy%cansto(j) / MAX( cansat(j), 0.01 ) ) )
         
         ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
         ! to avoid excessive evaporation:
         ccfevw(j) = MIN(canopy%cansto(j) * air%rlam(j) / dels, &
                       2.0 / (1440.0 / (dels/60.0)) * air%rlam(j) )
   
         canopy%fevw(j) = MIN( canopy%fwet(j) * ( air%dsatdk(j) *              &
                         ( sum_rad_rniso(j)- C%CAPP*C%rmair*( met%tvair(j)     &
                         - met%tk(j) ) * sum_rad_gradis(j) )                   &
                         + C%CAPP * C%rmair * met%dva(j) * ghrwet(j) )         &
                         / ( air%dsatdk(j)+air%psyc(j)*ghrwet(j) / gwwet(j) )  &
                         , ccfevw(j) )

         canopy%fevw_pot(j) = ( air%dsatdk(j)* (sum_rad_rniso(j) -             &
                              C%CAPP * C%rmair * ( met%tvair(j) - met%tk(j) )  &
                              *sum_rad_gradis(j) )                             &
                              + C%CAPP * C%rmair * met%dva(j) * ghrwet(j))     &
                              / (air%dsatdk(j)+air%psyc(j)*ghrwet(j)/gwwet(j) )
          
         ! calculate sens heat from wet canopy:
         canopy%fhvw(j) = canopy%fwet(j) * ( sum_rad_rniso(j) -C%CAPP * C%rmair&
                          * ( tlfy(j) - met%tk(j) ) * sum_rad_gradis(j) )      &
                           - canopy%fevw(j)

      ENDIF
       
   ENDDO 
           
END SUBROUTINE wetLeaf

! -----------------------------------------------------------------------------

END SUBROUTINE define_canopy

! -----------------------------------------------------------------------------

SUBROUTINE Surf_wetness_fact( cansat, canopy, ssnow,veg, met, soil, dels )

   USE cable_common_module
   USE cable_def_types_mod

   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
   TYPE (soil_snow_type), intent(inout):: ssnow
   TYPE (soil_parameter_type), intent(inout)   :: soil
   TYPE (canopy_type), INTENT(INOUT)   :: canopy
   TYPE (met_type), INTENT(INOUT)   :: met
   
   REAL, INTENT(IN) :: dels ! integration time setp (s)

   REAL,INTENT(IN), DIMENSION(:) :: cansat ! max canopy intercept. (mm)

   !local variables
   REAL, DIMENSION(mp)  :: lower_limit, upper_limit,ftemp
   
   INTEGER :: j

   ! Rainfall variable is limited so canopy interception is limited,
   ! used to stabilise latent fluxes.
   ! to avoid excessive direct canopy evaporation (EK nov2007, snow scheme)
   upper_limit = 4.0 * MIN(dels,1800.0) / (60.0 * 1440.0 ) 
   ftemp =MIN(met%precip-met%precip_sn, upper_limit )
   ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
   lower_limit = cansat - canopy%cansto
   upper_limit = max(lower_limit, 0.0) 
   canopy%wcint = MERGE( MIN( upper_limit, ftemp ), 0.0,                       &
                  ftemp > 0.0  .AND. met%tk > C%tfrz)  !EAK, 09/10

   ! Define canopy throughfall (100% of precip if temp < 0C, see above):
   canopy%through = met%precip_sn + MIN( met%precip - met%precip_sn ,          &
                    MAX( 0.0, met%precip - met%precip_sn - canopy%wcint) ) 

   ! Add canopy interception to canopy storage term:
   canopy%cansto = canopy%cansto + canopy%wcint

   ! Calculate fraction of canopy which is wet:
   canopy%fwet   = MAX( 0.0, MIN( 0.9, 0.8 * canopy%cansto /                   &
                   MAX( cansat, 0.01 ) ) )

   ssnow%wetfac = MAX( 1.e-6, MIN( 1.0,                                        &
                  ( REAL (ssnow%wb(:,1) ) - soil%swilt/ 2.0 )                  &
                  / ( soil%sfc - soil%swilt/2.0 ) ) )
  
   DO j=1,mp
   
      IF( ssnow%wbice(j,1) > 0. )                                              &
         ssnow%wetfac(j) = ssnow%wetfac(j) * MAX( 0.5, 1. - MIN( 0.2,          &
                           ( ssnow%wbice(j,1) / ssnow%wb(j,1) )**2 ) )

      IF( ssnow%snowd(j) > 0.1) ssnow%wetfac(j) = 0.9
      
      IF ( veg%iveg(j) == 16 .and. met%tk(j) >= C%tfrz + 5. )                  &
         ssnow%wetfac(j) = 1.0 ! lakes: hard-wired number to be removed
      
      IF( veg%iveg(j) == 16 .and. met%tk(j) < C%tfrz + 5. )                    &
         ssnow%wetfac(j) = 0.7 ! lakes: hard-wired number to be removed

   ENDDO 
      

   ! owetfac introduced to reduce sharp changes in dry regions,
   ! especially in offline runs in which there may be discrepancies b/n
   ! timing of precip and temperature change (EAK apr2009)
   ssnow%wetfac = 0.5*(ssnow%wetfac + ssnow%owetfac)

END SUBROUTINE Surf_wetness_fact

! -----------------------------------------------------------------------------
! dryleaf subr was here
! -----------------------------------------------------------------------------
!SUBROUTINE photosynthesis( csxz, cx1z, cx2z, gswminz,                          &
! ------------------------------------------------------------------------------

FUNCTION ej3x(parx,x) RESULT(z)
   
   REAL, INTENT(IN)     :: parx
   REAL, INTENT(IN)     :: x
   REAL                 :: z
   
   z = MAX(0.0,                                                                &
       0.25*((C%alpha3*parx+x-sqrt((C%alpha3*parx+x)**2 -                      &
       4.0*C%convx3*C%alpha3*parx*x)) /(2.0*C%convx3)) )

END FUNCTION ej3x

! ------------------------------------------------------------------------------

FUNCTION ej4x(parx,x) RESULT(z)
   
   REAL, INTENT(IN)     :: parx
   REAL, INTENT(IN)     :: x
   REAL                 :: z
 
   z = MAX(0.0,                                                                &
        (C%alpha4*parx+x-sqrt((C%alpha4*parx+x)**2 -                           &
        4.0*C%convx4*C%alpha4*parx*x))/(2.0*C%convx4))
 
END FUNCTION ej4x

! ------------------------------------------------------------------------------
SUBROUTINE fwsoil_calc_std(fwsoil, soil, ssnow, veg) 
   USE cable_def_types_mod
   TYPE (soil_snow_type), INTENT(INOUT):: ssnow
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
   REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
   REAL, DIMENSION(mp) :: rwater ! soil water availability

   rwater = MAX(1.0e-9,                                                    &
            SUM(veg%froot * MAX(1.0e-9,MIN(1.0_r_2,ssnow%wb -                   &
            SPREAD(soil%swilt, 2, ms))),2) /(soil%sfc-soil%swilt))
  
   fwsoil = MAX(1.0e-9,MIN(1.0, veg%vbeta * rwater))
      
END SUBROUTINE fwsoil_calc_std 

! ------------------------------------------------------------------------------

SUBROUTINE fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg) 
   USE cable_def_types_mod
   TYPE (soil_snow_type), INTENT(INOUT):: ssnow
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
   REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
   REAL, DIMENSION(mp) :: rwater ! soil water availability
   REAL, DIMENSION(mp,3)          :: xi, ti, si
   INTEGER :: j

   rwater = MAX(1.0e-9,                                                    &
            SUM(veg%froot * MAX(0.0,MIN(1.0_r_2,ssnow%wb -                   &
            SPREAD(soil%swilt, 2, ms))),2) /(soil%sfc-soil%swilt))

   fwsoil = 1.

   rwater = soil%swilt + rwater * (soil%sfc-soil%swilt)
   
   xi(:,1) = soil%swilt
   xi(:,2) = soil%swilt + (soil%sfc - soil%swilt)/2.0
   xi(:,3) = soil%sfc
   
   ti(:,1) = 0.
   ti(:,2) = 0.9
   ti(:,3) = 1.0
   
   si(:,1) = (rwater - xi(:,2)) / ( xi(:,1) - xi(:,2)) *                       &
             (rwater - xi(:,3)) / ( xi(:,1) - xi(:,3))

   si(:,2) = (rwater - xi(:,1)) / ( xi(:,2) - xi(:,1)) *                       &
             (rwater - xi(:,3)) / ( xi(:,2) - xi(:,3))

   si(:,3) = (rwater - xi(:,1)) / ( xi(:,3) - xi(:,1)) *                       &
             (rwater - xi(:,2)) / ( xi(:,3) - xi(:,2))

   DO j=1,mp
      IF (rwater(j) < soil%sfc(j) - 0.02)                                      &
         fwsoil(j) = max(0.,min(1., ti(j,1)*si(j,1) +                          &
                       ti(j,2)*si(j,2) + ti(j,3)*si(j,3)))

   ENDDO

END SUBROUTINE fwsoil_calc_non_linear 

! ------------------------------------------------------------------------------

! ypw 19/may/2010 soil water uptake efficiency (see Lai and Ktaul 2000)
SUBROUTINE fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssnow, veg) 
   USE cable_def_types_mod
   TYPE (soil_snow_type), INTENT(INOUT):: ssnow
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
   REAL, INTENT(OUT), DIMENSION(:):: fwsoil ! soil water modifier of stom. cond
   INTEGER   :: ns
   REAL, parameter ::rootgamma = 0.01   ! (19may2010)
   REAL, DIMENSION(mp)  :: dummy, normFac
   !--- local level dependent rwater 
   REAL, DIMENSION(mp,ms)  :: frwater

   fwsoil(:) = 0.0
   normFac(:) = 0.0

   DO ns=1,ms
     
      dummy(:) = rootgamma/max(1.0e-3,ssnow%wb(:,ns)-soil%swilt(:))

      frwater(:,ns) = MAX(1.0e-4,((ssnow%wb(:,ns)-soil%swilt(:))/soil%ssat(:)) &
                      ** dummy)
      
      fwsoil(:) = min(1.0,max(fwsoil(:),frwater(:,ns)))
      
      normFac(:) = normFac(:) + frwater(:,ns) * veg%froot(:,ns)

   ENDDO

END SUBROUTINE fwsoil_calc_Lai_Ktaul


    
END MODULE cable_canopy_module
