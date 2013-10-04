
MODULE soil_snow_module
   
   USE physical_constants
   USE define_dimensions       
   USE define_types       
   USE io_variables, ONLY: landpt

   PRIVATE

   REAL, PARAMETER ::                                                          &
      cgsnow = 2090.0,     & ! specific heat capacity for snow
      csice = 2.100e3,     & ! specific heat capacity for ice
      cswat = 4.218e3,     & ! specific heat capacity for water
      cp = capp,           & ! specific heat capacity for air
      rhowat = 1000.0,     & ! density of water
      snmin = 1.,          & ! for 3-layer;
      max_ssdn = 750.0,    & !
      max_sconds = 2.51,   & !
      frozen_limit = 0.85    ! EAK Feb2011 (could be 0.95)
   
   !jhan:make parameter
   REAL :: max_glacier_snowd
 
   ! This module contains the following subroutines:
   PUBLIC soil_snow ! must be available outside this module
   PRIVATE snowdensity, snow_melting, snowcheck, snowl_adjust 
   PRIVATE trimb, smoisturev, snow_accum, stempv
   PRIVATE soilfreeze, remove_trans

CONTAINS

! SUBROUTINE trimb
!
!      this routine solves the system
!	   a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kmax-1
!	   with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)	       for k=1
!	   and	 a(k)*u(k-1)+b(k)*u(k)=rhs(k)	       for k=kmax
!
!	 the Thomas algorithm is used for solving sets of linear equation
!	 rhs initially contains rhs; leaves with answer (jlm)
!	 n.b. this one does not assume b = 1-a-c
!
SUBROUTINE trimb (a, b, c, rhs, kmax)

   INTEGER, INTENT(IN)                  :: kmax ! no. of discrete layers    

   REAL(r_2), DIMENSION(:,:), INTENT(IN) ::                                    &
      a,    & ! coef "A" in finite diff eq
      b,    & ! coef "B" in finite diff eq
      c       ! coef "C" in finite diff eq
   
   REAL(r_2), DIMENSION(:,:), INTENT(INOUT)  :: rhs ! right hand side of eq
  
   REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) ::                                &
      e, temp, g 

   INTEGER :: k   ! do lloop counter
  

   e(:,1) = c(:,1) / b(:,1)
   DO k = 2, kmax - 1
     temp(:,k) = 1. / ( b(:,k) - a(:,k) * e(:,k-1) )
     e(:,k) = c(:,k) * temp(:,k)
   END DO

   g(:,1) = rhs(:,1) / b(:,1)
   DO k = 2, kmax - 1
     g(:,k) = ( rhs(:,k) - a(:,k) * g(:,k-1) ) * temp(:,k)
   END DO

   ! do back substitution to give answer now
   rhs(:,kmax) = ( rhs(:,kmax) - a(:,kmax) * g(:,kmax-1) )                   &
                 / ( b(:,kmax) - a(:,kmax) * e(:,kmax-1) )
                 
   DO k = kmax - 1, 1, - 1
     rhs(:,k) = g(:,k) - e(:,k) * rhs(:,k + 1)
   END DO
  
END SUBROUTINE trimb

! SUBROUTINE smoisturev (fwtop,dels,ssoil,soil)
!      Solves implicit soil moisture equation
!      Science development by Eva Kowalczyk and John McGregor, CMAR
!
SUBROUTINE smoisturev (dels,ssoil,soil,veg)
   
   USE cable_common_module
   
   REAL, INTENT(IN) :: dels    ! time step size (s)
   
   TYPE(soil_snow_type),      INTENT(INOUT) ::                                &
      ssoil ! soil and snow variables
 
   TYPE(soil_parameter_type), INTENT(INOUT) ::                                &
      soil  ! soil parameters
   
   TYPE(veg_parameter_type), INTENT(IN)  :: veg

   ! nmeth selects the solution method
   ! Values as follows:
   !  -1 for simple implicit D  - PREFERRED METHOD
   !   1 for fully implicit solution
   !   2 for simpler implicit
   !   3 for simple implicit D, explicit K 
   !   4 for simple implicit D, implicit K
   !   0 for simple implicit D, new jlm TVD K
   INTEGER, PARAMETER ::                                                       &
      nmeth = -1 

   REAL, DIMENSION(mp) ::                                                      &
      totwba,     & ! diagnostic
      totwbb,     & !
      totwbc,     & !
      totwblb,    & !
      totwblc,    & !
      wbficemx,   & ! 
      wblfmn,     & !
      wblfmx        !
 
   REAL(r_2), DIMENSION(mp) ::                                                 &
      fact,       & ! 
      fact2,      & !
      fluxhi,     & !
      fluxlo,     & !
      hydss,      & ! hydraulic conductivity adjusted for ice
      phi,        & !
      pwb,        & !
      rat,        & !
      speed_k,    & !
      ssatcurr_k, & !
      wbh_k,      & !
      wbl_k,      & !
      wbl_kp,     & !
      wh,         & !
      z3_k,       & !
      pwb_wbh       ! 
   
   
   REAL, DIMENSION(mp,ms+1) :: z1mult

   ! change dimension of at,bt,ct from 3*ms to ms (BP Jun2010)
   REAL(r_2), DIMENSION(mp,ms) ::                                              &
      at,      & ! coef "A" in finite diff eq
      bt,      & ! coef "B" in finite diff eq
      ct,      & ! coef "C" in finite diff eq
      ssatcurr   !
   
   REAL(r_2), DIMENSION(mp,ms+1) ::                                            &
      wbh,     & !
      z1,      & !
      z2,      & !
      z3         !
   
   REAL(r_2), DIMENSION(mp,0:ms) ::                                            &
      fluxh,      & !
      delt,       & !
      dtt           !
   
   LOGICAL :: is_open     ! Is file open?
   
   INTEGER ::                                                                  &
      u,    & ! I/O unit
      k
   


   at = 0.0
   bt = 1.0
   ct = 0.0
   z1mult(:,1) = 0.0       ! corresponds to 2b+3
   z1mult(:,ms+1) = 0.0    ! corresponds to 2b+3
   z1(:,1) = 0.0           ! i.e. K(.5),    value at surface
   z1(:,ms+1) = 0.0        ! i.e. K(ms+.5), value at bottom

   ! nmeth: equation solution technique
   IF (nmeth <= 0) THEN

      ! jlm split TVD version
      ! all land points
      delt(:,0) = 0.0
      fluxh(:,0) = 0.0
      fluxh(:,ms) = 0.0
      
      DO k = 1, ms-1
      
         ! Calculate amount of liquid soil water:
         wbl_k = MAX( 0.01_r_2, ssoil%wb(:,k) - ssoil%wbice(:,k) )
         wbl_kp = MAX( 0.01_r_2, ssoil%wb(:,k+1) - ssoil%wbice(:,k+1) )
         
         ! Calculate difference in liq soil water b/w consecutive layers:
         delt(:,k) = wbl_kp - wbl_k
         
         ! especially to allow for isolated frozen layers, use min speed
         wh = MIN( wbl_k, wbl_kp )
         WHERE( ssoil%wbice(:,k) > 0.05 .OR. ssoil%wbice(:,k+1) > 0.01 )       &
            wh = 0.9*wbl_k + 0.1*wbl_kp
         
         ! with 50% wbice, reduce hyds by 1.e-5
         ! Calculate hyd conductivity adjusted for ice:
         hydss = soil%hyds

         speed_k = hydss * (wh / soil%ssat )**( soil%i2bp3 - 1 )
         
         ! update wb by TVD method
         rat = delt(:,k - 1) / ( delt(:,k) + SIGN( REAL( 1.0e-20, r_2 ),       &
               delt(:,k) ) )

         phi = MAX( 0.0_r_2, MIN( 1.0_r_2, 2.0_r_2 * rat ),                    &
               MIN( 2.0_r_2, rat ) ) ! 0 for -ve rat
         fluxhi = wh
         fluxlo = wbl_k

         ! scale speed to grid lengths per dels & limit speed for stability
         ! 1. OK too for stability
         speed_k = MIN( speed_k, REAL( 0.5 * soil%zse(k) / dels , r_2 ) )
         fluxh(:,k) = speed_k * ( fluxlo + phi * ( fluxhi - fluxlo ) )
      
      END DO
      
      ! calculate drainage (this code replaces the code in the surfb)
      k = ms 

      WHERE( ssoil%wb(:,ms) > soil%sfc(:) )

         wbl_k = MAX( 0.001_r_2, ssoil%wb(:,ms) - ssoil%wbice(:,ms) )
         wbl_kp = MAX( 0.001_r_2, soil%ssat(:) - ssoil%wbice(:,ms) )
         
         wh = MIN( wbl_k, wbl_kp )
         
         WHERE( ssoil%wbice(:,ms) .GT. 0.05 ) wh = 0.9 * wbl_k + 0.1 * wbl_kp
       
         ! Calculate hyd conductivity adjusted for ice:
         hydss = soil%hyds
    
         speed_k = hydss * ( wh / soil%ssat )**( soil%i2bp3 - 1 )
         speed_k =  0.5 * speed_k / ( 1. - MIN( 0.5, 10. * ssoil%wbice(:,ms) ) )
         fluxlo = wbl_k
         
         ! scale speed to grid lengths per dt & limit speed for stability
         speed_k = MIN( 0.5 * speed_k, 0.5 * soil%zse(ms) / dels )
         fluxh(:,ms) = MAX( 0.0, speed_k * fluxlo )
     
      END WHERE

      ! update wb by TVD method
      DO k = ms, 1, -1
        
         IF(  nmeth == -1 ) THEN ! each new wb constrained by ssat
            fluxh(:,k-1) = MIN( fluxh (:,k-1), ( soil%ssat - ssoil%wb(:,k) )   &
                           * soil%zse(k) / dels + fluxh(:,k) )
         END IF

         ! fluxh (:,ms) is drainage
         ssoil%wb(:,k) = ssoil%wb(:,k) + dels * ( fluxh(:,k-1) - fluxh(:,k) )  &
                         / soil%zse(k)

         ! re-calculate wblf
         ssatcurr_k = soil%ssat - ssoil%wbice(:,k)
         dtt(:,k) = dels / ( soil%zse(k) * ssatcurr_k )

         ! this defn of wblf has different meaning from previous one in surfbv
         ! N.B. are imposing wbice<wb, so wblf <1
         ssoil%wblf(:,k) = ( ssoil%wb(:,k) - ssoil%wbice(:,k) ) / ssatcurr_k
      
      END DO
      
      ssoil%rnof2 = dels * REAL( fluxh(:,ms), r_1 ) * 1000.0

      ! wbh_k represents wblf(k-.5)
      DO k = 2, ms
         
         ssatcurr_k = REAL( soil%ssat, r_2 ) - ssoil%wbice(:,k)
         wbh_k = ( soil%zse(k) * ssoil%wblf(:,k-1) + soil%zse(k-1)             &
                 * ssoil%wblf(:,k) ) / ( soil%zse(k) + soil%zse(k-1) )
         ! i.e. wbh**(bch+1)
         fact = wbh_k**( soil%ibp2 - 1 )

         ! with 50% wbice, reduce hbsh by 1.e-5
         pwb_wbh = ( soil%hsbh * ( 1. - MIN( 2. * MIN (0.1_r_2, MAX(           &
                  ssoil%wbice(:,k-1) / MAX( 0.01_r_2, ssoil%wb(:,k-1) ),       &
                  ssoil%wbice(:,k)   / MAX( 0.01_r_2, ssoil%wb(:,k) ) ) )      &
                  , 0.1_r_2 ) ) )                                              &
                  * MAX( soil%pwb_min, wbh_k * fact )

         ! moisture diffusivity (D) is  wbh*pwb; hsbh includes b
         ! i.e. D(k-.5)/soil%zshh(k)
         z3_k = pwb_wbh / soil%zshh (k)

         ! where dtt=dels/(soil%zse(k)*ssatcurr_k)
         at (:,k) = - dtt(:,k) * z3_k
         ct (:,k-1) = - dtt(:,k-1) * z3_k
     
      END DO

      bt = 1. - at - ct
      ssoil%wblf(:,1) = ssoil%wblf(:,1) + dtt(:,1) * ssoil%fwtop1 / rhowat
      ssoil%wblf(:,2) = ssoil%wblf(:,2) + dtt(:,2) * ssoil%fwtop2 / rhowat
      ssoil%wblf(:,3) = ssoil%wblf(:,3) + dtt(:,3) * ssoil%fwtop3 / rhowat
    
   END IF
   ! END: IF (nmeth <= 0) THEN

   IF ( nmeth > 0 ) THEN
      
      wbficemx = 0.0
      
      DO k = 1, ms
         
         ssatcurr(:,k) = REAL(soil%ssat,r_2) - ssoil%wbice(:,k)
         
         ! this defn of wblf has different meaning from previous one in surfbv
         ! N.B. are imposing wbice<wb, so wblf <1
         ssoil%wblf(:,k) = ( ssoil%wb(:,k) - ssoil%wbice(:,k) ) / ssatcurr(:,k)
         
         ssoil%wbfice(:,k) = REAL( ssoil%wbice(:,k) ) / soil%ssat
         
         wbficemx = MAX( wbficemx, ssoil%wbfice(:,k) )
         dtt(:,k) = dels / ( soil%zse(k) * ssatcurr(:,k) )
      
      END DO

      IF( nmeth == 1 ) THEN ! full implicit method
        
         DO k = 2, ms
           
            wbh(:,k) = ( soil%zse(k) * ssoil%wblf(:,k-1) + soil%zse(k-1)       &
                       * ssoil%wblf(:,k) ) / ( soil%zse(k) + soil%zse(k-1) )
                       
            fact = wbh(:,k)**( soil%ibp2 - 1 ) ! i.e. wbh**(bch+1)
            fact2 = fact * fact
            pwb = soil%hsbh * fact
           
            ! moisture diffusivity (D) is  wbh*pwb
            ! other term (K) is wbh*soil%hyds*fact2
            z1(:,k) = wbh(:,k) * ( (soil%i2bp3 - 1 ) * soil%hyds * fact2       &
                      - soil%ibp2 * pwb *                                      &
                      ( ssoil%wblf(:,k) - ssoil%wblf(:,k-1 ) ) / soil%zshh (k) )
           
            z2(:,k) = - soil%i2bp3 * soil%hyds * fact2 + soil%ibp2 * pwb       &
                      * ( ssoil%wblf(:,k) - ssoil%wblf(:,k-1) ) / soil%zshh (k)
           
            z3(:,k) = pwb * wbh(:,k) / soil%zshh (k)
            
            at(:,k) = dtt(:,k) * (z2(:,k) * 0.5 * soil%zse(k) / soil%zshh (k)  &
                      - z3(:,k) )
        
         END DO
        
         DO k = 1, ms - 1
           
            ct(:,k) = dtt(:,k) * ( - z2(:,k+1) * 0.5 * soil%zse(k)             &
                      / soil%zshh (k+1) - z3(:,k+1) )

            bt(:,k) = 1.0 + dtt(:,k) * ( - z2(:,k+1) * 0.5 * soil%zse(k+1)     &
                      / soil%zshh (k+1) + z2(:,k) * 0.5 * soil%zse( MAX( k-1,  &
                      1 ) ) / soil%zshh (k) + z3(:,k+1) + z3(:,k) )

         END DO
         
         bt(:,ms) = 1.0 + dtt(:,ms) * ( z2(:,ms) * 0.5 * soil%zse(ms)          &
                    / soil%zshh (ms) + z3(:,ms) )

         DO k = 1, ms
            ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k) *                     &
                              ( z1(:,k+1) - z1(:,k) )
         END DO
     
      END IF ! (nmeth == 1)
     
      IF (nmeth >= 2) THEN ! part implicit method
         
         DO k = 2, ms
            z1mult(:,k) = soil%i2bp3 ! corresponds to 2b+3
         END DO

         DO k = 2, ms ! wbh(k) represents wblf(k-.5)
            wbh(:,k) = ( soil%zse(k) * ssoil%wblf(:,k-1) + soil%zse(k-1)       &
                       * ssoil%wblf(:,k) ) / ( soil%zse(k) + soil%zse(k-1) )
            
            fact = wbh(:,k)**( soil%ibp2 - 1 ) ! i.e. wbh**(bch+1)
            
            IF (nmeth == 2) pwb_wbh = soil%hsbh * wbh(:,k) * fact
            
            IF (nmeth >= 3)                                                    &
               pwb_wbh = soil%hsbh * MAX( soil%pwb_min,wbh(:,k) * fact)
               
            ! moisture diffusivity (D) is  wbh*pwb
            ! other term (K) is wbh*soil%hyds*fact2
            z1(:,k) = soil%hyds * fact2 !  i.e. K(k-.5)/wbh(:,k)
            z3(:,k) = pwb_wbh / soil%zshh(k) !  i.e. D(k-.5)/soil%zshh(k)
            at(:,k) = - dtt(:,k) * z3(:,k)
            ct(:,k-1) = - dtt(:,k-1) * z3(:,k)
         
         END DO
         
         bt = 1. - at - ct
        
         IF (nmeth == 4) THEN ! for simple implicit D, implicit K
            bt(:,1) = bt(:,1) + dtt(:,1) * z1mult(:,1+1) &
                * z1(:,1+1) * soil%zse(1+1) / (soil%zse(1) + soil%zse(1+1) )
           DO k = 2, ms
              at(:,k)   = at(:,k) - dtt(:,k) * z1mult(:,k) * z1(:,k)           &
                          * soil%zse(k) / (soil%zse(k) + soil%zse(k-1) )
             
              ct(:,k-1) = ct(:,k-1) + dtt(:,k-1) * z1mult(:,k) * z1(:,k)       &
                          * soil%zse(k-1) / ( soil%zse(k) + soil%zse(k-1) )

              bt(:,k)   = bt(:,k) - dtt(:,k) * z1mult(:,k) * z1(:,k)           &
                          * soil%zse(k-1) / ( soil%zse(k) + soil%zse(k-1) )    &
                          + dtt(:,k) * z1mult(:,k+1) * z1(:,k+1)               &
                          * soil%zse(k+1) / ( soil%zse(k) + soil%zse(k+1) )
           
            END DO

         END IF ! (nmeth == 4)
      
         DO k = 2, ms
            ! i.e. now K(k-.5)
            z1(:,k) = wbh(:,k) * z1(:,k)
         END DO

         ! the following top & bottom b.c.'s will preserve a uniform column
         !     z1(1) =z1(2)   ! simple dk/dz=0
         !     z1(ms+1)=z1(ms) ! simple dk/dz=0
         ! N.B. z1 are here +ve
         z1(:,1) = MIN( z1(:,2), z1(:,ms) )
         z1(:,ms + 1) = z1(:,1)
         
         ! no gravit. term if too much ice 11/12/00
         DO k = 1, ms
            
            IF (nmeth == 4) THEN
               
               WHERE( wbficemx < 0.75 )                                        &
                  ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k)                 &
                                    * ( ( z1mult(:,k+1) - 1.0 ) * z1(:,k+1)    &
                                    - (z1mult(:,k) - 1.0) * z1(:,k) )

            ELSE
               
               WHERE( wbficemx < 0.75 )                                        &
                  ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k)                 &
                                    * ( z1(:,k) - z1(:,k+1) )
              
            END IF

         END DO
         
      END IF

      ! Block below for testing purposes only:
      IF( ntest > 0 ) THEN
         wblfmx = MAXVAL( REAL( ssoil%wblf ), 2 )
         wblfmn = MINVAL( REAL( ssoil%wblf ), 2 )
        
         totwbb = 0.0
         totwblb = 0.0
         DO k = 1, ms
            totwbb = totwbb + soil%zse(k) * REAL(ssoil%wb(:,k),r_1)
            totwblb = totwblb + soil%zse(k) * REAL(ssoil%wblf(:,k),r_1)
         END DO
      
      END IF

     IF (nmeth == 3) THEN
         ! artificial fix applied here for safety (explicit nmeth only)
         DO k = 1, ms
            ssoil%wblf(:,k) = MAX( 0.0_r_2, MIN( ssoil%wblf(:,k), 1.0_r_2 ) )
         END DO
      END IF
   
      ssoil%wblf(:,1) = ssoil%wblf(:,1) + dtt(:,1) * ssoil%fwtop1 / rhowat
      ssoil%wblf(:,2) = ssoil%wblf(:,2) + dtt(:,2) * ssoil%fwtop2 / rhowat
      ssoil%wblf(:,3) = ssoil%wblf(:,3) + dtt(:,3) * ssoil%fwtop3 / rhowat
    
   END IF  ! IF (nmeth > 0)

   CALL trimb(at, bt, ct, ssoil%wblf, ms)
   
   DO k = 1, ms
      ssatcurr(:,k) = soil%ssat - ssoil%wbice(:,k)
      ssoil%wb(:,k) = ssoil%wblf(:,k) * ssatcurr(:,k) + ssoil%wbice(:,k)
      ssoil%wbice(:,k) = MIN( ssoil%wbice(:,k), frozen_limit * ssoil%wb(:,k) )
   END DO

   IF (ntest > 0) THEN
      totwbc = 0.
      totwblc = 0.
      DO k = 1, ms
         totwbc = totwbc + soil%zse(k) * REAL( ssoil%wb(:,k) )
         totwblc = totwblc + soil%zse(k) * REAL( ssoil%wblf(:,k) )
      END DO
   END IF

END SUBROUTINE smoisturev



SUBROUTINE snowdensity (dels, ssoil, soil)
   
   REAL, INTENT(IN) :: dels   ! integration time step (s)

   TYPE(soil_snow_type),      INTENT(INOUT) :: ssoil 
    
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil

   INTEGER, DIMENSION(mp,3) :: ssnow_isflag_ssdn 
   REAL, DIMENSION(mp) :: ssoil_tgg_min1
   REAL, DIMENSION(mp,3) :: dels_ssdn, ssoil_tgg_min
     
   ssnow_isflag_ssdn = SPREAD( ssoil%isflag,2,mp) 
   
   dels_ssdn = SPREAD( SPREAD( dels, 1, mp ), 2,  mp ) 
   ssoil_tgg_min1 = MIN( tfrz, ssoil%tgg(:,1) )
   
   !print *, "jhan: ", shape( ssnow_isflag_ssdn )
   !print *, "jhan: ", shape( dels_ssdn )
     
   WHERE( ssoil%snowd > 0.1 .AND. ssoil%isflag == 0 )
      
      ssoil%ssdn(:,1) = MIN( max_ssdn, MAX( 120.0, ssoil%ssdn(:,1) + dels      &
                        * ssoil%ssdn(:,1) * 3.1e-6 * EXP( -0.03 * ( 273.15 -     &
                        ssoil_tgg_min1 ) - MERGE( 0.046, 0.0,     &
                        ssoil%ssdn(:,1) >= 150.0 ) * ( ssoil%ssdn(:,1) - 150.0)&
                        ) ) )

      !ssoil%ssdn(:,1) = MIN( max_ssdn, ssoil%ssdn(:,1) + dels * 9.806          &
      !                  * ssoil%ssdn(:,1) * 0.75 * ssoil%snowd /               &
      !                  ( 3.0e7 * EXP( 0.021 * ssoil%ssdn(:,1) + 0.081 *( tfrz &
      !                  - ssoil_tgg_min1 ) ) ) )
      !!!!!jhan
            ssoil%ssdn(:,1) = MIN(max_ssdn,ssoil%ssdn(:,1) + dels * 9.806 &
          & * ssoil%ssdn(:,1) * 0.75 * ssoil%snowd &
          & / (3.0e7 * EXP(0.021 * ssoil%ssdn(:,1) + 0.081 &
          & * (273.15 - MIN(tfrz, ssoil%tgg(:,1))))))

      WHERE( soil%isoilm /= 9 ) ssoil%ssdn(:,1) = MIN( 450.0, ssoil%ssdn(:,1) )

      ssoil%sconds(:,1) = MAX( 0.2, MIN( 2.876e-6 * ssoil%ssdn(:,1)**2         &
                          + 0.074, max_sconds ) )

      ssoil%sconds(:,2) = ssoil%sconds(:,1) 
      ssoil%sconds(:,3) = ssoil%sconds(:,1) 
      
      ssoil%ssdnn = ssoil%ssdn(:,1)
      
      ssoil%ssdn(:,2) = ssoil%ssdn(:,1)
      ssoil%ssdn(:,3) = ssoil%ssdn(:,1)
    
   END WHERE
  
   !WHERE ( spread( ssoil%isflag,2,mp) == 1 )
      
   !   ssoil%ssdn = ssoil%ssdn + spread( spread( dels,1,mp),2, mp)  * ssoil%ssdn * 3.1e-6      &
    !                    * EXP( -0.03 * ( tfrz - MIN( tfrz, ssoil%tggsn ) )&
   !                     - MERGE( 0.046, 0.0, ssoil%ssdn >= 150.0)         &
    !                    * ( ssoil%ssdn - 150.0 ) )
     
     
     !ssoil%ssdn(:,1) = ssoil%ssdn(:,1) + dels * ssoil%ssdn(:,1) * 3.1e-6      &
     !                  * EXP( -0.03 * ( tfrz - MIN( tfrz, ssoil%tggsn(:,1) ) )&
     !                  - MERGE( 0.046, 0.0, ssoil%ssdn(:,1) >= 150.0)         &
     !                  * ( ssoil%ssdn(:,1) - 150.0 ) )

     !ssoil%ssdn(:,2) = ssoil%ssdn(:,2) + dels * ssoil%ssdn(:,2) * 3.1e-6      &
     !                  * EXP( -0.03 * ( tfrz - MIN( tfrz, ssoil%tggsn(:,2) ) )&
     !                  - MERGE( 0.046, 0.0, ssoil%ssdn(:,2) >= 150.0 )        &
     !                  * ( ssoil%ssdn(:,2) - 150.0 ) )

      ! ssoil%ssdn(:,3) = ssoil%ssdn(:,3) + dels * ssoil%ssdn(:,3) * 3.1e-6 &
      !      * EXP( -0.03 * (273.15 - MIN(tfrz, ssoil%tggsn(:,3))) &
      !      - MERGE(0.046, 0.0, ssoil%ssdn(:,3) >= 150.0) &
      !      * (ssoil%ssdn(:,3) - 150.0) )
       
    !   ssoil%ssdn(:,1) = ssoil%ssdn(:,1) + dels * 9.806 * ssoil%ssdn(:,1) &
    !        * ssoil%t_snwlr*ssoil%ssdn(:,1) &
    !        / (3.0e7 * EXP(.021 * ssoil%ssdn(:,1) + 0.081 &
    !        * (273.15 - MIN(tfrz, ssoil%tggsn(:,1)))))
    !   ssoil%ssdn(:,2) = ssoil%ssdn(:,2) + dels * 9.806 * ssoil%ssdn(:,2) &
    !        * (ssoil%t_snwlr * ssoil%ssdn(:,1) + 0.5 * ssoil%smass(:,2) ) &
    !        / (3.0e7 * EXP(.021 * ssoil%ssdn(:,2) + 0.081 &
    !        * (273.15 - MIN(tfrz, ssoil%tggsn(:,2)))))
    !   ssoil%ssdn(:,3) = ssoil%ssdn(:,3) + dels * 9.806 * ssoil%ssdn(:,3) &
    !        * (ssoil%t_snwlr*ssoil%ssdn(:,1) + ssoil%smass(:,2) &
    !        + 0.5*ssoil%smass(:,3)) &
    !        / (3.0e7 * EXP(.021 * ssoil%ssdn(:,3) + 0.081 &
    !        * (273.15 - MIN(tfrz, ssoil%tggsn(:,3)))))
    !  
    !  ssoil%sdepth =  ssoil%smass / ssoil%ssdn 
    !  
    !  !ssoil%sdepth(:,1) =  ssoil%smass(:,1) / ssoil%ssdn(:,1) 
    !  !ssoil%sdepth(:,2) =  ssoil%smass(:,2) / ssoil%ssdn(:,2) 
    !  !ssoil%sdepth(:,3) =  ssoil%smass(:,3) / ssoil%ssdn(:,3) 
    !  
    !  ssoil%ssdnn = (ssoil%ssdn(:,1) * ssoil%smass(:,1) + ssoil%ssdn(:,2) &
    !        * ssoil%smass(:,2) + ssoil%ssdn(:,3) * ssoil%smass(:,3) ) &
    !        / ssoil%snowd
    !  
    !  ssoil%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,1) ** 2 &
    !                                                    & + 0.074, max_sconds) )
    !  ssoil%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,2) ** 2 &
    !                                                    & + 0.074, max_sconds) )
    !  ssoil%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,3) ** 2 &
!                                                        & + 0.074, max_sconds) )
  ! END WHERE


   WHERE (ssoil%isflag == 1)
      
      ssoil%ssdn(:,1) = ssoil%ssdn(:,1) + dels * ssoil%ssdn(:,1) * 3.1e-6 &
            * EXP( -0.03 * (273.15 - MIN(tfrz, ssoil%tggsn(:,1))) &
            - MERGE(0.046, 0.0, ssoil%ssdn(:,1) >= 150.0) &
            * (ssoil%ssdn(:,1) - 150.0) )
      
      ssoil%ssdn(:,2) = ssoil%ssdn(:,2) + dels * ssoil%ssdn(:,2) * 3.1e-6 &
            * EXP( -0.03 * (273.15 - MIN(tfrz, ssoil%tggsn(:,2))) &
            - MERGE(0.046, 0.0, ssoil%ssdn(:,2) >= 150.0) &
            * (ssoil%ssdn(:,2) - 150.0) )
      
      ssoil%ssdn(:,3) = ssoil%ssdn(:,3) + dels * ssoil%ssdn(:,3) * 3.1e-6 &
            * EXP( -0.03 * (273.15 - MIN(tfrz, ssoil%tggsn(:,3))) &
            - MERGE(0.046, 0.0, ssoil%ssdn(:,3) >= 150.0) &
            * (ssoil%ssdn(:,3) - 150.0) )
      
      ssoil%ssdn(:,1) = ssoil%ssdn(:,1) + dels * 9.806 * ssoil%ssdn(:,1) &
            * ssoil%t_snwlr*ssoil%ssdn(:,1) &
            / (3.0e7 * EXP(.021 * ssoil%ssdn(:,1) + 0.081 &
            * (273.15 - MIN(tfrz, ssoil%tggsn(:,1)))))
      
      ssoil%ssdn(:,2) = ssoil%ssdn(:,2) + dels * 9.806 * ssoil%ssdn(:,2) &
            * (ssoil%t_snwlr * ssoil%ssdn(:,1) + 0.5 * ssoil%smass(:,2) ) &
            / (3.0e7 * EXP(.021 * ssoil%ssdn(:,2) + 0.081 &
            * (273.15 - MIN(tfrz, ssoil%tggsn(:,2)))))
      
      ssoil%ssdn(:,3) = ssoil%ssdn(:,3) + dels * 9.806 * ssoil%ssdn(:,3) &
            * (ssoil%t_snwlr*ssoil%ssdn(:,1) + ssoil%smass(:,2) &
            + 0.5*ssoil%smass(:,3)) &
            / (3.0e7 * EXP(.021 * ssoil%ssdn(:,3) + 0.081 &
            * (273.15 - MIN(tfrz, ssoil%tggsn(:,3)))))
      
      ssoil%sdepth(:,1) =  ssoil%smass(:,1) / ssoil%ssdn(:,1) 
      ssoil%sdepth(:,2) =  ssoil%smass(:,2) / ssoil%ssdn(:,2) 
      ssoil%sdepth(:,3) =  ssoil%smass(:,3) / ssoil%ssdn(:,3) 
      
      ssoil%ssdnn = (ssoil%ssdn(:,1) * ssoil%smass(:,1) + ssoil%ssdn(:,2) &
            * ssoil%smass(:,2) + ssoil%ssdn(:,3) * ssoil%smass(:,3) ) &
            / ssoil%snowd
      
      ssoil%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,1) ** 2 &
                                                        & + 0.074, max_sconds) )
      ssoil%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,2) ** 2 &
                                                        & + 0.074, max_sconds) )
      ssoil%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,3) ** 2 &
                                                        & + 0.074, max_sconds) )
   END WHERE

END SUBROUTINE snowdensity


SUBROUTINE snow_melting (dels, snowmlt, ssoil, soil )

   USE cable_common_module
   
   REAL, INTENT(IN) :: dels   ! integration time step (s)
   
   REAL, DIMENSION(mp), INTENT(OUT) :: snowmlt ! snow melt   
   
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil
   TYPE(soil_snow_type), INTENT(INOUT)   :: ssoil  ! soil+snow variables
  
   INTEGER                 :: k,j 
   
   REAL, DIMENSION(mp) ::                                                      &
      osm,     & !
      sgamm,   & !
      snowflx    !
   
   REAL, DIMENSION(mp,0:3) :: smelt1

   snowmlt= 0.0
   smelt1 = 0.0
    
   DO j=1,mp  
      
      IF( ssoil%snowd(j) > 0.0 .AND. ssoil%isflag(j) == 0                      &
          .AND. ssoil%tgg(j,1) >= tfrz ) THEN

         ! snow covered land
         ! following done in sflux  via  ga= ... +cls*egg + ...
         ! ** land,snow,melting
         snowflx(j) = REAL((ssoil%tgg(j,1) - tfrz) * ssoil%gammzz(j,1),r_1)
         
         ! prevent snow depth going negative
         snowmlt(j) = MIN(snowflx(j) / hlf, ssoil%snowd(j) )
       
         ssoil%dtmlt(j,1) = ssoil%dtmlt(j,1) + snowmlt(j) * hlf                &
                            / ssoil%gammzz(j,1)

         ssoil%snowd(j) = ssoil%snowd(j) - snowmlt(j)
         ssoil%tgg(j,1) = REAL(                                                &
                          ssoil%tgg(j,1) - snowmlt(j) * hlf / ssoil%gammzz(j,1))

      ENDIF
    
   END DO

   smelt1(:,0) = 0.0
   
   DO k = 1, 3
    
      !where there is snow 
      WHERE( ssoil%snowd > 0.0 .AND. ssoil%isflag > 0 )

         sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
       
         ! snow melt refreezing
         snowflx = smelt1(:,k-1) * hlf / dels

         ssoil%tggsn(:,k) = ssoil%tggsn(:,k) + ( snowflx * dels +              &
                            smelt1(:,k-1)*cswat *( tfrz-ssoil%tggsn(:,k) ) ) / &
                            ( sgamm + cswat*smelt1(:,k-1) )
         
         ! increase density due to snowmelt
         osm = ssoil%smass(:,k)
         ssoil%smass(:,k) = ssoil%smass(:,k) + smelt1(:,k-1)
         ssoil%ssdn(:,k) = MAX( 120.0, MIN( ssoil%ssdn(:,k) * osm /            &
                           ssoil%smass(:,k) + rhowat * ( 1.0 - osm /           &
                           ssoil%smass(:,k)), max_ssdn ) )

         WHERE( soil%isoilm /= 9 )                                             &
            ssoil%ssdn(:,k) = MIN( 450.0, ssoil%ssdn(:,k) )

         ssoil%sdepth(:,k) = ssoil%smass(:,k) / ssoil%ssdn(:,k)
         
         sgamm = ssoil%smass(:,k) * cgsnow
        
         smelt1(:,k) = 0.0
         
         ! snow melting
         WHERE (ssoil%tggsn(:,k) > tfrz)
           
            snowflx = ( ssoil%tggsn(:,k) - tfrz ) * sgamm
           
            smelt1(:,k) = MIN( snowflx / hlf, 0.6 * ssoil%smass(:,k) )
            
            ssoil%dtmlt(:,k) = ssoil%dtmlt(:,k) + smelt1(:,k) * hlf / sgamm
            
            osm = ssoil%smass(:,k)
            
            ssoil%smass(:,k) = ssoil%smass(:,k) - smelt1(:,k)

            ssoil%tggsn(:,k) = ssoil%tggsn(:,k) - smelt1(:,k) * hlf / sgamm

            ssoil%sdepth(:,k) = ssoil%smass(:,k) / ssoil%ssdn(:,k)
       
         END WHERE 
         ! END snow melting
     
      END WHERE    
      ! END where there is snow 
   
   END DO
   
   WHERE( ssoil%snowd > 0.0 .AND. ssoil%isflag > 0 )
      snowmlt = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)
      ssoil%snowd = ssoil%snowd - snowmlt
   END WHERE

END SUBROUTINE snow_melting



SUBROUTINE snow_accum ( dels,  canopy, met, ssoil, soil )

USE cable_common_module

   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(canopy_type), INTENT(INOUT)         :: canopy ! vegetation variables
   TYPE(met_type), INTENT(INOUT)            :: met   ! all met forcing
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil ! soil+snow variables
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
   
   REAL, DIMENSION(mp) ::                                                      &
      osm,     & !
      sgamm,   & !
      snowmlt, & !
      xxx        !

   INTEGER             :: k

   WHERE (canopy%precis > 0.0 .and. ssoil%isflag == 0)
      ! accumulate solid part
      ssoil%snowd = MAX( ssoil%snowd + met%precip_sn, 0.0 ) 
      
      canopy%precis = canopy%precis - met%precip_sn
      
      ssoil%ssdn(:,1) = MAX( 120.0, ssoil%ssdn(:,1)                            &
                        * ssoil%osnowd / MAX( 0.01, ssoil%snowd )              &
                        + 120.0 * met%precip_sn / MAX( 0.01, ssoil%snowd ) )
      
      ssoil%ssdnn = ssoil%ssdn(:,1)
      
      WHERE( canopy%precis > 0.0 .AND. ssoil%tgg(:,1) < tfrz )
         
         ssoil%snowd = MAX(ssoil%snowd + canopy%precis, 0.0)
        
         ssoil%tgg(:,1) = ssoil%tgg(:,1) + canopy%precis * hlf                 &
                          / ( REAL( ssoil%gammzz(:,1) ) + cswat *canopy%precis )  
         ! change density due to water being added 
         ssoil%ssdn(:,1) = MIN( max_ssdn, MAX( 120.0, ssoil%ssdn(:,1)          &
                           * ssoil%osnowd / MAX( 0.01, ssoil%snowd ) + rhowat  &
                           * canopy%precis / MAX( 0.01, ssoil%snowd )  ) )

         WHERE( soil%isoilm /= 9 )                                             &
            ssoil%ssdn(:,1) = MIN( 450.0, ssoil%ssdn(:,1) )

         canopy%precis = 0.0
         ssoil%ssdnn = ssoil%ssdn(:,1)
      
      END WHERE
   
   END WHERE ! (canopy%precis > 0. .and. ssoil%isflag == 0) 

   WHERE (canopy%precis > 0.0 .and.  ssoil%isflag > 0)
      
      ! add solid precip
      ssoil%snowd = MAX( ssoil%snowd + met%precip_sn, 0.0 )

      canopy%precis = canopy%precis - met%precip_sn  ! remaining liquid precip

      ! update top snow layer with fresh snow
      osm = ssoil%smass(:,1)
      ssoil%smass(:,1) = ssoil%smass(:,1) + met%precip_sn
      ssoil%ssdn(:,1) = MAX( 120.0,ssoil%ssdn(:,1) * osm / ssoil%smass(:,1)    &
                        + 120.0 * met%precip_sn / ssoil%smass(:,1) )

      ssoil%sdepth(:,1) = MAX( 0.02, ssoil%smass(:,1) / ssoil%ssdn(:,1) )

      ! add liquid precip
      WHERE( canopy%precis > 0.0 )
        
         ssoil%snowd = MAX( ssoil%snowd + canopy%precis, 0.0 )
         sgamm = ssoil%ssdn(:,1) * cgsnow * ssoil%sdepth(:,1)
         osm = ssoil%smass(:,1)
         
         ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + canopy%precis * hlf             &
                            * osm / (sgamm * ssoil%osnowd )
         ssoil%smass(:,1) = ssoil%smass(:,1) + canopy%precis                   &
                            * osm/ssoil%osnowd

         ssoil%ssdn(:,1) = MAX( 120.0, MIN( ssoil%ssdn(:,1) * osm /            &
                           ssoil%smass(:,1) +  rhowat *                        &
                           ( 1.0 - osm / ssoil%smass(:,1) ), max_ssdn ) )

         WHERE( soil%isoilm /= 9 )                                             &
            ssoil%ssdn(:,1) = MIN( 450.0, ssoil%ssdn(:,1) )

         ssoil%sdepth(:,1) = ssoil%smass(:,1)/ssoil%ssdn(:,1)

         !layer 2
         sgamm = ssoil%ssdn(:,2) * cgsnow * ssoil%sdepth(:,2)
         osm = ssoil%smass(:,2)
         ssoil%tggsn(:,2) = ssoil%tggsn(:,2) + canopy%precis * hlf             &
                            * osm / ( sgamm * ssoil%osnowd )
         ssoil%smass(:,2) = ssoil%smass(:,2) + canopy%precis                   &
                            * osm / ssoil%osnowd
         ssoil%ssdn(:,2) = MAX( 120.0, MIN( ssoil%ssdn(:,2) * osm /            &
                           ssoil%smass(:,2) + rhowat *                         &
                           ( 1.0 - osm / ssoil%smass(:,2) ), max_ssdn ) )

         WHERE( soil%isoilm /= 9 )                                             &
            ssoil%ssdn(:,2) = MIN( 450.0, ssoil%ssdn(:,2) )

         ssoil%sdepth(:,2) = ssoil%smass(:,2) / ssoil%ssdn(:,2)

         !layer 3        
         sgamm = ssoil%ssdn(:,3) * cgsnow * ssoil%sdepth(:,3)
         osm = ssoil%smass(:,3)
         ssoil%tggsn(:,3) = ssoil%tggsn(:,3) + canopy%precis * hlf             &
                            * osm / ( sgamm * ssoil%osnowd )
         ssoil%smass(:,3) = ssoil%smass(:,3) + canopy%precis                   &
                            * osm / ssoil%osnowd
        ssoil%ssdn(:,3) = MAX( 120.0, MIN( ssoil%ssdn(:,3) * osm /             &
                          ssoil%smass(:,3) + rhowat *                          &
                          ( 1.0 - osm / ssoil%smass(:,3) ), max_ssdn ) )

         WHERE( soil%isoilm /= 9 )                                             &
            ssoil%ssdn(:,3) = MIN(450.0,ssoil%ssdn(:,3))

         ssoil%sdepth(:,3) = ssoil%smass(:,3) / ssoil%ssdn(:,3)

         canopy%precis = 0.0
      
      END WHERE
   
   END WHERE


   ! 'fess' is for soil evap and 'fes' is for soil evap plus soil puddle evap
   canopy%segg = canopy%fess / hl
   canopy%segg = ( canopy%fess + canopy%fes_cor ) / hl

   ! Initialise snow evaporation:
   ssoil%evapsn = 0
    
   ! Snow evaporation and dew on snow
   WHERE( ssoil%snowd > 0.1 )
      
      ssoil%evapsn = dels * ( canopy%fess + canopy%fes_cor ) / ( hl + hlf )
      xxx = ssoil%evapsn

      WHERE( ssoil%isflag == 0 .AND. canopy%fess + canopy%fes_cor.GT. 0.0 )    &
         ssoil%evapsn = MIN( ssoil%snowd, xxx )
          
      WHERE( ssoil%isflag  > 0 .AND. canopy%fess + canopy%fes_cor .GT. 0.0 )   &
         ssoil%evapsn = MIN( 0.9 * ssoil%smass(:,1), xxx )

      ssoil%snowd = ssoil%snowd - ssoil%evapsn
      
      WHERE( ssoil%isflag > 0 )
         ssoil%smass(:,1) = ssoil%smass(:,1)  - ssoil%evapsn
         ssoil%sdepth(:,1) = MAX( 0.02, ssoil%smass(:,1) / ssoil%ssdn(:,1) )
      END WHERE

      canopy%segg = 0.0

   END WHERE

END SUBROUTINE snow_accum 


SUBROUTINE surfbv (dels, met, ssoil, soil, veg, canopy )

   USE cable_common_module

   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(canopy_type), INTENT(IN)       :: canopy
   
   TYPE(met_type),       INTENT(INOUT) :: met    ! all met forcing
   TYPE(soil_snow_type), INTENT(INOUT) :: ssoil  ! soil+snow variables
   
   TYPE(veg_parameter_type),  INTENT(IN)     :: veg
   TYPE(soil_parameter_type), INTENT(INOUT)  :: soil  ! soil parameters

!jhan:cable.nml
   INTEGER, PARAMETER      :: nglacier = 2 ! 0 original, 1 off, 2 new Eva

   REAL, DIMENSION(mp) ::                                                      &
      rnof5,      & !
      sfact,      & !
      sgamm,      & !
      smasstot,   & !
      talb,       & ! snow albedo
      tmp,        & ! temporary value
      xxx           !

   REAL, DIMENSION(mp,0:3) :: smelt1
    
   INTEGER :: k

   CALL smoisturev( dels, ssoil, soil, veg )

   DO k = 1, ms
      xxx = REAL( soil%ssat,r_2 )
      ssoil%rnof1 = ssoil%rnof1 + REAL( MAX( ssoil%wb(:,k) - xxx, 0.0_r_2 )  &
                    * 1000.0 )  * soil%zse(k)
      ssoil%wb(:,k) = MAX( 0.01, MIN( ssoil%wb(:,k), xxx ) )
   END DO

   ! for deep runoff use wb-sfc, but this value not to exceed .99*wb-wbice
   ! account for soil/ice cracking
   ! fracm = MIN(0.2, 1. - MIN(1., ssoil%wb(:,ms) / soil%sfc ) )
   ! ssoil%wb(:,ms) = ssoil%wb(:,ms) &
   !                  + fracm*ssoil%rnof1/(1000.0*soil%zse(ms))
   ! ssoil%rnof1 = (1. - fracm) * ssoil%rnof1 

   ! Scaling  runoff to kg/m^2/s to match rest of the model
   ssoil%sinfil = 0.0
   WHERE( veg%iveg == 16 )
      ssoil%sinfil = MIN( ssoil%rnof1, ssoil%wb_lake + MAX( 0.,canopy%segg ) )
      ssoil%rnof1 = MAX( 0.0, ssoil%rnof1 - ssoil%sinfil )
      ssoil%wb_lake = ssoil%wb_lake - ssoil%sinfil
      ssoil%rnof2 = MAX( 0.0, ssoil%rnof2 - ssoil%wb_lake )
   ENDWHERE

!jhan:replace nested wheres 

   !---  glacier formation
   rnof5= 0.

   IF (nglacier == 2) THEN
      
      smelt1=0.
      WHERE( ssoil%snowd > max_glacier_snowd )

         rnof5 = MIN( 0.1, ssoil%snowd - max_glacier_snowd )

         !---- change local tg to account for energy - clearly not best method
         WHERE( ssoil%isflag == 0 )
            smasstot = 0.0
            ssoil%tgg(:,1) = ssoil%tgg(:,1) - rnof5 * hlf                      &
                             / REAL( ssoil%gammzz(:,1) )
            ssoil%snowd = ssoil%snowd - rnof5
         ELSEWHERE
            smasstot = ssoil%smass(:,1) + ssoil%smass(:,2) + ssoil%smass(:,3)
         END WHERE

      END WHERE

      DO k = 1, 3
         
         WHERE( ssoil%snowd > max_glacier_snowd  .AND.  ssoil%isflag > 0 )
            sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
            smelt1(:,k) = MIN( rnof5 * ssoil%smass(:,k) / smasstot,            &
                          0.2 * ssoil%smass(:,k) )
            ssoil%smass(:,k) = ssoil%smass(:,k) - smelt1(:,k)
            ssoil%snowd = ssoil%snowd - smelt1(:,k)
         END WHERE
      
      END DO
   
      WHERE( ssoil%isflag > 0 ) rnof5 = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)
   
   END IF

   ssoil%rnof1 = ssoil%rnof1 / dels + rnof5/dels
   ssoil%rnof2 = ssoil%rnof2 / dels
   ssoil%runoff = ssoil%rnof1 + ssoil%rnof2 

END SUBROUTINE surfbv

  
! calculates temperatures of the soil
! tgg - new soil/snow temperature
! ga - heat flux from the atmosphere (ground heat flux)
! ccnsw - soil thermal conductivity, including water/ice
SUBROUTINE stempv(dels, canopy, ssoil, soil)
   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(canopy_type),    INTENT(INOUT) :: canopy
   TYPE(soil_snow_type), INTENT(INOUT) :: ssoil
   
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil
   
   REAL, DIMENSION(mp) ::                                                      &
      coefa, coefb,  & !
      sgamm            ! 

   REAL(r_2), DIMENSION(mp) ::                                                 &
      dtg,     & !
      ew,      & !
      xx,      & !
      wblfsp     ! 

   REAL(r_2), DIMENSION(mp,ms) ::                                              &
      ccnsw  ! soil thermal conductivity (incl water/ice)

   REAL(r_2), DIMENSION(mp, -2:ms) ::                                          &
      at, bt, ct, rhs !
   
   REAL(r_2), DIMENSION(mp,-2:ms+1) :: coeff
   
   REAL(r_2), DIMENSION(mp,ms+3)    :: tmp_mat ! temp. matrix for tggsn & tgg
   
   INTEGER :: j,k
   REAL :: snow_ccnsw, exp_arg
   LOGICAL :: direct2min = .FALSE.


   at = 0.0
   bt = 1.0
   ct = 0.0
   coeff = 0.0
   snow_ccnsw = 2.0

   DO k = 1, ms
      
      DO j = 1, mp
      
         IF( soil%isoilm(j) == 9 ) THEN
            ccnsw(j,k) = snow_ccnsw
         ELSE
            ew(j) = ssoil%wblf(j,k) * soil%ssat(j)
            exp_arg = ( ew(j) * LOG( 60.0 ) ) + ( ssoil%wbfice(j,k)            &
                      * soil%ssat(j) * LOG( 250.0 ) )

            IF( exp_arg > 30 ) direct2min = .TRUE.
            
            IF( direct2min) THEN
               
               ccnsw(j,k) = 1.5 * MAX( 1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5 *      &
                            soil%ssat(j) /                                     &
                            MIN( ew(j), 0.5_r_2 * soil%ssat(j) ) ) ) )

            ELSE         
               
               ccnsw(j,k) = MIN( soil%cnsd(j) * EXP( exp_arg ), 1.5_r_2 )      &
                            * MAX( 1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5 *          &
                            soil%ssat(j) /                                     &
                            MIN( ew(j), 0.5_r_2 * soil%ssat(j) ) ) ) )
            
            ENDIF          
           
            direct2min = .FALSE.
        
         ENDIF 
      
      END DO
    
   END DO
    
   xx = 0. 
    
   WHERE(ssoil%isflag == 0)
      xx = MAX( 0., ssoil%snowd / ssoil%ssdnn )
      ccnsw(:,1) = ( ccnsw(:,1) - 0.2 ) * ( soil%zse(1) / ( soil%zse(1) + xx ) &
                   ) + 0.2
   END WHERE
    
   DO k = 3, ms
      
      WHERE (ssoil%isflag == 0)
         coeff(:,k) = 2.0 / ( soil%zse(k-1) / ccnsw(:,k-1) + soil%zse(k) /     &
                      ccnsw(:,k) )
      END WHERE
   END DO

   k = 1
   WHERE( ssoil%isflag == 0 )
      coeff(:,2) = 2.0 / ( ( soil%zse(1) + xx ) / ccnsw(:,1) + soil%zse(2) /   &
                   ccnsw(:,2) )
      coefa = 0.0
      coefb = REAL( coeff(:,2) )

      wblfsp = ssoil%wblf(:,k)
      
      xx = soil%css * soil%rhosoil
      
      ssoil%gammzz(:,k) = MAX( ( 1.0 - soil%ssat ) * soil%css * soil%rhosoil   &
                          + soil%ssat * ( wblfsp * cswat * rhowat +            &
                          ssoil%wbfice(:,k) * csice * rhowat * 0.9 ), xx )     &
                          * soil%zse(k)

      ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd

      dtg = dels / ssoil%gammzz(:,k)
      
      at(:,k) = - dtg * coeff(:,k)
      ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
      bt(:,k) = 1.0 - at(:,k) - ct(:,k)

   END WHERE
   
   DO k = 2, ms
      
      WHERE( ssoil%isflag == 0 )
         
         wblfsp = ssoil%wblf(:,k)
         xx = soil%css * soil%rhosoil

         ssoil%gammzz(:,k) = MAX( ( 1.0 - soil%ssat ) * soil%css * soil%rhosoil&
                             + soil%ssat * ( wblfsp * cswat * rhowat +         &
                             ssoil%wbfice(:,k) * csice * rhowat * 0.9 ), xx )  &
                             * soil%zse(k)

         dtg = dels / ssoil%gammzz(:,k)
         at(:,k) = - dtg * coeff(:,k)
         ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
         bt(:,k) = 1.0 - at(:,k) - ct(:,k)
      
      END WHERE
    
   END DO
   
   WHERE( ssoil%isflag == 0 )
      bt(:,1) = bt(:,1) - canopy%dgdtg * dels / ssoil%gammzz(:,1)
      ssoil%tgg(:,1) = ssoil%tgg(:,1) + ( canopy%ga - ssoil%tgg(:,1)           &
                       * REAL( canopy%dgdtg ) ) * dels /                       &
                       REAL( ssoil%gammzz(:,1) )
   END WHERE
   
   coeff(:,1-3) = 0.0  ! SO DOES THIS MEAN coeff(:,-2) ??

!jhan:fix this into a loop over 3 lyers  
   ! 3-layer snow points done here
   WHERE( ssoil%isflag /= 0 )
      
      ssoil%sconds(:,1) = MAX( 0.2, MIN( 2.876e-6 * ssoil%ssdn(:,1)**2         &
                          + 0.074, max_sconds ) )
      ssoil%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,2)**2 &
                        & + 0.074, max_sconds) )
      ssoil%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,3)**2 &
                        & + 0.074, max_sconds) )
      coeff(:,-1) = 2.0 / (ssoil%sdepth(:,1) / ssoil%sconds(:,1) &
                       & + ssoil%sdepth(:,2) / ssoil%sconds(:,2) )
      coeff(:,0) = 2.0 / (ssoil%sdepth(:,2) / ssoil%sconds(:,2) &
                      & + ssoil%sdepth(:,3) / ssoil%sconds(:,3) )
      coeff(:,1) = 2.0 / (ssoil%sdepth(:,3) / ssoil%sconds(:,3) &
                      & + soil%zse(1) / ccnsw (:,1) )
   END WHERE
    
   DO k = 2, ms
      
      WHERE( ssoil%isflag /= 0 )                                               &
         coeff(:,k) = 2.0 / ( soil%zse(k-1) / ccnsw(:,k-1) + soil%zse(k) /     &
                      ccnsw(:,k) )

   END DO
   
   WHERE( ssoil%isflag /= 0 )
      coefa = REAL( coeff (:,-1) )
      coefb = REAL( coeff (:,1) )
   END WHERE
   
   DO k = 1, 3
      
      WHERE( ssoil%isflag /= 0 ) 
         sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
         dtg = dels / sgamm
         at(:,k-3) = - dtg * coeff(:,k-3)
         ct(:,k-3) = - dtg * coeff(:,k-2)
         bt(:,k-3) = 1.0 - at(:,k-3) - ct(:,k-3)
      END WHERE
   
   END DO
   
   DO k = 1, ms
     
      WHERE( ssoil%isflag /= 0 )
         wblfsp = ssoil%wblf(:,k)
         xx = soil%css * soil%rhosoil
         
         ssoil%gammzz(:,k) = MAX( ( 1.0 - soil%ssat ) * soil%css *             &
                             soil%rhosoil + soil%ssat * ( wblfsp * cswat *     &
                             rhowat + ssoil%wbfice(:,k) * csice * rhowat *     &
                             0.9) , xx ) * soil%zse(k)

         dtg = dels / ssoil%gammzz(:,k)
         at(:,k) = - dtg * coeff(:,k)
         ct(:,k) = - dtg * coeff(:,k + 1) ! c3(ms)=0 & not really used
         bt(:,k) = 1.0 - at(:,k) - ct(:,k)
      
      END WHERE
   
   END DO

   WHERE( ssoil%isflag /= 0 )
      sgamm = ssoil%ssdn(:,1) * cgsnow * ssoil%sdepth(:,1)
      
      bt(:,-2) = bt(:,-2) - canopy%dgdtg * dels / sgamm
      
      ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + ( canopy%ga - ssoil%tggsn(:,1 )    &
                         * REAL( canopy%dgdtg ) ) * dels / sgamm

      rhs(:,1-3) = ssoil%tggsn(:,1)
   END WHERE
 
   !     note in the following that tgg and tggsn are processed together
   tmp_mat(:,1:3) = REAL(ssoil%tggsn,r_2) 
   tmp_mat(:,4:(ms+3)) = REAL(ssoil%tgg,r_2)

   CALL trimb( at, bt, ct, tmp_mat, ms + 3 ) 
   
   ssoil%tggsn = REAL( tmp_mat(:,:3) )
   ssoil%tgg   = REAL( tmp_mat(:,4:(ms+3)) )
   canopy%sghflux = coefa * ( ssoil%tggsn(:,1) - ssoil%tggsn(:,2) )
   canopy%ghflux = coefb * ( ssoil%tgg(:,1) - ssoil%tgg(:,2) ) ! +ve downwards

END SUBROUTINE stempv




SUBROUTINE snowcheck(dels, ssoil, soil, met )
   
   USE cable_common_module
   
   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(soil_snow_type), INTENT(INOUT) :: ssoil
   TYPE(met_type),       INTENT(INOUT) :: met ! all met forcing
   
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
   
   INTEGER :: k,j
   
   !jhan:cable.nml
   LOGICAL :: cable_runtime_coupled = .FALSE.

   DO j=1,mp

      IF( ssoil%snowd(j) <= 0.0 ) THEN
         
         ssoil%isflag(j) = 0
         ssoil%ssdn(j,:) = 120.0
         ssoil%ssdnn(j) = 120.0
         ssoil%tggsn(j,:) = tfrz
         ssoil%sdepth(j,1) = ssoil%snowd(j) / ssoil%ssdn(j,1)
          
         ssoil%sdepth(j,2) = 0.
         ssoil%sdepth(j,3) = 0.

         ssoil%smass(j,1) = ssoil%snowd(j)
         ssoil%smass(j,2) = 0.0     ! EK to fix -ve sdepth 21Dec2007
         ssoil%smass(j,3) = 0.0     ! EK to fix -ve sdepth 21Dec2007
      
      ! in loop: IF( ssoil%snowd(j) <= 0.0 ) THEN
      ELSEIF( ssoil%snowd(j) < snmin * ssoil%ssdnn(j) ) THEN

         IF( ssoil%isflag(j) == 1 ) THEN
            ssoil%ssdn(j,1) = ssoil%ssdnn(j)
            ssoil%tgg(j,1) = ssoil%tggsn(j,1)
         ENDIF 

         ssoil%isflag(j) = 0
         ssoil%ssdnn(j) = MIN( 400.0, MAX( 120.0, ssoil%ssdn(j,1) ) ) 
     
         ssoil%tggsn(j,:) = MIN( tfrz,ssoil%tgg(j,1) )

         ssoil%sdepth(j,1) = ssoil%snowd(j) / ssoil%ssdn(j,1)
         ssoil%sdepth(j,2) = 0.0     
         ssoil%sdepth(j,3) = 0.0     

         ssoil%smass(j,1) = ssoil%snowd(j)     
         ssoil%smass(j,2) = 0.0     
         ssoil%smass(j,3) = 0.0     

         ssoil%ssdn(j,:) = ssoil%ssdnn(j)
         
         IF( .NOT.cable_runtime_coupled ) THEN
            IF( soil%isoilm(j) == 9 .AND. ktau_gl <= 2 )                       &
               ssoil%ssdnn(j) = 700.0
         ENDIF
      
      
      ELSE ! in loop: IF( ssoil%snowd(j) <= 0.0 ) THEN
           ! sufficient snow now for 3 layer snowpack
         
         IF( ssoil%isflag(j) == 0 ) THEN

            ssoil%tggsn(j,:) = MIN( tfrz, ssoil%tgg(j,1) )

            ssoil%ssdn(j,2) = ssoil%ssdn(j,1)
            ssoil%ssdn(j,3) = ssoil%ssdn(j,1)

            IF( .NOT. cable_runtime_coupled) THEN
               IF( soil%isoilm(j) == 9 .AND. ktau_gl <= 2 ) THEN
                  ssoil%ssdn(j,1)  = 450.0
                  ssoil%ssdn(j,2)  = 580.0
                  ssoil%ssdn(j,3)  = 600.0
               ENDIF
            ENDIF
               
            Ssoil%sdepth(j,1) = ssoil%t_snwlr(j)
            
            ssoil%smass(j,1)  =  ssoil%t_snwlr(j) * ssoil%ssdn(j,1)
            
            ssoil%smass(j,2)  = ( ssoil%snowd(j) - ssoil%smass(j,1) ) * 0.4
            ssoil%smass(j,3)  = ( ssoil%snowd(j) - ssoil%smass(j,1) ) * 0.6
            
            ssoil%sdepth(j,2) = ssoil%smass(j,2) / ssoil%ssdn(j,2)
            ssoil%sdepth(j,3) = ssoil%smass(j,3) / ssoil%ssdn(j,3)
            
            ssoil%ssdnn(j) = ( ssoil%ssdn(j,1) * ssoil%smass(j,1) +            &
                              ssoil%ssdn(j,2) * ssoil%smass(j,2) +             &
                              ssoil%ssdn(j,3) * ssoil%smass(j,3) )             &
                              / ssoil%snowd(j)
         
         ENDIF 
         
         ssoil%isflag(j) = 1
      
      ENDIF ! END: IF( ssoil%snowd(j) <= 0.0 ) THEN

               
   ENDDO ! END: DO j=1,mp

END SUBROUTINE snowcheck 



SUBROUTINE snowl_adjust(dels, ssoil, canopy )
   
   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(soil_snow_type), INTENT(INOUT) :: ssoil
   TYPE(canopy_type), INTENT(INOUT)    :: canopy
  
   INTEGER :: k
   
   REAL(r_2), DIMENSION(mp) ::                                                 &
      excd,    & !
      excm,    & !
      frac,    & ! 
      xfrac     ! 
   
   REAL, DIMENSION(mp) :: osm

   INTEGER :: api ! active patch counter

    
   ! adjust levels in the snowpack due to snow accumulation/melting,
   ! snow aging etc...
   WHERE( ssoil%isflag > 0 )

      WHERE( ssoil%sdepth(:,1) > ssoil%t_snwlr )

         excd = ssoil%sdepth(:,1) - ssoil%t_snwlr
         excm = excd * ssoil%ssdn(:,1)
         ssoil%sdepth(:,1) = ssoil%sdepth(:,1) - REAL(excd)
         osm = ssoil%smass(:,1)
         ssoil%smass(:,1) = ssoil%smass(:,1) - REAL(excm)
            
         osm = ssoil%smass(:,2)
         ssoil%smass(:,2) = MAX( 0.01, ssoil%smass(:,2) + REAL(excm) )
         
         ssoil%ssdn(:,2) = REAL( MAX( 120.0_r_2, MIN( REAL( max_ssdn, r_2 ),   &
                           ssoil%ssdn(:,2) * osmi / ssoil%smass(:,2) +         &
                           ssoil%ssdn(:,1) * excm / ssoil%smass(:,2) ) ) )

          ssoil%sdepth(:,2) =  ssoil%smass(:,2) / ssoil%ssdn(:,2) 
         
          ssoil%tggsn(:,2) = REAL( ssoil%tggsn(:,2) * osm / ssoil%smass(:,2)   &
                             + ssoil%tggsn(:,1) * excm / ssoil%smass(:,2) )

          ! following line changed to fix -ve sdepth (EK 21Dec2007)
          ssoil%smass(:,3) = MAX( 0.01, ssoil%snowd - ssoil%smass(:,1)         &
                             - ssoil%smass(:,2) )

      ELSEWHERE ! ssoil%sdepth(:,1) < ssoil%t_snwlr
      
         ! 1st layer
         excd = ssoil%t_snwlr - ssoil%sdepth(:,1)
         excm = excd * ssoil%ssdn(:,2)
         osm = ssoil%smass(:,1)
         ssoil%smass(:,1) = ssoil%smass(:,1) + REAL(excm)
         ssoil%sdepth(:,1) = ssoil%t_snwlr
         ssoil%ssdn(:,1) = REAL( MAX( 120.0_r_2, MIN( REAL( max_ssdn,r_2 ),    &
                           ssoil%ssdn(:,1) * osm / ssoil%smass(:,1)            &
                           + ssoil%ssdn(:,2) * excm / ssoil%smass(:,1) ) ) )

         ssoil%tggsn(:,1) = REAL( ssoil%tggsn(:,1) * osm / ssoil%smass(:,1)   &
                          + ssoil%tggsn(:,2) * excm / ssoil%smass(:,1) )
          
         ! 2nd layer
         ssoil%smass(:,2) = MAX( 0.01, ssoil%smass(:,2) - REAL(excm) )
         ssoil%sdepth(:,2) = ssoil%smass(:,2) / ssoil%ssdn(:,2)

         ! following line changed to fix -ve sdepth (EK 21Dec2007)
         ssoil%smass(:,3) = MAX( 0.01, ssoil%snowd - ssoil%smass(:,1)          &
                            - ssoil%smass(:,2) )

      END WHERE
   
   END WHERE 

   DO  api=1,mp
      
      IF( ssoil%isflag(api).GT.0 ) THEN
      
         frac(api) = ssoil%smass(api,2) / MAX( 0.02, ssoil%smass(api,3) )
         ! if frac > 0.6 or frac < 0.74 do nothing 
         ! HOW TO translate this to xfrac
         xfrac(api) = 2.0/3.0/ frac(api)
         
         IF( xfrac(api) > 1.0 ) THEN

            excm(api) = (xfrac(api) - 1.0) * ssoil%smass(api,2)
            osm(api) = ssoil%smass(api,2)
            
            ! changed 0.02 to 0.01 to fix -ve sdepth (EK 21Dec2007)
            ssoil%smass(api,2) = MAX( 0.01, ssoil%smass(api,2) +               &
                                 REAL( excm(api) ) )

            ssoil%tggsn(api,2) = ssoil%tggsn(api,2) * osm(api) /               &
                                 ssoil%smass(api,2) +  ssoil%tggsn(api,3)      &
                                 * REAL( excm(api) )/ ssoil%smass(api,2)

            ssoil%ssdn(api,2) = MAX( 120.0, MIN( max_ssdn, ssoil%ssdn(api,2) * &
                                osm(api) / ssoil%smass(api,2) +                &
                                ssoil%ssdn(api,3) * REAL( excm(api) )          &
                                / ssoil%smass(api,2) ) )

            ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
            ssoil%smass(api,3) = MAX( 0.01, ssoil%snowd(api) -                 &
                                 ssoil%smass(api,1) - ssoil%smass(api,2) )
            
            ssoil%sdepth(api,3) = MAX( 0.02, ssoil%smass(api,3) /              &
                                  ssoil%ssdn(api,3) )
            
         ELSE! xfrac < 1
           
            excm(api) = ( 1 - xfrac(api) ) * ssoil%smass(api,2)
            ssoil%smass(api,2) = MAX(0.01, ssoil%smass(api,2) - REAL(excm(api),r_1))
           ssoil%sdepth(api,2) = MAX(0.02, ssoil%smass(api,2) / ssoil%ssdn(api,2) )

           osm(api) = ssoil%smass(api,3)
           ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
           ssoil%smass(api,3) = MAX(0.01, &
                         & ssoil%snowd(api) - ssoil%smass(api,1) - ssoil%smass(api,2))


           ssoil%tggsn(api,3) = ssoil%tggsn(api,3) * osm(api) / ssoil%smass(api,3) +  &
                         & ssoil%tggsn(api,2) * REAL(excm(api),r_1)/ ssoil%smass(api,3)
           ssoil%ssdn(api,3) = MAX(120.0, MIN(max_ssdn, ssoil%ssdn(api,3)* &
                osm(api)/ssoil%smass(api,3) + ssoil%ssdn(api,2) &
                * REAL(excm(api),r_1) / ssoil%smass(api,3)) )
           ssoil%sdepth(api,3) = ssoil%smass(api,3) /  ssoil%ssdn(api,3)
        END IF
        ssoil%isflag(api) = 1
        ssoil%ssdnn(api) = (ssoil%ssdn(api,1) * ssoil%sdepth(api,1) + ssoil%ssdn(api,2) &
             & * ssoil%sdepth(api,2) + ssoil%ssdn(api,3) * ssoil%sdepth(api,3) ) &
             & / (ssoil%sdepth(api,1) + ssoil%sdepth(api,2) + ssoil%sdepth(api,3))
      END IF
   END DO

  END SUBROUTINE snowl_adjust



  SUBROUTINE soilfreeze(dels, soil, ssoil)
   use cable_common_module
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    REAL(r_2), DIMENSION(mp)           :: sicefreeze
    REAL(r_2), DIMENSION(mp)           :: sicemelt
    REAL, DIMENSION(mp)           :: xx
    INTEGER k

    xx = 0.
    DO k = 1, ms
       WHERE (ssoil%tgg(:,k) < tfrz &
            & .AND. frozen_limit * ssoil%wb(:,k) - ssoil%wbice(:,k) > .001)
          sicefreeze = MIN( MAX(0.0_r_2,(frozen_limit*ssoil%wb(:,k)-ssoil%wbice(:,k))) &
               & * soil%zse(k) * 1000.0, &
               & (tfrz - ssoil%tgg(:,k) ) * ssoil%gammzz(:,k) / hlf )
          ssoil%wbice(:,k) = MIN( ssoil%wbice(:,k) + sicefreeze / (soil%zse(k) &
               * 1000.0), frozen_limit * ssoil%wb(:,k) )
          xx=soil%css * soil%rhosoil
          ssoil%gammzz(:,k) = MAX( &
               REAL((1.0 - soil%ssat) * soil%css * soil%rhosoil ,r_2) &
               + (ssoil%wb(:,k) - ssoil%wbice(:,k)) * REAL(cswat * rhowat,r_2) &
               + ssoil%wbice(:,k) * REAL(csice * rhowat * 0.9,r_2), REAL(xx,r_2)) &
               * REAL(soil%zse(k),r_2)
          WHERE (k == 1 .AND. ssoil%isflag == 0)
             ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
          END WHERE
          ssoil%tgg(:,k) = ssoil%tgg(:,k) + REAL(sicefreeze,r_1) &
               * hlf / REAL(ssoil%gammzz(:,k),r_1)
       ELSEWHERE (ssoil%tgg(:,k) > tfrz .AND. ssoil%wbice(:,k) > 0.)
          sicemelt = MIN(ssoil%wbice(:,k) * soil%zse(k) * 1000.0, &
               & (ssoil%tgg(:,k) - tfrz) * ssoil%gammzz(:,k) / hlf)
          ssoil%wbice(:,k) = MAX(0.0_r_2, ssoil%wbice(:,k) - sicemelt &
               / (soil%zse(k) * 1000.0) )
          xx = soil%css * soil%rhosoil
          ssoil%gammzz(:,k) = MAX( &
               REAL((1.0-soil%ssat) * soil%css * soil%rhosoil,r_2) &
               + (ssoil%wb(:,k) - ssoil%wbice(:,k)) * REAL(cswat*rhowat,r_2) &
               + ssoil%wbice(:,k) * REAL(csice * rhowat * 0.9,r_2), &
               REAL(xx,r_2) ) * REAL(soil%zse(k),r_2)
          WHERE (k == 1 .AND. ssoil%isflag == 0)
             ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
          END WHERE
          ssoil%tgg(:,k) = ssoil%tgg(:,k) - REAL(sicemelt,r_1) &
               * hlf / REAL(ssoil%gammzz(:,k),r_1)
       END WHERE
    END DO
  END SUBROUTINE soilfreeze


  !******************************************************************
  SUBROUTINE remove_trans(dels, soil, ssoil, canopy, veg)
    ! Removes transpiration water from soil.
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    REAL(r_2), DIMENSION(mp,0:ms) :: diff 
    REAL(r_2), DIMENSION(mp)      :: xx,xxd,evap_cur
    INTEGER k


     xx = 0.; xxd = 0.; diff(:,:) = 0.
     DO k = 1,ms
     ! Removing transpiration from soil:
       WHERE (canopy%fevc > 0.0 )     ! convert to mm/dels
         ! Calculate the amount (perhaps moisture/ice limited)
         ! which can be removed:
         xx = canopy%fevc * dels / hl * veg%froot(:,k) + diff(:,k-1)   ! kg/m2
         diff(:,k) = MAX( 0.0, ssoil%wb(:,k) - soil%swilt) &      ! m3/m3
                   & * soil%zse(k)*1000.0
         xxd = xx - diff(:,k)
         
         WHERE ( xxd .gt. 0.0 )
           ssoil%wb(:,k) = ssoil%wb(:,k) - diff(:,k) / (soil%zse(k)*1000.0)
           diff(:,k) = xxd
         ELSEWHERE
           ssoil%wb(:,k) = ssoil%wb(:,k) - xx / (soil%zse(k)*1000.0)
           diff(:,k) = 0.0
         ENDWHERE
       END WHERE
     END DO

  END SUBROUTINE remove_trans 



  ! Inputs:
  !	 dt_in - time step in sec
  !	 ktau_in - time step no.
  !	 ga	 - ground heat flux W/m^2
  !	 dgdtg	 -
  !	 condxpr - total precip reaching the ground (liquid and solid)
  !	 scondxpr - precip (solid only)
  !	 fev   - transpiration (W/m2)
  !	 fes   - soil evaporation (W/m2)
  !	 isoil - soil type
  !	 ivegt - vegetation type
  ! Output
  !	 ssoil
  SUBROUTINE soil_snow(dels, soil, ssoil, canopy, met, bal, veg)
   use cable_common_module
!  use arraydiag_m
    REAL, INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
    TYPE (balances_type), INTENT(INOUT)      :: bal
    INTEGER, PARAMETER  :: ntest = 0 !  for prints
    INTEGER             :: k
    REAL, DIMENSION(mp) :: snowmlt
    REAL, DIMENSION(mp) :: totwet
    REAL, DIMENSION(mp) :: weting
    REAL, DIMENSION(mp) :: xxx, tgg_old, tggsn_old
    REAL(r_2), DIMENSION(mp) :: xx,deltat,sinfil1,sinfil2,sinfil3 
    REAL                :: zsetot
    logical :: cable_runtime_coupled = .false.
    integer, save :: ktau =0 

   ktau = ktau +1 

   if( cable_runtime%um) then
      max_glacier_snowd = 50000.0
   else
      max_glacier_snowd = 1100.0
   endif


    zsetot = sum(soil%zse) 
    ssoil%tggav = 0.
    DO k = 1, ms
     ssoil%tggav = ssoil%tggav  + soil%zse(k)*ssoil%tgg(:,k)/zsetot
    END DO


   if( cable_runtime%offline .or. cable_runtime%mk3l ) then
        ssoil%t_snwlr = 0.05
   endif

    ssoil%fwtop1 = 0.0
    ssoil%fwtop2 = 0.0
    ssoil%fwtop3 = 0.0
    ssoil%runoff = 0.0 ! initialise total runoff
    ssoil%rnof1 = 0.0 ! initialise surface runoff
    ssoil%rnof2 = 0.0 ! initialise deep drainage
    ssoil%smelt = 0.0 ! initialise snowmelt
    ssoil%dtmlt = 0.0 
    ssoil%osnowd = ssoil%snowd


if(.NOT.cable_runtime_coupled) then

    IF (ktau_gl <= 1) THEN
      canopy%dgdtg = 0.0
      ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
      ssoil%wbtot = 0.0
      DO k = 1, ms
        ssoil%wb(:,k)  = min( soil%ssat,max ( ssoil%wb(:,k), soil%swilt ))
      END DO

      ssoil%wb(:,ms-2)  = min( soil%ssat,max ( ssoil%wb(:,ms-2), 0.5*(soil%sfc+soil%swilt) ))
      ssoil%wb(:,ms-1)  = min( soil%ssat,max ( ssoil%wb(:,ms-1), 0.8*soil%sfc ))
      ssoil%wb(:,ms)    = min( soil%ssat,max ( ssoil%wb(:,ms),   soil%sfc) )
      
      DO k = 1, ms
        WHERE (ssoil%tgg(:,k) <= tfrz .and. ssoil%wbice(:,k) <= 0.01)
          ssoil%wbice(:,k) = 0.5 * ssoil%wb(:,k)
        END WHERE
        WHERE (ssoil%tgg(:,k) < tfrz)
          ssoil%wbice(:,k) = frozen_limit * ssoil%wb(:,k)
        END WHERE
      END DO

      WHERE (soil%isoilm == 9) 
        ssoil%snowd = max_glacier_snowd
        ssoil%osnowd = max_glacier_snowd
        ssoil%tgg(:,1) = ssoil%tgg(:,1) - 1.0
        ssoil%wb(:,1) = 0.95 * soil%ssat
        ssoil%wb(:,2) = 0.95 * soil%ssat
        ssoil%wb(:,3) = 0.95 * soil%ssat
        ssoil%wb(:,4) = 0.95 * soil%ssat
        ssoil%wb(:,5) = 0.95 * soil%ssat
        ssoil%wb(:,6) = 0.95 * soil%ssat
        ssoil%wbice(:,1) = 0.90 * ssoil%wb(:,1)
        ssoil%wbice(:,2) = 0.90 * ssoil%wb(:,2)
        ssoil%wbice(:,3) = 0.90 * ssoil%wb(:,3)
        ssoil%wbice(:,4) = 0.90 * ssoil%wb(:,4)
        ssoil%wbice(:,5) = 0.90 * ssoil%wb(:,5)
        ssoil%wbice(:,6) = 0.90 * ssoil%wb(:,6)
      ENDWHERE
      
      xx=soil%css * soil%rhosoil
      ssoil%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
           & + (ssoil%wb(:,1) - ssoil%wbice(:,1) ) * cswat * rhowat &
           & + ssoil%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1)
    END IF
endif  ! if(.NOT.cable_runtime_coupled)

    xx=soil%css * soil%rhosoil
    IF (ktau <= 1) ssoil%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
             & + (ssoil%wb(:,1) - ssoil%wbice(:,1) ) * cswat * rhowat &
             & + ssoil%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1) + &
             & (1. - ssoil%isflag) * cgsnow * ssoil%snowd



    DO k = 1, ms ! for stempv
       ! Set liquid soil water fraction (fraction of saturation value):
       ssoil%wblf(:,k) = MAX( 0.01_r_2, (ssoil%wb(:,k) - ssoil%wbice(:,k)) ) &
            & / REAL(soil%ssat,r_2)
       ! Set ice soil water fraction (fraction of saturation value):
       ssoil%wbfice(:,k) = REAL(ssoil%wbice(:,k),r_1) / soil%ssat
    END DO
  
    CALL snowcheck (dels, ssoil, soil, met )

    CALL snowdensity (dels, ssoil, soil)

    CALL snow_accum (dels, canopy, met, ssoil, soil )

    CALL snow_melting (dels, snowmlt, ssoil, soil )

    ! Add snow melt to global snow melt variable:
    ssoil%smelt = snowmlt

    ! Adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...
    CALL snowl_adjust(dels, ssoil, canopy )

    CALL stempv(dels, canopy, ssoil, soil)
    
    ssoil%tss =  (1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)

    CALL snow_melting (dels, snowmlt, ssoil, soil )
    
    ! Add new snow melt to global snow melt variable: 
    ssoil%smelt = ssoil%smelt + snowmlt

    CALL remove_trans(dels, soil, ssoil, canopy, veg)



    CALL  soilfreeze(dels, soil, ssoil)


    totwet = canopy%precis + ssoil%smelt 
    weting = totwet + max(0.,ssoil%pudsto - canopy%fesp/hl*dels) ! total available liquid including puddle
    xxx=soil%ssat - ssoil%wb(:,1)
   
    sinfil1 = MIN( 0.95*xxx*soil%zse(1)*rhowat, weting) !soil capacity
    xxx=soil%ssat - ssoil%wb(:,2)
    sinfil2 = MIN( 0.95*xxx*soil%zse(2)*rhowat, weting - sinfil1) !soil capacity
    xxx=soil%ssat - ssoil%wb(:,3)
    sinfil3 = MIN( 0.95*xxx*soil%zse(3)*rhowat,weting-sinfil1-sinfil2)
    ssoil%fwtop1 = sinfil1 / dels - canopy%segg          ! net water flux to the soil
    ssoil%fwtop2 = sinfil2 / dels           ! net water flux to the soil
    ssoil%fwtop3 = sinfil3 / dels           ! net water flux to the soil

!   Puddle for the next time step
    ssoil%pudsto = max( 0., weting - sinfil1 - sinfil2 - sinfil3 )
    ssoil%rnof1 = max(0.,ssoil%pudsto - ssoil%pudsmx)
    ssoil%pudsto = ssoil%pudsto - ssoil%rnof1

    CALL surfbv(dels, met, ssoil, soil, veg, canopy )


    canopy%fhs_cor = ssoil%dtmlt(:,1)*ssoil%dfh_dtg
    canopy%fes_cor = ssoil%dtmlt(:,1)*(ssoil%cls*ssoil%dfe_ddq * ssoil%ddq_dtg)

    canopy%fhs = canopy%fhs+canopy%fhs_cor
    canopy%fes = canopy%fes+canopy%fes_cor

    CALL hydraulic_redistribution(dels,soil,ssoil,canopy,veg, met)

    ssoil%smelt = ssoil%smelt/dels

    ! Set weighted soil/snow surface temperature
    ssoil%tss=(1-ssoil%isflag)*ssoil%tgg(:,1) + ssoil%isflag*ssoil%tggsn(:,1)

    ssoil%wbtot = 0.0
    DO k = 1, ms
       ssoil%wbtot = ssoil%wbtot + REAL(ssoil%wb(:,k)*1000.0*soil%zse(k),r_2)
    END DO

  END SUBROUTINE soil_snow

  !+++++++++++++++++++  Hydraulic Redistribution Section  ++++++++++++++++++++++
  ! Sciences from Ryel et al. Oecologia, 2002; Lee et al., 2005, PNAS
  ! Code by LiLH 16 Feb, 2011
  ! Fixed problem of negative wb in global run by BP Mar 2011
  SUBROUTINE hydraulic_redistribution(dels, soil, ssoil, canopy, veg, met)
    REAL,                 INTENT(IN) :: dels ! integration time step (s)
    TYPE(soil_parameter_type), INTENT(IN) :: soil
    TYPE(soil_snow_type),   INTENT(INOUT) :: ssoil
    TYPE(canopy_type),         INTENT(IN) :: canopy
    TYPE(veg_parameter_type),  INTENT(IN) :: veg
    TYPE(met_type), INTENT(INOUT)         :: met ! all met forcing
    INTEGER k
    INTEGER j
    INTEGER ii
    REAL, DIMENSION(mp,ms)    :: S_VG                ! --
    REAL, DIMENSION(mp,ms)    :: wpsy                ! MPa
    REAL, DIMENSION(mp)       :: frootX              ! --
    REAL, DIMENSION(mp,ms)    :: C_hr                ! --
!    REAL, DIMENSION(mp,ms)    :: hr                  ! cm/hour
    REAL, DIMENSION(mp,ms,ms) :: hr_term             ! cm/hour
    REAL, DIMENSION(mp)       :: Dtran               ! Swith for hr

    REAL, PARAMETER :: thetas=0.45  ! from Belk et al., 2007, WRR
    REAL, PARAMETER :: thetar=0.20  ! from Belk et al., 2007, WRR
    ! REAL, PARAMETER :: alpha_VG = 0.00045  ! from Belk et al., 2007, WRR,
                                                  ! cm^{-1} 1cmH2O=100Pa
    ! REAL, PARAMETER :: n_VG = 1.40         ! from Belk et al., 2007, WRR
    REAL, PARAMETER :: n_VG = 2.06           ! -- 2.06
    REAL, PARAMETER :: m_VG = 1.0-1.0/n_VG   ! --
    REAL, PARAMETER :: alpha_VG = 0.00423    ! cm^{-1} Note: 1cmH2O=100Pa
    REAL, PARAMETER :: n_hr = 3.22           ! --
    REAL, PARAMETER :: wpsy50 = -1.0         ! MPa
    REAL, PARAMETER :: CRT = 125.0           ! cm MPa^-1 h^-1, default value (0.097) from Ryel et al., 2002
    REAL, PARAMETER :: wiltParam = 0.5   
    REAL, PARAMETER :: satuParam = 0.8   
    REAL, DIMENSION(mp,ms,ms) :: hr_perTime
    REAL, DIMENSION(mp)       :: temp
    REAL, DIMENSION(mp)       :: available
    REAL, DIMENSION(mp)       :: accommodate
    REAL, DIMENSION(mp)       :: totalmoist,totalice
    REAL, DIMENSION(mp)       :: total2,zsetot

    zsetot = sum(soil%zse)
    totalmoist(:) = 0.0
    totalice(:) = 0.0
    DO k=1, ms
      totalmoist(:) = totalmoist(:) + ssoil%wb(:,k)*soil%zse(k)/zsetot
      totalice(:) = totalice(:) + ssoil%wbice(:,k)*soil%zse(k)/zsetot
    ENDDO

    Dtran=0.0
    WHERE ( canopy%fevc < 10.0 .and.  totalice  < 1.e-2)   Dtran=1.0
    
    DO k=1, ms
      S_VG(:,k) = MIN( 1.0, MAX(1.0E-4, ssoil%wb(:,k) - soil%swilt) &
                             / (soil%ssat - soil%swilt) )
      wpsy(:,k) = -1.0/alpha_VG*(S_VG(:,k)**(-1.0/m_VG)-1.0)**(1/n_VG) &
                   *100*1.0E-6     ! VG model, convert from cm to Pa by (*100),
                                   !                           to MPa (*1.0E-6)
      C_hr(:,k) = 1./(1+(wpsy(:,k)/wpsy50)**n_hr)
    ENDDO

    temp(:)        = 0.0
    hr_term(:,:,:) = 0.0    ! unit: cm h^{-1}
    hr_perTime(:,:,:) = 0.0
    ! setting hr_term=0 for top layer, follows Lee et al., 2005, PNAS
    DO k = ms, 3, -1
      DO j = k-1, 2, -1
        temp(:)        = 0.0
        available(:)   = 0.0
        accommodate(:) = 0.0
        frootX= max(0.01,max( veg%froot(:,k),veg%froot(:,j)))
        hr_term(:,k,j) = CRT*(wpsy(:,j)-wpsy(:,k))*MAX(C_hr(:,k),C_hr(:,j)) &
                       *(veg%froot(:,k)*veg%froot(:,j))/(1-frootX) * Dtran
        hr_perTime(:,k,j) = hr_term(:,k,j)*1.0E-2/3600.0*dels ! m per timestep
        hr_perTime(:,j,k) = -1.0 * hr_perTime(:,k,j)
        hr_perTime(:,k,j) = hr_perTime(:,k,j)/soil%zse(k)
        hr_perTime(:,j,k) = hr_perTime(:,j,k)/soil%zse(j)
        ! Restricting changes to all broadleaf forests, and
        ! other forests and woody savannas in the tropics
        ! Note that veg types here are based on IGBP classification (BP mar2011)
        WHERE (.NOT.(veg%iveg == 2 .OR. veg%iveg == 7 ))
          hr_perTime(:,k,j) = 0.0
          hr_perTime(:,j,k) = 0.0
        ENDWHERE
        WHERE (hr_perTime(:,k,j) < 0.0)

          available(:)   = MAX(0.0, ssoil%wb(:,k)-  &
                          ( soil%swilt(:) + (soil%sfc(:)-soil%swilt(:))/3.) )
          accommodate(:) = MAX(0.0, soil%ssat(:)-ssoil%wb(:,j))
          temp(:) = MAX(hr_perTime(:,k,j), &
                        -1.0*wiltParam*available(:), &
                        -1.0*satuParam*accommodate(:)*soil%zse(j)/soil%zse(k)) 
          hr_perTime(:,k,j) = temp(:)
          hr_perTime(:,j,k) = -1.0 * temp(:) * soil%zse(k) / soil%zse(j)
        ELSEWHERE (hr_perTime(:,j,k) < 0.0)

          available(:)   = MAX(0.0, ssoil%wb(:,j)-  &
                           ( soil%swilt(:) + (soil%sfc(:)-soil%swilt(:))/3.) )
          accommodate(:) = MAX(0.0, soil%ssat(:)-ssoil%wb(:,k))
          temp(:) = MAX(hr_perTime(:,j,k), &
                        -1.0*wiltParam*available(:), &
                        -1.0*satuParam*accommodate(:)*soil%zse(k)/soil%zse(j))
          hr_perTime(:,j,k) = temp(:)
          hr_perTime(:,k,j) = -1.0 * temp(:) * soil%zse(j) / soil%zse(k)
        ENDWHERE
        ssoil%wb(:,k) = ssoil%wb(:,k) + hr_perTime(:,k,j)
        ssoil%wb(:,j) = ssoil%wb(:,j) + hr_perTime(:,j,k)
      ENDDO 
    ENDDO

    WHERE ( met%tk < tfrz + 5.  ) Dtran=0.0
      DO k=1, ms
        S_VG(:,k) = MIN( 1.0, MAX(1.0E-4, ssoil%wb(:,k) - soil%swilt) &
                    / (soil%ssat - soil%swilt) )
        wpsy(:,k) = -1.0/alpha_VG*(S_VG(:,k)**(-1.0/m_VG)-1.0)**(1/n_VG) &
                     *100*1.0E-6     ! VG model, convert from cm to Pa by (*100),
           !                           to MPa (*1.0E-6)
        C_hr(:,k) = 1./(1+(wpsy(:,k)/wpsy50)**n_hr)
      ENDDO                                                                                                   

      hr_term(:,:,:) = 0.0    ! unit: cm h^{-1}
      hr_perTime(:,:,:) = 0.0  
      DO k = 1,ms-2
        DO j = k+1,ms-1
         temp(:)        = 0.0
         available(:)   = 0.0
         accommodate(:) = 0.0
         frootX= max(0.01,max( veg%froot(:,k),veg%froot(:,j)))
         hr_term(:,k,j) = CRT*(wpsy(:,j)-wpsy(:,k))*MAX(C_hr(:,k),C_hr(:,j)) &
            *(max(0.01,veg%froot(:,k))*max(0.01,veg%froot(:,j)))/(1-frootX)*Dtran
         hr_perTime(:,k,j) = hr_term(:,k,j)*1.0E-2/3600.0*dels ! m per timestep
         hr_perTime(:,j,k) = -1.0 * hr_perTime(:,k,j)
         hr_perTime(:,k,j) = hr_perTime(:,k,j)/soil%zse(k)
         hr_perTime(:,j,k) = hr_perTime(:,j,k)/soil%zse(j)
        ! Restricting changes to all broadleaf forests, and
        ! other forests and woody savannas in the tropics
        ! Note that veg types here are based on IGBP classification (BP mar2011)
        !        WHERE (.NOT.(veg%iveg == 1 .OR. veg%iveg == 6 ))
        WHERE (.NOT.(veg%iveg == 2 .OR. veg%iveg == 7 ))
           hr_perTime(:,k,j) = 0.0
           hr_perTime(:,j,k) = 0.0
        ENDWHERE
        WHERE (hr_perTime(:,k,j) < 0.0)
           available(:)   = MAX(0.0, ssoil%wb(:,k)- soil%sfc(:))
           accommodate(:) = MAX(0.0, soil%ssat(:)-ssoil%wb(:,j))
           temp(:) = MAX(hr_perTime(:,k,j), &
                         -1.0*wiltParam*available(:), &
                         -1.0*satuParam*accommodate(:)*soil%zse(j)/soil%zse(k))
           hr_perTime(:,k,j) = temp(:)
           hr_perTime(:,j,k) = -1.0 * temp(:) * soil%zse(k) / soil%zse(j)
        ELSEWHERE (hr_perTime(:,j,k) < 0.0)
           available(:)   = MAX(0.0, ssoil%wb(:,j)- soil%sfc(:))
           accommodate(:) = MAX(0.0, soil%ssat(:)-ssoil%wb(:,k))
           temp(:) = MAX(hr_perTime(:,j,k), &
                    -1.0*wiltParam*available(:), &
                    -1.0*satuParam*accommodate(:)*soil%zse(k)/soil%zse(j))
           hr_perTime(:,j,k) = temp(:)
           hr_perTime(:,k,j) = -1.0 * temp(:) * soil%zse(j) / soil%zse(k)
        ENDWHERE
        ssoil%wb(:,k) = ssoil%wb(:,k) + hr_perTime(:,k,j)
        ssoil%wb(:,j) = ssoil%wb(:,j) + hr_perTime(:,j,k)
      ENDDO
     ENDDO
                           

  END SUBROUTINE hydraulic_redistribution
  !+++++++++++++++++++  Hydraulic Redistribution Section  ++++++++++++++++++++++


END MODULE soil_snow_module
