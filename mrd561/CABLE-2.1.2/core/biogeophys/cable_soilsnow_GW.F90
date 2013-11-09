!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: All routines for calculating soil temperature and moisture
!          and snow calculations
!
! Contact: Eva.Kowalczyk@csiro.au
!
! History: v2.0 Tighter water budget
!          v2.0 Hydraulic redistribution subroutine (with namelist switch). 
!               NB Currently hard-wired to veg types 2 and 7 
!                  (usually evergreen broadleaf and c4 grass)
!          v2.0 ssoil variable renamed ssnow
!          Mark Decker - used ssnow as base for ssgw.  Could be part of same module
!
! ==============================================================================

MODULE cable_soil_snow_gw_module
   
   USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
                             veg_parameter_type, canopy_type, met_type,        &
                             balances_type, r_2, ms, mp           
   USE cable_data_module, ONLY : issnow_type, point2constants

   IMPLICIT NONE

   PRIVATE

   TYPE ( issnow_type ) :: C 
   
   REAL, PARAMETER ::                                                          &
      cgsnow = 2090.0,     & ! specific heat capacity for snow
      csice = 2.100e3,     & ! specific heat capacity for ice
      cswat = 4.218e3,     & ! specific heat capacity for water
      rhowat = 1000.0,     & ! density of water
      snmin = 1.,          & ! for 3-layer;
      max_ssdn = 750.0,    & !
      max_sconds = 2.51,   & !
      frozen_limit = 0.85    ! EAK Feb2011 (could be 0.95) 
      
   !MD GW params
   REAL, PARAMETER :: sucmin  = -10000000.0,  & ! minimum soil pressure head [mm]
                      qhmax   = 1e-6,         & ! max horizontal drainage [mm/s]
                      hkrz    = 2.0,          & ! GW_hksat e-folding depth [mm**-1]
                      volwatmin  = 0.05,         & !min soil water [mm]      
                      wtd_uncert = 5.0,       &  ! uncertaintiy in wtd calcultations [mm]
                      wtd_max = 100000.0,     & ! maximum wtd [mm]
                      wtd_min = 10.0,         & ! minimum wtd [mm]
                      denliq = 1000.0,        & ! denisty of liquid water [kg/m3]
                      denice = 917.0            !denisty of ice
                      
   INTEGER, PARAMETER :: mx_wtd_iterations = 25 ! maximum number of iterations to find the water table depth                    
  
   
   REAL :: cp    ! specific heat capacity for air
   
   !jhan:make parameter
   REAL :: max_glacier_snowd
 
   ! This module contains the following subroutines:
   PUBLIC soil_snow_gw ! must be available outside this module
   PRIVATE snowdensity, snow_melting, snowcheck, snowl_adjust 
   PRIVATE trimb, smoisturev, snow_accum, stempv
   PRIVATE soilfreeze, remove_trans
   PRIVATE smoistgw, ovrlndflx, solve_tridiag

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
   rhs(:,kmax) = ( rhs(:,kmax) - a(:,kmax) * g(:,kmax-1) )                     &
                 / ( b(:,kmax) - a(:,kmax) * e(:,kmax-1) )
                 
   DO k = kmax - 1, 1, - 1
     rhs(:,k) = g(:,k) - e(:,k) * rhs(:,k + 1)
   END DO
  
END SUBROUTINE trimb

! -----------------------------------------------------------------------------

! SUBROUTINE smoisturev (fwtop,dels,ssnow,soil)
!      Solves implicit soil moisture equation
!      Science development by Eva Kowalczyk and John McGregor, CMAR
!
SUBROUTINE smoisturev (dels,ssnow,soil,veg)
   
   USE cable_common_module
   
   REAL, INTENT(IN) :: dels    ! time step size (s)
   
   TYPE(soil_snow_type),      INTENT(INOUT) ::                                &
      ssnow ! soil and snow variables
 
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
         wbl_k = MAX( 0.01_r_2, ssnow%wb(:,k) - ssnow%wbice(:,k) )
         wbl_kp = MAX( 0.01_r_2, ssnow%wb(:,k+1) - ssnow%wbice(:,k+1) )
         
         ! Calculate difference in liq soil water b/w consecutive layers:
         delt(:,k) = wbl_kp - wbl_k
         
         ! especially to allow for isolated frozen layers, use min speed
         wh = MIN( wbl_k, wbl_kp )
         WHERE( ssnow%wbice(:,k) > 0.05 .OR. ssnow%wbice(:,k+1) > 0.01 )       &
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

      WHERE( ssnow%wb(:,ms) > soil%sfc(:) )

         wbl_k = MAX( 0.001_r_2, ssnow%wb(:,ms) - ssnow%wbice(:,ms) )
         wbl_kp = MAX( 0.001_r_2, soil%ssat(:) - ssnow%wbice(:,ms) )
         
         wh = MIN( wbl_k, wbl_kp )
         
         WHERE( ssnow%wbice(:,ms) .GT. 0.05 ) wh = 0.9 * wbl_k + 0.1 * wbl_kp
       
         ! Calculate hyd conductivity adjusted for ice:
         hydss = soil%hyds
    
         speed_k = hydss * ( wh / soil%ssat )**( soil%i2bp3 - 1 )
         speed_k =  0.5 * speed_k / ( 1. - MIN( 0.5, 10. * ssnow%wbice(:,ms) ) )
         fluxlo = wbl_k
         
         ! scale speed to grid lengths per dt & limit speed for stability
         speed_k = MIN( 0.5 * speed_k, 0.5 * soil%zse(ms) / dels )
         fluxh(:,ms) = MAX( 0.0, speed_k * fluxlo )
     
      END WHERE

      ! update wb by TVD method
      DO k = ms, 1, -1
        
         IF(  nmeth == -1 ) THEN ! each new wb constrained by ssat
            fluxh(:,k-1) = MIN( fluxh (:,k-1), ( soil%ssat - ssnow%wb(:,k) )   &
                           * soil%zse(k) / dels + fluxh(:,k) )
         END IF

         ! fluxh (:,ms) is drainage
         ssnow%wb(:,k) = ssnow%wb(:,k) + dels * ( fluxh(:,k-1) - fluxh(:,k) )  &
                         / soil%zse(k)

         ! re-calculate wblf
         ssatcurr_k = soil%ssat - ssnow%wbice(:,k)
         dtt(:,k) = dels / ( soil%zse(k) * ssatcurr_k )

         ! this defn of wblf has different meaning from previous one in surfbv
         ! N.B. are imposing wbice<wb, so wblf <1
         ssnow%wblf(:,k) = ( ssnow%wb(:,k) - ssnow%wbice(:,k) ) / ssatcurr_k
      
      END DO
      
      ssnow%rnof2 = dels * REAL( fluxh(:,ms) ) * 1000.0

      ! wbh_k represents wblf(k-.5)
      DO k = 2, ms
         
         ssatcurr_k = REAL( soil%ssat, r_2 ) - ssnow%wbice(:,k)
         wbh_k = ( soil%zse(k) * ssnow%wblf(:,k-1) + soil%zse(k-1)             &
                 * ssnow%wblf(:,k) ) / ( soil%zse(k) + soil%zse(k-1) )
         ! i.e. wbh**(bch+1)
         fact = wbh_k**( soil%ibp2 - 1 )

         ! with 50% wbice, reduce hbsh by 1.e-5
         pwb_wbh = ( soil%hsbh * ( 1. - MIN( 2. * MIN (0.1_r_2, MAX(           &
                  ssnow%wbice(:,k-1) / MAX( 0.01_r_2, ssnow%wb(:,k-1) ),       &
                  ssnow%wbice(:,k)   / MAX( 0.01_r_2, ssnow%wb(:,k) ) ) )      &
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
      ssnow%wblf(:,1) = ssnow%wblf(:,1) + dtt(:,1) * ssnow%fwtop1 / rhowat
      ssnow%wblf(:,2) = ssnow%wblf(:,2) + dtt(:,2) * ssnow%fwtop2 / rhowat
      ssnow%wblf(:,3) = ssnow%wblf(:,3) + dtt(:,3) * ssnow%fwtop3 / rhowat
    
   END IF
   ! END: IF (nmeth <= 0) THEN

   IF ( nmeth > 0 ) THEN
      
      wbficemx = 0.0
      
      DO k = 1, ms
         
         ssatcurr(:,k) = REAL(soil%ssat,r_2) - ssnow%wbice(:,k)
         
         ! this defn of wblf has different meaning from previous one in surfbv
         ! N.B. are imposing wbice<wb, so wblf <1
         ssnow%wblf(:,k) = ( ssnow%wb(:,k) - ssnow%wbice(:,k) ) / ssatcurr(:,k)
         
         ssnow%wbfice(:,k) = REAL( ssnow%wbice(:,k) ) / soil%ssat
         
         wbficemx = MAX( wbficemx, ssnow%wbfice(:,k) )
         dtt(:,k) = dels / ( soil%zse(k) * ssatcurr(:,k) )
      
      END DO

      IF( nmeth == 1 ) THEN ! full implicit method
        
         DO k = 2, ms
           
            wbh(:,k) = ( soil%zse(k) * ssnow%wblf(:,k-1) + soil%zse(k-1)       &
                       * ssnow%wblf(:,k) ) / ( soil%zse(k) + soil%zse(k-1) )
                       
            fact = wbh(:,k)**( soil%ibp2 - 1 ) ! i.e. wbh**(bch+1)
            fact2 = fact * fact
            pwb = soil%hsbh * fact
           
            ! moisture diffusivity (D) is  wbh*pwb
            ! other term (K) is wbh*soil%hyds*fact2
            z1(:,k) = wbh(:,k) * ( (soil%i2bp3 - 1 ) * soil%hyds * fact2       &
                      - soil%ibp2 * pwb *                                      &
                      ( ssnow%wblf(:,k) - ssnow%wblf(:,k-1 ) ) / soil%zshh (k) )
           
            z2(:,k) = - soil%i2bp3 * soil%hyds * fact2 + soil%ibp2 * pwb       &
                      * ( ssnow%wblf(:,k) - ssnow%wblf(:,k-1) ) / soil%zshh (k)
           
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
            ssnow%wblf(:,k) = ssnow%wblf(:,k) + dtt(:,k) *                     &
                              ( z1(:,k+1) - z1(:,k) )
         END DO
     
      END IF ! (nmeth == 1)
     
      IF (nmeth >= 2) THEN ! part implicit method
         
         DO k = 2, ms
            z1mult(:,k) = soil%i2bp3 ! corresponds to 2b+3
         END DO

         DO k = 2, ms ! wbh(k) represents wblf(k-.5)
            wbh(:,k) = ( soil%zse(k) * ssnow%wblf(:,k-1) + soil%zse(k-1)       &
                       * ssnow%wblf(:,k) ) / ( soil%zse(k) + soil%zse(k-1) )
            
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
                  ssnow%wblf(:,k) = ssnow%wblf(:,k) + dtt(:,k)                 &
                                    * ( ( z1mult(:,k+1) - 1.0 ) * z1(:,k+1)    &
                                    - (z1mult(:,k) - 1.0) * z1(:,k) )

            ELSE
               
               WHERE( wbficemx < 0.75 )                                        &
                  ssnow%wblf(:,k) = ssnow%wblf(:,k) + dtt(:,k)                 &
                                    * ( z1(:,k) - z1(:,k+1) )
              
            END IF

         END DO
         
      END IF

           IF (nmeth == 3) THEN
         ! artificial fix applied here for safety (explicit nmeth only)
         DO k = 1, ms
            ssnow%wblf(:,k) = MAX( 0.0_r_2, MIN( ssnow%wblf(:,k), 1.0_r_2 ) )
         END DO
      END IF
   
      ssnow%wblf(:,1) = ssnow%wblf(:,1) + dtt(:,1) * ssnow%fwtop1 / rhowat
      ssnow%wblf(:,2) = ssnow%wblf(:,2) + dtt(:,2) * ssnow%fwtop2 / rhowat
      ssnow%wblf(:,3) = ssnow%wblf(:,3) + dtt(:,3) * ssnow%fwtop3 / rhowat
    
   END IF  ! IF (nmeth > 0)

   CALL trimb(at, bt, ct, ssnow%wblf, ms)
   
   DO k = 1, ms
      ssatcurr(:,k) = soil%ssat - ssnow%wbice(:,k)
      ssnow%wb(:,k) = ssnow%wblf(:,k) * ssatcurr(:,k) + ssnow%wbice(:,k)
      ssnow%wbice(:,k) = MIN( ssnow%wbice(:,k), frozen_limit * ssnow%wb(:,k) )
   END DO

END SUBROUTINE smoisturev

! -----------------------------------------------------------------------------

SUBROUTINE snowdensity (dels, ssnow, soil)
   
   REAL, INTENT(IN) :: dels   ! integration time step (s)

   TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow 
    
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil

   INTEGER, DIMENSION(mp,3) :: ssnow_isflag_ssdn 
   REAL, DIMENSION(mp) :: ssnow_tgg_min1
   REAL, DIMENSION(mp,3) :: dels_ssdn, ssnow_tgg_min
     
   ssnow_isflag_ssdn = SPREAD( ssnow%isflag,2,mp) 
   
   dels_ssdn = SPREAD( SPREAD( dels, 1, mp ), 2,  mp ) 
   ssnow_tgg_min1 = MIN( C%TFRZ, ssnow%tgg(:,1) )
   
   WHERE( ssnow%snowd > 0.1 .AND. ssnow%isflag == 0 )
      
      ssnow%ssdn(:,1) = MIN( max_ssdn, MAX( 120.0, ssnow%ssdn(:,1) + dels      &
                        * ssnow%ssdn(:,1) * 3.1e-6 * EXP( -0.03 * ( 273.15 -   &
                        ssnow_tgg_min1 ) - MERGE( 0.046, 0.0,                  &
                        ssnow%ssdn(:,1) >= 150.0 ) * ( ssnow%ssdn(:,1) - 150.0)&
                        ) ) )

            ssnow%ssdn(:,1) = MIN(max_ssdn,ssnow%ssdn(:,1) + dels * 9.806      &
          & * ssnow%ssdn(:,1) * 0.75 * ssnow%snowd                             &
          & / (3.0e7 * EXP(0.021 * ssnow%ssdn(:,1) + 0.081                     &
          & * (273.15 - MIN(C%TFRZ, ssnow%tgg(:,1) ) ) ) ) )

      ! permanent ice: fix hard-wired number in next version
      WHERE( soil%isoilm /= 9 ) ssnow%ssdn(:,1) = MIN( 450.0, ssnow%ssdn(:,1) )

      ssnow%sconds(:,1) = MAX( 0.2, MIN( 2.876e-6 * ssnow%ssdn(:,1)**2         &
                          + 0.074, max_sconds ) )

      ssnow%sconds(:,2) = ssnow%sconds(:,1) 
      ssnow%sconds(:,3) = ssnow%sconds(:,1) 
      
      ssnow%ssdnn = ssnow%ssdn(:,1)
      
      ssnow%ssdn(:,2) = ssnow%ssdn(:,1)
      ssnow%ssdn(:,3) = ssnow%ssdn(:,1)
    
   END WHERE
  

   WHERE (ssnow%isflag == 1)
      
      ssnow%ssdn(:,1) = ssnow%ssdn(:,1) + dels * ssnow%ssdn(:,1) * 3.1e-6      &
            * EXP( -0.03 * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,1)))            &
            - MERGE(0.046, 0.0, ssnow%ssdn(:,1) >= 150.0)                      &
            * (ssnow%ssdn(:,1) - 150.0) )
      
      ssnow%ssdn(:,2) = ssnow%ssdn(:,2) + dels * ssnow%ssdn(:,2) * 3.1e-6      &
            * EXP( -0.03 * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,2)))            &
            - MERGE(0.046, 0.0, ssnow%ssdn(:,2) >= 150.0)                      &
            * (ssnow%ssdn(:,2) - 150.0) )
      
      ssnow%ssdn(:,3) = ssnow%ssdn(:,3) + dels * ssnow%ssdn(:,3) * 3.1e-6      &
            * EXP( -0.03 * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,3)))            &
            - MERGE(0.046, 0.0, ssnow%ssdn(:,3) >= 150.0)                      &
            * (ssnow%ssdn(:,3) - 150.0) )
      
      ssnow%ssdn(:,1) = ssnow%ssdn(:,1) + dels * 9.806 * ssnow%ssdn(:,1)       &
            * ssnow%t_snwlr*ssnow%ssdn(:,1)                                    &
            / (3.0e7 * EXP(.021 * ssnow%ssdn(:,1) + 0.081                      &
            * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,1)))))
      
      ssnow%ssdn(:,2) = ssnow%ssdn(:,2) + dels * 9.806 * ssnow%ssdn(:,2)       &
            * (ssnow%t_snwlr * ssnow%ssdn(:,1) + 0.5 * ssnow%smass(:,2) )      &
            / (3.0e7 * EXP(.021 * ssnow%ssdn(:,2) + 0.081                      &
            * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,2)))))
      
      ssnow%ssdn(:,3) = ssnow%ssdn(:,3) + dels * 9.806 * ssnow%ssdn(:,3)       &
            * (ssnow%t_snwlr*ssnow%ssdn(:,1) + ssnow%smass(:,2)                &
            + 0.5*ssnow%smass(:,3))                                            &
            / (3.0e7 * EXP(.021 * ssnow%ssdn(:,3) + 0.081                      &
            * (273.15 - MIN(C%TFRZ, ssnow%tggsn(:,3)))))
      
      ssnow%sdepth(:,1) =  ssnow%smass(:,1) / ssnow%ssdn(:,1) 
      ssnow%sdepth(:,2) =  ssnow%smass(:,2) / ssnow%ssdn(:,2) 
      ssnow%sdepth(:,3) =  ssnow%smass(:,3) / ssnow%ssdn(:,3) 
      
      ssnow%ssdnn = (ssnow%ssdn(:,1) * ssnow%smass(:,1) + ssnow%ssdn(:,2)      &
            * ssnow%smass(:,2) + ssnow%ssdn(:,3) * ssnow%smass(:,3) )          &
            / ssnow%snowd
      
      ssnow%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,1) ** 2         &
                                                        & + 0.074, max_sconds) )
      ssnow%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,2) ** 2 &
                                                        & + 0.074, max_sconds) )
      ssnow%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,3) ** 2 &
                                                        & + 0.074, max_sconds) )
   END WHERE

END SUBROUTINE snowdensity

! -----------------------------------------------------------------------------

SUBROUTINE snow_melting (dels, snowmlt, ssnow, soil )

   USE cable_common_module
   
   REAL, INTENT(IN) :: dels   ! integration time step (s)
   
   REAL, DIMENSION(mp), INTENT(OUT) :: snowmlt ! snow melt   
   
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil
   TYPE(soil_snow_type), INTENT(INOUT)   :: ssnow  ! soil+snow variables
  
   INTEGER                 :: k,j 
   
   REAL, DIMENSION(mp) ::                                                      &
      osm,     & !
      sgamm,   & !
      snowflx    !
   
   REAL, DIMENSION(mp,0:3) :: smelt1

   snowmlt= 0.0
   smelt1 = 0.0
    
   DO j=1,mp  
      
      IF( ssnow%snowd(j) > 0.0 .AND. ssnow%isflag(j) == 0                      &
          .AND. ssnow%tgg(j,1) >= C%TFRZ ) THEN

         ! snow covered land
         ! following done in sflux  via  ga= ... +cls*egg + ...
         ! ** land,snow,melting
         snowflx(j) = REAL((ssnow%tgg(j,1) - C%TFRZ) * ssnow%gammzz(j,1))
         
         ! prevent snow depth going negative
         snowmlt(j) = MIN(snowflx(j) / C%HLF, ssnow%snowd(j) )
       
         ssnow%dtmlt(j,1) = ssnow%dtmlt(j,1) + snowmlt(j) * C%HLF              &
                            / ssnow%gammzz(j,1)

         ssnow%snowd(j) = ssnow%snowd(j) - snowmlt(j)
         ssnow%tgg(j,1) = REAL( ssnow%tgg(j,1) - snowmlt(j) *                  &
                          C%HLF / ssnow%gammzz(j,1) )
      ENDIF
    
   END DO

   smelt1(:,0) = 0.0
   
   DO k = 1, 3
    
      !where there is snow 
      WHERE( ssnow%snowd > 0.0 .AND. ssnow%isflag > 0 )

         sgamm = ssnow%ssdn(:,k) * cgsnow * ssnow%sdepth(:,k)
       
         ! snow melt refreezing
         snowflx = smelt1(:,k-1) * C%HLF / dels

         ssnow%tggsn(:,k) = ssnow%tggsn(:,k) + ( snowflx * dels +              &
                            smelt1(:,k-1)*cswat *( C%TFRZ-ssnow%tggsn(:,k) ) ) &
                            / ( sgamm + cswat*smelt1(:,k-1) )
         
         ! increase density due to snowmelt
         osm = ssnow%smass(:,k)
         ssnow%smass(:,k) = ssnow%smass(:,k) + smelt1(:,k-1)
         ssnow%ssdn(:,k) = MAX( 120.0, MIN( ssnow%ssdn(:,k) * osm /            &
                           ssnow%smass(:,k) + rhowat * ( 1.0 - osm /           &
                           ssnow%smass(:,k)), max_ssdn ) )

         ! permanent ice: fix hard-wired number in next version
         WHERE( soil%isoilm /= 9 )                                             &
            ssnow%ssdn(:,k) = MIN( 450.0, ssnow%ssdn(:,k) )

         ssnow%sdepth(:,k) = ssnow%smass(:,k) / ssnow%ssdn(:,k)
         
         sgamm = ssnow%smass(:,k) * cgsnow
        
         smelt1(:,k) = 0.0
         
         ! snow melting
         WHERE (ssnow%tggsn(:,k) > C%TFRZ)
           
            snowflx = ( ssnow%tggsn(:,k) - C%TFRZ ) * sgamm
           
            smelt1(:,k) = MIN( snowflx / C%HLF, 0.6 * ssnow%smass(:,k) )
            
            ssnow%dtmlt(:,k) = ssnow%dtmlt(:,k) + smelt1(:,k) * C%HLF / sgamm
            
            osm = ssnow%smass(:,k)
            
            ssnow%smass(:,k) = ssnow%smass(:,k) - smelt1(:,k)

            ssnow%tggsn(:,k) = ssnow%tggsn(:,k) - smelt1(:,k) * C%HLF / sgamm

            ssnow%sdepth(:,k) = ssnow%smass(:,k) / ssnow%ssdn(:,k)
       
         END WHERE 
         ! END snow melting
     
      END WHERE    
      ! END where there is snow 
   
   END DO
   
   WHERE( ssnow%snowd > 0.0 .AND. ssnow%isflag > 0 )
      snowmlt = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)
      ssnow%snowd = ssnow%snowd - snowmlt
   END WHERE

END SUBROUTINE snow_melting

! -----------------------------------------------------------------------------

SUBROUTINE snow_accum ( dels,  canopy, met, ssnow, soil )

USE cable_common_module

   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(canopy_type), INTENT(INOUT)         :: canopy ! vegetation variables
   TYPE(met_type), INTENT(INOUT)            :: met   ! all met forcing
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow ! soil+snow variables
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
   
   REAL, DIMENSION(mp) ::                                                      &
      osm,     & !
      sgamm,   & !
      snowmlt, & !
      xxx        !

   INTEGER             :: k

   WHERE (canopy%precis > 0.0 .and. ssnow%isflag == 0)
      ! accumulate solid part
      ssnow%snowd = MAX( ssnow%snowd + met%precip_sn, 0.0 ) 
      
      canopy%precis = canopy%precis - met%precip_sn
      
      ssnow%ssdn(:,1) = MAX( 120.0, ssnow%ssdn(:,1)                            &
                        * ssnow%osnowd / MAX( 0.01, ssnow%snowd )              &
                        + 120.0 * met%precip_sn / MAX( 0.01, ssnow%snowd ) )
      
      ssnow%ssdnn = ssnow%ssdn(:,1)
      
      WHERE( canopy%precis > 0.0 .AND. ssnow%tgg(:,1) < C%TFRZ )
         
         ssnow%snowd = MAX(ssnow%snowd + canopy%precis, 0.0)
        
         ssnow%tgg(:,1) = ssnow%tgg(:,1) + canopy%precis * C%HLF               &
                          / ( REAL( ssnow%gammzz(:,1) ) + cswat *canopy%precis )  
         ! change density due to water being added 
         ssnow%ssdn(:,1) = MIN( max_ssdn, MAX( 120.0, ssnow%ssdn(:,1)          &
                           * ssnow%osnowd / MAX( 0.01, ssnow%snowd ) + rhowat  &
                           * canopy%precis / MAX( 0.01, ssnow%snowd )  ) )

         ! permanent ice: fix hard-wired number in next version
         WHERE( soil%isoilm /= 9 )                                             &
            ssnow%ssdn(:,1) = MIN( 450.0, ssnow%ssdn(:,1) )

         canopy%precis = 0.0
         ssnow%ssdnn = ssnow%ssdn(:,1)
      
      END WHERE
   
   END WHERE ! (canopy%precis > 0. .and. ssnow%isflag == 0) 

   WHERE (canopy%precis > 0.0 .and.  ssnow%isflag > 0)
      
      ! add solid precip
      ssnow%snowd = MAX( ssnow%snowd + met%precip_sn, 0.0 )

      canopy%precis = canopy%precis - met%precip_sn  ! remaining liquid precip

      ! update top snow layer with fresh snow
      osm = ssnow%smass(:,1)
      ssnow%smass(:,1) = ssnow%smass(:,1) + met%precip_sn
      ssnow%ssdn(:,1) = MAX( 120.0,ssnow%ssdn(:,1) * osm / ssnow%smass(:,1)    &
                        + 120.0 * met%precip_sn / ssnow%smass(:,1) )

      ssnow%sdepth(:,1) = MAX( 0.02, ssnow%smass(:,1) / ssnow%ssdn(:,1) )

      ! add liquid precip
      WHERE( canopy%precis > 0.0 )
        
         ssnow%snowd = MAX( ssnow%snowd + canopy%precis, 0.0 )
         sgamm = ssnow%ssdn(:,1) * cgsnow * ssnow%sdepth(:,1)
         osm = ssnow%smass(:,1)
         
         ssnow%tggsn(:,1) = ssnow%tggsn(:,1) + canopy%precis * C%HLF           &
                            * osm / (sgamm * ssnow%osnowd )
         ssnow%smass(:,1) = ssnow%smass(:,1) + canopy%precis                   &
                            * osm/ssnow%osnowd

         ssnow%ssdn(:,1) = MAX( 120.0, MIN( ssnow%ssdn(:,1) * osm /            &
                           ssnow%smass(:,1) +  rhowat *                        &
                           ( 1.0 - osm / ssnow%smass(:,1) ), max_ssdn ) )

         ! permanent ice: fix hard-wired number in next version
         WHERE( soil%isoilm /= 9 )                                             &
            ssnow%ssdn(:,1) = MIN( 450.0, ssnow%ssdn(:,1) )

         ssnow%sdepth(:,1) = ssnow%smass(:,1)/ssnow%ssdn(:,1)

         !layer 2
         sgamm = ssnow%ssdn(:,2) * cgsnow * ssnow%sdepth(:,2)
         osm = ssnow%smass(:,2)
         ssnow%tggsn(:,2) = ssnow%tggsn(:,2) + canopy%precis * C%HLF           &
                            * osm / ( sgamm * ssnow%osnowd )
         ssnow%smass(:,2) = ssnow%smass(:,2) + canopy%precis                   &
                            * osm / ssnow%osnowd
         ssnow%ssdn(:,2) = MAX( 120.0, MIN( ssnow%ssdn(:,2) * osm /            &
                           ssnow%smass(:,2) + rhowat *                         &
                           ( 1.0 - osm / ssnow%smass(:,2) ), max_ssdn ) )

         ! permanent ice: fix hard-wired number in next version
         WHERE( soil%isoilm /= 9 )                                             &
            ssnow%ssdn(:,2) = MIN( 450.0, ssnow%ssdn(:,2) )

         ssnow%sdepth(:,2) = ssnow%smass(:,2) / ssnow%ssdn(:,2)

         !layer 3        
         sgamm = ssnow%ssdn(:,3) * cgsnow * ssnow%sdepth(:,3)
         osm = ssnow%smass(:,3)
         ssnow%tggsn(:,3) = ssnow%tggsn(:,3) + canopy%precis * C%HLF           &
                            * osm / ( sgamm * ssnow%osnowd )
         ssnow%smass(:,3) = ssnow%smass(:,3) + canopy%precis                   &
                            * osm / ssnow%osnowd
        ssnow%ssdn(:,3) = MAX( 120.0, MIN( ssnow%ssdn(:,3) * osm /             &
                          ssnow%smass(:,3) + rhowat *                          &
                          ( 1.0 - osm / ssnow%smass(:,3) ), max_ssdn ) )

         ! permanent ice: fix hard-wired number in next version
         WHERE( soil%isoilm /= 9 )                                             &
            ssnow%ssdn(:,3) = MIN(450.0,ssnow%ssdn(:,3))

         ssnow%sdepth(:,3) = ssnow%smass(:,3) / ssnow%ssdn(:,3)

         canopy%precis = 0.0
      
      END WHERE
   
   END WHERE


   ! 'fess' is for soil evap and 'fes' is for soil evap plus soil puddle evap
   canopy%segg = canopy%fess / C%HL
   canopy%segg = ( canopy%fess + canopy%fes_cor ) / C%HL

   ! Initialise snow evaporation:
   ssnow%evapsn = 0
    
   ! Snow evaporation and dew on snow
   WHERE( ssnow%snowd > 0.1 )
      
      ssnow%evapsn = dels * ( canopy%fess + canopy%fes_cor ) / ( C%HL + C%HLF )
      xxx = ssnow%evapsn

      WHERE( ssnow%isflag == 0 .AND. canopy%fess + canopy%fes_cor.GT. 0.0 )    &
         ssnow%evapsn = MIN( ssnow%snowd, xxx )
          
      WHERE( ssnow%isflag  > 0 .AND. canopy%fess + canopy%fes_cor .GT. 0.0 )   &
         ssnow%evapsn = MIN( 0.9 * ssnow%smass(:,1), xxx )

      ssnow%snowd = ssnow%snowd - ssnow%evapsn
      
      WHERE( ssnow%isflag > 0 )
         ssnow%smass(:,1) = ssnow%smass(:,1)  - ssnow%evapsn
         ssnow%sdepth(:,1) = MAX( 0.02, ssnow%smass(:,1) / ssnow%ssdn(:,1) )
      END WHERE

      canopy%segg = 0.0

   END WHERE

END SUBROUTINE snow_accum 

! -----------------------------------------------------------------------------

SUBROUTINE surfbv (dels, met, ssnow, soil, veg, canopy )

   USE cable_common_module

   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(canopy_type), INTENT(IN)       :: canopy
   
   TYPE(met_type),       INTENT(INOUT) :: met    ! all met forcing
   TYPE(soil_snow_type), INTENT(INOUT) :: ssnow  ! soil+snow variables
   
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

   CALL smoisturev( dels, ssnow, soil, veg )

   DO k = 1, ms
      xxx = REAL( soil%ssat,r_2 )
      ssnow%rnof1 = ssnow%rnof1 + REAL( MAX( ssnow%wb(:,k) - xxx, 0.0_r_2 )  &
                    * 1000.0 )  * soil%zse(k)
      ssnow%wb(:,k) = MAX( 0.01, MIN( ssnow%wb(:,k), xxx ) )
   END DO

   ! for deep runoff use wb-sfc, but this value not to exceed .99*wb-wbice
   ! account for soil/ice cracking
   ! fracm = MIN(0.2, 1. - MIN(1., ssnow%wb(:,ms) / soil%sfc ) )
   ! ssnow%wb(:,ms) = ssnow%wb(:,ms) &
   !                  + fracm*ssnow%rnof1/(1000.0*soil%zse(ms))
   ! ssnow%rnof1 = (1. - fracm) * ssnow%rnof1 

   ! Scaling  runoff to kg/m^2/s to match rest of the model
   ssnow%sinfil = 0.0
   ! lakes: replace hard-wired vegetation number in next version
   WHERE( veg%iveg == 16 )
      ssnow%sinfil = MIN( ssnow%rnof1, ssnow%wb_lake + MAX( 0.,canopy%segg ) )
      ssnow%rnof1 = MAX( 0.0, ssnow%rnof1 - ssnow%sinfil )
      ssnow%wb_lake = ssnow%wb_lake - ssnow%sinfil
      ssnow%rnof2 = MAX( 0.0, ssnow%rnof2 - ssnow%wb_lake )
   ENDWHERE

!jhan:replace nested wheres 

   !---  glacier formation
   rnof5= 0.

   IF (nglacier == 2) THEN
      
      smelt1=0.
      WHERE( ssnow%snowd > max_glacier_snowd )

         rnof5 = MIN( 0.1, ssnow%snowd - max_glacier_snowd )

         !---- change local tg to account for energy - clearly not best method
         WHERE( ssnow%isflag == 0 )
            smasstot = 0.0
            ssnow%tgg(:,1) = ssnow%tgg(:,1) - rnof5 * C%HLF                    &
                             / REAL( ssnow%gammzz(:,1) )
            ssnow%snowd = ssnow%snowd - rnof5
         ELSEWHERE
            smasstot = ssnow%smass(:,1) + ssnow%smass(:,2) + ssnow%smass(:,3)
         END WHERE

      END WHERE

      DO k = 1, 3
         
         WHERE( ssnow%snowd > max_glacier_snowd  .AND.  ssnow%isflag > 0 )
            sgamm = ssnow%ssdn(:,k) * cgsnow * ssnow%sdepth(:,k)
            smelt1(:,k) = MIN( rnof5 * ssnow%smass(:,k) / smasstot,            &
                          0.2 * ssnow%smass(:,k) )
            ssnow%smass(:,k) = ssnow%smass(:,k) - smelt1(:,k)
            ssnow%snowd = ssnow%snowd - smelt1(:,k)
         END WHERE
      
      END DO
   
      WHERE( ssnow%isflag > 0 ) rnof5 = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)
   
   END IF

   ssnow%rnof1 = ssnow%rnof1 / dels + rnof5/dels
   ssnow%rnof2 = ssnow%rnof2 / dels
   ssnow%runoff = ssnow%rnof1 + ssnow%rnof2 

END SUBROUTINE surfbv

! -----------------------------------------------------------------------------
  
! calculates temperatures of the soil
! tgg - new soil/snow temperature
! ga - heat flux from the atmosphere (ground heat flux)
! ccnsw - soil thermal conductivity, including water/ice
SUBROUTINE stempv(dels, canopy, ssnow, soil)
   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(canopy_type),    INTENT(INOUT) :: canopy
   TYPE(soil_snow_type), INTENT(INOUT) :: ssnow
   
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
            ! permanent ice: fix hard-wired number in next version
            ccnsw(j,k) = snow_ccnsw
         ELSE
            ew(j) = ssnow%wblf(j,k) * soil%ssat(j)
            exp_arg = ( ew(j) * LOG( 60.0 ) ) + ( ssnow%wbfice(j,k)            &
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
    
   WHERE(ssnow%isflag == 0)
      xx = MAX( 0., ssnow%snowd / ssnow%ssdnn )
      ccnsw(:,1) = ( ccnsw(:,1) - 0.2 ) * ( soil%zse(1) / ( soil%zse(1) + xx ) &
                   ) + 0.2
   END WHERE
    
   DO k = 3, ms
      
      WHERE (ssnow%isflag == 0)
         coeff(:,k) = 2.0 / ( soil%zse(k-1) / ccnsw(:,k-1) + soil%zse(k) /     &
                      ccnsw(:,k) )
      END WHERE
   END DO

   k = 1
   WHERE( ssnow%isflag == 0 )
      coeff(:,2) = 2.0 / ( ( soil%zse(1) + xx ) / ccnsw(:,1) + soil%zse(2) /   &
                   ccnsw(:,2) )
      coefa = 0.0
      coefb = REAL( coeff(:,2) )

      wblfsp = ssnow%wblf(:,k)
      
      xx = soil%css * soil%rhosoil
      
      ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%ssat ) * soil%css * soil%rhosoil   &
                          + soil%ssat * ( wblfsp * cswat * rhowat +            &
                          ssnow%wbfice(:,k) * csice * rhowat * 0.9 ), xx )     &
                          * soil%zse(k)

      ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + cgsnow * ssnow%snowd

      dtg = dels / ssnow%gammzz(:,k)
      
      at(:,k) = - dtg * coeff(:,k)
      ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
      bt(:,k) = 1.0 - at(:,k) - ct(:,k)

   END WHERE
   
   DO k = 2, ms
      
      WHERE( ssnow%isflag == 0 )
         
         wblfsp = ssnow%wblf(:,k)
         xx = soil%css * soil%rhosoil

         ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%ssat ) * soil%css * soil%rhosoil&
                             + soil%ssat * ( wblfsp * cswat * rhowat +         &
                             ssnow%wbfice(:,k) * csice * rhowat * 0.9 ), xx )  &
                             * soil%zse(k)

         dtg = dels / ssnow%gammzz(:,k)
         at(:,k) = - dtg * coeff(:,k)
         ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
         bt(:,k) = 1.0 - at(:,k) - ct(:,k)
      
      END WHERE
    
   END DO
   
   WHERE( ssnow%isflag == 0 )
      bt(:,1) = bt(:,1) - canopy%dgdtg * dels / ssnow%gammzz(:,1)
      ssnow%tgg(:,1) = ssnow%tgg(:,1) + ( canopy%ga - ssnow%tgg(:,1)           &
                       * REAL( canopy%dgdtg ) ) * dels /                       &
                       REAL( ssnow%gammzz(:,1) )
   END WHERE
   
   coeff(:,1-3) = 0.0  ! SO DOES THIS MEAN coeff(:,-2) ??

   ! 3-layer snow points done here
   WHERE( ssnow%isflag /= 0 )
      
      ssnow%sconds(:,1) = MAX( 0.2, MIN( 2.876e-6 * ssnow%ssdn(:,1)**2         &
                          + 0.074, max_sconds ) )
      ssnow%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,2)**2 &
                        & + 0.074, max_sconds) )
      ssnow%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssnow%ssdn(:,3)**2 &
                        & + 0.074, max_sconds) )
      coeff(:,-1) = 2.0 / (ssnow%sdepth(:,1) / ssnow%sconds(:,1) &
                       & + ssnow%sdepth(:,2) / ssnow%sconds(:,2) )
      coeff(:,0) = 2.0 / (ssnow%sdepth(:,2) / ssnow%sconds(:,2) &
                      & + ssnow%sdepth(:,3) / ssnow%sconds(:,3) )
      coeff(:,1) = 2.0 / (ssnow%sdepth(:,3) / ssnow%sconds(:,3) &
                      & + soil%zse(1) / ccnsw (:,1) )
   END WHERE
    
   DO k = 2, ms
      
      WHERE( ssnow%isflag /= 0 )                                               &
         coeff(:,k) = 2.0 / ( soil%zse(k-1) / ccnsw(:,k-1) + soil%zse(k) /     &
                      ccnsw(:,k) )

   END DO
   
   WHERE( ssnow%isflag /= 0 )
      coefa = REAL( coeff (:,-1) )
      coefb = REAL( coeff (:,1) )
   END WHERE
   
   DO k = 1, 3
      
      WHERE( ssnow%isflag /= 0 ) 
         sgamm = ssnow%ssdn(:,k) * cgsnow * ssnow%sdepth(:,k)
         dtg = dels / sgamm
         at(:,k-3) = - dtg * coeff(:,k-3)
         ct(:,k-3) = - dtg * coeff(:,k-2)
         bt(:,k-3) = 1.0 - at(:,k-3) - ct(:,k-3)
      END WHERE
   
   END DO
   
   DO k = 1, ms
     
      WHERE( ssnow%isflag /= 0 )
         wblfsp = ssnow%wblf(:,k)
         xx = soil%css * soil%rhosoil
         
         ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%ssat ) * soil%css *             &
                             soil%rhosoil + soil%ssat * ( wblfsp * cswat *     &
                             rhowat + ssnow%wbfice(:,k) * csice * rhowat *     &
                             0.9) , xx ) * soil%zse(k)

         dtg = dels / ssnow%gammzz(:,k)
         at(:,k) = - dtg * coeff(:,k)
         ct(:,k) = - dtg * coeff(:,k + 1) ! c3(ms)=0 & not really used
         bt(:,k) = 1.0 - at(:,k) - ct(:,k)
      
      END WHERE
   
   END DO

   WHERE( ssnow%isflag /= 0 )
      sgamm = ssnow%ssdn(:,1) * cgsnow * ssnow%sdepth(:,1)
      
      bt(:,-2) = bt(:,-2) - canopy%dgdtg * dels / sgamm
      
      ssnow%tggsn(:,1) = ssnow%tggsn(:,1) + ( canopy%ga - ssnow%tggsn(:,1 )    &
                         * REAL( canopy%dgdtg ) ) * dels / sgamm

      rhs(:,1-3) = ssnow%tggsn(:,1)
   END WHERE
 
   !     note in the following that tgg and tggsn are processed together
   tmp_mat(:,1:3) = REAL(ssnow%tggsn,r_2) 
   tmp_mat(:,4:(ms+3)) = REAL(ssnow%tgg,r_2)

   CALL trimb( at, bt, ct, tmp_mat, ms + 3 ) 
   
   ssnow%tggsn = REAL( tmp_mat(:,:3) )
   ssnow%tgg   = REAL( tmp_mat(:,4:(ms+3)) )
   canopy%sghflux = coefa * ( ssnow%tggsn(:,1) - ssnow%tggsn(:,2) )
   canopy%ghflux = coefb * ( ssnow%tgg(:,1) - ssnow%tgg(:,2) ) ! +ve downwards

END SUBROUTINE stempv

! -----------------------------------------------------------------------------

SUBROUTINE snowcheck(dels, ssnow, soil, met )
   
   USE cable_common_module
   
   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE(met_type),       INTENT(INOUT) :: met ! all met forcing
   
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
   
   INTEGER :: k,j
   
   DO j=1,mp

      IF( ssnow%snowd(j) <= 0.0 ) THEN
         
         ssnow%isflag(j) = 0
         ssnow%ssdn(j,:) = 120.0
         ssnow%ssdnn(j) = 120.0
         ssnow%tggsn(j,:) = C%TFRZ
         ssnow%sdepth(j,1) = ssnow%snowd(j) / ssnow%ssdn(j,1)
          
         ssnow%sdepth(j,2) = 0.
         ssnow%sdepth(j,3) = 0.

         ssnow%smass(j,1) = ssnow%snowd(j)
         ssnow%smass(j,2) = 0.0     ! EK to fix -ve sdepth 21Dec2007
         ssnow%smass(j,3) = 0.0     ! EK to fix -ve sdepth 21Dec2007
      
      ! in loop: IF( ssnow%snowd(j) <= 0.0 ) THEN
      ELSEIF( ssnow%snowd(j) < snmin * ssnow%ssdnn(j) ) THEN

         IF( ssnow%isflag(j) == 1 ) THEN
            ssnow%ssdn(j,1) = ssnow%ssdnn(j)
            ssnow%tgg(j,1) = ssnow%tggsn(j,1)
         ENDIF 

         ssnow%isflag(j) = 0
         ssnow%ssdnn(j) = MIN( 400.0, MAX( 120.0, ssnow%ssdn(j,1) ) ) 
     
         ssnow%tggsn(j,:) = MIN( C%TFRZ,ssnow%tgg(j,1) )

         ssnow%sdepth(j,1) = ssnow%snowd(j) / ssnow%ssdn(j,1)
         ssnow%sdepth(j,2) = 0.0     
         ssnow%sdepth(j,3) = 0.0     

         ssnow%smass(j,1) = ssnow%snowd(j)     
         ssnow%smass(j,2) = 0.0     
         ssnow%smass(j,3) = 0.0     

         ssnow%ssdn(j,:) = ssnow%ssdnn(j)
         
         IF( .NOT.cable_user%CABLE_RUNTIME_COUPLED ) THEN
            IF( soil%isoilm(j) == 9 .AND. ktau_gl <= 2 )                       &
               ! permanent ice: fixed hard-wired number in next version
               ssnow%ssdnn(j) = 700.0
         ENDIF
      
      
      ELSE ! in loop: IF( ssnow%snowd(j) <= 0.0 ) THEN
           ! sufficient snow now for 3 layer snowpack
         
         IF( ssnow%isflag(j) == 0 ) THEN

            ssnow%tggsn(j,:) = MIN( C%TFRZ, ssnow%tgg(j,1) )

            ssnow%ssdn(j,2) = ssnow%ssdn(j,1)
            ssnow%ssdn(j,3) = ssnow%ssdn(j,1)

            IF( .NOT. cable_user%cable_runtime_coupled) THEN
               IF( soil%isoilm(j) == 9 .AND. ktau_gl <= 2 ) THEN
                  ! permanent ice: fix hard-wired number in next version
                  ssnow%ssdn(j,1)  = 450.0
                  ssnow%ssdn(j,2)  = 580.0
                  ssnow%ssdn(j,3)  = 600.0
               ENDIF
            ENDIF
               
            ssnow%sdepth(j,1) = ssnow%t_snwlr(j)
            
            ssnow%smass(j,1)  =  ssnow%t_snwlr(j) * ssnow%ssdn(j,1)
            
            ssnow%smass(j,2)  = ( ssnow%snowd(j) - ssnow%smass(j,1) ) * 0.4
            ssnow%smass(j,3)  = ( ssnow%snowd(j) - ssnow%smass(j,1) ) * 0.6
            
            ssnow%sdepth(j,2) = ssnow%smass(j,2) / ssnow%ssdn(j,2)
            ssnow%sdepth(j,3) = ssnow%smass(j,3) / ssnow%ssdn(j,3)
            
            ssnow%ssdnn(j) = ( ssnow%ssdn(j,1) * ssnow%smass(j,1) +            &
                              ssnow%ssdn(j,2) * ssnow%smass(j,2) +             &
                              ssnow%ssdn(j,3) * ssnow%smass(j,3) )             &
                              / ssnow%snowd(j)
         
         ENDIF 
         
         ssnow%isflag(j) = 1
      
      ENDIF ! END: IF( ssnow%snowd(j) <= 0.0 ) THEN

               
   ENDDO ! END: DO j=1,mp

END SUBROUTINE snowcheck 

! -----------------------------------------------------------------------------

SUBROUTINE snowl_adjust(dels, ssnow, canopy )
   
   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(soil_snow_type), INTENT(INOUT) :: ssnow
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
   WHERE( ssnow%isflag > 0 )

      WHERE( ssnow%sdepth(:,1) > ssnow%t_snwlr )

         excd = ssnow%sdepth(:,1) - ssnow%t_snwlr
         excm = excd * ssnow%ssdn(:,1)
         ssnow%sdepth(:,1) = ssnow%sdepth(:,1) - REAL(excd)
         osm = ssnow%smass(:,1)
         ssnow%smass(:,1) = ssnow%smass(:,1) - REAL(excm)
            
         osm = ssnow%smass(:,2)
         ssnow%smass(:,2) = MAX( 0.01, ssnow%smass(:,2) + REAL(excm) )
         
         ssnow%ssdn(:,2) = REAL( MAX( 120.0_r_2, MIN( REAL( max_ssdn, r_2 ),   &
                           ssnow%ssdn(:,2) * osm / ssnow%smass(:,2) +          &
                           ssnow%ssdn(:,1) * excm / ssnow%smass(:,2) ) ) )

          ssnow%sdepth(:,2) =  ssnow%smass(:,2) / ssnow%ssdn(:,2) 
         
          ssnow%tggsn(:,2) = REAL( ssnow%tggsn(:,2) * osm / ssnow%smass(:,2)   &
                             + ssnow%tggsn(:,1) * excm / ssnow%smass(:,2) )

          ! following line changed to fix -ve sdepth (EK 21Dec2007)
          ssnow%smass(:,3) = MAX( 0.01, ssnow%snowd - ssnow%smass(:,1)         &
                             - ssnow%smass(:,2) )

      ELSEWHERE ! ssnow%sdepth(:,1) < ssnow%t_snwlr
      
         ! 1st layer
         excd = ssnow%t_snwlr - ssnow%sdepth(:,1)
         excm = excd * ssnow%ssdn(:,2)
         osm = ssnow%smass(:,1)
         ssnow%smass(:,1) = ssnow%smass(:,1) + REAL(excm)
         ssnow%sdepth(:,1) = ssnow%t_snwlr
         ssnow%ssdn(:,1) = REAL( MAX( 120.0_r_2, MIN( REAL( max_ssdn,r_2 ),    &
                           ssnow%ssdn(:,1) * osm / ssnow%smass(:,1)            &
                           + ssnow%ssdn(:,2) * excm / ssnow%smass(:,1) ) ) )

         ssnow%tggsn(:,1) = REAL( ssnow%tggsn(:,1) * osm / ssnow%smass(:,1)   &
                          + ssnow%tggsn(:,2) * excm / ssnow%smass(:,1) )
          
         ! 2nd layer
         ssnow%smass(:,2) = MAX( 0.01, ssnow%smass(:,2) - REAL(excm) )
         ssnow%sdepth(:,2) = ssnow%smass(:,2) / ssnow%ssdn(:,2)

         ! following line changed to fix -ve sdepth (EK 21Dec2007)
         ssnow%smass(:,3) = MAX( 0.01, ssnow%snowd - ssnow%smass(:,1)          &
                            - ssnow%smass(:,2) )

      END WHERE
   
   END WHERE 

   DO  api=1,mp
      
      IF( ssnow%isflag(api).GT.0 ) THEN
      
         frac(api) = ssnow%smass(api,2) / MAX( 0.02, ssnow%smass(api,3) )
         ! if frac > 0.6 or frac < 0.74 do nothing 
         ! HOW TO translate this to xfrac
         xfrac(api) = 2.0/3.0/ frac(api)
         
         IF( xfrac(api) > 1.0 ) THEN

            excm(api) = (xfrac(api) - 1.0) * ssnow%smass(api,2)
            osm(api) = ssnow%smass(api,2)
            
            ! changed 0.02 to 0.01 to fix -ve sdepth (EK 21Dec2007)
            ssnow%smass(api,2) = MAX( 0.01, ssnow%smass(api,2) +               &
                                 REAL( excm(api) ) )

            ssnow%tggsn(api,2) = ssnow%tggsn(api,2) * osm(api) /               &
                                 ssnow%smass(api,2) +  ssnow%tggsn(api,3)      &
                                 * REAL( excm(api) )/ ssnow%smass(api,2)

            ssnow%ssdn(api,2) = MAX( 120.0, MIN( max_ssdn, ssnow%ssdn(api,2) * &
                                osm(api) / ssnow%smass(api,2) +                &
                                ssnow%ssdn(api,3) * REAL( excm(api) )          &
                                / ssnow%smass(api,2) ) )

            ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
            ssnow%smass(api,3) = MAX( 0.01, ssnow%snowd(api) -                 &
                                 ssnow%smass(api,1) - ssnow%smass(api,2) )
            
            ssnow%sdepth(api,3) = MAX( 0.02, ssnow%smass(api,3) /              &
                                  ssnow%ssdn(api,3) )
            
         ELSE! xfrac < 1
           
            excm(api) = ( 1 - xfrac(api) ) * ssnow%smass(api,2)
            ssnow%smass(api,2) = MAX(0.01, ssnow%smass(api,2) - REAL(excm(api)))
           ssnow%sdepth(api,2) = MAX(0.02, ssnow%smass(api,2) /                &
                                 ssnow%ssdn(api,2) )

           osm(api) = ssnow%smass(api,3)
           ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
           ssnow%smass(api,3) = MAX(0.01, &
                                ssnow%snowd(api) - ssnow%smass(api,1) -        &
                                ssnow%smass(api,2) )


           ssnow%tggsn(api,3) = ssnow%tggsn(api,3) * osm(api) /                &
                                ssnow%smass(api,3) +  ssnow%tggsn(api,2) *     &
                                REAL( excm(api) ) / ssnow%smass(api,3)
           ssnow%ssdn(api,3) = MAX(120.0, MIN( max_ssdn, ssnow%ssdn(api, 3 )*  &
                               osm(api) / ssnow%smass(api,3) +                 &
                               ssnow%ssdn(api,2) * REAL( excm(api) )           &
                               / ssnow%smass(api,3) ) )
           ssnow%sdepth(api,3) = ssnow%smass(api,3) /  ssnow%ssdn(api,3)

        END IF
        
        ssnow%isflag(api) = 1
        
        ssnow%ssdnn(api) = ( ssnow%ssdn(api,1) * ssnow%sdepth(api,1) +         &
                           ssnow%ssdn(api,2) * ssnow%sdepth(api,2) +           &
                           ssnow%ssdn(api,3) * ssnow%sdepth(api,3) )           &
                           / ( ssnow%sdepth(api,1) + ssnow%sdepth(api,2)       &
                           + ssnow%sdepth(api,3) )

      END IF

   END DO

END SUBROUTINE snowl_adjust

! -----------------------------------------------------------------------------

SUBROUTINE soilfreeze(dels, soil, ssnow)
   USE cable_common_module
   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil
   REAL(r_2), DIMENSION(mp)           :: sicefreeze
   REAL(r_2), DIMENSION(mp)           :: sicemelt
   REAL, DIMENSION(mp)           :: xx
   INTEGER k

   xx = 0.
   DO k = 1, ms
      
      WHERE (ssnow%tgg(:,k) < C%TFRZ &
          & .AND. frozen_limit * ssnow%wb(:,k) - ssnow%wbice(:,k) > .001)
         
         sicefreeze = MIN( MAX( 0.0_r_2, ( frozen_limit * ssnow%wb(:,k) -      &
                      ssnow%wbice(:,k) ) ) * soil%zse(k) * 1000.0,             &
                      ( C%TFRZ - ssnow%tgg(:,k) ) * ssnow%gammzz(:,k) / C%HLF )
         ssnow%wbice(:,k) = MIN( ssnow%wbice(:,k) + sicefreeze / (soil%zse(k)  &
                            * 1000.0), frozen_limit * ssnow%wb(:,k) )
         xx = soil%css * soil%rhosoil
         ssnow%gammzz(:,k) = MAX(                                              &
             REAL((1.0 - soil%ssat) * soil%css * soil%rhosoil ,r_2)            &
             + (ssnow%wb(:,k) - ssnow%wbice(:,k)) * REAL(cswat * rhowat,r_2)   &
             + ssnow%wbice(:,k) * REAL(csice * rhowat * 0.9,r_2),              &
             REAL(xx,r_2)) * REAL( soil%zse(k),r_2 )

         WHERE (k == 1 .AND. ssnow%isflag == 0)
            ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + cgsnow * ssnow%snowd
         END WHERE
         ssnow%tgg(:,k) = ssnow%tgg(:,k) + REAL(sicefreeze)                    &
                          * C%HLF / REAL(ssnow%gammzz(:,k) )
      
      ELSEWHERE( ssnow%tgg(:,k) > C%TFRZ .AND. ssnow%wbice(:,k) > 0. )
         
         sicemelt = MIN( ssnow%wbice(:,k) * soil%zse(k) * 1000.0,              &
                    ( ssnow%tgg(:,k) - C%TFRZ ) * ssnow%gammzz(:,k) / C%HLF )
         
         ssnow%wbice(:,k) = MAX( 0.0_r_2, ssnow%wbice(:,k) - sicemelt          &
                            / (soil%zse(k) * 1000.0) )
         xx = soil%css * soil%rhosoil
         ssnow%gammzz(:,k) = MAX(                                              &
              REAL((1.0-soil%ssat) * soil%css * soil%rhosoil,r_2)             &
              + (ssnow%wb(:,k) - ssnow%wbice(:,k)) * REAL(cswat*rhowat,r_2)   &
              + ssnow%wbice(:,k) * REAL(csice * rhowat * 0.9,r_2),            &
              REAL(xx,r_2) ) * REAL(soil%zse(k),r_2)
         WHERE (k == 1 .AND. ssnow%isflag == 0)
            ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + cgsnow * ssnow%snowd
         END WHERE
         ssnow%tgg(:,k) = ssnow%tgg(:,k) - REAL(sicemelt)                     &
                          * C%HLF / REAL(ssnow%gammzz(:,k))
       
      END WHERE
    
   END DO

END SUBROUTINE soilfreeze

! -----------------------------------------------------------------------------

SUBROUTINE remove_trans(dels, soil, ssnow, canopy, veg)
   
   USE cable_common_module, ONLY : redistrb

   ! Removes transpiration water from soil.
   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(canopy_type), INTENT(INOUT)         :: canopy
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
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
         xx = canopy%fevc * dels / C%HL * veg%froot(:,k) + diff(:,k-1)   ! kg/m2
         diff(:,k) = MAX( 0.0, ssnow%wb(:,k) - soil%swilt) &      ! m3/m3
                     * soil%zse(k)*1000.0
         xxd = xx - diff(:,k)
       
         WHERE ( xxd .GT. 0.0 )
            ssnow%wb(:,k) = ssnow%wb(:,k) - diff(:,k) / (soil%zse(k)*1000.0)
            diff(:,k) = xxd
         ELSEWHERE
            ssnow%wb(:,k) = ssnow%wb(:,k) - xx / (soil%zse(k)*1000.0)
            diff(:,k) = 0.0
         ENDWHERE
     
     END WHERE
   
   END DO

END SUBROUTINE remove_trans 

! -----------------------------------------------------------------------------

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
!	 ssnow
SUBROUTINE soil_snow(dels, soil, ssnow, canopy, met, bal, veg)
   USE cable_common_module
   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(canopy_type), INTENT(INOUT)         :: canopy
   TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
   TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
   TYPE (balances_type), INTENT(INOUT)      :: bal
   INTEGER             :: k
   REAL, DIMENSION(mp) :: snowmlt
   REAL, DIMENSION(mp) :: totwet
   REAL, DIMENSION(mp) :: weting
   REAL, DIMENSION(mp) :: xxx, tgg_old, tggsn_old
   REAL(r_2), DIMENSION(mp) :: xx,deltat,sinfil1,sinfil2,sinfil3 
   REAL                :: zsetot
   INTEGER, SAVE :: ktau =0 
   
   CALL point2constants( C ) 
   cp = C%CAPP
    
   ktau = ktau +1 
   
   !jhan - make switchable 
   ! appropriate for ACCESS1.0
   !max_glacier_snowd = 50000.0
   ! appropriate for ACCESS1.3 
   max_glacier_snowd = 1100.0

   zsetot = sum(soil%zse) 
   ssnow%tggav = 0.
   DO k = 1, ms
      ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot
   END DO


   IF( cable_runtime%offline .or. cable_runtime%mk3l ) THEN
        ssnow%t_snwlr = 0.05
   ENDIF

   ssnow%fwtop1 = 0.0
   ssnow%fwtop2 = 0.0
   ssnow%fwtop3 = 0.0
   ssnow%runoff = 0.0 ! initialise total runoff
   ssnow%rnof1 = 0.0 ! initialise surface runoff
   ssnow%rnof2 = 0.0 ! initialise deep drainage
   ssnow%smelt = 0.0 ! initialise snowmelt
   ssnow%dtmlt = 0.0 
   ssnow%osnowd = ssnow%snowd


   IF( .NOT.cable_user%cable_runtime_coupled ) THEN
   
      IF( ktau_gl <= 1 ) THEN
         
         IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
                                                  ! after discussion with BP
         ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
         ssnow%wbtot = 0.0
         DO k = 1, ms
            ssnow%wb(:,k)  = MIN( soil%ssat,MAX ( ssnow%wb(:,k), soil%swilt ))
         END DO
   
         ssnow%wb(:,ms-2)  = MIN( soil%ssat, MAX ( ssnow%wb(:,ms-2),           &
                             0.5 * ( soil%sfc + soil%swilt ) ) )
         ssnow%wb(:,ms-1)  = MIN( soil%ssat, MAX ( ssnow%wb(:,ms-1),           &
                             0.8 * soil%sfc ) )
         ssnow%wb(:,ms)    = MIN( soil%ssat, MAX ( ssnow%wb(:,ms), soil%sfc) )
         
         DO k = 1, ms
            
            WHERE( ssnow%tgg(:,k) <= C%TFRZ .AND. ssnow%wbice(:,k) <= 0.01 )   &
               ssnow%wbice(:,k) = 0.5 * ssnow%wb(:,k)
            
            WHERE( ssnow%tgg(:,k) < C%TFRZ)                                    &
               ssnow%wbice(:,k) = frozen_limit * ssnow%wb(:,k)
            
         END DO
   
         WHERE (soil%isoilm == 9) 
            ! permanent ice: fix hard-wired number in next version
            ssnow%snowd = max_glacier_snowd
            ssnow%osnowd = max_glacier_snowd
            ssnow%tgg(:,1) = ssnow%tgg(:,1) - 1.0
            ssnow%wb(:,1) = 0.95 * soil%ssat
            ssnow%wb(:,2) = 0.95 * soil%ssat
            ssnow%wb(:,3) = 0.95 * soil%ssat
            ssnow%wb(:,4) = 0.95 * soil%ssat
            ssnow%wb(:,5) = 0.95 * soil%ssat
            ssnow%wb(:,6) = 0.95 * soil%ssat
            ssnow%wbice(:,1) = 0.90 * ssnow%wb(:,1)
            ssnow%wbice(:,2) = 0.90 * ssnow%wb(:,2)
            ssnow%wbice(:,3) = 0.90 * ssnow%wb(:,3)
            ssnow%wbice(:,4) = 0.90 * ssnow%wb(:,4)
            ssnow%wbice(:,5) = 0.90 * ssnow%wb(:,5)
            ssnow%wbice(:,6) = 0.90 * ssnow%wb(:,6)
         ENDWHERE
         
         xx=soil%css * soil%rhosoil
         
         ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
              & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * cswat * rhowat &
              & + ssnow%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1)
      END IF
   ENDIF  ! if(.NOT.cable_runtime_coupled)

   xx=soil%css * soil%rhosoil
   IF (ktau <= 1)                                                              &
     ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil      &
            & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * cswat * rhowat           &
            & + ssnow%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1) +   &
            & (1. - ssnow%isflag) * cgsnow * ssnow%snowd



   DO k = 1, ms ! for stempv
      
      ! Set liquid soil water fraction (fraction of saturation value):
      ssnow%wblf(:,k) = MAX( 0.01_r_2, (ssnow%wb(:,k) - ssnow%wbice(:,k)) )    &
           & / REAL(soil%ssat,r_2)
      
      ! Set ice soil water fraction (fraction of saturation value):
      ssnow%wbfice(:,k) = REAL(ssnow%wbice(:,k)) / soil%ssat
   END DO
 
   CALL snowcheck (dels, ssnow, soil, met )

   CALL snowdensity (dels, ssnow, soil)

   CALL snow_accum (dels, canopy, met, ssnow, soil )

   CALL snow_melting (dels, snowmlt, ssnow, soil )

   ! Add snow melt to global snow melt variable:
   ssnow%smelt = snowmlt

   ! Adjust levels in the snowpack due to snow accumulation/melting,
   ! snow aging etc...
   CALL snowl_adjust(dels, ssnow, canopy )

   CALL stempv(dels, canopy, ssnow, soil)
   
   ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

   CALL snow_melting (dels, snowmlt, ssnow, soil )
   
   ! Add new snow melt to global snow melt variable: 
   ssnow%smelt = ssnow%smelt + snowmlt

   CALL remove_trans(dels, soil, ssnow, canopy, veg)

   CALL  soilfreeze(dels, soil, ssnow)


   totwet = canopy%precis + ssnow%smelt 
   
   ! total available liquid including puddle
   weting = totwet + max(0.,ssnow%pudsto - canopy%fesp/C%HL*dels) 
   xxx=soil%ssat - ssnow%wb(:,1)
  
   sinfil1 = MIN( 0.95*xxx*soil%zse(1)*rhowat, weting) !soil capacity
   xxx=soil%ssat - ssnow%wb(:,2)
   sinfil2 = MIN( 0.95*xxx*soil%zse(2)*rhowat, weting - sinfil1) !soil capacity
   xxx=soil%ssat - ssnow%wb(:,3)
   sinfil3 = MIN( 0.95*xxx*soil%zse(3)*rhowat,weting-sinfil1-sinfil2)
   
   ! net water flux to the soil
   ssnow%fwtop1 = sinfil1 / dels - canopy%segg          
   ssnow%fwtop2 = sinfil2 / dels           
   ssnow%fwtop3 = sinfil3 / dels           

   ! Puddle for the next time step
   ssnow%pudsto = max( 0., weting - sinfil1 - sinfil2 - sinfil3 )
   ssnow%rnof1 = max(0.,ssnow%pudsto - ssnow%pudsmx)
   ssnow%pudsto = ssnow%pudsto - ssnow%rnof1

   CALL surfbv(dels, met, ssnow, soil, veg, canopy )

   ! correction required for energy balance in online simulations 
   IF( cable_runtime%um) THEN
      canopy%fhs_cor = ssnow%dtmlt(:,1)*ssnow%dfh_dtg
      canopy%fes_cor = ssnow%dtmlt(:,1)*(ssnow%cls*ssnow%dfe_ddq * ssnow%ddq_dtg)

      canopy%fhs = canopy%fhs+canopy%fhs_cor
      canopy%fes = canopy%fes+canopy%fes_cor
   ENDIF

   ! redistrb (set in cable.nml) by default==.FALSE. 
   IF( redistrb )                                                              &
      CALL hydraulic_redistribution( dels, soil, ssnow, canopy, veg, met )

   ssnow%smelt = ssnow%smelt/dels

   ! Set weighted soil/snow surface temperature
   ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

   ssnow%wbtot = 0.0
   DO k = 1, ms
      ssnow%wbtot = ssnow%wbtot + REAL(ssnow%wb(:,k)*1000.0*soil%zse(k),r_2)
   END DO

END SUBROUTINE soil_snow

! -----------------------------------------------------------------------------

!+++++++++++++++++++  Hydraulic Redistribution Section  ++++++++++++++++++++++
! Science from Ryel et al. Oecologia, 2002; Lee et al., 2005, PNAS
! Code by LiLH 16 Feb, 2011
! Fixed problem of negative wb in global run by BP Mar 2011
SUBROUTINE hydraulic_redistribution(dels, soil, ssnow, canopy, veg, met)

   USE cable_common_module, ONLY : wiltParam, satuParam
   
   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(soil_parameter_type), INTENT(IN) :: soil
   TYPE(canopy_type),         INTENT(IN) :: canopy
   TYPE(veg_parameter_type),  INTENT(IN) :: veg
   
   TYPE(soil_snow_type),   INTENT(INOUT) :: ssnow
   TYPE(met_type),         INTENT(INOUT) :: met 
   
   REAL, PARAMETER ::                                                         &
      thetas=0.45,         & ! from Belk et al., 2007, WRR
      thetar=0.20 ,        & ! from Belk et al., 2007, WRR
      n_hr = 3.22,         & ! --
      wpsy50 = -1.0,       & ! MPa
      n_VG = 2.06,         & ! -- 2.06
      m_VG = 1.0-1.0/n_VG, & ! --
      alpha_VG = 0.00423,  & ! cm^{-1} Note: 1cmH2O=100Pa
      CRT = 125.0            ! cm MPa^-1 h^-1, default value (0.097) 
                             ! from Ryel et al., 2002
   REAL, DIMENSION(mp) ::                                                      &
      frootX,      & ! --
      Dtran,       & ! Swith for hr
      available,   &
      accommodate, &
      totalmoist,  &
      totalice,    &
      total2,      &
      zsetot,      &
      temp

   REAL, DIMENSION(mp,ms)::                                                    &
      S_VG, & ! --
      wpsy, & ! MPa
      C_hr    ! --
  
   REAL, DIMENSION(mp,ms,ms) ::                                                &
      hr_term,    & ! cm/hour
      hr_perTime    !

   INTEGER :: j, k

   zsetot = sum(soil%zse)
   totalmoist(:) = 0.0
   totalice(:) = 0.0
   DO k=1, ms
     totalmoist(:) = totalmoist(:) + ssnow%wb(:,k)*soil%zse(k)/zsetot
     totalice(:) = totalice(:) + ssnow%wbice(:,k)*soil%zse(k)/zsetot
   ENDDO

   Dtran=0.0
   WHERE( canopy%fevc < 10.0 .and.  totalice  < 1.e-2 )  Dtran=1.0
   
   DO k=1, ms
      S_VG(:,k) = MIN( 1.0, MAX( 1.0E-4, ssnow%wb(:,k) - soil%swilt )          &
                            / ( soil%ssat - soil%swilt ) )
      ! VG model, convert from cm to Pa by (*100), to MPa (*1.0E-6)
      wpsy(:,k) = -1.0 / alpha_VG * ( S_VG(:,k)**(-1.0/m_VG) - 1.0 )**(1/n_VG) &
                  * 100 * 1.0E-6  
      
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
         
         ! Overwrite to give zero redistribution for all types except
         ! evergreen broadleaf (2) and c4 grass (7)
         ! NB: Hard-wired numbers should be removed in future version
         WHERE( .NOT.(veg%iveg == 2 .OR. veg%iveg == 7 ) )
            hr_perTime(:,k,j) = 0.0
            hr_perTime(:,j,k) = 0.0
         ENDWHERE
         
         WHERE( hr_perTime(:,k,j) < 0.0 )

            available(:)   = MAX( 0.0, ssnow%wb(:,k) -                         &
                            ( soil%swilt(:) + ( soil%sfc(:) - soil%swilt(:) )  &
                             / 3. ) )
            accommodate(:) = MAX( 0.0, soil%ssat(:) - ssnow%wb(:,j) )
            
            temp(:) = MAX( hr_perTime(:,k,j),                                  &
                          -1.0 * wiltParam * available(:),                     &
                          -1.0 * satuParam * accommodate(:) * soil%zse(j) /    &
                          soil%zse(k) ) 
            
            hr_perTime(:,k,j) = temp(:)
            hr_perTime(:,j,k) = -1.0 * temp(:) * soil%zse(k) / soil%zse(j)
         
         ELSEWHERE (hr_perTime(:,j,k) < 0.0)

           available(:)   = MAX( 0.0, ssnow%wb(:,j) -                          &
                            ( soil%swilt(:) + ( soil%sfc(:) - soil%swilt(:) )  &
                            / 3. ) ) 
           
           accommodate(:) = MAX( 0.0, soil%ssat(:) - ssnow%wb(:,k) )
           
           temp(:) = MAX( hr_perTime(:,j,k),                                   &
                         - 1.0 * wiltParam * available(:),                     &
                         -1.0 * satuParam * accommodate(:) * soil%zse(k) /     &
                         soil%zse(j) )

           hr_perTime(:,j,k) = temp(:)
           hr_perTime(:,k,j) = -1.0 * temp(:) * soil%zse(j) / soil%zse(k)
         
         ENDWHERE
         
         ssnow%wb(:,k) = ssnow%wb(:,k) + hr_perTime(:,k,j)
         ssnow%wb(:,j) = ssnow%wb(:,j) + hr_perTime(:,j,k)
      
      ENDDO 
   
   ENDDO

   WHERE( met%tk < C%TFRZ + 5.  ) Dtran=0.0
     
   DO k=1, ms
      S_VG(:,k) = MIN( 1.0, MAX( 1.0E-4, ssnow%wb(:,k) - soil%swilt )          &
                  / ( soil%ssat - soil%swilt ) )
      
      ! VG model, convert from cm to Pa by (*100), to MPa (*1.0E-6)
      wpsy(:,k) = -1.0 / alpha_VG * ( S_VG(:,k)**(-1.0/m_VG) -1.0 )**(1/n_VG)  &
                  * 100 * 1.0E-6     
      
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
         
         ! Overwrite to give zero redistribution for all types except
         ! evergreen broadleaf (2) and c4 grass (7)
         ! NB: Hard-wired numbers should be removed in future version
         WHERE( .NOT.( veg%iveg == 2 .OR. veg%iveg == 7 ) )
            hr_perTime(:,k,j) = 0.0
            hr_perTime(:,j,k) = 0.0
         ENDWHERE
         
         WHERE( hr_perTime(:,k,j) < 0.0 )
            
            available(:)   = MAX( 0.0, ssnow%wb(:,k) - soil%sfc(:) )
            accommodate(:) = MAX( 0.0, soil%ssat(:) - ssnow%wb(:,j) )
            
            temp(:) = MAX(hr_perTime(:,k,j),                                   &
                          -1.0 * wiltParam*available(:),                       &
                          -1.0 * satuParam * accommodate(:) * soil%zse(j) /    &
                          soil%zse(k) )

            hr_perTime(:,k,j) = temp(:)
            hr_perTime(:,j,k) = -1.0 * temp(:) * soil%zse(k) / soil%zse(j)
         
         ELSEWHERE (hr_perTime(:,j,k) < 0.0)
            
            available(:)   = MAX( 0.0, ssnow%wb(:,j)- soil%sfc(:) )
            accommodate(:) = MAX( 0.0, soil%ssat(:)-ssnow%wb(:,k) )
            
            temp(:) = MAX(hr_perTime(:,j,k),                                   &
                     -1.0 * wiltParam*available(:),                            &
                     -1.0 * satuParam * accommodate(:) * soil%zse(k) /         &
                     soil%zse(j) )
            
            hr_perTime(:,j,k) = temp(:)
            hr_perTime(:,k,j) = -1.0 * temp(:) * soil%zse(j) / soil%zse(k)
         
         ENDWHERE
         
         ssnow%wb(:,k) = ssnow%wb(:,k) + hr_perTime(:,k,j)
         ssnow%wb(:,j) = ssnow%wb(:,j) + hr_perTime(:,j,k)
      ENDDO
   ENDDO
                          
END SUBROUTINE hydraulic_redistribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!MD GW code from here on!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------
! SUBROUTINE solve_tridiag
! solves tridiagonal set of linear equations.  returns the solution
!  
  subroutine solve_tridiag (at, bt, ct, rt, ut,n)

! !ARGUMENTS:
    implicit none
    real(r_2), intent(in)    :: at(:,:)    ! 1 left of diagonal 
    real(r_2), intent(in)    :: bt(:,:)    ! diagonal 
    real(r_2), intent(in)    :: ct(:,:)    ! 1 right of diagonal 
    real(r_2), intent(in)    :: rt(:,:)    ! right hand side
    real(r_2), intent(inout) :: ut(:,:)    ! solution to the system of eqs
    integer, intent(in) :: n
!   local variables
    integer  :: k
    REAL(r_2), DIMENSION(mp,ms+1) ::  gam  
    REAL(r_2), DIMENSION(mp)      ::  bet 
!-----------------------------------------------------------------------

    ! Solve the matrix
    bet(:) = bt(:,1)
       
    where (bet(:) .ne. 0.0_r_2) ut(:,1) = rt(:,1) / bet(:)
    where (bet(:) .eq. 0.0_r_2) ut(:,1) = rt(:,1) / (bet(:) + 0.000001_r_2)
    do k = 1,n
       gam(:,k) = ct(:,k-1) / bet(:)
       bet(:)   = max(bt(:,k) - at(:,k) * gam(:,k),0.00001_r_2)
       ut(:,k)   = (rt(:,k) - at(:,k)*ut(:,k-1)) / bet(:)
    end do

    do k = n-1,2,-1
       ut(:,k) = ut(:,k) - gam(:,k+1) * ut(:,k+1)
    end do

  end subroutine solve_tridiag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  !-------------------------------------------------------------------------
  SUBROUTINE ovrlndflx (dels, ktau, ssoil, soil,prin )
  IMPLICIT NONE
    REAL, INTENT(IN)                         :: dels ! integration time step (s)
    INTEGER, INTENT(IN)                      :: ktau ! integration step number
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    LOGICAL, INTENT(IN)                      :: prin
    INTEGER, PARAMETER                       :: ntest = 0 ! for snow diag prints
    INTEGER, PARAMETER                       :: nglacier = 2 ! 0 original, 1 off, 2 new Eva
    INTEGER                                  :: k
    REAL, DIMENSION(mp)                :: rnof5
    REAL, DIMENSION(mp)                :: sgamm
    REAL, DIMENSION(mp)                :: smasstot
    REAL, DIMENSION(mp,0:3)            :: smelt1
    REAL(r_2), DIMENSION(mp)           :: xxx,fice,icef,efpor
    REAL(r_2), DIMENSION(mp)           :: tmpa,tmpb,inflmx,fcov
    REAL(r_2), dimension(mp)           :: satfrac
    logical                                  :: prinall = .false.  !for debugging
    
    !For now assume there is no puddle?
    ssoil%pudsto = 0._r_2
    
    
    efpor(:) = max(soil%watsat(:,1) - ssoil%wbice(:,1),0.05_r_2)
    !srf frozen fraction.  should be based on topography
    icef(:) = ssoil%wb(:,1)*denliq/&
              (ssoil%wbice(:,1)*denice+ssoil%wb(:,1)*denliq)
    fice(:) = min(max(exp(-3.0*(1.0-ssoil%icefrac(:,1)))- exp(-3.0)&
              ,0.0_r_2),0.99_r_2)
    ! Saturated fraction
    satfrac(:) = (1.0-fice(:))*0.3*exp(-0.2_r_2*ssoil%wtd(:)/1000.0_r_2)+fice(:)

    ! Maximum infiltration capacity
    tmpa(:)  = max(ssoil%wb(:,1)/max(efpor(:),0.075_r_2),0.01_r_2)
    tmpb(:)  = max((tmpa(:)-satfrac(:)) / (1.0_r_2-satfrac(:)),0.0_r_2)
    inflmx(:)= ((1.0-soil%clappB(:,1)*soil%smpsat(:,1)/(0.5*soil%zse(1)*1000.0)&
                    *(tmpb-1.0))*soil%hksat(:,1))
     ! Surface runoff
    where (ssoil%fwtop(:) .gt. inflmx*dels)
       ssoil%rnof1(:) =  (satfrac(:) * ssoil%fwtop(:)/dels + &
                       (1.0-satfrac(:)) * ssoil%fwtop(:)/dels-inflmx)  !in mm/s
    elsewhere
       ssoil%fwtop(:) = ssoil%fwtop(:).dels - ssoil%rnof1(:)  !in mm/s
    end where
           
    !in soil_snow_gw subroutine of this module
    !ssoil%runoff = ssoil%rnof1! + ssoil%rnof2    

    !---  glacier formation
    IF (nglacier == 2) THEN
       WHERE (ssoil%snowd > 1000.0)
          rnof5 = ssoil%snowd - 1000.0
          ssoil%runoff = ssoil%runoff + rnof5
          !---- change local tg to account for energy - clearly not best method
          WHERE (ssoil%isflag == 0)
             smasstot = 0.0
             ssoil%tgg(:,1) = ssoil%tgg(:,1) - rnof5 * C%HLF &
                  & / ssoil%gammzz(:,1)
             ssoil%snowd = 1000.0
          ELSEWHERE
             smasstot = ssoil%smass(:,1) + ssoil%smass(:,2) + ssoil%smass(:,3)
          END WHERE
       END WHERE
       
       DO k = 1, 3
          WHERE (ssoil%snowd > 1000.0 .AND. ssoil%isflag > 0)
             sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
             smelt1(:,k) = MIN(rnof5 * ssoil%smass(:,k) / smasstot, &
                  & 0.9 * ssoil%smass(:,k) )
             ssoil%smass(:,k) = ssoil%smass(:,k) - smelt1(:,k)
             ssoil%snowd = ssoil%snowd - smelt1(:,k)
             ssoil%tggsn(:,k) = ssoil%tggsn(:,k) - smelt1(:,k) * C%HLF / sgamm
          END WHERE
       END DO
    END IF

  END SUBROUTINE ovrlndflx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  !----------------------------------------------------------------------
  ! SUBROUTINE calcwtd
  !
  ! Iteratively calcs the water table depth by equating the mass of water in the
  ! soil column to the mass of a hydrostatic column inegrated from the surface to the 
  ! water table depth
  !  
  SUBROUTINE calcwtd (ssoil, soil,prin)
  IMPLICIT NONE
   TYPE (soil_snow_type), INTENT(INOUT)      :: ssoil ! soil and snow variables
   TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
   LOGICAL, INTENT(IN)                       :: prin  !print info?
  !Local vars 
  REAL(r_2), DIMENSION(mp,ms)          :: dzmm
  REAL(r_2), DIMENSION(0:ms+1)                 :: zimm
  REAL(r_2), DIMENSION(ms)                   :: zmm
  REAL(r_2), DIMENSION(mp)             :: GWzimm,temp
  REAL(r_2), DIMENSION(mp)             :: def,defc     

  REAL(r_2)                                  :: deffunc,tempa,tempb,derv,calc,tmpc
  REAL(r_2), DIMENSION(mp)             :: invB,Nsmpsat  !inverse of C&H B,Nsmpsat
  !local loop counters
  INTEGER :: k,i,wttd,jlp
  LOGICAL :: keeplooping     !used to exit iterations on wtd
     
   
  !make code cleaner define these here 
  invB       = soil%clappB(:,ms)                            !1 over C&H B
  Nsmpsat(:) = soil%smpsat(:,ms)                            !psi_saturated mm
  dzmm(:,:)  = spread((soil%zse(:)) * 1000.0,1,mp)    !layer thickness mm
  zimm(0)    = 0.0_r_2                                      !depth of layer interfaces mm
  zimm(1:ms) = zimm(0:ms-1) + real(dzmm(1,1:ms),r_2)
  
  defc(:) = (soil%watsat(:,ms))*(zimm(ms)+Nsmpsat(:)/(1-invB(:))* &
    (1.0_r_2-((Nsmpsat(:)+zimm(ms))/Nsmpsat(:))**(1.0_r_2-invB(:))))             !def if wtd=zimm(ms)
  where (defc(:) .le. 0.0_r_2) defc(:) = 0.1_r_2
  def(:) = sum((soil%watsat(:,:)-ssoil%wb(:,:))*dzmm(:,:),2)
     
  do i=1,mp
    keeplooping = .TRUE.
    if (defc(i) > def(i)) then                 !iterate tfor wtd
       jlp=0
       mainloop: DO WHILE (keeplooping)
          tempa   = 1.0_r_2
          tempb   = (1.0_r_2+ssoil%wtd(i)/Nsmpsat(i))**(-invB(i))
          derv    = (soil%watsat(i,ms))*(tempa-tempb) + &
                                          soil%watsat(i,ms)
          if (abs(derv) .lt. 1e-5) derv = sign(1e-5,derv)
          tempa   = 1.0_r_2
          tempb   = (1.0_r_2+ssoil%wtd(i)/Nsmpsat(i))**(1.0_r_2-invB(i))
          deffunc = (soil%watsat(i,ms))*(ssoil%wtd(i) +&
                     Nsmpsat(i)/(1.0_r_2-invB(i))* &
                        (tempa-tempb)) - def(i)
          calc    = ssoil%wtd(i) - deffunc/derv
          IF ((abs(calc-ssoil%wtd(i))) .le. wtd_uncert) THEN
              ssoil%wtd(i) = calc
              keeplooping = .FALSE.
          ELSEIF (jlp==mx_wtd_iterations) THEN
              keeplooping = .FALSE.
          ELSE
              jlp=jlp+1
              ssoil%wtd(i) = calc
          END IF
       END DO mainloop
    elseif (defc(i) .lt. def(i)) then
       jlp=0
       mainloop2: DO WHILE (keeplooping)
          tmpc     = Nsmpsat(i)+ssoil%wtd(i)-zimm(ms)
          if (abs(tmpc) .lt. 1e-4) tmpc = sign(1e-4,tmpc)
          tempa    = (abs(tmpc/Nsmpsat(i)))**(-invB(i))
          
          derv     = (soil%watsat(i,ms))*(tempa-tempb)
          if (abs(derv) .lt. 1e-1) derv = sign(1e-1,derv)
          
          tempa    = ((Nsmpsat(i)+ssoil%wtd(i)-zimm(ms))/Nsmpsat(i))**(1.0_r_2-invB(i))
          tempb    = (1.0_r_2+ssoil%wtd(i)/Nsmpsat(i))**(1.0-invB(i))
          deffunc  = (soil%watsat(i,ms))*(zimm(ms) +&
                      Nsmpsat(i)/(1.0_r_2-invB(i))*(tempa-tempb))-def(i)
          !calc     = ssoil%wtd(i) - deffunc/derv
          calc     = max(ssoil%wtd(i) - deffunc/derv,1.0)   !prevent wtd < 0
          IF ((abs(calc-ssoil%wtd(i))) .le. wtd_uncert) THEN
             ssoil%wtd(i) = calc
             keeplooping = .FALSE.
          ELSEIF (jlp==mx_wtd_iterations) THEN
             keeplooping = .FALSE.
          ELSE
             jlp=jlp+1
             ssoil%wtd(i) = calc
          END IF
       END DO mainloop2
    else
       ssoil%wtd(i) = zimm(ms)
    endif
    
  end do   !mp patches loop

  where (ssoil%wtd(:) .gt. wtd_max) ssoil%wtd(:) = wtd_max
  where (ssoil%wtd(:) .lt. wtd_min) ssoil%wtd(:) = wtd_min


  END SUBROUTINE calcwtd
    

  !-------------------------------------------------------------------------
  ! SUBROUTINE smoistgw (fwtop,dt,ktau,ssoil,soil,prin)
  ! solves the modified richards equation (Zeng and Decker 2009) to find
  ! vertical mocement of soil water.  Bottom boundary condition is determined
  ! using a single layer groundwater module
  !
  SUBROUTINE smoistgw (dels,ktau,ssoil,soil,prin)
  IMPLICIT NONE
  
    REAL, INTENT(IN)                          :: dels  ! time step size (s)
    INTEGER, INTENT(IN)                       :: ktau  ! integration step number
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssoil ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
    LOGICAL, INTENT(IN)                       :: prin
    
    !Local variables.  
    REAL(r_2), DIMENSION(mp,ms+1)       :: at     ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: bt     ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: ct     ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: rt
!    REAL(r_2), DIMENSION(mp)      :: fact
!    REAL(r_2), DIMENSION(mp)      :: fact2
!    REAL(r_2), DIMENSION(mp)      :: fluxhi
!    REAL(r_2), DIMENSION(mp)      :: fluxlo
    REAL(r_2), DIMENSION(mp)            :: hydss  ! hydraulic
                                                ! conductivity adjusted for ice
	
    INTEGER                                   :: k,kk
    REAL(r_2), DIMENSION(mp,ms)         :: eff_por,old_wb  !effective porosity and wb at start
    REAL(r_2), DIMENSION(mp)            :: den
    REAL(r_2), DIMENSION(mp)            :: dne
    REAL(r_2), DIMENSION(mp)            :: num
    REAL(r_2), DIMENSION(mp)            :: qin
    REAL(r_2), DIMENSION(mp)            :: qout
    REAL(r_2), DIMENSION(mp)            :: dqidw0
    REAL(r_2), DIMENSION(mp)            :: dqidw1
    REAL(r_2), DIMENSION(mp)            :: dqodw0
    REAL(r_2), DIMENSION(mp)            :: dqodw1,dqodw2
    REAL(r_2), DIMENSION(mp)            :: s1,s2,tmpi,temp0,voleq1,tempi
    REAL(r_2), DIMENSION(ms)                  :: dzmm
    REAL(r_2), DIMENSION(0:ms+1)                :: zimm
    REAL(r_2), DIMENSION(ms)                  :: zmm
    REAL(r_2), DIMENSION(mp)            :: GWzimm,xs,zaq,fice_avg,s_mid,GWdzmm
    REAL(r_2), DIMENSION(mp)            :: xsi,xs1,GWmsliq
    REAL(r_2), DIMENSION(mp,ms+1)       :: qhlev,del_wb
    INTEGER, DIMENSION(mp)              :: idlev
    REAL(r_2),DIMENSION(mp)             :: watmin
    REAL(r_2)                                 :: dri  !ratio of density of ice to density of water
                                                !ensures smp is independent of liq->ice conversion
    REAL(r_2), DIMENSION(mp,ms)         :: msliq,msice
    logical                                   :: prinall = .false.   !another debug flag
    character (len=30)                        :: fmt  !format to output some debug info
    
    fmt='(A6,6(1X,F8.6))'
 

    !make code cleaner define these here
    dzmm(:)    = 1000.0_r_2 * real(soil%zse(:),r_2)
    zimm(0)    = 0.0
    zmm(1)     = 0.5*dzmm(1)
    do k=1,ms
       zimm(k) = zimm(k-1) + dzmm(k)
    end do
    do k=2,ms
       zmm(k)  = zimm(k) + 0.5*dzmm(k)
    end do 
    GWdzmm(:) = real(1000.0*soil%GWdz(:),r_2)   
    GWzimm(:) = zimm(ms)+GWdzmm(:)
    zaq(:)    = real(1000.0*soil%GWz(:),r_2)
    
    
    watmin(:) = volwatmin * denliq

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! preset to allow for non-land & snow points in trimb
    at = 0.0
    bt = 1.0
    ct = 0.0
    old_wb(:,:) = ssoil%wb(:,:)
    
    !equilibrium water content
    do k=1,ms
       WHERE ((ssoil%wtd(:) .le. zimm(k-1)))          !fully saturated
          ssoil%wbeq(:,k) = real(soil%watsat(:,k)-soil%watr(:,k),r_2)
       ELSEWHERE
          !use the weighted average of sat & unsat
          WHERE ((ssoil%wtd(:) .le. zimm(k)) .and. (ssoil%wtd(:) .gt. zimm(k-1)))
             tempi = 1.0
             temp0 = (((soil%smpsat(:,k)+ssoil%wtd(:)-zimm(k-1))/soil%smpsat(:,k)))**(1.0-1.0/soil%clappB(:,k))               
             voleq1 = -soil%smpsat(:,k)*(soil%watsat(:,k)-soil%watr(:,k))/&
		         (1.0-1.0/soil%clappB(:,k))/(ssoil%wtd(:)-zimm(k-1))*(tempi-temp0)
             ssoil%wbeq(:,k) = (voleq1*(ssoil%wtd(:)-zimm(k-1)) + (soil%watsat(:,k)-soil%watr(:,k))&
		                  *(zimm(k)-ssoil%wtd(:)))/(zimm(k)-zimm(k-1)) + soil%watr(:,k)
             ssoil%wbeq(:,k) = min(soil%watsat(:,k),ssoil%wbeq(:,k))
             ssoil%wbeq(:,k) = max(ssoil%wbeq(:,k),0.01_r_2)
          ELSEWHERE 
             tempi = (((soil%smpsat(:,k)+ssoil%wtd(:)-zimm(k))/soil%smpsat(:,k)))**(1-1/soil%clappB(:,k))
             temp0 = (((soil%smpsat(:,k)+ssoil%wtd(:)-zimm(k-1))/soil%smpsat(:,k)))**(1-1/soil%clappB(:,k))   
             ssoil%wbeq(:,k) = -soil%smpsat(:,k)*(soil%watsat(:,k)-soil%watr(:,k))/&
		                 (1-1/soil%clappB(:,k))/(zimm(k)-zimm(k-1))*(tempi-temp0)+soil%watr(:,k)
             ssoil%wbeq(:,k) = max(ssoil%wbeq(:,k),0.01_r_2)
             ssoil%wbeq(:,k) = min(real(soil%watsat(:,k),r_2),ssoil%wbeq(:,k))
          END WHERE
       END WHERE
       ssoil%zq(:,k) = -soil%smpsat(:,k)*(max((ssoil%wbeq(:,k)-soil%watr(:,k))/&
	                real(soil%watsat(:,k)-soil%watr(:,k),r_2),0.01_r_2))**(-soil%clappB(:,k))
       ssoil%zq(:,k) = max(sucmin, ssoil%zq(:,k))
    end do
       
    !Aquifer Equilibrium water content
    WHERE (ssoil%wtd(:) .lt. zimm(ms))                                          !fully saturated
       ssoil%GWwbeq(:) = real(soil%GWwatsat(:),r_2)
    END WHERE
    WHERE ((ssoil%wtd(:) .gt. GWzimm(:)))                                       !fully unsaturated
       tempi = (((soil%GWsmpsat(:)+ssoil%wtd(:)-GWzimm(:))/soil%GWsmpsat(:)))**(1-1/soil%GWclappB(:))
       temp0 = (((soil%GWsmpsat(:)+ssoil%wtd(:)-zimm(ms))/soil%GWsmpsat(:)))**(1-1/soil%GWclappB(:))   
       ssoil%GWwbeq(:) = -soil%GWsmpsat(:)*soil%GWwatsat(:)/&
	                 (1-1/soil%GWclappB(:))/(GWzimm(:)-zimm(ms))*(tempi-temp0) + soil%GWwatr(:)
       ssoil%GWwbeq(:) = max(real(ssoil%GWwbeq(:),r_2),0.01_r_2)
       ssoil%GWwbeq(:) = min(real(soil%GWwatsat(:),r_2),ssoil%GWwbeq(:))	 
    END WHERE           
    WHERE ((ssoil%wtd(:) .lt. GWzimm(:)) .and. (ssoil%wtd(:) .gt. zimm(ms)))    !partially saturated
       tempi  = 1.0_r_2
       temp0  = (((soil%GWsmpsat(:)+ssoil%wtd(:)-zimm(ms))/soil%GWsmpsat(:)))**(1.0-1.0/soil%GWclappB(:))               
       voleq1 = -soil%GWsmpsat(:)*(soil%GWwatsat(:)-soil%GWwatr(:))/&
	          (1.0-1.0/soil%GWclappB(:))/(ssoil%wtd(:)-zimm(ms))*(tempi-temp0) + soil%GWwatr(:)
       ssoil%GWwbeq(:) = (voleq1*(ssoil%wtd(:)-zimm(ms)) + (soil%GWwatsat(:)-soil%GWwatr(:))*&
	                   (GWzimm(:)-ssoil%wtd(:)))/(GWzimm(:)-zimm(ms)) + soil%GWwatr(:)
       ssoil%GWwbeq(:) = min(real(soil%GWwatsat(:),r_2),ssoil%GWwbeq(:))
       ssoil%GWwbeq(:) = max(ssoil%GWwbeq(:),0.01_r_2)
    END WHERE
    ssoil%GWzq(:) = -soil%GWsmpsat(:)*(max((ssoil%GWwbeq(:)-soil%GWwatr(:))/&
                      (soil%GWwatsat(:)-soil%GWwatr(:)),0.01_r_2))**(-soil%GWclappB(:))
    ssoil%GWzq(:) = max(sucmin, ssoil%GWzq(:))
              
    !soil matric potential, hydraulic conductivity, and derivatives of each with respect to water
    do k=1,ms
       ssoil%icefrac(:,k) = ssoil%wbice(:,k)/(max(ssoil%wbice(:,k)+ssoil%wb(:,k),0.01_r_2))
       ssoil%fracice(:,k) = (ssoil%icefrac(:,k)- exp(-3.0))/(1.0-exp(-3.0))
       ssoil%fracice(:,k) = max(min(ssoil%fracice(:,k),1.0_r_2),0.0_r_2)
    end do
    fice_avg(:)  = sum(ssoil%fracice(:,:)*spread(dzmm(:),1,mp),2) / sum(dzmm(:))
    fice_avg(:)  = min(max(fice_avg(:) / ms,0.0_r_2),0.95_r_2)
    where(fice_avg(:) < ssoil%fracice(:,ms)) fice_avg(:) = ssoil%fracice(:,ms)         !frozen ms limits qh
       
    do k=1,ms   
       kk=min(k+1,ms)
       if (k < ms) then
          s1 = min(0.5*(ssoil%wb(:,k)-soil%watr(:,k) + ssoil%wb(:,kk)-soil%watr(:,kk)) / &
               (0.5*(soil%watsat(:,k)-soil%watr(:,k) + soil%watsat(:,kk)-soil%watr(:,kk))),1.0_r_2)
          s2 = soil%hksat(:,k)*s1**(2.0*soil%clappB(:,k)+2.0)
          ssoil%hk(:,k)    = (1.0-0.5*(ssoil%fracice(:,k)+ssoil%fracice(:,kk)))*s1*s2
          ssoil%dhkdw(:,k) = (1.0-0.5*(ssoil%fracice(:,k)+ssoil%fracice(:,kk)))* &
                    (2.0*soil%clappB(:,k)+3.0)*s2*0.5/(soil%watsat(:,k)-soil%watr(:,k))
       else
          s1 = min(0.5*(ssoil%wb(:,k)-soil%watr(:,k) + ssoil%GWwb(:)-soil%GWwatr(:)) / &
               (0.5*(soil%watsat(:,k)-soil%watr(:,k) + soil%GWwatsat(:)-soil%GWwatr(:))),1.0_r_2)
          s2 = soil%hksat(:,k)*s1**(2.0*soil%clappB(:,k)+2.0)
          ssoil%hk(:,k)    = s1*s2*(1.0-ssoil%fracice(:,k))* exp (-hkrz*zimm(ms)/1000.0_r_2)
          ssoil%dhkdw(:,k) = (1.0-ssoil%fracice(:,k))* (2.0*soil%clappB(:,k)+3.0)*&
	                         s2*0.5/(soil%watsat(:,k)-soil%watr(:,k)*exp (hkrz*zimm(ms)/1000.0_r_2))
       end if
	  
       s_mid = max(min((ssoil%wb(:,k)+dri*ssoil%wbice(:,k)-soil%watr(:,k))/&
              (soil%watsat(:,k)-soil%watr(:,k)),1.0_r_2),0.01_r_2)
       ssoil%smp(:,k) = max(min(-soil%smpsat(:,k)*s_mid**(-soil%clappB(:,k)),&
	                         -soil%smpsat(:,k)),sucmin)
       ssoil%dsmpdw(:,k) = -soil%clappB(:,k)*ssoil%smp(:,k)/&
	                    (max(s_mid*(soil%watsat(:,k)-soil%watr(:,k)),0.001_r_2))          
    end do
    !Aquifer properties
    s_mid = (ssoil%GWwb(:)-soil%GWwatr(:))/(soil%GWwatsat(:)-soil%GWwatr(:))
    s_mid = max(min(s_mid,1.0_r_2),0.01_r_2)   
    s2    = soil%GWhksat(:)*s_mid**(2.0*soil%clappB(:,ms)+2.0)
    ssoil%GWhk(:)     = s_mid*s2*(1.0-ssoil%fracice(:,k))

    ssoil%GWdhkdw(:)  = (1.0-ssoil%fracice(:,ms))* (2.0*soil%clappB(:,ms)+3.0)*&
                         s2*0.5/(soil%GWwatsat(:)-soil%GWwatr(:))
    ssoil%GWsmp(:)    = -soil%smpsat(:,ms)*s_mid**(-soil%clappB(:,ms))
    ssoil%GWdsmpdw(:) = -soil%clappB(:,ms)*ssoil%GWsmp(:)/(s_mid*(soil%GWwatsat(:)-soil%GWwatr(:)))
       
    !Note: temporaary parameteriation of horizontal drainage
    !too be replaced with explivit treatment of subgrid scale, topographically
    !based subsurface flux convergence flowing to river channels 
       
    ssoil%qhz(:)  = qhmax *exp(-8.0_r_2*ssoil%wtd(:)/1000.0)*((1.0_r_2 - fice_avg(:))**3.0)
    !find index of soil layer with the water table
    qhlev(:,:)   = 0.0  !set to zero except for layer that contains the wtd
    idlev(:)     = ms
    do k=1,ms
       WHERE ((ssoil%wtd(:) > zimm(k-1)) .and. (ssoil%wtd(:) <= zimm(k)))
         idlev(:) = k
         qhlev(:,k) = ssoil%qhz(:)
       END WHERE
    end do  
    WHERE (ssoil%wtd(:) > zimm(ms))
       idlev(:) = ms+1
       qhlev(:,ms+1) = ssoil%qhz(:)
    END WHERE
       
    rt(:,:) = 0.0_r_2; at(:,:) = 0.0_r_2     !ensure input to tridiag is valid
    bt(:,:) = 1.0_r_2; ct(:,:) = 0.0_r_2
	
    k = 1
       qin    = ssoil%fwtop(:)/dels
       den    = (zmm(k+1)-zmm(k))
       dne    = (ssoil%zq(:,k+1)-ssoil%zq(:,k))
       num    = (ssoil%smp(:,k+1)-ssoil%smp(:,k)) - dne
       qout   = -ssoil%hk(:,k)*num/den
       dqodw1 = -(-ssoil%hk(:,k)*ssoil%dsmpdw(:,k)   + num*ssoil%dhkdw(:,k))/den
       dqodw2 = -( ssoil%hk(:,k)*ssoil%dsmpdw(:,k+1) + num*ssoil%dhkdw(:,k))/den
       rt(:,k) =  qin - qout  - qhlev(:,k)! - ssoil%rex(:,k)
       at(:,k) =  0.0_r_2
       bt(:,k) =  dzmm(k)/dels + dqodw1
       ct(:,k) =  dqodw2      

    do k = 2, ms - 1
       den    = (zmm(k) - zmm(k-1))
       dne    = (ssoil%zq(:,k)-ssoil%zq(:,k-1))
       num    = (ssoil%smp(:,k)-ssoil%smp(:,k-1)) - dne
       qin    = -ssoil%hk(:,k-1)*num/den
       dqidw0 = -(-ssoil%hk(:,k-1)*ssoil%dsmpdw(:,k-1) + num*ssoil%dhkdw(:,k-1))/den
       dqidw1 = -( ssoil%hk(:,k-1)*ssoil%dsmpdw(:,k)   + num*ssoil%dhkdw(:,k-1))/den
       den    = (zmm(k+1)-zmm(k))
       dne    = (ssoil%zq(:,k+1)-ssoil%zq(:,k))
       num    = (ssoil%smp(:,k+1)-ssoil%smp(:,k)) - dne
       qout   = -ssoil%hk(:,k)*num/den
       dqodw1 = -(-ssoil%hk(:,k)*ssoil%dsmpdw(:,k)   + num*ssoil%dhkdw(:,k))/den
       dqodw2 = -( ssoil%hk(:,k)*ssoil%dsmpdw(:,k+1) + num*ssoil%dhkdw(:,k))/den
       rt(:,k) =  qin - qout  - qhlev(:,k)! - ssoil%rex(:,k)
       at(:,k) = -dqidw0
       bt(:,k) =  dzmm(k)/dels - dqidw1 + dqodw1
       ct(:,k) =  dqodw2
    end do
       
    k = ms
       den    = (zmm(k) - zmm(k-1))
       dne    = (ssoil%zq(:,k)-ssoil%zq(:,k-1))
       num    = (ssoil%smp(:,k)-ssoil%smp(:,k-1)) - dne
       qin    = -ssoil%hk(:,k-1)*num/den
       dqidw0 = -(-ssoil%hk(:,k-1)*ssoil%dsmpdw(:,k-1) + num*ssoil%dhkdw(:,k-1))/den
       dqidw1 = -( ssoil%hk(:,k-1)*ssoil%dsmpdw(:,k)   + num*ssoil%dhkdw(:,k-1))/den
       den    = zaq(:) - zmm(k)!dzmm(ms)
       dne    = (ssoil%GWzq(:)-ssoil%zq(:,k))
       num    =  (ssoil%GWsmp(:)-ssoil%smp(:,k)) - dne
       qout   = -ssoil%hk(:,k)*num/den
       dqodw1 = -(-ssoil%hk(:,k)*ssoil%dsmpdw(:,k)   + num*ssoil%dhkdw(:,k))/den
       dqodw2 = -( ssoil%hk(:,k)*ssoil%GWdsmpdw(:) + num*ssoil%dhkdw(:,k))/den
       rt(:,k) =  qin - qout  - qhlev(:,k) !- ssoil%rex(:,k)
       at(:,k) = -dqidw0
       bt(:,k) =  dzmm(k)/dels - dqidw1 + dqodw1
       ct(:,k) =  dqodw2
          
    k = ms+1
       den    = (zaq(:) - zmm(k-1))
       dne    = (ssoil%GWzq(:)-ssoil%zq(:,k-1))
       num    = (ssoil%GWsmp(:)-ssoil%smp(:,k-1)) - dne
       qin    = -ssoil%hk(:,k-1)*num/den
       dqidw0 = -(-ssoil%hk(:,k-1)*ssoil%dsmpdw(:,k-1) + num*ssoil%dhkdw(:,k-1))/den
       dqidw1 = -( ssoil%hk(:,k-1)*ssoil%GWdsmpdw(:)   + num*ssoil%dhkdw(:,k-1))/den
       den    = zaq(:) - zmm(k-1)!dzmm(ms)
       dne    = (ssoil%GWzq(:)-ssoil%zq(:,k-1))
       num    =  (ssoil%GWsmp(:)-ssoil%smp(:,k-1)) - dne
       qout   = 0.0_r_2
       dqodw1 = 0.0_r_2
       dqodw2 = 0.0_r_2
       rt(:,k) =  qin - qout  - qhlev(:,k)
       at(:,k) = -dqidw0
       bt(:,k) =  GWdzmm(:)/dels - dqidw1
       ct(:,k) =  0.0_r_2
          
    CALL solve_tridiag(at, bt, ct, rt, del_wb,ms+1)                      !solve system of eqns

    ssoil%GWwb(:) = ssoil%GWwb(:) + del_wb(:,ms+1)                       !add del-h2o to soil  
    do k=1,ms 
       ssoil%wb(:,k) = old_wb(:,k) + del_wb(:,k)
       msliq(:,k)    = (ssoil%wb(:,k))*dzmm(k)                           !mass of soil liq [mm] from volumetric
       msice(:,k)    = denice*ssoil%wbice(:,k)*dzmm(k)/denliq            !mass of soil ice
       eff_por(:,k)  = (soil%watsat(:,k)) - ssoil%wbice(:,k)
    end do
    GWmsliq(:) = ssoil%GWwb(:)*GWdzmm                                    !mass aquifer liq 
          
    xsi(:)      = GWmsliq(:) - soil%GWwatsat(:)*GWdzmm
    where (xsi(:) .lt. 0.0_r_2) xsi(:) = 0.0_r_2
    where (xsi(:) .gt. 0.0_r_2) GWmsliq(:) = soil%GWwatsat(:)*GWdzmm(:) 	  
    msliq(:,ms) = msliq(:,ms) + xsi(:)
 
    do k = ms,2,-1
       where ((xsi(:) .lt. (msliq(:,k)-eff_por(:,k)*dzmm(k))) .and. (xsi(:) .gt. 0.0_r_2))
          msliq(:,k) = msliq(:,k) + xsi(:)
       elsewhere (xsi(:) .gt. 0.0_r_2)
          xsi(:)     = xsi(:) - (eff_por(:,k)*dzmm(k)-msliq(:,k))
          msliq(:,k) = eff_por(:,k)*dzmm(k)
       end where
       !xsi(:)       = max(msliq(:,k)-eff_por(:,k)*dzmm(k),0._r_2)
       !msliq(:,k)   = min(eff_por(:,k)*dzmm(k), msliq(:,k))
       !msliq(:,k-1) = msliq(:,k-1) + xsi(:)
    end do
    xs1(:) = msliq(:,1)-eff_por(:,1)*dzmm(1)
    where(xs1(:) .lt. 0.0_r_2) xs1(:) = 0.0_r_2

    where(xs1(:) .gt. 0.0_r_2) 
       msliq(:,1) = eff_por(:,1)*dzmm(1)
       ssoil%qhz(:) = ssoil%qhz(:) + xs1(:) / dt
       xs1(:) = 0.0_r_2
   end where
 
    do k = 1,ms-1                                                           !ensure liq < liq_minimum (using mm)
       xs(:) = 0.0_r_2
       where (msliq(:,k) .lt. watmin)
          xs(:) = watmin - msliq(:,k)
       end where
       msliq(:,k  ) = msliq(:,k  ) + xs(:)
       if (k .le. (ms-1)) then
          msliq(:,k+1) = msliq(:,k+1) - xs(:)
       else
          GWmsliq(:) = GWmsliq(:) - xs(:)
       endif
    end do

    ssoil%wb(:,:) = msliq(:,:) / spread(dzmm(:),1,mp_patch)                 !convert from mm to volumetric
    ssoil%GWwb(:) = GWmsliq(:) / GWdzmm(:)  

    ssoil%rnof2(:) = ssoil%qhz(:)          
          
          


 END SUBROUTINE smoistgw



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
!	 ssnow
SUBROUTINE soil_snow_gw(dels, soil, ssnow, canopy, met, bal, veg)
   USE cable_common_module
   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(soil_parameter_type), INTENT(INOUT) :: soil
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(canopy_type), INTENT(INOUT)         :: canopy
   TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
   TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
   TYPE (balances_type), INTENT(INOUT)      :: bal
   INTEGER             :: k
   REAL, DIMENSION(mp) :: snowmlt
   REAL, DIMENSION(mp) :: totwet
   REAL, DIMENSION(mp) :: weting
   REAL, DIMENSION(mp) :: xxx, tgg_old, tggsn_old
   REAL(r_2), DIMENSION(mp) :: xx,deltat,sinfil1,sinfil2,sinfil3 
   REAL                :: zsetot
   INTEGER, SAVE :: ktau =0 
   LOGICAL :: prin
   
   prin = .FALSE.
   
   CALL point2constants( C ) 
   cp = C%CAPP
    
   ktau = ktau +1 
   
   !jhan - make switchable 
   ! appropriate for ACCESS1.0
   !max_glacier_snowd = 50000.0
   ! appropriate for ACCESS1.3 
   max_glacier_snowd = 1100.0

   zsetot = sum(soil%zse) 
   ssnow%tggav = 0.
   DO k = 1, ms
      ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot
   END DO


   IF( cable_runtime%offline .or. cable_runtime%mk3l ) THEN
        ssnow%t_snwlr = 0.05
   ENDIF

   ssnow%fwtop1 = 0.0
   ssnow%fwtop2 = 0.0
   ssnow%fwtop3 = 0.0
   ssnow%runoff = 0.0 ! initialise total runoff
   ssnow%rnof1 = 0.0 ! initialise surface runoff
   ssnow%rnof2 = 0.0 ! initialise deep drainage
   ssnow%smelt = 0.0 ! initialise snowmelt
   ssnow%dtmlt = 0.0 
   ssnow%osnowd = ssnow%snowd


   IF( .NOT.cable_user%cable_runtime_coupled ) THEN
   
      IF( ktau_gl <= 1 ) THEN
         
         IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
                                                  ! after discussion with BP
         ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
         ssnow%wbtot = 0.0
         DO k = 1, ms
            ssnow%wb(:,k)  = MIN( soil%ssat,MAX ( ssnow%wb(:,k), soil%swilt ))
         END DO
   
         ssnow%wb(:,ms-2)  = MIN( soil%ssat, MAX ( ssnow%wb(:,ms-2),           &
                             0.5 * ( soil%sfc + soil%swilt ) ) )
         ssnow%wb(:,ms-1)  = MIN( soil%ssat, MAX ( ssnow%wb(:,ms-1),           &
                             0.8 * soil%sfc ) )
         ssnow%wb(:,ms)    = MIN( soil%ssat, MAX ( ssnow%wb(:,ms), soil%sfc) )
         
         DO k = 1, ms
            
            WHERE( ssnow%tgg(:,k) <= C%TFRZ .AND. ssnow%wbice(:,k) <= 0.01 )   &
               ssnow%wbice(:,k) = 0.5 * ssnow%wb(:,k)
            
            WHERE( ssnow%tgg(:,k) < C%TFRZ)                                    &
               ssnow%wbice(:,k) = frozen_limit * ssnow%wb(:,k)
            
         END DO
   
         WHERE (soil%isoilm == 9) 
            ! permanent ice: fix hard-wired number in next version
            ssnow%snowd = max_glacier_snowd
            ssnow%osnowd = max_glacier_snowd
            ssnow%tgg(:,1) = ssnow%tgg(:,1) - 1.0
            ssnow%wb(:,1) = 0.95 * soil%ssat
            ssnow%wb(:,2) = 0.95 * soil%ssat
            ssnow%wb(:,3) = 0.95 * soil%ssat
            ssnow%wb(:,4) = 0.95 * soil%ssat
            ssnow%wb(:,5) = 0.95 * soil%ssat
            ssnow%wb(:,6) = 0.95 * soil%ssat
            ssnow%wbice(:,1) = 0.90 * ssnow%wb(:,1)
            ssnow%wbice(:,2) = 0.90 * ssnow%wb(:,2)
            ssnow%wbice(:,3) = 0.90 * ssnow%wb(:,3)
            ssnow%wbice(:,4) = 0.90 * ssnow%wb(:,4)
            ssnow%wbice(:,5) = 0.90 * ssnow%wb(:,5)
            ssnow%wbice(:,6) = 0.90 * ssnow%wb(:,6)
         ENDWHERE
         
         xx=soil%css * soil%rhosoil
         
         ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
              & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * cswat * rhowat &
              & + ssnow%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1)
      END IF
   ENDIF  ! if(.NOT.cable_runtime_coupled)

   xx=soil%css * soil%rhosoil
   IF (ktau <= 1)                                                              &
     ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil      &
            & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * cswat * rhowat           &
            & + ssnow%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1) +   &
            & (1. - ssnow%isflag) * cgsnow * ssnow%snowd



   DO k = 1, ms ! for stempv
      
      ! Set liquid soil water fraction (fraction of saturation value):
      ssnow%wblf(:,k) = MAX( 0.01_r_2, (ssnow%wb(:,k) - ssnow%wbice(:,k)) )    &
           & / REAL(soil%ssat,r_2)
      
      ! Set ice soil water fraction (fraction of saturation value):
      ssnow%wbfice(:,k) = REAL(ssnow%wbice(:,k)) / soil%ssat
   END DO
 
   CALL snowcheck (dels, ssnow, soil, met )

   CALL snowdensity (dels, ssnow, soil)

   CALL snow_accum (dels, canopy, met, ssnow, soil )

   CALL snow_melting (dels, snowmlt, ssnow, soil )

   ! Add snow melt to global snow melt variable:
   ssnow%smelt = snowmlt

   ! Adjust levels in the snowpack due to snow accumulation/melting,
   ! snow aging etc...
   CALL snowl_adjust(dels, ssnow, canopy )

   !!!CALL stempv(dels, canopy, ssnow, soil)
   !do the soil and snow melting, freezing prior to water movement
   
   ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

   CALL snow_melting (dels, snowmlt, ssnow, soil )
   
   ! Add new snow melt to global snow melt variable: 
   ssnow%smelt = ssnow%smelt + snowmlt

   !CALL remove_trans(dels, soil, ssnow, canopy, veg)

   CALL  soilfreeze(dels, soil, ssnow)


   ssnow%fwtop = canopy%precis + ssnow%smelt   !water for infiltration   
   

    CALL calcwtd (ssnow, soil,prin)                         !update the wtd
    CALL ovrlndflx (dels, ktau, ssnow, soil, prin )         !surface runoff, incorporate ssnow%pudsto?

    CALL smoistgw (dels,ktau,ssnow,soil,prin)               !vertical soil moisture movement.  LIS verions trans is rex removed here
    !note: canopy%segg appears to be soil evap.
    !canopy%fesp/C%HL*dels is the puddle evaporation
    
    ssnow%runoff = (ssnow%rnof1 + ssnow%rnof2)*dels          !total runoff
  
    ! Scaling  runoff to kg/m^2/s (mm/s) to match rest of the model
    ssnow%sinfil = 0.0
    ! lakes: replace hard-wired vegetation number in next version
    WHERE( veg%iveg == 16 )
      ssnow%sinfil = MIN( ssnow%rnof1, ssnow%wb_lake + MAX( 0.,canopy%segg ) )
      ssnow%rnof1 = MAX( 0.0, ssnow%rnof1 - ssnow%sinfil )
      ssnow%wb_lake = ssnow%wb_lake - ssnow%sinfil
      ssnow%rnof2 = MAX( 0.0, ssnow%rnof2 - ssnow%wb_lake )
    ENDWHERE    
    
    CALL remove_trans(dels, soil, ssnow, canopy, veg)        !transpiration loss per soil layer
    
    !test is top layer fully sat then put into pond?
    

!    ! total available liquid including puddle
!    weting = totwet + max(0.,ssnow%pudsto - canopy%fesp/C%HL*dels) 
!    xxx=soil%ssat - ssnow%wb(:,1)
!   
!    sinfil1 = MIN( 0.95*xxx*soil%zse(1)*rhowat, weting) !soil capacity
!    xxx=soil%ssat - ssnow%wb(:,2)
!    sinfil2 = MIN( 0.95*xxx*soil%zse(2)*rhowat, weting - sinfil1) !soil capacity
!    xxx=soil%ssat - ssnow%wb(:,3)
!    sinfil3 = MIN( 0.95*xxx*soil%zse(3)*rhowat,weting-sinfil1-sinfil2)
!    
!    ! net water flux to the soil
!    ssnow%fwtop1 = sinfil1 / dels - canopy%segg          
!    ssnow%fwtop2 = sinfil2 / dels           
!    ssnow%fwtop3 = sinfil3 / dels           
! 
!    ! Puddle for the next time step
!    ssnow%pudsto = max( 0., weting - sinfil1 - sinfil2 - sinfil3 )
!    ssnow%rnof1 = max(0.,ssnow%pudsto - ssnow%pudsmx)
!    ssnow%pudsto = ssnow%pudsto - ssnow%rnof1
! 
!    CALL surfbv(dels, met, ssnow, soil, veg, canopy )

   ! correction required for energy balance in online simulations 
   IF( cable_runtime%um) THEN
      canopy%fhs_cor = ssnow%dtmlt(:,1)*ssnow%dfh_dtg
      canopy%fes_cor = ssnow%dtmlt(:,1)*(ssnow%cls*ssnow%dfe_ddq * ssnow%ddq_dtg)

      canopy%fhs = canopy%fhs+canopy%fhs_cor
      canopy%fes = canopy%fes+canopy%fes_cor
   ENDIF

   ! redistrb (set in cable.nml) by default==.FALSE. 
   !IF( redistrb )                                                              &
   !   CALL hydraulic_redistribution( dels, soil, ssnow, canopy, veg, met )

   ssnow%smelt = ssnow%smelt/dels

   ! Set weighted soil/snow surface temperature
   ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

   ssnow%wbtot = 0._r_2
   DO k = 1, ms
      ssnow%wbtot = ssnow%wbtot + REAL(ssnow%wb(:,k)*1000.0*soil%zse(k),r_2)
   END DO

END SUBROUTINE soil_snow_gw







END MODULE cable_soil_snow_gw_module
