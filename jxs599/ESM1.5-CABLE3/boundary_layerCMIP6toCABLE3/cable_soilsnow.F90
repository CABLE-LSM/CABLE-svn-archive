!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
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
!
! ==============================================================================

MODULE cbl_soil_snow_main_module
   
   
   USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
                             veg_parameter_type, canopy_type, met_type,        &
                             balances_type, r_2, ms, mp           
   USE cable_data_module, ONLY : issnow_type, point2constants
USE trimb_mod,                    ONLY:trimb 

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
   
   REAL :: cp    ! specific heat capacity for air
   
   !jhan:make parameter
   REAL :: max_glacier_snowd
 
   ! This module contains the following subroutines:
   PUBLIC soil_snow ! must be available outside this module
   PRIVATE smoisturev, stempv

CONTAINS

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
         IF (.not. cable_user%l_new_runoff_speed) THEN
            wbl_k = MAX( 0.01_r_2, ssnow%wb(:,k) - ssnow%wbice(:,k) )
            wbl_kp = MAX( 0.01_r_2, ssnow%wb(:,k+1) - ssnow%wbice(:,k+1) )
         ELSE
            wbl_k = MAX( 0.001_r_2, ssnow%wb(:,k) - ssnow%wbice(:,k) )
            wbl_kp = MAX( 0.001_r_2, ssnow%wb(:,k+1) - ssnow%wbice(:,k+1) )
         ENDIF
         
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

      IF (.not. cable_user%l_new_runoff_speed) then

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

      ELSE

      WHERE( ssnow%wb(:,ms) > soil%sfc(:) )

         wbl_k = MAX( 0.001_r_2, ssnow%wb(:,ms) - ssnow%wbice(:,ms) )
         wbl_kp = MAX( 0.001_r_2, soil%ssat(:) - ssnow%wbice(:,ms) )

         wh = MIN( wbl_k, wbl_kp )

         WHERE( ssnow%wbice(:,ms) .GT. 0.05 ) wh = 0.9 * wbl_k + 0.1 * wbl_kp

         ! Calculate hyd conductivity adjusted for ice:
         hydss = soil%hyds

         speed_k = hydss * ( wh / soil%ssat )**( soil%i2bp3 - 1 )
         speed_k =  speed_k / ( 1. - MIN( 0.5, 10. * ssnow%wbice(:,ms) ) )
         fluxlo = wbl_k

         ! scale speed to grid lengths per dt & limit speed for stability
         speed_k = MIN( speed_k, 0.5 * soil%zse(ms) / dels )
         fluxh(:,ms) = MAX( 0.0, speed_k * fluxlo )

      END WHERE

      ENDIF

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
    
   REAL :: wb_lake_T, rnof2_T, ratio
   INTEGER :: k,j

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

!  Rescale drainage to remove water added to lakes (wb_lake) 
   ssnow%sinfil = 0.0
   WHERE( veg%iveg == 16 )
      ssnow%sinfil  = MIN( ssnow%rnof1, ssnow%wb_lake ) ! water that can be extracted from the rnof1
      ssnow%rnof1   = MAX( 0.0, ssnow%rnof1 - ssnow%sinfil )
      ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
      ssnow%sinfil  = MIN( ssnow%rnof2, ssnow%wb_lake ) ! water that can be extracted from the rnof2
      ssnow%rnof2   = MAX( 0.0, ssnow%rnof2 - ssnow%sinfil )
      ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
      xxx = MAX(0.0, (ssnow%wb(:,ms) - soil%sfc(:))*soil%zse(ms)*1000.0)
      ssnow%sinfil  = MIN( xxx, ssnow%wb_lake )
      ssnow%wb(:,ms) = ssnow%wb(:,ms) - ssnow%sinfil / (soil%zse(ms)*1000.0)
      ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
      xxx = MAX(0.0, (ssnow%wb(:,ms) - .5*(soil%sfc + soil%swilt))*soil%zse(ms)*1000.0)
      ssnow%sinfil  = MIN( xxx, ssnow%wb_lake )
      ssnow%wb(:,ms) = ssnow%wb(:,ms) - ssnow%sinfil / (soil%zse(ms)*1000.0)
      ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ssnow%sinfil)
   ENDWHERE

   !wb_lake_T = sum( ssnow%wb_lake )
   !rnof2_T = sum( ssnow%rnof2 )
   !ratio = min( 1., wb_lake_T/max(rnof2_T,1.))
   !ssnow%rnof2 = ssnow%rnof2 - ratio*ssnow%rnof2
   !ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ratio*ssnow%rnof2)

!  Rescale drainage to remove water added to lakes (wb_lake)
   !wb_lake_T = 0.0
   !rnof2_T = 0.
   !DO j=1,mp
   !   IF( ssnow%wb_lake(j) >  0.0 ) wb_lake_T = wb_lake_T + ssnow%wb_lake(j)
   !   rnof2_T = rnof2_T + ssnow%rnof2(j)
   !END DO
   !ratio = min( 1., wb_lake_T/max(rnof2_T,1.))
   !ssnow%rnof2 = ssnow%rnof2 - ratio*ssnow%rnof2
   !ssnow%wb_lake = MAX( 0.0, ssnow%wb_lake - ratio*ssnow%rnof2)

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
                      * REAL( canopy%dgdtg ) ) * dels / REAL( ssnow%gammzz(:,1) )
   END WHERE
   
   coeff(:,1-3) = 0.0  ! coeff(:,-2)

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
   
   ssnow%tggsn = REAL( tmp_mat(:,1:3) )
   ssnow%tgg   = REAL( tmp_mat(:,4:(ms+3)) )
   canopy%sghflux = coefa * ( ssnow%tggsn(:,1) - ssnow%tggsn(:,2) )
   canopy%ghflux = coefb * ( ssnow%tgg(:,1) - ssnow%tgg(:,2) ) ! +ve downwards

END SUBROUTINE stempv

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
 !called subrs
USE hydraulic_redistribution_mod, ONLY: hydraulic_redistribution
USE soilfreeze_mod,               ONLY: soilfreeze
USE remove_trans_mod,             ONLY: remove_trans
USE snowl_adjust_mod,             ONLY: snowl_adjust
USE snowCheck_mod,                ONLY: snowCheck
USE snow_melting_mod,             ONLY: snow_melting
USE snow_accum_mod,               ONLY: snow_accum
USE snowdensity_mod,              ONLY: snowDensity
USE trimb_mod,                    ONLY:trimb 
!USE stempv_mod,                   ONLY: stempv
!USE surfbv_mod,                   ONLY: surfbv
 
USE cbl_ssnow_data_mod, ONLY : heat_cap_lower_limit

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

   CALL  soilfreeze(dels, soil, ssnow,heat_cap_lower_limit)


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
   IF( cable_runtime%um ) THEN
      canopy%fhs_cor = ssnow%dtmlt(:,1)*ssnow%dfh_dtg
      canopy%fes_cor = ssnow%dtmlt(:,1)*(ssnow%dfe_ddq * ssnow%ddq_dtg)

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
END MODULE cbl_soil_snow_main_module
