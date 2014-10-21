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

   USE cable_common_module, ONLY : gw_params

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
      
   !mrd561 GW params
   !Should read some in from namelist
   REAL(r_2), PARAMETER :: sucmin  = -10000000.0, &! minimum soil pressure head [mm]
                      hkrz         = 1.0,         &! GW_hksat e-folding depth [m**-1]
                      volwatmin    = 0.00001,        &!min soil water [mm]      
                      wtd_uncert   = 0.1,         &! uncertaintiy in wtd calcultations [mm]
                      wtd_max      = 100000.0,    &! maximum wtd [mm]
                      wtd_min      = 100.0,       &! minimum wtd [mm]
                      dri          = 1.0             !ratio of density of ice to density of liquid [unitless]

   INTEGER, PARAMETER :: wtd_iter_max = 10 ! maximum number of iterations to find the water table depth                    
   
   REAL :: cp    ! specific heat capacity for air
   
   !jhan:make parameter
   REAL :: max_glacier_snowd
 
   ! This module contains the following subroutines:
   PUBLIC soil_snow_gw,calc_srf_wet_fraction,calc_equilibrium_water_content ! must be available outside this module
   PRIVATE snowdensity, snow_melting, snowcheck, snowl_adjust 
   PRIVATE trimb,snow_accum, stempv
   PRIVATE soilfreeze, remove_trans,calcwtd
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
   REAL, DIMENSION(mp)                :: xx, ice_mass,liq_mass,tot_mass
   INTEGER k

   xx = 0.
   DO k = 1, ms

      ice_mass = ssnow%wbice(:,k)*real(soil%zse(k)*C%denice,r_2)
      liq_mass = ssnow%wbliq(:,k)*real(soil%zse(k)*C%denliq,r_2)
      tot_mass = liq_mass + ice_mass
      
      WHERE (ssnow%tgg(:,k) < C%TFRZ &
          & .AND. frozen_limit * ssnow%wb(:,k) - ssnow%wbice(:,k) > .001)
         
         sicefreeze = MIN( MAX( 0.0_r_2, ( frozen_limit * ssnow%wb(:,k) -      &
                      ssnow%wbice(:,k) ) ) * soil%zse(k) * C%denice,             &
                      ( C%TFRZ - ssnow%tgg(:,k) ) * ssnow%gammzz(:,k) / C%HLF )
         ssnow%wbice(:,k) = MIN( ssnow%wbice(:,k) + sicefreeze / (soil%zse(k)  &
                            * 1000.0), frozen_limit * ssnow%wb(:,k) )
         xx = soil%css * soil%rhosoil
         ssnow%gammzz(:,k) = MAX(                                              &
             REAL((1.0 - soil%ssat) * soil%css * soil%rhosoil ,r_2)            &
             + (ssnow%wb(:,k) - ssnow%wbice(:,k)) * REAL(cswat * rhowat,r_2)   &
             + ssnow%wbice(:,k) * REAL(csice * C%denice,r_2),              &
             REAL(xx,r_2)) * REAL( soil%zse(k),r_2 )

         WHERE (k == 1 .AND. ssnow%isflag == 0)
            ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + cgsnow * ssnow%snowd
         END WHERE
         ssnow%tgg(:,k) = ssnow%tgg(:,k) + REAL(sicefreeze)                    &
                          * C%HLF / REAL(ssnow%gammzz(:,k) )
      
      ELSEWHERE( ssnow%tgg(:,k) > C%TFRZ .AND. ssnow%wbice(:,k) > 0. )
         
         sicemelt = MIN( ssnow%wbice(:,k) * soil%zse(k) * C%denice,              &
                    ( ssnow%tgg(:,k) - C%TFRZ ) * ssnow%gammzz(:,k) / C%HLF )
         
         ssnow%wbice(:,k) = MAX( 0.0_r_2, ssnow%wbice(:,k) - sicemelt          &
                            / (soil%zse(k) * C%denice) )
         xx = soil%css * soil%rhosoil
         ssnow%gammzz(:,k) = MAX(                                              &
              REAL((1.0-soil%ssat) * soil%css * soil%rhosoil,r_2)             &
              + (ssnow%wb(:,k) - ssnow%wbice(:,k)) * REAL(cswat*rhowat,r_2)   &
              + ssnow%wbice(:,k) * REAL(csice * C%denice,r_2),            &
              REAL(xx,r_2) ) * REAL(soil%zse(k),r_2)
         WHERE (k == 1 .AND. ssnow%isflag == 0)
            ssnow%gammzz(:,k) = ssnow%gammzz(:,k) + cgsnow * ssnow%snowd
         END WHERE
         ssnow%tgg(:,k) = ssnow%tgg(:,k) - REAL(sicemelt)                     &
                          * C%HLF / REAL(ssnow%gammzz(:,k))
       
      END WHERE
      !update the liq and ice volume and mass
      ice_mass = ssnow%wbice(:,k)*real(soil%zse(k)*C%denice,r_2)
      liq_mass = tot_mass - ice_mass
      ssnow%wbliq(:,k) = liq_mass / real(soil%zse(k)*C%denliq,r_2)
      ssnow%wbice(:,k) = ice_mass / real(soil%zse(k)*C%denice,r_2)
      ssnow%wb(:,k)    = ssnow%wbliq(:,k) + ssnow%wbice(:,k)
    
   END DO

END SUBROUTINE soilfreeze

!-----------------------------------------------------------------------------------------

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
 
 
   xx = 0._r_2; xxd = 0._r_2; diff(:,:) = 0._r_2
   DO k = 1,ms
   
      ! Removing transpiration from soil:
      WHERE (canopy%fevc > 0._r_2 )     ! convert to mm/dels
      
         ! Calculate the amount (perhaps moisture/ice limited)
         ! which can be removed:
         xx = canopy%fevc * dels / C%HL * veg%froot(:,k) + diff(:,k-1)   ! kg/m2
         diff(:,k) = MAX( 0._r_2, ssnow%wbliq(:,k) - soil%swilt) &      ! m3/m3  
                     * real(soil%zse(k)*C%denliq,r_2)
         xxd = xx - diff(:,k)
       
         WHERE ( xxd .GT. 0._r_2 )
            ssnow%wbliq(:,k) = ssnow%wbliq(:,k) - diff(:,k)/real(soil%zse(k)*C%denliq,r_2)  !volume
            diff(:,k) = xxd
         ELSEWHERE
            ssnow%wbliq(:,k) = ssnow%wbliq(:,k) - xx/real(soil%zse(k)*C%denliq,r_2)  !volume
            diff(:,k) = 0._r_2
         END WHERE

         ssnow%wmliq(:,k) = ssnow%wbliq(:,k)*real(soil%zse(k)*C%denliq,r_2)  !mass
         ssnow%wmtot(:,k) = ssnow%wmliq(:,k) + ssnow%wmice(:,k)  !mass
         ssnow%wb(:,k)    = ssnow%wbliq(:,k) + ssnow%wbice(:,k)  !volume
     
     END WHERE
   
   END DO

END SUBROUTINE remove_trans 


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

    ut(:,1) = rt(:,1) / bet(:)

    do k = 2,n
       gam(:,k) = ct(:,k-1) / bet(:)
       bet(:)   = max(bt(:,k) - at(:,k) * gam(:,k),0.00001_r_2)
       ut(:,k)   = (rt(:,k) - at(:,k)*ut(:,k-1)) / bet(:)
    end do

    do k = n-1,1,-1
       ut(:,k) = ut(:,k) - gam(:,k+1) * ut(:,k+1)
    end do

  end subroutine solve_tridiag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  !-------------------------------------------------------------------------
  SUBROUTINE ovrlndflx (dels, ktau, ssnow, soil,md_prin )
USE cable_common_module

  IMPLICIT NONE
    REAL, INTENT(IN)                         :: dels ! integration time step (s)
    INTEGER, INTENT(IN)                      :: ktau ! integration step number
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    LOGICAL, INTENT(IN)                      :: md_prin
    INTEGER, PARAMETER                       :: ntest = 0 ! for snow diag prints
    INTEGER, PARAMETER                       :: nglacier = 2 ! 0 original, 1 off, 2 new Eva
    INTEGER                                  :: k
    REAL, DIMENSION(mp)                :: rnof5
    REAL, DIMENSION(mp)                :: sgamm
    REAL, DIMENSION(mp)                :: smasstot
    REAL, DIMENSION(mp,0:3)            :: smelt1                   !snow melt
    REAL(r_2), DIMENSION(mp)           :: xxx,fice,icef,efpor      !tmp vars, fraction of ice in gridcell
    REAL(r_2), DIMENSION(mp)           :: tmpa,tmpb,qinmax         !tmp vars, maximum infiltration [mm/s]
    REAL(r_2), DIMENSION(mp)           :: satfrac,wtd_meters       !saturated fraction of cell, wtd in m
    REAL(r_2), DIMENSION(mp,ms)        :: liqmass,icemass,totmass  !liquid mass,ice mass, total mass [mm]
    REAL(r_2), DIMENSION(mp,ms)        :: dzmm_mp                  !soil layer thickness [mm] spread over mp
    logical                                  :: prinall = .false.  !for debugging
    
   if (md_prin) write(*,*) 'inside ovrlndflux '   !MDeck

    !For now assume there is no puddle?
    !ssnow%pudsto = 0.0!1e-5

    dzmm_mp  = 1000._r_2 * spread(soil%zse,1,mp)
    
    icemass  = ssnow%wbice(:,:) * dzmm_mp * dri   !dri = denice/denliq
    liqmass  = (ssnow%wb-ssnow%wbice) * dzmm_mp
    totmass  = icemass + liqmass

    where (totmass .lt. real(1e-5,r_2)) totmass = real(1e-5,r_2) 

    efpor(:) = soil%watsat(:,1) - ssnow%wbice(:,1)

    where (efpor .lt. 0.01_r_2) efpor = 0.01_r_2

    !srf frozen fraction.  should be based on topography
    icef(:) = max(0._r_2,min(1._r_2,icemass(:,1) / totmass(:,1)))
    fice(:) = (exp(-3._r_2*(1._r_2-icef(:)))- exp(-3._r_2))/(1.0-exp(-3.0)) !this normalizes of from 0-1

    where (fice(:) .lt. 0.0_r_2) fice(:) = 0.0_r_2
    where (fice(:) .gt. 1.0_r_2) fice(:) = 1.0_r_2

    ! Saturated fraction
    wtd_meters = ssnow%wtd / 1000.0_r_2

    satfrac(:) = (1._r_2-fice(:))*gw_params%MaxSatFraction*exp(-wtd_meters/gw_params%EfoldMaxSatFrac)+fice(:)

    ! Maximum infiltration capacity
    tmpa = ssnow%wbliq(:,1) / efpor(:)

    where (satfrac .lt. 0.99_r_2)

       tmpb = (tmpa-satfrac) / (1.0_r_2 - satfrac)

    elsewhere

       tmpb = (tmpa - satfrac) / 0.01_r_2

    end where

    where ( tmpb .lt. 0.0_r_2) tmpb = 0.0_r_2

    tmpa = -2._r_2*soil%clappB(:,1)*soil%smpsat(:,1)/dzmm_mp(:,1)

    qinmax = (1._r_2 + tmpa*(tmpb-1._r_2))*soil%hksat(:,1)

     ! Surface runoff
    where (ssnow%fwtop(:) .gt. qinmax)

       ssnow%rnof1(:) =  satfrac(:) * ssnow%fwtop(:)  + &
                       (1._r_2-satfrac(:)) * (ssnow%fwtop(:)-qinmax)  !in mm/s

    elsewhere

       ssnow%rnof1(:) = satfrac(:) * ssnow%fwtop(:)  !in mm/s

    end where

   ssnow%fwtop(:) = ssnow%fwtop(:) - ssnow%rnof1(:)

   !ssnow%pudsto(:) = ssnow%pudsto(:) + ssnow%rnof1(:)*dels
   !ssnow%rnof1(:) = 0._r_2

   !where (ssnow%pudsto(:) .gt. satfrac(:)*ssnow%pudsmx)

   !  ssnow%rnof1 = (ssnow%pudsto(:) - satfrac(:)*ssnow%pudsmx)/dels
   !  ssnow%pudsto(:)  =  satfrac(:)*ssnow%pudsmx

   !end where
           
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

      ssnow%rnof1 = ssnow%rnof1 + rnof5/dels   !include this runoff in suface runoff term
   
   END IF

  END SUBROUTINE ovrlndflx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SUBROUTINE simple_wtd(ssnow, soil, veg, ktau, md_prin)
  IMPLICIT NONE
  TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
  TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
  TYPE (veg_parameter_type), INTENT(IN)     :: veg
  INTEGER, INTENT(IN)                       :: ktau  ! integration step number
  LOGICAL, INTENT(IN)                       :: md_prin  !print info?

  REAL(r_2), DIMENSION(mp)            :: fz, wmean
  REAL(r_2), DIMENSION(mp,ms)         :: stot
  REAL(r_2), DIMENSION(mp)            :: ztot
  INTEGER                             :: k

  wmean(:) = 0._r_2
  fz(:)    = 2._r_2
  ztot(:)  = 0._r_2

  stot(:,:) = (ssnow%wb(:,:)-soil%watr(:,:)) / (soil%watsat(:,:)-soil%watr(:,:))

  do k  = 1, ms
     wmean(:) = wmean(:) + stot(:,k)*soil%zse(k)*1000._r_2
     ztot(:) = ztot(:) + soil%zse(k)*1000._r_2
  enddo
  wmean(:) = wmean(:) + ssnow%GWwb(:)/soil%GWwatsat(:) * soil%GWdz(:)*1000._r_2
  ztot(:)     = ztot(:) + soil%GWdz(:)*1000._r_2

  ssnow%wtd(:) = fz(:) * (ztot(:) - wmean(:))

  END SUBROUTINE simple_wtd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


  !----------------------------------------------------------------------
  ! SUBROUTINE calcwtd
  !
  ! Iteratively calcs the water table depth by equating the mass of water in the
  ! soil column to the mass of a hydrostatic column inegrated from the surface to the 
  ! water table depth
  !  
  SUBROUTINE liscalcwtd (ssnow, soil, veg, ktau, md_prin)
  IMPLICIT NONE
  TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
  TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
  TYPE (veg_parameter_type), INTENT(IN)     :: veg
  INTEGER, INTENT(IN)                       :: ktau  ! integration step number
  LOGICAL, INTENT(IN)                       :: md_prin  !print info?

 
  !Local vars 
  REAL(r_2), DIMENSION(mp,ms)   :: dzmm_mp,tmp_def
  REAL(r_2), DIMENSION(0:ms)    :: zimm
  REAL(r_2), DIMENSION(ms)      :: zmm
  REAL(r_2), DIMENSION(mp)      :: GWzimm,temp
  REAL(r_2), DIMENSION(mp)      :: def,defc     

  REAL(r_2)                     :: deffunc,tempa,tempb,derv,calc,tmpc
  REAL(r_2), DIMENSION(mp)      :: invB,Nsmpsat  !inverse of C&H B,Nsmpsat
  INTEGER :: k,i,wttd,jlp
  LOGICAL :: empwtd
   

  empwtd = .false.

  !make code cleaner define these here 
  invB     = 1._r_2/soil%clappB(:,ms)                                !1 over C&H B
  Nsmpsat  = soil%smpsat(:,ms)                                !psi_saturated mm
  dzmm_mp  = real(spread((soil%zse(:)) * 1000.0,1,mp),r_2)    !layer thickness mm
  zimm(0)  = 0.0_r_2                                          !depth of layer interfaces mm

  do k=1,ms

    zimm(k) = zimm(k-1) + soil%zse(k)*1000._r_2

  end do

  zimm(ms) = zimm(ms) + soil%GWdz(1)*1000._r_2
  
  !find the deficit if the water table is at the bottom of the soil column
  defc(:) = (soil%watsat(:,ms))*(zimm(ms)+Nsmpsat(:)/(1._r_2-invB(:))*            &
             (1._r_2-((Nsmpsat(:)+zimm(ms))/Nsmpsat(:))**(1._r_2-invB(:)))) 

  
  where (defc(:) .le. 0._r_2) defc(:) = 0.1_r_2

  where (soil%watsat .gt. ssnow%wb)

    tmp_def = (soil%watsat(:,:)-(ssnow%wbliq+dri*ssnow%wbice))  !prevent freezing from changing wtd

  elsewhere

    tmp_def = 0._r_2

  end where

  def(:) = sum(tmp_def*dzmm_mp,2)
  def(:) = def(:) + max(soil%GWwatsat(:) - ssnow%GWwb(:),0._r_2)*soil%GWdz*1000._r_2

  if (empwtd) then

     ssnow%wtd(:) = zimm(ms)*def(:)/defc(:)

  else

     if (md_prin) write(*,*) 'start wtd iterations'

     ssnow%wtd(:) = zimm(ms)*def(:)/defc(:)
     
     do i=1,mp

      if ((veg%iveg(i) .lt. 16) .and. (soil%isoilm(i) .ne. 9)) then      

       if (defc(i) > def(i)) then                 !iterate tfor wtd

         jlp=0

         mainloop: DO

           tempa   = 1.0_r_2
           tempb   = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(-invB(i))
           derv    = (soil%watsat(i,ms))*(tempa-tempb) + &
                                          soil%watsat(i,ms)

           if (abs(derv) .lt. real(1e-8,r_2)) derv = sign(real(1e-8,r_2),derv)

           tempa   = 1.0_r_2
           tempb   = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(1._r_2-invB(i))
           deffunc = (soil%watsat(i,ms))*(ssnow%wtd(i) +&
                              Nsmpsat(i)/(1-invB(i))* &
                        (tempa-tempb)) - def(i)
           calc    = ssnow%wtd(i) - deffunc/derv

           IF ((abs(calc-ssnow%wtd(i))) .le. wtd_uncert) THEN

             ssnow%wtd(i) = calc
             EXIT mainloop

           ELSEIF (jlp .ge. wtd_iter_max) THEN

              EXIT mainloop

           ELSE

              jlp=jlp+1
              ssnow%wtd(i) = calc

           END IF

         END DO mainloop

       elseif (defc(i) .lt. def(i)) then

         jlp=0

         mainloop2: DO

           tmpc     = Nsmpsat(i)+ssnow%wtd(i)-zimm(ms)
           tempa    = (abs(tmpc/Nsmpsat(i)))**(-invB(i))
           tempb    = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(-invB(i))
           derv     = (soil%watsat(i,ms))*(tempa-tempb)
           if (abs(derv) .lt. real(1e-8,r_2)) derv = sign(real(1e-8,r_2),derv)

           tempa    = (abs((Nsmpsat(i)+ssnow%wtd(i)-zimm(ms))/Nsmpsat(i)))**(1._r_2-invB(i))
           tempb    = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(1._r_2-invB(i))
           deffunc  = (soil%watsat(i,ms))*(zimm(ms) +&
                      Nsmpsat(i)/(1._r_2-invB(i))*(tempa-tempb))-def(i)
           calc     = ssnow%wtd(i) - deffunc/derv

           IF ((abs(calc-ssnow%wtd(i))) .le. wtd_uncert) THEN

             ssnow%wtd(i) = calc
             EXIT mainloop2

           ELSEIF (jlp==wtd_iter_max) THEN

             EXIT mainloop2

           ELSE

             jlp=jlp+1
             ssnow%wtd(i) = calc

           END IF

         END DO mainloop2

       else  !water table depth is exactly on bottom boundary

         ssnow%wtd(i) = zimm(ms)

       endif

      endif  !check veg and soils

     end do   !Loop over all mp tiles

  end if  !debug by using empirical wtd

  !limit wtd to be within a psecified range
  where (ssnow%wtd(:) .gt. wtd_max) ssnow%wtd(:) = wtd_max
  where (ssnow%wtd(:) .lt. wtd_min) ssnow%wtd(:) = wtd_min

  if (md_prin) write(*,*) 'done iterating for wtd'

  END SUBROUTINE liscalcwtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


  !----------------------------------------------------------------------
  ! SUBROUTINE calcwtd
  !
  ! Iteratively calcs the water table depth by equating the mass of water in the
  ! soil column to the mass of a hydrostatic column inegrated from the surface to the 
  ! water table depth
  !  
  SUBROUTINE calcwtd (ssnow, soil, veg, ktau, md_prin)
  IMPLICIT NONE
  TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
  TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
  TYPE (veg_parameter_type), INTENT(IN)     :: veg
  INTEGER, INTENT(IN)                       :: ktau  ! integration step number
  LOGICAL, INTENT(IN)                       :: md_prin  !print info?

 
  !Local vars 
  REAL(r_2), DIMENSION(mp,ms)   :: dzmm_mp,tmp_def
  REAL(r_2), DIMENSION(0:ms)    :: zimm
  REAL(r_2), DIMENSION(ms)      :: zmm
  REAL(r_2), DIMENSION(mp)      :: GWzimm,temp
  REAL(r_2), DIMENSION(mp)      :: def,defc     

  REAL(r_2)                     :: deffunc,tempa,tempb,derv,calc,tmpc,delzb,mx_wtd
  REAL(r_2), DIMENSION(mp)      :: invB,Nsmpsat  !inverse of C&H B,Nsmpsat
  INTEGER :: k,i,wttd,jlp
  LOGICAL :: empwtd
   

  empwtd = .false.

  !make code cleaner define these here 
  invB     = 1._r_2/soil%clappB(:,ms)                         !1 over C&H B
  Nsmpsat  = soil%smpsat(:,ms)                                !psi_saturated mm
  dzmm_mp  = real(spread((soil%zse(:)) * 1000.0,1,mp),r_2)    !layer thickness mm
  zimm(0)  = 0.0_r_2                                          !depth of layer interfaces mm

  do k=1,ms

    zimm(k) = zimm(k-1) + soil%zse(k)*1000._r_2

  end do

  zimm(ms) = zimm(ms) + soil%GWdz(1)*1000._r_2
  
  !find the deficit if the water table is at the bottom of the soil column
  defc(:) = (soil%watsat(:,ms))*(zimm(ms)+Nsmpsat(:)/(1._r_2-invB(:))*            &
             (1._r_2-((Nsmpsat(:)+zimm(ms))/Nsmpsat(:))**(1._r_2-invB(:)))) 
  
  where (defc(:) .le. 0._r_2) defc(:) = 0.1_r_2

  where (soil%watsat .gt. ssnow%wb)

    tmp_def = (soil%watsat(:,:)-(ssnow%wbliq+dri*ssnow%wbice))  !prevent freezing from changing wtd

  elsewhere

    tmp_def = 0._r_2

  end where

  def(:) = sum(tmp_def*dzmm_mp,2)
  def(:) = def(:) + max(soil%GWwatsat(:) - ssnow%GWwb(:),0._r_2)*soil%GWdz*1000._r_2


  if (empwtd) then

    ssnow%wtd(:) = zimm(ms)*def(:)/defc(:)

  else

    if (md_prin) write(*,*) 'start wtd iterations'

    ssnow%wtd(:) = zimm(ms)*def(:)/defc(:)
     
    do i=1,mp

       !if (def(i) .lt. defc(i)) then                 !iterate tfor wtd


      if ((veg%iveg(i) .lt. 16) .and. (soil%isoilm(i) .ne. 9)) then


      jlp=0

      wtd_loop: DO while (jlp .le. wtd_iter_max)

        if (ssnow%wtd(i) .lt. zimm(ms)) then

          tempa   = 1.0_r_2
          tempb   = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(-invB(i))

          derv    = (soil%watsat(i,ms))*(tempa-tempb) + &
                                     soil%watsat(i,ms)

          if (abs(derv) .lt. real(1e-6,r_2)) derv = sign(real(1e-6,r_2),derv)

          tempa   = 1.0_r_2
          tempb   = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(1._r_2-invB(i))
          deffunc = (soil%watsat(i,ms))*(ssnow%wtd(i) +&
                            Nsmpsat(i)/(1-invB(i))* &
                        (tempa-tempb)) - def(i)

          calc    = ssnow%wtd(i) - deffunc/derv

          IF ((abs(calc-ssnow%wtd(i))) .le. wtd_uncert) THEN
            ssnow%wtd(i) = calc
            jlp = wtd_iter_max + 1
          ELSE
            ssnow%wtd(i) = calc
            jlp=jlp+1
          END IF

        elseif (ssnow%wtd(i) .gt. zimm(ms)) then

          tmpc  = Nsmpsat(i)+ssnow%wtd(i)-zimm(ms)
          tempa = (abs(tmpc/Nsmpsat(i)))**(-invB(i))
          tempb = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(-invB(i))
          derv  = (soil%watsat(i,ms))*(tempa-tempb)

          if (abs(derv) .lt. real(1e-8,r_2)) derv = sign(real(1e-8,r_2),derv)

          tempa = (abs((Nsmpsat(i)+ssnow%wtd(i)-zimm(ms))/Nsmpsat(i)))**(1._r_2-invB(i))  
          tempb = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(1._r_2-invB(i))
          deffunc = (soil%watsat(i,ms))*(zimm(ms) +&                     
                      Nsmpsat(i)/(1._r_2-invB(i))*(tempa-tempb))-def(i)

          calc = ssnow%wtd(i) - deffunc/derv

          IF ((abs(calc-ssnow%wtd(i))) .le. wtd_uncert) THEN
            ssnow%wtd(i) = calc
            jlp = wtd_iter_max + 1
          ELSE
            ssnow%wtd(i) = calc
            jlp=jlp+1
          END IF

        else  !water table depth is exactly on bottom boundary

          ssnow%wtd(i) = zimm(ms)

        endif

      end do wtd_loop

      end if  !check if isoil is 9 or iveg is .ge. 16  

    end do   !Loop over all mp tiles

  end if  !debug by using empirical wtd

  !limit wtd to be within a psecified range
  where (ssnow%wtd(:) .gt. wtd_max) ssnow%wtd(:) = wtd_max
  where (ssnow%wtd(:) .lt. wtd_min) ssnow%wtd(:) = wtd_min

  if (md_prin) write(*,*) 'done iterating for wtd'

  END SUBROUTINE calcwtd

  !-------------------------------------------------------------------------
  ! SUBROUTINE smoistgw (fwtop,dt,ktau,ssnow,soil,prin)
  ! solves the modified richards equation (Zeng and Decker 2009) to find
  ! vertical mocement of soil water.  Bottom boundary condition is determined
  ! using a single layer groundwater module
  !
  SUBROUTINE smoistgw (dels,ktau,ssnow,soil,md_prin)
  USE cable_common_module

  IMPLICIT NONE
  
    REAL, INTENT(IN)                          :: dels  ! time step size (s)
    INTEGER, INTENT(IN)                       :: ktau  ! integration step number
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
    LOGICAL, INTENT(IN)                       :: md_prin
    
    !Local variables.  
    REAL(r_2), DIMENSION(mp,ms+1)       :: at     ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: bt     ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: ct     ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: rt

    INTEGER                             :: k,kk,i
    REAL(r_2), DIMENSION(mp,ms)         :: eff_por,old_wb,mss_por  !effective porosity (mm3/mm3),wb(mm3/mm3),mass (mm) of eff_por
    REAL(r_2), DIMENSION(mp,ms)         :: msliq,msice             !mass of the soil liquid and ice water    
    REAL(r_2), DIMENSION(mp,ms)         :: wbrat                   !ratio of volumetric water - watr / saturated - watr
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
    REAL(r_2), DIMENSION(ms)            :: dzmm
    REAL(r_2), DIMENSION(mp,ms)         :: dzmm_mp
    REAL(r_2), DIMENSION(0:ms+1)        :: zimm
    REAL(r_2), DIMENSION(ms)            :: zmm
    REAL(r_2), DIMENSION(mp)            :: GWzimm,xs,zaq,fice_avg,s_mid,GWdzmm
    REAL(r_2), DIMENSION(mp)            :: xsi,xs1,GWmsliq          !mass (mm) of liquid over/under saturation, mass of aquifer water
    REAL(r_2), DIMENSION(mp,ms+1)       :: qhlev,del_wb
    INTEGER, DIMENSION(mp)              :: idlev
    REAL(r_2),DIMENSION(mp,ms+1)           :: masswatmin
    logical                             :: prinall = .false.   !another debug flag
    character (len=30)                  :: fmt  !format to output some debug info
    !MD DEBUG VARS
    INTEGER :: imp,ims

    fmt='(A6,6(1X,F8.6))'       !not needed.  was used to nicely format debug output
    !make code cleaner define these here
    dzmm    = 1000.0_r_2 * real(soil%zse(:),r_2)
    dzmm_mp = spread(dzmm,1,mp)
    zimm(0) = 0._r_2
    zmm(1)  = 0.5_r_2*dzmm(1)

    do k=1,ms

       zimm(k) = zimm(k-1) + dzmm(k)

    end do

    do k=2,ms

       zmm(k)  = zimm(k) + 0.5_r_2*dzmm(k)

    end do 

    GWdzmm(:) = real(soil%GWdz(:),r_2)*1000._r_2
    GWzimm(:) = zimm(ms)+GWdzmm(:)
    zaq(:)    = zimm(ms) + 0.5_r_2*GWdzmm(:)
    
    masswatmin(:,1:ms) = volwatmin * dzmm_mp        !soil must retain this much liquid (mm)
    masswatmin(:,ms+1) = volwatmin * real(C%denliq,r_2) * soil%GWdz

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! preset to allow for non-land & snow points in trimb
    old_wb(:,:) = ssnow%wb(:,:)
    
    !equilibrium water content
    CALL calc_equilibrium_water_content(ssnow,soil)

    !if (ktau .le. 1) then
    !   ssnow%wbliq = ssnow%wbeq
    !   ssnow%wb    = ssnow%wbliq
    !   ssnow%GWwb  = ssnow%GWwbeq
    !end if

    !soil matric potential, hydraulic conductivity, and derivatives of each with respect to water (calculated using total (not liquid))
    do k=1,ms

       ssnow%icefrac(:,k) = ssnow%wbice(:,k)/(max(ssnow%wb(:,k),0.01_r_2))
       ssnow%fracice(:,k) = ssnow%icefrac(:,k)- exp(-3._r_2)! / (1._r_2 - exp(-3._r_2))  !normalizing allows range of 0-1
       ssnow%fracice(:,k) = max(min(ssnow%fracice(:,k),1.0_r_2),0.0_r_2)

    end do

    fice_avg(:)  = sum(ssnow%fracice*dzmm_mp,2) / sum(dzmm)
    fice_avg(:)  = min(max(fice_avg(:),0.0_r_2),1._r_2)  

    where(fice_avg(:) < ssnow%fracice(:,ms)) fice_avg(:) = ssnow%fracice(:,ms)         !frozen ms limits qh

       
    do k=1,ms   
       kk=min(k+1,ms)

       if (k .lt. ms) then

          s1 = 0.5_r_2*(ssnow%wb(:,k)-soil%watr(:,k) + ssnow%wb(:,kk)-soil%watr(:,kk)) / &
               (0.5_r_2*(soil%watsat(:,k)-soil%watr(:,k) + soil%watsat(:,kk)-soil%watr(:,kk)))

          where (s1 .gt. 1._r_2)   s1 = 1._r_2
          where (s1 .lt. 0.01_r_2) s1 = 0.01_r_2

          s2 = soil%hksat(:,k)*s1**(2._r_2*soil%clappB(:,k)+2._r_2)
          ssnow%hk(:,k)    = (1._r_2-0.5_r_2*(ssnow%fracice(:,k)+ssnow%fracice(:,kk)))*s1*s2*&
                              exp(-hkrz*zimm(k)/1000.0_r_2)
          ssnow%dhkdw(:,k) = (1._r_2-0.5_r_2*(ssnow%fracice(:,k)+ssnow%fracice(:,kk)))* &
                    (2._r_2*soil%clappB(:,k)+3._r_2)*s2*0.5_r_2/(soil%watsat(:,k)-soil%watr(:,k))*&
                     exp(-hkrz*zimm(k)/1000.0_r_2)

       else

          s1 = 0.5_r_2*(ssnow%wb(:,k)-soil%watr(:,k) + ssnow%GWwb(:)-soil%GWwatr(:)) / &
               (0.5_r_2*(soil%watsat(:,k)-soil%watr(:,k) + soil%GWwatsat(:)-soil%GWwatr(:)))

          where (s1 .gt. 1._r_2)   s1 = 1._r_2
          where (s1 .lt. 0.01_r_2) s1 = 0.01_r_2

          s2 = soil%hksat(:,k)*s1**(2._r_2*soil%clappB(:,k)+2._r_2)
          ssnow%hk(:,k)    = s1*s2*(1._r_2-ssnow%fracice(:,k))* exp(-hkrz*zimm(ms)/1000._r_2)
          ssnow%dhkdw(:,k) = (1._r_2-ssnow%fracice(:,k))* (2._r_2*soil%clappB(:,k)+3._r_2)*&
                         s2*0.5_r_2/(soil%watsat(:,k)-soil%watr(:,k))*exp(-hkrz*zimm(ms)/1000._r_2)

       end if
  
       s_mid = (ssnow%wb(:,k)-soil%watr(:,k))/&  !+dri*ssnow%wbice(:,k)
              (soil%watsat(:,k)-soil%watr(:,k))

       where (s_mid .gt. 1._r_2)   s_mid = 1._r_2
       where (s_mid .lt. 0.01_r_2) s_mid = 0.01_r_2

       ssnow%smp(:,k) = -soil%smpsat(:,k)*s_mid**(-soil%clappB(:,k))

       where (ssnow%smp(:,k) .gt. -soil%smpsat(:,k)) ssnow%smp(:,k) = -soil%smpsat(:,k)
       where (ssnow%smp(:,k) .lt. sucmin)            ssnow%smp(:,k) = sucmin

       ssnow%dsmpdw(:,k) = -soil%clappB(:,k)*ssnow%smp(:,k)/&
                    (max(s_mid*(soil%watsat(:,k)-soil%watr(:,k)),0.001_r_2))          
    end do


    !Aquifer properties
    s_mid = (ssnow%GWwb(:)-soil%GWwatr(:))/(soil%GWwatsat(:)-soil%GWwatr(:))
    where (s_mid .gt. 1._r_2) s_mid = 1._r_2
    where (s_mid .lt. 0.01_r_2) s_mid = 0.01_r_2
    s2    = soil%GWhksat(:)*s_mid**(2._r_2*soil%clappB(:,ms)+2._r_2)

    ssnow%GWhk(:)     = s_mid*s2*(1._r_2-ssnow%fracice(:,ms))

    ssnow%GWdhkdw(:)  = (1._r_2-ssnow%fracice(:,ms))* (2._r_2*soil%clappB(:,ms)+3._r_2)*&
                         s2*0.5_r_2/(soil%GWwatsat(:)-soil%GWwatr(:))
    ssnow%GWsmp(:)    = -soil%smpsat(:,ms)*s_mid**(-soil%clappB(:,ms))

    where (ssnow%GWsmp .lt. sucmin) &
           ssnow%GWsmp = sucmin

    where (ssnow%GWsmp .gt. -soil%GWsmpsat) &
           ssnow%GWsmp = -soil%GWsmpsat

    ssnow%GWdsmpdw(:) = -soil%GWclappB(:)*ssnow%GWsmp(:)/(s_mid*(soil%GWwatsat(:)-soil%GWwatr(:)))
       
    !Note: temporaary parameteriation of horizontal drainage
    !too be replaced with explivit treatment of subgrid scale, topographically
    !based subsurface flux convergence flowing to river channels 
       
    ssnow%qhz(:)  = gw_params%MaxHorzDrainRate *(1._r_2 - fice_avg(:)) * &
                    exp(-ssnow%wtd(:)/(1000._r_2*gw_params%EfoldHorzDrainRate))
    !find index of soil layer with the water table
    qhlev(:,:)   = 0._r_2  !set to zero except for layer that contains the wtd
    idlev(:)     = ms+1    !assume the water table is in the aquifer
    do k=ms,1,-1

       where (ssnow%wtd(:) <= zimm(k) .and. ssnow%wtd(:) > zimm(k-1))

         idlev(:) = k

       endwhere

    end do

    do i=1,mp
         !If the layer is frozen there shouldn't be any drainage
         ssnow%qhz(i)      = ssnow%qhz(i) * (1._r_2 - ssnow%fracice(i,min(idlev(i),ms)))
         qhlev(i,idlev(i)) = ssnow%qhz(i)

    end do  


    rt(:,:) = 0._r_2; at(:,:) = 0._r_2     !ensure input to tridiag is valid
    bt(:,:) = 1._r_2; ct(:,:) = 0._r_2

    k = 1     !top soil layer
       qin     = ssnow%sinfil(:)
       den     = (zmm(k+1)-zmm(k))
       dne     = (ssnow%zq(:,k+1)-ssnow%zq(:,k))
       num     = (ssnow%smp(:,k+1)-ssnow%smp(:,k)) - dne
       qout    = -ssnow%hk(:,k)*num/den
       dqodw1  = -(-ssnow%hk(:,k)*ssnow%dsmpdw(:,k)   + num*ssnow%dhkdw(:,k))/den
       dqodw2  = -( ssnow%hk(:,k)*ssnow%dsmpdw(:,k+1) + num*ssnow%dhkdw(:,k))/den
       rt(:,k) =  qin - qout !- qhlev(:,k)
       at(:,k) =  0._r_2
       bt(:,k) =  dzmm(k)/dels + dqodw1
       ct(:,k) =  dqodw2      

    do k = 2, ms - 1     !middle soil layers

       den     = (zmm(k) - zmm(k-1))
       dne     = (ssnow%zq(:,k)-ssnow%zq(:,k-1))
       num     = (ssnow%smp(:,k)-ssnow%smp(:,k-1)) - dne
       qin     = -ssnow%hk(:,k-1)*num/den
       dqidw0  = -(-ssnow%hk(:,k-1)*ssnow%dsmpdw(:,k-1) + num*ssnow%dhkdw(:,k-1))/den
       dqidw1  = -( ssnow%hk(:,k-1)*ssnow%dsmpdw(:,k)   + num*ssnow%dhkdw(:,k-1))/den
       den     = (zmm(k+1)-zmm(k))
       dne     = (ssnow%zq(:,k+1)-ssnow%zq(:,k))
       num     = (ssnow%smp(:,k+1)-ssnow%smp(:,k)) - dne
       qout    = -ssnow%hk(:,k)*num/den
       dqodw1  = -(-ssnow%hk(:,k)*ssnow%dsmpdw(:,k)   + num*ssnow%dhkdw(:,k))/den
       dqodw2  = -( ssnow%hk(:,k)*ssnow%dsmpdw(:,k+1) + num*ssnow%dhkdw(:,k))/den
       rt(:,k) =  qin - qout !- qhlev(:,k)
       at(:,k) = -dqidw0
       bt(:,k) =  dzmm(k)/dels - dqidw1 + dqodw1
       ct(:,k) =  dqodw2
       
    end do
       
    k = ms   !Bottom soil layer
       den     = (zmm(k) - zmm(k-1))
       dne     = (ssnow%zq(:,k)-ssnow%zq(:,k-1))
       num     = (ssnow%smp(:,k)-ssnow%smp(:,k-1)) - dne
       qin     = -ssnow%hk(:,k-1)*num/den
       dqidw0  = -(-ssnow%hk(:,k-1)*ssnow%dsmpdw(:,k-1) + num*ssnow%dhkdw(:,k-1))/den
       dqidw1  = -( ssnow%hk(:,k-1)*ssnow%dsmpdw(:,k)   + num*ssnow%dhkdw(:,k-1))/den
       den     = zaq(:) - zmm(k)
       dne     = (ssnow%GWzq(:)-ssnow%zq(:,k))
       num     =  (ssnow%GWsmp(:)-ssnow%smp(:,k)) - dne
       qout    = -ssnow%hk(:,k)*num/den
       dqodw1  = -(-ssnow%hk(:,k)*ssnow%dsmpdw(:,k)   + num*ssnow%dhkdw(:,k))/den
       dqodw2  = -( ssnow%hk(:,k)*ssnow%GWdsmpdw(:) + num*ssnow%dhkdw(:,k))/den
       rt(:,k) =  qin - qout !- qhlev(:,k)
       at(:,k) = -dqidw0
       bt(:,k) =  dzmm(k)/dels - dqidw1 + dqodw1
       ct(:,k) =  dqodw2 

       !den = dzmm(k)
       !dne = -ssnow%zq(:,k)
       !num = -ssnow%smp(:,k) - dne
       !qout = -ssnow%hk(:,k)*num/den
       !dqodw1 = -(-ssnow%hk(:,k)*ssnow%dsmpdw(:,k) + num*ssnow%dhkdw(:,k))/den
       !dqodw2 = 0._r_2
       !rt(:,k) = qin-qout
       !at(:,k) = -dqidw0
       !bt(:,k) = dzmm(k)/dels - dqidw1 + dqodw1
       !ct(:,k) = dqodw2
         
    k = ms+1   !Aquifer layer
       den     = (zaq(:) - zmm(k-1))
       dne     = (ssnow%GWzq(:)-ssnow%zq(:,k-1))
       num     = (ssnow%GWsmp(:)-ssnow%smp(:,k-1)) - dne
       qin     = -ssnow%hk(:,k-1)*num/den
       dqidw0  = -(-ssnow%hk(:,k-1)*ssnow%dsmpdw(:,k-1) + num*ssnow%dhkdw(:,k-1))/den
       dqidw1  = -( ssnow%hk(:,k-1)*ssnow%GWdsmpdw(:)   + num*ssnow%dhkdw(:,k-1))/den
       den     = zaq(:) - zmm(k-1)
       dne     = (ssnow%GWzq(:)-ssnow%zq(:,k-1))
       num     = (ssnow%GWsmp(:)-ssnow%smp(:,k-1)) - dne
       qout    = 0._r_2
       dqodw1  = 0._r_2
       dqodw2  = 0._r_2
       rt(:,k) =  qin - qout ! - qhlev(:,k)
       at(:,k) = -dqidw0
       bt(:,k) =  GWdzmm(:)/dels - dqidw1
       ct(:,k) =  0._r_2

    !CALL solve_tridiag(at, bt, ct, rt, del_wb,ms+1)                      !solve system of eqns
    CALL trimb(at,bt,ct,rt,ms+1)                       !use the defulat cable tridiag solution

    ssnow%wbliq = ssnow%wbliq + rt(:,1:ms) - qhlev(:,1:ms)*dels/spread(dzmm,1,mp)   !volutermic liquid
    ssnow%GWwb  = ssnow%GWwb  + rt(:,ms+1) - qhlev(:,ms+1)*dels/GWdzmm
    !ssnow%GWwb = ssnow%GWwb + qout + dqodw1*rt(:,ms) - qhlev(:,ms+1)*dels/GWdzmm  !use if solve 1:ms but not if solve 1:ms+1

    msliq       = ssnow%wbliq * dzmm_mp                    !liquid mass
    msice       = ssnow%wbice * dzmm_mp * dri              !ice mass
    GWmsliq     = ssnow%GWwb  * GWdzmm

    !determine the available pore space
    !volumetric
    eff_por     = soil%watsat - ssnow%wbice
    !mass of liquid that can fit in eff_por
    mss_por     = eff_por * dzmm_mp
    
    xsi(:) = 0._r_2

    where (GWmsliq .gt. soil%GWwatsat*GWdzmm)  !if GW oversat put into bottom soil layer

      xsi(:)  = GWmsliq - soil%GWwatsat*GWdzmm
      GWmsliq =  soil%GWwatsat*GWdzmm

    end where

    do k=1,ms

       where(msliq(:,k) .gt. mss_por(:,k))

          xsi(:) = xsi(:) + msliq(:,k) - mss_por(:,k)
          msliq(:,k) = mss_por(:,k)

       end where

    end do
     
    do k = ms,1,-1  !loop from bottom to top adding extra water to each layer

       where ((xsi(:) .lt. (mss_por(:,k)-msliq(:,k))) .and. (xsi(:) .gt. 0._r_2))

          msliq(:,k) = msliq(:,k) + xsi(:)
          xsi(:)     = 0._r_2

       elsewhere ((xsi(:) .gt. 0.0_r_2) .and. (xsi(:) .gt. (mss_por(:,k)-msliq(:,k))))

          xsi(:)     = xsi(:) - (mss_por(:,k)-msliq(:,k))
          msliq(:,k) = mss_por(:,k)

       end where

    end do

    where (xsi(:) .gt. 0._r_2) &
          ssnow%qhz(:) = ssnow%qhz(:) + xsi(:)/dels

    xsi(:) = 0._r_2
    
   
    do k = 1,ms

       do i=1,mp  !ensure liq > liq_minimum (using mm)

          xs(i) = 0._r_2             !should be a single float (array not needed)

          if (msliq(i,k) .lt. masswatmin(i,k)) then

             xs(i) = masswatmin(i,k) - msliq(i,k)
             msliq(i,k) = masswatmin(i,k)

             if (k .lt. ms) then
                msliq(i,k+1) = msliq(i,k+1) - xs(i)
             else
                GWmsliq(i) = GWmsliq(i) - xs(i)
                !ssnow%qhz(i)  = ssnow%qhz(i) - xs(i)/dels  !take from runnof may make runoff < 0
             end if

          end if


       end do

    end do
 
     xs(:) = 0._r_2
     do i=1,mp

        if (GWmsliq(i) .lt. masswatmin(i,ms+1)) then

           xs(i)        = masswatmin(i,ms+1) - GWmsliq(i)
           GWmsliq(i)   = masswatmin(i,ms+1)
           ssnow%qhz(i) = ssnow%qhz(i) - xs(i)/dels

        end if
     end do


    !update all variabes
    ssnow%wmliq(:,:) = msliq(:,:)
    ssnow%wmice(:,:) = msice(:,:)
    ssnow%wbliq(:,:) = msliq(:,:) /  (dzmm_mp)     !convert from mm to volumetric
    ssnow%wbice(:,:) = msice(:,:) /  (dzmm_mp * dri)     !convert from mm to volumetric
    ssnow%GWwb(:)    = GWmsliq(:) / GWdzmm(:)  
    ssnow%wb         = ssnow%wbliq(:,:) + ssnow%wbice(:,:)
    ssnow%wmtot      = ssnow%wmliq(:,:) + ssnow%wmice(:,:)
    ssnow%rnof2(:)   = ssnow%qhz(:)                                       !rnof2 is default cable deep runoff var 
       

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

   INTEGER             :: k,i
   REAL, DIMENSION(mp) :: snowmlt
   REAL, DIMENSION(mp) :: totwet
   REAL, DIMENSION(mp) :: weting,GWwb_ic,wberr
   REAL, DIMENSION(mp) :: xxx, tgg_old, tggsn_old,wbtot_ic,del_wbtot
   REAL(r_2), DIMENSION(mp) :: xx,deltat,sinfil1,sinfil2,sinfil3 
   REAL                :: zsetot
   INTEGER, SAVE :: ktau =0 
   LOGICAL :: prin,md_prin
   REAL(r_2) :: wb_lake_T, rnof2_T, ratio

   
   prin = .FALSE.
   md_prin = .false.
   
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


   IF( cable_runtime%offline .or. cable_runtime%mk3l ) ssnow%t_snwlr = 0.05_r_2

   ssnow%fwtop1 = 0.0
   ssnow%fwtop2 = 0.0
   ssnow%fwtop3 = 0.0
   ssnow%runoff = 0.0 ! initialise total runoff
   ssnow%rnof1  = 0.0 ! initialise surface runoff
   ssnow%rnof2  = 0.0 ! initialise deep drainage
   ssnow%smelt  = 0.0 ! initialise snowmelt
   ssnow%dtmlt  = 0.0 
   ssnow%osnowd = ssnow%snowd
   ! Scaling  runoff to kg/m^2/s (mm/s) to match rest of the model
   ssnow%sinfil = 0.0   

   IF( .NOT.cable_user%cable_runtime_coupled ) THEN
   
      IF( ktau_gl <= 1 ) THEN
         
         IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
                                                  ! after discussion with BP
         ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
         ssnow%wbtot = 0.0
         ssnow%wb(:,:)  = MIN( soil%watsat(:,:), MAX ( ssnow%wb(:,:), spread(soil%swilt(:),2,ms) ) )   

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

         END WHERE

         WHERE (spread(soil%isoilm,2,ms) == 9)

              ssnow%wb    = 0.95 * soil%watsat
              ssnow%wbice = 0.90 * ssnow%wb

         END WHERE
         
         xx=soil%css * soil%rhosoil

         ssnow%gammzz(:,1) = MAX( (1.0 - soil%watsat(:,1)) * soil%css * soil%rhosoil &
              & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * cswat * C%denliq &
              & + ssnow%wbice(:,1) * csice * C%denice, xx ) * soil%zse(1)

      END IF

   ENDIF  ! if(.NOT.cable_runtime_coupled)

   !Start with wb and wbice.  Need wbliq, wmliq,wmice,wmtot
   !find the mass of ice and liq from the prognostic volumetric values
   ssnow%wbliq = ssnow%wb - ssnow%wbice                     !liquid volume
   ssnow%wmice = ssnow%wbice*C%denice*spread(soil%zse,1,mp) !ice mass
   ssnow%wmliq = ssnow%wbliq*C%denliq*spread(soil%zse,1,mp) !liquid mass
   ssnow%wmtot = ssnow%wmice + ssnow%wmliq                  !liq+ice mass

   xx=soil%css * soil%rhosoil

   IF (ktau <= 1)                                                              &
     ssnow%gammzz(:,1) = MAX( (1.0 - soil%watsat(:,1)) * soil%css * soil%rhosoil      &
            & + ssnow%wbliq(:,1) * cswat * C%denliq           &
            & + ssnow%wbice(:,1) * csice * C%denice, xx ) * soil%zse(1) +   &
            & (1. - ssnow%isflag) * cgsnow * ssnow%snowd

   ssnow%wblf   = max(0.01_r_2,ssnow%wbliq/soil%watsat)
   ssnow%wbfice = max(0.01_r_2,ssnow%wbice/soil%watsat)  

   !initial water in the soil column
   wbtot_ic  = sum(ssnow%wbliq(:,:)*C%denliq*spread(soil%zse,1,mp),2) + &
               sum(ssnow%wbice(:,:)*C%denice*spread(soil%zse,1,mp),2) + &
               ssnow%GWwb(:)*soil%GWdz*C%denliq
               
   GWwb_ic = ssnow%GWwb

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

   !do the soil and snow melting, freezing prior to water movement
   ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

   CALL snow_melting (dels, snowmlt, ssnow, soil )
   
   ! Add new snow melt to global snow melt variable: 
   ssnow%smelt = ssnow%smelt + snowmlt

   CALL remove_trans(dels, soil, ssnow, canopy, veg)        !transpiration loss per soil layer

   CALL  soilfreeze(dels, soil, ssnow)

   ssnow%fwtop = canopy%precis/dels + ssnow%smelt/dels   !water from canopy and snowmelt [mm/s]   
   !ssnow%rnof1 = ssnow%rnof1 + ssnow%smelt / dels          !adding snow melt directly to the runoff

   !CALL calcwtd (ssnow, soil, veg, ktau, md_prin)                  !update the wtd
   !CALL liscalcwtd (ssnow, soil, veg, ktau, md_prin)            !test the calcwtd from lis to bug hunt
   CALL simple_wtd (ssnow, soil, veg, ktau, md_prin)

   CALL ovrlndflx (dels, ktau, ssnow, soil, md_prin )         !surface runoff, incorporate ssnow%pudsto?
   
   ssnow%sinfil = ssnow%fwtop - canopy%segg  !canopy%fes/C%HL               !remove soil evap from throughfall
   !ssnow%pudsto = max(ssnow%pudsto - canopy%fesp/C%HL*dels,0._r_2)  !currently pudsto = 0.0 always

   CALL smoistgw (dels,ktau,ssnow,soil,md_prin)               !vertical soil moisture movement. 
  
   ! lakes: replace hard-wired vegetation number in next version
   WHERE( veg%iveg == 16 )

      ssnow%sinfil = MIN( ssnow%rnof1, ssnow%wb_lake&
                   + MAX( 0.,canopy%segg ) )  !segg not in other version.  segg is mm/s?????
      ssnow%rnof1 = MAX( 0.0, ssnow%rnof1 - ssnow%sinfil )
      ssnow%wb_lake = ssnow%wb_lake - ssnow%sinfil
      ssnow%rnof2 = MAX( 0.0, ssnow%rnof2 - ssnow%wb_lake )

   ENDWHERE    

   ! correction required for energy balance in online simulations 
   IF( cable_runtime%um) THEN

      canopy%fhs_cor = ssnow%dtmlt(:,1)*ssnow%dfh_dtg
      canopy%fes_cor = ssnow%dtmlt(:,1)*(ssnow%cls*ssnow%dfe_ddq * ssnow%ddq_dtg)

      canopy%fhs = canopy%fhs+canopy%fhs_cor
      canopy%fes = canopy%fes+canopy%fes_cor

   ENDIF

   ssnow%smelt  = ssnow%smelt/dels    !change units to mm/s.  cable_driver then reverts back to mm
   ssnow%runoff = (ssnow%rnof1 + ssnow%rnof2)!*dels  !cable_driver converts from mm/s to mm
                                                     !rnof1 and rnof2 are already in mm/s

   ! Set weighted soil/snow surface temperature
   ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

   !total water mass at the end of the soilsnow_GW routine
   ssnow%wbtot = sum(ssnow%wbliq(:,:)*real(C%denliq,r_2)*real(spread(soil%zse,1,mp),r_2),2) + &
                 sum(ssnow%wbice(:,:)*real(C%denice,r_2)*real(spread(soil%zse,1,mp),r_2),2) + &
                 ssnow%GWwb(:)*soil%GWdz*real(C%denliq,r_2)
                 
   !for debug water balance.  del_wbtot = fluxes = infiltration [though-evap] - trans - qhorz drainage
   del_wbtot   = dels * (ssnow%sinfil - ssnow%rnof2 - canopy%fevc / C%HL)
   !set below to keep track of water imbalance within the GW module explicitly.  also must change cable_checks
   ssnow%wbtot = ssnow%wbtot-wbtot_ic

END SUBROUTINE soil_snow_gw



SUBROUTINE calc_srf_wet_fraction(ssnow,soil)

  IMPLICIT NONE
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters

    !local variables
    REAL(r_2), DIMENSION(mp)           :: xxx,fice,icef,efpor,xx
    REAL(r_2), DIMENSION(mp)           :: satfrac,wtd_meters
    REAL(r_2), DIMENSION(mp,ms)        :: liqmass,icemass,totmass
    REAL(r_2), DIMENSION(mp,ms)        :: dzmm_mp

    xxx(:)   = 0._r_2
    dzmm_mp  = 1000._r_2 * real(spread(soil%zse,1,mp),r_2)

    icemass  = ssnow%wbice(:,:) * dzmm_mp * dri
    liqmass  = (ssnow%wb-ssnow%wbice) * dzmm_mp
    totmass  = icemass + liqmass

    where (totmass .lt. real(1e-2,r_2)) totmass = real(1e-2,r_2)

    efpor(:) = soil%watsat(:,1) - ssnow%wbice(:,1)!-soil%watr(:,1)
    where (efpor .lt. 0.05_r_2) efpor = 0.05_r_2

    !srf frozen fraction.  should be based on topography
    icef(:) = max(0._r_2,min(1._r_2,icemass(:,1) / totmass(:,1)))
    fice(:) = (exp(-3._r_2*(1._r_2-icef(:)))- exp(-3._r_2))/(1._r_2-exp(-3._r_2))
    where (fice(:) .lt. 0._r_2) fice(:) = 0._r_2
    where (fice(:) .gt. 1._r_2) fice(:) = 1._r_2

    ! Saturated fraction
    wtd_meters = ssnow%wtd / 1000._r_2


    xx(:) = max(real(1e-6,r_2), min(1._r_2, (ssnow%wbliq(:,1)-0.5_r_2*real(soil%swilt(:),r_2))/&
                                    (real(soil%sfc(:)-0.5*soil%swilt(:),r_2))))

    !xxx(:) = xx(:)! (exp(-2._r_2*xx(:))-1._r_2) / (exp(-2._r_2)-1._r_2)

    satfrac(:) = gw_params%MaxSatFraction*exp(-wtd_meters/gw_params%EfoldMaxSatFrac) &
                 +(1._r_2 - gw_params%MaxSatFraction*exp(-wtd_meters/gw_params%EfoldMaxSatFrac))* xx(:)

    ssnow%wetfac(:) = fice(:) + ( 1._r_2 - fice(:) )*satfrac(:)



END SUBROUTINE calc_srf_wet_fraction

SUBROUTINE calc_equilibrium_water_content(ssnow,soil)

  IMPLICIT NONE
  
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
    !local variables
    REAL(r_2), dimension(mp)    :: zaq      !node depth of the aquifer
    REAL(r_2), dimension(mp,ms) :: dzmm_mp  !layer thicknesses at each land tile
    REAL(r_2), dimension(ms)    :: dzmm     !layer thickness for single tile
    REAL(r_2), dimension(mp)    :: GWdzmm   !aquifer thickness at each tile
    REAL(r_2), dimension(mp)    :: GWzimm   !aquifer layer interface depth
    REAL(r_2), dimension(0:ms)  :: zimm     !layer interface depth in mm  
    REAL(r_2), dimension(ms)    :: zmm      !node depths in mm

    REAL(r_2), dimension(mp,ms) :: wbrat    !temporary variable.  eq water content prior to checks

    REAL(r_2), dimension(mp)    :: voleq1   !temp variable
    REAL(r_2), dimension(mp)    :: tempi, temp0

    INTEGER :: k

    !make code cleaner define these here
    dzmm    = 1000.0_r_2 * real(soil%zse(:),r_2)
    dzmm_mp = spread(dzmm,1,mp)
    zimm(0) = 0._r_2
    zmm(1)  = 0.5_r_2*dzmm(1)

    do k=1,ms
       zimm(k) = zimm(k-1) + dzmm(k)
    end do

    do k=2,ms
       zmm(k)  = zimm(k) + 0.5_r_2*dzmm(k)
    end do 

    GWdzmm(:) = real(soil%GWdz(:),r_2)*1000._r_2
    GWzimm(:) = zimm(ms)+GWdzmm(:)
    zaq(:)    = zimm(ms) + 0.5_r_2*GWdzmm(:)
    
    !equilibrium water content
    do k=1,ms

       WHERE ((ssnow%wtd(:) .le. zimm(k-1)))          !fully saturated

          ssnow%wbeq(:,k) = real(soil%watsat(:,k)-soil%watr(:,k),r_2)

       END WHERE

       WHERE ((ssnow%wtd(:) .le. zimm(k)) .and. (ssnow%wtd(:) .gt. zimm(k-1)))

          tempi = 1._r_2
          temp0 = (((soil%smpsat(:,k)+ssnow%wtd(:)-zimm(k-1))/soil%smpsat(:,k)))**(1._r_2-1._r_2/soil%clappB(:,k))               
          voleq1 = -soil%smpsat(:,k)*(soil%watsat(:,k)-soil%watr(:,k))/&
                       (1._r_2-1._r_2/soil%clappB(:,k))/(ssnow%wtd(:)-zimm(k-1))*(tempi-temp0)
          ssnow%wbeq(:,k) = (voleq1*(ssnow%wtd(:)-zimm(k-1)) + (soil%watsat(:,k)-soil%watr(:,k))&
                            *(zimm(k)-ssnow%wtd(:)))/(zimm(k)-zimm(k-1)) + soil%watr(:,k)

       END WHERE

       WHERE (ssnow%wtd .ge. zimm(k))

          tempi = (((soil%smpsat(:,k)+ssnow%wtd(:)-zimm(k))/soil%smpsat(:,k)))**(1._r_2-1._r_2/soil%clappB(:,k))
          temp0 = (((soil%smpsat(:,k)+ssnow%wtd(:)-zimm(k-1))/soil%smpsat(:,k)))**(1._r_2-1._r_2/soil%clappB(:,k))   
          ssnow%wbeq(:,k) = -soil%smpsat(:,k)*(soil%watsat(:,k)-soil%watr(:,k))/&
                               (1._r_2-1._r_2/soil%clappB(:,k))/(zimm(k)-zimm(k-1))*(tempi-temp0)+soil%watr(:,k)

       END WHERE
    end do
 
    where (ssnow%wbeq(:,:) .gt. soil%watsat(:,:)) ssnow%wbeq(:,:) = soil%watsat(:,:)

    where (ssnow%wbeq(:,:) .le. soil%watr(:,:)) ssnow%wbeq(:,:) = soil%watr(:,:)+0.01_r_2 !0.01 is arbitrary

    wbrat(:,:) = (ssnow%wbeq(:,:) - soil%watr(:,:))/(soil%watsat(:,:) - soil%watr(:,:))

    where (wbrat .lt. 0.1_r_2)  &
           wbrat = 0.1_r_2

    where (wbrat .gt. 1._r_2)   &
           wbrat = 1._r_2

    ssnow%zq(:,:) = -soil%smpsat(:,:)*(wbrat(:,:)**(-soil%clappB(:,:)))

    where (ssnow%zq(:,:) .lt. sucmin) ssnow%zq(:,:) = sucmin

    where (ssnow%zq(:,:) .gt. -soil%smpsat) ssnow%zq(:,:) = -soil%smpsat

       
    !Aquifer Equilibrium water content
    WHERE (ssnow%wtd(:) .le. zimm(ms))                                          !fully saturated

       ssnow%GWwbeq(:) = real(soil%GWwatsat(:)-soil%GWwatr(:),r_2)

    END WHERE

    WHERE ((ssnow%wtd(:) .gt. GWzimm(:)))                                       !fully unsaturated

       tempi = (((soil%GWsmpsat(:)+ssnow%wtd(:)-GWzimm(:))/soil%GWsmpsat(:)))**(1._r_2-1._r_2/soil%GWclappB(:))
       temp0 = (((soil%GWsmpsat(:)+ssnow%wtd(:)-zimm(ms))/soil%GWsmpsat(:)))**(1._r_2-1._r_2/soil%GWclappB(:))   
       ssnow%GWwbeq(:) = -soil%GWsmpsat(:)*soil%GWwatsat(:)/&
                          (1._r_2-1._r_2/soil%GWclappB(:))/(GWzimm(:)-zimm(ms))*(tempi-temp0) + soil%GWwatr(:)   

    END WHERE           

    WHERE ((ssnow%wtd(:) .le. GWzimm(:)) .and. (ssnow%wtd(:) .gt. zimm(ms)))    !partially saturated

       tempi  = 1._r_2
       temp0  = (((soil%GWsmpsat(:)+ssnow%wtd(:)-zimm(ms))/soil%GWsmpsat(:)))**(1._r_2-1._r_2/soil%GWclappB(:))               
       voleq1 = -soil%GWsmpsat(:)*(soil%GWwatsat(:)-soil%GWwatr(:))/&
                (1._r_2-1._r_2/soil%GWclappB(:))/(ssnow%wtd(:)-zimm(ms))*(tempi-temp0) + soil%GWwatr(:)
       ssnow%GWwbeq(:) = (voleq1*(ssnow%wtd(:)-zimm(ms)) + (soil%GWwatsat(:)-soil%GWwatr(:))*&
                         (GWzimm(:)-ssnow%wtd(:)))/(GWzimm(:)-zimm(ms)) + soil%GWwatr(:)

    END WHERE

    ssnow%GWwbeq(:) = min(real(soil%GWwatsat(:),r_2),ssnow%GWwbeq(:))
    ssnow%GWwbeq(:) = max(ssnow%GWwbeq(:),soil%GWwatr(:)+0.01_r_2) 

    ssnow%GWzq(:) = -soil%GWsmpsat(:)*(max((ssnow%GWwbeq(:)-soil%GWwatr(:))/           &
                    (soil%GWwatsat(:)-soil%GWwatr(:)),0.01_r_2))**(-soil%GWclappB(:))
    ssnow%GWzq(:) = max(sucmin, ssnow%GWzq(:))

END SUBROUTINE calc_equilibrium_water_content



END MODULE cable_soil_snow_gw_module
