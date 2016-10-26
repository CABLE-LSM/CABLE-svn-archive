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

   USE cable_common_module, ONLY : gw_params,cable_user

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
      max_sconds = 2.51!,   & !
   !   frozen_limit = 0.85    ! EAK Feb2011 (could be 0.95) 
      
   !mrd561 GW params
   !Should read some in from namelist
   REAL(r_2), PARAMETER :: sucmin  = -1.0e8, &! minimum soil pressure head [mm]
                      volwatmin    = 1e-4,        &!min soil water [mm]      
                      wtd_uncert   = 0.1,         &! uncertaintiy in wtd calcultations [mm]
                      wtd_max      = 1000000.0,    &! maximum wtd [mm]
                      wtd_min      = 100.0,       &! minimum wtd [mm]
                      dri          = 1.0           !ratio of density of ice to density of liquid [unitless]

   INTEGER, PARAMETER :: wtd_iter_max = 20 ! maximum number of iterations to find the water table depth                    
   
   REAL :: cp    ! specific heat capacity for air
   
   !jhan:make parameter
   REAL :: max_glacier_snowd
 
   ! This module contains the following subroutines:
   PUBLIC soil_snow_gw,calc_srf_wet_fraction,calc_soil_hydraulic_props ! must be available outside this module
   PRIVATE snowdensity, snow_melting, snowcheck, snowl_adjust 
   PRIVATE trimb,snow_accum, stempv,calc_equilibrium_water_content
   PRIVATE GWsoilfreeze, remove_trans,iterative_wtd,simple_wtd
   PRIVATE smoistgw, ovrlndflx

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
    
   TYPE(soil_parameter_type), INTENT(IN)    :: soil

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
      WHERE( soil%isoilm .ne. 9 ) ssnow%ssdn(:,1) = MIN( 450.0, ssnow%ssdn(:,1) )

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
   
   TYPE(soil_parameter_type), INTENT(IN)    :: soil
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
         WHERE( soil%isoilm .eq. 9)                                             &
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
   TYPE(soil_parameter_type), INTENT(IN)    :: soil ! soil parameters
   
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
         WHERE( soil%isoilm .eq. 9 )                                             &
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
         WHERE( soil%isoilm .ne. 9 )                                             &
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
         WHERE( soil%isoilm .ne. 9 )                                             &
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
         WHERE( soil%isoilm .ne. 9 )                                             &
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
  
! calculates temperatures of the soil
! tgg - new soil/snow temperature
! ga - heat flux from the atmosphere (ground heat flux)
! ccnsw - soil thermal conductivity, including water/ice
SUBROUTINE GWstempv(dels, canopy, ssnow, soil)
   REAL, INTENT(IN) :: dels ! integration time step (s)
   
   TYPE(canopy_type),    INTENT(INOUT) :: canopy
   TYPE(soil_snow_type), INTENT(INOUT) :: ssnow
   
   TYPE(soil_parameter_type), INTENT(IN)    :: soil
   
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
   real(r_2) :: tmp,tmpsat
   real(r_2) :: om_tkm, tkm(mp,ms),tkmg(mp,ms),liq_frac
   LOGICAL :: direct2min = .FALSE.

   om_tkm       = 0.25_r_2
   tkm(:,:) = (1._r_2-soil%Forg(:,:))*(8.80_r_2*soil%Fsand(:,:)+2.92_r_2*soil%Fclay(:,:))/(soil%Fsand(:,:)+2.92_r_2*soil%Fclay(:,:))+om_tkm*soil%Forg(:,:)

   tkmg(:,:) = tkm(:,:) ** (1._r_2- soil%watsat(:,:))

   at = 0.0
   bt = 1.0
   ct = 0.0
   coeff = 0.0
   snow_ccnsw = 2.0

   DO k = 1, ms
      
      DO j = 1, mp
      
         IF( soil%isoilm(j) .eq. 9 ) THEN
            ! permanent ice: fix hard-wired number in next version
            ccnsw(j,k) = snow_ccnsw
         ELSE  !test alt formation for thermal conductivity
            ew(j) =min(1.0_r_2, (ssnow%wmliq(j,k)/rhowat + ssnow%wmice(j,k)/(rhowat*dri))/(1000.*soil%zse(k)*soil%watsat(j,k)))
            if (ew(j) .gt. 1.0e-7) then
               liq_frac = ssnow%wmliq(j,k)/(ssnow%wmliq(j,k)+ssnow%wmice(j,k))
               if (ssnow%tgg(j,k) .gt. 273.15) then
                  tmp = max(0._r_2, log10(ew(j)) + 1.0_r_2)
                  tmpsat = tkmg(j,k)*0.57_r_2**soil%watsat(j,k)
                else                               ! Frozen soil
                  tmp =  ew(j)
                  tmpsat = tkmg(j,k)*0.249_r_2**(liq_frac*soil%watsat(k,j))*2.29_r_2**soil%watsat(k,j)
                endif
                ccnsw(j,k) = tmp*tmpsat + (1._r_2-tmp)*soil%cnsd(j)
             else 
                ccnsw(j,k) = soil%cnsd(j)
             endif

        
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
      
      ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%watsat(:,k) ) * soil%css * soil%rhosoil   &
                          + soil%watsat(:,k) * ( wblfsp * cswat * rhowat +            &
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

         ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%watsat(:,k) ) * soil%css * soil%rhosoil&
                             + soil%watsat(:,k) * ( wblfsp * cswat * rhowat +         &
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
         
         ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%watsat(:,k) ) * soil%css *             &
                             soil%rhosoil + soil%watsat(:,k) * ( wblfsp * cswat *     &
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

END SUBROUTINE GWstempv

! -----------------------------------------------------------------------------
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
   
   TYPE(soil_parameter_type), INTENT(IN)    :: soil
   
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
      
         IF( soil%isoilm(j) .eq. 9 ) THEN
            ! permanent ice: fix hard-wired number in next version
            ccnsw(j,k) = snow_ccnsw
         ELSE
            ew(j) = ssnow%wblf(j,k) * soil%watsat(j,k)
            exp_arg = ( ew(j) * LOG( 60.0 ) ) + ( ssnow%wbfice(j,k)            &
                      * soil%watsat(j,k) * LOG( 250.0 ) )

            IF( exp_arg > 30 ) direct2min = .TRUE.
            
            IF( direct2min) THEN
               
               ccnsw(j,k) = 1.5 * MAX( 1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5 *      &
                            soil%watsat(j,k) /                                     &
                            MIN( ew(j), 0.5_r_2 * soil%watsat(j,k) ) ) ) )

            ELSE         
               
               ccnsw(j,k) = MIN( soil%cnsd(j) * EXP( exp_arg ), 1.5_r_2 )      &
                            * MAX( 1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5 *          &
                            soil%watsat(j,k) /                                     &
                            MIN( ew(j), 0.5_r_2 * soil%watsat(j,k) ) ) ) )
            
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
      
      ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%watsat(:,k) ) * soil%css * soil%rhosoil   &
                          + soil%watsat(:,k) * ( wblfsp * cswat * rhowat +            &
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

         ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%watsat(:,k) ) * soil%css * soil%rhosoil&
                             + soil%watsat(:,k) * ( wblfsp * cswat * rhowat +         &
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
         
         ssnow%gammzz(:,k) = MAX( ( 1.0 - soil%watsat(:,k) ) * soil%css *             &
                             soil%rhosoil + soil%watsat(:,k) * ( wblfsp * cswat *     &
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
   
   TYPE(soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
   
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
            IF( (soil%isoilm(j) .eq. 9) .AND. ktau_gl <= 2 )                       &
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
               IF( (soil%isoilm(j) .eq. 9) .AND. ktau_gl <= 2 ) THEN
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

SUBROUTINE GWsoilfreeze(dels, soil, ssnow)
   USE cable_common_module
   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(soil_parameter_type), INTENT(IN)    :: soil
   REAL(r_2), DIMENSION(mp)           :: sicefreeze
   REAL(r_2), DIMENSION(mp)           :: sicemelt
   REAL, DIMENSION(mp)                :: xx, ice_mass,liq_mass,tot_mass
   INTEGER k
   REAL(r_2),DIMENSION(mp,ms) :: frozen_limit,iceF  !Decker and Zeng 2009

   where (ssnow%tgg .le. C%TFRZ)
      frozen_limit(:,:) = (1. - exp(-2.*(ssnow%wb(:,:)/soil%watsat(:,:))**4.0 *&
                       (ssnow%tgg(:,:)-273.16)))/exp(1. - ssnow%wb(:,:)/soil%watsat(:,:))
      frozen_limit(:,:) = max(0.5,frozen_limit(:,:))
   elsewhere (ssnow%tgg .gt. C%TFRZ)
      frozen_limit(:,:) = 0.
   endwhere

   !allow more freezing for permenant glacier ice regions
   where ( spread(soil%isoilm,2,ms) .eq. 9 ) frozen_limit(:,:) = 0.95

   xx = 0.
   DO k = 1, ms

      ice_mass = ssnow%wbice(:,k)*real(soil%zse(k)*C%denice,r_2)
      liq_mass = ssnow%wbliq(:,k)*real(soil%zse(k)*C%denliq,r_2)
      tot_mass = liq_mass + ice_mass
      
      WHERE (ssnow%tgg(:,k) .lt. C%TFRZ .and. frozen_limit(:,k) * ssnow%wb(:,k) - ssnow%wbice(:,k) > .001)
         
         sicefreeze = MIN( MAX( 0.0_r_2, ( frozen_limit(:,k) * ssnow%wb(:,k) -      &
                      ssnow%wbice(:,k) ) ) * soil%zse(k) * C%denice,             &
                      ( C%TFRZ - ssnow%tgg(:,k) ) * ssnow%gammzz(:,k) / C%HLF )
         ssnow%wbice(:,k) = MIN( ssnow%wbice(:,k) + sicefreeze / (soil%zse(k)  &
                            * 1000.0), frozen_limit(:,k) * ssnow%wb(:,k) )
         xx = soil%css * soil%rhosoil
         ssnow%gammzz(:,k) = MAX(                                              &
             REAL((1.0 - soil%watsat(:,k)) * soil%css * soil%rhosoil ,r_2)            &
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
              REAL((1.0-soil%watsat(:,k)) * soil%css * soil%rhosoil,r_2)             &
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

END SUBROUTINE GWsoilfreeze

!-----------------------------------------------------------------------------------------

! -----------------------------------------------------------------------------

SUBROUTINE remove_trans(dels, soil, ssnow, canopy, veg)
   
   USE cable_common_module, ONLY : redistrb

   ! Removes transpiration water from soil.
   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(canopy_type), INTENT(INOUT)         :: canopy
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(soil_parameter_type), INTENT(IN)    :: soil
   TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
   REAL(r_2), DIMENSION(mp,0:ms+1) :: diff 
   REAL(r_2), DIMENSION(mp)      :: xx,xxd,evap_cur
   INTEGER :: k,i
 
 
   xx = 0._r_2; xxd = 0._r_2; diff(:,:) = 0._r_2
   DO k = ms,1,-1  !I like the idea of removing from the bottom first.shouldn't matter
                   !should layers that are wetter get more water removed?
                   !logical but is it supported by evidence?
      DO i=1,mp

         if (canopy%fevc(i) .gt. 0._r_2) then

            xx(i) = canopy%fevc(i) * dels / C%HL * veg%froot(i,k) + diff(i,k+1)
            diff(i,k) = max(0._r_2,ssnow%wbliq(i,k)-soil%wiltp(i,k)) &
                       * real(soil%zse(k)*C%denliq,r_2)
            xxd(i) = xx(i) - diff(i,k)

            if (xxd(i) .gt. 0._r_2) then
               ssnow%wbliq(i,k) = ssnow%wbliq(i,k) - diff(i,k)/real(soil%zse(k)*C%denliq,r_2)  !volume
               diff(i,k) = xxd(i)
            else
               ssnow%wbliq(i,k) = ssnow%wbliq(i,k) - xx(i)/real(soil%zse(k)*C%denliq,r_2)  !volume
               diff(i,k) = 0._r_2
            end if

            ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*real(soil%zse(k)*C%denliq,r_2)  !mass
            ssnow%wmtot(i,k) = ssnow%wmliq(i,k) + ssnow%wmice(i,k)  !mass
            ssnow%wb(i,k)    = ssnow%wbliq(i,k) + ssnow%wbice(i,k)  !volume

          end if  !fvec > 0

      END DO  !mp
   END DO     !ms

END SUBROUTINE remove_trans 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!MD GW code from here on!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  !-------------------------------------------------------------------------
  SUBROUTINE ovrlndflx (dels, ktau, ssnow, soil,veg, md_prin )
  USE cable_common_module

  IMPLICIT NONE
    REAL, INTENT(IN)                         :: dels ! integration time step (s)
    INTEGER, INTENT(IN)                      :: ktau ! integration step number
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(IN)    :: soil ! soil parameters
    TYPE(veg_parameter_type) , INTENT(IN)    :: veg  ! veg parameters
    LOGICAL, INTENT(IN)                      :: md_prin
    INTEGER, PARAMETER                       :: ntest = 0 ! for snow diag prints
    INTEGER, PARAMETER                       :: nglacier = 2 ! 0 original, 1 off, 2 new Eva
    INTEGER                                  :: k, i, j
    REAL, DIMENSION(mp)                :: rnof5
    REAL, DIMENSION(mp)                :: sgamm
    REAL, DIMENSION(mp)                :: smasstot
    REAL, DIMENSION(mp,0:3)            :: smelt1                   !snow melt
    REAL(r_2), DIMENSION(mp)           :: icef,efpor               !tmp vars, fraction of ice in gridcell
    REAL(r_2)                          :: tmpa,tmpb,qinmax         !tmp vars, maximum infiltration [mm/s]
    REAL(r_2), DIMENSION(mp)           :: satfrac_liqice,S       !saturated fraction of cell, wtd in m
    REAL(r_2)                          :: liqmass,icemass,totmass  !liquid mass,ice mass, total mass [mm]
    REAL(r_2), parameter               :: pi=3.1415926535898
    REAL(r_2)                          :: fice
    REAL(r_2)                          :: dzmm,slopeSTDmm
    logical                            :: prinall = .false.  !for debugging
    
   !For now assume there is no puddle
   dzmm = 1000._r_2 * soil%zse(1)

   do i = 1,mp
      efpor(i) = max(0.001_r_2, soil%watsat(i,1) - ssnow%wbice(i,1))
      icemass  = ssnow%wbice(i,1) * dzmm * dri
      liqmass  = (ssnow%wb(i,1)-ssnow%wbice(i,1)) * dzmm
      totmass  = max(liqmass+icemass,real(1e-2,r_2))
      !icef(i)     = max(0._r_2,min(1._r_2,1.25_r_2*icemass / totmass))
      icef(i)     = max(0._r_2,min(1._r_2,gw_params%IceBeta*icemass / totmass))
   end do
   S(:) = 0._r_2
   do k=1,ms
     S(:) = S(:) + max(0.01,min(1.0, (ssnow%wb(:,k)-ssnow%wbice(:,k)-soil%watr(:,k))/max(0.001,soil%watsat(:,k)-soil%watr(:,k)) ) )*soil%zse(k)
   end do
   S(:) = S(:)/sum(soil%zse(1:ms),dim=1)
   !srf frozen fraction.  should be based on topography
   do i = 1,mp
      !fice = (exp(-3._r_2*(1._r_2-icef(i)))-exp(-3._r_2))/(1._r_2-exp(-3._r_2))
      fice = (exp(gw_params%IceAlpha*(1._r_2-icef(i)))-exp(gw_params%IceAlpha))/(1._r_2-exp(gw_params%IceAlpha))
      fice  = min(max(fice,0._r_2),1._r_2)
      !Saturated fraction
       if (gw_params%MaxSatFraction .gt. 1e-7) then 
          slopeSTDmm = sqrt(max(gw_params%MaxSatFraction*soil%slope_std(i),1e-5)) ! ensure some variability
          ssnow%satfrac(i)    = max(1e-6,min(0.95,1._r_2 - erf( slopeSTDmm / sqrt(2.0* S(i)) ) ) )  
       else
          ssnow%satfrac(i) = 0. 
       end if
       satfrac_liqice(i)   = max(1e-6,min(0.95,fice + (1._r_2-fice)*ssnow%satfrac(i) ) )
   end do

   do i=1,mp
      tmpa = ssnow%wbliq(i,1) / efpor(i)
      tmpb = max( (tmpa-satfrac_liqice(i))/max(0.01_r_2,(1._r_2-satfrac_liqice(i))), 0._r_2)
      tmpa = -2._r_2*soil%clappB(i,1)*soil%smpsat(i,1)/dzmm
      qinmax = (1._r_2 + tmpa*(tmpb-1._r_2))*soil%hksat(i,1)*exp(-gw_params%hkrz*(0.5*dzmm/1000.0_r_2-gw_params%zdepth))

      ssnow%rnof1(i) = satfrac_liqice(i) * ssnow%fwtop(i) + &
                         (1._r_2-satfrac_liqice(i))*max((ssnow%fwtop(i)-qinmax) , 0._r_2)

      ssnow%fwtop(i) = ssnow%fwtop(i) - ssnow%rnof1(i)

   end do  !mp

  !add back to the lakes to keep saturated instead of drying
  where (veg%iveg .eq. 16)
     ssnow%fwtop(:) = ssnow%fwtop(:) + ssnow%rnof1(:)
     ssnow%rnof1(:) = 0._r_2
  endwhere
           
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

  REAL(r_2), DIMENSION(mp)            :: fz, wmean,ztot
  REAL(r_2), DIMENSION(mp,ms)         :: stot
  INTEGER                             :: k,i

  do i=1,mp
     wmean(i) = 0._r_2
     fz(i)    = 5._r_2
     ztot(i)  = 0._r_2
     stot(i,:) = (ssnow%wb(i,:)-soil%watr(i,:)) / (soil%watsat(i,:)-soil%watr(i,:))
  end do
  do k  = 1, ms
     do i=1,mp
        wmean(i) = wmean(i) + stot(i,k)*soil%zse(k)*1000._r_2
        ztot(i)  = ztot(i) + soil%zse(k)*1000._r_2
     end do
  end do

  do i=1,mp
     wmean(i) = wmean(i) + ssnow%GWwb(i)/soil%GWwatsat(i) * soil%GWdz(i)*1000._r_2
     ztot(i)  = ztot(i) + soil%GWdz(i)*1000._r_2

     ssnow%wtd(i) = min(200000._r_2, fz(i) * (ztot(i) - wmean(i)))
  end do

  END SUBROUTINE simple_wtd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


  !----------------------------------------------------------------------
  ! SUBROUTINE iterative_wtd
  !
  ! Iteratively calcs the water table depth by equating the mass of water in the
  ! soil column to the mass of a hydrostatic column inegrated from the surface to the 
  ! water table depth
  !  
  SUBROUTINE iterative_wtd (ssnow, soil, veg, ktau, md_prin)
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
  do i=1,mp
     defc(i) = (soil%watsat(i,ms))*(zimm(ms)+Nsmpsat(i)/(1._r_2-invB(i))*            &
             (1._r_2-((Nsmpsat(i)+zimm(ms))/Nsmpsat(i))**(1._r_2-invB(i)))) 
     defc(i) = max(0.1_r_2,defc(i)) 
  end do

  def(:) = 0._r_2
  do k=1,ms
     do i=1,mp

        if (soil%watsat(i,k) .gt. ssnow%wb(i,k)) then
          def(i) = def(i) +                                                           &
                   (soil%watsat(i,k)-(ssnow%wbliq(i,k)+dri*ssnow%wbice(i,k)))*dzmm_mp(i,k)
        end if
      end do  !mp
  end do  !ms

  do i=1,mp
    def(i) = def(i) + max(0._r_2,soil%GWwatsat(i)-ssnow%GWwb(i))*soil%GWdz(i)*1000._r_2
  end do   

  ssnow%wtd(:) = zimm(ms)*def(:)/defc(:)

  do i=1,mp

    if ((veg%iveg(i) .ne. 16) .and. (soil%isoilm(i) .ne. 9)) then      

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

        END DO mainloop  !defc .gt. def

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

        END DO mainloop2  !defc .lt. def

      else  !water table depth is exactly on bottom boundary

        ssnow%wtd(i) = zimm(ms)

      endif

    endif  !check veg and soils

  end do   !mp loop


  !limit wtd to be within a psecified range
  where (ssnow%wtd(:) .gt. wtd_max) ssnow%wtd(:) = wtd_max
  where (ssnow%wtd(:) .lt. wtd_min) ssnow%wtd(:) = wtd_min

  END SUBROUTINE iterative_wtd

  !-------------------------------------------------------------------------
  ! SUBROUTINE smoistgw (fwtop,dt,ktau,ssnow,soil,prin)
  ! solves the modified richards equation (Zeng and Decker 2009) to find
  ! vertical mocement of soil water.  Bottom boundary condition is determined
  ! using a single layer groundwater module
  !
  SUBROUTINE smoistgw (dels,ktau,ssnow,soil,veg,canopy,md_prin)
  USE cable_common_module

  IMPLICIT NONE
  
    REAL, INTENT(IN)                          :: dels  ! time step size (s)
    INTEGER, INTENT(IN)                       :: ktau  ! integration step number
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)     :: veg
    TYPE(canopy_type), INTENT(INOUT)          :: canopy ! vegetation variables
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
    REAL(r_2), DIMENSION(mp)            :: GWzimm,xs,zaq,s_mid,GWdzmm
    REAL(r_2), DIMENSION(mp)            :: xs1,GWmsliq!xsi    !mass (mm) of liquid over/under saturation, mass of aquifer water
    REAL(r_2)                           :: xsi
    REAL(r_2), DIMENSION(mp,ms+1)       :: qhlev,del_wb
    REAL(r_2), DIMENSION(mp)            :: sm_tot,drainmod  !total column soil water available for drainage
    INTEGER, DIMENSION(mp)              :: idlev
    logical                             :: prinall = .false.   !another debug flag
    character (len=30)                  :: fmt  !format to output some debug info
    !MD DEBUG VARS
    INTEGER :: imp,ims,k_drain

    drainmod(:) = 1._r_2  !parameter to modify qhrz params by basin or veg type
    !where(veg%iveg .eq. 2) drainmod(:) = 0.1_r_2*drainmod(:)
    !drainmod(:) = (1._r_2 + soil%FOrg(:,1))*drainmod(:)


    fmt='(A6,6(1X,F8.6))'       !not needed.  was used to nicely format debug output
    !make code cleaner define these here
    dzmm    = 1000.0_r_2 * real(soil%zse(:),r_2)
    dzmm_mp = spread(dzmm,1,mp)
    zimm(0) = 0._r_2

    zimm(1:ms) = zimm(0:(ms-1)) + dzmm(1:ms)
    zmm(1:ms)  = zimm(0:(ms-1)) + 0.5_r_2*dzmm(1:ms)

    GWdzmm(:) = real(soil%GWdz(:),r_2)*1000._r_2
    GWzimm(:) = zimm(ms)+GWdzmm(:)
    zaq(:)    = zimm(ms) + 0.5_r_2*GWdzmm(:)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! preset to allow for non-land & snow points in trimb
    old_wb(:,:) = ssnow%wb(:,:)
    
    !equilibrium water content
    CALL calc_equilibrium_water_content(ssnow,soil)

    CALL calc_soil_hydraulic_props(ssnow,soil,veg)

    do i=1,mp

       !Note: future revision will have interaction with river here. nned to
       !work on router and add river type cells
       ssnow%qhz(i)  = max(tan(soil%slope(i)),0.001) * drainmod(i)*gw_params%MaxHorzDrainRate* &
                    exp(-ssnow%wtd(i)/(1000._r_2*(gw_params%EfoldHorzDrainRate)))

       !Keep "lakes" saturated forcing qhz = 0.  runoff only from lakes
       !overflowing
       if (soil%isoilm(i) .eq. 9 .or. veg%iveg(i) .eq. 16) ssnow%qhz(i) = 0._r_2
 
       !identify first no frozen layer.  drinage from that layer and below
       !drain from sat layers
       k_drain = ms+1
       do k=ms,2,-1
           !below what was in paper
          !if ( ssnow%wbliq(i,k+1) .le. gw_params%frozen_frac*ssnow%wb(i,k+1)  .and.&
          !     ssnow%wbliq(i,k  ) .gt. gw_params%frozen_frac*ssnow%wb(i,k  ) ) then
          if (ssnow%wtd(i) .le. sum(dzmm(1:k),dim=1)) then
             k_drain = k
          end if
       end do
       k_drain = min(k_drain,2)

!       k_drain = ms + 1
!       do k=ms,2,-1
!          if ((all(ssnow%wbliq(i,k:min(k_drain,ms)) .ge. 0.95*soil%watsat(i,k:min(k_drain,ms)))) .and. &
!                 (ssnow%wbliq(i,k-1) .lt. 0.95*soil%watsat(i,k-1)))) then
!             k_drain = k
!          end if
!       end do

       qhlev(i,:) = 0._r_2
       sm_tot(i) = 0._r_2
       if (k_drain .le. ms) then
          do k=k_drain,ms
             sm_tot(i) = sm_tot(i) + max(ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)!*dzmm(k)
          end do
          sm_tot(i) = max(sm_tot(i),0.01_r_2)

         do k=k_drain,ms
             qhlev(i,k) = ssnow%qhz(i)*max(ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)/sm_tot(i)!*dzmm(k)/sm_tot(i)
          end do

       else
          qhlev(i,ms+1) = max(1._r_2-ssnow%fracice(i,ms),0._r_2)*ssnow%qhz(i)!*max(ssnow%GWwb(i)-soil%watr(i,ms),0._r_2)/sm_tot(i)
       end if

       !incase every layer is frozen very dry
       ssnow%qhz(i) = qhlev(i,ms+1)
       do k=k_drain,ms
          ssnow%qhz(i) = ssnow%qhz(i) +qhlev(i,k)
       end do

    end do  


    rt(:,:) = 0._r_2; at(:,:) = 0._r_2     !ensure input to tridiag is valid
    bt(:,:) = 1._r_2; ct(:,:) = 0._r_2

    k = 1     !top soil layer
    do i=1,mp
       qin(i)     = ssnow%sinfil(i)
       den(i)     = (zmm(k+1)-zmm(k))
       dne(i)     = (ssnow%zq(i,k+1)-ssnow%zq(i,k))
       num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
       dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k))/den(i)
       dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
       rt(i,k) =  qin(i) - qout(i)! - qhlev(i,k) - ssnow%rex(i,k)
       at(i,k) =  0._r_2
       bt(i,k) =  dzmm(k)/dels + dqodw1(i)
       ct(i,k) =  dqodw2(i)      
    end do
    do k = 2, ms - 1     !middle soil layers
       do i=1,mp
          den(i)     = (zmm(k) - zmm(k-1))
          dne(i)     = (ssnow%zq(i,k)-ssnow%zq(i,k-1))
          num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
          qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
          dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          den(i)     = (zmm(k+1)-zmm(k))
          dne(i)     = (ssnow%zq(i,k+1)-ssnow%zq(i,k))
          num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
          qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
          dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k))/den(i)
          dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
          rt(i,k) =  qin(i) - qout(i)! - qhlev(i,k) - ssnow%rex(i,k)
          at(i,k) = -dqidw0(i)
          bt(i,k) =  dzmm(k)/dels - dqidw1(i) + dqodw1(i)
          ct(i,k) =  dqodw2(i)
       end do
    end do
       
    k = ms   !Bottom soil layer
    do i=1,mp
       den(i)     = (zmm(k) - zmm(k-1))
       dne(i)     = (ssnow%zq(i,k)-ssnow%zq(i,k-1))
       num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
       qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
       dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       den(i)     = zaq(i) - zmm(k)
       dne(i)     = (ssnow%GWzq(i)-ssnow%zq(i,k))
       num(i)     =  (ssnow%GWsmp(i)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = 0._r_2
       dqodw1(i)  = 0._r_2
       dqodw2(i)  = 0._r_2
       rt(i,k) =  qin(i) - qout(i)
       at(i,k) = -dqidw0(i)
       bt(i,k) =  dzmm(k)/dels - dqidw1(i) + dqodw1(i)
       ct(i,k) =  dqodw2(i) 
    end do
       
    !Doing the recharge outside of the soln of Richards Equation makes it easier to track total recharge amount.
    !Add to ssnow at some point 
    do i=1,mp
       if ((ssnow%wtd(i) .ge. sum(dzmm,dim=1)) .and. (veg%iveg(i) .ne. 17) .and. (soil%isoilm(i) .ne. 9))  then
          ssnow%Qrecharge(i) = -ssnow%hk(i,ms)*((ssnow%GWsmp(i)-ssnow%smp(i,k-1)) - (ssnow%GWzq(i)-ssnow%zq(i,k)))/(zaq(i) - zmm(k))
       else
          ssnow%Qrecharge(i) = 0._r_2
       end if
    end do

    CALL trimb(at,bt,ct,rt,ms)                       !use the defulat cable tridiag solution

    do k=1,ms
       do i=1,mp
          ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + rt(i,k) - qhlev(i,k)*dels/dzmm(k)   !volutermic liquid
       end do
    end do

    do i=1,mp
       ssnow%wbliq(i,ms) = ssnow%wbliq(i,ms) - ssnow%Qrecharge(i)*dels/dzmm(ms)
    end do
    do i=1,mp
       ssnow%GWwb(i) = ssnow%GWwb(i)  +  (ssnow%Qrecharge(i)-qhlev(i,ms+1))*dels/GWdzmm(i)
    end do

    !determine the available pore space
    !volumetric
    eff_por(:,:)  = soil%watsat(:,:) - ssnow%wbice(:,:)
    do i=1,mp
       xsi = 0._r_2

       if (ssnow%GWwb(i) .gt. soil%GWwatsat(i)) then
          xsi = (ssnow%GWwb(i) - soil%GWwatsat(i))*GWdzmm(i)
          ssnow%GWwb(i) = soil%GWwatsat(i)
       end if

       do k=1,ms
          if (ssnow%wbliq(i,k) .gt. eff_por(i,k)) then
             xsi = xsi + (ssnow%wbliq(i,k) - eff_por(i,k))*dzmm(k)
             ssnow%wbliq(i,k) = eff_por(i,k)
           end if
       end do
     
       do k = ms,1,-1  !loop from bottom to top adding extra water to each layer
          if (xsi .gt. 0._r_2) then
             if (xsi .lt. (eff_por(i,k)-ssnow%wbliq(i,k))*dzmm(k)) then
                ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + xsi/dzmm(k)
                xsi = 0._r_2
             else
                xsi = xsi - (eff_por(i,k) - ssnow%wbliq(i,k))*dzmm(k)
                ssnow%wbliq(i,k) = eff_por(i,k)
             end if
          end if
       end do  !ms loop

       if (xsi .gt. 0._r_2) then
          ssnow%qhz(i) = ssnow%qhz(i) + xsi/dels
          xsi = 0._r_2
       end if

       do k = 1,ms
          xsi = 0._r_2             !should be a single float (array not needed)
          if (ssnow%wbliq(i,k) .lt. volwatmin) then
             xsi = (volwatmin - ssnow%wbliq(i,k))*dzmm(k)  !in mm
             ssnow%wbliq(i,k) = volwatmin
             if (k .lt. ms) then
                ssnow%wbliq(i,k+1) = ssnow%wbliq(i,k+1) - xsi/dzmm(k+1)
             else
                ssnow%GWwb(i) = ssnow%GWwb(i) - xsi / GWdzmm(i)
             end if
          end if
       end do  !ms loop
 
       if ( (ssnow%GWwb(i) .lt. volwatmin) .and. (soil%isoilm(i) .ne. 9) ) then
          xsi = (volwatmin - ssnow%GWwb(i)) / GWdzmm(i)  !mm
          ssnow%GWwb(i) = volwatmin
          ssnow%qhz(i) = ssnow%qhz(i) - xsi / dels
       end if
         
       !update mass variables
       ssnow%wmliq(i,:)      = ssnow%wbliq(i,:) * dzmm_mp(i,:)
       ssnow%wmice(i,:)      = ssnow%wbice(i,:) * dzmm_mp(i,:) * dri
       ssnow%wb(i,:)         = ssnow%wbliq(i,:) + ssnow%wbice(i,:)
       ssnow%wmtot(i,:)      = ssnow%wmliq(i,:) + ssnow%wmice(i,:)
       ssnow%rnof2(i)        = ssnow%qhz(i)               !rnof2 is default cable deep runoff var 

   end do  !mp loop
       

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
   REAL                     , INTENT(IN)     :: dels ! integration time step (s)
   TYPE(soil_parameter_type), INTENT(INOUT)  :: soil
   TYPE(soil_snow_type)     , INTENT(INOUT)  :: ssnow
   TYPE(canopy_type)        , INTENT(INOUT)  :: canopy
   TYPE(veg_parameter_type) , INTENT(INOUT)  :: veg
   TYPE(met_type)           , INTENT(INOUT)  :: met ! all met forcing
   TYPE (balances_type)     , INTENT(INOUT)  :: bal

   INTEGER             :: k,i
   REAL, DIMENSION(mp) :: snowmlt
   REAL, DIMENSION(mp) :: totwet
   REAL, DIMENSION(mp) :: weting,GWwb_ic,wberr
   REAL, DIMENSION(mp) :: tgg_old, tggsn_old,wbtot_ic,del_wbtot
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
         ssnow%wb(:,:)  = MIN( soil%watsat(:,:), MAX ( ssnow%wb(:,:), soil%wiltp(:,:) ) )   

         DO k = 1, ms
            
            WHERE( ssnow%tgg(:,k) <= C%TFRZ .AND. ssnow%wbice(:,k) <= 0.01 )   &
               ssnow%wbice(:,k) = 0.5 * ssnow%wb(:,k)

            WHERE( ssnow%tgg(:,k) < C%TFRZ)                                    &
               ssnow%wbice(:,k) = 0.8 * ssnow%wb(:,k)
            
         END DO

         WHERE ( soil%isoilm .eq. 9 ) 

            ! permanent ice: fix hard-wired number in next version
            ssnow%snowd = max_glacier_snowd
            ssnow%osnowd = max_glacier_snowd
            ssnow%tgg(:,1) = ssnow%tgg(:,1) - 1.0

         END WHERE

         WHERE ( spread(soil%isoilm,2,ms) .eq. 9 )

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

   ssnow%wblf   = max(ssnow%wbliq/soil%watsat,0.00001_r_2)
   ssnow%wbfice = max(ssnow%wbice/soil%watsat,0.00001_r_2)

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
   !CALL GWstempv(dels, canopy, ssnow, soil)

   !do the soil and snow melting, freezing prior to water movement
   ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

   CALL snow_melting (dels, snowmlt, ssnow, soil )
   
   ! Add new snow melt to global snow melt variable: 
   ssnow%smelt = ssnow%smelt + snowmlt

   CALL remove_trans(dels, soil, ssnow, canopy, veg)        !transpiration loss per soil layer

   CALL  GWsoilfreeze(dels, soil, ssnow)

   ssnow%fwtop = canopy%precis/dels + ssnow%smelt/dels   !water from canopy and snowmelt [mm/s]   
   !ssnow%rnof1 = ssnow%rnof1 + ssnow%smelt / dels          !adding snow melt directly to the runoff

   CALL iterative_wtd (ssnow, soil, veg, ktau, md_prin)  
   !CALL simple_wtd(ssnow, soil, veg, ktau, md_prin)

   CALL ovrlndflx (dels, ktau, ssnow, soil, veg, md_prin )         !surface runoff, incorporate ssnow%pudsto?

   ssnow%sinfil = ssnow%fwtop - canopy%segg  !canopy%fes/C%HL               !remove soil evap from throughfall
   !ssnow%pudsto = max(ssnow%pudsto - canopy%fesp/C%HL*dels,0._r_2)  !currently pudsto = 0.0 always

   CALL smoistgw (dels,ktau,ssnow,soil,veg,canopy, md_prin)               !vertical soil moisture movement. 
  
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


SUBROUTINE calc_equilibrium_water_content(ssnow,soil)

  IMPLICIT NONE
  
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    !local variables
    REAL(r_2), dimension(mp)    :: zaq      !node depth of the aquifer
    REAL(r_2), dimension(ms)    :: dzmm     !layer thickness for single tile
    REAL(r_2), dimension(mp)    :: GWdzmm   !aquifer thickness at each tile
    REAL(r_2), dimension(mp)    :: GWzimm   !aquifer layer interface depth
    REAL(r_2), dimension(0:ms)  :: zimm     !layer interface depth in mm  
    REAL(r_2), dimension(ms)    :: zmm      !node depths in mm
    REAL(r_2)                   :: tempi, temp0,voleq1,wbrat

    INTEGER :: k,i

    !make code cleaner define these here
    dzmm    = 1000.0_r_2 * real(soil%zse(:),r_2)
    zimm(0) = 0._r_2

    do k=1,ms
       zimm(k) = zimm(k-1) + dzmm(k)
       zmm(k)  = zimm(k-1) + 0.5_r_2*dzmm(k)
    end do 

    do i=1,mp
       GWdzmm(i) = real(soil%GWdz(i),r_2)*1000._r_2
       GWzimm(i) = zimm(ms)+GWdzmm(i)
       zaq(i)    = zimm(ms) + 0.5_r_2*GWdzmm(i)
    end do
    !equilibrium water content
    do k=1,ms
       do i=1,mp

          if ((ssnow%wtd(i) .le. zimm(k-1))) then         !fully saturated

             ssnow%wbeq(i,k) = soil%watsat(i,k)
             
          elseif ((ssnow%wtd(i) .le. zimm(k)) .and. (ssnow%wtd(i) .gt. zimm(k-1))) then

             tempi = 1._r_2
             temp0 = (((soil%smpsat(i,k)+ssnow%wtd(i)-zimm(k-1))/soil%smpsat(i,k)))**(1._r_2-1._r_2/soil%clappB(i,k))               
             voleq1 = -soil%smpsat(i,k)*(soil%watsat(i,k)-soil%watr(i,k))/&
                       (1._r_2-1._r_2/soil%clappB(i,k))/(ssnow%wtd(i)-zimm(k-1))*(tempi-temp0)
             ssnow%wbeq(i,k) = (voleq1*(ssnow%wtd(i)-zimm(k-1)) + (soil%watsat(i,k)-soil%watr(i,k))&
                            *(zimm(k)-ssnow%wtd(i)))/(zimm(k)-zimm(k-1)) + soil%watr(i,k)

          else

             tempi = (((soil%smpsat(i,k)+ssnow%wtd(i)-zimm(k))/soil%smpsat(i,k)))**(1._r_2-1._r_2/soil%clappB(i,k))
             temp0 = (((soil%smpsat(i,k)+ssnow%wtd(i)-zimm(k-1))/soil%smpsat(i,k)))**(1._r_2-1._r_2/soil%clappB(i,k))   
             ssnow%wbeq(i,k) = -soil%smpsat(i,k)*(soil%watsat(i,k)-soil%watr(i,k))/&
                               (1._r_2-1._r_2/soil%clappB(i,k))/(zimm(k)-zimm(k-1))*(tempi-temp0)+soil%watr(i,k)

          end if
          
          ssnow%wbeq(i,k) = min(max(ssnow%wbeq(i,k),soil%watr(i,k)),soil%watsat(i,k))
          
          wbrat = min(max((ssnow%wbeq(i,k) - soil%watr(i,k))/(soil%watsat(i,k)-soil%watr(i,k)),&
                          0.001_r_2),1._r_2)

          ssnow%zq(i,k) = max(-soil%smpsat(i,k)*(wbrat**(-soil%clappB(i,k))),sucmin)
          
          
       end do  !mp
    end do  !ms
 
    do i=1,mp
    !Aquifer Equilibrium water content
       if (ssnow%wtd(i) .le. zimm(ms)) then      !fully saturated

          ssnow%GWwbeq(i) = soil%GWwatsat(i)-soil%GWwatr(i)

       elseif ((ssnow%wtd(i) .gt. GWzimm(i)))   then     !fully unsaturated

          tempi = (((soil%GWsmpsat(i)+ssnow%wtd(i)-GWzimm(i))/soil%GWsmpsat(i)))**(1._r_2-1._r_2/soil%GWclappB(i))
          temp0 = (((soil%GWsmpsat(i)+ssnow%wtd(i)-zimm(ms))/soil%GWsmpsat(i)))**(1._r_2-1._r_2/soil%GWclappB(i))   
          ssnow%GWwbeq(i) = -soil%GWsmpsat(i)*soil%GWwatsat(i)/&
                          (1._r_2-1._r_2/soil%GWclappB(i))/(GWzimm(i)-zimm(ms))*(tempi-temp0) + soil%GWwatr(i)   

       else           

          tempi  = 1._r_2
          temp0  = (((soil%GWsmpsat(i)+ssnow%wtd(i)-zimm(ms))/soil%GWsmpsat(i)))**(1._r_2-1._r_2/soil%GWclappB(i))               
          voleq1 = -soil%GWsmpsat(i)*(soil%GWwatsat(i)-soil%GWwatr(i))/&
                (1._r_2-1._r_2/soil%GWclappB(i))/(ssnow%wtd(i)-zimm(ms))*(tempi-temp0) + soil%GWwatr(i)
          ssnow%GWwbeq(i) = (voleq1*(ssnow%wtd(i)-zimm(ms)) + (soil%GWwatsat(i)-soil%GWwatr(i))*&
                         (GWzimm(i)-ssnow%wtd(i)))/(GWzimm(i)-zimm(ms)) + soil%GWwatr(i)

       end if

       ssnow%GWwbeq(i) = min(max(ssnow%GWwbeq(i),soil%GWwatr(i)),soil%GWwatsat(i))

       ssnow%GWzq(i) = -soil%GWsmpsat(i)*(max((ssnow%GWwbeq(i)-soil%GWwatr(i))/           &
                    (soil%GWwatsat(i)-soil%GWwatr(i)),0.001_r_2))**(-soil%GWclappB(i))
       ssnow%GWzq(i) = max(sucmin, ssnow%GWzq(i))
       
    end do

END SUBROUTINE calc_equilibrium_water_content

SUBROUTINE calc_srf_wet_fraction(ssnow,soil)

  IMPLICIT NONE
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(IN)    :: soil ! soil parameters

    !local variables
    REAL(r_2), DIMENSION(mp)           :: icef,satfrac_liqice,S
    REAL(r_2)                          :: fice,xx
    REAL(r_2)                          :: dzmm_one,liqmass,icemass,totmass
    INTEGER                            :: i,j,k
    REAL(r_2), parameter               :: pi=3.1415926535898
    INTEGER,   parameter               :: max_iter = 10
    REAL(r_2)                          :: wb_unsat,wb_lin,funcval,derv,slopeSTDmm,func_step
    REAL(r_2)                          :: wb_evap_threshold

    dzmm_one  = 1000._r_2 * real(soil%zse(1),r_2)

    do i = 1,mp
       icemass  = ssnow%wbice(i,1) * dzmm_one * dri
       liqmass  = (ssnow%wb(i,1)-ssnow%wbice(i,1)) * dzmm_one
       totmass  = max(liqmass+icemass,real(1e-2,r_2))
       icef(i)     = max(0._r_2,min(1._r_2, gw_params%IceBeta*icemass / totmass))
   end do

   !S(:) = 0._r_2
   !do k=1,2
   !   S(:) = S(:) + max(0.01,min(1.0, (ssnow%wb(:,k)-ssnow%wbice(:,k)-soil%watr(:,k))/max(0.001,soil%watsat(:,k)-soil%watr(:,k)) ) )
   !end do
   !S(:) = S(:)/2._r_2

   S(:) = 0._r_2
   do k=1,ms
     S(:) = S(:) + max(0.01,min(1.0,(ssnow%wb(:,k)-ssnow%wbice(:,k)-soil%watr(:,k))/max(0.001,soil%watsat(:,k)-soil%watr(:,k))) )*soil%zse(k)
   end do
   S(:) = S(:)/sum(soil%zse(1:ms),dim=1)

   !srf frozen fraction.  should be based on topography
   do i = 1,mp
      fice = (exp(gw_params%IceAlpha*(1._r_2-icef(i)))-exp(gw_params%IceAlpha))/(1._r_2-exp(gw_params%IceAlpha))
      fice = min(1._r_2,max(0._r_2,fice))

      !Saturated fraction
       if (gw_params%MaxSatFraction .gt. 1e-7) then 
         slopeSTDmm = sqrt(max(gw_params%MaxSatFraction*soil%slope_std(i),1e-2)) ! ensure some variability
         ssnow%satfrac(i)    = 1._r_2 - erf( slopeSTDmm / sqrt(2.0* S(i)) )
      else
         ssnow%satfrac(i)  = 0._r_2
      end if
      ssnow%satfrac(i)    = max(1e-6,min(0.95,ssnow%satfrac(i)))
      satfrac_liqice(i) = fice + (1._r_2-fice)*ssnow%satfrac(i)

      wb_unsat = ((ssnow%wb(i,1)-ssnow%wbice(i,1)) - ssnow%satfrac(i)*soil%watsat(i,1))/(1.-ssnow%satfrac(i))
      wb_unsat = min(soil%watsat(i,1),max(0.,wb_unsat))

      wb_evap_threshold = min( max( gw_params%SoilEvapAlpha*soil%fldcap(i,1), soil%wiltp(i,1) ), soil%watsat(i,1) )

      !Sakguchi and Zeng 2009
      if (wb_unsat .ge. wb_evap_threshold) then
         xx = 1.
      else
         xx = 0.25 * (1._r_2 - cos(pi*wb_unsat/(wb_evap_threshold)))**2.0
      end if

      if (.not.cable_user%or_evap) then
         ssnow%wetfac(i) = max(0.0,min(1.0,satfrac_liqice(i) +&
                                     (1. - satfrac_liqice(i))*xx ) )
      else
         ssnow%wetfac(i) = ssnow%satfrac(i)
      end if

   end do


END SUBROUTINE calc_srf_wet_fraction

SUBROUTINE calc_soil_hydraulic_props(ssnow,soil,veg)
   USE cable_common_module
   TYPE(soil_parameter_type), INTENT(IN)     :: soil 
   TYPE(soil_snow_type)     , INTENT(INOUT)  :: ssnow
   TYPE(veg_parameter_type) , INTENT(IN)     :: veg

   INTEGER :: i,k,kk

   REAL(r_2), DIMENSION(mp) :: s1, &  !temporary variables for calculating hydraulic properties
                               s2, &
                               s_mid

   REAL(r_2), DIMENSION(0:ms) :: zimm  !depths at interface between layers

    !soil matric potential, hydraulic conductivity, and derivatives of each with respect to water (calculated using total (not liquid))

   zimm(0)  = 0.0_r_2                                          !depth of layer interfaces mm
   do k=1,ms
     zimm(k) = zimm(k-1) + soil%zse(k)*1000._r_2
   end do
   zimm(ms) = zimm(ms) + soil%GWdz(1)*1000._r_2


    do k=1,ms
       do i=1,mp
          ssnow%icefrac(i,k) = ssnow%wbice(i,k)/(max(ssnow%wb(i,k),0.01_r_2))
          ssnow%fracice(i,k) = (exp(gw_params%IceAlpha*(1._r_2-ssnow%icefrac(i,k)))-exp(gw_params%IceAlpha))/(1._r_2-exp(gw_params%IceAlpha))
       end do
    end do

    ssnow%fracice(:,:) = max( min( ssnow%fracice, 1._r_2), 0._r_2)

    do k=1,ms-1
       kk=k+1
       do i=1,mp

          s1(i) = 0.5_r_2*(max(ssnow%wb(i,k)-soil%watr(i,k),0.) + max(ssnow%wb(i,kk)-soil%watr(i,kk),0.)) / &
            (0.5_r_2*((soil%watsat(i,k)-soil%watr(i,k)) + (soil%watsat(i,kk)-soil%watr(i,kk))))

          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)

          s2(i) = soil%hksat(i,k)*s1(i)**(2._r_2*soil%clappB(i,k)+2._r_2)
          ssnow%hk(i,k)    =  (1.0 - 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))*s1(i)*s2(i)*&
                               exp(-gw_params%hkrz*(zimm(k)/1000.0_r_2-gw_params%zdepth)) 
          ssnow%dhkdw(i,k) = (1.0 - 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))* &
                             (2._r_2*soil%clappB(i,k)+3._r_2)*s2(i)*0.5_r_2/(soil%watsat(i,k)-soil%watr(i,k))*&
                             exp(-gw_params%hkrz*(zimm(k)/1000.0_r_2-gw_params%zdepth))
       end do
    end do

    k = ms 
       do i=1,mp
          s1(i) = 0.5_r_2*(max(ssnow%wb(i,k)-soil%watr(i,k),0.) + max(ssnow%GWwb(i)-soil%GWwatr(i),0.)) / &
                  (0.5_r_2*(soil%watsat(i,k)-soil%watr(i,k) + soil%GWwatsat(i)-soil%GWwatr(i)))

          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)

          s2(i) = soil%hksat(i,k)*s1(i)**(2._r_2*soil%clappB(i,k)+2._r_2)
          ssnow%hk(i,k)    = s1(i)*s2(i)*(1.-ssnow%fracice(i,k))*&
                             exp(-gw_params%hkrz*(zimm(k)/1000._r_2-gw_params%zdepth))
          ssnow%dhkdw(i,k) = (1.-ssnow%fracice(i,k))* (2._r_2*soil%clappB(i,k)+3._r_2)*&
                             s2(i)*0.5_r_2/(soil%watsat(i,k)-soil%watr(i,k))*&
                             exp(-gw_params%hkrz*(zimm(k)/1000._r_2-gw_params%zdepth))
       end do
 
    do k=1,ms 
       do i=1,mp
          s_mid(i) = (ssnow%wb(i,k)-soil%watr(i,k))/&  !+dri*ssnow%wbice(:,k)
              (soil%watsat(i,k)-soil%watr(i,k))

          s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)

          ssnow%smp(i,k) = -soil%smpsat(i,k)*s_mid(i)**(-soil%clappB(i,k))

          ssnow%smp(i,k) = max(min(ssnow%smp(i,k),-soil%smpsat(i,k)),sucmin)

          ssnow%dsmpdw(i,k) = -soil%clappB(i,k)*ssnow%smp(i,k)/&
                    (max(s_mid(i)*(soil%watsat(i,k)-soil%watr(i,k)),0.001_r_2))       
       end do   
    end do

    do i=1,mp
       !Aquifer properties
       s_mid(i) = (ssnow%GWwb(i)-soil%GWwatr(i))/(soil%GWwatsat(i)-soil%GWwatr(i))
       s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)
       s2(i)    = soil%GWhksat(i)*s_mid(i)**(2._r_2*soil%clappB(i,ms)+2._r_2)

       ssnow%GWhk(i)     =s_mid(i)*s2(i)*(1._r_2-ssnow%fracice(i,ms))*&
                          exp(-gw_params%hkrz*(zimm(ms)/1000._r_2-gw_params%zdepth))

       ssnow%GWdhkdw(i)  = (1._r_2-ssnow%fracice(i,ms))* (2._r_2*soil%clappB(i,ms)+3._r_2)*&
                           s2(i)*0.5_r_2/(soil%GWwatsat(i)-soil%GWwatr(i))*&
                           exp(-gw_params%hkrz*(zimm(ms)/1000._r_2-gw_params%zdepth))

       ssnow%GWsmp(i)    = -soil%smpsat(i,ms)*s_mid(i)**(-soil%clappB(i,ms))
       ssnow%GWsmp(i)    = max(min(ssnow%GWsmp(i),-soil%GWsmpsat(i)),sucmin)
       ssnow%GWdsmpdw(i) = -soil%GWclappB(i)*ssnow%GWsmp(i)/(s_mid(i)*(soil%GWwatsat(i)-soil%GWwatr(i)))
    end do

END SUBROUTINE calc_soil_hydraulic_props


END MODULE cable_soil_snow_gw_module
