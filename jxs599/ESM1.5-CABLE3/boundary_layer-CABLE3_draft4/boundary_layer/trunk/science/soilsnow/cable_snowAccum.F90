MODULE snow_accum_mod

USE cable_soil_snow_data_mod

CONTAINS 


SUBROUTINE snow_accum ( dels,  canopy, met, ssnow, soil )
USE cable_common_module
IMPLICIT NONE

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
      
      WHERE( canopy%precis > 0.0 .AND. ssnow%tgg(:,1) < CTFRZ )
         
         ssnow%snowd = MAX(ssnow%snowd + canopy%precis, 0.0)
        
         ssnow%tgg(:,1) = ssnow%tgg(:,1) + canopy%precis * CHLF               &
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
         
         ssnow%tggsn(:,1) = ssnow%tggsn(:,1) + canopy%precis * CHLF           &
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
         ssnow%tggsn(:,2) = ssnow%tggsn(:,2) + canopy%precis * CHLF           &
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
         ssnow%tggsn(:,3) = ssnow%tggsn(:,3) + canopy%precis * CHLF           &
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
   canopy%segg = canopy%fess / CHL
   canopy%segg = ( canopy%fess + canopy%fes_cor ) / CHL

   ! Initialise snow evaporation:
   ssnow%evapsn = 0
    
   ! Snow evaporation and dew on snow
   WHERE( ssnow%snowd > 0.1 )
      
      ssnow%evapsn = dels * ( canopy%fess + canopy%fes_cor ) / ( CHL + CHLF )
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


END MODULE snow_accum_mod
