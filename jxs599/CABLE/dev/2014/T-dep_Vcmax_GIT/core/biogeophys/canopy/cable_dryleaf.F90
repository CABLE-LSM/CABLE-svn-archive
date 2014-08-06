MODULE cable_canopy_dryleaf_mod

contains
   
SUBROUTINE dryLeaf( dels, rad, rough, air, met,                                &
                    veg, canopy, soil, ssnow, dsx,                             &
                    fwsoil, tlfx,  tlfy,  ecy, hcy,                            &
                    rny, gbhu, gbhf, csx,                                      &
                    cansat, ghwet, iter )

   USE cable_def_types_mod
   USE cable_common_module
   use cable_photosynthesis_mod

   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (air_type),       INTENT(INOUT) :: air
   TYPE (met_type),       INTENT(INOUT) :: met
   TYPE (canopy_type),    INTENT(INOUT) :: canopy
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow

   TYPE (veg_parameter_type),  INTENT(INOUT)   :: veg
   TYPE (soil_parameter_type), INTENT(inout)   :: soil
   
   REAL, INTENT(INOUT), DIMENSION(:) ::                                        &
      dsx,        & ! leaf surface vpd
      fwsoil,     & ! soil water modifier of stom. cond
      tlfx,       & ! leaf temp prev. iter (K)
      tlfy          ! leaf temp (K)
   
   REAL(R_2),INTENT(INOUT), DIMENSION(:) ::                                    &
      ecy,        & ! lat heat fl dry big leaf
      hcy,        & ! veg. sens heat
      rny         !& !

   REAL(R_2),INTENT(INOUT), DIMENSION(:,:) ::                                  &
      gbhu,       & ! forcedConvectionBndryLayerCond
      gbhf,       & ! freeConvectionBndryLayerCond
      csx           ! leaf surface CO2 concentration

   REAL,INTENT(IN), DIMENSION(:) :: cansat

   REAL(r_2), INTENT(OUT), DIMENSION(:) ::                                     &
      ghwet  ! cond for heat for a wet canopy

   INTEGER,INTENT(IN) :: iter
 
   REAL, INTENT(IN)     :: dels ! integration time step (s)
   
   !local variables
   REAL, PARAMETER  ::                                                         &
      co2cp3 = 0.0,  & ! CO2 compensation pt C3
      jtomol = 4.6e-6  ! Convert from J to Mol for light

   REAL, DIMENSION(mp) ::                                                      &
      conkct,        & ! Michaelis Menton const.
      conkot,        & ! Michaelis Menton const.
      cx1,           & ! "d_{3}" in Wang and Leuning,
      cx2,           & !     1998, appendix E
      tdiff,         & ! leaf air temp diff.
      tlfxx,         & ! leaf temp of current iteration (K)
      abs_deltlf,    & ! ABS(deltlf)
      deltlf,        & ! deltlfy of prev iter.
      deltlfy,       & ! del temp successive iter.
      gras,          & ! Grashof coeff
      evapfb,        & !
      sum_rad_rniso, & !
      sum_rad_gradis,& ! 
      gwwet,         & ! cond for water for a wet canopy
      ghrwet,        & ! wet canopy cond: heat & thermal rad
      sum_gbh,       & !
      ccfevw,        & ! limitation term for
                       ! wet canopy evaporation rate
      temp             !

   REAL(r_2), DIMENSION(mp)  ::                                                &
      ecx,        & ! lat. hflux big leaf
      ecx_t,      & ! lat. hflux big leaf
      hcx,        & ! sens heat fl big leaf prev iteration
      rnx,        & ! net rad prev timestep
      fwsoil_coef   !

   REAL, DIMENSION(mp,ms)  :: oldevapfbl
      
   REAL, DIMENSION(mp,mf)  ::                                                  &
      gw,         & ! cond for water for a dry canopy
      gh,         & ! cond for heat for a dry canopy
      ghr,        & ! dry canopy cond for heat & thermal rad
      anx,        & ! net photos. prev iteration
      an_y,       & ! net photosynthesis soln
      rdx,        & ! daytime leaf resp rate, prev iteration
      rdy,        & ! daytime leaf resp rate
      ejmax2,     & ! jmax of big leaf
      ejmxt3,     & ! jmax big leaf C3 plants
      vcmxt3,     & ! vcmax big leaf C3
      vcmxt4,     & ! vcmax big leaf C4
      vx3,        & ! carboxylation C3 plants
      vx4,        & ! carboxylation C4 plants
      xleuning,   & ! leuning stomatal coeff
      psycst,     & ! modified pych. constant
      frac42,     & ! 2D frac4
      temp2

   REAL, DIMENSION(:,:), POINTER :: gswmin ! min stomatal conductance
   
   REAL, DIMENSION(mp,2) ::  gsw_term, lower_limit2  ! local temp var 

   INTEGER :: i, j, k, kk  ! iteration count
   
   ! END header


   ALLOCATE( gswmin(mp,mf ))

   ! Soil water limitation on stomatal conductance:
   IF( iter ==1) THEN
   
      IF(cable_user%FWSOIL_SWITCH == 'standard') THEN
         CALL fwsoil_calc_std( fwsoil, soil, ssnow, veg) 
      ELSEIf (cable_user%FWSOIL_SWITCH == 'non-linear extrapolation') THEN
         !EAK, 09/10 - replace linear approx by polynomial fitting
         CALL fwsoil_calc_non_linear(fwsoil, soil, ssnow, veg) 
      ELSEIF(cable_user%FWSOIL_SWITCH == 'Lai and Ktaul 2000') THEN
         CALL fwsoil_calc_Lai_Ktaul(fwsoil, soil, ssnow, veg) 
      ELSE
         STOP 'fwsoil_switch failed.'
      ENDIF

   ENDIF

   ! weight min stomatal conductance by C3 an C4 plant fractions
   frac42 = SPREAD(veg%frac4, 2, mf) ! frac C4 plants

   gsw_term = C%gsw03 * (1. - frac42) + C%gsw04 * frac42
   lower_limit2 = rad%scalex * (C%gsw03 * (1. - frac42) + C%gsw04 * frac42)
   gswmin = max(1.e-6,lower_limit2)
         

   gw = 1.0e-3 ! default values of conductance
   gh = 1.0e-3
   ghr= 1.0e-3
   rdx = 0.0
   anx = 0.0
   rnx = SUM(rad%rniso,2)
   abs_deltlf = 999.0


   gras = 1.0e-6
   an_y= 0.0
   hcx = 0.0              ! init sens heat iteration memory variable
   hcy = 0.0
   rdy = 0.0
   ecx = SUM(rad%rniso,2) ! init lat heat iteration memory variable
   tlfxx = tlfx
   psycst(:,:) = SPREAD(air%psyc,2,mf)
   canopy%fevc = 0.0
   ssnow%evapfbl = 0.0

   ghwet = 1.0e-3
   gwwet = 1.0e-3
   ghrwet= 1.0e-3
   canopy%fevw = 0.0
   canopy%fhvw = 0.0
   sum_gbh = SUM((gbhu+gbhf),2)
   sum_rad_rniso = SUM(rad%rniso,2)
   sum_rad_gradis = SUM(rad%gradis,2)

   DO kk=1,mp

      IF(canopy%vlaiw(kk) <= C%LAI_THRESH) THEN
         rnx(kk) = 0.0 ! intialise
         ecx(kk) = 0.0 ! intialise
         ecy(kk) = ecx(kk) ! store initial values
         abs_deltlf(kk)=0.0
         rny(kk) = rnx(kk) ! store initial values
         ! calculate total thermal resistance, rthv in s/m
      END IF

   ENDDO
   
   deltlfy = abs_deltlf
   k = 0

   !kdcorbin, 08/10 - doing all points all the time
   DO WHILE (k < C%MAXITER)
      k = k + 1

      DO i=1,mp
         
         IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. abs_deltlf(i) > 0.1) THEN
         
            ghwet(i) = 2.0   * sum_gbh(i)
            gwwet(i) = 1.075 * sum_gbh(i)
            ghrwet(i) = sum_rad_gradis(i) + ghwet(i)
       
            ! Calculate fraction of canopy which is wet:
            canopy%fwet(i) = MAX( 0.0, MIN( 1.0, 0.8 * canopy%cansto(i)/ MAX(  &
                             cansat(i),0.01 ) ) )

            ! Calculate lat heat from wet canopy, may be neg.                  
            ! if dew on wet canopy to avoid excessive evaporation:
            ccfevw(i) = MIN(canopy%cansto(i) * air%rlam(i) / dels,             &
                        2.0 / (1440.0 / (dels/60.0)) * air%rlam(i) )
   
            ! Grashof number (Leuning et al, 1995) eq E4:
            gras(i) = MAX(1.0e-6,                                              &
                      1.595E8* ABS( tlfx(i)-met%tvair(i))* (veg%dleaf(i)**3.0) )

            ! See Appendix E in (Leuning et al, 1995):
            gbhf(i,1) = rad%fvlai(i,1) * air%cmolar(i) * 0.5*C%dheat           &
                       * ( gras(i)**0.25 ) / veg%dleaf(i)
            gbhf(i,2) = rad%fvlai(i,2) * air%cmolar(i) * 0.5 * C%dheat         &
                        * ( gras(i)**0.25 ) / veg%dleaf(i)
            gbhf(i,:) = MAX( 1.e-6, gbhf(i,:) )
      
            ! Conductance for heat:
            gh(i,:) = 2.0 * (gbhu(i,:) + gbhf(i,:))
      
            ! Conductance for heat and longwave radiation:
            ghr(i,:) = rad%gradis(i,:)+gh(i,:)
!Vcmax      
            call Temperature_dependence_Vcmax( tlfx(i), veg%vcmax(i), &
                     veg%frac4(i), vcmxt3(i,:), rad%scalex(i,:) )
   
            ! Temperature response Vcmax, C4 plants (Collatz et al 1989):
            temp(i) = xvcmxt4(tlfx(i)-C%tfrz) * veg%vcmax(i) * veg%frac4(i)
            vcmxt4(i,1) = rad%scalex(i,1) * temp(i)
            vcmxt4(i,2) = rad%scalex(i,2) * temp(i)
    
!Jmax
            ! Leuning 2002 (P C & E) equation for temperature response
            ! used for Jmax for C3 plants:
            temp(i) = xejmxt3(tlfx(i)) * veg%ejmax(i) * (1.0-veg%frac4(i))
            ejmxt3(i,1) = rad%scalex(i,1) * temp(i)
            ejmxt3(i,2) = rad%scalex(i,2) * temp(i)
!!!!            
            ! Difference between leaf temperature and reference temperature:
            tdiff(i) = tlfx(i) - C%TREFK
            
            ! Michaelis menten constant of Rubisco for CO2:
            conkct(i) = C%conkc0 * EXP( (C%ekc / ( C%rgas*C%trefk) ) *         &
                        ( 1.0 - C%trefk/tlfx(i) ) )

            ! Michaelis menten constant of Rubisco for oxygen:
            conkot(i) = C%conko0 * EXP( ( C%eko / (C%rgas*C%trefk) ) *         &
                        ( 1.0 - C%trefk/tlfx(i) ) )
   
            ! Store leaf temperature
            tlfxx(i) = tlfx(i)
   
            ! "d_{3}" in Wang and Leuning, 1998, appendix E:
            cx1(i) = conkct(i) * (1.0+0.21/conkot(i))
            cx2(i) = 2.0 * C%gam0 * ( 1.0 + C%gam1 * tdiff(i) +                    &
                     C%gam2 * tdiff(i) * tdiff(i ))
    
            ! All equations below in appendix E in Wang and Leuning 1998 are
            ! for calculating anx, csx and gswx for Rubisco limited,
            ! RuBP limited, sink limited
            temp2(i,1) = rad%qcan(i,1,1) * jtomol * (1.0-veg%frac4(i))
            temp2(i,2) = rad%qcan(i,2,1) * jtomol * (1.0-veg%frac4(i))
            vx3(i,1)  = ej3x(temp2(i,1),ejmxt3(i,1))
            vx3(i,2)  = ej3x(temp2(i,2),ejmxt3(i,2))
    
            temp2(i,1) = rad%qcan(i,1,1) * jtomol * veg%frac4(i)
            temp2(i,2) = rad%qcan(i,2,1) * jtomol * veg%frac4(i)
            vx4(i,1)  = ej4x(temp2(i,1),vcmxt4(i,1))
            vx4(i,2)  = ej4x(temp2(i,2),vcmxt4(i,2))
    
            rdx(i,1) = (C%cfrd3*vcmxt3(i,1) + C%cfrd4*vcmxt4(i,1))
            rdx(i,2) = (C%cfrd3*vcmxt3(i,2) + C%cfrd4*vcmxt4(i,2))
            
            xleuning(i,1) = ( fwsoil(i) / ( csx(i,1) - co2cp3 ) )              &
                          * ( ( 1.0 - veg%frac4(i) ) * C%A1C3 / ( 1.0 + dsx(i) &
                          / C%d0c3 ) + veg%frac4(i)    * C%A1C4 / (1.0+dsx(i)/ &
                          C%d0c4) )

            xleuning(i,2) = ( fwsoil(i) / ( csx(i,2)-co2cp3 ) )                &
                            * ( (1.0-veg%frac4(i) ) * C%A1C3 / ( 1.0 + dsx(i) /&
                            C%d0c3 ) + veg%frac4(i)    * C%A1C4 / (1.0+ dsx(i)/&
                            C%d0c4) )
    
         
         ENDIF
         
      ENDDO !i=1,mp
   
      CALL photosynthesis( csx(:,:),                                           &
                           SPREAD( cx1(:), 2, mf ),                            &
                           SPREAD( cx2(:), 2, mf ),                            &
                           gswmin(:,:), rdx(:,:), vcmxt3(:,:),                 &
                           vcmxt4(:,:), vx3(:,:), vx4(:,:),                    &
                           xleuning(:,:), rad%fvlai(:,:),                      &
                           SPREAD( abs_deltlf, 2, mf ),                        &
                           anx(:,:), fwsoil(:) )

      DO i=1,mp
         
         IF (canopy%vlaiw(i) > C%LAI_THRESH .AND. abs_deltlf(i) > 0.1) Then
      
            DO kk=1,mf
               
               IF(rad%fvlai(i,kk)>C%LAI_THRESH) THEN

                  csx(i,kk) = met%ca(i) - C%RGBWC*anx(i,kk) / (                &
                              gbhu(i,kk) + gbhf(i,kk) )
                  csx(i,kk) = MAX( 1.0e-4, csx(i,kk) )

                  canopy%gswx(i,kk) = MAX( 1.e-3, gswmin(i,kk)*fwsoil(i) +     &
                                      MAX( 0.0, C%RGSWC * xleuning(i,kk) *     &
                                      anx(i,kk) ) )

                  !Recalculate conductance for water:
                  gw(i,kk) = 1.0 / ( 1.0 / canopy%gswx(i,kk) +                 &
                             1.0 / ( 1.075 * ( gbhu(i,kk) + gbhf(i,kk) ) ) )

                  gw(i,kk) = MAX( gw(i,kk), 0.00001 )

                  ! Modified psychrometric constant 
                  ! (Monteith and Unsworth, 1990)
                  psycst(i,kk) = air%psyc(i) * REAL( ghr(i,kk) / gw(i,kk) )
           
               ENDIF
            
            ENDDO

            ecx(i) = ( air%dsatdk(i) * ( rad%rniso(i,1) - C%capp * C%rmair     &
                     * ( met%tvair(i) - met%tk(i) ) * rad%gradis(i,1) )        &
                     + C%capp * C%rmair * met%dva(i) * ghr(i,1) )              &
                     / ( air%dsatdk(i) + psycst(i,1) ) + ( air%dsatdk(i)       &
                     * ( rad%rniso(i,2) - C%capp * C%rmair * ( met%tvair(i) -  &
                     met%tk(i) ) * rad%gradis(i,2) ) + C%capp * C%rmair *      &
                     met%dva(i) * ghr(i,2) ) /                                 &
                     ( air%dsatdk(i) + psycst(i,2) ) 


            IF (ecx(i) > 0.0 .AND. canopy%fwet(i) < 1.0) Then
               evapfb(i) = ( 1.0 - canopy%fwet(i)) * REAL( ecx(i) ) *dels      &
                           / air%rlam(i)

               DO kk = 1,ms
                  
                  ssnow%evapfbl(i,kk) = MIN( evapfb(i) * veg%froot(i,kk),      &
                                        MAX( 0.0, REAL( ssnow%wb(i,kk) ) -     &
                                        1.1 * soil%swilt(i) ) *                &
                                        soil%zse(kk) * 1000.0 )

               ENDDO

               canopy%fevc(i) = SUM(ssnow%evapfbl(i,:))*air%rlam(i)/dels
    
               ecx(i) = canopy%fevc(i) / (1.0-canopy%fwet(i))

            ENDIF

            ! Update canopy sensible heat flux:
            hcx(i) = (SUM(rad%rniso(i,:))-ecx(i)                               &
               - C%capp*C%rmair*(met%tvair(i)-met%tk(i))                       &
               * SUM(rad%gradis(i,:)))                                         &
               * SUM(gh(i,:))/ SUM(ghr(i,:))

            ! Update leaf temperature:
            tlfx(i)=met%tvair(i)+REAL(hcx(i))/(C%capp*C%rmair*SUM(gh(i,:)))
      
            ! Update net radiation for canopy:
            rnx(i) = SUM( rad%rniso(i,:)) -                                    &
                     C%CAPP * C%rmair *( tlfx(i)-met%tk(i) ) *                 &
                     SUM( rad%gradis(i,:) )

            ! Update leaf surface vapour pressure deficit:
            dsx(i) = met%dva(i) + air%dsatdk(i) * (tlfx(i)-met%tvair(i))

            ! Store change in leaf temperature between successive iterations:
            deltlf(i) = tlfxx(i)-tlfx(i)
            abs_deltlf(i) = ABS(deltlf(i))

         ENDIF !lai/abs_deltlf

      ENDDO !i=1,mp

      ! Whhere leaf temp change b/w iterations is significant, and
      ! difference is smaller than the previous iteration, store results:
      DO i=1,mp
      
         IF ( abs_deltlf(i) < ABS( deltlfy(i) ) ) THEN

            deltlfy(i) = deltlf(i)
            tlfy(i) = tlfx(i)
            rny(i) = rnx(i)
            hcy(i) = hcx(i)
            ecy(i) = ecx(i)
            rdy(i,1) = rdx(i,1)
            rdy(i,2) = rdx(i,2)
            an_y(i,1) = anx(i,1)
            an_y(i,2) = anx(i,2)
            
            ! save last values calculated for ssnow%evapfbl
            oldevapfbl(i,1) = ssnow%evapfbl(i,1)
            oldevapfbl(i,2) = ssnow%evapfbl(i,2)
            oldevapfbl(i,3) = ssnow%evapfbl(i,3)
            oldevapfbl(i,4) = ssnow%evapfbl(i,4)
            oldevapfbl(i,5) = ssnow%evapfbl(i,5)
            oldevapfbl(i,6) = ssnow%evapfbl(i,6)

         ENDIF
          
         IF( abs_deltlf(i) > 0.1 )                                             &
            
            ! after 4 iterations, take mean of current & previous estimates
            ! as the next estimate of leaf temperature, to avoid oscillation
            tlfx(i) = ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) *tlfxx(i) + &
                      ( 1.0 - ( 0.5 * ( MAX( 0, k-5 ) / ( k - 4.9999 ) ) ) )   &
                      * tlfx(i)
     
         IF(k==1) THEN
         
            ! take the first iterated estimates as the defaults
            tlfy(i) = tlfx(i)
            rny(i) = rnx(i)
            hcy(i) = hcx(i)
            ecy(i) = ecx(i)
            rdy(i,:) = rdx(i,:)
            an_y(i,:) = anx(i,:)
            ! save last values calculated for ssnow%evapfbl
            oldevapfbl(i,:) = ssnow%evapfbl(i,:)
         
         END IF
      
      END DO !over mp 

   END DO  ! DO WHILE (ANY(abs_deltlf > 0.1) .AND.  k < C%MAXITER)


   ! dry canopy flux
   canopy%fevc = (1.0-canopy%fwet) * ecy

   ! Recalculate ssnow%evapfbl as ecy may not be updated with the ecx
   ! calculated in the last iteration.
   ! DO NOT use simple scaling as there are times that ssnow%evapfbl is zero.
   ! ** ssnow%evapfbl(i,:) = ssnow%evapfbl(i,:) * ecy(i) / ecx(i) **
   DO i = 1, mp
      
       IF( ecy(i) > 0.0 .AND. canopy%fwet(i) < 1.0 ) THEN
         
         IF( ABS( ecy(i) - ecx(i) ) > 1.0e-6 ) THEN
            
            IF( ABS( canopy%fevc(i) - ( SUM( oldevapfbl(i,:)) * air%rlam(i)    &
                /dels ) ) > 1.0e-4 ) THEN
               
               PRINT *, 'Error! oldevapfbl not right.', ktau_gl, i
               PRINT *, 'ecx, ecy = ', ecx(i), ecy(i)
               PRINT *, 'or in mm = ', ecx(i) * ( 1.0 - canopy%fwet(i) )       &
                                       / air%rlam(i) * dels,                   &
                                       ecy(i) * ( 1.0 - canopy%fwet(i) ) /     &
                                       air%rlam(i) * dels

               PRINT *,'fevc = ', canopy%fevc(i), SUM( oldevapfbl(i,:) ) *     &
                                  air%rlam(i) / dels
               PRINT *, 'fwet = ', canopy%fwet(i)
               PRINT *, 'oldevapfbl = ', oldevapfbl(i,:)
               PRINT *, 'ssnow%evapfbl before rescaling: ',                    &
                                                           ssnow%evapfbl(i,:)
               STOP
            
            ELSE
            
               ssnow%evapfbl(i,:) = oldevapfbl(i,:)
            
            END IF
         
         END IF

      END IF
   
   END DO

   canopy%frday = 12.0 * SUM(rdy, 2)
   canopy%fpn = -12.0 * SUM(an_y, 2)
   canopy%evapfbl = ssnow%evapfbl
   
   DEALLOCATE( gswmin )

END SUBROUTINE dryLeaf

! -----------------------------------------------------------------------------
   
END MODULE cable_canopy_dryleaf_mod
