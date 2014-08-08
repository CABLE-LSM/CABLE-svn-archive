module cable_photosynthesis_module

contains

! ------------------------------------------------------------------------------

SUBROUTINE photosynthesis( csxz, cx1z, cx2z, gswminz,                          &
                           rdxz, vcmxt3z, vcmxt4z, vx3z,                       &
                           vx4z, xleuningz, vlaiz, deltlfz, anxz, fwsoilz )
   USE cable_def_types_mod, only : mp, mf, r_2
   
   REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: csxz
   
   REAL, DIMENSION(mp,mf), INTENT(IN) ::                                       &
      cx1z,       & !
      cx2z,       & !     
      gswminz,    & !
      rdxz,       & !
      vcmxt3z,    & !
      vcmxt4z,    & !
      vx4z,       & !
      vx3z,       & !
      xleuningz,  & !
      vlaiz,      & !
      deltlfz       ! 

   REAL, DIMENSION(mp,mf), INTENT(INOUT) :: anxz
   
   ! local variables
   REAL(r_2), DIMENSION(mp,mf) ::                                              &
      coef0z,coef1z,coef2z, ciz,delcxz,                                        &
      anrubiscoz,anrubpz,ansinkz

   REAL, DIMENSION(mp) :: fwsoilz  
 
   REAL, PARAMETER  :: effc4 = 4000.0  ! Vc=effc4*Ci*Vcmax (see
                                       ! Bonan,LSM version 1.0, p106)

   INTEGER :: i,j   
  
   
   DO i=1,mp
      
      IF (sum(vlaiz(i,:)) .GT. C%LAI_THRESH) THEN
      
         DO j=1,mf
            
            IF( vlaiz(i,j) .GT. C%LAI_THRESH .AND. deltlfz(i,j) .GT. 0.1) THEN

               ! Rubisco limited:
               coef2z(i,j) = gswminz(i,j)*fwsoilz(i) / C%RGSWC + xleuningz(i,j) * &
                             ( vcmxt3z(i,j) - ( rdxz(i,j)-vcmxt4z(i,j) ) )

               coef1z(i,j) = (1.0-csxz(i,j)*xleuningz(i,j)) *                  &
                             (vcmxt3z(i,j)+vcmxt4z(i,j)-rdxz(i,j))             &
                             + (gswminz(i,j)*fwsoilz(i)/C%RGSWC)*(cx1z(i,j)-csxz(i,j)) &
                             - xleuningz(i,j)*(vcmxt3z(i,j)*cx2z(i,j)/2.0      &
                             + cx1z(i,j)*(rdxz(i,j)-vcmxt4z(i,j) ) )
               
                
               coef0z(i,j) = -(1.0-csxz(i,j)*xleuningz(i,j)) *                 &    
                             (vcmxt3z(i,j)*cx2z(i,j)/2.0                       &
                             + cx1z(i,j)*( rdxz(i,j)-vcmxt4z(i,j ) ) )         &
                             -( gswminz(i,j)*fwsoilz(i)/C%RGSWC ) * cx1z(i,j)*csxz(i,j)

               ! kdcorbin,09/10 - new calculations
               IF( ABS(coef2z(i,j)) .GT. 1.0e-9 .AND. &
                   ABS(coef1z(i,j)) .LT. 1.0e-9) THEN
                  
                  ! no solution, give it a huge number as 
                  ! quadratic below cannot handle zero denominator
                  ciz(i,j) = 99999.0        
                  
                  anrubiscoz(i,j) = 99999.0 ! OR do ciz=0 and calc anrubiscoz
                   
               ENDIF

               ! solve linearly
               IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                   ABS( coef1z(i,j) ) >= 1e-9 ) THEN
                  
                  ! same reason as above
                  ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)  
          
                  ciz(i,j) = MAX( 0.0_r_2, ciz(i,j) )

                  anrubiscoz(i,j) = vcmxt3z(i,j)*(ciz(i,j)-cx2z(i,j) / 2.0 ) / &
                                    ( ciz(i,j) + cx1z(i,j)) + vcmxt4z(i,j) -   &
                                    rdxz(i,j)
   
               ENDIF
   
               ! solve quadratic (only take the more positive solution)
               IF( ABS( coef2z(i,j) ) >= 1.e-9 ) THEN
                  
                  delcxz(i,j) = coef1z(i,j)**2 -4.0 * coef0z(i,j)              &
                                * coef2z(i,j)
          
                  ciz(i,j) = ( -coef1z(i,j) + SQRT( MAX( 0.0_r_2 ,             &
                             delcxz(i,j) ) ) ) / ( 2.0*coef2z(i,j) )
                             
                  ciz(i,j) = MAX( 0.0_r_2, ciz(i,j) )   ! must be positive, why?
   
                  anrubiscoz(i,j) = vcmxt3z(i,j) * ( ciz(i,j) - cx2z(i,j)      &
                                    / 2.0)  / ( ciz(i,j) + cx1z(i,j) ) +       &
                                    vcmxt4z(i,j) - rdxz(i,j)
               
               ENDIF
   
               ! RuBP limited:
               coef2z(i,j) = gswminz(i,j)*fwsoilz(i) / C%RGSWC + xleuningz(i,j) &
                             * ( vx3z(i,j) - ( rdxz(i,j) - vx4z(i,j) ) )
   
               coef1z(i,j) = ( 1.0 - csxz(i,j) * xleuningz(i,j) ) *            &
                             ( vx3z(i,j) + vx4z(i,j) - rdxz(i,j) )             &
                             + ( gswminz(i,j)*fwsoilz(i) / C%RGSWC ) *          &
                             ( cx2z(i,j) - csxz(i,j) ) - xleuningz(i,j)        &
                             * ( vx3z(i,j) * cx2z(i,j) / 2.0 + cx2z(i,j) *     &
                             ( rdxz(i,j) - vx4z(i,j) ) )                          
                             
                             coef0z(i,j) = -(1.0-csxz(i,j)*xleuningz(i,j)) *   &
                             (vx3z(i,j)*cx2z(i,j)/2.0                          &
                             + cx2z(i,j)*(rdxz(i,j)-vx4z(i,j)))                &
                         - (gswminz(i,j)*fwsoilz(i)/C%RGSWC)*cx2z(i,j)*csxz(i,j)
   
   
               !kdcorbin, 09/10 - new calculations
               ! no solution, give it a huge number
               IF( ABS( coef2z(i,j) ) < 1.0e-9 .AND.                           &
                   ABS( coef1z(i,j) ) < 1.0e-9 ) THEN

                  ciz(i,j) = 99999.0
                  anrubpz(i,j)  = 99999.0

               ENDIF
   
               ! solve linearly
               IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                   ABS( coef1z(i,j) ) >= 1.e-9) THEN
                  
                  ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)
                  
                  ciz(i,j) = MAX(0.0_r_2,ciz(i,j))
                  
                  anrubpz(i,j) = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) /          &
                                 (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)
               
               ENDIF
   
               ! solve quadratic (only take the more positive solution)
               IF ( ABS( coef2z(i,j)) >= 1.e-9 ) THEN
                  
                  delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)
                  
                  ciz(i,j) = (-coef1z(i,j)+SQRT(MAX(0.0_r_2,delcxz(i,j))))     &
                             /(2.0*coef2z(i,j))
                   
                  ciz(i,j) = MAX(0.0_r_2,ciz(i,j)) 
                  
                  anrubpz(i,j)  = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) /         &
                                  (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)
   
               ENDIF
                 
               ! Sink limited:
               coef2z(i,j) = xleuningz(i,j)
               
               coef1z(i,j) = gswminz(i,j)*fwsoilz(i)/C%RGSWC + xleuningz(i,j)   &
                             * (rdxz(i,j) - 0.5*vcmxt3z(i,j))                  &
                             + effc4 * vcmxt4z(i,j) - xleuningz(i,j)           &
                             * csxz(i,j) * effc4 * vcmxt4z(i,j)  
                                            
               coef0z(i,j) = -( gswminz(i,j)*fwsoilz(i)/C%RGSWC )*csxz(i,j)*effc4 &
                             * vcmxt4z(i,j) + ( rdxz(i,j)                      &
                           - 0.5 * vcmxt3z(i,j)) * gswminz(i,j)*fwsoilz(i)/C%RGSWC
          
               ! no solution, give it a huge number
               IF( ABS( coef2z(i,j) ) < 1.0e-9 .AND.                           &
                   ABS( coef1z(i,j)) < 1.0e-9 ) THEN

                  ciz(i,j) = 99999.0
                  ansinkz(i,j)  = 99999.0

               ENDIF
          
               ! solve linearly
               IF( ABS( coef2z(i,j) ) < 1.e-9 .AND.                            &
                   ABS( coef1z(i,j) ) >= 1.e-9 ) THEN

                  ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)
                  ansinkz(i,j)  = ciz(i,j)

               ENDIF
               
               ! solve quadratic (only take the more positive solution)
               IF( ABS( coef2z(i,j) ) >= 1.e-9 ) THEN
                  
                  delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)
                  
                  ciz(i,j) = (-coef1z(i,j)+SQRT (MAX(0.0_r_2,delcxz(i,j)) ) )  &
                             / ( 2.0 * coef2z(i,j) )
   
                  ansinkz(i,j) = ciz(i,j)
   
               ENDIF
          
               ! minimal of three limited rates
               anxz(i,j) = MIN(anrubiscoz(i,j),anrubpz(i,j),ansinkz(i,j))
        
            ENDIF
      
         ENDDO
    
      ENDIF

   ENDDO
     
END SUBROUTINE photosynthesis



! ------------------------------------------------------------------------------

subroutine Temperature_dependence_Vcmax( tlfx, vcmax, &
               frac4, vcmxt3, scalex )
   real ::       
            ! Leuning 2002 (P C & E) equation for temperature response
            ! used for Vcmax for C3 plants:
            ! jhan: frac4 = frac of c4 plnts - but why isnt this 0 here 
            temp(i) =  xvcmxt3(tlfx(i)) * veg%vcmax(i) * (1.0-veg%frac4(i))
            
            ! jhan: 1$2<-mf=sunlkit,shaded. scalex = scaling param from ?  
            vcmxt3(1) = rad%scalex(i,1) * temp(i)
            vcmxt3(2) = rad%scalex(i,2) * temp(i)
 
End subroutine Temperature_dependence_Vcmax 

! ------------------------------------------------------------------------------

FUNCTION xvcmxt3(x) RESULT(z)
   
   !  leuning 2002 (p c & e) equation for temperature response
   !  used for vcmax for c3 plants
   REAL, INTENT(IN) :: x
   REAL :: xvcnum,xvcden,z
 
   REAL, PARAMETER  :: EHaVc  = 73637.0  ! J/mol (Leuning 2002)
   REAL, PARAMETER  :: EHdVc  = 149252.0 ! J/mol (Leuning 2002)
   REAL, PARAMETER  :: EntropVc = 486.0  ! J/mol/K (Leuning 2002)
   REAL, PARAMETER  :: xVccoef = 1.17461 ! derived parameter
                     ! xVccoef=1.0+exp((EntropJx*C%TREFK-EHdJx)/(Rconst*C%TREFK))
 
   xvcnum=xvccoef*exp( ( ehavc / ( C%rgas*C%TREFK ) )* ( 1.-C%TREFK/x ) )
   xvcden=1.0+exp( ( entropvc*x-ehdvc ) / ( C%rgas*x ) )
   z = max( 0.0,xvcnum / xvcden )

END FUNCTION xvcmxt3

! ------------------------------------------------------------------------------

! Explicit array dimensions as temporary work around for NEC inlining problem
FUNCTION xvcmxt4(x) RESULT(z)
   
   REAL, PARAMETER      :: q10c4 = 2.0
   REAL, INTENT(IN) :: x
   REAL :: z
 
   z = q10c4 ** (0.1 * x - 2.5) /                                              &
        ((1.0 + exp(0.3 * (13.0 - x))) * (1.0 + exp(0.3 * (x - 36.0))))
 
END FUNCTION xvcmxt4

! ------------------------------------------------------------------------------

FUNCTION xejmxt3(x) RESULT(z)
   
   !  leuning 2002 (p c & e) equation for temperature response
   !  used for jmax for c3 plants
 
   REAL, INTENT(IN) :: x
   REAL :: xjxnum,xjxden,z   
 
   REAL, PARAMETER  :: EHaJx  = 50300.0  ! J/mol (Leuning 2002)
   REAL, PARAMETER  :: EHdJx  = 152044.0 ! J/mol (Leuning 2002)
   REAL, PARAMETER  :: EntropJx = 495.0  ! J/mol/K (Leuning 2002)
   REAL, PARAMETER  :: xjxcoef = 1.16715 ! derived parameter
 
   xjxnum = xjxcoef*exp( ( ehajx / ( C%rgas*C%TREFK ) ) * ( 1.-C%TREFK / x ) )
   xjxden=1.0+exp( ( entropjx*x-ehdjx) / ( C%rgas*x ) )
   z = max(0.0, xjxnum/xjxden)

END FUNCTION xejmxt3

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

end module cable_photosynthesis_module
