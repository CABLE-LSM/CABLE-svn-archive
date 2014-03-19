#if !defined(UM_JULES)
! *****************************COPYRIGHT******************************* 
! (c) CROWN COPYRIGHT 2000, Met Office, All Rights Reserved.            
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details   
! *****************************COPYRIGHT******************************* 
! SUBROUTINE FREEZE_SOIL -----------------------------------------   
                                                                        
!                                                                   
! Subroutine Interface:                                                 
  SUBROUTINE FREEZE_SOIL (NPNTS,NSHYD,B,DZ,SATHH,SMCL,TSOIL,V_SAT,STHU,STHF)         

                                                                        
  USE c_0_dg_c, ONLY : zeroDegC
  USE c_perma
  USE c_densty

  IMPLICIT NONE                                                     
!                                                                   
! Description:                                                          
! Calculates the unfrozen and frozen water within a soil layer      
! as a fraction of saturation.                          (Cox, 6/95) 
!                                                                   
! Documentation : UM Documentation Paper 25                             
!                                                                   
! Current Code Owner : David Gregory                                    
!                                                                   
! History:                                                              
! Version   Date     Comment                                            
! -------   ----     -------                                            
!  3.4      6/95     Original code    Peter Cox                         
!  4.4  18/09/97     Rename from FREEZE to FREEZE_SOIL. D. Robinson.    
!                                                                   
! Code Description:                                                     
!   Language: FORTRAN 77 + common extensions.                           
!                                                                   
! System component covered: P25                                         
! System Task: P25                                                      
!                                                                   
                                                                    
                                                                    
! Subroutine arguments                                                  
! Scalar arguments with intent(IN) :                                  
  INTEGER  ::  &
     & NPNTS               &! IN Number of gridpoints.                  
     &,NSHYD                ! IN Number of soil layers.                 
                                                                        
                                                                    
!   Array arguments with intent(IN) :                                   
                                                                        
  REAL ::  &
     & B(NPNTS,NSHYD)      &! IN Clapp-Hornberger exponent.             
     &,DZ(NSHYD)           &! IN Thicknesses of the soil layers (m).
     &,SATHH(NPNTS,NSHYD)  &! IN Saturated soil water pressure (m).     
     &,SMCL(NPNTS,NSHYD)   &! IN Soil moisture content of               
                            !    layers (kg/m2).                        
     &,TSOIL(NPNTS,NSHYD)  &! IN Sub-surface temperatures (K).          
     &,V_SAT(NPNTS,NSHYD)   ! IN Volumetric soil moisture               
                            !    concentration at saturation            
                            !    (m3 H2O/m3 soil).                      
                                                                        
!   Array arguments with intent(OUT) :                                  
      REAL ::  &
     & STHF(NPNTS,NSHYD)   &! OUT Frozen soil moisture content of       
                            !     the layers as a fraction of           
                            !     saturation.                           
     &,STHU(NPNTS,NSHYD)    ! OUT Unfrozen soil moisture content of     
                            !     the layers as a fraction of           
                            !     saturation.                           
                                                                        
! Local scalars:                                                        
  INTEGER                                                           &
     & I,N                  ! WORK Loop counters.

  REAL                                                              &
     small_value,tiny_0     ! WORK
                                                                        
! Local arrays:                                                         
  REAL                                                              &
     & SMCLF(NPNTS,NSHYD)  &! WORK Frozen moisture content of the       
                            !      soil layers (kg/m2).                 
     &,SMCLU(NPNTS,NSHYD)  &! WORK Unfrozen moisture content of the     
                            !      soil layers (kg/m2).                 
     &,SMCLSAT(NPNTS,NSHYD)&! WORK The saturation moisture content of   
                            !      the layers (kg/m2).                  
     &,TMAX(NPNTS)         &! WORK Temperature above which all water is 
                            !      unfrozen (Celsius)                   
     &,TSL(NPNTS,NSHYD)    &! WORK Soil layer temperatures (Celsius).   
     &,WORK1(NPNTS)         ! WORK

!-----------------------------------------------------------------------
  tiny_0=TINY(0.0)
  small_value=EPSILON(0.0)

  DO N=1,NSHYD                                                      
                                                                    
    DO I=1,NPNTS                                                    
!-----------------------------------------------------------------------
! Calculate TMAX, the temperature above which all soil water is         
! unfrozen                                                              
!-----------------------------------------------------------------------
      SMCLSAT(I,N)=RHO_WATER*DZ(N)*V_SAT(I,N)
      TSL(I,N)=TSOIL(I,N)-ZERODEGC                                  
      IF ( (V_SAT(I,N) > tiny_0) .and. (SMCL(I,N) > small_value)) THEN
        WORK1(I)=(SMCL(I,N)/SMCLSAT(I,N))**(B(I,N))
        IF ( WORK1(I) > small_value ) THEN
          TMAX(I)=-SATHH(I,N)/(DPSIDT*WORK1(I))
          TMAX(I)=MAX(TMAX(I),-ZERODEGC)                                           
        ELSE
          TMAX(I)=-ZERODEGC
        ENDIF  
      ELSE
        TMAX(I)=-ZERODEGC
      ENDIF                                                        
             

!--------------------------------------------------------------------   
! Diagnose unfrozen and frozen water contents                           
!--------------------------------------------------------------------   
      IF (TSL(I,N).GE.TMAX(I)) THEN                                 

        SMCLU(I,N)=SMCL(I,N)                                        
        SMCLF(I,N)=0.0              
                                
      ELSE                                                          

!-----------------------------------------------------------------      
! For ice points (V_SAT=0) set SMCLU=0.0 and SMCLF=0.0                  
!-----------------------------------------------------------------      
        IF (V_SAT(I,N).EQ.0.0) THEN                                 
          SMCLU(I,N)=0.0                                            
          SMCLF(I,N)=0.0                                            
        ELSE                                                        
          SMCLU(I,N)=SMCLSAT(I,N)*(-DPSIDT*TSL(I,N)/SATHH(I,N))**(-1.0/B(I,N)) 
          SMCLF(I,N)=SMCL(I,N)-SMCLU(I,N)                           
        ENDIF                                                       
      ENDIF                                                         

      IF (SMCLSAT(I,N).GT.0.0) THEN                                 
        STHF(I,N)=SMCLF(I,N)/SMCLSAT(I,N)                           
        STHU(I,N)=SMCLU(I,N)/SMCLSAT(I,N)                           
      ELSE                                                          
        STHF(I,N)=0.0                                               
        STHU(I,N)=0.0                                               
      ENDIF                                                         

    ENDDO  !  i (points)                                                           
                                                                    
  ENDDO  !  n (layers)                                                           
                                                                    
  RETURN                                                            
  END SUBROUTINE freeze_soil
#endif
