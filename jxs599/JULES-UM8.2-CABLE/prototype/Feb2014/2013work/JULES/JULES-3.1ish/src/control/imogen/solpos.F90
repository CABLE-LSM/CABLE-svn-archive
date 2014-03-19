#if !defined(UM_JULES)
!---------------------------------------------------------------------  
!  Unified model deck SOLPOS, containing only routine SOLPOS.         
!    This is part of logical component P233, performing the           
!  calculations of the earth's orbit described in the first page of   
!  the "Calculation of incoming insolation" section of UMDP 23, i.e.  
!  from the day of the year (and, in forecast mode, whether it is a   
!  leap year) and the orbital "constants" (which vary over            
!  "Milankovitch" timescales) it calculates the sin of the solar      
!  declination and the inverse-square scaling factor for the solar    
!  "constant".  It is thus intrinsically scalar.  The FORTRAN code    
!  present depends on whether *DEF CAL360 is set during UPDATE: this  
!  replaces the Julian calendar with the climate-mode 360-day calendar
!    Written in FORTRAN 77, with the addition of "!" comments and     
!  underscores in variable names.                                     
!    Written to comply with 12/9/89 version of UMDP 4 (meteorological 
!  standard).                                                         
!   Author:    William Ingram  22/3/89                                
!                      Reviewer: Clive Wilson Winter 1989/90          
!  First version.                                                     
!
  SUBROUTINE SOLPOS(DAY,YEAR,SINDEC,SCS)

    USE c_pi

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::                                        &
      DAY,                                                        &
             !  Day-number in the year
      YEAR   !  Calendar year

    REAL, INTENT(OUT) ::                                          &
      SINDEC,                                                     &
             !  sin(solar declination)
      SCS    !  solar constant scaling factor

! This routine has no dynamically allocated work areas and no          
!  significant structure.  It calls the intrinsic functions FLOAT, SIN 
!  & COS, but no user functions or subroutines.                        
!
    REAL ::                                                       &
! Basic orbital constants    
      GAMMA,                                                      &
             ! Gamma
      E,                                                          &
             ! e
      TAU0,                                                       &
             ! True date of perihelion
      SINOBL,                                                     &
             ! Sin (obliquity)
! Derived orbital constants
      TAU1,                                                       &
      E1,E2,                                                      &
             ! Coefficients for 3.1.2
      E3,                                                         &
      E4,                                                         &
             ! Constant for 3.1.4
      TWOPI  ! 2pi
      
    REAL :: DINY      ! Number of days in the year 
    REAL :: M,V       ! Mean & true anomaly

!-----------------------------------------------------------------
! Set up parameters
!-----------------------------------------------------------------
    PARAMETER( TWOPI = 2. * PI )
    PARAMETER(                                                    &
      GAMMA  = 1.352631,                                          &
      E      = .0167,                                             &
      TAU0   = 2.5,                                               &
      SINOBL = .397789                                            &
    )
    PARAMETER (                                                   &
      E1 = E * (2.-.25*E*E),                                      &
      E2 = 1.25 * E*E,                                            &
      E3 = E*E*E * 13./12.,                                       &
      E4 = ((1.+E*E*.5) / (1.-E*E))**2                            &
    )

!  In climate mode, DINY=360 always, and as well as applying 3.3.1,     
!  TAU1 is modified so as to include the conversion of day-ordinal into 
!  fractional-number-of-days-into-the-year-at-12-Z-on-this-day.         
    PARAMETER(                                                    &
      DINY = 360.,                                                &
      TAU1 = TAU0 * DINY / 365.25 + 0.71 + .5                     &
    )



    M = TWOPI * (FLOAT(DAY)-TAU1) / DINY                    ! Eq 3.1.1
    V = M + E1*SIN(M) + E2*SIN(2.*M) + E3*SIN(3.*M)         ! Eq 3.1.2
    SCS = E4 * ( 1. + E * COS(V) )**2                       ! Eq 3.1.4
    SINDEC = SINOBL * SIN (V - GAMMA)                       ! Eq 3.1.6

    RETURN                                                            

  END SUBROUTINE SOLPOS
#endif
