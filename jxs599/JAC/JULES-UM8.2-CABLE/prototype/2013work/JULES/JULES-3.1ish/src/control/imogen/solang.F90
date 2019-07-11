#if !defined(UM_JULES)
!  Unified model deck SOLANG, containing only routine SOLANG.         
!    This is part of logical component P233, performing the           
!  calculations of the earth's orbit described in the second page of  
!  the "Calculation of incoming insolation" section of UMDP 23, i.e.  
!  from the sin of the solar  declination, the position of each point 
!  and the time limits it calculates how much sunlight, if any, it    
!  receives.                                                          
!    Written in FORTRAN 77 with the addition of "!" comments and      
!  underscores in variable names.                                     
!    Written to comply with 12/9/89 version of UMDP 4 (meteorological 
!  standard).                                                         
!    Author:    William Ingram  21/3/89                               
!                      Reviewer: Clive Wilson Winter 1989/90          
!  First version.                                                     


!*********************************************************************  
! A01_2A selected for SCM use equivalent to frozen version              
! physics 1/3/93 J.Lean                                                 
!*********************************************************************  
  SUBROUTINE SOLANG(SINDEC,T,DT,SINLAT,LONGIT,K,LIT,COSZ)

    USE c_pi

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::                                        &
      K     ! Number of points

    REAL, INTENT(IN) ::                                           &
      SINDEC,                                                     &
            ! Sin(solar declination)
      T,DT,                                                       &
            ! Start time (GMT) & timestep
      SINLAT(K),LONGIT(K)
            ! sin(latitude) & longitude of each point

    REAL, INTENT(OUT) ::                                          &
      LIT(K),                                                     &
            ! Sunlit fraction of the timestep
      COSZ(K)
            ! Mean cos(solar zenith angle) during the sunlit
            ! fraction

! This routine has no dynamically allocated work areas.  It calls the  
! intrinsic functions SQRT, ACOS & SIN, but no user functions or       
! subroutines.  The only structure is a loop over all the points to be 
! dealt with, with IF blocks nested inside to cover the various        
! possibilities.                                                       

    INTEGER :: J  ! Loop counter over points

    REAL ::                                                       &
      TWOPI,                                                      &
            ! 2*pi
      S2R,                                                        &
            ! Seconds-to-radians converter
      SINSIN,COSCOS,                                              &
            ! Products of the sines and of the cosines
            ! of solar declination and of latitude.
      HLD,COSHLD,                                                 &
            ! Half-length of the day in radians (equal to the
            ! hour-angle of sunset, and minus the hour-angle of
            ! sunrise) & its cosine.
      HAT,                                                        &
            ! Local hour angle at the start time.
      OMEGAB,OMEGAE,OMEGA1,OMEGA2,OMEGAS,                         &
            ! Beginning and end of the timestep and of the
            ! period over which cosz is integrated, and sunset
            ! - all measured in radians after local sunrise,
            ! not from local noon as the true hour angle is.
      DIFSIN,DIFTIM,                                              &
            ! A difference-of-sines intermediate value and the
            ! corresponding time period
      TRAD, DTRAD
            ! These are the start-time and length of the
            ! timestep (T & DT) converted to radians after
            ! midday GMT, or equivalently, hour angle of the
            ! sun on the Greenwich meridian.

!-----------------------------------------------------------------
! Set up parameter values
!-----------------------------------------------------------------
    PARAMETER( TWOPI = 2. * PI, S2R = PI / 43200.)



    TRAD = T * S2R - PI
    DTRAD = DT * S2R

!DIR$ IVDEP
    DO J = 1,K                          ! Loop over points

! Logically unnecessary statement without which the CRAY compiler will not vectorize this code
      HLD = 0.0

      SINSIN = SINDEC * SINLAT(J)
      COSCOS = SQRT( (1.-SINDEC**2) * (1.-SINLAT(J)**2) )
      COSHLD = SINSIN / COSCOS
      IF(COSHLD < -1.0) THEN             ! Perpetual night
        LIT(J)  = 0.0
        COSZ(J) = 0.0
      ELSE
        HAT = LONGIT(J) + TRAD           ! (3.2.2)                
        IF(COSHLD > 1.0) THEN            !   Perpetual day - hour 
          OMEGA1 = HAT                   ! angles for (3.2.3) are 
          OMEGA2 = HAT + DTRAD           ! start & end of timestep
        ELSE                                !   At this latitude some
! points are sunlit, some not.  Different ones need different treatment.
          HLD = ACOS(-COSHLD)               ! (3.2.4)
! The logic seems simplest if one takes all "times" - actually hour     
! angles - relative to sunrise (or sunset), but they must be kept in the
! range 0 to 2pi for the tests on their orders to work.                 
          OMEGAB = HAT + HLD                                         
          IF(OMEGAB < 0.0)    OMEGAB = OMEGAB + TWOPI
          IF(OMEGAB >= TWOPI) OMEGAB = OMEGAB - TWOPI
          IF(OMEGAB >= TWOPI) OMEGAB = OMEGAB - TWOPI
!            !  Line repeated - otherwise could have failure if         
!            !  longitudes W are > pi rather than < 0.                  
          OMEGAE = OMEGAB + DTRAD                                    
          IF(OMEGAE > TWOPI)  OMEGAE = OMEGAE - TWOPI
          OMEGAS = 2. * HLD

! Now that the start-time, end-time and sunset are set in terms of hour 
! angle, can set the two hour-angles for (3.2.3).  The simple cases are 
! start-to-end-of-timestep, start-to-sunset, sunrise-to-end and sunrise-
! -to-sunset, but two other cases exist and need special treatment.     
          IF(OMEGAB <= OMEGAS .OR. OMEGAB < OMEGAE) THEN
            OMEGA1 = OMEGAB - HLD
          ELSE
            OMEGA1 = - HLD
          ENDIF

          IF(OMEGAE <= OMEGAS) THEN
            OMEGA2 = OMEGAE - HLD
          ELSE
            OMEGA2 = OMEGAS - HLD
          ENDIF

!  Put in an arbitrary marker for the case when the sun does not rise   
!  during the timestep (though it is up elsewhere at this latitude).    
!  (Cannot set COSZ & LIT within the ELSE ( COSHLD < 1 ) block          
!  because 3.2.3 is done outside this block.)                           
          IF(OMEGAE > OMEGAB .AND. OMEGAB > OMEGAS)               &
            OMEGA2=OMEGA1

        ENDIF           ! This finishes the ELSE (perpetual day) block

        DIFSIN = SIN(OMEGA2) - SIN(OMEGA1)             ! Begin (3.2.3)
        DIFTIM = OMEGA2 - OMEGA1

! Next, deal with the case where the sun sets and then rises again      
! within the timestep.  There the integration has actually been done    
! backwards over the night, and the resulting negative DIFSIN and DIFTIM
! must be combined with positive values representing the whole of the   
! timestep to get the right answer, which combines contributions from   
! the two separate daylit periods.  A simple analytic expression for the
! total sun throughout the day is used.  (This could of course be used  
! alone at points where the sun rises and then sets within the timestep)
        IF(DIFTIM < 0.0) THEN
          DIFSIN = DIFSIN + 2. * SQRT(1.-COSHLD**2)                   
          DIFTIM = DIFTIM + 2. * HLD                                  
        ENDIF                                                         

        IF (DIFTIM == 0.0) THEN
! Pick up the arbitrary marker for night points at a partly-lit latitude
          COSZ(J) = 0.0
          LIT(J)  = 0.0
        ELSE
          COSZ(J) = DIFSIN*COSCOS/DIFTIM + SINSIN     ! (3.2.3)      
          LIT(J) = DIFTIM / DTRAD
        ENDIF
      ENDIF            ! This finishes the ELSE (perpetual night) block
    ENDDO

    RETURN                                                            

 END SUBROUTINE SOLANG
#endif
