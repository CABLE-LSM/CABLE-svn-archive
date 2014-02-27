#if !defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine CALC_FIT_FSAT-------------------------------------------
!
!   Purpose: To speed up the large scale hydrology code (LTOP=TRUE)
!            dramatically. This is done by fitting exponential
!            functions to the incomplete gamma function for each grid
!            box and the complete range of possible "water table"
!            (top_crit) cases - see documentation.
!            Estimates the fitted parameters for Fsat=function(ZW)
!            and  Fwet=function(ZW) for each land grid point.
!            (Calculating the incomplete gamma function for each grid
!            box at each time step was very time consuming).
!                                                             !
! Documentation: UNIFIED MODEL DOCUMENTATION PAPER NO 25
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  6.4    10/01/07   New Deck         Nic Gedney
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.

SUBROUTINE CALC_FIT_FSAT()

  USE ancil_info, ONLY : land_pts, sm_levels, soil_pts, soil_index

  USE top_pdm, ONLY : a_fsat, a_fwet, c_fsat, c_fwet, fexp, gamtot, ti_mean,  &
                      ti_sig

  USE c_topog, ONLY : zw_max

  USE soil_param, ONLY : dzsoil

  USE SWITCHES, ONLY : l_top

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Subroutine arguments - NONE.
!-----------------------------------------------------------------------------

! Work variables
  REAL :: ZDEPTH           ! Standard Soil model DEPTH.

! Local scalars:
  INTEGER, PARAMETER :: NZW = 100  ! Number of ZW values used in fit.
                                   ! Maximum value for a significant
                                   ! improvement in the fit.

  INTEGER :: I,J,IZ,N,IFITA  ! Loop counters.
  INTEGER, PARAMETER :: NFITA = 20  ! Number of loops for fitting.

  REAL ::                                                                     &
    DZW,                                                                      &
                   ! WORK ZW increment ; defined by ZW_MAX and NZW.
    RMS,                                                                      &
                   ! WORK RMS errors for given fsat fit values.
    RMSW,                                                                     &
                   ! WORK RMS errors for given fwet fit values.
    RMSOLD,                                                                   &
                   ! WORK RMS errors for given fsat fit values.
!                  !      for best fit so far.
    RMSWOLD,                                                                  &
                   ! WORK RMS errors for given fwet fit values
!                  !      for best fit so far.
     CFIT
                   ! WORK CFit value for given loop.

  REAL, PARAMETER ::                                                          &
    THR_ERR = 5.0E-3,                                                         &
                   ! Error threshold value

    CFITMIN = 0.0,                                                            &
                   ! Minimum possible value for Cfit.
    CFITMAX = 3.0
                   ! Maximum possible value for Cfit.

! Local arrays:
  REAL ::                                                                     &
    FSAT_CALC(LAND_PTS,NZW),                                                  &
                              ! WORK Surface saturation fraction.
    FSAT_FIT(NZW),                                                            &
                              ! WORK Fitted surface saturation fraction.
    FWET_CALC(LAND_PTS,NZW),                                                  &
                              ! WORK Wetland fraction.
    FWET_FIT(NZW),                                                            &
                              ! WORK Fitted wetland fraction.
    DUMZW(NZW),                                                               &
                              ! WORK Dummy water table depth (m).
    DUMFSAT(LAND_PTS),                                                        &
                              ! WORK Dummy surface saturation fraction.
    DUMFWETL(LAND_PTS)
                              ! WORK Dummy wetland fraction.
! DBC We could use this local version of gamtot in this version,
!     but using the module version to match other versions of code.
!     &,GAMTOT(LAND_PTS)       ! WORK Integrated complete Gamma function

  REAL ::                                                                     &
    TOP_CRIT(LAND_PTS,NZW),                                                   &
                              ! WORK LOG(QBASE_MAX/QBASE) -see document.
    TOP_CRIT1Z(LAND_PTS),                                                     &
                              ! WORK As above but for an individual zw.
    TOP_MIN(LAND_PTS),                                                        &
                              ! WORK value for when zw=zw_max.
    WUTOT(LAND_PTS)           ! WORK Dummy (set to 1.0).

  INTEGER :: ERRORSTATUS
  CHARACTER (LEN=80) :: CMESSAGE
  CHARACTER (LEN=13) :: ROUTINENAME


!-----------------------------------------------------------------------------


  ROUTINENAME='CALC_FIT_FSAT'

  IF ( L_TOP ) THEN

! Define the water table depths to be used in the fitting process:
    DZW = 1.0 / REAL(NZW) * ZW_MAX
    DO IZ = 1,NZW
      DUMZW(IZ) = REAL(IZ - 1) * DZW
    END DO

    DO I = 1,LAND_PTS        ! initialise to zero
      WUTOT(I)    = 1.0
      DUMFSAT(I)  = 0.0
      DUMFWETL(I) = 0.0
    END DO
    ZDEPTH = 0
         
! Calculate total soil depth
    DO N = 1,SM_LEVELS
      ZDEPTH = ZDEPTH + DZSOIL(N)
    END DO
         
! Calculate Gamtot
    gamtot(:) = 0.0
! DEPENDS ON: calc_fsat
    CALL CALC_FSAT(.TRUE., SOIL_PTS, SOIL_INDEX, LAND_PTS, TI_MEAN, TI_SIG,   &
                   WUTOT, TOP_CRIT1Z, GAMTOT, DUMFSAT, DUMFWETL)
         
! Calculate TOP_CRIT for the water table depths:
    DO IZ = 1,NZW
      DO J = 1,SOIL_PTS
        I = SOIL_INDEX(J)
               
        FSAT_CALC(I,IZ) = 0.0
        FWET_CALC(I,IZ) = 0.0
        TOP_CRIT(I,IZ)  = 0.0
               
        IF ( TI_MEAN(I) > 0.0 .AND. TI_SIG(I) > 0.0 ) THEN
          TOP_MIN(I) = 1.0 / FEXP(I) * EXP(-FEXP(I) * (ZW_MAX - ZDEPTH))
                  
          IF ( DUMZW(IZ) <= ZDEPTH )                                          &
            TOP_CRIT1Z(I) = LOG(ZDEPTH + 1.0 / FEXP(I) - TOP_MIN(I))          &
                          - LOG(ZDEPTH - DUMZW(IZ) + 1.0 / FEXP(I) -          &
                                TOP_MIN(I))

          IF ( DUMZW(IZ) >  ZDEPTH )                                          &
            TOP_CRIT1Z(I) = LOG(ZDEPTH + 1.0 / FEXP(I) - TOP_MIN(I))          &
                          - LOG(1 / FEXP(I) *                                 &
                                EXP(-FEXP(I) * (DUMZW(IZ) - ZDEPTH)) -        &
                                TOP_MIN(I))
        ENDIF
      ENDDO
            
! Calculate FSAT and FWET for one ZW at all soil land_pts:
! DEPENDS ON: calc_fsat
      CALL CALC_FSAT(.FALSE., SOIL_PTS, SOIL_INDEX, LAND_PTS, TI_MEAN,        &
                     TI_SIG, WUTOT, TOP_CRIT1Z, GAMTOT, DUMFSAT, DUMFWETL)
            
      DO J = 1,SOIL_PTS
        I = SOIL_INDEX(J)

        IF ( TI_MEAN(I) > 0.0 .AND. TI_SIG(I) > 0.0 ) THEN
          FSAT_CALC(I,IZ) = DUMFSAT(I)
          FWET_CALC(I,IZ) = DUMFWETL(I)
          TOP_CRIT(I,IZ)  = TOP_CRIT1Z(I)

          IF ( DUMZW(IZ) <  DZW ) THEN ! Values at zw=0m
            A_FSAT(I) = FSAT_CALC(I,IZ)
            A_FWET(I) = FWET_CALC(I,IZ)
          ENDIF
        ENDIF
      ENDDO
    ENDDO                  !ZW calc_fsat loop
         
         
! Now carry out fit for FSAT, where FSAT=function(ZW). (Likewise FWET)
    DO J = 1,SOIL_PTS
      I = SOIL_INDEX(J)

      IF ( TI_MEAN(I) > 0.0 .AND. TI_SIG(I) > 0.0 ) THEN
        RMSOLD = 1.0E10
        RMSWOLD = 1.0E10
               
        DO IFITA = 1,NFITA
          CFIT = CFITMAX * (IFITA) / FLOAT(NFITA)
                  
          RMS  = 0.0
          RMSW = 0.0
! TOP_CRIT=TI_MAX when zw=zw_max
          DO IZ = 1,NZW
            FSAT_FIT(IZ) = A_FSAT(I) * EXP(-CFIT * TOP_CRIT(I,IZ))
            FWET_FIT(IZ) = A_FWET(I) * EXP(-CFIT * TOP_CRIT(I,IZ))
            RMS  = RMS + (FSAT_CALC(I,IZ) - FSAT_FIT(IZ))**2
            RMSW = RMSW + (FWET_CALC(I,IZ) - FWET_FIT(IZ))**2
          ENDDO         !ZW

          RMS  = SQRT(RMS) / FLOAT(NZW)
          RMSW = SQRT(RMSW) / FLOAT(NZW)
                  
          IF ( RMS < RMSOLD ) THEN
            RMSOLD = RMS
            C_FSAT(I) = CFIT
          ENDIF
          IF ( RMSW < RMSWOLD ) THEN
            RMSWOLD = RMSW
            C_FWET(I) = CFIT
          ENDIF
        ENDDO
               
        DO IZ = 1,NZW
          FSAT_FIT(IZ) = A_FSAT(I) * EXP(-C_FSAT(I) * TOP_CRIT(I,IZ))
          FWET_FIT(IZ) = A_FWET(I) * EXP(-C_FWET(I) * TOP_CRIT(I,IZ))
        ENDDO            !ZW
               
        IF ( RMSOLD >= THR_ERR ) THEN
          IF ( C_FSAT(I) <= CFITMIN .OR. C_FSAT(I) >= CFITMAX ) THEN
            WRITE(6,*) 'ERROR CFIT FSAT', I, C_FSAT(I), CFITMIN, CFITMAX

            WRITE(6,*) 'fsat_calc=',                                          &
                        FSAT_CALC(I,1), FSAT_CALC(I,3), FSAT_CALC(I,5)
            WRITE(6,*) 'fsat_fit=',                                           &
                        FSAT_FIT(1), FSAT_FIT(3), FSAT_FIT(5)
            WRITE(6,*) 'RMS=', RMSOLD
            ErrorStatus = 35
            WRITE(CMessage, *) 'Error in CFIT FSAT in LSH model setup'
          ENDIF
        ENDIF
               
        IF ( RMSWOLD >= THR_ERR ) THEN
          IF ( C_FWET(I) <= CFITMIN .OR. C_FWET(I) >= CFITMAX ) THEN
            WRITE(6,*) 'ERROR CFIT FWET', I, C_FWET(I), CFITMIN, CFITMAX

            WRITE(6,*) 'fwet_calc=',                                          &
                        FWET_CALC(I,1), FWET_CALC(I,3), FWET_CALC(I,5)
            WRITE(6,*) 'fwet_fit=',                                           &
                        FWET_FIT(1), FWET_FIT(3), FWET_FIT(5)
            WRITE(6,*) 'RMSW=', RMSWOLD
            WRITE(6,*) '(fsat_calc=)',                                        &
                         FSAT_CALC(I,1), FSAT_CALC(I,3), FSAT_CALC(I,5)
            WRITE(6,*) '(fsat_fit=)',                                         &
                         FSAT_FIT(1), FSAT_FIT(3), FSAT_FIT(5)
            WRITE(6,*) '(RMS=)', RMSOLD
            ErrorStatus = 35
            ErrorStatus = 40
            WRITE(CMessage,*) 'Error in CFIT FWET in LSH model setup'
          ENDIF
        ENDIF

        IF( RMSOLD > THR_ERR )                                                &
          WRITE(6,*) 'Warning LSH RMS Error in fit:', RMSOLD, RMSWOLD
               
      ENDIF
            
    ENDDO                  ! LAND_PTS
  ELSE
    DO J = 1,SOIL_PTS
      I = SOIL_INDEX(J)

      A_FSAT(I) = 0.0
      C_FSAT(I) = 0.0
      A_FWET(I) = 0.0
      C_FWET(I) = 0.0
      GAMTOT(I) = 0.0
    ENDDO
  ENDIF

  RETURN

END SUBROUTINE CALC_FIT_FSAT
#endif
